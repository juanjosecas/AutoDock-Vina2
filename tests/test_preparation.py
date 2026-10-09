"""Real RDKit/Meeko preparation tests; native Vina is stubbed, never built."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'build/python/vina'
spec = importlib.util.spec_from_file_location('prep_under_test', PACKAGE / 'preparation.py')
prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prep)
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import meeko
    AVAILABLE = meeko.__version__ == '0.7.1'
except ImportError:
    AVAILABLE = False


@unittest.skipUnless(AVAILABLE, 'RDKit and Meeko 0.7.1 required')
class PreparationTests(unittest.TestCase):
    def peptide(self, directory):
        mol = Chem.AddHs(Chem.MolFromFASTA('AG'))
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=42), 0)
        filename = Path(directory)/'peptide.pdb'
        filename.write_text(Chem.MolToPDBBlock(Chem.RemoveHs(mol)))
        return filename

    def test_smiles_reproducible_charge_stereo_and_messages(self):
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            first = prep.prepare_ligand('C[C@H](O)C(=O)[O-]', input_format='smiles', verbosity=2)
        second = prep.prepare_ligand('C[C@H](O)C(=O)[O-]', input_format='smiles', verbosity=0)
        self.assertEqual(first['pdbqt_string'], second['pdbqt_string'])
        self.assertEqual(first['report']['formal_charge'], -1)
        self.assertIn('@', first['report']['input_smiles'])
        self.assertEqual(first['report']['optimization']['status'], 0)
        self.assertIn('Gasteiger', stream.getvalue())
        self.assertIn('no pH optimization', stream.getvalue())
        self.assertIn('REMARK SMILES', first['pdbqt_string'])

    def test_existing_3d_coordinates_and_output_protection(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
        AllChem.EmbedMolecule(mol, randomSeed=17)
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'ethanol.sdf'
            with Chem.SDWriter(str(source)) as writer:
                writer.write(mol)
            output = Path(directory)/'ethanol.pdbqt'
            result = prep.prepare_ligand(source, output, verbosity=0)
            self.assertFalse(result['report']['generated_3d'])
            records = [line for line in result['pdbqt_string'].splitlines() if line.startswith('ATOM')]
            # Meeko preserves heavy-atom coordinates, subject to PDBQT rounding.
            position = mol.GetConformer().GetAtomPosition(0)
            xyz = [float(records[0][start:start+8]) for start in (30,38,46)]
            for actual, expected in zip(xyz, position):
                self.assertLessEqual(abs(actual-expected), 0.00051)
            self.assertEqual(json.loads(output.with_suffix('.preparation.json').read_text())['heavy_atoms'], 3)
            with self.assertRaises(FileExistsError):
                prep.prepare_ligand(source, output, verbosity=0)

    def test_invalid_fragments_multirecord_and_geometry_failures(self):
        with self.assertRaises(ValueError):
            prep.prepare_ligand('CCO.[Na+]', input_format='smiles', verbosity=0)
        with self.assertRaises(ValueError):
            prep.prepare_ligand('not a molecule', input_format='smiles', verbosity=0)
        with patch.object(AllChem, 'EmbedMolecule', return_value=-1):
            with self.assertRaisesRegex(RuntimeError, '3D conformer'):
                prep.prepare_ligand('CCO', input_format='smiles', verbosity=0)
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'multi.sdf'
            with Chem.SDWriter(str(source)) as writer:
                writer.write(Chem.MolFromSmiles('CCO'))
                writer.write(Chem.MolFromSmiles('CC'))
            with self.assertRaisesRegex(ValueError, 'exactly one'):
                prep.prepare_ligand(source, verbosity=0)

    def test_real_receptor_templates_hydrogens_and_saved_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.peptide(directory)
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                result = prep.prepare_receptor(source, set_template='A:1=NALA', verbosity=2)
            report = result['report']
            self.assertEqual(report['residue_templates']['A:1'], 'NALA')
            self.assertEqual(report['residue_templates']['A:2'], 'CGLY')
            self.assertGreater(report['residue_changes']['A:1']['prepared']['hydrogens'], 0)
            self.assertEqual(report['residue_changes']['A:1']['input']['hydrogens'], 0)
            self.assertIn('no pKa prediction', stream.getvalue())
            self.assertIn('Template assignments', stream.getvalue())
            output = Path(result['pdbqt_filename'])
            saved = json.loads(output.with_suffix('.preparation.json').read_text())
            self.assertEqual(saved['residue_templates'], report['residue_templates'])
            self.assertEqual(saved['input_sha256'], report['input_sha256'])
            with self.assertRaises(FileExistsError):
                prep.prepare_receptor(source, verbosity=0)

    def test_incomplete_receptor_fails_without_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.peptide(directory)
            lines = source.read_text().splitlines()
            source.write_text('\n'.join(line for line in lines
                                        if not (line.startswith('ATOM') and line[12:16].strip() == 'CA')) + '\n')
            output = Path(directory)/'failed.pdbqt'
            with self.assertRaisesRegex(RuntimeError, 'Meeko receptor preparation failed'):
                prep.prepare_receptor(source, output, verbosity=0)
            self.assertFalse(output.exists())
            self.assertFalse(output.with_suffix('.preparation.json').exists())

    def test_api_automatic_loading_without_native_extension(self):
        package = types.ModuleType('prep_test_vina')
        package.__path__ = [str(PACKAGE)]
        native = types.ModuleType('prep_test_vina.vina_wrapper')
        native.Vina = lambda *args: Mock(seed=Mock(return_value=42))
        with patch.dict(sys.modules, {'prep_test_vina': package, 'prep_test_vina.vina_wrapper': native}):
            spec = importlib.util.spec_from_file_location('prep_test_vina.vina', PACKAGE/'vina.py')
            wrapper = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(wrapper)
            instance = wrapper.Vina(verbosity=0)
            report = instance.set_ligand_from_smiles('CCO')
            self.assertEqual(report['heavy_atoms'], 3)
            instance._vina.set_ligand_from_string.assert_called_once()
            with tempfile.TemporaryDirectory() as directory:
                source = self.peptide(directory)
                instance.set_receptor(source, set_template='A:1=NALA')
                instance._vina.set_receptor.assert_called_once_with(str(source.with_suffix('.prepared.pdbqt')))
                ligand = Path(directory)/'ligand.sdf'
                with Chem.SDWriter(str(ligand)) as writer:
                    writer.write(Chem.MolFromSmiles('CCO'))
                instance.set_ligand_from_file(ligand)
                self.assertEqual(instance._vina.set_ligand_from_string.call_count, 2)
            with tempfile.TemporaryDirectory() as directory:
                existing = Path(directory)/'existing.pdbqt'
                existing.write_text(instance._ligand_pdbqt[0])
                with patch.object(prep, '_dependencies', side_effect=AssertionError('Unexpected dependency import')):
                    with patch.dict(sys.modules, {'prep_test_vina.preparation': prep}):
                        instance.set_ligand_from_file(existing)
                        instance.set_receptor(existing)
                instance._vina.set_ligand_from_file.assert_called_once_with(str(existing))
            self.assertEqual(len(instance.preparation_reports()), 3)
            copy = instance.preparation_reports()
            copy[0]['heavy_atoms'] = 99
            self.assertEqual(instance.preparation_reports()[0]['heavy_atoms'], 3)


if __name__ == '__main__':
    unittest.main()
