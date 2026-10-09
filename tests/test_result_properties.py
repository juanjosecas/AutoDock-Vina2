"""Reports tested without loading or building the native Vina extension."""
import csv
import importlib.util
import math
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'build/python/vina'
spec = importlib.util.spec_from_file_location('result_properties', PACKAGE / 'results.py')
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def pdbqt(types_, smiles=None):
    text = '\n'.join('ATOM  %5d C LIG 1 0.0 0.0 0.0 0.0 %s' % (i, t)
                     for i, t in enumerate(types_, 1))
    return text + '\nTORSDOF 1\nBRANCH 1 2\n' + ('REMARK SMILES %s\n' % smiles if smiles else '')


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.input = pdbqt(['C', 'C', 'OA', 'HD', 'CG0', 'G0', 'W'])
        self.energies = [[-7.123456789, -8, 2, 0.876543211, 2], [-6, -7, 3, 0, 2]]
        self.details = [[-7.123456789, -6, -2, 0, 0, 2, 0.876543211, 2, 0, 0],
                        [-6, -5, -2, 0, 0, 3, 0, 2, 1.2, 1.5]]

    def rows(self, **kwargs):
        return report.result_rows(self.energies, self.details, [self.input], 'vina', **kwargs)

    def test_count_real_atoms_and_missing_descriptors(self):
        rows = self.rows()
        self.assertEqual(rows[0]['heavy_atom_count'], 4)
        self.assertIsNone(rows[0]['clogp'])
        self.assertEqual(rows[0]['docking_branches'], 1)
        self.assertEqual(rows[0]['score'], -7.123456789)
        self.assertAlmostEqual(rows[0]['docking_efficiency'], 7.123456789/4)
        self.assertAlmostEqual(rows[1]['delta_score'], 1.123456789)
        self.assertEqual(rows[1]['intra_delta'], 1)
        for row in rows:
            self.assertAlmostEqual(row['score'], row['energy_inter'] + row['energy_intra']
                                   + row['torsion_term'] - row['unbound_reference'])

    def test_experimental_efficiency_and_validation(self):
        row = self.rows(experimental_pvalue=7, activity_type='Kd')[0]
        self.assertAlmostEqual(row['experimental_le'], report.R_KCAL*298.15*math.log(10)*7/4)
        self.assertIsNone(row['experimental_lle'])
        for kwargs in ({'temperature': 0}, {'temperature': float('nan')},
                       {'experimental_pvalue': 7}, {'activity_type': 'Ki'}):
            with self.assertRaises(ValueError):
                self.rows(**kwargs)

    def test_ad4_reference_sign(self):
        row = report.result_rows([[-7, -8, 3, 1, 3]],
                                 [[-7, -8, 0, 0, 0, 3, 1, 3, 0, 0]],
                                 [self.input], 'ad4')[0]
        self.assertEqual(row['unbound_reference'], 3)
        self.assertEqual(row['score'], -7)

    def test_rdkit_properties_and_mismatch(self):
        try:
            import rdkit
        except ImportError:
            self.skipTest('Optional RDKit unavailable')
        props = report.ligand_properties([pdbqt(['C', 'C', 'OA'], 'CCO')])
        self.assertAlmostEqual(props['molecular_weight'], 46.069, places=3)
        self.assertEqual(props['hbd'], 1)
        self.assertEqual(props['hba'], 1)
        self.assertEqual(props['lipinski_violations'], 0)
        self.assertEqual(props['chemical_source'], 'pdbqt_smiles')
        with self.assertRaises(ValueError):
            report.ligand_properties([self.input], 'CCO')
        rows = report.result_rows(self.energies, self.details, [pdbqt(['C','C','OA'], 'CCO')],
                                  'vina', experimental_pvalue=7, activity_type='IC50')
        self.assertAlmostEqual(rows[0]['experimental_lle'], 7 - props['clogp'])

    def test_python_api_export_without_native_extension(self):
        package = types.ModuleType('report_test_vina')
        package.__path__ = [str(PACKAGE)]
        native = types.ModuleType('report_test_vina.vina_wrapper')
        native.Vina = object
        with patch.dict(sys.modules, {'report_test_vina': package,
                                     'report_test_vina.vina_wrapper': native}):
            spec = importlib.util.spec_from_file_location('report_test_vina.vina', PACKAGE/'vina.py')
            wrapper = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(wrapper)
            instance = wrapper.Vina.__new__(wrapper.Vina)
            instance._ligand_pdbqt = [self.input]
            instance._sf_name = 'vina'
            instance._vina = types.SimpleNamespace(get_poses_energies=lambda *args: self.energies,
                                                   get_poses_components=lambda *args: self.details)
            self.assertEqual(instance.results()[0]['score'], -7.123456789)
            with tempfile.TemporaryDirectory() as directory:
                file = Path(directory)/'results.tsv'
                instance.write_results(file)
                with file.open() as stream:
                    rows = list(csv.DictReader(stream, delimiter='\t'))
                self.assertEqual(rows[0]['score'], '-7.123457')
                self.assertEqual(rows[0]['clogp'], '')
                with self.assertRaises(FileExistsError):
                    instance.write_results(file)


if __name__ == '__main__':
    unittest.main()
