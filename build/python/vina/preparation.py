"""Optional, auditable Meeko preparation. No native Vina dependency here."""
import argparse
from collections import Counter
from contextlib import redirect_stdout, redirect_stderr
import io
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def _dependencies():
    if sys.version_info < (3, 10):
        raise RuntimeError("Automatic preparation requires Python 3.10 or newer.")
    try:
        import meeko
        from rdkit import Chem
    except ImportError as error:
        raise ImportError("Preparation requires: pip install 'vina[preparation]' "
                          "(Meeko 0.7.1, RDKit, gemmi and scipy).") from error
    if meeko.__version__ != '0.7.1':
        raise RuntimeError('Preparation supports Meeko 0.7.1; found %s.' % meeko.__version__)
    return meeko, Chem


def _message(report, text, verbosity, detailed=False):
    report['messages'].append(text)
    if verbosity > (1 if detailed else 0):
        print('[preparation] ' + text, flush=True)


def _stats(text):
    atoms = [line for line in text.splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    types = Counter(line.split()[-1] for line in atoms)
    return {'atom_records': len(atoms), 'atom_types': dict(sorted(types.items())),
            'partial_charge_sum': sum(float(line.split()[-2]) for line in atoms),
            'torsdof': next((int(line.split()[1]) for line in text.splitlines()
                            if line.startswith('TORSDOF')), None)}


def _outputs(filename, overwrite):
    if filename is None:
        return None, None
    path = Path(filename).absolute()
    if path.suffix.lower() != '.pdbqt':
        raise ValueError('Prepared output must have a .pdbqt extension.')
    sidecar = path.with_suffix('.preparation.json')
    if not overwrite:
        for candidate in (path, sidecar):
            if candidate.exists():
                raise FileExistsError('Refusing to overwrite %s.' % candidate)
    return path, sidecar


def _save(text, report, path, sidecar, overwrite):
    if path is not None:
        # Exclusive creation also protects against another process creating the output.
        with path.open('w' if overwrite else 'x') as stream:
            stream.write(text)
        with sidecar.open('w' if overwrite else 'x') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')


def prepare_ligand(source, output_filename=None, input_format=None, seed=42,
                   verbosity=1, overwrite=False):
    """Prepare one SDF/MOL/MOL2 file, SMILES (input_format='smiles'), or RDKit Mol.

    Preserve formal charge/stereochemistry; add H and generate 3D only if needed.
    Return a dictionary with pdbqt_string and report. No pH/tautomer enumeration.
    """
    path, sidecar = _outputs(output_filename, overwrite)
    meeko, Chem = _dependencies()
    from rdkit.Chem import AllChem
    from rdkit import __version__ as rdkit_version
    report = {'kind': 'ligand', 'meeko_version': meeko.__version__, 'messages': [],
              'output': str(path) if path else None, 'charge_model': 'gasteiger',
              'ph_selection': None, 'seed': seed, 'rdkit_version': rdkit_version}
    if input_format == 'smiles':
        mol = Chem.MolFromSmiles(source)
        report['input'] = source
    elif isinstance(source, Chem.Mol):
        mol = Chem.Mol(source)
        report['input'] = 'RDKit Mol'
    else:
        filename = Path(source)
        report['input'] = str(filename.absolute())
        report['input_sha256'] = hashlib.sha256(filename.read_bytes()).hexdigest()
        suffix = filename.suffix.lower()
        if suffix == '.sdf':
            supplier = Chem.SDMolSupplier(str(filename), removeHs=False)
            if len(supplier) != 1:
                raise ValueError('Supply exactly one molecule per SDF; multiple records are not docked implicitly.')
            mol = supplier[0]
        elif suffix == '.mol':
            mol = Chem.MolFromMolFile(str(filename), removeHs=False)
        elif suffix == '.mol2':
            mol = Chem.MolFromMol2File(str(filename), removeHs=False)
        else:
            raise ValueError('Ligand preparation accepts SDF, MOL, MOL2, SMILES or RDKit Mol; use SDF instead of PDB.')
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError('Could not read a valid ligand.')
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('Disconnected ligand fragments: select the intended compound/salt explicitly.')
    report['input_smiles'] = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    report['input_atoms'] = mol.GetNumAtoms()
    report['heavy_atoms'] = mol.GetNumHeavyAtoms()
    report['formal_charge'] = Chem.GetFormalCharge(mol)
    _message(report, 'Ligand: %s | heavy atoms: %d | formal charge: %+d' %
             (report['input'], report['heavy_atoms'], report['formal_charge']), verbosity)
    _message(report, 'Protonation/tautomer state is taken from the input; no pH optimization or enumeration.', verbosity)
    has_3d = mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()
    before = mol.GetNumAtoms()
    mol = Chem.AddHs(mol, addCoords=has_3d)
    report['hydrogens_added'] = mol.GetNumAtoms() - before
    _message(report, 'Added %d explicit hydrogens.' % report['hydrogens_added'], verbosity)
    report['generated_3d'] = not has_3d
    report['optimization'] = None
    if not has_3d:
        mol.RemoveAllConformers()
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = seed
        if AllChem.EmbedMolecule(mol, parameters) != 0:
            raise RuntimeError('RDKit could not generate a 3D conformer.')
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            raise ValueError('MMFF94 parameters unavailable; supply a prepared 3D structure.')
        status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        report['optimization'] = {'force_field': 'MMFF94', 'max_iterations': 500, 'status': status}
        _message(report, 'Generated 3D with ETKDGv3 (seed %d); MMFF94 optimization status: %d (0=converged).' %
                 (seed, status), verbosity)
        if status != 0:
            raise RuntimeError('MMFF94 optimization did not converge; supply a prepared 3D structure.')
    else:
        _message(report, 'Preserved input 3D heavy-atom coordinates; no geometry minimization.', verbosity)
    report['unassigned_stereocenters'] = [index for index, tag in
                                        Chem.FindMolChiralCenters(mol, includeUnassigned=True) if tag == '?']
    if report['unassigned_stereocenters']:
        _message(report, 'WARNING: unassigned stereocenters: %s; generated geometry is not a stereochemical assignment.' %
                 report['unassigned_stereocenters'], max(verbosity, 1))
    log = io.StringIO()
    try:
        with redirect_stdout(log), redirect_stderr(log):
            setups = meeko.MoleculePreparation(charge_model='gasteiger').prepare(mol)
            if len(setups) != 1:
                raise ValueError('Meeko produced multiple setups; choose one explicitly.')
            text, success, error = meeko.PDBQTWriterLegacy.write_string(setups[0])
        if not success:
            raise RuntimeError(error)
    except Exception as error:
        raise RuntimeError('Meeko ligand preparation failed: %s\n%s' % (error, log.getvalue())) from error
    report['meeko_log'] = log.getvalue()
    if report['meeko_log']:
        _message(report, report['meeko_log'].rstrip(), verbosity)
    report['pdbqt'] = _stats(text)
    _message(report, 'Gasteiger charges; Meeko atom typing; nonpolar H merged into parent atoms. '
             'PDBQT records: %d | TORSDOF: %s | charge sum: %.6f' %
             (report['pdbqt']['atom_records'], report['pdbqt']['torsdof'], report['pdbqt']['partial_charge_sum']), verbosity)
    _message(report, 'PDBQT atom types: %s' % report['pdbqt']['atom_types'], verbosity, detailed=True)
    _save(text, report, path, sidecar, overwrite)
    if path:
        _message(report, 'Saved %s and %s' % (path, sidecar), verbosity)
    return {'pdbqt_string': text, 'report': report}


def prepare_receptor(source, output_filename=None, delete_residues=None,
                     set_template=None, verbosity=1, overwrite=False):
    """Prepare a rigid PDB/PQR receptor via Meeko's CLI; never silently drop bad residues.

    set_template uses Meeko syntax, e.g. 'A:42=HID'. No pKa calculation is run.
    Returns pdbqt_filename and an auditable report, also saved beside the PDBQT.
    """
    filename = Path(source).absolute()
    if filename.suffix.lower() not in ('.pdb', '.pqr'):
        raise ValueError('Receptor preparation accepts PDB/PQR; existing PDBQT needs no conversion.')
    if not filename.is_file():
        raise FileNotFoundError(str(filename))
    if output_filename is None:
        output_filename = filename.with_name(filename.stem + '.prepared.pdbqt')
    path, sidecar = _outputs(output_filename, overwrite)
    meeko, Chem = _dependencies()
    from rdkit.Chem import rdMolInterchange
    from rdkit import __version__ as rdkit_version
    report = {'kind': 'receptor', 'input': str(filename), 'output': str(path),
              'meeko_version': meeko.__version__, 'messages': [], 'charge_model': 'gasteiger',
              'delete_residues': delete_residues, 'set_template': set_template,
              'ph_selection': None, 'rdkit_version': rdkit_version,
              'input_sha256': hashlib.sha256(filename.read_bytes()).hexdigest()}
    records = [line for line in filename.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    report['input_atom_records'] = len(records)
    _message(report, 'Receptor: %s | input atom records: %d' % (filename, len(records)), verbosity)
    _message(report, 'Residue deletion requested: %s | template overrides: %s' %
             (delete_residues or 'none', set_template or 'none'), verbosity)
    _message(report, 'Meeko matches residue templates and adds missing H; no pKa prediction or pH optimization. '
             'Incomplete/unsupported residues cause an error; no automatic water/cofactor removal.', verbosity)
    with tempfile.TemporaryDirectory(prefix='vina-receptor-') as directory:
        base = Path(directory) / 'receptor'
        command = [sys.executable, '-u', '-m', 'meeko.cli.mk_prepare_receptor',
                   '--read_pqr' if filename.suffix.lower() == '.pqr' else '--read_pdb',
                   str(filename), '-o', str(base), '-p', '-j', '--charge_model', 'gasteiger']
        if delete_residues:
            command.extend(['--delete_residues', delete_residues])
        if set_template:
            command.extend(['--set_template', set_template])
        report['command'] = command
        _message(report, 'Running Meeko receptor preparation...', verbosity)
        lines = []
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              universal_newlines=True) as process:
            for line in process.stdout:
                lines.append(line)
                if verbosity > 0:
                    print('[meeko] ' + line.rstrip(), flush=True)
            code = process.wait()
        report['meeko_log'] = ''.join(lines)
        if code != 0 or not base.with_suffix('.pdbqt').is_file():
            raise RuntimeError('Meeko receptor preparation failed (exit %s).\n%s' %
                               (code, report['meeko_log']))
        text = base.with_suffix('.pdbqt').read_text()
        parameterization = json.loads(base.with_suffix('.json').read_text())
        report['residue_templates'] = {key: value.get('residue_template_key')
                                      for key, value in parameterization['monomers'].items()}
        report['template_matching_log'] = parameterization['log']
        report['residue_changes'] = {}
        for key, value in parameterization['monomers'].items():
            changes = {}
            for name in ('raw_rdkit_mol', 'rdkit_mol'):
                molecule = rdMolInterchange.JSONToMols(value[name])[0] if value[name] else None
                changes['input' if name == 'raw_rdkit_mol' else 'prepared'] = (
                    {'atoms': molecule.GetNumAtoms(), 'heavy_atoms': molecule.GetNumHeavyAtoms(),
                     'hydrogens': sum(atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms()),
                     'formal_charge': Chem.GetFormalCharge(molecule)} if molecule else None)
            report['residue_changes'][key] = changes
    report['pdbqt'] = _stats(text)
    if not report['pdbqt']['atom_records']:
        raise RuntimeError('Meeko produced an empty receptor.')
    _message(report, 'Matched residue templates: %s' % dict(Counter(report['residue_templates'].values())), verbosity)
    _message(report, 'Template assignments and changes (input -> prepared):', verbosity, detailed=True)
    _message(report, '%-12s %-10s %-13s %-13s %s' %
             ('Residue', 'Template', 'Atoms', 'Hydrogens', 'Formal charge'), verbosity, detailed=True)
    for residue, changes in report['residue_changes'].items():
        before, after = changes['input'], changes['prepared']
        def transition(key):
            return '%s -> %s' % (before[key] if before else 'NA', after[key] if after else 'deleted')
        _message(report, '%-12s %-10s %-13s %-13s %s' %
                 (residue, report['residue_templates'][residue] or 'deleted', transition('atoms'),
                  transition('hydrogens'), transition('formal_charge')), verbosity, detailed=True)
    _message(report, 'Gasteiger charge model; PDBQT records: %d | partial-charge sum: %.6f' %
             (report['pdbqt']['atom_records'], report['pdbqt']['partial_charge_sum']), verbosity)
    _message(report, 'PDBQT atom types: %s' % report['pdbqt']['atom_types'], verbosity, detailed=True)
    _save(text, report, path, sidecar, overwrite)
    _message(report, 'Saved %s and %s' % (path, sidecar), verbosity)
    return {'pdbqt_filename': str(path), 'report': report}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=['ligand', 'receptor'])
    parser.add_argument('source')
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('--smiles', action='store_true', help='Interpret ligand source as SMILES')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--delete-residues')
    parser.add_argument('--set-template')
    parser.add_argument('--verbosity', type=int, choices=[0, 1, 2], default=1)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if args.kind == 'ligand':
        if args.delete_residues or args.set_template:
            parser.error('Residue options apply only to receptors.')
        prepare_ligand(args.source, args.output, 'smiles' if args.smiles else None,
                       args.seed, args.verbosity, args.overwrite)
    else:
        if args.smiles:
            parser.error('--smiles applies only to ligands.')
        prepare_receptor(args.source, args.output, args.delete_residues,
                         args.set_template, args.verbosity, args.overwrite)


if __name__ == '__main__':
    main()
