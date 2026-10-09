"""Post-docking reports. These calculations never alter energies or pose ranking."""

import math
from collections import Counter

R_KCAL = 8.31446261815324 / 4184.0
CHEMICAL_FIELDS = ('molecular_weight', 'clogp', 'hbd', 'hba', 'tpsa',
                   'rotatable_bonds', 'formal_charge', 'rings', 'fraction_csp3',
                   'lipinski_violations')
ATOM_ELEMENTS = {'A': 'C', 'NA': 'N', 'NS': 'N', 'OA': 'O', 'OS': 'O',
                 'SA': 'S', 'CG0': 'C', 'CG1': 'C', 'CG2': 'C', 'CG3': 'C'}
PSEUDO_TYPES = {'G0', 'G1', 'G2', 'G3', 'W'}


def pdbqt_properties(text):
    """Read real ligand atoms; CG carbons are real, G closure sites are not."""
    elements = Counter()
    branches = 0
    torsdof = None
    smiles = None
    if sum(line.startswith('MODEL') for line in text.splitlines()) > 1:
        raise ValueError('Use the prepared ligand input, not a multi-pose PDBQT.')
    for line in text.splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            fields = line.split()
            atom_type = fields[-1]
            if atom_type in PSEUDO_TYPES or atom_type in ('H', 'HD', 'HS'):
                continue
            elements[ATOM_ELEMENTS.get(atom_type, atom_type)] += 1
        elif line.startswith('BRANCH '):
            branches += 1
        elif line.startswith('TORSDOF '):
            torsdof = int(line.split()[1])
        elif line.startswith('REMARK SMILES '):
            fields = line.split()
            if len(fields) == 3:
                smiles = fields[2]
    if not elements:
        raise ValueError('No real heavy atoms found in ligand PDBQT.')
    return {'heavy_atom_count': sum(elements.values()), 'docking_branches': branches,
            'docking_torsdof': torsdof, 'elements': elements, 'smiles': smiles}


def ligand_properties(pdbqt_strings, molecule=None):
    inputs = [pdbqt_properties(text) for text in pdbqt_strings]
    if not inputs:
        raise ValueError('No ligand supplied.')
    props = {key: None for key in CHEMICAL_FIELDS}
    props.update(heavy_atom_count=sum(item['heavy_atom_count'] for item in inputs),
                 ligand_count=len(inputs),
                 docking_branches=sum(item['docking_branches'] for item in inputs),
                 docking_torsdof=(sum(item['docking_torsdof'] for item in inputs)
                                 if all(item['docking_torsdof'] is not None for item in inputs)
                                 else None),
                 chemical_source='unavailable')
    if len(inputs) != 1:
        if molecule is not None:
            raise ValueError('Chemical descriptors require single-ligand docking.')
        props['chemical_source'] = 'multiple_ligands'
        return props
    source = molecule if molecule is not None else inputs[0]['smiles']
    if source is None:
        return props
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    except ImportError:
        if molecule is not None:
            raise ImportError('Chemical descriptors require RDKit: install vina[properties].')
        props['chemical_source'] = 'rdkit_unavailable'
        return props
    mol = Chem.MolFromSmiles(source) if isinstance(source, str) else Chem.Mol(source)
    if mol is None:
        raise ValueError('Invalid ligand SMILES.')
    Chem.SanitizeMol(mol)
    elements = Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if elements != inputs[0]['elements']:
        raise ValueError('Chemical structure and PDBQT heavy-atom composition differ.')
    props.update(molecular_weight=float(Descriptors.MolWt(mol)),
                 clogp=float(Crippen.MolLogP(mol)),
                 hbd=int(Lipinski.NumHDonors(mol)), hba=int(Lipinski.NumHAcceptors(mol)),
                 tpsa=float(rdMolDescriptors.CalcTPSA(mol)),
                 rotatable_bonds=int(Lipinski.NumRotatableBonds(mol)),
                 formal_charge=int(Chem.GetFormalCharge(mol)),
                 rings=int(rdMolDescriptors.CalcNumRings(mol)),
                 fraction_csp3=float(rdMolDescriptors.CalcFractionCSP3(mol)),
                 chemical_source='provided_structure' if molecule is not None else 'pdbqt_smiles')
    props['lipinski_violations'] = sum((props['molecular_weight'] > 500,
                                       props['clogp'] > 5, props['hbd'] > 5, props['hba'] > 10))
    return props


def result_rows(energies, components, pdbqt_strings, scoring_function, molecule=None,
                experimental_pvalue=None, activity_type=None, temperature=298.15):
    """Build rows; reference energy is a positive stored term that is subtracted."""
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError('Temperature must be finite and positive (kelvin).')
    if experimental_pvalue is not None:
        if not math.isfinite(experimental_pvalue):
            raise ValueError('Experimental p-value must be finite.')
        if activity_type not in ('Kd', 'Ki', 'IC50', 'EC50'):
            raise ValueError('Specify experimental activity_type: Kd, Ki, IC50 or EC50.')
        if len(pdbqt_strings) != 1:
            raise ValueError('Experimental efficiencies require single-ligand docking.')
    elif activity_type is not None:
        raise ValueError('activity_type requires experimental_pvalue.')
    props = ligand_properties(pdbqt_strings, molecule)
    if len(energies) != len(components):
        raise ValueError('Pose energies and components have different lengths.')
    rows = []
    for i, (energy, detail) in enumerate(zip(energies, components)):
        if len(energy) != 5 or len(detail) != 10:
            raise ValueError('Unexpected energy/component layout; rebuild the updated Vina extension.')
        score, inter, intra, torsion, reference = map(float, energy)
        d = list(map(float, detail))
        row = dict(pose=i + 1, scoring_function=scoring_function, score=score,
                   energy_inter=inter, energy_intra=intra, torsion_term=torsion,
                   unbound_reference=reference, ligand_rigid=d[1], ligand_flex=d[2],
                   flex_rigid=d[3], flex_internal=d[4], ligand_internal=d[5],
                   rmsd_lb=d[8], rmsd_ub=d[9], delta_score=score-float(energies[0][0]),
                   intra_delta=intra-float(energies[0][2]),
                   docking_efficiency=-score/props['heavy_atom_count'],
                   inter_efficiency=-inter/props['heavy_atom_count'])
        row.update(props)
        row.update(experimental_pvalue=experimental_pvalue, activity_type=activity_type,
                   temperature_kelvin=temperature if experimental_pvalue is not None else None,
                   experimental_le=None, experimental_lle=None)
        if experimental_pvalue is not None:
            row['experimental_le'] = R_KCAL * temperature * math.log(10) * experimental_pvalue / props['heavy_atom_count']
            if props['clogp'] is not None:
                row['experimental_lle'] = experimental_pvalue - props['clogp']
        rows.append(row)
    return rows
