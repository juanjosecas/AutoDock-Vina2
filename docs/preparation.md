# Automatic preparation with Meeko

Preparation is an optional Python layer requiring Python 3.10 or newer. The C++ `vina` executable still takes PDBQT; its scoring/search code is unchanged. Install the dependencies into the same environment as this fork's Python package:

```bash
python -m pip install 'meeko==0.7.1' rdkit gemmi scipy
```

The package also declares a `preparation` extra (`vina[preparation]`). Install the fork itself using your usual build procedure; installing upstream `vina` does not supply these new methods. Meeko is imported only when conversion is needed, and the supported version is 0.7.1.

## Python: prepare and load automatically

```python
from vina import Vina

v = Vina(verbosity=2)
v.set_receptor('protein.pdb', prepared_filename='protein.prepared.pdbqt')
v.set_ligand_from_file('ligand.sdf')

# Or start from a SMILES, with deterministic 3D generation:
v.set_ligand_from_smiles('CCO', seed=42)

reports = v.preparation_reports()
```

`set_ligand_from_file` accepts PDBQT unchanged, or prepares SDF/MOL/MOL2 in memory before loading it. One molecule is required per file. A list of files continues to mean simultaneous ligands, not independent screening. Multi-record SDF files and disconnected salts/fragments are rejected rather than silently selecting or combining compounds. For existing PDBQT, no Meeko dependency is needed.

`set_receptor` converts a rigid PDB/PQR input into `<input_stem>.prepared.pdbqt` by default. Supply `prepared_filename` for a different location. Existing files are protected; use `overwrite=True` to explicitly replace generated files. Flexible-residue files must already be PDBQT. mmCIF conversion and automatic flexible-residue generation are not exposed by this first integration.

Explicit receptor preparation choices use Meeko syntax:

```python
v.set_receptor(
    'protein.pdb',
    prepared_filename='protein.prepared.pdbqt',
    delete_residues='A:501,B:601',  # only these residues are explicitly requested for deletion
    set_template='A:42=HID',       # choose the appropriate histidine state for your system
    preparation_verbosity=2,
)
```

Preparation inherits the Vina instance's verbosity. `preparation_verbosity=0/1/2` overrides it. Level 1 prints the main steps and streams Meeko's receptor diagnostics; level 2 adds atom types, per-residue templates, atom/hydrogen counts and formal charges. Errors include the Meeko diagnostic even in quiet mode. Unassigned ligand stereocenters produce a warning.

## Standalone helpers and CLI

```python
from vina.preparation import prepare_ligand, prepare_receptor

ligand = prepare_ligand('ligand.sdf', 'ligand.pdbqt', verbosity=2)
receptor = prepare_receptor('protein.pdb', 'protein.pdbqt', verbosity=2)
# ligand['pdbqt_string']; ligand['report']
# receptor['pdbqt_filename']; receptor['report']
```

A ligand helper also accepts an RDKit Mol or `input_format='smiles'`.

```bash
python -m vina.preparation ligand ligand.sdf -o ligand.pdbqt --verbosity 2
python -m vina.preparation ligand 'CCO' --smiles -o ethanol.pdbqt --seed 42
python -m vina.preparation receptor protein.pdb -o protein.pdbqt --verbosity 2
```

Whenever an output PDBQT is written, a `<output_stem>.preparation.json` audit accompanies it. Reports include software versions, input file SHA-256 where applicable, preparation settings, operations and output atom types/counts, partial-charge sums and ligand TORSDOF. Receptor reports retain the Meeko log, selected templates and per-residue input/prepared atom counts and formal charges. The stored command may mention a temporary working directory which has already been removed. Reports from automatic API conversions can also be saved explicitly:

```python
import json
with open('preparation_reports.json', 'w') as stream:
    json.dump(v.preparation_reports(), stream, indent=2)
```

## Chemical behavior

Ligands retain their input formal charge, protonation/tautomer state and specified stereochemistry. RDKit adds explicit hydrogens. Existing 3D heavy-atom coordinates are retained; no minimization is performed. For 2D/SMILES input, RDKit generates one ETKDGv3 conformer with the requested seed and minimizes it with MMFF94 (500 iterations maximum). Missing force-field parameters, embedding failures or nonconvergence stop preparation. Meeko applies Gasteiger charges, AutoDock atom typing, its default torsion handling and merging of nonpolar hydrogens. PDBQT coordinates and charges have their format's normal rounding; their charge sum can differ slightly from the integer formal charge.

Receptors are matched to Meeko residue templates, which define chemistry/protonation and add missing hydrogens. This does **not** run PROPKA or select states at a requested pH. Input PDB formal charges may be unspecified; the report's input formal charges reflect RDKit's interpretation. PQR input is also parameterized using Gasteiger, not its original PQR partial charges. Waters, metals, cofactors and incomplete residues are not silently deleted. Unmatched residues or interrupted residue records stop conversion with the original Meeko diagnostics. Resolve them or request specific deletions/template overrides explicitly. Preparation helpers should run serially: ligand library-output capture temporarily redirects Python stdout/stderr.

## Validation

`tests/test_preparation.py` uses real RDKit/Meeko preparation to check reproducibility, preservation of charge/stereochemistry and 3D coordinates, output protection, invalid/multirecord inputs, receptor templates/hydrogen addition and rejection of incomplete receptors. Native Vina loading is stubbed to verify API routing without compiling or executing docking.

References: [Meeko repository](https://github.com/forlilab/Meeko), [ligand preparation](https://meeko.readthedocs.io/en/develop/lig_prep_basic.html), [receptor templates](https://meeko.readthedocs.io/en/develop/rec_overview.html).
