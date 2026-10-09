# Docking results, energy components and ligand properties

The docking score and pose ranking are unchanged. The terminal table now reports
score, intermolecular and internal energies, effective torsional contribution,
unbound reference, score difference from the best pose, real heavy-atom count,
`-score/HAC`, and the existing RMSD bounds. Energies have six decimal places.
This changes the console table format; parsers should use the Python report API
instead of parsing terminal output. PDBQT coordinate and RMSD precision is unchanged.

## Python reports

After `v.dock(...)`:

```python
rows = v.results(n_poses=20, energy_range=5.0)
v.write_results('results.tsv', n_poses=20, energy_range=5.0)
```

`results()` returns a list of dictionaries. CSV/TSV export needs no pandas.
Missing properties are `None` in Python and empty cells in exports. Existing files
are protected unless `overwrite=True`. The filters are identical to `energies()`.
Numerical calculations use the full engine precision; exported floats use six decimals.

The updated native extension is required for `get_poses_components()`; an older
installed Vina binary or wheel will not gain these functions from a source update.
This change was checked with Python tests and source inspection without building
or running the native docking engine.

### Energy columns

| Column | Meaning |
| --- | --- |
| `score` | Final docking score |
| `energy_inter` | Ligand interaction with rigid receptor and flexible residues |
| `energy_intra` | Internal energy of the mobile system |
| `torsion_term` | Effective contribution from the conformation-independent term |
| `unbound_reference` | Stored reference energy, **subtracted** from the total |
| `ligand_rigid` | Ligand–rigid receptor contribution |
| `ligand_flex` | Ligand–flexible residue contribution |
| `flex_rigid` | Flexible residue–rigid receptor contribution |
| `flex_internal` | Within/between flexible residue contributions |
| `ligand_internal` | Internal ligand contribution |
| `rmsd_lb`, `rmsd_ub` | Existing lower/upper RMSD bounds, in angstrom |

Energy units are kcal/mol. The identity is:

`score = energy_inter + energy_intra + torsion_term - unbound_reference`

In Vina/Vinardo, the reference comes from the best pose's internal system energy.
In AD4, it equals the current pose's internal system energy. Earlier Python
API documentation called the last AD4 column `-intra`, but the stored value is
positive `intra` and subtraction occurs in the total. The docstrings now reflect this.

The Vina/Vinardo torsional contribution is not a constant times RDKit's rotatable
bond count. Internal energy and `intra_delta` are not ligand strain energy relative
to a separately optimized free ligand. These components do not separate Gaussian,
repulsion, hydrophobic and hydrogen-bond potentials.

### Derived metrics

| Column | Definition |
| --- | --- |
| `heavy_atom_count` | Real non-hydrogen ligand atoms; G0–G3 and W excluded, CG0–CG3 included |
| `docking_efficiency` | `-score / heavy_atom_count` |
| `inter_efficiency` | `-energy_inter / heavy_atom_count` |
| `delta_score` | Current score minus best score |
| `intra_delta` | Current internal energy minus best pose internal energy |
| `docking_torsdof` | Prepared PDBQT TORSDOF, if available |
| `docking_branches` | PDBQT BRANCH count |

Efficiencies are score normalizations in kcal/mol/heavy atom, not individual atomic
energy assignments. Score gaps are not confidence estimates or pose probabilities.
For simultaneous docking, energies and HAC refer to the complete ligand set;
per-ligand chemical descriptors and experimental efficiencies are not assigned.

## Optional chemical descriptors

RDKit is optional (`pip install 'vina[properties]'`). With RDKit installed,
`REMARK SMILES` from Meeko-prepared PDBQT inputs is used automatically. Otherwise
provide the original SMILES or RDKit Mol:

```python
from rdkit import Chem

mol = Chem.SDMolSupplier('ligand.sdf', removeHs=False)[0]
rows = v.results(molecule=mol)
v.write_results('results.csv', molecule=mol)
# A SMILES string is also accepted: molecule='CCO' for an ethanol ligand.
```

The structure must match the loaded ligand, including its intended protonation
state. Heavy-atom elemental composition is checked; this check cannot establish
connectivity identity. Bond orders are never guessed from PDBQT coordinates.
Explicitly supplying a structure without RDKit raises an informative ImportError.

Available descriptors: molecular weight, Wildman–Crippen cLogP, HBD, HBA, TPSA,
RDKit rotatable bonds, formal charge, ring count and fraction Csp3. Lipinski
violations count MW > 500, cLogP > 5, HBD > 5 and HBA > 10. These flags are not
predictions of absorption or clinical suitability. `chemical_source` records whether
properties came from a provided structure, a PDBQT SMILES, or are unavailable.
Descriptors are calculated once per report and repeated for each pose.

## Experimental efficiencies

```python
v.write_results('results.tsv', molecule=mol,
                experimental_pvalue=7.0, activity_type='Kd', temperature=298.15)
```

The p-value is `-log10(activity in mol/L)`. Supported labels are Kd, Ki, IC50 and
EC50; the label is mandatory and retained. IC50/EC50-based values are activity
normalizations, not thermodynamic affinity measurements.

`experimental_le = R * T * ln(10) * experimental_pvalue / HAC`

`experimental_lle = experimental_pvalue - cLogP`

R = 8.31446261815324/4184 kcal/mol/K. LLE is unavailable without chemical descriptors.
No experimental affinity, pKd or LLE is inferred from a docking score. All poses
of a ligand repeat the same experimental efficiency values.

## Validation without native compilation

```bash
python tests/test_result_properties.py
```

The tests exercise the reporting module and actual Python wrapper with a stub
native extension. RDKit-dependent checks skip when the optional package is absent.
Native compilation and docking regression were not performed for this change.
