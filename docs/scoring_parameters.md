# Scoring parameters, physical constants and output precision

Audit date: 2026-10-08. The default scoring equations have not been changed.

## Output precision

CLI energy tables, score components and PDBQT energy remarks use six decimal
places. RMSDs and PDBQT coordinates retain three decimal places. The Python
`energies()`, `score()` and `optimize()` methods return the engine's floating-point
values without rounding. This changes presentation and API output precision,
not the computed energy, the scoring model or its predictive uncertainty.

PDBQT energy remarks remain whitespace-separated; consumers that assume fixed
energy column widths or exactly three decimals may need updating. Atom records
retain their original fixed-width layout. Serialized affinity maps retain their
existing precision: changing map serialization is a separate numerical change.

## Published model parameters

Scoring weights, radii and shape parameters are calibrated parts of a complete
model. There is no universal table of more accurate replacement weights. Mixing
parameters across force fields or inserting generic van der Waals radii changes
the fitted model. Additional digits cannot be inferred from published rounding.

The following defaults are present in `src/lib/vina.h` and `src/main/main.cpp`:

| Term | Vina | Vinardo | AD4 |
| --- | ---: | ---: | ---: |
| Gaussian 1 | -0.035579 | -0.045 | — |
| Gaussian 2 | -0.005156 | — | — |
| Repulsion | 0.840245 | 0.8 | — |
| Hydrophobic | -0.035069 | -0.035 | — |
| Hydrogen bonding | -0.587439 | -0.6 | 0.1209 |
| van der Waals | — | — | 0.1662 |
| Electrostatic | — | — | 0.1406 |
| Desolvation | — | — | 0.1322 |
| Torsional parameter | 0.05846 | 0.05846 | 0.2983 |
| Macrocycle glue weight (Vina implementation) | 50 | 50 | 50 |

These numbers multiply different terms and are not interchangeable energies.
Vina/Vinardo use a conformation-independent normalization with a torsional
parameter; AD4 uses an additive torsional term.

`src/lib/scoring_function.h` additionally fixes the Vina Gaussian widths at
0.5 and 2.0 Å, offsets at 0 and 3 Å, and its regular interaction cutoff at 8 Å.
Vinardo uses one Gaussian of width 0.8 Å and its own radii and interaction
thresholds. The implementation has a maximum cutoff of 20 Å to accommodate
macrocycle closure. AD4 desolvation uses sigma=3.6 Å and QASP=0.01097 in its
charge-dependent term. These are model parameters, not CODATA constants.

The Python metadata previously omitted the torsional weight in the default
Vinardo and AD4 tuples. They now include 0.05846 and 0.2983, respectively,
matching the C++ defaults and the six arguments accepted by `set_weights()`.
This corrects reported metadata; engine defaults were already using these values.

Vinardo is an established alternative already available through
`--scoring vinardo` or `Vina(sf_name='vinardo')`. Its publication reports improved
results on the datasets evaluated there; that is not evidence of universal
superiority on a new target. Any replacement of default weights requires separate
validation of pose recovery, affinity ranking and screening enrichment.

Sources:

- Vina: https://doi.org/10.1002/jcc.21334
- Vina 1.2: https://doi.org/10.1021/acs.jcim.1c00203
- Vinardo, equations and parameter tables:
  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155183
- AD4 parameter documentation:
  https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf
- Authoritative AD4 dielectric implementation:
  https://github.com/ccsb-scripps/AutoDock4/blob/master/distdepdiel.cc

## Physical constants for affinity-derived metrics

Use the SI defining constants, rather than the rounded factor 1.37:

| Quantity | Value | Unit |
| --- | ---: | --- |
| Boltzmann constant, kB | 1.380649e-23 | J/K |
| Avogadro constant, NA | 6.02214076e23 | 1/mol |
| Gas constant, R = NA × kB | 8.31446261815324 | J/(mol K) |
| R in kcal units, R / 4184 | 0.0019872042586408316 | kcal/(mol K) |

The thermochemical conversion is 1 kcal = 4184 J. Temperature is an explicit
input in kelvin, not a fitted universal constant.

For experimental dissociation affinity with a 1 M standard state:

    pKd = -log10(Kd / 1 M)
    delta_G_standard = -R * T * ln(10) * pKd
    LE = R * T * ln(10) * pKd / heavy_atom_count

With R in kcal units, the conversion factor is:

| T (K) | R × T × ln(10), kcal/mol per pKd unit |
| ---: | ---: |
| 298.00 | 1.363560656999 |
| 298.15 | 1.364247013034 |
| 310.00 | 1.418469139831 |

IC50 is assay-dependent and must not silently be interpreted as Kd. LLE is
pAffinity minus logP (or a separately identified logD-based variant); it has no
additional physical conversion constant. If using cLogP, report its calculation
method and the molecular protonation/tautomer policy. Dividing the docking score
by heavy-atom count is a score normalization, not an atom-by-atom energy
contribution or a measurement of experimental LE. Converting a Vina score to
pKd with RT does not establish a validated prediction of affinity.

Source: NIST/CODATA 2022 table:
https://physics.nist.gov/cuu/Constants/Table/allascii.txt

## AD4 electrostatic conversion candidate

`ad4_electrostatic::eval()` in `src/lib/potentials.h` uses 332.0 to convert
charges in elementary-charge units and distances in Å to kcal/mol.
Using NA=6.02214076e23, e=1.602176634e-19 C, and CODATA 2022 vacuum permittivity
8.8541878188e-12 F/m gives:

    NA * e**2 / (4*pi*epsilon0 * 1e-10 * 4184) = 332.063713074171

This is approximately 0.01919% larger than 332.0. It is a better-resolved physical
conversion, not demonstrated evidence of better docking. It is deliberately
not substituted in the default model: AD4 receptor interactions depend on
external AutoGrid maps, while other terms are computed locally. Updating only
one side would create inconsistent electrostatics. A future opt-in model would
need matching regenerated maps and recorded charge/dielectric conventions.

The AD4 distance-dependent dielectric uses 78.4, -8.5525, 0.003627 and 7.7839.
These define the existing Mehler–Solmajer model. The 78.4 value represents a
water dielectric reference, not vacuum permittivity and not a universal value
independent of temperature. It cannot simply be replaced with CODATA epsilon0.
