# Precalculation performance and regression checks

The candidate transfers completed interaction tables with C++11 move assignment
instead of copying their nested vectors. Both sources are local temporaries and
are never read after the transfer. Scoring, interpolation, search and RNG are
unchanged. This removes two deep copies from each ligand-loading path and one
from atom-type map preparation; the dense quadratic table remains in place.

## Numerical regression

Build the reference at commit `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645` and the candidate
with the same compiler, dependencies and flags. From the repository root:

```bash
python tests/compare_binaries.py /path/to/reference/vina build/linux/release/vina
```

The script uses bundled rigid, flexible, macrocycle and simultaneous-ligand
examples with Vina, Vinardo and AD4. It requires identical printed energy
components for scoring and byte-identical pose files for docking. Docking uses
CPU=1, seed=42, exhaustiveness=1, max_evals=1000 and at most three poses. These
are regression checks, not evidence of docking accuracy or convergence. The two
prepared simultaneous ligands are centered at the origin; score-only contains
their supplied coordinates, while docking uses the supplied receptor box center.
AD4 uses the bundled maps, whose grid dimensions supersede the CLI dimensions.

## Loading benchmark

After building each version, link the harness against its library objects:

```bash
g++ -O3 -DNDEBUG -std=c++11 -pthread -Isrc/lib \
    tests/benchmark_ligand_loading.cpp \
    $(find build/linux/release -name '*.o' ! -name main.o ! -name split.o) \
    -lboost_system -lboost_thread -lboost_serialization \
    -lboost_filesystem -lboost_program_options \
    -o build/linux/release/benchmark_loading

python tests/measure_loading.py /path/to/reference/benchmark_loading \
    build/linux/release/benchmark_loading \
    example/basic_docking/solution/1iep_ligand.pdbqt \
    --output benchmark_results.json
```

Linux only: peak resident memory is collected with `wait4`. Each binary performs
one warm-up load, then three timed replacements. Five runs alternate reference
and candidate order. The reported time is the median seconds per timed load;
peak RSS covers the entire harness process, including its warm-up. The vector
path contains one ligand to isolate overload behavior. Numerical regression
separately covers actual simultaneous ligands.

Recorded results in `benchmark_results.json` were measured in this workspace
with GCC 13.3.0, Boost 1.83, `-O3 -DNDEBUG -std=c++11`, dynamic linking and the
bundled 1IEP ligand. Shared-machine timings are preliminary and apply only to
ligand loading, not total docking or scientific accuracy. Raw repetitions and
per-run memory are included so variability can be inspected.

| Scoring | Loading API | Reference ms | Candidate ms | Time reduction |
| --- | --- | ---: | ---: | ---: |
| Vina | string | 333.7 | 260.9 | 21.8% |
| Vina | vector | 339.0 | 280.6 | 17.2% |
| Vinardo | string | 394.5 | 281.2 | 28.7% |
| Vinardo | vector | 321.2 | 234.6 | 27.0% |
| AD4 | string | 520.6 | 444.3 | 14.7% |
| AD4 | vector | 537.6 | 453.8 | 15.6% |

Peak RSS was approximately 724 to 487 MiB for Vina/Vinardo, and 759 to 509 MiB
for AD4. All 24 numerical regression comparisons passed locally.

## Report validation without compiling Vina

Run `python tests/test_result_properties.py` to check named energy components,
reference-energy signs, heavy-atom counting, derived metrics, optional RDKit
properties, and CSV/TSV export. The tests use a stub native extension and do not
build or run the docking engine.
