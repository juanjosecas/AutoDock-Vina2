# Console output

Use `--verbosity 0` for quiet output, `--verbosity 1` (default) for progress and ranked results, or `--verbosity 2` for search settings and detailed pose energies. Warnings and errors remain visible in quiet mode. The Python `Vina(verbosity=...)` option uses the same levels.

The docking search prints the actual random seed, number of trajectories and active workers. An interactive terminal shows a single updating bar, step-budget percentage, elapsed seconds and approximate remaining search time. Rendering is limited to five updates per second; progress callbacks retain their existing frequency and serialization. The estimate excludes subsequent refinement and rescoring, and is disabled when an evaluation limit can stop trajectories early. Search completion may therefore show less than 100% of the step budget.

Redirected output and CI logs use plain lines at approximately 10% intervals, plus a final completion line. They contain no carriage-return animations or ANSI colors. `TERM=dumb` also selects plain output.

Compatible terminals use cyan headings and green completion messages and best-pose scores. Set `NO_COLOR=1` to disable colors. Windows consoles enable virtual-terminal processing only when supported; otherwise output stays uncolored.

Results appear in two aligned tables: ranked poses (score, difference from the best score, negative score per ligand heavy atom, and RMSD bounds), followed by intermolecular, internal, torsional and reference energies. Energies retain six decimal places and RMSD retains three. Formatting uses local streams so it does not change precision settings of subsequent output. Docking scores per atom are descriptive metrics, not experimental ligand efficiencies.

The changes affect console presentation only: scoring equations, random seeds, search trajectories, pose ordering and machine-readable outputs are unchanged. No build or docking run was performed for this change; validation was limited to source review and whitespace checks.
