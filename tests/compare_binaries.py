"""Check scores and seeded poses against a reference executable, using bundled examples."""
import argparse
import re
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("reference", type=Path)
parser.add_argument("candidate", type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
executables = [args.reference.resolve(), args.candidate.resolve()]

basic = root / "example/basic_docking/solution"
flex = root / "example/flexible_docking/solution"
macro = root / "example/docking_with_macrocycles/solution"
multi = root / "example/mulitple_ligands_docking/solution"
cases = [
    ("rigid", basic / "1iep_receptor.pdbqt", [basic / "1iep_ligand.pdbqt"], None,
     (15.190, 53.903, 16.917)),
    ("flex", flex / "1fpu_receptor_rigid.pdbqt", [flex / "1iep_ligand.pdbqt"],
     flex / "1fpu_receptor_flex.pdbqt", (15.190, 53.903, 16.917)),
    ("macrocycle", macro / "BACE_1_receptor.pdbqt", [macro / "BACE_1_ligand.pdbqt"], None,
     (0, 0, 0)),
    ("multiple", multi / "5x72_receptor.pdbqt",
     [multi / "5x72_ligand_p59.pdbqt", multi / "5x72_ligand_p69.pdbqt"], None, (0, 0, 0)),
]
# Read the supplied grid centers for the larger examples.
for idx, box in [(2, macro / "BACE_1_receptor_vina_box.txt"),
                 (3, multi / "5x72_receptor.box.txt")]:
    centers = dict(re.findall(r"center_([xyz])\s*=\s*([-\d.]+)", box.read_text()))
    case = cases[idx]
    cases[idx] = (*case[:4], tuple(float(centers[axis]) for axis in "xyz"))

with tempfile.TemporaryDirectory() as tmp:
    for name, receptor, ligands, flexible, center in cases:
        for scoring in ("vina", "vinardo", "ad4"):
            for mode in ("score_only", "docking"):
                results = []
                for idx, executable in enumerate(executables):
                    output = Path(tmp) / f"{name}-{scoring}-{mode}-{idx}.pdbqt"
                    receptor_input = ["--maps", str(receptor.with_suffix(""))] if scoring == "ad4" else ["--receptor", str(receptor)]
                    command = [str(executable), *receptor_input, "--ligand",
                               *map(str, ligands), "--scoring", scoring, "--cpu", "1",
                               "--seed", "42", "--verbosity", "1"]
                    if flexible:
                        command += ["--flex", str(flexible)]
                    # These two prepared ligands are centered at the origin, not
                    # at their receptor's docking site; scoring must contain them.
                    box_center = (0, 0, 0) if name == "multiple" and mode == "score_only" else center
                    for axis, value in zip("xyz", box_center):
                        command += [f"--center_{axis}", str(value), f"--size_{axis}", "30"]
                    if mode == "score_only":
                        command += ["--score_only"]
                    else:
                        command += ["--exhaustiveness", "1", "--max_evals", "1000",
                                    "--num_modes", "3", "--out", str(output)]
                    run = subprocess.run(command, text=True, capture_output=True, timeout=120)
                    if run.returncode:
                        raise RuntimeError(f"{name}/{scoring}/{mode}: {run.stderr}\n{run.stdout}")
                    if mode == "score_only":
                        result = [line for line in run.stdout.splitlines()
                                  if "kcal/mol" in line]
                        if not result:
                            raise RuntimeError("No scores found in output")
                    else:
                        result = output.read_text()
                        if "ATOM" not in result and "HETATM" not in result:
                            raise RuntimeError("No pose atoms found in output")
                    results.append(result)
                if results[0] != results[1]:
                    raise AssertionError(f"Output mismatch: {name}/{scoring}/{mode}")
                print(f"PASS {name}/{scoring}/{mode}", flush=True)
