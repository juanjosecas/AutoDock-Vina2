"""Linux benchmark: paired runs, median loading time and per-process peak RSS."""
import argparse
import json
import os
import statistics
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("reference", type=Path)
parser.add_argument("candidate", type=Path)
parser.add_argument("ligand", type=Path)
parser.add_argument("--repeats", type=int, default=5)
parser.add_argument("--iterations", type=int, default=3)
parser.add_argument("--output", type=Path, default=Path("benchmark_results.json"))
args = parser.parse_args()
if args.repeats < 1 or args.iterations < 1:
    parser.error("repeats and iterations must be positive")
binaries = {"baseline": args.reference.resolve(), "candidate": args.candidate.resolve()}
rows = []
for scoring in ("vina", "vinardo", "ad4"):
    for mode in ("single", "vector"):
        times = {name: [] for name in binaries}
        memory = {name: [] for name in binaries}
        for repeat in range(args.repeats):
            order = list(binaries) if repeat % 2 == 0 else list(reversed(binaries))
            for name in order:
                with tempfile.TemporaryFile(mode="w+") as out:
                    proc = subprocess.Popen([
                        str(binaries[name]), str(args.ligand.resolve()),
                        str(args.iterations), scoring, mode,
                    ], stdout=out)
                    _, status, usage = os.wait4(proc.pid, 0)
                    proc.returncode = os.waitstatus_to_exitcode(status)
                    if proc.returncode:
                        raise RuntimeError(f"{name} exited with {proc.returncode}")
                    out.seek(0)
                    times[name].append(float(out.read().strip()))
                    memory[name].append(usage.ru_maxrss)
        row = {
            "scoring": scoring, "mode": mode,
            "baseline_seconds": statistics.median(times["baseline"]),
            "candidate_seconds": statistics.median(times["candidate"]),
            "samples": times, "peak_rss_kib": memory,
        }
        row["reduction_percent"] = 100 * (1 - row["candidate_seconds"] / row["baseline_seconds"])
        rows.append(row)
        print(f"{scoring}/{mode}: {row['reduction_percent']:.1f}% shorter loading time", flush=True)
args.output.write_text(json.dumps(rows, indent=2) + "\n")
