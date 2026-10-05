#!/usr/bin/env python3
"""Benchmarks claude_vrp on a set of instances against best known solutions.

Every instance is solved once by build/claude_vrp with the given time limit
and number of threads; several instances can be solved at the same time. The
results go to a CSV file and are summarized per instance type (C1, C2, R1,
R2, RC1, RC2), counting vehicles over the best known solution and the distance
gap on the instances where the vehicle count matches it.

Examples (from the repository root):
    python3 run_heuristic_benchmark.py
    python3 run_heuristic_benchmark.py --pattern 'RC2*' --seconds 60
    python3 run_heuristic_benchmark.py --solutions solutions/
    python3 run_heuristic_benchmark.py --threads 16 --jobs 1 -- --set alns.max_remove=40
"""

import argparse
import csv
import re
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
COLUMNS = ["instance", "type", "vehicles", "distance", "bks_vehicles",
           "bks_distance", "extra_vehicles", "gap_percent", "valid",
           "wall_seconds"]
# The solver's final line, e.g.
# "Best: seed 3 (alns), vehicles 18, distance 23658.12 (validated), 30.0 s".
RESULT_LINE = re.compile(
    r"^Best: .*vehicles (?P<vehicles>\d+), distance (?P<distance>[\d.]+)"
    r"(?P<validated> \(validated\))?", re.MULTILINE)


def read_best_known(path):
    """Maps instance name -> (vehicles, distance) from 'name veh dist' lines."""
    best_known = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if len(fields) == 3 and not line.startswith("#"):
            best_known[fields[0].upper()] = (int(fields[1]), float(fields[2]))
    return best_known


def instance_type(name):
    """'RC2_10_9' -> 'RC2'."""
    return name.split("_")[0].upper()


def solve(solver, instance, seconds, threads, solver_args, solutions_dir):
    """Runs the solver once; returns (vehicles, distance, valid, wall_seconds).

    With `solutions_dir`, the routes are written to <dir>/<instance>.sol.
    """
    command = [str(solver), str(instance), str(seconds), "--threads",
               str(threads), *solver_args]
    if solutions_dir is not None:
        command += ["--output", str(solutions_dir / f"{instance.stem}.sol")]
    started = time.monotonic()
    completed = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.monotonic() - started
    match = RESULT_LINE.search(completed.stdout)
    if completed.returncode != 0 or match is None:
        return None, None, False, elapsed
    return (int(match["vehicles"]), float(match["distance"]),
            match["validated"] is not None, elapsed)


def make_row(instance, result, best_known):
    vehicles, distance, valid, elapsed = result
    name = instance.stem
    bks_vehicles, bks_distance = best_known.get(name.upper(), (None, None))
    row = {"instance": name, "type": instance_type(name),
           "vehicles": vehicles,
           "distance": None if distance is None else round(distance, 2),
           "bks_vehicles": bks_vehicles, "bks_distance": bks_distance,
           "extra_vehicles": None, "gap_percent": None,
           "valid": "yes" if valid else "no",
           "wall_seconds": round(elapsed, 1)}
    if vehicles is not None and bks_vehicles is not None:
        row["extra_vehicles"] = vehicles - bks_vehicles
        row["gap_percent"] = round(100.0 * (distance / bks_distance - 1.0), 2)
    return row


def describe(row):
    if row["vehicles"] is None:
        return "no solution"
    text = f"{row['vehicles']:3d} veh {row['distance']:10.2f}"
    if row["bks_vehicles"] is not None:
        text += (f"   (BKS {row['bks_vehicles']} / {row['bks_distance']:.2f}:"
                 f" {row['extra_vehicles']:+d} veh, {row['gap_percent']:+.2f}%)")
    if row["valid"] != "yes":
        text += "   NOT VALIDATED"
    return text


def print_summary(rows):
    by_type = defaultdict(list)
    for row in rows:
        by_type[row["type"]].append(row)
    print()
    print(f"{'type':<6}{'solved':>8}{'extra veh':>11}{'at BKS veh':>12}"
          f"{'gap at BKS veh':>16}")
    all_gaps = []
    for type_name in sorted(by_type):
        group = [r for r in by_type[type_name]
                 if r["valid"] == "yes" and r["extra_vehicles"] is not None]
        at_bks = [r["gap_percent"] for r in group if r["extra_vehicles"] == 0]
        all_gaps += at_bks
        gap = f"{statistics.mean(at_bks):+.2f}%" if at_bks else "-"
        print(f"{type_name:<6}{len(group):>5}/{len(by_type[type_name]):<2}"
              f"{sum(r['extra_vehicles'] for r in group):>+11d}"
              f"{len(at_bks):>9}/{len(group):<2}{gap:>16}")
    solved = [r for r in rows if r["valid"] == "yes"]
    extra = sum(r["extra_vehicles"] for r in solved
                if r["extra_vehicles"] is not None)
    print(f"\n{len(solved)} of {len(rows)} instances solved and validated; "
          f"{extra:+d} vehicles over the best known solutions in total"
          + (f"; mean gap {statistics.mean(all_gaps):+.2f}% on the "
             f"{len(all_gaps)} instances at the best known vehicle count"
             if all_gaps else ""))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog="Arguments after '--' are passed to the solver.")
    parser.add_argument("--instances", default=ROOT / "data", type=Path,
                        help="directory with instance files (default: data)")
    parser.add_argument("--pattern", default="*",
                        help="glob for instance names, e.g. 'C1*' (default: all)")
    parser.add_argument("--seconds", type=float, default=120,
                        help="time limit per instance (default: 120)")
    parser.add_argument("--threads", type=int, default=8,
                        help="solver threads per instance (default: 8)")
    parser.add_argument("--jobs", type=int, default=2,
                        help="instances solved at the same time (default: 2); "
                             "keep threads x jobs within the CPU cores")
    parser.add_argument("--bks", default=ROOT / "experiments" / "bks.txt",
                        type=Path, help="best known solutions file")
    parser.add_argument("--solver", default=ROOT / "build" / "claude_vrp",
                        type=Path, help="solver executable")
    parser.add_argument("--csv", type=Path, default=None,
                        help="output file (default: benchmark_results/"
                             "run_<date>_<time>.csv)")
    parser.add_argument("--solutions", type=Path, default=None,
                        help="directory for the routes of each solution "
                             "(<instance>.sol, SINTEF format); default: none")
    parser.add_argument("solver_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    solver_args = args.solver_args
    if solver_args[:1] == ["--"]:
        solver_args = solver_args[1:]

    if not args.solver.is_file():
        sys.exit(f"{args.solver} not found; build it first "
                 "(cmake -S . -B build && cmake --build build)")
    instances = sorted(p for p in args.instances.glob(args.pattern)
                       if p.is_file() and p.suffix.lower() == ".txt")
    if not instances:
        sys.exit(f"no instances match {args.pattern!r} in {args.instances}")
    best_known = read_best_known(args.bks)
    if args.solutions is not None:
        args.solutions.mkdir(parents=True, exist_ok=True)
    csv_path = args.csv or (ROOT / "benchmark_results" /
                            f"run_{datetime.now():%Y%m%d_%H%M%S}.csv")

    print(f"{len(instances)} instances, {args.seconds:g} s and "
          f"{args.threads} threads each, {args.jobs} at a time"
          + (f", solver arguments: {' '.join(solver_args)}" if solver_args
             else ""))
    rows = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(solve, args.solver, instance, args.seconds,
                               args.threads, solver_args,
                               args.solutions): instance
                   for instance in instances}
        for done, future in enumerate(as_completed(futures), start=1):
            row = make_row(futures[future], future.result(), best_known)
            rows.append(row)
            print(f"[{done:>{len(str(len(instances)))}}/{len(instances)}] "
                  f"{row['instance']:<10} {describe(row)}", flush=True)

    rows.sort(key=lambda r: r["instance"])
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print_summary(rows)
    print(f"results: {csv_path}")


if __name__ == "__main__":
    main()
