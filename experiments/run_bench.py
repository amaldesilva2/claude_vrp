#!/usr/bin/env python3
"""Benchmarks the solver on a fixed instance set and reports gaps to the BKS.

Each (instance, seed) pair is one single-threaded run of build/claude_vrp; runs are
executed in parallel. Results are appended to experiments/results.csv.

usage: run_bench.py NAME [--seconds S] [--seeds K] [--jobs J] [--instances SET]
                    [-- vrp args...]
example: run_bench.py ils-default -- --search ils
"""
import argparse
import csv
import os
import re
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
VRP = os.path.join(ROOT, "build", "claude_vrp")

INSTANCE_SETS = {
    # Two instances per type, plus RC2_10_9 (the original target instance).
    "core": ["C1_10_2", "C1_10_6", "C2_10_3", "C2_10_6", "R1_10_1", "R1_10_5",
             "R2_10_1", "R2_10_5", "RC1_10_1", "RC1_10_5", "RC2_10_1",
             "RC2_10_5", "RC2_10_9"],
    "quick": ["C1_10_2", "R1_10_1", "R2_10_1", "RC1_10_1", "RC2_10_1",
              "RC2_10_9"],
    # Core instances where the vehicle count is above the BKS.
    "vehicles": ["C1_10_2", "C1_10_6", "C2_10_3", "C2_10_6", "RC2_10_1"],
    "vehicles+": ["C1_10_2", "C1_10_6", "C2_10_3", "C2_10_6", "RC2_10_1",
                  "R1_10_1"],
    "all": None,  # filled from bks.txt
}


def load_bks():
    bks = {}
    with open(os.path.join(HERE, "bks.txt")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            name, vehicles, distance = line.split()
            bks[name] = (int(vehicles), float(distance))
    return bks


def run_one(instance, seed, seconds, extra, threads=1):
    # Seeds of multi-threaded runs are spaced so that runs do not share seeds.
    cmd = [VRP, os.path.join(DATA, instance + ".TXT"), str(seconds),
           "--threads", str(threads), "--seed", str(1 + (seed - 1) * threads)]
    cmd += extra
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    m = re.search(r"^Best: .*vehicles (\d+), distance ([\d.]+)", out, re.M)
    if not m:
        return instance, seed, None, None
    return instance, seed, int(m.group(1)), float(m.group(2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("name")
    parser.add_argument("--seconds", type=float, default=30)
    parser.add_argument("--seeds", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--instances", default="core")
    parser.add_argument("--run-threads", type=int, default=1,
                        help="threads per run (jobs is divided by this)")
    parser.add_argument("extra", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    extra = args.extra[1:] if args.extra[:1] == ["--"] else args.extra

    bks = load_bks()
    instances = INSTANCE_SETS[args.instances] or sorted(bks)
    tasks = [(i, s) for i in instances for s in range(1, args.seeds + 1)]
    jobs = max(1, args.jobs // args.run_threads)
    with ThreadPoolExecutor(jobs) as pool:
        results = list(pool.map(
            lambda t: run_one(t[0], t[1], args.seconds, extra,
                              args.run_threads), tasks))

    csv_path = os.path.join(HERE, "results.csv")
    new_file = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="") as f:
        writer = csv.writer(f)
        if new_file:
            writer.writerow(["config", "args", "seconds", "instance", "seed",
                             "vehicles", "distance", "bks_vehicles",
                             "bks_distance"])
        for inst, seed, v, d in results:
            writer.writerow([args.name, " ".join(extra), args.seconds, inst,
                             seed, v, d, bks[inst][0], bks[inst][1]])

    print(f"config {args.name}: {' '.join(extra)}  ({args.seconds:g} s, "
          f"{args.seeds} seeds)")
    gaps, extra_vehicles, failures = [], 0.0, 0
    for inst in instances:
        rows = [r for r in results if r[0] == inst]
        ok = [r for r in rows if r[2] is not None]
        failures += len(rows) - len(ok)
        if not ok:
            print(f"  {inst:9} FAILED")
            continue
        bv, bd = bks[inst]
        veh = statistics.mean(r[2] for r in ok)
        dist = statistics.mean(r[3] for r in ok)
        gap = 100.0 * (dist / bd - 1.0)
        extra_vehicles += veh - bv
        gaps.append(gap)
        print(f"  {inst:9} veh {veh:6.2f} (bks {bv:3d})  dist {dist:10.2f}"
              f"  gap {gap:+6.2f}%   seeds: "
              + " ".join(f"{r[2]}/{r[3]:.0f}" for r in ok))
    print(f"  SUMMARY {args.name}: vehicles over BKS {extra_vehicles:+.2f}, "
          f"mean distance gap {statistics.mean(gaps):+.2f}%"
          + (f", {failures} FAILED runs" if failures else ""))


if __name__ == "__main__":
    sys.exit(main())
