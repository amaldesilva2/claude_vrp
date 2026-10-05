#!/usr/bin/env python3
"""Compares benchmark configs from experiments/results.csv, instance by instance.

For each instance the configs are ranked lexicographically (mean vehicles,
then mean distance); the summary counts wins and gives each config's mean
rank, total vehicles over the BKS, and mean distance gap on the instances
where it matches the best vehicle count.

usage: compare.py CONFIG [CONFIG ...]
"""
import csv
import os
import statistics
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    configs = sys.argv[1:]
    rows = defaultdict(list)
    bks = {}
    with open(os.path.join(HERE, "results.csv")) as f:
        for r in csv.DictReader(f):
            if r["config"] in configs and r["vehicles"]:
                rows[(r["config"], r["instance"])].append(
                    (int(r["vehicles"]), float(r["distance"])))
                bks[r["instance"]] = (int(r["bks_vehicles"]),
                                      float(r["bks_distance"]))
    instances = sorted({i for (_, i) in rows})
    ranks = defaultdict(list)
    wins = defaultdict(int)
    over = defaultdict(float)
    gaps = defaultdict(list)
    width = max(len(c) for c in configs)
    print(f"{'instance':9} " + " ".join(f"{c:>{max(width, 18)}}"
                                       for c in configs))
    for inst in instances:
        stats = {}
        for c in configs:
            v = rows.get((c, inst))
            if v:
                stats[c] = (statistics.mean(x[0] for x in v),
                            statistics.mean(x[1] for x in v))
        order = sorted(stats, key=lambda c: stats[c])
        for rank, c in enumerate(order):
            ranks[c].append(rank + 1)
        wins[order[0]] += 1
        best_veh = stats[order[0]][0]
        cells = []
        for c in configs:
            if c not in stats:
                cells.append(f"{'-':>{max(width, 18)}}")
                continue
            veh, dist = stats[c]
            over[c] += veh - bks[inst][0]
            if veh == best_veh:
                gaps[c].append(100 * (dist / bks[inst][1] - 1))
            mark = "*" if c == order[0] else " "
            cells.append(f"{veh:6.2f}/{dist:9.0f}{mark}".rjust(max(width, 18)))
        print(f"{inst:9} " + " ".join(cells))
    print()
    for c in configs:
        print(f"{c:{width}}  wins {wins[c]:2d}  mean rank "
              f"{statistics.mean(ranks[c]):.2f}  vehicles over BKS "
              f"{over[c]:+6.2f}  gap where best-vehicles "
              f"{statistics.mean(gaps[c]) if gaps[c] else float('nan'):+.2f}%"
              f" ({len(gaps[c])} inst)")


if __name__ == "__main__":
    main()
