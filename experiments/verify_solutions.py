#!/usr/bin/env python3
"""Checks VRPTW solution files independently of the C++ solver.

Reads solution files in the SINTEF format written by `claude_vrp --output`
("Instance name : X", then "Route i : c1 c2 ..." lines), loads the matching
instance from the data directory, and checks that

  * every customer is visited exactly once,
  * every route respects the vehicle capacity,
  * service starts within each customer's time window (arriving early means
    waiting; travel time equals the Euclidean distance),
  * every vehicle is back at the depot before it closes,
  * no more routes are used than the instance has vehicles.

It prints the number of routes and the recomputed distance (plain double
Euclidean). With --results (a CSV from run_heuristic_benchmark.py) the
recomputed values are also compared with the ones the solver reported.
The exit code is 1 if any check fails.

Examples (from the repository root):
    python3 experiments/verify_solutions.py solutions/
    python3 experiments/verify_solutions.py solutions/RC2_10_9.sol --data data
    python3 experiments/verify_solutions.py solutions/ --results benchmark_results/run_x.csv
"""

import argparse
import csv
import math
import sys
from pathlib import Path

TOLERANCE = 1e-6  # for floating-point time comparisons


def read_instance(path):
    """Returns (vehicles, capacity, nodes); nodes[i] = (x, y, demand, ready,
    due, service), node 0 is the depot."""
    lines = Path(path).read_text().splitlines()
    vehicles = capacity = None
    nodes = {}
    for i, line in enumerate(lines):
        fields = line.split()
        if fields[:2] == ["NUMBER", "CAPACITY"]:
            vehicles, capacity = (int(x) for x in lines[i + 1].split())
        elif len(fields) == 7 and fields[0].isdigit():
            nodes[int(fields[0])] = tuple(float(x) for x in fields[1:])
    if vehicles is None or 0 not in nodes:
        raise ValueError(f"{path}: not a Solomon-format instance")
    return vehicles, capacity, nodes


def read_solution(path):
    """Returns (instance name, list of routes as customer lists)."""
    name = None
    routes = []
    for line in Path(path).read_text().splitlines():
        key, _, value = line.partition(":")
        if key.strip() == "Instance name":
            name = value.strip()
        elif key.strip().startswith("Route"):
            routes.append([int(c) for c in value.split()])
    return name, routes


def check(instance_path, routes):
    """Returns (problems, distance) for the routes on the instance."""
    vehicles, capacity, nodes = read_instance(instance_path)

    def travel(a, b):
        return math.hypot(nodes[a][0] - nodes[b][0], nodes[a][1] - nodes[b][1])

    problems = []
    visits = {}
    total = 0.0
    for number, route in enumerate(routes, start=1):
        time = nodes[0][3]  # depot ready time
        load = 0
        previous = 0
        for customer in route:
            if customer not in nodes or customer == 0:
                problems.append(f"route {number}: unknown customer {customer}")
                continue
            visits[customer] = visits.get(customer, 0) + 1
            load += nodes[customer][2]
            total += travel(previous, customer)
            arrival = time + nodes[previous][5] + travel(previous, customer)
            time = max(arrival, nodes[customer][3])
            if time > nodes[customer][4] + TOLERANCE:
                problems.append(f"route {number}: customer {customer} starts at "
                                f"{time:.2f}, due {nodes[customer][4]:g}")
            previous = customer
        total += travel(previous, 0)
        back = time + nodes[previous][5] + travel(previous, 0)
        if back > nodes[0][4] + TOLERANCE:
            problems.append(f"route {number}: back at the depot at {back:.2f}, "
                            f"closes at {nodes[0][4]:g}")
        if load > capacity:
            problems.append(f"route {number}: load {load:g} > capacity "
                            f"{capacity}")

    for customer in nodes:
        if customer != 0 and visits.get(customer, 0) != 1:
            problems.append(f"customer {customer} visited "
                            f"{visits.get(customer, 0)} times")
    if len(routes) > vehicles:
        problems.append(f"{len(routes)} routes, only {vehicles} vehicles")
    return problems, total


def find_instance(data_dir, name):
    for path in Path(data_dir).iterdir():
        if path.stem.upper() == name.upper() and path.suffix.lower() == ".txt":
            return path
    raise FileNotFoundError(f"instance {name} not found in {data_dir}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("solutions", nargs="+", type=Path,
                        help="solution files or directories of .sol files")
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).resolve().parent.parent / "data",
                        help="directory with the instances (default: data)")
    parser.add_argument("--results", type=Path, default=None,
                        help="CSV from run_heuristic_benchmark.py to compare "
                             "the reported vehicles and distances with")
    args = parser.parse_args()

    files = []
    for path in args.solutions:
        files += sorted(path.glob("*.sol")) if path.is_dir() else [path]
    reported = {}
    if args.results is not None:
        for row in csv.DictReader(open(args.results)):
            if row["vehicles"]:
                reported[row["instance"].upper()] = (int(row["vehicles"]),
                                                    float(row["distance"]))

    failed = 0
    for path in files:
        name, routes = read_solution(path)
        try:
            if name is None:
                raise ValueError("no 'Instance name' line")
            problems, distance = check(find_instance(args.data, name), routes)
        except (OSError, ValueError) as error:
            problems, distance = [str(error)], math.nan
        if name is not None and name.upper() in reported:
            vehicles, reported_distance = reported[name.upper()]
            if vehicles != len(routes):
                problems.append(f"{len(routes)} routes, reported {vehicles}")
            if abs(round(distance, 2) - reported_distance) > 0.0051:
                problems.append(f"distance {distance:.2f}, reported "
                                f"{reported_distance:.2f}")
        status = "OK  " if not problems else "FAIL"
        print(f"{status} {path.name:<16} {len(routes):4d} routes  "
              f"distance {distance:12.4f}")
        for problem in problems[:10]:
            print(f"       {problem}")
        if len(problems) > 10:
            print(f"       ... and {len(problems) - 10} more")
        failed += bool(problems)

    print(f"\n{len(files) - failed} of {len(files)} solutions passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
