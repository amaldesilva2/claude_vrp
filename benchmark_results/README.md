# Benchmark results

## `homberger1000_120s_8threads.csv`

All 60 Gehring & Homberger 1000-customer instances (`data/`), one run each,
compared with the best known solutions in `experiments/bks.txt` (SINTEF).

| type | vehicles over BKS | instances at BKS vehicle count | distance gap on those |
|---|---|---|---|
| C1 | +6 | 6/10 | +0.53% |
| C2 | +5 | 5/10 | +1.00% |
| R1 | 0 | 10/10 | +4.12% |
| R2 | 0 | 10/10 | +1.85% |
| RC1 | 0 | 10/10 | +2.11% |
| RC2 | 0 | 10/10 | +2.28% |
| **all** | **+11** | **51/60** | **+2.19%** |

All 60 solutions were validated by the solver (each customer served exactly
once, vehicle capacity, time windows, fleet size).

### Run settings

* Solver: default configuration (ALNS, vehicle minimization with the
  elimination portfolio and capacity-tight rule, thread cooperation), no
  `--set` overrides.
* 120 s per instance, 8 threads per instance, 2 instances solved at the same
  time (16 threads on 18 cores):

  ```bash
  python3 run_heuristic_benchmark.py --seconds 120 --threads 8 --jobs 2
  ```

* Seeds: 1-8 (one per thread).

### Hardware and software

| | |
|---|---|
| Machine | MacBook Pro (Mac17,8) |
| CPU | Apple M5 Pro, 18 cores (6 super + 12 performance) |
| Memory | 24 GB |
| OS | macOS 26.5.1 |
| Compiler | Apple clang 21.0.0, `-O2`, C++17 |
| Build | CMake 4.2.3, `CMAKE_BUILD_TYPE=Release` |

### Notes

* The time limit is wall-clock time and includes reading the instance and
  the construction heuristic, so results depend on the machine's speed and
  load; a faster or less loaded machine gets more iterations in the same time.
* The search is randomized and parallel, so repeated runs give slightly
  different results even with the same seeds.
* Distances are plain double Euclidean distances (no rounding), as in the
  SINTEF best known solutions.

## Columns

| column | meaning |
|---|---|
| `instance`, `type` | instance name and its type (C1, C2, R1, R2, RC1, RC2) |
| `vehicles`, `distance` | our solution |
| `bks_vehicles`, `bks_distance` | best known solution |
| `extra_vehicles` | `vehicles - bks_vehicles` |
| `gap_percent` | `100 * (distance / bks_distance - 1)`; only comparable when `extra_vehicles` is 0 (more vehicles usually allow a shorter distance) |
| `valid` | `yes` if the solver validated the solution |
| `wall_seconds` | wall-clock time of the run |
