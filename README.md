# claude_vrp

A heuristic solver for the **Vehicle Routing Problem with Time Windows
(VRPTW)** in C++17. The objective is hierarchical: first the number of
vehicles, then the total travel distance (plain double Euclidean distances).

It is tuned for the large **Gehring & Homberger 1000-customer** benchmark
instances (included in `data/`).

## Who wrote it

This code was written by **Claude** (Anthropic's AI model, Claude Opus 5.5,
working in Claude Code), with guidance from **Amal de Silva**.

* **Claude** designed and implemented the solver, read the research
  literature and implemented the methods (SISR, Nagata & Braysy route
  minimization, ALNS, hybrid genetic search), built the benchmark tools, and
  ran and evaluated the experiments that set the defaults (see
  [`experiments/RESEARCH_LOG.md`](experiments/RESEARCH_LOG.md)).
* **Amal de Silva** set the goals and direction: the problem and benchmark
  instances, the coding style, which ideas to pursue next, and reviewed the
  results along the way.

*claude_vrp is an independent project. It is not affiliated with, sponsored
by, or endorsed by Anthropic. Claude is a trademark of Anthropic, PBC.*

## Results

120 s per instance, 8 threads, all 60 Gehring & Homberger 1000-customer
instances, compared with the best known solutions (BKS) from
[SINTEF](https://www.sintef.no/projectweb/top/vrptw/homberger-benchmark/1000-customers/)
(full table: `benchmark_results/homberger1000_120s_8threads.csv`):

| type | vehicles over BKS | instances with BKS vehicle count | distance gap on those |
|---|---|---|---|
| C1 | +6 | 6/10 | +0.53% |
| C2 | +5 | 5/10 | +1.00% |
| R1 | 0 | 10/10 | +4.12% |
| R2 | 0 | 10/10 | +1.85% |
| RC1 | 0 | 10/10 | +2.11% |
| RC2 | 0 | 10/10 | +2.28% |
| **all** | **+11** | **51/60** | **+2.19%** |

## Build

Requires CMake >= 3.16 and a C++17 compiler.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Usage

```bash
./build/claude_vrp [instance] [seconds] [--seed N] [--threads N] \
    [--search alns|ils|hgs|mixed] [--output FILE] [--set key=value ...]

./build/claude_vrp data/RC2_10_9.TXT             # defaults: 120 s, 8 threads, ALNS
./build/claude_vrp data/C1_10_2.TXT 300 --threads 16
./build/claude_vrp data/R1_10_1.TXT 60 --threads 1 --seed 3
./build/claude_vrp data/RC2_10_9.TXT --output RC2_10_9.sol   # write the routes
```

Every final solution is validated (each customer exactly once, capacity, time
windows, return to the depot, fleet size, and a fresh recomputation of the
distance). `--output FILE` writes the routes in the format of the SINTEF best
known solutions (`Route i : c1 c2 ...`), so they can be checked
independently, e.g. with `experiments/verify_solutions.py` (pure Python, no
code shared with the solver).

`--set key=value` changes any algorithm parameter without recompiling; the
keys are listed in `applySetting()` in `src/main.cc`.

## How it works

1. **Construction:** Solomon's I1 insertion heuristic over a small parameter grid.
2. **Local search:** granular neighborhoods (relocate, swap, 2-opt, 2-opt*,
   moves of strings of customers); feasibility and cost of every move in O(1)
   via segment concatenation (Vidal et al. 2013).
3. **Vehicle minimization** (30% of the time): ejection-based route
   elimination, then SISR fleet minimization with absence counters
   (Christiaens & Vanden Berghe 2020) or, for capacity-tight fleets, the route
   minimization heuristic of Nagata & Braysy (2009). With several threads,
   both strategies run side by side.
4. **Distance minimization:** adaptive large neighborhood search without local
   search in the inner loop (SISR string removal, random / worst / related /
   route / radial removal, greedy and regret insertion with blinks, simulated
   annealing), with periodic local search and route elimination. Threads share
   their best solution every 10% of the run.

Also available: iterated local search (`--search ils`), a hybrid genetic search
with SREX crossover and penalized local search (`--search hgs`).

The design decisions and all experiments behind the defaults are documented
in [`experiments/RESEARCH_LOG.md`](experiments/RESEARCH_LOG.md).

## Benchmarking

```bash
python3 run_heuristic_benchmark.py                 # all instances, 120 s, 8 threads
python3 run_heuristic_benchmark.py --pattern 'RC2*' --seconds 60
python3 run_heuristic_benchmark.py --solutions solutions/   # also save the routes
python3 experiments/verify_solutions.py solutions/          # check them independently
python3 experiments/run_bench.py my-config -- --set alns.max_remove=40
python3 experiments/compare.py my-config other-config
```

## Repository layout

```
include/, src/             solver (C++17)
data/                      Gehring & Homberger 1000-customer instances
experiments/               tuning tools, experiment results, research log,
                           best known solutions, independent solution checker
benchmark_results/         benchmark output (CSV)
run_heuristic_benchmark.py solve instances and compare with the best known
                           solutions (CSV + per-type summary)
```

## References

* M. M. Solomon (1987). Algorithms for the vehicle routing and scheduling
  problems with time window constraints. *Operations Research* 35(2).
* H. Gehring, J. Homberger (1999). A parallel hybrid evolutionary metaheuristic
  for the vehicle routing problem with time windows.
* T. Vidal, T. G. Crainic, M. Gendreau, C. Prins (2013). A hybrid genetic
  algorithm with adaptive diversity management for a large class of vehicle
  routing problems with time-windows. *Computers & Operations Research* 40(1).
* Y. Nagata, O. Braysy (2009). A powerful route minimization heuristic for the
  vehicle routing problem with time windows. *Operations Research Letters* 37(5).
* S. Ropke, D. Pisinger (2006). An adaptive large neighborhood search heuristic
  for the pickup and delivery problem with time windows. *Transportation
  Science* 40(4).
* J. Christiaens, G. Vanden Berghe (2020). Slack induction by string removals
  for vehicle routing problems. *Transportation Science* 54(2).
* W. Kool et al. (2022). Hybrid genetic search for the vehicle routing problem
  with time windows: a high-performance implementation. 12th DIMACS
  Implementation Challenge.

## License

[MIT](LICENSE). The benchmark instances in `data/` are the public
Gehring & Homberger instances and are included for convenience.
