# CLAUDE.md

Guidance for Claude Code (and humans) working on this repository.

## What this is

A C++17 heuristic for the vehicle routing problem with time windows (VRPTW).
Objective is **lexicographic: fewest vehicles first, then shortest distance**
(plain double Euclidean distances, travel time = distance). Tuned for the
Gehring & Homberger 1000-customer instances in `data/`; best known solutions
are in `experiments/bks.txt`. Read `README.md` for the algorithm overview and
`experiments/RESEARCH_LOG.md` for what has been tried and why the defaults are
what they are **before proposing changes**.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
ctest --test-dir build --output-on-failure      # ~10 s, must pass
./build/claude_vrp data/RC2_10_9.TXT 30         # default is 120 s, 8 threads
```

* `tests/solver_tests.cc`: unit tests (no framework; `CHECK`, one ctest entry
  per test, `./build/claude_vrp_tests <name>` runs one). They compare the O(1)
  move/insertion evaluations with full recomputation and run every search
  phase with `LocalSearchParams::check_moves`, which throws if an applied move
  was evaluated wrongly.
* `end_to_end_*` ctest entries solve an instance with `--output` and check it
  with `experiments/verify_solutions.py` (independent Python checker).
* Debug build: `cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug`.
* Thread safety (after touching threading or shared state):
  `clang++ -std=c++17 -O1 -g -fsanitize=thread -Iinclude src/*.cc -o /tmp/vrp_tsan`
  then run it with `--threads 4` long enough to reach the ALNS phase (TSan is
  10-20x slower).

## Layout

* `include/`, `src/`: solver. `main.cc` (CLI, pipeline, threads, `--set`
  parameters in `applySetting()`), `instance` (parsing), `solution` (routes,
  cached schedules, `validate()`), `segment.h` (O(1) concatenation),
  `local_search`, `construction` (Solomon I1), `destroy_repair` (ruin and
  recreate operators), `adaptive_lns` (default search), `route_elimination`,
  `fleet_minimization` (SISR), `ejection_search` (Nagata & Braysy),
  `iterated_local_search`, `genetic_search` (HGS), `solution_exchange`.
* `experiments/`: `run_bench.py` (tuning runs, appends to `results.csv`),
  `compare.py`, `verify_solutions.py`, `bks.txt`, `RESEARCH_LOG.md`.
* `run_heuristic_benchmark.py` + `benchmark_results/`: full benchmark runs.

## Coding conventions

* Google C++ style formatting (`.clang-format`; run
  `clang-format -i include/*.h src/*.cc tests/*.cc`) with these naming rules,
  enforced by `.clang-tidy`:
  * functions, methods **and getters**: camelCase (`addRoute`, `numCustomers`);
  * types: PascalCase; variables, parameters, public struct fields:
    snake_case; private members: snake_case with trailing `_`;
    constants: `kName`; files: snake_case `.h` / `.cc`;
  * include guards `CLAUDE_VRP_INCLUDE_<FILE>_H_`.
* C++ **exceptions are allowed** (`std::runtime_error`, `std::invalid_argument`).
* Comments explain intent and cite the paper an algorithm comes from; keep the
  density of the surrounding code. No new dependencies beyond the standard
  library and CMake.
* Build with no warnings (`-Wall -Wextra -Wpedantic`) and no clang-tidy
  findings.

## Correctness rules

* Every best solution is checked with `Solution::validate()` (visits,
  capacity, time windows, depot horizon, fleet size, distance). Keep it that
  way when adding search code; never report or keep an unvalidated solution.
* New local search moves go through the generic piece/plan evaluation in
  `local_search.cc`; add them to the tests' coverage and run with
  `check_moves`.
* Route caches (`start_`, `latest_`, segments) are rebuilt from the node list;
  change routes only through `Route` methods (`insert`, `remove`, `assign`).

## Experiments and benchmarking

* Tuning: `python3 experiments/run_bench.py [options] NAME -- <solver args>`
  (options **before** NAME), then `python3 experiments/compare.py A B ...`.
  `--set key=value` changes parameters without recompiling.
* Rank configurations per instance **lexicographically (vehicles, then
  distance)**; a distance gap is only meaningful at equal vehicle counts.
* Results depend on the machine: time limits are wall-clock. On a new
  machine, **re-run the baseline first** and only compare runs made under the
  same conditions. Don't run other CPU-heavy work during a benchmark.
* Noise: with 2-3 seeds, differences below ~0.3% distance or one vehicle on a
  single seed are not significant. Use `--seeds 3` before adopting a change.
* Record every experiment (idea, result, kept or not) in
  `experiments/RESEARCH_LOG.md`.

## Lessons so far (details in the research log)

* Many cheap iterations beat few expensive ones: ALNS **without** local search
  in the inner loop (SISR style) is far better than ILS / ALNS with local
  search after every repair (~+4% vs ~+10% gap at 30 s).
* Vehicles are the primary objective and the main remaining gap (C1/C2).
  SISR fleet minimization and the Nagata & Braysy method are complementary;
  the defaults combine them (capacity-tight rule, thread portfolio).
* Did not help at <= 2 min: soft-constraint (penalized) local search in ILS,
  hybrid genetic search (too few offspring, infeasible offspring with a tight
  fleet), `-O3 -march=native`, larger elimination time share.
* Untried ideas: HGS with an unlimited fleet and a vehicle penalty for long
  runs, EAX crossover for distance on R1/RC1.
