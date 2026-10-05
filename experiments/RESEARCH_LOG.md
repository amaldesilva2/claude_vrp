# Research log: improving the VRPTW heuristic

All experiments use `experiments/run_bench.py` (results appended to
`experiments/results.csv`, compare with `experiments/compare.py CONFIG...`).

Protocol unless stated otherwise: **core** set of 13 Gehring & Homberger
1000-customer instances (2 per type + RC2_10_9), single thread per run,
30 s per run, 2 or 3 seeds, 16 runs in parallel.
Ranking is lexicographic per instance (vehicles, then distance); "gap" is the
distance gap to the best known solution (BKS, `experiments/bks.txt`, SINTEF).

## Starting point

| config | vehicles over BKS | mean gap |
|---|---|---|
| ILS (local search + ruin/recreate), old defaults | +9.5 | +6.9% |
| ALNS (with local search after every repair) | +11.0 | +6.8% |

Weakest spots: vehicles on C1/C2 (C1_10_2: 94 vs 90), distance on RC1 (+22-26%)
and R1_10_5 (+15%).

## What was tried

| # | idea | result | kept? |
|---|---|---|---|
| 1 | `-O3 -march=native` | no measurable speedup | no |
| 2 | smaller local search move set (runs <= 2, no reversed runs, K=20 neighbors) | 2-3x faster LS, better ILS (mean rank 1.77 vs 4.15) | yes |
| 3 | soft constraints in local search (penalized time warp / excess load, HGS-style adaptive penalties, x10/x100 repair) | ILS clearly worse (11 of 13 instances lost): candidates often stay infeasible, repair is expensive | option only (`--set soft.tw=...`) |
| 4 | **ALNS without local search** (SISR philosophy: many cheap ruin & recreate iterations + SA) | **gap +10.4% -> +4.0%**, 2000 -> 344 000 iterations in 30 s | yes |
| 5 | lazy route segment summaries | +14% iterations | yes |
| 6 | insertion restricted to routes near the customer (fallback: all routes) | +60% iterations; 40 neighbors better than 20 | yes (40) |
| 7 | SA start temperature / removal size / end temperature variations | defaults were best | - |
| 8 | **SISR fleet minimization** (absence counters) after ejection-based route elimination | vehicles +8.0 -> **+6.5**; RC2_10_1 reaches BKS 20 vehicles | yes |
| 9 | ALNS as default instead of ILS (both with fleet minimization) | gap +4.4% vs +11.0% | yes |
| 10 | final local-search polish of the best ALNS solution | small (0.01-0.26%) but free | yes |
| 11 | local search on the current ALNS solution every 2000 iterations | slight gain (mean rank 2.15 vs 2.46) | yes |
| 12 | SISR blinks (1%) and split strings in ALNS | neutral within noise; kept for fidelity to SISR | yes |
| 13 | elimination share 0.5 instead of 0.3 | equal | no |
| 14 | Nagata & Braysy route minimization (random insertion, squeeze with penalized LS, k<=5 ejection DFS with counters, 1000 random perturbation moves) | at 30 s: worse (eats the elimination budget); alone on C1_10_2 for ~108 s: **90 vehicles = BKS** | see below |
| 15 | ALNS leftover ejections restricted to nearby routes | better distance, but C instances lost vehicles | yes, combined with 16 |
| 16 | ALNS periodic route elimination every 5000 iterations | with 15: vehicles +6.67, gap **+2.91%** (from +4.10%) | yes |
| 17 | cooperation between threads (share best solution every 10% of the time) | 8-thread runs: wins 5 of 6 instances vs independent threads | yes |
| 18 | Nagata & Braysy search at 120 s (alone, 30% / 50% of time) | complementary: better on C1_10_2 (92.0-92.5 vs 93.5 veh) but loses RC2_10_1 (21 vs 20) and distance elsewhere | option (`--set nb.enabled=1`), off |
| 19 | fleet minimization with stagnation stop + Nagata & Braysy after it | worse (RC2_10_1 needs long fleet-minimization runs) | option (`fleet.max_stagnation`), off |
| 20 | **Hybrid genetic search** (`--search hgs`): feasible/infeasible populations, biased fitness with broken-pairs diversity, SREX crossover, penalized local search education, repair, adaptive penalties, fixed fleet | 30 s: far worse (+21% gap, ~1400 offspring); 120 s: RC2_10_9 23 874 vs ALNS 23 567, R1_10_1 57 139 vs 54 493 (only 5% of offspring feasible, time-warp penalty rose to 1300 with the tight fleet) | option, off |
| 21 | ALNS for 70% of the time, then HGS from its result (`--set hgs.alns_share=0.7`) | ties ALNS (quick set, 30 s: 3 wins each, gap 3.14% vs 3.01%); 120 s RC2_10_9 23 540 vs 23 567, R1_10_1 no feasible offspring | option, off |
| 22 | fleet-minimization stagnation streaks measured: RC2_10_1 needs up to ~97 000 iterations without progress before a success | a stagnation-based switch to Nagata & Braysy would cut it off | - |
| 23 | **elimination portfolio across threads**: odd threads use Nagata & Braysy instead of ejection elimination + fleet minimization; fleets spread via the exchange | 8 threads, 120 s: C1_10_2 **91** vs 92 vehicles, others equal; 30 s: neutral on vehicles | yes (default, `elim.portfolio`) |
| 24 | **capacity-tight rule**: if the fleet after ejection elimination is within 10% of the capacity lower bound, use Nagata & Braysy instead of fleet minimization (5% never triggered: C1_10_2 is at 96 vs 90 at that point) | 1 thread, 120 s: C1_10_2 **91.5** vs 93 vehicles, others equal; 30 s: neutral (same vehicles, 4 of 6 slightly better distance) | yes (default, `elim.auto_gap=0.1`) |

## Notes

* Noise: with 2-3 seeds per instance, differences below ~0.3% gap or one
  vehicle on one seed are not significant.
* Vehicles dominate the ranking; distance gaps are only comparable at equal
  vehicle counts (`compare.py` reports the gap on instances where a config
  matches the best vehicle count).

* Results are only compared between runs made under the same machine load;
  a fleet-minimization stagnation limit of 100 000 iterations was tested
  under different load conditions and is inconclusive (off by default).

## Results of the default configuration

`./build/claude_vrp <instance> [seconds] [--threads N]` runs: I1 -> local
search -> ejection-based route elimination -> SISR fleet minimization or, for
capacity-tight fleets and on odd threads, Nagata & Braysy route minimization
(30% of the time together) -> ALNS without local search (SISR ruin/recreate,
regret repair, SA, periodic local search and route elimination, thread
cooperation) -> final local search.

Results on the 13 core instances (before experiments 23-24):

| setting | vehicles over BKS (13 core instances) | mean gap |
|---|---|---|
| starting point: ILS, 30 s, 1 thread | +9.5 | +6.9% |
| final, 30 s, 1 thread (3 seeds) | **+6.67** | **+2.91%** |
| final, 120 s, 1 thread (2 seeds) | +6.5 | **+1.84%** |
| final, 30 s, 8 threads with exchange (6 instances) | +3.0 on 6 inst. | +2.22% |

Per instance (30 s, 1 thread, mean of 3 seeds) the remaining gaps are mostly
vehicles on C1/C2 (C1_10_2 93.3 vs 90, C1_10_6 100 vs 99, C2 +1) and 3-6%
distance on RC1 / R1_10_5. RC2_10_9: 18 vehicles, ~23 600-23 900 (BKS 22 919,
gap ~3-4%).

All 60 instances with the final defaults (120 s, 8 threads): 51 of 60 at the
best known vehicle count (+11 vehicles in total, all on C1/C2), mean distance
gap +2.19% on those 51 instances; see `benchmark_results/` and the README.

## Ideas not tried yet

* HGS needs much longer runs than 2 minutes on 1000 customers and struggles
  with fleets fixed at the minimized size; it might pay off with an unlimited
  fleet and a separate vehicle penalty, or for runs of 10+ minutes.
* EAX crossover (Nagata) for distance on R1/RC1.
