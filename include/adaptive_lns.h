#ifndef CLAUDE_VRP_INCLUDE_ADAPTIVE_LNS_H_
#define CLAUDE_VRP_INCLUDE_ADAPTIVE_LNS_H_

#include <chrono>
#include <string>
#include <vector>

#include "local_search.h"
#include "solution.h"
#include "solution_exchange.h"

struct AlnsParams {
  // Customers removed per iteration, uniform in [min_remove, max_remove].
  int min_remove = 10;
  int max_remove = 30;
  // Rank exponent for worst and related removal (larger = greedier).
  double removal_randomness = 6.0;
  // SISR string removal settings (see destroy_repair.h).
  double avg_removed = 15.0;
  int max_string = 10;
  double split_probability = 0.5;  // SISR split-string variant
  double blink_rate = 0.01;        // SISR blinks in repair

  // Repair tries each customer only in routes containing one of its this
  // many nearest customers (falling back to all routes); 0 = all routes.
  int insertion_neighbors = 40;

  // Run local search on the changed routes after every repair.
  bool local_search = false;  // off: many fast iterations work better
  // Every this many iterations, run local search on the current solution
  // (0 disables it). Independent of local_search.
  int ls_period = 2000;

  // Leftover customers are inserted again after local search with up to this
  // many single ejections before they get new routes.
  int max_pool_ejections = 100;

  // Adaptive weights (Ropke & Pisinger, 2006): every segment_length
  // iterations, weight = (1 - reaction) * weight + reaction * average score,
  // with the scores below per use, and never below min_weight.
  int segment_length = 100;
  double reaction = 0.1;
  double score_new_best = 33.0;
  double score_better = 9.0;  // better than the current solution
  double score_accepted = 13.0;
  double min_weight = 0.1;

  // Simulated annealing: the start temperature accepts a solution
  // start_worse (relative) worse with probability 0.5; it decays
  // exponentially over the run to end_temperature_ratio times that.
  double start_worse = 0.01;
  double end_temperature_ratio = 0.01;

  // Every this many iterations, try to eliminate a route (0 disables it).
  int elimination_period = 5000;
  int max_elimination_ejections = 2000;
  int max_elimination_ejected = 2;  // escalates from 1 up to this

  // Cooperation between parallel searches: every exchange_seconds the best
  // solution is offered to `exchange`, and a better one from another search
  // replaces both the best and the current solution (nullptr or <= 0: off).
  SolutionExchange* exchange = nullptr;
  double exchange_seconds = 0.0;

  unsigned seed = 1;
};

struct AlnsOperatorStats {
  std::string name;
  int uses = 0;
  int new_best = 0;
  double weight = 1.0;  // at the end of the run
};

struct AlnsStats {
  int iterations = 0;
  int accepted = 0;
  int improvements = 0;  // new best solutions
  int elimination_attempts = 0;
  int routes_eliminated = 0;
  int adopted = 0;  // better solutions taken from other searches
  std::vector<AlnsOperatorStats> destroy;
  std::vector<AlnsOperatorStats> repair;
};

// Adaptive large neighborhood search. Each iteration picks a destroy and a
// repair operator by roulette-wheel selection on adaptive weights:
//   destroy: random, worst, related (Shaw), route, radial, string (SISR);
//   repair:  greedy (random order), greedy (earliest deadline first),
//            regret-2, regret-3.
// It removes min_remove..max_remove customers, repairs, optionally runs
// `local_search` on the changed routes, reinserts leftovers with ejections
// (new routes only as a last resort), and accepts by simulated annealing on
// distance, after comparing vehicles: fewer is always accepted, more never.
// Periodically it tries to eliminate a route. `solution` receives the best
// solution found; with local_search on, it must be a local optimum on entry.
AlnsStats adaptiveLns(const AlnsParams& params,
                      std::chrono::steady_clock::time_point deadline,
                      LocalSearch* local_search, Solution* solution);

#endif  // CLAUDE_VRP_INCLUDE_ADAPTIVE_LNS_H_
