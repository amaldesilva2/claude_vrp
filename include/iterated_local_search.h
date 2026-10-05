#ifndef CLAUDE_VRP_INCLUDE_ITERATED_LOCAL_SEARCH_H_
#define CLAUDE_VRP_INCLUDE_ITERATED_LOCAL_SEARCH_H_

#include <chrono>

#include "local_search.h"
#include "solution.h"

struct IlsParams {
  // Radial ruin: a random customer and its nearest neighbors.
  int min_ruin = 10;
  int max_ruin = 40;
  // String ruin (SISR, Christiaens & Vanden Berghe 2020): remove strings of
  // consecutive customers from routes near a random customer.
  double string_ruin_probability = 0.5;
  double avg_removed = 15.0;  // c-bar: average number of customers removed
  int max_string = 10;        // L-max: longest string removed from one route

  // Customers that recreate cannot place are inserted again after local
  // search, with up to this many ejection steps, each ejecting a string of up
  // to max_pool_ejected customers.
  int max_pool_ejections = 100;
  int max_pool_ejected = 1;

  // If true, customers that still do not fit stay unplaced in a pool and cost
  // unplaced_penalty * 2 * d(0, c) each; otherwise each gets a new route.
  bool penalize_unplaced = false;  // off: hurt RC1 by ~4% in tests
  double unplaced_penalty = 1.0;
  // While the current state has unplaced customers the weight is multiplied
  // by this every iteration; it is reset once the state is complete.
  double penalty_growth = 1.05;
  // If the current state has had unplaced customers for this many
  // consecutive iterations (e.g. after a soft elimination), it is reset to
  // the best solution.
  int max_incomplete_iterations = 50;
  // If true (and penalize_unplaced), a route that cannot be eliminated with
  // ejections is moved into the pool instead (soft elimination).
  bool soft_elimination = false;

  // Every this many iterations, try to eliminate a route of the current
  // solution (0 disables it).
  int elimination_period = 100;
  int max_elimination_ejections = 2000;
  int max_elimination_ejected = 2;  // escalates from 1 up to this

  double initial_threshold = 0.01;  // accept up to +1% over best, decays to 0
  unsigned seed = 1;
};

struct IlsStats {
  int iterations = 0;
  int accepted = 0;      // candidates that became the current solution
  int improvements = 0;  // candidates that became the best solution
  int string_ruins = 0;  // iterations that used string removal
  int pool_rescues = 0;  // leftover customers placed after local search
  int new_routes = 0;    // leftover customers that needed a new route
  int elimination_attempts = 0;
  int routes_eliminated = 0;    // removed directly with ejections
  int soft_eliminations = 0;    // routes moved into the pool
  int incomplete_restarts = 0;  // current reset to best: pool not emptied
};

// Iterated local search. The search state is a solution plus a pool of
// unplaced customers. Each iteration copies the current state and:
//   1. ruin:     removes customers, either radially (a random customer and its
//                nearest neighbors) or as SISR strings;
//   2. recreate: reinserts them and the pool at their cheapest feasible
//                positions, in a random order chosen from random / farthest
//                from depot first / largest demand first / earliest due time
//                first;
//   3. runs `local_search` on the changed routes, then reinserts leftovers
//      with ejections (route_elimination.h); customers that still do not fit
//      stay in the pool at a penalty that grows while the state is
//      incomplete (or, if penalize_unplaced is false, get
//      new routes), followed by another local search.
// Periodically it tries to eliminate a route of the current state: directly
// with ejections, or else by moving its customers into the pool. A state
// that stays incomplete for max_incomplete_iterations is reset to the best.
// Acceptance is lexicographic (vehicles, then penalized cost) with a
// record-to-record threshold: a candidate with the same vehicle count is
// accepted if its cost is below the current one or within threshold * best
// cost of the best, the threshold decaying linearly to 0 at the deadline.
// Only complete states (empty pool) can become the best. `solution` must be
// a local optimum on entry and receives the best solution found.
IlsStats iteratedLocalSearch(const IlsParams& params,
                             std::chrono::steady_clock::time_point deadline,
                             LocalSearch* local_search, Solution* solution);

#endif  // CLAUDE_VRP_INCLUDE_ITERATED_LOCAL_SEARCH_H_
