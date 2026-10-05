#ifndef CLAUDE_VRP_INCLUDE_FLEET_MINIMIZATION_H_
#define CLAUDE_VRP_INCLUDE_FLEET_MINIMIZATION_H_

#include <chrono>

#include "solution.h"

struct FleetMinimizationParams {
  // String removal (SISR) settings for the ruin step.
  double avg_removed = 10.0;
  int max_string = 10;
  // Repair tries each customer only in routes containing one of its this
  // many nearest customers (falling back to all routes); 0 = all routes.
  int insertion_neighbors = 40;
  // Give up after this many iterations without fewer absent customers or a
  // removed route (0: run until the deadline).
  int max_stagnation = 0;
  unsigned seed = 1;
};

struct FleetMinimizationStats {
  int iterations = 0;
  int routes_removed = 0;  // vehicles saved
  // Longest run of iterations without fewer absent customers that still
  // ended in a removed route (helps choosing max_stagnation).
  int longest_successful_stagnation = 0;
  bool stagnated = false;  // stopped by max_stagnation
};

// SISR fleet minimization (Christiaens & Vanden Berghe, 2020). Removes the
// smallest route and puts its customers in an "absent" set. Then repeatedly:
// ruin (string removal) a copy, recreate the removed and absent customers
// greedily without opening routes; the customers that still do not fit form
// the new absent set. The copy is accepted if it has fewer absent customers
// or a smaller sum of absence counters, where every iteration increments the
// counter of each customer that is absent. When no customer is absent, the
// fleet is one vehicle smaller and the next route is removed. Stops at the
// deadline or the capacity lower bound; `solution` receives the complete
// solution with the fewest vehicles.
FleetMinimizationStats minimizeFleet(
    const FleetMinimizationParams& params,
    std::chrono::steady_clock::time_point deadline, Solution* solution);

#endif  // CLAUDE_VRP_INCLUDE_FLEET_MINIMIZATION_H_
