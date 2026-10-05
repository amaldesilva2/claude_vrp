#ifndef CLAUDE_VRP_INCLUDE_EJECTION_SEARCH_H_
#define CLAUDE_VRP_INCLUDE_EJECTION_SEARCH_H_

#include <chrono>

#include "solution.h"

struct EjectionSearchParams {
  int max_ejected = 5;            // k_max: customers ejected per insertion
  int perturbation_moves = 1000;  // I_rand: random moves after an ejection
  int neighbors = 30;             // candidate partners for random moves
  bool squeeze = true;            // try a penalized insertion + repair first
  // Penalty weight of the squeeze local search (per unit of time warp and
  // excess load); large so that distance hardly matters.
  double squeeze_penalty = 1000.0;
  int max_dfs_nodes = 2000;  // search budget per (route, position)
  unsigned seed = 1;
};

struct EjectionSearchStats {
  int routes_removed = 0;
  int attempts = 0;
  long insertions = 0;  // direct (random feasible) insertions
  long squeezes = 0;    // successful squeezes
  long ejections = 0;   // insertions with ejections
};

// Route minimization heuristic of Nagata & Braysy (2009), "A powerful route
// minimization heuristic for the vehicle routing problem with time windows".
// Repeatedly deletes a random route: its customers go to an ejection pool
// (LIFO). Each pool customer v is inserted at a random feasible position; if
// there is none, a "squeeze" inserts v at the least infeasible position and
// repairs with penalized local search; if that fails, v is inserted with up
// to max_ejected customers ejected from the same route, minimizing the sum
// of the ejected customers' counters p (p[v] is incremented each time v
// needs ejections), followed by perturbation_moves random feasible moves
// (relocate, swap, 2-opt*). A deletion succeeds when the pool is empty; at
// the deadline an unfinished deletion is undone. Stops at the deadline or at
// the capacity lower bound.
EjectionSearchStats minimizeRoutesEjectionSearch(
    const EjectionSearchParams& params,
    std::chrono::steady_clock::time_point deadline, Solution* solution);

#endif  // CLAUDE_VRP_INCLUDE_EJECTION_SEARCH_H_
