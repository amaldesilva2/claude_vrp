#ifndef CLAUDE_VRP_INCLUDE_ROUTE_ELIMINATION_H_
#define CLAUDE_VRP_INCLUDE_ROUTE_ELIMINATION_H_

#include <chrono>
#include <vector>

#include "destroy_repair.h"
#include "local_search.h"
#include "solution.h"

struct RouteEliminationParams {
  int max_ejections = 2000;  // ejection steps per attempt before giving up
  int max_ejected = 2;       // customers one insertion may eject (1 or more)
};

struct RouteEliminationStats {
  int attempts = 0;  // routes we tried to eliminate
  int removed = 0;   // routes successfully eliminated
};

// Inserts the customers of `pool` (taken LIFO) into existing non-empty routes.
// A customer is inserted at its cheapest feasible position; if it fits
// nowhere it replaces a string of 1..`max_ejected` consecutive customers,
// choosing the lowest total ejection count (ties: smallest distance increase),
// and the replaced customers go back into the pool. Stops when the pool is
// empty (returns true), after `max_ejections` ejection steps, or when a
// customer fits nowhere even by replacement (returns false; the remaining
// customers stay in `pool`). If `touched` is given, the flag of every route
// that changed is set (one flag per route). With a restricted `scope`,
// replacements are first searched only in the customer's nearby routes.
bool insertWithEjections(int max_ejections, int max_ejected,
                         std::vector<int>* pool, std::vector<bool>* touched,
                         Solution* solution, const InsertionScope& scope = {});

// Tries to remove route `r` by moving its customers into the other routes
// with insertWithEjections. On failure `solution` is left unchanged.
bool eliminateRoute(int r, const RouteEliminationParams& params,
                    Solution* solution);

// Reduces the number of vehicles by repeatedly calling eliminateRoute (a
// simplified version of Nagata & Braysy, 2009). Routes are tried smallest
// first, all of them with single ejections before any with up to
// params.max_ejected; after each success
// `local_search` is run. Stops at the deadline, when every route fails, or
// when the capacity lower bound is reached.
RouteEliminationStats minimizeRoutes(
    const RouteEliminationParams& params,
    std::chrono::steady_clock::time_point deadline, LocalSearch* local_search,
    Solution* solution);

#endif  // CLAUDE_VRP_INCLUDE_ROUTE_ELIMINATION_H_
