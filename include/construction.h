#ifndef CLAUDE_VRP_INCLUDE_CONSTRUCTION_H_
#define CLAUDE_VRP_INCLUDE_CONSTRUCTION_H_

#include "instance.h"
#include "solution.h"

// How a new route picks its first customer.
enum class SeedRule {
  kFarthest,          // unrouted customer farthest from the depot
  kEarliestDeadline,  // unrouted customer with the smallest due time
};

// Parameters of Solomon's I1 insertion heuristic.
struct I1Params {
  double mu = 1.0;      // weight of the saved direct arc i -> j in c11
  double lambda = 1.0;  // weight of the depot distance in c2
  double alpha1 = 1.0;  // weight of distance (c11); time (c12) gets 1 - alpha1
  SeedRule seed = SeedRule::kFarthest;
};

// Builds a solution with Solomon's (1987) I1 insertion heuristic: routes are
// opened one at a time with a seed customer and filled by repeatedly inserting
// the customer with the best c2 score until none fits. Throws
// std::runtime_error if a seed customer cannot be served on its own route.
Solution solomonI1(const Instance& inst, const I1Params& params);

// Inserts `customer` at its cheapest feasible position in any non-empty route
// of `solution`. Returns the route index, or -1 if it fits nowhere.
int insertCheapest(int customer, Solution* solution);

#endif  // CLAUDE_VRP_INCLUDE_CONSTRUCTION_H_
