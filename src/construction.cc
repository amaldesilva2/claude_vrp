#include "construction.h"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Index into `unrouted` of the seed customer for a new route.
std::size_t pickSeed(const Instance& inst, const std::vector<int>& unrouted,
                     SeedRule rule) {
  std::size_t best = 0;
  for (std::size_t k = 1; k < unrouted.size(); ++k) {
    const int u = unrouted[k];
    const int b = unrouted[best];
    const bool better = rule == SeedRule::kFarthest
                            ? inst.distance(0, u) > inst.distance(0, b)
                            : inst.nodes[u].due < inst.nodes[b].due;
    if (better) best = k;
  }
  return best;
}

// Finds the feasible position in `route` with the lowest Solomon c1 score for
// `customer`:
//   c11 = d(i,u) + d(u,j) - mu * d(i,j)   (extra distance)
//   c12 = pushForward                     (extra time at j)
//   c1  = alpha1 * c11 + (1 - alpha1) * c12
// Returns false if the customer fits nowhere in the route.
bool bestInsertion(const Instance& inst, const Route& route, int customer,
                   const I1Params& params, int* best_pos, double* best_c1) {
  const std::vector<int>& nodes = route.nodes();
  *best_pos = -1;
  *best_c1 = kInf;
  for (int pos = 0; pos <= route.numCustomers(); ++pos) {
    if (!route.canInsert(customer, pos)) continue;
    const int i = nodes[pos];
    const int j = nodes[pos + 1];
    const double c11 = inst.distance(i, customer) + inst.distance(customer, j) -
                       params.mu * inst.distance(i, j);
    const double c12 = route.pushForward(customer, pos);
    const double c1 = params.alpha1 * c11 + (1.0 - params.alpha1) * c12;
    if (c1 < *best_c1) {
      *best_c1 = c1;
      *best_pos = pos;
    }
  }
  return *best_pos >= 0;
}

// Removes unrouted[k] in O(1); order of `unrouted` does not matter.
void removeAt(std::vector<int>* unrouted, std::size_t k) {
  (*unrouted)[k] = unrouted->back();
  unrouted->pop_back();
}

}  // namespace

Solution solomonI1(const Instance& inst, const I1Params& params) {
  Solution solution(inst);

  std::vector<int> unrouted;
  unrouted.reserve(inst.numCustomers());
  for (int c = 1; c < inst.numNodes(); ++c) unrouted.push_back(c);

  while (!unrouted.empty()) {
    // Open a new route with a seed customer.
    const std::size_t seed_index = pickSeed(inst, unrouted, params.seed);
    const int seed = unrouted[seed_index];
    Route& route = solution.addRoute();
    if (!route.canInsert(seed, 0)) {
      throw std::runtime_error("customer " + std::to_string(seed) +
                               " cannot be served on its own route");
    }
    route.insert(seed, 0);
    removeAt(&unrouted, seed_index);

    // Fill it: insert the customer with the highest c2 until none fits.
    //   c2 = lambda * d(0,u) - c1(u)
    while (true) {
      std::size_t best_index = 0;
      int best_pos = -1;
      double best_c2 = -kInf;
      for (std::size_t k = 0; k < unrouted.size(); ++k) {
        const int u = unrouted[k];
        int pos;
        double c1;
        if (!bestInsertion(inst, route, u, params, &pos, &c1)) continue;
        const double c2 = params.lambda * inst.distance(0, u) - c1;
        if (c2 > best_c2) {
          best_c2 = c2;
          best_index = k;
          best_pos = pos;
        }
      }
      if (best_pos < 0) break;
      route.insert(unrouted[best_index], best_pos);
      removeAt(&unrouted, best_index);
    }
  }

  return solution;
}

int insertCheapest(int customer, Solution* solution) {
  std::vector<Route>& routes = *solution->mutableRoutes();
  int best_route = -1;
  int best_pos = -1;
  double best_cost = kInf;
  for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
    const Route& route = routes[r];
    if (route.empty()) continue;
    for (int pos = 0; pos <= route.numCustomers(); ++pos) {
      if (!route.canInsert(customer, pos)) continue;
      const double cost = route.insertionCost(customer, pos);
      if (cost < best_cost) {
        best_cost = cost;
        best_route = r;
        best_pos = pos;
      }
    }
  }
  if (best_route >= 0) routes[best_route].insert(customer, best_pos);
  return best_route;
}
