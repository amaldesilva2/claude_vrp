#ifndef CLAUDE_VRP_INCLUDE_DESTROY_REPAIR_H_
#define CLAUDE_VRP_INCLUDE_DESTROY_REPAIR_H_

#include <random>
#include <vector>

#include "instance.h"
#include "solution.h"

// Destroy (ruin) and repair (recreate) operators shared by the iterated local
// search and the adaptive large neighborhood search.
//
// Destroy operators only choose customers; removeCustomers() takes them out.
// Repair operators insert customers into existing non-empty routes and return
// the ones that fit nowhere. `touched` holds one flag per route and is set for
// every route that changed.

// Customer -> (route, node index) for the current state of a solution;
// -1 for customers that are not in any route.
struct Positions {
  std::vector<int> route_of;
  std::vector<int> index_of;

  explicit Positions(const Solution& solution);
};

// ---------------------------------------------------------------------------
// Destroy operators
// ---------------------------------------------------------------------------

// A random customer and its count - 1 nearest neighbors (`nearest` from
// nearestCustomers()). May include customers that are not routed.
std::vector<int> radialRemoval(const Instance& inst,
                               const std::vector<std::vector<int>>& nearest,
                               int count, std::mt19937* rng);

// SISR string removal (Christiaens & Vanden Berghe, 2020): from the routes of
// the customers nearest to a random seed customer, remove one string of
// consecutive customers each. About `avg_removed` customers in total, strings
// at most `max_string` long.
//
// With probability `split_probability` a route gets the "split string"
// variant: a window of len + m customers around the chosen one, of which a
// random block of m consecutive customers is kept (m grows by one while a
// uniform draw exceeds 0.01, up to the route size).
std::vector<int> stringRemoval(const Solution& solution,
                               const std::vector<std::vector<int>>& nearest,
                               double avg_removed, int max_string,
                               std::mt19937* rng,
                               double split_probability = 0.0);

// `count` routed customers chosen uniformly at random.
std::vector<int> randomRemoval(const Solution& solution, int count,
                               std::mt19937* rng);

// Customers whose removal saves the most distance (d(a,u) + d(u,b) - d(a,b)).
// The k-th pick is taken at rank floor(y^randomness * size) of the remaining
// candidates, y uniform in [0, 1): larger `randomness` is greedier.
std::vector<int> worstRemoval(const Solution& solution, int count,
                              double randomness, std::mt19937* rng);

// Shaw removal (Ropke & Pisinger, 2006): starting from a random customer,
// repeatedly adds a customer related to an already removed one, where
// relatedness combines distance, service start time and demand. Candidates
// come from `nearest`; ranks are randomized as in worstRemoval.
std::vector<int> relatedRemoval(const Solution& solution,
                                const std::vector<std::vector<int>>& nearest,
                                int count, double randomness,
                                std::mt19937* rng);

// All customers of randomly chosen routes, smaller routes more likely, until
// at least `count` customers are selected (at least one route).
std::vector<int> routeRemoval(const Solution& solution, int count,
                              std::mt19937* rng);

// Removes `removed` (customers that are not routed are ignored).
void removeCustomers(const std::vector<int>& removed,
                     std::vector<bool>* touched, Solution* solution);

// ---------------------------------------------------------------------------
// Repair operators
// ---------------------------------------------------------------------------

// Optional restriction of repair to "nearby" routes: a customer is only
// tried in routes that contain one of its `neighbors` nearest customers
// (from `nearest`), falling back to all routes if none of those fits it.
// Default: all routes.
struct InsertionScope {
  const std::vector<std::vector<int>>* nearest = nullptr;
  int neighbors = 0;
  // SISR "blinks": each candidate position in nearby routes is skipped with
  // this probability (needs `rng`).
  double blink_rate = 0.0;
  std::mt19937* rng = nullptr;

  bool restricted() const { return nearest != nullptr && neighbors > 0; }
};

// Distinct routes (per `route_of`, -1 = unrouted) of the first
// `scope.neighbors` nearest customers of `customer`.
void nearbyRoutes(int customer, const InsertionScope& scope,
                  const std::vector<int>& route_of, std::vector<int>* routes);

enum class InsertionOrder {
  kRandom,
  kFarthestFirst,  // farthest from the depot first
  kDemandFirst,    // largest demand first
  kDeadlineFirst,  // earliest due time first
};

// Cheapest feasible insertion of each customer, in the given order (ties
// broken by a random shuffle). Returns the customers that fit nowhere.
std::vector<int> greedyRepair(std::vector<int> customers, InsertionOrder order,
                              std::mt19937* rng, std::vector<bool>* touched,
                              Solution* solution,
                              const InsertionScope& scope = {});

// Regret-k insertion: repeatedly inserts the customer with the largest
// regret, sum over h = 2..k of (h-th best - best) insertion cost over
// different routes, at its best position. Customers with fewer than k
// feasible routes get priority (fewest options first). Returns the customers
// that fit nowhere.
std::vector<int> regretRepair(std::vector<int> customers, int k,
                              std::vector<bool>* touched, Solution* solution,
                              const InsertionScope& scope = {});

#endif  // CLAUDE_VRP_INCLUDE_DESTROY_REPAIR_H_
