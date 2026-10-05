#include "route_elimination.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

#include "construction.h"
#include "segment.h"

namespace {

constexpr double kEps = 1e-6;

// Replaces a string of 1..`max_ejected` consecutive customers of some route
// by `customer`, choosing the feasible replacement with the lowest total
// ejection count of the removed customers (ties: smallest distance increase).
// Appends the ejected customers to `ejected` and returns the changed route,
// or returns -1 (and changes nothing) if no replacement is feasible.
// Only routes in `candidates` are searched if it is non-null.
int replaceWithEjection(int customer, int max_ejected,
                        const std::vector<int>& ejections,
                        std::vector<int>* ejected, Solution* solution,
                        const std::vector<int>* candidates = nullptr) {
  const Instance& inst = solution->instance();
  std::vector<Route>& routes = *solution->mutableRoutes();
  const Segment single = Segment::single(inst, customer);
  const int demand = inst.nodes[customer].demand;

  int best_route = -1;
  int best_k = -1;
  int best_len = 0;
  int best_ejections = std::numeric_limits<int>::max();
  double best_delta = std::numeric_limits<double>::infinity();

  const int num_routes =
      static_cast<int>(candidates ? candidates->size() : routes.size());
  for (int i = 0; i < num_routes; ++i) {
    const int r = candidates ? (*candidates)[i] : i;
    const Route& route = routes[r];
    const int last = route.numNodes() - 2;  // last customer index
    for (int k = 1; k <= last; ++k) {
      int removed_ejections = 0;
      int removed_demand = 0;
      for (int len = 1; len <= max_ejected && k + len - 1 <= last; ++len) {
        const int w = route.node(k + len - 1);
        removed_ejections += ejections[w];
        removed_demand += inst.nodes[w].demand;
        if (removed_ejections > best_ejections) break;
        if (route.load() - removed_demand + demand > inst.capacity) continue;
        const Segment seg =
            concat(inst, concat(inst, route.forward(k - 1), single),
                   route.backward(k + len));
        if (seg.time_warp > kEps) continue;
        const double delta = seg.distance - route.distance();
        if (removed_ejections < best_ejections || delta < best_delta) {
          best_route = r;
          best_k = k;
          best_len = len;
          best_ejections = removed_ejections;
          best_delta = delta;
        }
      }
    }
  }
  if (best_route < 0) return -1;

  Route& route = routes[best_route];
  std::vector<int> nodes = route.nodes();
  ejected->insert(ejected->end(), nodes.begin() + best_k,
                  nodes.begin() + best_k + best_len);
  nodes.erase(nodes.begin() + best_k + 1, nodes.begin() + best_k + best_len);
  nodes[best_k] = customer;
  route.assign(std::move(nodes));
  return best_route;
}

}  // namespace

bool insertWithEjections(int max_ejections, int max_ejected,
                         std::vector<int>* pool, std::vector<bool>* touched,
                         Solution* solution, const InsertionScope& scope) {
  std::vector<int> nearby;
  const auto mark = [touched](int r) {
    if (touched) (*touched)[r] = true;
  };
  std::vector<int> ejections(solution->instance().numNodes(), 0);
  int num_ejections = 0;
  while (!pool->empty()) {
    const int u = pool->back();
    pool->pop_back();
    const int r = insertCheapest(u, solution);
    if (r >= 0) {
      mark(r);
      continue;
    }
    if (num_ejections == max_ejections) {
      pool->push_back(u);
      return false;
    }

    ++ejections[u];
    int changed = -1;
    if (scope.restricted()) {
      nearbyRoutes(u, scope, Positions(*solution).route_of, &nearby);
      changed = replaceWithEjection(u, max_ejected, ejections, pool, solution,
                                    &nearby);
    }
    if (changed < 0) {
      changed = replaceWithEjection(u, max_ejected, ejections, pool, solution);
    }
    if (changed < 0) {
      pool->push_back(u);  // u fits nowhere, even by replacement: give up
      return false;
    }
    ++num_ejections;
    mark(changed);
  }
  return true;
}

bool eliminateRoute(int r, const RouteEliminationParams& params,
                    Solution* solution) {
  const Solution backup = *solution;
  std::vector<Route>& routes = *solution->mutableRoutes();

  std::vector<int> pool(routes[r].nodes().begin() + 1,
                        routes[r].nodes().end() - 1);
  routes.erase(routes.begin() + r);

  if (insertWithEjections(params.max_ejections, params.max_ejected, &pool,
                          nullptr, solution)) {
    return true;
  }
  *solution = backup;
  return false;
}

RouteEliminationStats minimizeRoutes(
    const RouteEliminationParams& params,
    std::chrono::steady_clock::time_point deadline, LocalSearch* local_search,
    Solution* solution) {
  const int lower_bound = capacityLowerBound(solution->instance());

  RouteEliminationStats stats;
  bool progress = true;
  while (progress && solution->numUsedRoutes() > lower_bound &&
         std::chrono::steady_clock::now() < deadline) {
    progress = false;
    solution->removeEmptyRoutes();

    // Smallest routes first: they have the fewest customers to place.
    const std::vector<Route>& routes = solution->routes();
    std::vector<int> order(routes.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      return routes[a].numCustomers() < routes[b].numCustomers();
    });

    // Escalate: try every route with single ejections before allowing
    // larger ones, which disturb the solution more.
    for (int ejected = 1; ejected <= params.max_ejected && !progress;
         ++ejected) {
      RouteEliminationParams level = params;
      level.max_ejected = ejected;
      for (int r : order) {
        if (std::chrono::steady_clock::now() >= deadline) break;
        ++stats.attempts;
        if (eliminateRoute(r, level, solution)) {
          ++stats.removed;
          local_search->run(solution);
          progress = true;
          break;  // route indices changed; recompute the order
        }
      }
    }
  }
  return stats;
}
