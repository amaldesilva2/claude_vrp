#include "destroy_repair.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "construction.h"

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
// Stands in for a missing h-th best insertion in regret repair, so that
// customers with fewer feasible routes are inserted first.
constexpr double kMissingRegret = 1e9;

// Index in [0, size) drawn as floor(y^randomness * size), y uniform in [0, 1).
int randomizedRank(int size, double randomness, std::mt19937* rng) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  const int rank = static_cast<int>(std::pow(unit(*rng), randomness) * size);
  return std::min(rank, size - 1);
}

std::vector<int> routedCustomers(const Solution& solution) {
  std::vector<int> customers;
  for (const Route& route : solution.routes()) {
    for (int k = 1; k + 1 < route.numNodes(); ++k) {
      customers.push_back(route.node(k));
    }
  }
  return customers;
}

// Cheapest feasible position of `customer` in `route`; false if none.
bool bestPosition(const Route& route, int customer, double* cost, int* pos,
                  const InsertionScope* blinks = nullptr) {
  *cost = kInf;
  *pos = -1;
  const bool blink =
      blinks != nullptr && blinks->blink_rate > 0 && blinks->rng != nullptr;
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (int p = 0; p <= route.numCustomers(); ++p) {
    if (blink && unit(*blinks->rng) < blinks->blink_rate) continue;
    if (!route.canInsert(customer, p)) continue;
    const double c = route.insertionCost(customer, p);
    if (c < *cost) {
      *cost = c;
      *pos = p;
    }
  }
  return *pos >= 0;
}

}  // namespace

Positions::Positions(const Solution& solution)
    : route_of(solution.instance().numNodes(), -1),
      index_of(solution.instance().numNodes(), -1) {
  const std::vector<Route>& routes = solution.routes();
  for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
    for (int k = 1; k + 1 < routes[r].numNodes(); ++k) {
      route_of[routes[r].node(k)] = r;
      index_of[routes[r].node(k)] = k;
    }
  }
}

void nearbyRoutes(int customer, const InsertionScope& scope,
                  const std::vector<int>& route_of, std::vector<int>* routes) {
  routes->clear();
  const std::vector<int>& near = (*scope.nearest)[customer];
  const int k = std::min<int>(scope.neighbors, static_cast<int>(near.size()));
  for (int i = 0; i < k; ++i) {
    const int r = route_of[near[i]];
    if (r >= 0 &&
        std::find(routes->begin(), routes->end(), r) == routes->end()) {
      routes->push_back(r);
    }
  }
}

// ---------------------------------------------------------------------------
// Destroy operators
// ---------------------------------------------------------------------------

std::vector<int> radialRemoval(const Instance& inst,
                               const std::vector<std::vector<int>>& nearest,
                               int count, std::mt19937* rng) {
  std::uniform_int_distribution<int> pick_customer(1, inst.numCustomers());
  const int center = pick_customer(*rng);
  const int size = std::min(count, inst.numCustomers());
  std::vector<int> removed{center};
  for (int i = 0; i + 1 < size && i < static_cast<int>(nearest[center].size());
       ++i) {
    removed.push_back(nearest[center][i]);
  }
  return removed;
}

std::vector<int> stringRemoval(const Solution& solution,
                               const std::vector<std::vector<int>>& nearest,
                               double avg_removed, int max_string,
                               std::mt19937* rng, double split_probability) {
  const Instance& inst = solution.instance();
  const std::vector<Route>& routes = solution.routes();
  const Positions positions(solution);

  double avg_route_size = 0.0;
  int used = 0;
  for (const Route& route : routes) {
    if (route.empty()) continue;
    avg_route_size += route.numCustomers();
    ++used;
  }
  avg_route_size /= std::max(1, used);

  const double max_len =
      std::min(static_cast<double>(max_string), avg_route_size);
  const double max_strings =
      std::max(1.0, 4.0 * avg_removed / (1.0 + max_len) - 1.0);
  std::uniform_real_distribution<double> pick_strings(1.0, max_strings + 1.0);
  const int num_strings = static_cast<int>(pick_strings(*rng));

  std::uniform_int_distribution<int> pick_customer(1, inst.numCustomers());
  const int seed = pick_customer(*rng);
  std::vector<int> candidates{seed};
  candidates.insert(candidates.end(), nearest[seed].begin(),
                    nearest[seed].end());

  std::vector<int> removed;
  std::vector<bool> route_ruined(routes.size(), false);
  int strings = 0;
  for (int c : candidates) {
    if (strings == num_strings) break;
    const int r = positions.route_of[c];
    if (r < 0 || route_ruined[r]) continue;
    route_ruined[r] = true;
    ++strings;

    // String of length `len` containing c, at a random offset.
    const Route& route = routes[r];
    const int route_len = route.numCustomers();
    const double len_cap = std::min(static_cast<double>(route_len), max_len);
    std::uniform_real_distribution<double> pick_len(1.0, len_cap + 1.0);
    const int len = std::min(route_len, static_cast<int>(pick_len(*rng)));
    const int k = positions.index_of[c];

    // Split string: keep m consecutive customers inside a window of len + m.
    int keep = 0;
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    if (unit(*rng) < split_probability && len < route_len) {
      keep = 1;
      while (len + keep < route_len && unit(*rng) > 0.01) ++keep;
    }
    const int window = len + keep;
    const int lo = std::max(1, k - window + 1);
    const int hi = std::min(k, route_len - window + 1);
    std::uniform_int_distribution<int> pick_start(lo, hi);
    const int first = pick_start(*rng);
    std::uniform_int_distribution<int> pick_kept(first, first + len);
    const int kept_first = keep > 0 ? pick_kept(*rng) : first + window;
    for (int i = first; i < first + window; ++i) {
      if (i >= kept_first && i < kept_first + keep) continue;
      removed.push_back(route.node(i));
    }
  }
  return removed;
}

std::vector<int> randomRemoval(const Solution& solution, int count,
                               std::mt19937* rng) {
  std::vector<int> customers = routedCustomers(solution);
  std::shuffle(customers.begin(), customers.end(), *rng);
  customers.resize(std::min<std::size_t>(customers.size(), count));
  return customers;
}

std::vector<int> worstRemoval(const Solution& solution, int count,
                              double randomness, std::mt19937* rng) {
  const Instance& inst = solution.instance();
  std::vector<std::pair<double, int>> gains;  // (saving, customer)
  for (const Route& route : solution.routes()) {
    for (int k = 1; k + 1 < route.numNodes(); ++k) {
      const int prev = route.node(k - 1);
      const int u = route.node(k);
      const int next = route.node(k + 1);
      gains.emplace_back(inst.distance(prev, u) + inst.distance(u, next) -
                             inst.distance(prev, next),
                         u);
    }
  }
  std::sort(gains.begin(), gains.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

  // Savings are computed once on the original solution.
  std::vector<int> removed;
  while (static_cast<int>(removed.size()) < count && !gains.empty()) {
    const int rank =
        randomizedRank(static_cast<int>(gains.size()), randomness, rng);
    removed.push_back(gains[rank].second);
    gains.erase(gains.begin() + rank);
  }
  return removed;
}

std::vector<int> relatedRemoval(const Solution& solution,
                                const std::vector<std::vector<int>>& nearest,
                                int count, double randomness,
                                std::mt19937* rng) {
  const Instance& inst = solution.instance();
  const std::vector<Route>& routes = solution.routes();
  const Positions positions(solution);
  std::vector<int> routed = routedCustomers(solution);
  if (routed.empty()) return {};

  // Service start time of each routed customer.
  std::vector<double> start(inst.numNodes(), 0.0);
  for (int c : routed) {
    start[c] =
        routes[positions.route_of[c]].startTime(positions.index_of[c] - 1);
  }
  // Relatedness (lower = more related), with weights from Ropke & Pisinger.
  // Distances and times share units, so both are scaled by the horizon.
  const double horizon = inst.horizon();
  const auto relatedness = [&](int i, int j) {
    return 9.0 * inst.distance(i, j) / horizon +
           3.0 * std::fabs(start[i] - start[j]) / horizon +
           2.0 * std::abs(inst.nodes[i].demand - inst.nodes[j].demand) /
               static_cast<double>(inst.capacity);
  };

  std::uniform_int_distribution<std::size_t> pick_routed(0, routed.size() - 1);
  std::vector<bool> is_removed(inst.numNodes(), false);
  std::vector<int> removed{routed[pick_routed(*rng)]};
  is_removed[removed[0]] = true;
  const int target = std::min<int>(count, static_cast<int>(routed.size()));

  std::vector<std::pair<double, int>> candidates;
  while (static_cast<int>(removed.size()) < target) {
    std::uniform_int_distribution<std::size_t> pick_removed(0,
                                                            removed.size() - 1);
    const int r = removed[pick_removed(*rng)];
    candidates.clear();
    for (int c : nearest[r]) {
      if (positions.route_of[c] >= 0 && !is_removed[c]) {
        candidates.emplace_back(relatedness(r, c), c);
      }
    }
    int chosen;
    if (candidates.empty()) {
      // No unremoved neighbor: fall back to any unremoved routed customer.
      do {
        chosen = routed[pick_routed(*rng)];
      } while (is_removed[chosen]);
    } else {
      std::sort(candidates.begin(), candidates.end());
      chosen = candidates[randomizedRank(static_cast<int>(candidates.size()),
                                         randomness, rng)]
                   .second;
    }
    is_removed[chosen] = true;
    removed.push_back(chosen);
  }
  return removed;
}

std::vector<int> routeRemoval(const Solution& solution, int count,
                              std::mt19937* rng) {
  const std::vector<Route>& routes = solution.routes();
  std::vector<int> candidates;
  std::vector<double> weights;
  for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
    if (routes[r].empty()) continue;
    candidates.push_back(r);
    weights.push_back(1.0 / routes[r].numCustomers());
  }

  std::vector<int> removed;
  while (!candidates.empty() &&
         (removed.empty() || static_cast<int>(removed.size()) < count)) {
    std::discrete_distribution<int> pick(weights.begin(), weights.end());
    const int i = pick(*rng);
    const Route& route = routes[candidates[i]];
    removed.insert(removed.end(), route.nodes().begin() + 1,
                   route.nodes().end() - 1);
    candidates.erase(candidates.begin() + i);
    weights.erase(weights.begin() + i);
  }
  return removed;
}

void removeCustomers(const std::vector<int>& removed,
                     std::vector<bool>* touched, Solution* solution) {
  std::vector<bool> is_removed(solution->instance().numNodes(), false);
  for (int c : removed) is_removed[c] = true;
  std::vector<Route>& routes = *solution->mutableRoutes();
  for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
    const std::vector<int>& old_nodes = routes[r].nodes();
    std::vector<int> kept;
    kept.reserve(old_nodes.size());
    for (int node : old_nodes) {
      if (!is_removed[node]) kept.push_back(node);
    }
    if (kept.size() != old_nodes.size()) {
      routes[r].assign(std::move(kept));
      (*touched)[r] = true;
    }
  }
}

// ---------------------------------------------------------------------------
// Repair operators
// ---------------------------------------------------------------------------

std::vector<int> greedyRepair(std::vector<int> customers, InsertionOrder order,
                              std::mt19937* rng, std::vector<bool>* touched,
                              Solution* solution, const InsertionScope& scope) {
  const Instance& inst = solution->instance();
  const std::vector<Node>& nodes = inst.nodes;
  std::shuffle(customers.begin(), customers.end(), *rng);
  switch (order) {
    case InsertionOrder::kRandom:
      break;
    case InsertionOrder::kFarthestFirst:
      std::stable_sort(customers.begin(), customers.end(), [&](int a, int b) {
        return inst.distance(0, a) > inst.distance(0, b);
      });
      break;
    case InsertionOrder::kDemandFirst:
      std::stable_sort(customers.begin(), customers.end(), [&](int a, int b) {
        return nodes[a].demand > nodes[b].demand;
      });
      break;
    case InsertionOrder::kDeadlineFirst:
      std::stable_sort(customers.begin(), customers.end(), [&](int a, int b) {
        return nodes[a].due < nodes[b].due;
      });
      break;
  }

  std::vector<Route>& routes = *solution->mutableRoutes();
  std::vector<int> route_of;
  if (scope.restricted()) route_of = Positions(*solution).route_of;
  std::vector<int> nearby;

  std::vector<int> leftover;
  for (int c : customers) {
    int r = -1;
    if (scope.restricted()) {
      nearbyRoutes(c, scope, route_of, &nearby);
      double best_cost = kInf;
      int best_pos = -1;
      for (int candidate : nearby) {
        double cost;
        int pos;
        if (bestPosition(routes[candidate], c, &cost, &pos, &scope) &&
            cost < best_cost) {
          best_cost = cost;
          best_pos = pos;
          r = candidate;
        }
      }
      if (r >= 0) routes[r].insert(c, best_pos);
    }
    if (r < 0) r = insertCheapest(c, solution);  // all routes
    if (r < 0) {
      leftover.push_back(c);
    } else {
      (*touched)[r] = true;
      if (scope.restricted()) route_of[c] = r;
    }
  }
  return leftover;
}

std::vector<int> regretRepair(std::vector<int> customers, int k,
                              std::vector<bool>* touched, Solution* solution,
                              const InsertionScope& scope) {
  std::vector<Route>& routes = *solution->mutableRoutes();
  std::vector<int> active;  // non-empty routes; they stay non-empty
  for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
    if (!routes[r].empty()) active.push_back(r);
  }
  const int num_routes = static_cast<int>(active.size());

  std::vector<int> column_of(routes.size(), -1);
  for (int j = 0; j < num_routes; ++j) column_of[active[j]] = j;

  // cost[i][j], pos[i][j]: best insertion of customers[i] into active[j];
  // considered[i][j]: whether that route is evaluated for customers[i] (all
  // routes, or with a restricted scope its nearby routes, or all routes if
  // none of those fits).
  std::vector<std::vector<double>> cost(customers.size(),
                                        std::vector<double>(num_routes, kInf));
  std::vector<std::vector<int>> pos(customers.size(),
                                    std::vector<int>(num_routes, -1));
  std::vector<std::vector<char>> considered(
      customers.size(), std::vector<char>(num_routes, !scope.restricted()));
  std::vector<int> route_of;
  if (scope.restricted()) route_of = Positions(*solution).route_of;
  std::vector<int> nearby;
  for (std::size_t i = 0; i < customers.size(); ++i) {
    bool any = false;
    if (scope.restricted()) {
      nearbyRoutes(customers[i], scope, route_of, &nearby);
      for (int r : nearby) {
        const int j = column_of[r];
        if (j < 0) continue;
        considered[i][j] = 1;
        any |= bestPosition(routes[r], customers[i], &cost[i][j], &pos[i][j],
                            &scope);
      }
      if (!any) std::fill(considered[i].begin(), considered[i].end(), 1);
    }
    if (!any) {
      for (int j = 0; j < num_routes; ++j) {
        bestPosition(routes[active[j]], customers[i], &cost[i][j], &pos[i][j]);
      }
    }
  }

  std::vector<double> smallest(k);
  while (!customers.empty()) {
    int best_i = -1;
    int best_j = -1;
    double best_regret = -kInf;
    double best_cost = kInf;
    for (std::size_t i = 0; i < customers.size(); ++i) {
      // The k smallest costs over routes (kInf where fewer exist).
      std::fill(smallest.begin(), smallest.end(), kInf);
      int argmin = -1;
      for (int j = 0; j < num_routes; ++j) {
        const double c = cost[i][j];
        if (c >= smallest[k - 1]) continue;
        if (c < smallest[0]) argmin = j;
        int h = k - 1;
        while (h > 0 && smallest[h - 1] > c) {
          smallest[h] = smallest[h - 1];
          --h;
        }
        smallest[h] = c;
      }
      if (argmin < 0) continue;  // fits nowhere

      double regret = 0.0;
      for (int h = 1; h < k; ++h) {
        regret +=
            (smallest[h] == kInf ? kMissingRegret : smallest[h]) - smallest[0];
      }
      if (regret > best_regret ||
          (regret == best_regret && smallest[0] < best_cost)) {
        best_regret = regret;
        best_cost = smallest[0];
        best_i = static_cast<int>(i);
        best_j = argmin;
      }
    }
    if (best_i < 0) break;  // nobody fits anywhere

    const int r = active[best_j];
    routes[r].insert(customers[best_i], pos[best_i][best_j]);
    (*touched)[r] = true;

    // Drop the inserted customer and refresh the changed route's column.
    // (Swap-and-pop; guard against self-move when best_i is the last row.)
    const std::size_t last = customers.size() - 1;
    if (static_cast<std::size_t>(best_i) != last) {
      customers[best_i] = customers[last];
      cost[best_i] = std::move(cost[last]);
      pos[best_i] = std::move(pos[last]);
      considered[best_i] = std::move(considered[last]);
    }
    customers.pop_back();
    cost.pop_back();
    pos.pop_back();
    considered.pop_back();
    for (std::size_t i = 0; i < customers.size(); ++i) {
      if (!considered[i][best_j]) continue;
      bestPosition(routes[r], customers[i], &cost[i][best_j], &pos[i][best_j],
                   &scope);
    }
  }

  // With a restricted scope, a customer whose nearby routes filled up may
  // still fit elsewhere.
  std::vector<int> leftover;
  for (int c : customers) {
    const int r = scope.restricted() ? insertCheapest(c, solution) : -1;
    if (r >= 0) {
      (*touched)[r] = true;
    } else {
      leftover.push_back(c);
    }
  }
  return leftover;
}
