#include "iterated_local_search.h"

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

#include "construction.h"
#include "destroy_repair.h"
#include "route_elimination.h"

namespace {

// Neighbor lists used by the ruin operators.
constexpr int kRuinNeighbors = 100;

class IteratedLocalSearch {
 public:
  IteratedLocalSearch(const IlsParams& params, LocalSearch* local_search,
                      const Instance& inst)
      : params_(params),
        inst_(inst),
        local_search_(local_search),
        nearest_(nearestCustomers(inst, kRuinNeighbors)),
        lower_bound_(capacityLowerBound(inst)),
        rng_(params.seed),
        penalty_weight_(params.unplaced_penalty) {}

  IlsStats run(std::chrono::steady_clock::time_point deadline,
               Solution* solution) {
    using Clock = std::chrono::steady_clock;
    const Clock::time_point start = Clock::now();
    const double total =
        std::chrono::duration<double>(deadline - start).count();

    // `best` is always complete (empty pool).
    State best{*solution, {}, solution->numUsedRoutes()};
    State current = best;
    int incomplete_iterations = 0;
    bool elimination_due = false;

    while (Clock::now() < deadline) {
      ++stats_.iterations;

      // A current state that stays incomplete for too long (e.g. after a
      // soft elimination) is replaced by the best (complete) one.
      if (current.pool.empty()) {
        incomplete_iterations = 0;
        penalty_weight_ = params_.unplaced_penalty;
      } else if (++incomplete_iterations > params_.max_incomplete_iterations) {
        ++stats_.incomplete_restarts;
        incomplete_iterations = 0;
        current = best;
      }

      // Every elimination_period iterations an elimination becomes due; it
      // runs at the next iteration where the current state is complete.
      if (params_.elimination_period > 0 &&
          stats_.iterations % params_.elimination_period == 0) {
        elimination_due = true;
      }
      if (elimination_due && current.pool.empty()) {
        elimination_due = false;
        if (current.solution.numUsedRoutes() > lower_bound_ &&
            tryElimination(&current)) {
          updateBest(current, &best);
          continue;
        }
      }

      // Staying incomplete gets more expensive every iteration.
      if (!current.pool.empty()) penalty_weight_ *= params_.penalty_growth;

      const double elapsed =
          std::chrono::duration<double>(Clock::now() - start).count();
      const double threshold =
          params_.initial_threshold * std::max(0.0, 1.0 - elapsed / total);

      State candidate = current;
      if (!perturb(&candidate)) continue;

      // Fewer vehicles only counts once every customer is placed; a route
      // emptied into the pool is judged by penalized cost like any other
      // candidate.
      const int vehicles = candidate.solution.numUsedRoutes();
      if (vehicles > candidate.fleet) continue;
      const double cost = penalizedCost(candidate);
      const bool fewer_vehicles =
          candidate.pool.empty() && vehicles < candidate.fleet;
      const bool accept = fewer_vehicles || cost < penalizedCost(current) ||
                          cost < penalizedCost(best) * (1.0 + threshold);
      if (!accept) continue;
      if (fewer_vehicles) candidate.fleet = vehicles;

      ++stats_.accepted;
      updateBest(candidate, &best);
      current = std::move(candidate);
    }

    *solution = std::move(best.solution);
    return stats_;
  }

 private:
  // A solution plus the customers that are currently not in any route.
  struct State {
    Solution solution;
    std::vector<int> pool;
    int fleet;  // vehicles the state may use; below it only with empty pool
  };

  // Distance plus, per unplaced customer c, unplaced_penalty * 2 * d(0, c),
  // roughly the cost of serving c on a route of its own.
  double penalizedCost(const State& state) const {
    double cost = state.solution.totalDistance();
    for (int c : state.pool) {
      cost += penalty_weight_ * 2.0 * inst_.distance(0, c);
    }
    return cost;
  }

  // Validation is O(n) and new bests are rare, so every one is checked: a
  // bookkeeping bug fails here instead of corrupting the result.
  void updateBest(const State& state, State* best) {
    if (state.pool.empty() && state.solution.isBetterThan(best->solution)) {
      state.solution.validate();
      *best = state;
      ++stats_.improvements;
    }
  }

  // Ruin, recreate, local search, and placement of leftover customers.
  // Returns false if local search could not make the candidate feasible
  // (soft constraints only); the candidate must then be discarded.
  bool perturb(State* state) {
    Solution* solution = &state->solution;
    std::vector<bool> touched(solution->routes().size(), false);
    std::vector<int> removed;
    std::bernoulli_distribution use_strings(params_.string_ruin_probability);
    if (use_strings(rng_)) {
      ++stats_.string_ruins;
      removed = stringRuin(*solution);
    } else {
      removed = radialRuin();
    }
    removeCustomers(removed, &touched, solution);

    // Customers already waiting in the pool get another chance as well. The
    // radial ruin may have picked some of them; keep each customer once.
    std::vector<bool> listed(inst_.numNodes(), false);
    for (int c : removed) listed[c] = true;
    for (int c : state->pool) {
      if (!listed[c]) removed.push_back(c);
    }
    std::vector<int> pool = recreate(&removed, &touched, solution);
    if (!local_search_->runAndRepair(solution, &touched)) return false;
    if (pool.empty()) {
      state->pool.clear();
      return true;
    }

    // Leftovers: retry after local search, with ejections.
    const std::vector<int> leftover = pool;
    std::vector<bool> touched_again(solution->routes().size(), false);
    insertWithEjections(params_.max_pool_ejections, params_.max_pool_ejected,
                        &pool, &touched_again, solution);
    std::vector<bool> still_unplaced(inst_.numNodes(), false);
    for (int c : pool) still_unplaced[c] = true;
    for (int c : leftover) {
      if (!still_unplaced[c]) ++stats_.pool_rescues;
    }
    if (!params_.penalize_unplaced) {
      for (int c : pool) {
        solution->addRoute().insert(c, 0);
        touched_again.push_back(true);
        ++stats_.new_routes;
      }
      pool.clear();
    }
    state->pool = std::move(pool);
    return local_search_->runAndRepair(solution, &touched_again);
  }

  std::vector<int> radialRuin() {
    std::uniform_int_distribution<int> pick_size(params_.min_ruin,
                                                 params_.max_ruin);
    return radialRemoval(inst_, nearest_, pick_size(rng_), &rng_);
  }

  std::vector<int> stringRuin(const Solution& solution) {
    return stringRemoval(solution, nearest_, params_.avg_removed,
                         params_.max_string, &rng_);
  }

  // Cheapest feasible insertion in a randomly chosen order. Returns the
  // customers that fit in no existing route.
  std::vector<int> recreate(std::vector<int>* removed,
                            std::vector<bool>* touched, Solution* solution) {
    std::uniform_int_distribution<int> pick_order(0, 3);
    const InsertionOrder order = static_cast<InsertionOrder>(pick_order(rng_));
    return greedyRepair(std::move(*removed), order, &rng_, touched, solution);
  }

  // Tries to remove one of the smallest routes of the (complete) state,
  // rotating among the five smallest so that repeated calls try different
  // routes. First with ejections; if that fails and unplaced customers are
  // allowed, the route's customers are moved to the pool (soft elimination)
  // and later iterations try to place them.
  bool tryElimination(State* state) {
    ++stats_.elimination_attempts;
    Solution* solution = &state->solution;
    solution->removeEmptyRoutes();
    const std::vector<Route>& routes = solution->routes();
    std::vector<int> order(routes.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      return routes[a].numCustomers() < routes[b].numCustomers();
    });
    const int r = order[stats_.elimination_attempts %
                        std::min<int>(5, static_cast<int>(order.size()))];

    // Single ejections first; larger ones disturb the solution more.
    RouteEliminationParams params;
    params.max_ejections = params_.max_elimination_ejections;
    for (int ejected = 1; ejected <= params_.max_elimination_ejected;
         ++ejected) {
      params.max_ejected = ejected;
      if (eliminateRoute(r, params, solution)) {
        ++stats_.routes_eliminated;
        const Solution eliminated = *solution;
        if (!local_search_->runAndRepair(solution)) *solution = eliminated;
        state->fleet = solution->numUsedRoutes();
        return true;
      }
    }
    if (!params_.penalize_unplaced || !params_.soft_elimination) return false;

    std::vector<Route>& mutable_routes = *solution->mutableRoutes();
    state->pool.assign(mutable_routes[r].nodes().begin() + 1,
                       mutable_routes[r].nodes().end() - 1);
    mutable_routes.erase(mutable_routes.begin() + r);
    ++stats_.soft_eliminations;
    state->fleet = solution->numUsedRoutes();
    return true;
  }

  const IlsParams& params_;
  const Instance& inst_;
  LocalSearch* local_search_;
  std::vector<std::vector<int>> nearest_;
  int lower_bound_;
  std::mt19937 rng_;
  double penalty_weight_ = 0.0;  // current weight of unplaced customers
  IlsStats stats_;
};

}  // namespace

IlsStats iteratedLocalSearch(const IlsParams& params,
                             std::chrono::steady_clock::time_point deadline,
                             LocalSearch* local_search, Solution* solution) {
  return IteratedLocalSearch(params, local_search, solution->instance())
      .run(deadline, solution);
}
