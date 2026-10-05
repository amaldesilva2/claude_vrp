#include "fleet_minimization.h"

#include <algorithm>
#include <random>
#include <utility>
#include <vector>

#include "destroy_repair.h"

namespace {

// Neighbor lists for string removal and insertion.
constexpr int kNeighbors = 100;

// Removes the route with the fewest customers (ties: first) and returns its
// customers.
std::vector<int> removeSmallestRoute(Solution* solution) {
  solution->removeEmptyRoutes();
  std::vector<Route>& routes = *solution->mutableRoutes();
  const auto smallest = std::min_element(
      routes.begin(), routes.end(), [](const Route& a, const Route& b) {
        return a.numCustomers() < b.numCustomers();
      });
  std::vector<int> customers(smallest->nodes().begin() + 1,
                             smallest->nodes().end() - 1);
  routes.erase(smallest);
  return customers;
}

}  // namespace

FleetMinimizationStats minimizeFleet(
    const FleetMinimizationParams& params,
    std::chrono::steady_clock::time_point deadline, Solution* solution) {
  const Instance& inst = solution->instance();
  const int lower_bound = capacityLowerBound(inst);
  const std::vector<std::vector<int>> nearest =
      nearestCustomers(inst, kNeighbors);
  InsertionScope scope;
  scope.nearest = &nearest;
  scope.neighbors = params.insertion_neighbors;
  std::mt19937 rng(params.seed);
  std::uniform_int_distribution<int> pick_order(0, 3);

  FleetMinimizationStats stats;
  Solution best = *solution;  // complete, fewest vehicles so far
  if (best.numUsedRoutes() <= lower_bound) return stats;

  Solution current = best;
  std::vector<int> absent = removeSmallestRoute(&current);
  std::vector<long> absences(inst.numNodes(), 0);
  const auto sum_absences = [&](const std::vector<int>& customers) {
    long sum = 0;
    for (int c : customers) sum += absences[c];
    return sum;
  };

  std::size_t fewest_absent = absent.size();
  int stagnation = 0;
  int longest_streak = 0;  // within the current route removal
  while (std::chrono::steady_clock::now() < deadline) {
    ++stats.iterations;
    ++stagnation;
    if (params.max_stagnation > 0 && stagnation > params.max_stagnation) {
      stats.stagnated = true;
      break;
    }

    Solution candidate = current;
    std::vector<int> removed = stringRemoval(
        candidate, nearest, params.avg_removed, params.max_string, &rng);
    std::vector<bool> touched(candidate.routes().size(), false);
    removeCustomers(removed, &touched, &candidate);
    removed.insert(removed.end(), absent.begin(), absent.end());
    const InsertionOrder order = static_cast<InsertionOrder>(pick_order(rng));
    std::vector<int> candidate_absent = greedyRepair(
        std::move(removed), order, &rng, &touched, &candidate, scope);

    if (candidate_absent.size() < absent.size() ||
        sum_absences(candidate_absent) < sum_absences(absent)) {
      current = std::move(candidate);
      absent = std::move(candidate_absent);
    }
    for (int c : absent) ++absences[c];
    if (absent.size() < fewest_absent) {
      fewest_absent = absent.size();
      longest_streak = std::max(longest_streak, stagnation);
      stagnation = 0;
    }

    if (absent.empty()) {
      // One vehicle fewer.
      current.removeEmptyRoutes();
      current.validate();
      best = current;
      ++stats.routes_removed;
      stats.longest_successful_stagnation =
          std::max(stats.longest_successful_stagnation, longest_streak);
      longest_streak = 0;
      if (best.numUsedRoutes() <= lower_bound) break;
      absent = removeSmallestRoute(&current);
      fewest_absent = absent.size();
      stagnation = 0;
    }
  }

  *solution = std::move(best);
  return stats;
}
