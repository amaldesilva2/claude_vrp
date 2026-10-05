#include "adaptive_lns.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <utility>

#include "destroy_repair.h"
#include "route_elimination.h"

namespace {

// Neighbor lists used by the radial, string and related removal operators.
constexpr int kRemovalNeighbors = 100;

enum class Destroy { kRandom, kWorst, kRelated, kRoute, kRadial, kString };
enum class Repair { kGreedy, kGreedyDeadline, kRegret2, kRegret3 };

constexpr Destroy kDestroyOps[] = {Destroy::kRandom,  Destroy::kWorst,
                                   Destroy::kRelated, Destroy::kRoute,
                                   Destroy::kRadial,  Destroy::kString};
constexpr const char* kDestroyNames[] = {"random", "worst",  "related",
                                         "route",  "radial", "string"};
constexpr Repair kRepairOps[] = {Repair::kGreedy, Repair::kGreedyDeadline,
                                 Repair::kRegret2, Repair::kRegret3};
constexpr const char* kRepairNames[] = {"greedy", "greedy-deadline", "regret-2",
                                        "regret-3"};

// Roulette-wheel operator selection with adaptive weights.
class OperatorWeights {
 public:
  OperatorWeights(const char* const* names, int count) : ops_(count) {
    for (int i = 0; i < count; ++i) ops_[i].stats.name = names[i];
  }

  int pick(std::mt19937* rng) const {
    std::vector<double> weights;
    for (const Op& op : ops_) weights.push_back(op.stats.weight);
    std::discrete_distribution<int> pick(weights.begin(), weights.end());
    return pick(*rng);
  }

  void reward(int i, double score, bool new_best) {
    ops_[i].segment_score += score;
    ++ops_[i].segment_uses;
    ++ops_[i].stats.uses;
    if (new_best) ++ops_[i].stats.new_best;
  }

  // End of a segment: blend the average score into each used weight.
  void update(const AlnsParams& params) {
    for (Op& op : ops_) {
      if (op.segment_uses > 0) {
        const double average = op.segment_score / op.segment_uses;
        op.stats.weight = std::max(params.min_weight,
                                   (1.0 - params.reaction) * op.stats.weight +
                                       params.reaction * average);
      }
      op.segment_score = 0.0;
      op.segment_uses = 0;
    }
  }

  std::vector<AlnsOperatorStats> stats() const {
    std::vector<AlnsOperatorStats> result;
    for (const Op& op : ops_) result.push_back(op.stats);
    return result;
  }

 private:
  struct Op {
    AlnsOperatorStats stats;
    double segment_score = 0.0;
    int segment_uses = 0;
  };
  std::vector<Op> ops_;
};

class AdaptiveLns {
 public:
  AdaptiveLns(const AlnsParams& params, LocalSearch* local_search,
              const Instance& inst)
      : params_(params),
        inst_(inst),
        local_search_(local_search),
        nearest_(nearestCustomers(inst, kRemovalNeighbors)),
        lower_bound_(capacityLowerBound(inst)),
        destroy_(kDestroyNames, std::size(kDestroyOps)),
        repair_(kRepairNames, std::size(kRepairOps)),
        rng_(params.seed) {}

  AlnsStats run(std::chrono::steady_clock::time_point deadline,
                Solution* solution) {
    using Clock = std::chrono::steady_clock;
    const Clock::time_point start = Clock::now();
    const double total =
        std::chrono::duration<double>(deadline - start).count();

    Solution best = *solution;
    Solution current = *solution;
    // Accept start_worse * distance worse with probability 0.5 at the start.
    const double start_temperature =
        params_.start_worse * best.totalDistance() / std::log(2.0);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    const bool exchanging =
        params_.exchange != nullptr && params_.exchange_seconds > 0;
    const auto exchange_interval = std::chrono::duration_cast<Clock::duration>(
        std::chrono::duration<double>(params_.exchange_seconds));
    Clock::time_point next_exchange = start + exchange_interval;

    while (Clock::now() < deadline) {
      ++stats_.iterations;
      if (exchanging && stats_.iterations % 256 == 0 &&
          Clock::now() >= next_exchange) {
        next_exchange += exchange_interval;
        Solution global = params_.exchange->offer(best);
        if (global.isBetterThan(best)) {
          best = global;
          current = std::move(global);
          ++stats_.adopted;
        }
      }
      if (stats_.iterations % params_.segment_length == 0) {
        destroy_.update(params_);
        repair_.update(params_);
      }

      if (params_.elimination_period > 0 &&
          stats_.iterations % params_.elimination_period == 0 &&
          current.numUsedRoutes() > lower_bound_ && tryElimination(&current)) {
        if (current.isBetterThan(best)) {
          current.validate();
          best = current;
          ++stats_.improvements;
        }
        continue;
      }

      if (params_.ls_period > 0 && stats_.iterations % params_.ls_period == 0) {
        Solution improved = current;
        if (local_search_->runAndRepair(&improved)) {
          current = std::move(improved);
          if (current.isBetterThan(best)) {
            current.validate();
            best = current;
            ++stats_.improvements;
          }
        }
        continue;
      }

      const double progress =
          std::chrono::duration<double>(Clock::now() - start).count() / total;
      const double temperature =
          start_temperature *
          std::pow(params_.end_temperature_ratio, std::min(1.0, progress));

      const int d = destroy_.pick(&rng_);
      const int r = repair_.pick(&rng_);
      Solution candidate = current;
      const bool feasible = perturb(kDestroyOps[d], kRepairOps[r], &candidate);

      // Vehicles first: fewer is always accepted, more never.
      const int vehicles = candidate.numUsedRoutes();
      const int current_vehicles = current.numUsedRoutes();
      if (!feasible || vehicles > current_vehicles) {
        destroy_.reward(d, 0.0, false);
        repair_.reward(r, 0.0, false);
        continue;
      }
      const double delta = candidate.totalDistance() - current.totalDistance();
      const bool better = vehicles < current_vehicles || delta < -1e-9;
      const bool accept = better || unit(rng_) < std::exp(-delta / temperature);

      double score = 0.0;
      bool new_best = false;
      if (accept) {
        ++stats_.accepted;
        if (candidate.isBetterThan(best)) {
          candidate.validate();  // O(n); new bests are rare
          best = candidate;
          ++stats_.improvements;
          new_best = true;
          score = params_.score_new_best;
        } else {
          score = better ? params_.score_better : params_.score_accepted;
        }
        current = std::move(candidate);
      }
      destroy_.reward(d, score, new_best);
      repair_.reward(r, score, new_best);
    }

    stats_.destroy = destroy_.stats();
    stats_.repair = repair_.stats();
    *solution = std::move(best);
    return stats_;
  }

 private:
  // Returns false if local search could not make the candidate feasible
  // (soft constraints only); the candidate must then be discarded.
  bool perturb(Destroy destroy, Repair repair, Solution* solution) {
    std::uniform_int_distribution<int> pick_count(params_.min_remove,
                                                  params_.max_remove);
    const int count = pick_count(rng_);
    std::vector<int> removed;
    switch (destroy) {
      case Destroy::kRandom:
        removed = randomRemoval(*solution, count, &rng_);
        break;
      case Destroy::kWorst:
        removed =
            worstRemoval(*solution, count, params_.removal_randomness, &rng_);
        break;
      case Destroy::kRelated:
        removed = relatedRemoval(*solution, nearest_, count,
                                 params_.removal_randomness, &rng_);
        break;
      case Destroy::kRoute:
        removed = routeRemoval(*solution, count, &rng_);
        break;
      case Destroy::kRadial:
        removed = radialRemoval(inst_, nearest_, count, &rng_);
        break;
      case Destroy::kString:
        removed =
            stringRemoval(*solution, nearest_, params_.avg_removed,
                          params_.max_string, &rng_, params_.split_probability);
        break;
    }

    std::vector<bool> touched(solution->routes().size(), false);
    removeCustomers(removed, &touched, solution);

    InsertionScope scope;
    scope.nearest = &nearest_;
    scope.neighbors = params_.insertion_neighbors;
    scope.blink_rate = params_.blink_rate;
    scope.rng = &rng_;
    std::vector<int> leftover;
    switch (repair) {
      case Repair::kGreedy:
        leftover = greedyRepair(std::move(removed), InsertionOrder::kRandom,
                                &rng_, &touched, solution, scope);
        break;
      case Repair::kGreedyDeadline:
        leftover =
            greedyRepair(std::move(removed), InsertionOrder::kDeadlineFirst,
                         &rng_, &touched, solution, scope);
        break;
      case Repair::kRegret2:
        leftover =
            regretRepair(std::move(removed), 2, &touched, solution, scope);
        break;
      case Repair::kRegret3:
        leftover =
            regretRepair(std::move(removed), 3, &touched, solution, scope);
        break;
    }

    if (params_.local_search &&
        !local_search_->runAndRepair(solution, &touched)) {
      return false;
    }
    if (leftover.empty()) return true;

    // Leftovers: retry after local search with ejections, then new routes.
    std::vector<bool> touched_again(solution->routes().size(), false);
    insertWithEjections(params_.max_pool_ejections, 1, &leftover,
                        &touched_again, solution, scope);
    for (int c : leftover) {
      solution->addRoute().insert(c, 0);
      touched_again.push_back(true);
    }
    return !params_.local_search ||
           local_search_->runAndRepair(solution, &touched_again);
  }

  // Tries to remove one of the five smallest routes (rotating), with single
  // ejections first and then up to max_elimination_ejected.
  bool tryElimination(Solution* solution) {
    ++stats_.elimination_attempts;
    solution->removeEmptyRoutes();
    const std::vector<Route>& routes = solution->routes();
    std::vector<int> order(routes.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      return routes[a].numCustomers() < routes[b].numCustomers();
    });
    const int r = order[stats_.elimination_attempts %
                        std::min<int>(5, static_cast<int>(order.size()))];

    RouteEliminationParams params;
    params.max_ejections = params_.max_elimination_ejections;
    for (int ejected = 1; ejected <= params_.max_elimination_ejected;
         ++ejected) {
      params.max_ejected = ejected;
      if (eliminateRoute(r, params, solution)) {
        ++stats_.routes_eliminated;
        const Solution eliminated = *solution;
        if (params_.local_search && !local_search_->runAndRepair(solution)) {
          *solution = eliminated;
        }
        return true;
      }
    }
    return false;
  }

  const AlnsParams& params_;
  const Instance& inst_;
  LocalSearch* local_search_;
  std::vector<std::vector<int>> nearest_;
  int lower_bound_;
  OperatorWeights destroy_;
  OperatorWeights repair_;
  std::mt19937 rng_;
  AlnsStats stats_;
};

}  // namespace

AlnsStats adaptiveLns(const AlnsParams& params,
                      std::chrono::steady_clock::time_point deadline,
                      LocalSearch* local_search, Solution* solution) {
  return AdaptiveLns(params, local_search, solution->instance())
      .run(deadline, solution);
}
