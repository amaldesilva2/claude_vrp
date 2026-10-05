#include "ejection_search.h"

#include <algorithm>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "destroy_repair.h"
#include "local_search.h"
#include "segment.h"

namespace {

constexpr double kEps = 1e-6;
// Deletions in a row that may get stuck before the search gives up.
constexpr int kMaxFailures = 10;

class EjectionSearch {
 public:
  EjectionSearch(const EjectionSearchParams& params, const Instance& inst)
      : params_(params),
        inst_(inst),
        nearest_(nearestCustomers(inst, params.neighbors)),
        squeezer_(inst, squeezeParams(params)),
        rng_(params.seed),
        counters_(inst.numNodes(), 1) {}

  EjectionSearchStats run(std::chrono::steady_clock::time_point deadline,
                          Solution* solution) {
    deadline_ = deadline;
    const int lower_bound = capacityLowerBound(inst_);
    solution->removeEmptyRoutes();
    Solution best = *solution;
    int failures = 0;  // consecutive deletions that got stuck
    while (std::chrono::steady_clock::now() < deadline_ &&
           best.numUsedRoutes() > lower_bound && failures < kMaxFailures) {
      ++stats_.attempts;
      Solution candidate = best;
      const Outcome outcome = deleteRoute(&candidate);
      if (outcome == Outcome::kTimeout) break;
      if (outcome == Outcome::kStuck) {
        ++failures;
        continue;
      }
      failures = 0;
      candidate.removeEmptyRoutes();
      candidate.validate();
      best = std::move(candidate);
      ++stats_.routes_removed;
    }
    *solution = std::move(best);
    return stats_;
  }

 private:
  static LocalSearchParams squeezeParams(const EjectionSearchParams& params) {
    LocalSearchParams ls;
    ls.tw_penalty = params.squeeze_penalty;
    ls.load_penalty = params.squeeze_penalty;
    ls.adapt_period = 0;
    ls.seed = params.seed;
    return ls;
  }

  bool pastDeadline() const {
    return std::chrono::steady_clock::now() >= deadline_;
  }

  enum class Outcome { kSuccess, kTimeout, kStuck };

  // Deletes a random route and reinserts its customers. Unless the result is
  // kSuccess, `solution` is left in an unspecified state.
  Outcome deleteRoute(Solution* solution) {
    solution->removeEmptyRoutes();
    std::vector<Route>& routes = *solution->mutableRoutes();
    std::uniform_int_distribution<std::size_t> pick_route(0, routes.size() - 1);
    const std::size_t r = pick_route(rng_);
    std::vector<int> pool(routes[r].nodes().begin() + 1,
                          routes[r].nodes().end() - 1);
    routes.erase(routes.begin() + r);
    std::fill(counters_.begin(), counters_.end(), 1);

    while (!pool.empty()) {
      if (pastDeadline()) return Outcome::kTimeout;
      const int v = pool.back();
      pool.pop_back();
      if (insertRandomFeasible(v, solution)) {
        ++stats_.insertions;
        continue;
      }
      if (params_.squeeze && squeeze(v, solution)) {
        ++stats_.squeezes;
        continue;
      }
      ++counters_[v];
      // No ejection set of size <= max_ejected makes room for v.
      if (!insertWithEjections(v, solution, &pool)) return Outcome::kStuck;
      ++stats_.ejections;
      perturb(solution);
    }
    return Outcome::kSuccess;
  }

  // Inserts v at a uniformly random feasible position, if any.
  bool insertRandomFeasible(int v, Solution* solution) {
    std::vector<Route>& routes = *solution->mutableRoutes();
    int chosen_route = -1;
    int chosen_pos = -1;
    int count = 0;
    for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
      for (int pos = 0; pos <= routes[r].numCustomers(); ++pos) {
        if (!routes[r].canInsert(v, pos)) continue;
        // Reservoir sampling: each feasible position equally likely.
        if (std::uniform_int_distribution<int>(0, count++)(rng_) == 0) {
          chosen_route = r;
          chosen_pos = pos;
        }
      }
    }
    if (chosen_route < 0) return false;
    routes[chosen_route].insert(v, chosen_pos);
    return true;
  }

  // Inserts v at the position with the least excess load + time warp, then
  // repairs with penalized local search. Undone if still infeasible.
  bool squeeze(int v, Solution* solution) {
    std::vector<Route>& routes = *solution->mutableRoutes();
    const Segment single = Segment::single(inst_, v);
    int best_route = -1;
    int best_pos = -1;
    double best_penalty = std::numeric_limits<double>::infinity();
    for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
      const Route& route = routes[r];
      for (int k = 0; k + 1 < route.numNodes(); ++k) {
        const Segment seg =
            concat(inst_, concat(inst_, route.forward(k), single),
                   route.backward(k + 1));
        const double penalty =
            seg.time_warp + std::max(0, seg.load - inst_.capacity);
        if (penalty < best_penalty) {
          best_penalty = penalty;
          best_route = r;
          best_pos = k;  // between node k and k + 1
        }
      }
    }
    if (best_route < 0) return false;

    const Solution backup = *solution;
    routes[best_route].insert(v, best_pos);
    std::vector<bool> touched(routes.size(), false);
    touched[best_route] = true;
    squeezer_.run(solution, &touched);
    if (solution->isFeasible()) return true;
    *solution = backup;
    return false;
  }

  // Depth-first search over ejection sets for one route with v inserted
  // between nodes `pos` and `pos + 1`.
  struct Search {
    const Instance* inst;
    const Route* route;
    const std::vector<int>* counters;
    std::vector<int> seq;   // nodes after the start depot, v included
    std::vector<int> orig;  // original node index of seq[i], -1 for v
    int v_index = 0;        // index of v in seq
    int extra_load = 0;     // route load + demand of v
    int max_ejected = 0;
    int budget = 0;
    std::vector<int> path;  // seq indices ejected so far
    // Best found (shared across routes and positions).
    int* best_psum;
    int* ties;
    std::vector<int>* best_ejected;  // seq indices
    bool* improved;
    std::mt19937* rng;

    void record(int psum) {
      if (psum < *best_psum) {
        *best_psum = psum;
        *ties = 1;
      } else if (psum == *best_psum) {
        ++*ties;
        if (std::uniform_int_distribution<int>(1, *ties)(*rng) != 1) return;
      } else {
        return;
      }
      *best_ejected = path;
      *improved = true;
    }

    void dfs(int i, int last, double depart, int psum, int ejected_load) {
      if (--budget < 0) return;
      const int w = seq[i];
      const double arrival = depart + inst->distance(last, w);
      const bool load_ok = extra_load - ejected_load <= inst->capacity;
      if (w == 0) {  // end depot
        if (load_ok && arrival <= inst->horizon() + kEps) record(psum);
        return;
      }
      // After v the rest of the original route is unchanged: if no further
      // ejection is needed, the latest start times decide feasibility.
      if (i > v_index && load_ok &&
          arrival <= route->latestStart(orig[i]) + kEps) {
        record(psum);
      }
      const Node& node = inst->nodes[w];
      const double start = std::max(arrival, node.ready);
      if (start <= node.due + kEps) {
        dfs(i + 1, w, start + node.service, psum, ejected_load);
      }
      if (i != v_index && static_cast<int>(path.size()) < max_ejected &&
          psum + (*counters)[w] <= *best_psum) {
        path.push_back(i);
        dfs(i + 1, last, depart, psum + (*counters)[w],
            ejected_load + node.demand);
        path.pop_back();
      }
    }
  };

  // Inserts v with ejections minimizing the sum of ejected counters; the
  // ejected customers are pushed onto the pool.
  bool insertWithEjections(int v, Solution* solution, std::vector<int>* pool) {
    std::vector<Route>& routes = *solution->mutableRoutes();
    int best_psum = std::numeric_limits<int>::max();
    int ties = 0;
    int best_route = -1;
    std::vector<int> best_ejected;
    std::vector<int> best_seq;

    for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
      const Route& route = routes[r];
      const int n = route.numNodes();
      for (int pos = 0; pos + 1 < n; ++pos) {
        Search search;
        search.inst = &inst_;
        search.route = &route;
        search.counters = &counters_;
        for (int k = 1; k <= pos; ++k) {
          search.seq.push_back(route.node(k));
          search.orig.push_back(k);
        }
        search.v_index = static_cast<int>(search.seq.size());
        search.seq.push_back(v);
        search.orig.push_back(-1);
        for (int k = pos + 1; k < n; ++k) {
          search.seq.push_back(route.node(k));
          search.orig.push_back(k);
        }
        search.extra_load = route.load() + inst_.nodes[v].demand;
        search.max_ejected = params_.max_ejected;
        search.budget = params_.max_dfs_nodes;
        search.best_psum = &best_psum;
        search.ties = &ties;
        std::vector<int> ejected;
        search.best_ejected = &ejected;
        bool improved = false;
        search.improved = &improved;
        search.rng = &rng_;
        search.dfs(0, 0, inst_.nodes[0].ready, 0, 0);
        if (improved) {
          best_route = r;
          best_ejected = ejected;
          best_seq = search.seq;
        }
      }
    }
    if (best_route < 0) return false;

    std::vector<bool> drop(best_seq.size(), false);
    for (int i : best_ejected) {
      drop[i] = true;
      pool->push_back(best_seq[i]);
    }
    std::vector<int> nodes{0};
    for (std::size_t i = 0; i < best_seq.size(); ++i) {
      if (!drop[i]) nodes.push_back(best_seq[i]);
    }
    routes[best_route].assign(std::move(nodes));
    return true;
  }

  // Random feasible relocate / swap / 2-opt* moves between neighbors.
  void perturb(Solution* solution) {
    std::vector<Route>& routes = *solution->mutableRoutes();
    Positions positions(*solution);
    const auto refresh = [&](int r) {
      for (int k = 1; k + 1 < routes[r].numNodes(); ++k) {
        positions.route_of[routes[r].node(k)] = r;
        positions.index_of[routes[r].node(k)] = k;
      }
    };
    std::uniform_int_distribution<int> pick_customer(1, inst_.numCustomers());
    std::uniform_int_distribution<int> pick_move(0, 2);

    for (int t = 0; t < params_.perturbation_moves; ++t) {
      const int u = pick_customer(rng_);
      const int ra = positions.route_of[u];
      if (ra < 0 || nearest_[u].empty()) continue;
      std::uniform_int_distribution<std::size_t> pick_near(
          0, nearest_[u].size() - 1);
      const int v = nearest_[u][pick_near(rng_)];
      const int rb = positions.route_of[v];
      if (rb < 0 || rb == ra) continue;
      Route& a = routes[ra];
      Route& b = routes[rb];
      const int ia = positions.index_of[u];
      const int ib = positions.index_of[v];

      switch (pick_move(rng_)) {
        case 0: {  // relocate u right after v
          if (!b.canInsert(u, ib)) break;
          b.insert(u, ib);
          a.remove(ia - 1);
          positions.route_of[u] = rb;
          refresh(ra);
          refresh(rb);
          break;
        }
        case 1: {  // swap u and v
          const Segment sa = concat(
              inst_,
              concat(inst_, a.forward(ia - 1), Segment::single(inst_, v)),
              a.backward(ia + 1));
          const Segment sb = concat(
              inst_,
              concat(inst_, b.forward(ib - 1), Segment::single(inst_, u)),
              b.backward(ib + 1));
          if (sa.time_warp > kEps || sb.time_warp > kEps ||
              sa.load > inst_.capacity || sb.load > inst_.capacity) {
            break;
          }
          std::vector<int> na = a.nodes();
          std::vector<int> nb = b.nodes();
          na[ia] = v;
          nb[ib] = u;
          a.assign(std::move(na));
          b.assign(std::move(nb));
          refresh(ra);
          refresh(rb);
          break;
        }
        default: {  // 2-opt*: exchange the tails after u and after v
          const Segment sa = concat(inst_, a.forward(ia), b.backward(ib + 1));
          const Segment sb = concat(inst_, b.forward(ib), a.backward(ia + 1));
          if (sa.time_warp > kEps || sb.time_warp > kEps ||
              sa.load > inst_.capacity || sb.load > inst_.capacity) {
            break;
          }
          std::vector<int> na(a.nodes().begin(), a.nodes().begin() + ia + 1);
          na.insert(na.end(), b.nodes().begin() + ib + 1, b.nodes().end());
          std::vector<int> nb(b.nodes().begin(), b.nodes().begin() + ib + 1);
          nb.insert(nb.end(), a.nodes().begin() + ia + 1, a.nodes().end());
          a.assign(std::move(na));
          b.assign(std::move(nb));
          refresh(ra);
          refresh(rb);
          break;
        }
      }
    }
  }

  const EjectionSearchParams& params_;
  const Instance& inst_;
  std::vector<std::vector<int>> nearest_;
  LocalSearch squeezer_;
  std::mt19937 rng_;
  std::vector<int> counters_;  // p[v]: how often v needed ejections
  std::chrono::steady_clock::time_point deadline_;
  EjectionSearchStats stats_;
};

}  // namespace

EjectionSearchStats minimizeRoutesEjectionSearch(
    const EjectionSearchParams& params,
    std::chrono::steady_clock::time_point deadline, Solution* solution) {
  return EjectionSearch(params, solution->instance()).run(deadline, solution);
}
