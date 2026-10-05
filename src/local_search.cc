#include "local_search.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "segment.h"

namespace {

constexpr double kEps = 1e-6;

// Nodes route.node(from..to) inclusive, possibly reversed. Empty if from > to.
struct Piece {
  int route = 0;
  int from = 0;
  int to = -1;
  bool reversed = false;

  bool empty() const { return from > to; }
};

// A candidate route described as a concatenation of pieces of current routes.
struct Plan {
  std::array<Piece, 5> pieces;
  int size = 0;

  Plan& add(int route, int from, int to, bool reversed = false) {
    pieces[size++] = Piece{route, from, to, reversed};
    return *this;
  }
};

// Identifies a move variant for statistics; formatted only when applied.
struct MoveName {
  const char* kind;
  bool intra = false;
  int a = 0;
  bool rev_a = false;
  int b = 0;
  bool rev_b = false;

  std::string str() const {
    std::string s = intra ? "intra " : "inter ";
    s += kind;
    if (a > 0) {
      s += "(" + std::to_string(a) + (rev_a ? "r" : "");
      if (b > 0) s += "," + std::to_string(b) + (rev_b ? "r" : "");
      s += ")";
    }
    return s;
  }
};

// State of one local search run on one solution.
class LocalSearchRun {
 public:
  LocalSearchRun(const LocalSearchParams& params,
                 const std::vector<std::vector<int>>& neighbors,
                 std::mt19937* rng, Solution* solution)
      : params_(params),
        inst_(solution->instance()),
        solution_(solution),
        routes_(*solution->mutableRoutes()),
        neighbors_(neighbors),
        route_of_(inst_.numNodes(), -1),
        index_of_(inst_.numNodes(), -1),
        last_tested_(inst_.numNodes(), 0),
        route_modified_(routes_.size(), 1),
        soft_(params.tw_penalty > 0 || params.load_penalty > 0),
        rng_(*rng) {}

  LocalSearchStats run(const std::vector<bool>* touched_routes) {
    for (int r = 0; r < static_cast<int>(routes_.size()); ++r) {
      updatePositions(r);
      // Untouched routes count as unmodified since time 0, so pairs between
      // them are skipped until one of them changes.
      if (touched_routes && !(*touched_routes)[r]) route_modified_[r] = 0;
    }

    std::vector<int> order(inst_.numCustomers());
    std::iota(order.begin(), order.end(), 1);

    bool improved = true;
    while (improved) {
      improved = false;
      ++stats_.rounds;
      std::shuffle(order.begin(), order.end(), rng_);
      for (int u : order) {
        if (route_of_[u] < 0) continue;  // not routed (e.g. ILS pool)
        const int last_tested = last_tested_[u];
        last_tested_[u] = move_count_;
        for (int v : neighbors_[u]) {
          if (route_of_[v] < 0) continue;
          const int changed = std::max(route_modified_[route_of_[u]],
                                       route_modified_[route_of_[v]]);
          if (changed <= last_tested) continue;
          if (tryAllMoves(u, v)) improved = true;
        }
      }
    }

    solution_->removeEmptyRoutes();
    return stats_;
  }

 private:
  // ---------------------------------------------------------------------
  // Setup and bookkeeping
  // ---------------------------------------------------------------------

  void updatePositions(int r) {
    const Route& route = routes_[r];
    for (int k = 1; k + 1 < route.numNodes(); ++k) {
      route_of_[route.node(k)] = r;
      index_of_[route.node(k)] = k;
    }
  }

  int depotEnd(int r) const { return routes_[r].numNodes() - 1; }

  // ---------------------------------------------------------------------
  // Piece and plan evaluation
  // ---------------------------------------------------------------------

  int firstNode(const Piece& p) const {
    return routes_[p.route].node(p.reversed ? p.to : p.from);
  }
  int lastNode(const Piece& p) const {
    return routes_[p.route].node(p.reversed ? p.from : p.to);
  }

  // Distances are symmetric, so a reversed piece has the same length.
  double pieceDistance(const Piece& p) const {
    const Route& route = routes_[p.route];
    return route.forward(p.to).distance - route.forward(p.from).distance;
  }

  int pieceLoad(const Piece& p) const {
    const Route& route = routes_[p.route];
    return route.forward(p.to).load -
           (p.from > 0 ? route.forward(p.from - 1).load : 0);
  }

  // O(1) for route prefixes/suffixes; otherwise O(piece length).
  Segment pieceSegment(const Piece& p) const {
    const Route& route = routes_[p.route];
    if (!p.reversed) {
      if (p.from == 0) return route.forward(p.to);
      if (p.to == route.numNodes() - 1) return route.backward(p.from);
      Segment s = Segment::single(inst_, route.node(p.from));
      for (int k = p.from + 1; k <= p.to; ++k) {
        s = concat(inst_, s, Segment::single(inst_, route.node(k)));
      }
      return s;
    }
    Segment s = Segment::single(inst_, route.node(p.to));
    for (int k = p.to - 1; k >= p.from; --k) {
      s = concat(inst_, s, Segment::single(inst_, route.node(k)));
    }
    return s;
  }

  double planDistance(const Plan& plan) const {
    double total = 0.0;
    int prev = -1;
    for (int i = 0; i < plan.size; ++i) {
      const Piece& p = plan.pieces[i];
      if (p.empty()) continue;
      total += pieceDistance(p);
      if (prev >= 0) total += inst_.distance(prev, firstNode(p));
      prev = lastNode(p);
    }
    return total;
  }

  int planLoad(const Plan& plan) const {
    int load = 0;
    for (int i = 0; i < plan.size; ++i) {
      if (!plan.pieces[i].empty()) load += pieceLoad(plan.pieces[i]);
    }
    return load;
  }

  // Summary of the whole candidate route.
  Segment planSegment(const Plan& plan) const {
    Segment total;
    bool started = false;
    for (int i = 0; i < plan.size; ++i) {
      const Piece& p = plan.pieces[i];
      if (p.empty()) continue;
      const Segment s = pieceSegment(p);
      total = started ? concat(inst_, total, s) : s;
      started = true;
    }
    return total;
  }

  // Penalty part of the soft-constraint cost of a route or segment.
  double penalty(double time_warp, int load) const {
    return params_.tw_penalty * time_warp +
           params_.load_penalty * std::max(0, load - inst_.capacity);
  }
  double routePenalty(int r) const {
    return penalty(routes_[r].timeWarp(), routes_[r].load());
  }

  bool planTimeFeasible(const Plan& plan) const {
    Segment total;
    bool started = false;
    for (int i = 0; i < plan.size; ++i) {
      const Piece& p = plan.pieces[i];
      if (p.empty()) continue;
      const Segment s = pieceSegment(p);
      total = started ? concat(inst_, total, s) : s;
      started = true;
      // Time warp never decreases as pieces are appended.
      if (total.time_warp > kEps) return false;
    }
    return true;
  }

  std::vector<int> planNodes(const Plan& plan) const {
    std::vector<int> nodes;
    for (int i = 0; i < plan.size; ++i) {
      const Piece& p = plan.pieces[i];
      if (p.empty()) continue;
      const Route& route = routes_[p.route];
      if (p.reversed) {
        for (int k = p.to; k >= p.from; --k) nodes.push_back(route.node(k));
      } else {
        for (int k = p.from; k <= p.to; ++k) nodes.push_back(route.node(k));
      }
    }
    return nodes;
  }

  // ---------------------------------------------------------------------
  // Generic move: replace route ra by plan_a and, for inter-route moves,
  // route rb by plan_b. Applied if feasible and distance decreases.
  // ---------------------------------------------------------------------

  bool tryMove(const MoveName& name, int ra, const Plan& plan_a, int rb = -1,
               const Plan* plan_b = nullptr) {
    // 1. Distance (cheapest; rejects most candidates).
    const double new_a = planDistance(plan_a);
    const double new_b = plan_b ? planDistance(*plan_b) : 0.0;
    double delta = new_a - routes_[ra].distance();
    if (plan_b) delta += new_b - routes_[rb].distance();

    double expected_a = new_a;
    double expected_b = new_b;
    if (!soft_) {
      if (delta > -kEps) return false;

      // 2. Capacity (an intra-route move keeps the load unchanged).
      if (plan_b && (planLoad(plan_a) > inst_.capacity ||
                     planLoad(*plan_b) > inst_.capacity)) {
        return false;
      }

      // 3. Time windows.
      if (!planTimeFeasible(plan_a)) return false;
      if (plan_b && !planTimeFeasible(*plan_b)) return false;
    } else {
      // Soft constraints: new penalties are >= 0, so the distance change
      // minus the old penalties bounds the cost change from below.
      const double old_penalty =
          routePenalty(ra) + (plan_b ? routePenalty(rb) : 0.0);
      if (delta - old_penalty > -kEps) return false;
      const Segment seg_a = planSegment(plan_a);
      expected_a = seg_a.distance + penalty(seg_a.time_warp, seg_a.load);
      double cost_delta =
          expected_a - routes_[ra].distance() - routePenalty(ra);
      if (plan_b) {
        const Segment seg_b = planSegment(*plan_b);
        expected_b = seg_b.distance + penalty(seg_b.time_warp, seg_b.load);
        cost_delta += expected_b - routes_[rb].distance() - routePenalty(rb);
      }
      if (cost_delta > -kEps) return false;
    }

    // Apply. Build both node lists before touching either route.
    std::vector<int> nodes_a = planNodes(plan_a);
    std::vector<int> nodes_b = plan_b ? planNodes(*plan_b) : std::vector<int>();
    ++move_count_;
    routes_[ra].assign(std::move(nodes_a));
    route_modified_[ra] = move_count_;
    updatePositions(ra);
    if (plan_b) {
      routes_[rb].assign(std::move(nodes_b));
      route_modified_[rb] = move_count_;
      updatePositions(rb);
    }

    if (params_.check_moves) {
      checkApplied(name, ra, expected_a);
      if (plan_b) checkApplied(name, rb, expected_b);
    }
    ++stats_.moves_applied;
    ++stats_.by_move[name.str()];
    return true;
  }

  // `expected` is the distance (hard) or penalized cost (soft) of route r.
  void checkApplied(const MoveName& name, int r, double expected) const {
    const Route& route = routes_[r];
    const double actual = route.distance() + (soft_ ? routePenalty(r) : 0.0);
    if (std::fabs(actual - expected) > 1e-6 ||
        (!soft_ && !route.isFeasible())) {
      throw std::logic_error("move " + name.str() +
                             " evaluated incorrectly on route " +
                             std::to_string(r));
    }
  }

  // ---------------------------------------------------------------------
  // Move generators for customer u and neighbor v
  // ---------------------------------------------------------------------

  bool tryAllMoves(int u, int v) {
    const bool intra = route_of_[u] == route_of_[v];

    const int max_run = params_.max_run;
    const bool reversed = params_.reversed_runs;
    for (int a = 1; a <= max_run; ++a) {
      for (bool rev : {false, true}) {
        if (rev && (a == 1 || !reversed)) continue;
        for (bool after : {true, false}) {
          if (intra ? relocateIntra(u, v, a, rev, after)
                    : relocateInter(u, v, a, rev, after)) {
            return true;
          }
        }
      }
    }

    for (int a = 1; a <= max_run; ++a) {
      for (bool rev_a : {false, true}) {
        if (rev_a && (a == 1 || !reversed)) continue;
        for (int b = 1; b <= max_run; ++b) {
          for (bool rev_b : {false, true}) {
            if (rev_b && (b == 1 || !reversed)) continue;
            if (intra ? swapIntra(u, v, a, rev_a, b, rev_b)
                      : swapInter(u, v, a, rev_a, b, rev_b)) {
              return true;
            }
          }
        }
      }
    }

    if (intra) return twoOpt(u, v);
    return twoOptStar(u, v) || twoOptStarReversed(u, v);
  }

  // Moves the run of `a` customers starting at u to just after (or before) v
  // in another route.
  bool relocateInter(int u, int v, int a, bool rev, bool after) {
    const int ra = route_of_[u], pu = index_of_[u], ea = depotEnd(ra);
    const int rb = route_of_[v], pv = index_of_[v], eb = depotEnd(rb);
    if (pu + a - 1 > ea - 1) return false;
    const int k = after ? pv : pv - 1;  // insert between node k and k + 1

    Plan plan_a, plan_b;
    plan_a.add(ra, 0, pu - 1).add(ra, pu + a, ea);
    plan_b.add(rb, 0, k).add(ra, pu, pu + a - 1, rev).add(rb, k + 1, eb);
    return tryMove(
        {after ? "relocate-after" : "relocate-before", false, a, rev}, ra,
        plan_a, rb, &plan_b);
  }

  bool relocateIntra(int u, int v, int a, bool rev, bool after) {
    const int r = route_of_[u], pu = index_of_[u], e = depotEnd(r);
    const int pv = index_of_[v];
    if (pu + a - 1 > e - 1) return false;
    const int k = after ? pv : pv - 1;  // insert between node k and k + 1
    const int last = pu + a - 1;

    Plan plan;
    if (k < pu - 1) {
      plan.add(r, 0, k)
          .add(r, pu, last, rev)
          .add(r, k + 1, pu - 1)
          .add(r, last + 1, e);
    } else if (k > last) {
      plan.add(r, 0, pu - 1)
          .add(r, last + 1, k)
          .add(r, pu, last, rev)
          .add(r, k + 1, e);
    } else {
      return false;  // v inside the run, or the run would not move
    }
    return tryMove({after ? "relocate-after" : "relocate-before", true, a, rev},
                   r, plan);
  }

  // Exchanges the run of `a` customers starting at u with the run of `b`
  // customers starting at v (different routes).
  bool swapInter(int u, int v, int a, bool rev_a, int b, bool rev_b) {
    const int ra = route_of_[u], pu = index_of_[u], ea = depotEnd(ra);
    const int rb = route_of_[v], pv = index_of_[v], eb = depotEnd(rb);
    if (pu + a - 1 > ea - 1 || pv + b - 1 > eb - 1) return false;

    Plan plan_a, plan_b;
    plan_a.add(ra, 0, pu - 1)
        .add(rb, pv, pv + b - 1, rev_b)
        .add(ra, pu + a, ea);
    plan_b.add(rb, 0, pv - 1)
        .add(ra, pu, pu + a - 1, rev_a)
        .add(rb, pv + b, eb);
    return tryMove({"swap", false, a, rev_a, b, rev_b}, ra, plan_a, rb,
                   &plan_b);
  }

  bool swapIntra(int u, int v, int a, bool rev_a, int b, bool rev_b) {
    const int r = route_of_[u], e = depotEnd(r);
    const int pu = index_of_[u], pv = index_of_[v];
    if (pu + a - 1 > e - 1 || pv + b - 1 > e - 1) return false;

    // Order the two runs along the route; they must not overlap.
    int p1 = pu, n1 = a, p2 = pv, n2 = b;
    bool rev1 = rev_a, rev2 = rev_b;
    if (pv < pu) {
      std::swap(p1, p2);
      std::swap(n1, n2);
      std::swap(rev1, rev2);
    }
    if (p1 + n1 - 1 >= p2) return false;

    Plan plan;
    plan.add(r, 0, p1 - 1)
        .add(r, p2, p2 + n2 - 1, rev2)
        .add(r, p1 + n1, p2 - 1)
        .add(r, p1, p1 + n1 - 1, rev1)
        .add(r, p2 + n2, e);
    return tryMove({"swap", true, a, rev_a, b, rev_b}, r, plan);
  }

  // A: 0..u | x..0, B: 0..v | y..0  ->  A': 0..u y..0, B': 0..v x..0
  bool twoOptStar(int u, int v) {
    const int ra = route_of_[u], pu = index_of_[u], ea = depotEnd(ra);
    const int rb = route_of_[v], pv = index_of_[v], eb = depotEnd(rb);

    Plan plan_a, plan_b;
    plan_a.add(ra, 0, pu).add(rb, pv + 1, eb);
    plan_b.add(rb, 0, pv).add(ra, pu + 1, ea);
    return tryMove({"2-opt*"}, ra, plan_a, rb, &plan_b);
  }

  // A: 0..u | x..0, B: 0..v | y..0
  //   ->  A': 0..u v..(B start) 0,  B': 0 (A end)..x y..0
  bool twoOptStarReversed(int u, int v) {
    const int ra = route_of_[u], pu = index_of_[u], ea = depotEnd(ra);
    const int rb = route_of_[v], pv = index_of_[v], eb = depotEnd(rb);

    Plan plan_a, plan_b;
    plan_a.add(ra, 0, pu).add(rb, 1, pv, true).add(ra, ea, ea);
    plan_b.add(rb, 0, 0).add(ra, pu + 1, ea - 1, true).add(rb, pv + 1, eb);
    return tryMove({"2-opt*rev"}, ra, plan_a, rb, &plan_b);
  }

  // Reverses the path between u and v so that they become adjacent:
  // 0..i | i+1..j | j+1..0  ->  0..i j..i+1 j+1..0
  bool twoOpt(int u, int v) {
    const int r = route_of_[u], e = depotEnd(r);
    const int i = std::min(index_of_[u], index_of_[v]);
    const int j = std::max(index_of_[u], index_of_[v]);
    if (j - i < 2) return false;

    Plan plan;
    plan.add(r, 0, i).add(r, i + 1, j, true).add(r, j + 1, e);
    return tryMove({"2-opt", true}, r, plan);
  }

  const LocalSearchParams& params_;
  const Instance& inst_;
  Solution* solution_;
  std::vector<Route>& routes_;
  const std::vector<std::vector<int>>& neighbors_;  // K nearest customers
  std::vector<int> route_of_;                       // customer -> route
  std::vector<int> index_of_;                       // customer -> node index
  std::vector<int> last_tested_;     // customer -> move_count_ when examined
  std::vector<int> route_modified_;  // route -> move_count_ when changed
  int move_count_ = 1;
  bool soft_;  // soft constraints (penalized time warp and excess load)
  std::mt19937& rng_;
  LocalSearchStats stats_;
};

}  // namespace

LocalSearch::LocalSearch(const Instance& inst, const LocalSearchParams& params)
    : inst_(inst),
      params_(params),
      neighbors_(nearestCustomers(inst, params.num_neighbors)),
      rng_(params.seed) {}

LocalSearchStats LocalSearch::run(Solution* solution,
                                  const std::vector<bool>* touched_routes) {
  if (&solution->instance() != &inst_) {
    throw std::invalid_argument("solution belongs to a different instance");
  }
  if (touched_routes && touched_routes->size() != solution->routes().size()) {
    throw std::invalid_argument("touched_routes size != number of routes");
  }
  return LocalSearchRun(params_, neighbors_, &rng_, solution)
      .run(touched_routes);
}

void LocalSearch::setPenalties(double tw_penalty, double load_penalty) {
  params_.tw_penalty = tw_penalty;
  params_.load_penalty = load_penalty;
}

bool LocalSearch::runAndRepair(Solution* solution,
                               const std::vector<bool>* touched_routes) {
  run(solution, touched_routes);
  if (!softConstraints()) return true;

  // Penalty adaptation (Vidal et al.): steer towards target_feasible.
  bool tw_ok = true;
  bool load_ok = true;
  for (const Route& route : solution->routes()) {
    if (route.timeWarp() > 1e-6) tw_ok = false;
    if (route.excessLoad() > 0) load_ok = false;
  }
  if (params_.adapt_period > 0) {
    ++adapt_calls_;
    adapt_tw_feasible_ += tw_ok;
    adapt_load_feasible_ += load_ok;
    if (adapt_calls_ == params_.adapt_period) {
      const auto adapt = [&](int feasible, double* penalty) {
        const double fraction = static_cast<double>(feasible) / adapt_calls_;
        if (fraction < params_.target_feasible - 0.05) {
          *penalty = std::min(1e5, *penalty * 1.2);
        } else if (fraction > params_.target_feasible + 0.05) {
          *penalty = std::max(0.01, *penalty * 0.85);
        }
      };
      if (params_.tw_penalty > 0)
        adapt(adapt_tw_feasible_, &params_.tw_penalty);
      if (params_.load_penalty > 0) {
        adapt(adapt_load_feasible_, &params_.load_penalty);
      }
      adapt_calls_ = adapt_tw_feasible_ = adapt_load_feasible_ = 0;
    }
  }
  if (tw_ok && load_ok) return true;

  // Repair with stronger penalties.
  const double tw = params_.tw_penalty;
  const double load = params_.load_penalty;
  bool feasible = false;
  for (double factor : {10.0, 100.0}) {
    setPenalties(tw * factor, load * factor);
    run(solution);
    if (solution->isFeasible()) {
      feasible = true;
      break;
    }
  }
  setPenalties(tw, load);
  return feasible;
}
