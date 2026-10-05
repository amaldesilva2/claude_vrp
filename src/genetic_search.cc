#include "genetic_search.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include "destroy_repair.h"
#include "local_search.h"
#include "segment.h"

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

struct Individual {
  Solution solution;
  double distance = 0.0;
  double time_warp = 0.0;
  int excess_load = 0;
  bool feasible = false;
  std::vector<int> succ;  // successor of each customer (0 = depot)
  std::vector<int> pred;  // predecessor of each customer (0 = depot)
  // (broken-pairs distance, other individual), ascending.
  std::vector<std::pair<double, Individual*>> proximity;
  double fitness = 0.0;  // biased fitness, lower is better

  explicit Individual(Solution s) : solution(std::move(s)) {
    const Instance& inst = solution.instance();
    succ.assign(inst.numNodes(), 0);
    pred.assign(inst.numNodes(), 0);
    for (const Route& route : solution.routes()) {
      distance += route.distance();
      time_warp += route.timeWarp();
      excess_load += route.excessLoad();
      for (int k = 1; k + 1 < route.numNodes(); ++k) {
        pred[route.node(k)] = route.node(k - 1);
        succ[route.node(k)] = route.node(k + 1);
      }
    }
    feasible = time_warp < 1e-6 && excess_load == 0;
  }

  double cost(double tw_penalty, double load_penalty) const {
    return distance + tw_penalty * time_warp + load_penalty * excess_load;
  }
};

// Broken-pairs distance (Vidal et al.): share of customers whose neighbors
// differ between the two solutions.
double brokenPairsDistance(const Individual& a, const Individual& b) {
  const int n = static_cast<int>(a.succ.size());
  int differences = 0;
  for (int c = 1; c < n; ++c) {
    if (a.succ[c] != b.succ[c] && a.succ[c] != b.pred[c]) ++differences;
    if (a.pred[c] == 0 && b.pred[c] != 0 && b.succ[c] != 0) ++differences;
  }
  return static_cast<double>(differences) / (n - 1);
}

using Population = std::vector<std::unique_ptr<Individual>>;

class GeneticSearch {
 public:
  GeneticSearch(const GeneticSearchParams& params, const Instance& inst)
      : params_(params),
        inst_(inst),
        local_search_(inst, localSearchParams(params)),
        rng_(params.seed) {
    // Initial penalties: 1 per unit of time warp, and distance per unit of
    // demand for excess load.
    double demand = 0.0;
    double distance = 0.0;
    for (int c = 1; c < inst.numNodes(); ++c) {
      demand += inst.nodes[c].demand;
      distance += inst.distance(0, c);
    }
    tw_penalty_ = 1.0;
    load_penalty_ = std::max(0.1, distance / std::max(1.0, demand));
  }

  GeneticSearchStats run(std::chrono::steady_clock::time_point deadline,
                         Solution* solution) {
    using Clock = std::chrono::steady_clock;
    best_ = std::make_unique<Solution>(*solution);

    // Initial population: the input and perturbed copies of it.
    addIndividual(*solution);
    const int removal = std::max(
        1, static_cast<int>(params_.initial_removal * inst_.numCustomers()));
    for (int i = 1; i < params_.initial_population && Clock::now() < deadline;
         ++i) {
      Solution s = *solution;
      std::vector<int> removed = randomRemoval(s, removal, &rng_);
      std::vector<bool> touched(s.routes().size(), false);
      removeCustomers(removed, &touched, &s);
      std::shuffle(removed.begin(), removed.end(), rng_);
      for (int c : removed) insertPenalized(c, &s);
      educateAndAdd(std::move(s));
    }

    while (Clock::now() < deadline) {
      const Individual* a = tournament();
      const Individual* b = tournament();
      for (int tries = 0; tries < 5 && b == a; ++tries) b = tournament();
      educateAndAdd(srex(a->solution, b->solution));
    }

    stats_.final_tw_penalty = tw_penalty_;
    stats_.final_load_penalty = load_penalty_;
    *solution = std::move(*best_);
    return stats_;
  }

 private:
  static LocalSearchParams localSearchParams(
      const GeneticSearchParams& params) {
    LocalSearchParams ls;
    ls.num_neighbors = params.num_neighbors;
    ls.tw_penalty = 1.0;  // set before every use
    ls.load_penalty = 1.0;
    ls.adapt_period = 0;  // adapted here instead
    ls.seed = params.seed;
    return ls;
  }

  double routeCost(const Route& route) const {
    return route.distance() + tw_penalty_ * route.timeWarp() +
           load_penalty_ * route.excessLoad();
  }

  double solutionCost(const Solution& s) const {
    double cost = 0.0;
    for (const Route& route : s.routes()) cost += routeCost(route);
    return cost;
  }

  // Inserts c where distance + penalties increase least (may be infeasible).
  void insertPenalized(int c, Solution* s) {
    std::vector<Route>& routes = *s->mutableRoutes();
    const Segment single = Segment::single(inst_, c);
    int best_route = -1;
    int best_pos = -1;
    double best_delta = kInf;
    for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
      const Route& route = routes[r];
      const double old_cost = routeCost(route);
      for (int k = 0; k + 1 < route.numNodes(); ++k) {
        const Segment seg =
            concat(inst_, concat(inst_, route.forward(k), single),
                   route.backward(k + 1));
        const double delta =
            seg.distance + tw_penalty_ * seg.time_warp +
            load_penalty_ * std::max(0, seg.load - inst_.capacity) - old_cost;
        if (delta < best_delta) {
          best_delta = delta;
          best_route = r;
          best_pos = k;
        }
      }
    }
    routes[best_route].insert(c, best_pos);
  }

  // Penalized local search, penalty adaptation, insertion into the
  // population and (with repair_probability) a repaired copy.
  void educateAndAdd(Solution s) {
    local_search_.setPenalties(tw_penalty_, load_penalty_);
    local_search_.run(&s);
    ++stats_.generations;

    bool tw_ok = true;
    bool load_ok = true;
    for (const Route& route : s.routes()) {
      if (route.timeWarp() > 1e-6) tw_ok = false;
      if (route.excessLoad() > 0) load_ok = false;
    }
    adaptPenalties(tw_ok, load_ok);
    if (tw_ok && load_ok) ++stats_.feasible_offspring;

    const bool feasible = tw_ok && load_ok;
    std::bernoulli_distribution repair(params_.repair_probability);
    Solution repaired = feasible ? Solution(inst_) : s;
    addIndividual(std::move(s));
    if (!feasible && repair(rng_)) {
      for (double factor : {10.0, 100.0}) {
        local_search_.setPenalties(tw_penalty_ * factor,
                                   load_penalty_ * factor);
        local_search_.run(&repaired);
        if (repaired.isFeasible()) break;
      }
      if (repaired.isFeasible()) {
        ++stats_.repaired;
        addIndividual(std::move(repaired));
      }
    }
  }

  void adaptPenalties(bool tw_ok, bool load_ok) {
    ++adapt_count_;
    adapt_tw_ok_ += tw_ok;
    adapt_load_ok_ += load_ok;
    if (adapt_count_ < params_.adapt_period) return;
    const auto adapt = [&](int ok, double* penalty) {
      const double fraction = static_cast<double>(ok) / adapt_count_;
      if (fraction < params_.target_feasible - 0.05) {
        *penalty = std::min(1e5, *penalty * 1.2);
      } else if (fraction > params_.target_feasible + 0.05) {
        *penalty = std::max(0.1, *penalty * 0.85);
      }
    };
    adapt(adapt_tw_ok_, &tw_penalty_);
    adapt(adapt_load_ok_, &load_penalty_);
    adapt_count_ = adapt_tw_ok_ = adapt_load_ok_ = 0;
  }

  void addIndividual(Solution s) {
    s.removeEmptyRoutes();
    auto individual = std::make_unique<Individual>(std::move(s));
    Individual* added = individual.get();
    if (added->feasible && added->solution.isBetterThan(*best_)) {
      added->solution.validate();  // O(n); new bests are rare
      *best_ = added->solution;
      ++stats_.improvements;
    }
    Population& population = added->feasible ? feasible_ : infeasible_;
    for (const auto& other : population) {
      const double d = brokenPairsDistance(*added, *other);
      insertSorted(&added->proximity, {d, other.get()});
      insertSorted(&other->proximity, {d, added});
    }
    population.push_back(std::move(individual));
    if (static_cast<int>(population.size()) >=
        params_.min_population + params_.generation_size) {
      selectSurvivors(&population);
    }
  }

  static void insertSorted(std::vector<std::pair<double, Individual*>>* list,
                           std::pair<double, Individual*> entry) {
    list->insert(std::upper_bound(list->begin(), list->end(), entry,
                                  [](const auto& x, const auto& y) {
                                    return x.first < y.first;
                                  }),
                 entry);
  }

  double averageClosest(const Individual& individual) const {
    const int k = std::min<int>(params_.num_close,
                                static_cast<int>(individual.proximity.size()));
    if (k == 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < k; ++i) sum += individual.proximity[i].first;
    return sum / k;
  }

  // Biased fitness = cost rank + (1 - elite/size) * diversity rank.
  void updateFitness(Population* population) {
    const int size = static_cast<int>(population->size());
    if (size == 0) return;
    if (size == 1) {
      (*population)[0]->fitness = 0.0;
      return;
    }
    std::vector<int> by_cost(size);
    std::iota(by_cost.begin(), by_cost.end(), 0);
    std::sort(by_cost.begin(), by_cost.end(), [&](int a, int b) {
      return (*population)[a]->cost(tw_penalty_, load_penalty_) <
             (*population)[b]->cost(tw_penalty_, load_penalty_);
    });
    std::vector<int> by_diversity(size);
    std::iota(by_diversity.begin(), by_diversity.end(), 0);
    std::vector<double> diversity(size);
    for (int i = 0; i < size; ++i) {
      diversity[i] = averageClosest(*(*population)[i]);
    }
    std::sort(by_diversity.begin(), by_diversity.end(),
              [&](int a, int b) { return diversity[a] > diversity[b]; });

    std::vector<double> cost_rank(size);
    std::vector<double> diversity_rank(size);
    for (int i = 0; i < size; ++i) {
      cost_rank[by_cost[i]] = static_cast<double>(i) / (size - 1);
      diversity_rank[by_diversity[i]] = static_cast<double>(i) / (size - 1);
    }
    const double diversity_weight =
        size <= params_.num_elite
            ? 0.0
            : 1.0 - static_cast<double>(params_.num_elite) / size;
    for (int i = 0; i < size; ++i) {
      (*population)[i]->fitness =
          cost_rank[i] + diversity_weight * diversity_rank[i];
    }
  }

  // Removes the worst individuals (clones first) down to min_population.
  void selectSurvivors(Population* population) {
    while (static_cast<int>(population->size()) > params_.min_population) {
      updateFitness(population);
      int worst = -1;
      bool worst_is_clone = false;
      for (int i = 0; i < static_cast<int>(population->size()); ++i) {
        const Individual& ind = *(*population)[i];
        const bool clone =
            !ind.proximity.empty() && ind.proximity.front().first < 1e-9;
        if (worst < 0 || (clone && !worst_is_clone) ||
            (clone == worst_is_clone &&
             ind.fitness > (*population)[worst]->fitness)) {
          worst = i;
          worst_is_clone = clone;
        }
      }
      Individual* removed = (*population)[worst].get();
      for (const auto& other : *population) {
        auto& list = other->proximity;
        list.erase(
            std::remove_if(list.begin(), list.end(),
                           [&](const auto& e) { return e.second == removed; }),
            list.end());
      }
      population->erase(population->begin() + worst);
    }
  }

  // Binary tournament on biased fitness over both sub-populations.
  const Individual* tournament() {
    updateFitness(&feasible_);
    updateFitness(&infeasible_);
    const int total = static_cast<int>(feasible_.size() + infeasible_.size());
    std::uniform_int_distribution<int> pick(0, total - 1);
    const auto get = [&](int i) -> const Individual* {
      return i < static_cast<int>(feasible_.size())
                 ? feasible_[i].get()
                 : infeasible_[i - feasible_.size()].get();
    };
    const Individual* x = get(pick(rng_));
    const Individual* y = get(pick(rng_));
    return x->fitness <= y->fitness ? x : y;
  }

  // Routes of `s` (non-empty) ordered by the polar angle of their centroid
  // around the depot.
  std::vector<int> routesByAngle(const Solution& s) const {
    std::vector<std::pair<double, int>> angles;
    const std::vector<Route>& routes = s.routes();
    const Node& depot = inst_.nodes[0];
    for (int r = 0; r < static_cast<int>(routes.size()); ++r) {
      if (routes[r].empty()) continue;
      double x = 0.0;
      double y = 0.0;
      for (int k = 1; k + 1 < routes[r].numNodes(); ++k) {
        x += inst_.nodes[routes[r].node(k)].x;
        y += inst_.nodes[routes[r].node(k)].y;
      }
      const int n = routes[r].numCustomers();
      angles.emplace_back(std::atan2(y / n - depot.y, x / n - depot.x), r);
    }
    std::sort(angles.begin(), angles.end());
    std::vector<int> order;
    for (const auto& [angle, r] : angles) order.push_back(r);
    return order;
  }

  // Selective route exchange crossover (Nagata & Kobayashi; as in PyVRP).
  Solution srex(const Solution& a, const Solution& b) {
    const std::vector<int> routes_a = routesByAngle(a);
    const std::vector<int> routes_b = routesByAngle(b);
    const int na = static_cast<int>(routes_a.size());
    const int nb = static_cast<int>(routes_b.size());
    std::uniform_int_distribution<int> pick_moved(1, std::min(na, nb));
    const int moved = pick_moved(rng_);
    const int start_a = std::uniform_int_distribution<int>(0, na - 1)(rng_);

    const int n = inst_.numNodes();
    std::vector<char> in_a(n, 0);  // customer in a selected route of a
    for (int k = 0; k < moved; ++k) {
      for (int node : a.routes()[routes_a[(start_a + k) % na]].nodes()) {
        if (node != 0) in_a[node] = 1;
      }
    }
    // Block of b's routes that overlaps most with a's selected customers.
    std::vector<int> overlap(nb, 0);
    for (int j = 0; j < nb; ++j) {
      for (int node : b.routes()[routes_b[j]].nodes()) overlap[j] += in_a[node];
    }
    int start_b = 0;
    int best_overlap = -1;
    for (int s = 0; s < nb; ++s) {
      int sum = 0;
      for (int k = 0; k < moved; ++k) sum += overlap[(s + k) % nb];
      if (sum > best_overlap) {
        best_overlap = sum;
        start_b = s;
      }
    }
    std::vector<char> selected_b(nb, 0);
    std::vector<char> in_b(n, 0);  // customer in a selected route of b
    for (int k = 0; k < moved; ++k) {
      const int j = (start_b + k) % nb;
      selected_b[j] = 1;
      for (int node : b.routes()[routes_b[j]].nodes()) {
        if (node != 0) in_b[node] = 1;
      }
    }

    // Customers that neither the kept routes of b nor the selected routes of
    // a will contain.
    std::vector<int> missing;
    for (int c = 1; c < n; ++c) {
      if (in_b[c] && !in_a[c]) missing.push_back(c);
    }

    // Offspring 1 removes a's customers from b's kept routes; offspring 2
    // removes b's kept customers from a's selected routes.
    Solution best(inst_);
    double best_cost = kInf;
    for (int variant = 0; variant < 2; ++variant) {
      Solution child(inst_);
      for (int j = 0; j < nb; ++j) {
        if (selected_b[j]) continue;
        std::vector<int> nodes;
        for (int node : b.routes()[routes_b[j]].nodes()) {
          if (variant == 0 && node != 0 && in_a[node]) continue;
          nodes.push_back(node);
        }
        child.addRoute().assign(std::move(nodes));
      }
      for (int k = 0; k < moved; ++k) {
        std::vector<int> nodes;
        for (int node : a.routes()[routes_a[(start_a + k) % na]].nodes()) {
          // b's kept routes still hold customers that are not in in_b.
          if (variant == 1 && node != 0 && !in_b[node]) continue;
          nodes.push_back(node);
        }
        child.addRoute().assign(std::move(nodes));
      }
      std::vector<int> to_insert = missing;
      std::shuffle(to_insert.begin(), to_insert.end(), rng_);
      for (int c : to_insert) insertPenalized(c, &child);
      const double cost = solutionCost(child);
      if (cost < best_cost) {
        best_cost = cost;
        best = std::move(child);
      }
    }
    return best;
  }

  const GeneticSearchParams& params_;
  const Instance& inst_;
  LocalSearch local_search_;
  std::mt19937 rng_;
  double tw_penalty_ = 1.0;
  double load_penalty_ = 1.0;
  int adapt_count_ = 0;
  int adapt_tw_ok_ = 0;
  int adapt_load_ok_ = 0;
  Population feasible_;
  Population infeasible_;
  std::unique_ptr<Solution> best_;
  GeneticSearchStats stats_;
};

}  // namespace

GeneticSearchStats geneticSearch(const GeneticSearchParams& params,
                                 std::chrono::steady_clock::time_point deadline,
                                 Solution* solution) {
  return GeneticSearch(params, solution->instance()).run(deadline, solution);
}
