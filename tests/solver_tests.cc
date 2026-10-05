// Tests for the claude_vrp solver. No framework: each test is a function that
// records failed CHECKs; `claude_vrp_tests NAME` runs one test (as ctest
// does), `claude_vrp_tests` runs all of them.
//
// Most tests use the first customers of a benchmark instance (subInstance)
// so that the suite runs in seconds. Several tests compare the O(1) move and
// insertion evaluations against a full recomputation, and the search tests
// run with LocalSearchParams::check_moves, which throws if any applied move
// was evaluated incorrectly.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "adaptive_lns.h"
#include "construction.h"
#include "destroy_repair.h"
#include "ejection_search.h"
#include "fleet_minimization.h"
#include "genetic_search.h"
#include "instance.h"
#include "iterated_local_search.h"
#include "local_search.h"
#include "route_elimination.h"
#include "solution.h"

namespace {

int g_failures = 0;

#define CHECK(condition)                                     \
  do {                                                       \
    if (!(condition)) {                                      \
      ++g_failures;                                          \
      std::cerr << __FILE__ << ":" << __LINE__               \
                << ": CHECK failed: " << #condition << '\n'; \
    }                                                        \
  } while (false)

using Clock = std::chrono::steady_clock;

Clock::time_point after(double seconds) {
  return Clock::now() + std::chrono::duration_cast<Clock::duration>(
                            std::chrono::duration<double>(seconds));
}

std::string dataPath(const std::string& name) {
  return std::string(CLAUDE_VRP_DATA_DIR) + "/" + name;
}

// The depot and the first `customers` customers of a benchmark instance.
Instance subInstance(const std::string& file, int customers) {
  const Instance full = loadInstance(dataPath(file));
  Instance inst;
  inst.name = full.name + "_first" + std::to_string(customers);
  inst.num_vehicles = full.num_vehicles;
  inst.capacity = full.capacity;
  inst.nodes.assign(full.nodes.begin(), full.nodes.begin() + customers + 1);
  const std::size_t n = inst.nodes.size();
  inst.dist.resize(n * n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      inst.dist[i * n + j] = computeDistance(inst.nodes[i], inst.nodes[j]);
    }
  }
  return inst;
}

bool throwsOnValidate(const Solution& solution) {
  try {
    solution.validate();
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}

// A local optimum: a fresh full local search finds nothing to improve.
bool isLocalOptimum(const Solution& solution) {
  Solution copy = solution;
  LocalSearchParams params;
  params.check_moves = true;
  LocalSearch local_search(solution.instance(), params);
  return local_search.run(&copy).moves_applied == 0;
}

// I1 followed by hard-constraint local search (every move checked).
Solution startSolution(const Instance& inst) {
  I1Params i1;
  i1.alpha1 = 0.5;
  Solution solution = solomonI1(inst, i1);
  LocalSearchParams params;
  params.check_moves = true;
  LocalSearch local_search(inst, params);
  local_search.run(&solution);
  return solution;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

void testInstance() {
  const Instance inst = loadInstance(dataPath("RC2_10_9.TXT"));
  CHECK(inst.numCustomers() == 1000);
  CHECK(inst.num_vehicles == 250);
  CHECK(inst.capacity == 1000);
  CHECK(inst.horizon() == 7284);
  CHECK(inst.distance(3, 3) == 0.0);
  CHECK(inst.distance(1, 2) == inst.distance(2, 1));
  CHECK(capacityLowerBound(inst) == 18);

  const std::vector<std::vector<int>> nearest = nearestCustomers(inst, 10);
  CHECK(nearest[0].empty());
  for (int c = 1; c <= 50; ++c) {
    CHECK(nearest[c].size() == 10);
    CHECK(std::find(nearest[c].begin(), nearest[c].end(), c) ==
          nearest[c].end());
    CHECK(
        std::is_sorted(nearest[c].begin(), nearest[c].end(), [&](int a, int b) {
          return inst.distance(c, a) < inst.distance(c, b);
        }));
  }

  bool threw = false;
  try {
    loadInstance(dataPath("NO_SUCH_FILE.TXT"));
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CHECK(threw);
}

// Segment concatenation agrees with the explicit schedule on random routes
// (many of them infeasible).
void testSegments() {
  const Instance inst = subInstance("R1_10_1.TXT", 200);
  std::mt19937 rng(3);
  std::uniform_int_distribution<int> customer(1, inst.numCustomers());
  for (int t = 0; t < 2000; ++t) {
    std::vector<int> nodes{0};
    const int length = 1 + t % 15;
    for (int i = 0; i < length; ++i) nodes.push_back(customer(rng));
    nodes.push_back(0);
    Route route(inst);
    route.assign(nodes);
    const Segment& all = route.forward(route.numNodes() - 1);
    const bool segment_feasible =
        all.time_warp <= 1e-6 && all.load <= inst.capacity;
    CHECK(segment_feasible == route.isFeasible());
    CHECK(std::fabs(all.distance - route.distance()) < 1e-9);
    CHECK(std::fabs(route.backward(0).time_warp - all.time_warp) < 1e-6);
  }
}

// canInsert and pushForward agree with inserting and recomputing.
void testInsertion() {
  const Instance inst = subInstance("RC1_10_1.TXT", 200);
  std::mt19937 rng(7);
  std::uniform_int_distribution<int> customer(1, inst.numCustomers());
  int feasible = 0;
  for (int t = 0; t < 200; ++t) {
    Route route(inst);
    for (int step = 0; step < 30; ++step) {
      const int c = customer(rng);
      const int pos = static_cast<int>(rng() % (route.numCustomers() + 1));
      Route copy = route;
      copy.insert(c, pos);
      CHECK(route.canInsert(c, pos) == copy.isFeasible());
      if (!route.canInsert(c, pos)) continue;
      ++feasible;
      const double before =
          pos < route.numCustomers() ? route.startTime(pos) : route.endTime();
      const double after_insert = pos + 1 < copy.numCustomers()
                                      ? copy.startTime(pos + 1)
                                      : copy.endTime();
      CHECK(std::fabs(route.pushForward(c, pos) - (after_insert - before)) <
            1e-6);
      CHECK(std::fabs(route.insertionCost(c, pos) -
                      (copy.distance() - route.distance())) < 1e-6);
      route.insert(c, pos);
    }
  }
  CHECK(feasible > 100);  // the test exercised feasible insertions
}

void testConstruction() {
  for (const char* file : {"C1_10_1.TXT", "R2_10_1.TXT", "RC2_10_9.TXT"}) {
    const Instance inst = subInstance(file, 150);
    for (SeedRule rule : {SeedRule::kFarthest, SeedRule::kEarliestDeadline}) {
      I1Params params;
      params.seed = rule;
      const Solution solution = solomonI1(inst, params);
      CHECK(!throwsOnValidate(solution));
      CHECK(solution.unroutedCustomers().empty());
    }
  }
}

// validate() rejects broken solutions.
void testValidateDetectsErrors() {
  const Instance inst = subInstance("R1_10_1.TXT", 150);
  const Solution good = startSolution(inst);
  CHECK(!throwsOnValidate(good));

  const auto modified =
      [&](const std::function<void(std::vector<Route>*)>& change) {
        Solution s = good;
        change(s.mutableRoutes());
        return s;
      };
  // A customer visited twice.
  CHECK(throwsOnValidate(modified([](std::vector<Route>* routes) {
    std::vector<int> nodes = (*routes)[1].nodes();
    nodes.insert(nodes.end() - 1, (*routes)[0].customer(0));
    (*routes)[1].assign(nodes);
  })));
  // A customer missing.
  CHECK(throwsOnValidate(
      modified([](std::vector<Route>* routes) { (*routes)[0].remove(0); })));
  // A route in reverse order (R1 time windows are tight).
  CHECK(throwsOnValidate(modified([](std::vector<Route>* routes) {
    for (Route& route : *routes) {
      if (route.numCustomers() < 4) continue;
      std::vector<int> nodes = route.nodes();
      std::reverse(nodes.begin() + 1, nodes.end() - 1);
      route.assign(nodes);
      return;
    }
  })));
  // All customers in one route: capacity (and time windows) violated.
  CHECK(throwsOnValidate(modified([&](std::vector<Route>* routes) {
    std::vector<int> nodes{0};
    for (int c = 1; c <= inst.numCustomers(); ++c) nodes.push_back(c);
    nodes.push_back(0);
    routes->clear();
    routes->emplace_back(inst);
    routes->back().assign(nodes);
  })));
}

// Hard-constraint local search: every applied move is checked against a
// recomputation, the result is valid and a true local optimum.
void testLocalSearch() {
  for (const char* file : {"C1_10_1.TXT", "R1_10_1.TXT", "RC2_10_9.TXT"}) {
    const Instance inst = subInstance(file, 200);
    I1Params i1;
    Solution solution = solomonI1(inst, i1);
    const double before = solution.totalDistance();
    LocalSearchParams params;
    params.check_moves = true;
    params.max_run = 3;
    params.reversed_runs = true;
    LocalSearch local_search(inst, params);
    const LocalSearchStats stats = local_search.run(&solution);  // may throw
    CHECK(stats.moves_applied > 0);
    CHECK(solution.totalDistance() < before);
    CHECK(!throwsOnValidate(solution));
    CHECK(isLocalOptimum(solution));
  }
}

// Soft-constraint local search: penalized costs are evaluated correctly and
// repair returns a feasible solution when it reports success.
void testSoftLocalSearch() {
  const Instance inst = subInstance("RC1_10_1.TXT", 200);
  Solution solution = startSolution(inst);
  LocalSearchParams params;
  params.check_moves = true;
  params.tw_penalty = 0.5;
  params.load_penalty = 0.5;
  LocalSearch local_search(inst, params);
  const bool feasible = local_search.runAndRepair(&solution);  // may throw
  CHECK(feasible == solution.isFeasible());
  if (feasible) CHECK(!throwsOnValidate(solution));
}

// Search phases that pass "touched routes" to the local search end in true
// local optima (the skipping of unchanged pairs misses nothing).
void testIteratedLocalSearch() {
  const Instance inst = subInstance("RC2_10_9.TXT", 200);
  Solution solution = startSolution(inst);
  LocalSearchParams params;
  params.check_moves = true;
  LocalSearch local_search(inst, params);
  for (unsigned seed = 1; seed <= 3; ++seed) {
    IlsParams ils;
    ils.seed = seed;
    ils.elimination_period = 20;
    iteratedLocalSearch(ils, after(0.3), &local_search, &solution);
    CHECK(!throwsOnValidate(solution));
    CHECK(isLocalOptimum(solution));
  }
}

void testDestroyRepair() {
  const Instance inst = subInstance("R1_10_1.TXT", 200);
  const Solution solution = startSolution(inst);
  const std::vector<std::vector<int>> nearest = nearestCustomers(inst, 50);
  std::mt19937 rng(5);

  const auto routed_and_distinct = [&](const std::vector<int>& customers) {
    const std::set<int> distinct(customers.begin(), customers.end());
    const Positions positions(solution);
    return distinct.size() == customers.size() &&
           std::all_of(customers.begin(), customers.end(),
                       [&](int c) { return positions.route_of[c] >= 0; });
  };
  for (int count : {5, 20, 40}) {
    CHECK(static_cast<int>(randomRemoval(solution, count, &rng).size()) ==
          count);
    const std::vector<std::vector<int>> removals = {
        randomRemoval(solution, count, &rng),
        worstRemoval(solution, count, 6.0, &rng),
        relatedRemoval(solution, nearest, count, 6.0, &rng),
        routeRemoval(solution, count, &rng),
        radialRemoval(inst, nearest, count, &rng),
        stringRemoval(solution, nearest, 15.0, 10, &rng, 0.5)};
    for (const std::vector<int>& removed : removals) {
      CHECK(!removed.empty());
      CHECK(routed_and_distinct(removed));
    }
    CHECK(static_cast<int>(routeRemoval(solution, count, &rng).size()) >=
          count);
  }

  // Every repair yields a valid solution once leftovers get their own route.
  InsertionScope scope;
  scope.nearest = &nearest;
  scope.neighbors = 20;
  scope.blink_rate = 0.01;
  scope.rng = &rng;
  for (int repair = 0; repair < 4; ++repair) {
    for (bool restricted : {false, true}) {
      Solution s = solution;
      std::vector<int> removed = relatedRemoval(s, nearest, 30, 6.0, &rng);
      std::vector<bool> touched(s.routes().size(), false);
      removeCustomers(removed, &touched, &s);
      const InsertionScope used = restricted ? scope : InsertionScope();
      std::vector<int> leftover =
          repair < 2
              ? greedyRepair(removed,
                             repair == 0 ? InsertionOrder::kRandom
                                         : InsertionOrder::kDeadlineFirst,
                             &rng, &touched, &s, used)
              : regretRepair(removed, repair, &touched, &s, used);
      for (int c : leftover) s.addRoute().insert(c, 0);
      CHECK(!throwsOnValidate(s));
    }
  }
}

// Short runs of every search phase end with valid solutions.
void testSearchPhases() {
  const Instance inst = subInstance("C1_10_2.TXT", 200);
  const Solution start = startSolution(inst);
  LocalSearchParams ls_params;
  ls_params.check_moves = true;
  LocalSearch local_search(inst, ls_params);

  Solution eliminated = start;
  minimizeRoutes(RouteEliminationParams(), after(0.5), &local_search,
                 &eliminated);
  CHECK(!throwsOnValidate(eliminated));
  CHECK(eliminated.numUsedRoutes() <= start.numUsedRoutes());

  Solution fleet = start;
  minimizeFleet(FleetMinimizationParams(), after(0.5), &fleet);
  CHECK(!throwsOnValidate(fleet));
  CHECK(fleet.numUsedRoutes() <= start.numUsedRoutes());

  Solution ejection = start;
  minimizeRoutesEjectionSearch(EjectionSearchParams(), after(0.5), &ejection);
  CHECK(!throwsOnValidate(ejection));
  CHECK(ejection.numUsedRoutes() <= start.numUsedRoutes());

  for (bool with_local_search : {false, true}) {
    Solution alns = start;
    AlnsParams params;
    params.local_search = with_local_search;
    params.elimination_period = 200;
    adaptiveLns(params, after(0.5), &local_search, &alns);
    CHECK(!throwsOnValidate(alns));
    CHECK(!start.isBetterThan(alns));
  }

  Solution hgs = start;
  geneticSearch(GeneticSearchParams(), after(1.0), &hgs);
  CHECK(!throwsOnValidate(hgs));
  CHECK(!start.isBetterThan(hgs));
}

// writeRoutes produces the SINTEF format with every route on one line.
void testWriteRoutes() {
  const Instance inst = subInstance("RC2_10_9.TXT", 100);
  const Solution solution = startSolution(inst);
  std::ostringstream out;
  solution.writeRoutes("RC2_10_9", "2026-01-01", out);
  std::istringstream in(out.str());
  std::string line;
  std::getline(in, line);
  CHECK(line == "Instance name : RC2_10_9");
  int routes = 0;
  std::multiset<int> customers;
  while (std::getline(in, line)) {
    if (line.rfind("Route ", 0) != 0) continue;
    ++routes;
    std::istringstream fields(line.substr(line.find(':') + 1));
    for (int c; fields >> c;) customers.insert(c);
  }
  CHECK(routes == solution.numUsedRoutes());
  CHECK(static_cast<int>(customers.size()) == inst.numCustomers());
  CHECK(static_cast<int>(
            std::set<int>(customers.begin(), customers.end()).size()) ==
        inst.numCustomers());
}

const std::map<std::string, std::function<void()>>& allTests() {
  static const std::map<std::string, std::function<void()>> tests = {
      {"instance", testInstance},
      {"segments", testSegments},
      {"insertion", testInsertion},
      {"construction", testConstruction},
      {"validate", testValidateDetectsErrors},
      {"local_search", testLocalSearch},
      {"soft_local_search", testSoftLocalSearch},
      {"iterated_local_search", testIteratedLocalSearch},
      {"destroy_repair", testDestroyRepair},
      {"search_phases", testSearchPhases},
      {"write_routes", testWriteRoutes},
  };
  return tests;
}

}  // namespace

int main(int argc, char** argv) {
  const auto& tests = allTests();
  std::vector<std::string> selected;
  if (argc > 1) {
    selected.assign(argv + 1, argv + argc);
  } else {
    for (const auto& [name, test] : tests) selected.push_back(name);
  }
  for (const std::string& name : selected) {
    const auto it = tests.find(name);
    if (it == tests.end()) {
      std::cerr << "unknown test " << name << '\n';
      return 2;
    }
    const int failures_before = g_failures;
    try {
      it->second();
    } catch (const std::exception& e) {
      ++g_failures;
      std::cerr << name << ": unexpected exception: " << e.what() << '\n';
    }
    std::cout << (g_failures == failures_before ? "PASS " : "FAIL ") << name
              << '\n';
  }
  return g_failures == 0 ? 0 : 1;
}
