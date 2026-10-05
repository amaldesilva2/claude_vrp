#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "adaptive_lns.h"
#include "construction.h"
#include "ejection_search.h"
#include "fleet_minimization.h"
#include "genetic_search.h"
#include "instance.h"
#include "iterated_local_search.h"
#include "local_search.h"
#include "route_elimination.h"
#include "solution.h"
#include "solution_exchange.h"

namespace {

constexpr char kDefaultInstancePath[] = "data/RC2_10_9.TXT";
constexpr double kDefaultTimeLimitSeconds = 120.0;
// Share of the time limit that route elimination may use; ILS gets the rest.
constexpr double kRouteEliminationShare = 0.3;

using Clock = std::chrono::steady_clock;

double secondsSince(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now() - t).count();
}

constexpr unsigned kDefaultSeed = 1;
constexpr int kMaxDefaultThreads = 8;

enum class SearchMethod { kIls, kAlns, kHgs, kMixed };

const char* methodName(SearchMethod method) {
  switch (method) {
    case SearchMethod::kIls:
      return "ils";
    case SearchMethod::kAlns:
      return "alns";
    case SearchMethod::kHgs:
      return "hgs";
    case SearchMethod::kMixed:
      return "mixed";
  }
  return "?";
}

// All algorithm parameters, adjustable with --set key=value (for
// experiments without recompiling).
struct Tuning {
  LocalSearchParams ls;
  RouteEliminationParams elimination;
  double elimination_share = kRouteEliminationShare;
  // SISR fleet minimization for the rest of the elimination share.
  bool fleet_minimization = true;
  // Nagata-Braysy ejection search before the other elimination methods.
  bool ejection_search = false;
  EjectionSearchParams ejection;
  bool final_polish = true;  // hard-constraint local search on the result
  // Multi-threaded runs: odd threads replace ejection-based elimination and
  // fleet minimization by the Nagata-Braysy ejection search, which is
  // stronger on capacity-tight instances (C1) but weaker elsewhere; better
  // fleets then spread to all threads through the solution exchange.
  bool elimination_portfolio = true;
  // If the fleet after ejection-based elimination is within this relative
  // gap above the capacity lower bound, use the ejection search instead of
  // fleet minimization (0: never).
  double auto_ejection_gap = 0.1;
  // ALNS threads share their best solution every this many seconds (0: off;
  // negative: a tenth of the time limit, which won 5 of 6 instances against
  // independent threads in 8-thread benchmarks).
  double exchange_seconds = -1.0;
  FleetMinimizationParams fleet;
  // Soft-constraint penalties for local search during ILS/ALNS (0 = hard
  // constraints). A load penalty <= 0 with a positive time-warp penalty is
  // replaced by a scale-based default.
  double soft_tw_penalty = 0.0;
  double soft_load_penalty = 0.0;
  IlsParams ils;
  AlnsParams alns;
  GeneticSearchParams hgs;
  // With --search hgs, ALNS first runs for this share of the remaining time
  // and HGS starts from its result.
  double hgs_alns_share = 0.0;
};

struct Options {
  std::string path = kDefaultInstancePath;
  std::string output;  // file for the best solution's routes (empty: none)
  double time_limit = kDefaultTimeLimitSeconds;
  unsigned seed = kDefaultSeed;
  int threads = 0;  // 0: min(kMaxDefaultThreads, hardware threads)
  // Mixed alternates ILS (even threads) and ALNS (odd threads).
  // Default ALNS: in benchmarks it beat ILS by ~6% distance.
  SearchMethod search = SearchMethod::kAlns;
  Tuning tuning;
};

void printUsage(const char* program) {
  std::cerr << "usage: " << program
            << " [instance] [seconds] [--seed N] [--threads N]"
               " [--search ils|alns|hgs|mixed] [--output FILE]"
               " [--set key=value ...]\n"
               "  Runs independent searches in parallel (seeds N, N+1, ...) and"
               " keeps the best;\n"
               "  mixed alternates ILS and ALNS across threads.\n";
}

// Parses all of `text` as a number; throws std::invalid_argument naming
// `what` otherwise (std::stod/stoi alone accept trailing garbage like "5x").
template <typename T>
T parseNumber(const std::string& text, const std::string& what) {
  std::size_t used = 0;
  T value{};
  try {
    if constexpr (std::is_floating_point_v<T>) {
      value = static_cast<T>(std::stod(text, &used));
    } else {
      value = static_cast<T>(std::stoll(text, &used));
    }
  } catch (const std::exception&) {
    used = 0;
  }
  if (used == 0 || used != text.size()) {
    throw std::invalid_argument("invalid number for " + what + ": " + text);
  }
  return value;
}

// Applies one "key=value" setting; throws std::invalid_argument for unknown
// keys or bad values. Keys are "<group>.<field>", e.g. "ils.max_ruin".
void applySetting(const std::string& setting, Tuning* tuning) {
  const std::size_t eq = setting.find('=');
  if (eq == std::string::npos) {
    throw std::invalid_argument("--set expects key=value, got " + setting);
  }
  const std::string key = setting.substr(0, eq);
  const std::string text = setting.substr(eq + 1);
  const auto real = [&] { return parseNumber<double>(text, key); };
  const auto integer = [&] {
    return static_cast<int>(parseNumber<int64_t>(text, key));
  };
  const auto flag = [&] { return integer() != 0; };

  if (key == "ls.num_neighbors") {
    tuning->ls.num_neighbors = integer();
  } else if (key == "ls.max_run") {
    tuning->ls.max_run = integer();
  } else if (key == "ls.reversed_runs") {
    tuning->ls.reversed_runs = flag();
  } else if (key == "soft.tw") {
    tuning->soft_tw_penalty = real();
  } else if (key == "soft.load") {
    tuning->soft_load_penalty = real();
  } else if (key == "soft.target") {
    tuning->ls.target_feasible = real();
  } else if (key == "soft.adapt_period") {
    tuning->ls.adapt_period = integer();
  } else if (key == "hgs.alns_share") {
    tuning->hgs_alns_share = real();
  } else if (key == "hgs.min_population") {
    tuning->hgs.min_population = integer();
  } else if (key == "hgs.generation_size") {
    tuning->hgs.generation_size = integer();
  } else if (key == "hgs.target_feasible") {
    tuning->hgs.target_feasible = real();
  } else if (key == "hgs.initial_population") {
    tuning->hgs.initial_population = integer();
  } else if (key == "hgs.num_neighbors") {
    tuning->hgs.num_neighbors = integer();
  } else if (key == "elim.auto_gap") {
    tuning->auto_ejection_gap = real();
  } else if (key == "elim.portfolio") {
    tuning->elimination_portfolio = flag();
  } else if (key == "exchange.seconds") {
    tuning->exchange_seconds = real();
  } else if (key == "polish") {
    tuning->final_polish = flag();
  } else if (key == "nb.enabled") {
    tuning->ejection_search = flag();
  } else if (key == "nb.max_ejected") {
    tuning->ejection.max_ejected = integer();
  } else if (key == "nb.perturbation_moves") {
    tuning->ejection.perturbation_moves = integer();
  } else if (key == "nb.squeeze") {
    tuning->ejection.squeeze = flag();
  } else if (key == "fleet.enabled") {
    tuning->fleet_minimization = flag();
  } else if (key == "fleet.avg_removed") {
    tuning->fleet.avg_removed = real();
  } else if (key == "fleet.max_stagnation") {
    tuning->fleet.max_stagnation = integer();
  } else if (key == "fleet.insertion_neighbors") {
    tuning->fleet.insertion_neighbors = integer();
  } else if (key == "elim.share") {
    tuning->elimination_share = real();
  } else if (key == "elim.max_ejections") {
    tuning->elimination.max_ejections = integer();
  } else if (key == "elim.max_ejected") {
    tuning->elimination.max_ejected = integer();
  } else if (key == "ils.min_ruin") {
    tuning->ils.min_ruin = integer();
  } else if (key == "ils.max_ruin") {
    tuning->ils.max_ruin = integer();
  } else if (key == "ils.string_ruin_probability") {
    tuning->ils.string_ruin_probability = real();
  } else if (key == "ils.avg_removed") {
    tuning->ils.avg_removed = real();
  } else if (key == "ils.max_string") {
    tuning->ils.max_string = integer();
  } else if (key == "ils.max_pool_ejections") {
    tuning->ils.max_pool_ejections = integer();
  } else if (key == "ils.penalize_unplaced") {
    tuning->ils.penalize_unplaced = flag();
  } else if (key == "ils.elimination_period") {
    tuning->ils.elimination_period = integer();
  } else if (key == "ils.initial_threshold") {
    tuning->ils.initial_threshold = real();
  } else if (key == "alns.min_remove") {
    tuning->alns.min_remove = integer();
  } else if (key == "alns.max_remove") {
    tuning->alns.max_remove = integer();
  } else if (key == "alns.insertion_neighbors") {
    tuning->alns.insertion_neighbors = integer();
  } else if (key == "alns.ls_period") {
    tuning->alns.ls_period = integer();
  } else if (key == "alns.blink_rate") {
    tuning->alns.blink_rate = real();
  } else if (key == "alns.split_probability") {
    tuning->alns.split_probability = real();
  } else if (key == "alns.max_pool_ejections") {
    tuning->alns.max_pool_ejections = integer();
  } else if (key == "alns.local_search") {
    tuning->alns.local_search = flag();
  } else if (key == "alns.start_worse") {
    tuning->alns.start_worse = real();
  } else if (key == "alns.end_temperature_ratio") {
    tuning->alns.end_temperature_ratio = real();
  } else if (key == "alns.elimination_period") {
    tuning->alns.elimination_period = integer();
  } else {
    throw std::invalid_argument("unknown setting " + key);
  }
}

// Parses positional [instance] [seconds] and the --seed/--threads flags.
// Throws std::invalid_argument on bad input.
Options parseOptions(int argc, char** argv) {
  Options options;
  int positional = 0;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--output") {
      if (i + 1 >= argc) throw std::invalid_argument(arg + " needs a value");
      options.output = argv[++i];
    } else if (arg == "--set") {
      if (i + 1 >= argc) throw std::invalid_argument(arg + " needs a value");
      applySetting(argv[++i], &options.tuning);
    } else if (arg == "--search") {
      if (i + 1 >= argc) throw std::invalid_argument(arg + " needs a value");
      const std::string value = argv[++i];
      if (value == "ils") {
        options.search = SearchMethod::kIls;
      } else if (value == "alns") {
        options.search = SearchMethod::kAlns;
      } else if (value == "hgs") {
        options.search = SearchMethod::kHgs;
      } else if (value == "mixed") {
        options.search = SearchMethod::kMixed;
      } else {
        throw std::invalid_argument("unknown search method " + value);
      }
    } else if (arg == "--seed" || arg == "--threads") {
      if (i + 1 >= argc) throw std::invalid_argument(arg + " needs a value");
      const int64_t value = parseNumber<int64_t>(argv[++i], arg);
      if (value < (arg == "--seed" ? 0 : 1) || value > 1'000'000) {
        throw std::invalid_argument("invalid value for " + arg);
      }
      if (arg == "--seed") {
        options.seed = static_cast<unsigned>(value);
      } else {
        options.threads = static_cast<int>(value);
      }
    } else if (!arg.empty() && arg[0] == '-') {
      throw std::invalid_argument("unknown option " + arg);
    } else if (positional == 0) {
      options.path = arg;
      ++positional;
    } else if (positional == 1) {
      options.time_limit = parseNumber<double>(arg, "seconds");
      if (options.time_limit <= 0) {
        throw std::invalid_argument("time limit must be positive");
      }
      ++positional;
    } else {
      throw std::invalid_argument("unexpected argument " + arg);
    }
  }
  if (options.threads == 0) {
    const int hardware = static_cast<int>(std::thread::hardware_concurrency());
    options.threads = std::max(1, std::min(kMaxDefaultThreads, hardware));
  }
  return options;
}

// Result of one search (local search -> route elimination -> ILS or ALNS).
struct RunResult {
  unsigned seed = 0;
  SearchMethod method = SearchMethod::kIls;
  std::optional<Solution> after_local_search;
  std::optional<Solution> after_elimination;
  std::optional<Solution> solution;  // final, validated
  LocalSearchStats ls_stats;
  RouteEliminationStats re_stats;
  FleetMinimizationStats fleet_stats;
  EjectionSearchStats nb_stats;
  double before_polish = 0.0;
  bool auto_ejection_search = false;  // capacity-tight rule triggered    //
                                      // distance before the final local search
  IlsStats ils_stats;                 // if method is ILS
  AlnsStats alns_stats;               // if method is ALNS
  GeneticSearchStats hgs_stats;       // if method is HGS
  std::string error;                  // non-empty if the run failed
};

// Runs the improvement pipeline from `initial` with the given seed. Safe to
// call concurrently: it only reads `inst` and `initial`.
RunResult runPipeline(const Instance& inst, const Solution& initial,
                      unsigned seed, SearchMethod method, const Tuning& tuning,
                      SolutionExchange* exchange, Clock::time_point start,
                      double time_limit) {
  RunResult result;
  result.seed = seed;
  result.method = method;
  try {
    const Clock::time_point deadline =
        start + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(time_limit));
    Solution solution = initial;

    // 1. Local search.
    LocalSearchParams ls_params = tuning.ls;
    ls_params.seed = seed;
    LocalSearch local_search(inst, ls_params);
    result.ls_stats = local_search.run(&solution);
    result.after_local_search = solution;

    // 2. Route elimination.
    const Clock::time_point elimination_deadline = std::min(
        deadline, start + std::chrono::duration_cast<Clock::duration>(
                              std::chrono::duration<double>(
                                  time_limit * tuning.elimination_share)));
    result.re_stats = minimizeRoutes(tuning.elimination, elimination_deadline,
                                     &local_search, &solution);
    // Capacity-tight fleets (close to the capacity lower bound) are left to
    // the Nagata-Braysy ejection search instead of fleet minimization.
    bool ejection_search = tuning.ejection_search;
    bool fleet_minimization = tuning.fleet_minimization;
    const int fleet = solution.numUsedRoutes();
    const int lower_bound = capacityLowerBound(inst);
    if (tuning.auto_ejection_gap > 0 && !ejection_search &&
        fleet > lower_bound &&
        fleet <= lower_bound * (1.0 + tuning.auto_ejection_gap)) {
      ejection_search = true;
      fleet_minimization = false;
      result.auto_ejection_search = true;
    }
    // With the ejection search enabled, fleet minimization gets half of the
    // remaining elimination time (or, with a stagnation limit, runs until it
    // stagnates) and the ejection search the rest.
    const Clock::time_point now = Clock::now();
    const Clock::time_point fleet_deadline =
        ejection_search && elimination_deadline > now &&
                tuning.fleet.max_stagnation == 0
            ? now + (elimination_deadline - now) / 2
            : elimination_deadline;
    if (fleet_minimization) {
      FleetMinimizationParams fleet_params = tuning.fleet;
      fleet_params.seed = seed;
      result.fleet_stats =
          minimizeFleet(fleet_params, fleet_deadline, &solution);
      local_search.run(&solution);
    }
    if (ejection_search) {
      EjectionSearchParams ejection_params = tuning.ejection;
      ejection_params.seed = seed;
      result.nb_stats = minimizeRoutesEjectionSearch(
          ejection_params, elimination_deadline, &solution);
      local_search.run(&solution);
    }
    result.after_elimination = solution;

    // 3. ILS or ALNS for the remaining time, optionally with soft
    // constraints in local search.
    if (tuning.soft_tw_penalty > 0) {
      double load_penalty = tuning.soft_load_penalty;
      if (load_penalty <= 0) {
        // Distance per unit of demand: average demand vs. average distance
        // to the depot.
        double demand = 0.0;
        double distance = 0.0;
        for (int c = 1; c < inst.numNodes(); ++c) {
          demand += inst.nodes[c].demand;
          distance += inst.distance(0, c);
        }
        load_penalty = std::max(0.1, distance / std::max(1.0, demand));
      }
      local_search.setPenalties(tuning.soft_tw_penalty, load_penalty);
    }
    if (method == SearchMethod::kHgs) {
      if (tuning.hgs_alns_share > 0) {
        const Clock::time_point alns_deadline =
            Clock::now() +
            std::chrono::duration_cast<Clock::duration>(
                (deadline - Clock::now()) * tuning.hgs_alns_share);
        AlnsParams alns_params = tuning.alns;
        alns_params.seed = seed;
        result.alns_stats =
            adaptiveLns(alns_params, alns_deadline, &local_search, &solution);
        local_search.run(&solution);
      }
      GeneticSearchParams hgs_params = tuning.hgs;
      hgs_params.seed = seed;
      result.hgs_stats = geneticSearch(hgs_params, deadline, &solution);
    } else if (method == SearchMethod::kAlns) {
      AlnsParams alns_params = tuning.alns;
      alns_params.seed = seed;
      alns_params.exchange = exchange;
      alns_params.exchange_seconds = tuning.exchange_seconds < 0
                                         ? 0.1 * time_limit
                                         : tuning.exchange_seconds;
      result.alns_stats =
          adaptiveLns(alns_params, deadline, &local_search, &solution);
    } else {
      IlsParams ils_params = tuning.ils;
      ils_params.seed = seed;
      result.ils_stats =
          iteratedLocalSearch(ils_params, deadline, &local_search, &solution);
    }

    // Final polish: ALNS without local search leaves easy improvements.
    if (tuning.final_polish) {
      local_search.setPenalties(0.0, 0.0);
      result.before_polish = solution.totalDistance();
      local_search.run(&solution);
    }

    solution.validate();
    result.solution = std::move(solution);
  } catch (const std::exception& e) {
    result.error = e.what();
  }
  return result;
}

std::string summary(const std::optional<Solution>& solution) {
  if (!solution) return "-";
  std::ostringstream os;
  os << std::fixed << std::setprecision(2) << std::setw(4)
     << solution->numUsedRoutes() << " / " << std::setw(10)
     << solution->totalDistance();
  return os.str();
}

const char* seedRuleName(SeedRule rule) {
  return rule == SeedRule::kFarthest ? "farthest" : "deadline";
}

}  // namespace

int main(int argc, char** argv) {
  Options options;
  try {
    options = parseOptions(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    printUsage(argv[0]);
    return 1;
  }
  const Clock::time_point start = Clock::now();

  Instance inst;
  try {
    inst = loadInstance(options.path);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }

  int64_t total_demand = 0;
  for (int i = 1; i < inst.numNodes(); ++i) {
    total_demand += inst.nodes[i].demand;
  }
  const int64_t min_vehicles =
      (total_demand + inst.capacity - 1) / inst.capacity;

  std::cout << "Instance       : " << inst.name << '\n'
            << "Vehicles       : " << inst.num_vehicles << '\n'
            << "Capacity       : " << inst.capacity << '\n'
            << "Customers      : " << inst.numCustomers() << '\n'
            << "Horizon        : " << inst.horizon() << '\n'
            << "Total demand   : " << total_demand << '\n'
            << "Min vehicles LB: " << min_vehicles << '\n';

  // Customers that cannot be served even by a dedicated depot -> i -> depot
  // route.
  int num_unservable = 0;
  for (int i = 1; i < inst.numNodes(); ++i) {
    const Node& c = inst.nodes[i];
    const double earliest = std::max(inst.distance(0, i), c.ready);
    if (earliest > c.due ||
        earliest + c.service + inst.distance(i, 0) > inst.horizon() ||
        c.demand > inst.capacity) {
      std::cout << "  unservable customer " << i << '\n';
      ++num_unservable;
    }
  }
  std::cout << "Unservable     : " << num_unservable << '\n';
  if (num_unservable > 0) return 2;

  // Solomon I1 over a small parameter grid; keep the best solution.
  std::cout << "\nSolomon I1\n"
            << "  lambda alpha1 seed      vehicles   distance   ms\n";
  std::optional<Solution> best;
  I1Params best_params;
  for (SeedRule seed : {SeedRule::kFarthest, SeedRule::kEarliestDeadline}) {
    for (double lambda : {1.0, 1.5, 2.0}) {
      for (double alpha1 : {0.0, 0.5, 1.0}) {
        I1Params params;
        params.lambda = lambda;
        params.alpha1 = alpha1;
        params.seed = seed;

        const auto t0 = std::chrono::steady_clock::now();
        Solution solution = solomonI1(inst, params);
        const auto t1 = std::chrono::steady_clock::now();
        const auto ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                .count();

        std::cout << std::fixed << std::setprecision(1) << "  " << std::setw(6)
                  << lambda << ' ' << std::setw(6) << alpha1 << ' ' << std::left
                  << std::setw(9) << seedRuleName(seed) << std::right
                  << std::setw(9) << solution.numUsedRoutes() << std::setw(11)
                  << std::setprecision(2) << solution.totalDistance()
                  << std::setw(5) << ms << '\n';

        if (!best || solution.isBetterThan(*best)) {
          best = std::move(solution);
          best_params = params;
        }
      }
    }
  }

  try {
    best->validate();
  } catch (const std::exception& e) {
    std::cerr << "error: best solution is invalid: " << e.what() << '\n';
    return 1;
  }

  std::cout << std::fixed << std::setprecision(1)
            << "\nBest I1: lambda=" << best_params.lambda
            << " alpha1=" << best_params.alpha1
            << " seed=" << seedRuleName(best_params.seed) << '\n'
            << std::setprecision(2) << "  vehicles " << best->numUsedRoutes()
            << ", distance " << best->totalDistance() << " (validated)\n";

  // Independent searches in parallel, one seed per thread.
  const unsigned last_seed = options.seed + options.threads - 1;
  std::cout << "\nSearch: " << methodName(options.search) << ", "
            << options.threads << " parallel runs, seeds " << options.seed
            << ".." << last_seed << ", time limit " << std::setprecision(0)
            << options.time_limit << " s\n";

  std::vector<RunResult> results(options.threads);
  SolutionExchange exchange;
  std::vector<std::thread> workers;
  workers.reserve(options.threads);
  for (int t = 0; t < options.threads; ++t) {
    SearchMethod method = options.search;
    if (method == SearchMethod::kMixed) {
      method = t % 2 == 0 ? SearchMethod::kIls : SearchMethod::kAlns;
    }
    Tuning tuning = options.tuning;
    if (tuning.elimination_portfolio && t % 2 == 1) {
      tuning.ejection_search = true;
      tuning.fleet_minimization = false;
      tuning.elimination.max_ejections = 0;
    }
    workers.emplace_back([&, t, method, tuning] {
      results[t] = runPipeline(inst, *best, options.seed + t, method, tuning,
                               &exchange, start, options.time_limit);
    });
  }
  for (std::thread& worker : workers) worker.join();

  std::cout << "  seed  method   local search (veh / dist)   route "
               "elimination        final (veh / dist)   iterations\n";
  const RunResult* winner = nullptr;
  for (const RunResult& r : results) {
    const int iterations =
        r.method == SearchMethod::kAlns  ? r.alns_stats.iterations
        : r.method == SearchMethod::kHgs ? r.hgs_stats.generations
                                         : r.ils_stats.iterations;
    std::cout << "  " << std::setw(4) << r.seed << "  " << std::left
              << std::setw(6) << methodName(r.method) << std::right << "   "
              << std::setw(22) << summary(r.after_local_search) << "   "
              << std::setw(22) << summary(r.after_elimination) << "   "
              << std::setw(22) << summary(r.solution) << "   " << std::setw(6)
              << iterations;
    if (!r.error.empty()) std::cout << "   FAILED: " << r.error;
    std::cout << '\n';
    if (r.solution &&
        (!winner || r.solution->isBetterThan(*winner->solution))) {
      winner = &r;
    }
  }
  if (!winner) {
    std::cerr << "error: every run failed\n";
    return 1;
  }

  const Solution& final_solution = *winner->solution;
  try {
    final_solution.validate();
  } catch (const std::exception& e) {
    std::cerr << "error: best solution is invalid: " << e.what() << '\n';
    return 1;
  }

  if (!options.output.empty()) {
    // Instance name from the file name, e.g. "data/RC2_10_9.TXT" -> RC2_10_9.
    std::string name = options.path.substr(options.path.find_last_of('/') + 1);
    name = name.substr(0, name.find('.'));
    const std::time_t now = std::time(nullptr);
    char date[11];
    std::strftime(date, sizeof(date), "%Y-%m-%d", std::localtime(&now));
    std::ofstream file(options.output);
    final_solution.writeRoutes(name, date, file);
    if (!file) {
      std::cerr << "error: cannot write " << options.output << '\n';
      return 1;
    }
    std::cout << "Routes written to " << options.output << '\n';
  }

  std::cout << std::fixed << std::setprecision(2) << "\nBest: seed "
            << winner->seed << " (" << methodName(winner->method)
            << "), vehicles " << final_solution.numUsedRoutes() << ", distance "
            << final_solution.totalDistance() << " (validated), "
            << std::setprecision(1) << secondsSince(start) << " s\n"
            << "  local search: " << winner->ls_stats.moves_applied
            << " moves in " << winner->ls_stats.rounds << " rounds\n"
            << "  route elimination: " << winner->re_stats.removed << " of "
            << winner->re_stats.attempts << " attempts succeeded\n"
            << "  ejection search"
            << (winner->auto_ejection_search ? " (capacity-tight)" : "") << ": "
            << winner->nb_stats.routes_removed << " routes removed in "
            << winner->nb_stats.attempts << " attempts ("
            << winner->nb_stats.insertions << " insertions, "
            << winner->nb_stats.squeezes << " squeezes, "
            << winner->nb_stats.ejections << " ejections)\n"
            << "  fleet minimization: " << winner->fleet_stats.routes_removed
            << " routes removed in " << winner->fleet_stats.iterations
            << " iterations (longest "
            << winner->fleet_stats.longest_successful_stagnation
            << " without progress before a success"
            << (winner->fleet_stats.stagnated ? ", stopped: stagnation" : "")
            << ")\n"
            << "  final local search: " << winner->before_polish << " -> "
            << final_solution.totalDistance() << '\n';
  if (winner->method == SearchMethod::kHgs) {
    const GeneticSearchStats& hgs = winner->hgs_stats;
    std::cout << "  HGS: " << hgs.generations << " offspring ("
              << hgs.feasible_offspring << " feasible after local search, "
              << hgs.repaired << " repaired), " << hgs.improvements
              << " new best, final penalties: time warp "
              << hgs.final_tw_penalty << ", load " << hgs.final_load_penalty
              << '\n';
  } else if (winner->method == SearchMethod::kAlns) {
    const AlnsStats& alns = winner->alns_stats;
    std::cout << "  ALNS: " << alns.iterations << " iterations, "
              << alns.accepted << " accepted, " << alns.improvements
              << " new best, " << alns.routes_eliminated << " of "
              << alns.elimination_attempts << " route eliminations, "
              << alns.adopted << " solutions adopted from other threads\n";
    for (const auto* group : {&alns.destroy, &alns.repair}) {
      std::cout << (group == &alns.destroy ? "  destroy" : "  repair ")
                << " operators (final weight / uses / new best):";
      for (const AlnsOperatorStats& op : *group) {
        std::cout << "  " << op.name << ' ' << std::setprecision(1) << op.weight
                  << '/' << op.uses << '/' << op.new_best;
      }
      std::cout << '\n';
    }
  } else {
    const IlsStats& ils = winner->ils_stats;
    std::cout << "  ILS: " << ils.iterations << " iterations, " << ils.accepted
              << " accepted, " << ils.improvements << " new best, "
              << ils.routes_eliminated << " of " << ils.elimination_attempts
              << " route eliminations\n";
  }
  return 0;
}
