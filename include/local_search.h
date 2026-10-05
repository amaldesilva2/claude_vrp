#ifndef CLAUDE_VRP_INCLUDE_LOCAL_SEARCH_H_
#define CLAUDE_VRP_INCLUDE_LOCAL_SEARCH_H_

#include <map>
#include <random>
#include <string>
#include <vector>

#include "solution.h"

struct LocalSearchParams {
  int num_neighbors = 20;      // moves pair each customer with its K nearest
  int max_run = 2;             // longest run of customers relocated/swapped
  bool reversed_runs = false;  // also try runs inserted in reverse order
  unsigned seed = 1;           // random order in which customers are visited
  // Soft constraints: if either penalty is positive, routes may violate time
  // windows and capacity while searching, and moves minimize
  //   distance + tw_penalty * time warp + load_penalty * excess load.
  double tw_penalty = 0.0;
  double load_penalty = 0.0;
  // With soft constraints, runAndRepair() adapts both penalties every
  // adapt_period calls so that about target_feasible of the first-pass
  // results satisfy each constraint (0 disables adaptation).
  int adapt_period = 100;
  double target_feasible = 0.3;
  // Debug: after every applied move, recompute the changed routes and throw
  // std::logic_error if distance or feasibility differ from the evaluation.
  bool check_moves = false;
};

struct LocalSearchStats {
  int rounds = 0;                      // passes over all customers
  int moves_applied = 0;               // total improving moves applied
  std::map<std::string, int> by_move;  // applied count per move variant
};

// Improves `solution` in place with first-improvement local search over
// granular neighborhoods (each customer u with its K nearest neighbors v).
// Moves, all checked in O(1) per piece with Segment concatenation:
//   relocate   a run of 1-3 customers starting at u to just after/before v,
//              optionally reversed (inter- and intra-route)
//   swap       a run of 1-3 starting at u with a run of 1-3 starting at v,
//              each optionally reversed (inter- and intra-route)
//   2-opt*     exchange route tails after u and after v (inter-route)
//   2-opt*rev  connect u -> v, reversing the other parts (inter-route)
//   2-opt      reverse the path between u and v (intra-route)
// Only moves that keep every route feasible and reduce total distance are
// applied. Routes that become empty are removed at the end. Customers that
// are not in any route are ignored.
//
// A (u, v) pair is only re-evaluated if the route of u or of v has changed
// since u was last examined. The object keeps the neighbor lists and the
// random generator, so it is cheap to run many times (e.g. inside ILS).
class LocalSearch {
 public:
  LocalSearch(const Instance& inst, const LocalSearchParams& params);

  // Improves `solution` in place until no improving move exists.
  //
  // If `touched_routes` is given (one flag per route), `solution` must be a
  // local optimum apart from the flagged routes: in the first round, pairs
  // where neither route is flagged are skipped.
  LocalSearchStats run(Solution* solution,
                       const std::vector<bool>* touched_routes = nullptr);

  // run(), and with soft constraints, if the result is infeasible, run again
  // on the whole solution with both penalties x10 and then x100. Returns
  // whether `solution` is feasible at the end (always true without soft
  // constraints, given a feasible input).
  bool runAndRepair(Solution* solution,
                    const std::vector<bool>* touched_routes = nullptr);

  void setPenalties(double tw_penalty, double load_penalty);
  bool softConstraints() const {
    return params_.tw_penalty > 0 || params_.load_penalty > 0;
  }
  const LocalSearchParams& params() const { return params_; }

 private:
  const Instance& inst_;
  LocalSearchParams params_;
  std::vector<std::vector<int>> neighbors_;  // K nearest customers
  std::mt19937 rng_;
  // Penalty adaptation counters (see runAndRepair).
  int adapt_calls_ = 0;
  int adapt_tw_feasible_ = 0;
  int adapt_load_feasible_ = 0;
};

#endif  // CLAUDE_VRP_INCLUDE_LOCAL_SEARCH_H_
