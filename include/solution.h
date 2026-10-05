#ifndef CLAUDE_VRP_INCLUDE_SOLUTION_H_
#define CLAUDE_VRP_INCLUDE_SOLUTION_H_

#include <ostream>
#include <string>
#include <vector>

#include "instance.h"
#include "segment.h"

// A single vehicle route: depot -> customers -> depot.
//
// Caches load, distance, service start times and latest feasible start times
// so that insertion feasibility can be checked in O(1), plus prefix/suffix
// Segment summaries for local search.
//
// Two index conventions are used: "pos" is a 0-based customer position that
// excludes the depot; "k" is a node index into nodes(), where 0 and
// numNodes() - 1 are the depot.
class Route {
 public:
  explicit Route(const Instance& inst);

  int numCustomers() const { return static_cast<int>(nodes_.size()) - 2; }
  bool empty() const { return numCustomers() == 0; }
  int customer(int pos) const { return nodes_[pos + 1]; }
  int load() const { return load_; }
  double distance() const { return distance_; }

  // Node sequence including the depot at both ends.
  const std::vector<int>& nodes() const { return nodes_; }
  int numNodes() const { return static_cast<int>(nodes_.size()); }
  int node(int k) const { return nodes_[k]; }

  // Segment summary of nodes 0..k (forward) and k..numNodes()-1 (backward).
  // Computed lazily (only local search needs them): the first call after a
  // change is O(route length). Not safe to call concurrently on one Route.
  const Segment& forward(int k) const {
    if (!segments_valid_) computeSegments();
    return forward_[k];
  }
  const Segment& backward(int k) const {
    if (!segments_valid_) computeSegments();
    return backward_[k];
  }

  // Latest service start at node k that keeps the rest of the route
  // feasible (node index, 0 and numNodes() - 1 are the depot).
  double latestStart(int k) const { return latest_[k]; }

  // Service start time of the customer at `pos`.
  double startTime(int pos) const { return start_[pos + 1]; }

  // Time the vehicle arrives back at the depot.
  double endTime() const;

  // True if `customer` can be inserted before position `pos`
  // (0..numCustomers()) without violating capacity or any time window.
  // Assumes the route is currently feasible.
  bool canInsert(int customer, int pos) const;

  // Change in route distance if `customer` is inserted before `pos`.
  double insertionCost(int customer, int pos) const;

  // How much later service would start at the node currently at `pos` (or at
  // the depot, if `pos` == numCustomers()) if `customer` were inserted before
  // it.
  double pushForward(int customer, int pos) const;

  // Inserts `customer` before `pos`. Does not check feasibility.
  void insert(int customer, int pos);

  // Removes and returns the customer at `pos`.
  int remove(int pos);

  // Replaces the whole node sequence. `nodes` must start and end with the
  // depot. Does not check feasibility.
  void assign(std::vector<int> nodes);

  // Capacity and time-window feasibility, based on the cached schedule.
  bool isFeasible() const;

  // Total lateness forced by due times (0 iff time-window feasible) and load
  // above the vehicle capacity.
  double timeWarp() const {
    return forward(static_cast<int>(nodes_.size()) - 1).time_warp;
  }
  int excessLoad() const;

 private:
  // Recomputes load, distance, start_ and latest_ from nodes_ from scratch,
  // and invalidates the segment summaries. O(route length).
  void recomputeSchedule();

  const Instance* inst_;
  std::vector<int> nodes_;      // 0, c1, ..., ck, 0
  std::vector<double> start_;   // service start time at each node
  std::vector<double> latest_;  // latest start keeping the suffix feasible
  // Fills forward_ and backward_ from nodes_.
  void computeSegments() const;

  mutable std::vector<Segment> forward_;   // summary of nodes 0..k
  mutable std::vector<Segment> backward_;  // summary of nodes k..end
  mutable bool segments_valid_ = false;
  int load_ = 0;
  double distance_ = 0.0;
};

// A set of routes for an instance. Objective is lexicographic: first the
// number of vehicles, then total distance.
class Solution {
 public:
  explicit Solution(const Instance& inst);

  const Instance& instance() const { return *inst_; }
  const std::vector<Route>& routes() const { return routes_; }
  std::vector<Route>* mutableRoutes() { return &routes_; }

  // Appends an empty route and returns it. The reference is invalidated by
  // any later change to the route list.
  Route& addRoute();
  void removeEmptyRoutes();

  int numUsedRoutes() const;
  double totalDistance() const;
  // True if every route respects capacity and time windows.
  bool isFeasible() const;
  std::vector<int> unroutedCustomers() const;

  // True if this solution uses fewer vehicles, or the same number of vehicles
  // and less distance.
  bool isBetterThan(const Solution& other) const;

  // Recomputes the schedule from scratch (ignoring route caches) and checks
  // that every customer is visited exactly once, that capacity, time windows
  // and the fleet size are respected, and that the cached route distances
  // and totalDistance() match a fresh computation. Throws std::runtime_error
  // describing the first violation found.
  void validate() const;

  void print(std::ostream& os) const;

  // Writes the routes in the format of the SINTEF best known solutions:
  // header lines ("Instance name", "Authors", "Date", "Reference"),
  // "Solution", then one "Route i : c1 c2 ..." line per non-empty route.
  void writeRoutes(const std::string& instance_name, const std::string& date,
                   std::ostream& os) const;

 private:
  const Instance* inst_;
  std::vector<Route> routes_;
};

#endif  // CLAUDE_VRP_INCLUDE_SOLUTION_H_
