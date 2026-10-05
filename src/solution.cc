#include "solution.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

// Tolerance for floating-point time and distance comparisons.
constexpr double kEps = 1e-6;

}  // namespace

// ----------------------------------------------------------------------------
// Route
// ----------------------------------------------------------------------------

Route::Route(const Instance& inst) : inst_(&inst), nodes_{0, 0} {
  recomputeSchedule();
}

double Route::endTime() const { return start_.back(); }

bool Route::canInsert(int customer, int pos) const {
  const std::vector<Node>& nodes = inst_->nodes;
  const Node& c = nodes[customer];
  if (load_ + c.demand > inst_->capacity) return false;

  const int prev = nodes_[pos];
  const int next = nodes_[pos + 1];
  const double start = std::max(
      start_[pos] + nodes[prev].service + inst_->distance(prev, customer),
      c.ready);
  if (start > c.due + kEps) return false;

  // Waiting is allowed, so the suffix stays feasible as long as we reach
  // `next` no later than its latest feasible start time.
  return start + c.service + inst_->distance(customer, next) <=
         latest_[pos + 1] + kEps;
}

double Route::insertionCost(int customer, int pos) const {
  const int prev = nodes_[pos];
  const int next = nodes_[pos + 1];
  return inst_->distance(prev, customer) + inst_->distance(customer, next) -
         inst_->distance(prev, next);
}

double Route::pushForward(int customer, int pos) const {
  const std::vector<Node>& nodes = inst_->nodes;
  const int prev = nodes_[pos];
  const int next = nodes_[pos + 1];
  const double start = std::max(
      start_[pos] + nodes[prev].service + inst_->distance(prev, customer),
      nodes[customer].ready);
  const double next_start = std::max(
      start + nodes[customer].service + inst_->distance(customer, next),
      nodes[next].ready);
  return next_start - start_[pos + 1];
}

void Route::insert(int customer, int pos) {
  nodes_.insert(nodes_.begin() + pos + 1, customer);
  recomputeSchedule();
}

void Route::assign(std::vector<int> nodes) {
  if (nodes.size() < 2 || nodes.front() != 0 || nodes.back() != 0) {
    throw std::invalid_argument("route must start and end at the depot");
  }
  nodes_ = std::move(nodes);
  recomputeSchedule();
}

int Route::remove(int pos) {
  const int customer = nodes_[pos + 1];
  nodes_.erase(nodes_.begin() + pos + 1);
  recomputeSchedule();
  return customer;
}

int Route::excessLoad() const { return std::max(0, load_ - inst_->capacity); }

bool Route::isFeasible() const {
  if (load_ > inst_->capacity) return false;
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    if (start_[i] > inst_->nodes[nodes_[i]].due + kEps) return false;
  }
  return true;
}

void Route::recomputeSchedule() {
  const std::vector<Node>& nodes = inst_->nodes;
  const std::size_t n = nodes_.size();
  start_.assign(n, 0.0);
  latest_.assign(n, 0.0);
  load_ = 0;
  distance_ = 0.0;

  // Forward pass: earliest service start times.
  start_[0] = nodes[0].ready;
  for (std::size_t i = 1; i < n; ++i) {
    const int prev = nodes_[i - 1];
    const int cur = nodes_[i];
    distance_ += inst_->distance(prev, cur);
    load_ += nodes[cur].demand;
    start_[i] = std::max(
        start_[i - 1] + nodes[prev].service + inst_->distance(prev, cur),
        nodes[cur].ready);
  }

  // Backward pass: latest start times that keep the rest of the route
  // feasible.
  latest_[n - 1] = nodes[0].due;
  for (std::size_t i = n - 1; i-- > 0;) {
    const int cur = nodes_[i];
    const int next = nodes_[i + 1];
    latest_[i] =
        std::min(nodes[cur].due, latest_[i + 1] - inst_->distance(cur, next) -
                                     nodes[cur].service);
  }

  segments_valid_ = false;
}

void Route::computeSegments() const {
  const std::size_t n = nodes_.size();
  forward_.resize(n);
  backward_.resize(n);
  forward_[0] = Segment::single(*inst_, nodes_[0]);
  for (std::size_t i = 1; i < n; ++i) {
    forward_[i] =
        concat(*inst_, forward_[i - 1], Segment::single(*inst_, nodes_[i]));
  }
  backward_[n - 1] = Segment::single(*inst_, nodes_[n - 1]);
  for (std::size_t i = n - 1; i-- > 0;) {
    backward_[i] =
        concat(*inst_, Segment::single(*inst_, nodes_[i]), backward_[i + 1]);
  }
  segments_valid_ = true;
}

// ----------------------------------------------------------------------------
// Solution
// ----------------------------------------------------------------------------

Solution::Solution(const Instance& inst) : inst_(&inst) {}

Route& Solution::addRoute() {
  routes_.emplace_back(*inst_);
  return routes_.back();
}

void Solution::removeEmptyRoutes() {
  routes_.erase(std::remove_if(routes_.begin(), routes_.end(),
                               [](const Route& r) { return r.empty(); }),
                routes_.end());
}

int Solution::numUsedRoutes() const {
  return static_cast<int>(
      std::count_if(routes_.begin(), routes_.end(),
                    [](const Route& r) { return !r.empty(); }));
}

bool Solution::isFeasible() const {
  return std::all_of(routes_.begin(), routes_.end(),
                     [](const Route& r) { return r.isFeasible(); });
}

double Solution::totalDistance() const {
  double total = 0.0;
  for (const Route& route : routes_) total += route.distance();
  return total;
}

std::vector<int> Solution::unroutedCustomers() const {
  std::vector<bool> routed(inst_->numNodes(), false);
  for (const Route& route : routes_) {
    for (int pos = 0; pos < route.numCustomers(); ++pos) {
      routed[route.customer(pos)] = true;
    }
  }
  std::vector<int> unrouted;
  for (int c = 1; c < inst_->numNodes(); ++c) {
    if (!routed[c]) unrouted.push_back(c);
  }
  return unrouted;
}

bool Solution::isBetterThan(const Solution& other) const {
  const int vehicles = numUsedRoutes();
  const int other_vehicles = other.numUsedRoutes();
  if (vehicles != other_vehicles) return vehicles < other_vehicles;
  return totalDistance() < other.totalDistance() - kEps;
}

void Solution::validate() const {
  const std::vector<Node>& nodes = inst_->nodes;
  std::vector<int> visits(inst_->numNodes(), 0);

  double total_distance = 0.0;
  for (std::size_t r = 0; r < routes_.size(); ++r) {
    const Route& route = routes_[r];
    double time = nodes[0].ready;
    double distance = 0.0;
    int load = 0;
    int prev = 0;

    for (int pos = 0; pos < route.numCustomers(); ++pos) {
      const int c = route.customer(pos);
      if (c < 1 || c >= inst_->numNodes()) {
        std::ostringstream msg;
        msg << "route " << r << ": invalid customer id " << c;
        throw std::runtime_error(msg.str());
      }
      ++visits[c];
      load += nodes[c].demand;
      distance += inst_->distance(prev, c);
      time = std::max(time + nodes[prev].service + inst_->distance(prev, c),
                      nodes[c].ready);
      if (time > nodes[c].due + kEps) {
        std::ostringstream msg;
        msg << "route " << r << ": customer " << c << " starts at " << time
            << ", after due time " << nodes[c].due;
        throw std::runtime_error(msg.str());
      }
      prev = c;
    }

    distance += inst_->distance(prev, 0);
    total_distance += distance;
    // The cached distance is summed in the same order, so any difference
    // beyond rounding means the cache is stale.
    if (std::fabs(route.distance() - distance) >
        kEps * std::max(1.0, distance)) {
      std::ostringstream msg;
      msg << "route " << r << ": cached distance " << route.distance()
          << " differs from recomputed " << distance;
      throw std::runtime_error(msg.str());
    }

    time += nodes[prev].service + inst_->distance(prev, 0);
    if (time > inst_->horizon() + kEps) {
      std::ostringstream msg;
      msg << "route " << r << ": returns to depot at " << time
          << ", after horizon " << inst_->horizon();
      throw std::runtime_error(msg.str());
    }
    if (load > inst_->capacity) {
      std::ostringstream msg;
      msg << "route " << r << ": load " << load << " exceeds capacity "
          << inst_->capacity;
      throw std::runtime_error(msg.str());
    }
  }

  for (int c = 1; c < inst_->numNodes(); ++c) {
    if (visits[c] != 1) {
      std::ostringstream msg;
      msg << "customer " << c << " visited " << visits[c] << " times";
      throw std::runtime_error(msg.str());
    }
  }

  if (std::fabs(totalDistance() - total_distance) >
      kEps * std::max(1.0, total_distance)) {
    std::ostringstream msg;
    msg << "total distance " << totalDistance() << " differs from recomputed "
        << total_distance;
    throw std::runtime_error(msg.str());
  }

  if (numUsedRoutes() > inst_->num_vehicles) {
    std::ostringstream msg;
    msg << numUsedRoutes() << " routes used, fleet size is "
        << inst_->num_vehicles;
    throw std::runtime_error(msg.str());
  }
}

void Solution::writeRoutes(const std::string& instance_name,
                           const std::string& date, std::ostream& os) const {
  os << "Instance name : " << instance_name << '\n'
     << "Authors       : claude_vrp\n"
     << "Date          : " << date << '\n'
     << "Reference     : claude_vrp "
        "(https://github.com/amaldesilva2/claude_vrp)\n"
     << "Solution\n";
  int index = 1;
  for (const Route& route : routes_) {
    if (route.empty()) continue;
    os << "Route " << index++ << " :";
    for (int pos = 0; pos < route.numCustomers(); ++pos) {
      os << ' ' << route.customer(pos);
    }
    os << '\n';
  }
}

void Solution::print(std::ostream& os) const {
  int index = 1;
  for (const Route& route : routes_) {
    if (route.empty()) continue;
    os << "Route " << index++ << ":";
    for (int pos = 0; pos < route.numCustomers(); ++pos) {
      os << ' ' << route.customer(pos);
    }
    os << '\n';
  }
  os << "Vehicles: " << numUsedRoutes() << '\n'
     << "Distance: " << std::fixed << std::setprecision(2) << totalDistance()
     << '\n';
}
