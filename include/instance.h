#ifndef CLAUDE_VRP_INCLUDE_INSTANCE_H_
#define CLAUDE_VRP_INCLUDE_INSTANCE_H_

#include <cstddef>
#include <string>
#include <vector>

// A depot or customer from a Solomon / Gehring & Homberger VRPTW file.
struct Node {
  int id = 0;
  double x = 0.0;
  double y = 0.0;
  int demand = 0;
  double ready = 0.0;    // earliest service start
  double due = 0.0;      // latest service start
  double service = 0.0;  // service duration
};

struct Instance {
  std::string name;
  int num_vehicles = 0;
  int capacity = 0;
  std::vector<Node> nodes;   // nodes[0] = depot, nodes[1..n] = customers
  std::vector<double> dist;  // flat (n+1) x (n+1) distance matrix

  int numNodes() const { return static_cast<int>(nodes.size()); }
  int numCustomers() const { return numNodes() - 1; }
  double horizon() const { return nodes[0].due; }
  double distance(int i, int j) const {
    return dist[static_cast<std::size_t>(i) * nodes.size() + j];
  }
};

// Euclidean distance, plain double (no rounding/truncation).
double computeDistance(const Node& a, const Node& b);

// Minimum number of vehicles needed for the total demand: ceil(sum q / Q).
int capacityLowerBound(const Instance& inst);

// For every node, its `k` nearest customers (excluding itself and the depot),
// closest first. Entry 0 (the depot) is left empty.
std::vector<std::vector<int>> nearestCustomers(const Instance& inst, int k);

// Parses a Solomon / G&H formatted instance file. Throws std::runtime_error on
// bad input.
Instance loadInstance(const std::string& path);

#endif  // CLAUDE_VRP_INCLUDE_INSTANCE_H_
