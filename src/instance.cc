#include "instance.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace {

// Advances the stream until `keyword` has been read. Throws if EOF is hit
// first.
void skipTo(std::istream& in, const std::string& keyword,
            const std::string& path) {
  std::string token;
  while (in >> token) {
    if (token == keyword) return;
  }
  throw std::runtime_error(path + ": expected '" + keyword + "'");
}

}  // namespace

double computeDistance(const Node& a, const Node& b) {
  return std::hypot(a.x - b.x, a.y - b.y);
}

int capacityLowerBound(const Instance& inst) {
  int64_t total_demand = 0;
  for (int c = 1; c < inst.numNodes(); ++c)
    total_demand += inst.nodes[c].demand;
  return static_cast<int>((total_demand + inst.capacity - 1) / inst.capacity);
}

std::vector<std::vector<int>> nearestCustomers(const Instance& inst, int k) {
  const int n = inst.numNodes();
  k = std::max(0, std::min(k, inst.numCustomers() - 1));
  std::vector<std::vector<int>> nearest(n);
  std::vector<int> others;
  for (int u = 1; u < n; ++u) {
    others.clear();
    for (int v = 1; v < n; ++v) {
      if (v != u) others.push_back(v);
    }
    std::partial_sort(others.begin(), others.begin() + k, others.end(),
                      [&](int a, int b) {
                        return inst.distance(u, a) < inst.distance(u, b);
                      });
    nearest[u].assign(others.begin(), others.begin() + k);
  }
  return nearest;
}

Instance loadInstance(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);

  Instance inst;

  // Token-based parsing: operator>> skips spaces, '\r' and blank lines,
  // so CRLF files and whitespace-only lines need no special handling.
  if (!(in >> inst.name)) throw std::runtime_error(path + ": empty file");

  skipTo(in, "VEHICLE", path);
  skipTo(in, "CAPACITY", path);  // header: NUMBER CAPACITY
  if (!(in >> inst.num_vehicles >> inst.capacity))
    throw std::runtime_error(path + ": bad vehicle section");

  skipTo(in, "CUSTOMER", path);
  // Header: CUST NO. XCOORD. YCOORD. DEMAND READY TIME DUE DATE SERVICE TIME
  skipTo(in, "TIME", path);  // READY TIME
  skipTo(in, "TIME", path);  // SERVICE TIME

  Node node;
  while (in >> node.id >> node.x >> node.y >> node.demand >> node.ready >>
         node.due >> node.service) {
    if (node.id != inst.numNodes())
      throw std::runtime_error(path + ": unexpected customer id " +
                               std::to_string(node.id) + " (expected " +
                               std::to_string(inst.numNodes()) + ")");
    inst.nodes.push_back(node);
  }
  if (!in.eof()) throw std::runtime_error(path + ": malformed customer line");
  if (inst.nodes.size() < 2) throw std::runtime_error(path + ": no customers");

  const std::size_t n = inst.nodes.size();
  inst.dist.resize(n * n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      inst.dist[i * n + j] = computeDistance(inst.nodes[i], inst.nodes[j]);

  return inst;
}
