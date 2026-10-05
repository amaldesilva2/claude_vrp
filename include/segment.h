#ifndef CLAUDE_VRP_INCLUDE_SEGMENT_H_
#define CLAUDE_VRP_INCLUDE_SEGMENT_H_

#include <algorithm>

#include "instance.h"

// Summary of a sequence of consecutive visits, following Vidal et al. (2013),
// "A hybrid genetic algorithm with adaptive diversity management for a large
// class of vehicle routing problems with time windows". Two segments can be
// concatenated in O(1), which makes any move that rearranges whole pieces of
// routes checkable in O(number of pieces).
//
// A route (depot ... depot) is time-feasible iff its time_warp is zero.
struct Segment {
  int first = 0;           // first node of the sequence
  int last = 0;            // last node of the sequence
  double duration = 0.0;   // minimal travel + service + waiting time
  double time_warp = 0.0;  // total lateness forced by the due times
  double earliest = 0.0;   // earliest service start at `first`
  double latest = 0.0;     // latest service start at `first` without warp
  double distance = 0.0;   // travel distance inside the sequence
  int load = 0;            // total demand

  static Segment single(const Instance& inst, int node) {
    const Node& n = inst.nodes[node];
    Segment s;
    s.first = node;
    s.last = node;
    s.duration = n.service;
    s.earliest = n.ready;
    s.latest = n.due;
    s.load = n.demand;
    return s;
  }
};

// Summary of `a` followed directly by `b`.
inline Segment concat(const Instance& inst, const Segment& a,
                      const Segment& b) {
  const double travel = inst.distance(a.last, b.first);
  const double delta = a.duration - a.time_warp + travel;
  const double wait = std::max(b.earliest - delta - a.latest, 0.0);
  const double warp = std::max(a.earliest + delta - b.latest, 0.0);

  Segment s;
  s.first = a.first;
  s.last = b.last;
  s.duration = a.duration + b.duration + travel + wait;
  s.time_warp = a.time_warp + b.time_warp + warp;
  s.earliest = std::max(b.earliest - delta, a.earliest) - wait;
  s.latest = std::min(b.latest - delta, a.latest) + warp;
  s.distance = a.distance + b.distance + travel;
  s.load = a.load + b.load;
  return s;
}

#endif  // CLAUDE_VRP_INCLUDE_SEGMENT_H_
