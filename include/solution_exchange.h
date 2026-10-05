#ifndef CLAUDE_VRP_INCLUDE_SOLUTION_EXCHANGE_H_
#define CLAUDE_VRP_INCLUDE_SOLUTION_EXCHANGE_H_

#include <mutex>
#include <optional>

#include "solution.h"

// Thread-safe board where parallel searches share their best solutions.
class SolutionExchange {
 public:
  // Records `solution` if it is the best offered so far and returns a copy
  // of the best solution offered by any thread (possibly `solution`).
  Solution offer(const Solution& solution) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!best_ || solution.isBetterThan(*best_)) best_ = solution;
    return *best_;
  }

 private:
  std::mutex mutex_;
  std::optional<Solution> best_;
};

#endif  // CLAUDE_VRP_INCLUDE_SOLUTION_EXCHANGE_H_
