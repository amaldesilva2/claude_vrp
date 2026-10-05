#ifndef CLAUDE_VRP_INCLUDE_GENETIC_SEARCH_H_
#define CLAUDE_VRP_INCLUDE_GENETIC_SEARCH_H_

#include <chrono>

#include "solution.h"

struct GeneticSearchParams {
  int min_population = 25;   // mu: sub-population size after survivor selection
  int generation_size = 40;  // lambda: offspring before survivor selection
  int num_elite = 4;         // individuals protected from diversity pressure
  int num_close = 5;         // neighbors used for the diversity contribution
  int initial_population = 20;
  double initial_removal = 0.2;  // share of customers re-inserted per initial
                                 // individual
  double repair_probability = 0.5;
  // Penalty adaptation: every adapt_period offspring, move each penalty
  // towards target_feasible (fraction of offspring satisfying it after
  // local search).
  int adapt_period = 100;
  double target_feasible = 0.3;
  int num_neighbors = 20;  // local search granularity
  unsigned seed = 1;
};

struct GeneticSearchStats {
  int generations = 0;   // offspring created
  int improvements = 0;  // new best feasible solutions
  int feasible_offspring = 0;
  int repaired = 0;
  double final_tw_penalty = 0.0;
  double final_load_penalty = 0.0;
};

// Hybrid genetic search for the VRPTW (Vidal et al. 2013; time-window
// version as in Kool et al. 2022 / PyVRP), with a fixed fleet: it never uses
// more routes than `solution` has.
//
// Feasible and infeasible sub-populations are ranked by biased fitness
// (penalized cost rank + diversity rank, diversity = mean broken-pairs
// distance to the num_close closest individuals). Parents are chosen by
// binary tournaments, combined with SREX (selective route exchange: a block
// of routes of one parent replaces overlapping routes of the other, missing
// customers are inserted at the cheapest penalized position) and educated
// by local search with soft time-window and capacity constraints. Infeasible
// offspring are repaired with stronger penalties with repair_probability.
// `solution` must be feasible; it receives the best feasible solution found.
GeneticSearchStats geneticSearch(const GeneticSearchParams& params,
                                 std::chrono::steady_clock::time_point deadline,
                                 Solution* solution);

#endif  // CLAUDE_VRP_INCLUDE_GENETIC_SEARCH_H_
