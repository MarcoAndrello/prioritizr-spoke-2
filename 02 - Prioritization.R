# Prioritization
rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(prioritizr)

# p1: No cost
# p2: with cost from Mazor
# p3: No cost + boundary penalty
# p4: No cost + connectivity penalty
# p5: No cost + representation targets proportional to species distributon range
# p6: like p1 but maximal coverage
# p7: like p1 but locking in existing PAs


load("Planning_units_SoS.RData")
feature_names
pus$unitary_cost <- 1

# Limit to conservation features present in the pilot area
features_present <- feature_names[which(local_area > 0)]

# p1: No cost
p1 <-
  problem(pus, features = features_present, cost_column="unitary_cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = T)
s1 <- solve(p1)
plot(s1["solution_1"],border=NA, main="p1: No cost")

# p2: with cost from Mazor
p2 <-
  problem(pus, features = features_present, cost_column="cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = T)
s2 <- solve(p2)
plot(s2["solution_1"],border=NA, main="p2: using cost from Mazor et al 2013")

# p3: No cost + boundary penalty
p3 <-
  problem(pus, features = features_present, cost_column="unitary_cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_boundary_penalties(0.0001) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.25, verbose = T)
s3 <- solve(p3)
plot(s3["solution_1"],border=NA,  main="p3: No cost + boundary penalties")
eval_target_coverage_summary(p3,s3["solution_1"]) %>% print(n=400)

# p4: No cost + connectivity penalty
p4 <-
  problem(pus, features = features_present, cost_column="unitary_cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_boundary_penalties(0.001) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.6, verbose = T)
s4 <- solve(p4)
plot(s4["solution_1"],border=NA,  main="p4: No cost + connectivity penalty")



# p5: No cost + representation targets proportional to species distributon range
# p6: like p1 but maximal coverage
# p7: like p1 but locking in existing PAs