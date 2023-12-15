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
features_present <- feature_names[which(local_area > 0)] # 398 species

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
# 
# # p4: No cost + connectivity penalty
# p4 <-
#   problem(pus, features = features_present, cost_column="unitary_cost") %>%
#   add_min_set_objective() %>%
#   add_relative_targets(0.2) %>%
#   add_boundary_penalties(0.001) %>%
#   add_binary_decisions() %>%
#   add_gurobi_solver(gap = 0.6, verbose = T)
# s4 <- solve(p4)
# plot(s4["solution_1"],border=NA,  main="p4: No cost + connectivity penalty")

save(p1, s1,
     p2, s2,
     p3, s3,
     file="Solutions_SoS_temp.RData")

# p5: No cost + representation targets proportional to species distribution range
# Rodrigues et al 2004 Bioscience
# More demanding representation targets (a larger percentage of the range)
# were set for species with more restricted ranges 
relative_local_area <- local_area/nrow(pus) ## 635 species
relative_targets<- c(1, 1, 1, 1, 0.9, 0.7, 0.5, 0.3, 0.1, 0.1)
hist(relative_local_area)
plot(cbind(seq(0.1,1,0.1),relative_targets),type="l", xlab="Distribution range (relative)", xlab="Representation target (relative)")
relative_targets_species <- relative_targets[cut(relative_local_area,breaks=seq(0,1,0.1))]
relative_targets_species <- relative_targets_species[which(local_area > 0)]

p5l <-
  problem(pus, features = features_present, cost_column="unitary_cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(relative_targets_species) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = T)
s5l <- solve(p5l)
eval_target_coverage_summary(p5l,s5l["solution_1"]) %>% print(n=400)


relative_global_area <- med_area/27095 ## 635 species; 27095 is approx number of sea cells
# Each cell of FishMed grid is 0.1 * 0.1 degree so approx 11 km * 6 km = 70 km2
absolute_global_area = 70 * med_area 
hist(absolute_global_area)
### now find range area to target relationship 
## there was a publication for setting targets
# https://www.sciencedirect.com/science/article/pii/S0006320723003191

plot(cbind(seq(0.1,1,0.1),relative_targets),type="l", xlab="Distribution range (relative)", ylab="Representation target (relative)")
relative_targets_species <- relative_targets[cut(relative_global_area,breaks=seq(0,1,0.1))]
relative_targets_species <- relative_targets_species[which(local_area > 0)]

p5g <-
  problem(pus, features = features_present, cost_column="unitary_cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(relative_targets_species) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = T)
s5g <- solve(p5g)
eval_target_coverage_summary(p5g,s5g["solution_1"]) %>% print(n=400)

save(p1, s1,
     p2, s2,
     p3, s3,
     p5l, s5l,
     file="Solutions_SoS_temp.RData")

###

# p6: like p1 but maximal coverage



# p7: like p1 but locking in existing PAs