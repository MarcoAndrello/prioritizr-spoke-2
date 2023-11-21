# Prioritization
rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(prioritizr)

load("Planning_units.RData")
feature_names <- names(pus)[4:ncol(st_drop_geometry(pus))]

p1 <-
  problem(pus, features = feature_names, cost_column="cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = T)

s1 <- solve(p1)

plot(s1["solution_1"],border=NA)

## Prossimi passi:
# shapefile delle tree aree pilota per zoom
# climate change
# irreplaceability
# connectivity
# targets diversi
