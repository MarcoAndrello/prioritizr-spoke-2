# Prioritization
rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(prioritizr)
library(tmap)
library(rnaturalearth)

# p1: with cost from Mazor
# p3: No cost + boundary penalty
# p4: No cost + connectivity penalty
# p5: No cost + representation targets proportional to species distributon range
# p6: like p1 but maximal coverage
# p7: like p1 but locking in existing PAs


load("Planning_units_SoS.RData")
feature_names

# Limit to conservation features present in the pilot area
id.features.present <- which(local_area > 0)
features_present <- feature_names[id.features.present] # 398 species

# Maps
ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(pus)) -> countries
## Planning units
tm_shape(pus) +
  tm_borders(col="black") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray")
## Costi di conservazione
tm_shape(pus) +
  tm_fill(col="cost", title="Cost (k€)") +
  tm_legend(legend.outside = T, legend.outside.position = "right") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray")
## Biodiversity features
pus_richness <- pus
pus %>% st_drop_geometry %>% select(all_of(features_present)) %>% rowSums %>% mutate(pus_richness, richness=.) -> pus_richness
tm_shape(pus_richness) +
  tm_fill(col="richness", title="Number of species", palette="Greens") +
  tm_legend(legend.outside = T, legend.outside.position = "right") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray")
rm(pus_richness)

# p1: Cost from Mazor et al, and 20% target
p1 <-
  problem(pus, features = features_present, cost_column="cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = F)
s1 <- solve(p1)
s1$sol_1 <- factor(s1$solution_1)
tm_shape(s1) +
  tm_fill(col="sol_1", title="Solution") +
  tm_legend(legend.outside = T, legend.outside.position = "right") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray")

eval_n_summary(p1,s1["solution_1"])
eval_cost_summary(p1,s1["solution_1"])
eval_boundary_summary(p1,s1["solution_1"])
eval_target_coverage_summary(p1,s1["solution_1"])$met %>% which %>% length


# p2: Cost from Mazor et al, and variable targets
# Approccio di Rodrigues et al 2004 Bioscience
# Each cell of FishMed grid is 0.1 * 0.1 degree so approx 11 km * 8.5 km ~= 90 km2
absolute_med_area <- med_area * 90
hist(absolute_med_area)
# Imposta target (valori arbitrari a intervalli arbitrari)
relative_targets<- c(0.5, 0.5, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# Grafico area di distribuzzione e target di rappresentazione
plot(cbind(seq(250000,2500000,250000),relative_targets),type="l", xlab="Distribution range (km^2)", ylab="Representation target")
# Assegna a ogni specie un target sulla base dell'area di distribuzione
relative_targets_species <- relative_targets[cut(absolute_med_area,breaks=seq(0,2500000,250000))]
# Distribuzione dei targets
table(relative_targets_species) # su tutte le 635 specie
relative_targets_species_present <- relative_targets_species[id.features.present]
table(relative_targets_species_present) # sulle 398 specie presenti nell'area
# Calcolo targets in termini assoluti (numero di PUs)
absolute_targets_species <- relative_targets_species * local_area
absolute_targets_species_present <- absolute_targets_species[id.features.present]
## there was a publication for setting targets
## https://www.sciencedirect.com/science/article/pii/S0006320723003191

p2 <-
  problem(pus, features = features_present, cost_column="cost") %>%
  add_min_set_objective() %>%
  add_absolute_targets(absolute_targets_species_present) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = T)
s2 <- solve(p2)
s2$sol_1 <- factor(s2$solution_1)
tm_shape(s2) +
  tm_fill(col="sol_1", title="Solution") +
  tm_legend(legend.outside = T, legend.outside.position = "right") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray")

# Replacement importance
rc1 <- eval_replacement_importance(p1, s1["solution_1"],run_checks=F)
tm_shape(rc1) +
  tm_fill(col="rc", title="Irreplaceability") +
  tm_legend(legend.outside = T, legend.outside.position = "right") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray")


# p3: Boundary penalty
p3 <-
  problem(pus, features = features_present, cost_column="cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_boundary_penalties(0.1) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.25, verbose = T)
s3 <- solve(p3)
s3$sol_1 <- factor(s3$solution_1)
##Plot map
tm_shape(s3) +
  tm_fill(col="sol_1", title="Solution") +
  tm_legend(legend.outside = T, legend.outside.position = "right") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray")
## Summary statistics
eval_n_summary(p3,s3["solution_1"])
eval_cost_summary(p3,s3["solution_1"])
eval_boundary_summary(p3,s3["solution_1"])
#


# p4: Maximal feature
p4 <-
  problem(pus, features = features_present, cost_column="cost") %>%
  add_max_features_objective(1000000) %>%
  add_relative_targets(0.2) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = T)
s4 <- solve(p4)
s4$sol_1 <- factor(s4$solution_1)
##Plot map
tm_shape(s4) +
  tm_fill(col="sol_1", title="Solution") +
  tm_legend(legend.outside = T, legend.outside.position = "right") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray")
## Summary statistics
eval_n_summary(p4,s4["solution_1"])
eval_cost_summary(p4,s4["solution_1"])
eval_target_coverage_summary(p4,s4["solution_1"])$met %>% which %>% length
#


# p5: Locked-in
## Download MPAs with WDPAr
library(wdpar)
pa <- wdpa_fetch("ITA")
pa %>% filter(MARINE == 2) %>% st_transform(st_crs(pus)) %>% st_crop(pus) -> mpas_sos
a <- st_intersects(pus, mpas_sos)
pus$is.protected <- 0
pus$is.protected[which(lengths(a)>0)] <- 1
plot(pus["is.protected"])

p5 <-
  problem(pus, features = features_present, cost_column="cost") %>%
  add_min_set_objective() %>%
  add_relative_targets(0.2) %>%
  add_locked_in_constraints(which(lengths(a)>0)) %>%
  add_binary_decisions() %>%
  add_gurobi_solver(gap = 0.1, verbose = F)
s5 <- solve(p5)
s5$sol_1 <- factor(s5$solution_1)
tm_shape(s5) +
  tm_fill(col="sol_1", title="Solution") +
  tm_legend(legend.outside = T, legend.outside.position = "right") +
  tm_shape(countries, bbox = res) +
  tm_polygons(col="lightgray") +
  tm_shape(mpas_sos) + 
  tm_borders(col="red",lwd=2)
## Summary statistics
eval_n_summary(p5,s5["solution_1"])
eval_cost_summary(p5,s5["solution_1"])
eval_target_coverage_summary(p5,s5["solution_1"])$met %>% which %>% length
#


save(p1, s1,
     p2, s2,
     p3, s3,
     p4, s4,
     p5, s5,
     file="Solutions_SoS_temp.RData")

