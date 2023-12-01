# 01 - Create planning units

rm(list=ls())

library(terra)
library(sf)
library(tidyverse)
library(units)

# Read shapefile of pilot area
pa <- st_read(paste0(getwd(),"/data/SoS_pilot_area.shp"))

# Make planning unit geometries
st_make_grid(st_bbox(pa),cellsize=5000) %>% # rectangular grid
  vect %>% crop(vect(pa)) %>% # convert to SpatVect to use terra::crop
  st_as_sf -> pus # convert back to sf
plot(pus)

# Caluclate area
pus$AREA <- st_area(pus)
plot(pus["AREA"])


# Read ETOPO1 depth data
etopo1 <- rast("C:/Users/marco/Il mio Drive/Maps/ETOPO1-ice-Mediterranean.tif")
# Extract depth from ETOPO1
pus %>% st_centroid() %>% st_transform(crs=st_crs(4326)) %>% st_coordinates() -> pu_coord_4326
terra::extract(etopo1, pu_coord_4326) -> depth
# Add depth to pus
pus$depth <- depth$`ETOPO1-ice-Mediterranean`
# remove land cells
pus %>% filter(depth < 0) -> pus_temp
# st_write(pus_temp, dsn="pus_temp.shp")
# # In QGis, delete cells out of domain
# pus <- st_read("pus_temp.shp")
pus$ID <- 1:nrow(pus)
rm(etopo1, pu_coord_4326, pus_temp, depth)




################################################################################
# Import conservation costs
################################################################################

# Read shapefile and values
# The orginal files are available upon request from Tessa Mazor: t.mazor@uq.edu.au or Sylvaine Giakoumi: sylvaine.giakoumi@szn.it
grid <- st_read(paste0(getwd(),"/../Genetic SCP Diplodus Mullus/data/Mazor et al 2014 Opportunity costs/GridMEDSEA_Clip_etrs89.shp"))
values <- read.csv(paste0(getwd(),"/../Genetic SCP Diplodus Mullus/data/Mazor et al 2014 Opportunity costs/CombinedCostSenario7.csv"),sep=";")

# Join values to the sf
grid %>% 
    left_join(values,by="id") %>% 
    select(id,cost) %>% 
    st_transform(st_crs(pus)) -> 
    grid

# Intersect the two datasets and calculate the area of each intersection
st_intersection(pus, grid) %>%
    mutate(AREA = st_area(.)) %>%
    select(ID, AREA, cost) -> 
    intersection

# Calculate the mean cost as weighted means, the weights are the areas
intersection %>%
    mutate(AREA_NUM = as.numeric(AREA)) %>%
    group_by(ID) %>%
    st_drop_geometry() %>% 
    summarise(mean_cost = weighted.mean(cost,AREA_NUM)) ->
    mean_cost

# Join the mean costs to the PUs
left_join(pus, mean_cost, by = "ID") -> pus
pus %>% mutate(cost=mean_cost) %>% select(-mean_cost) -> pus
# Rescale mean cost to avoid numerical issues
pus$cost <- pus$cost / 1000

# Calculate coordinates of the centroid of each PU
pus %>%
    st_geometry() %>%
    st_centroid() %>%
    st_coordinates() -> pus_centroid



################################################################################
# Extract species presence from FishMed
################################################################################
# Read Fishmed data
FishMed_grid <- read.csv(paste0(getwd(),"/../Genetic SCP Diplodus Mullus/data/Albouy et al 2015 FishMed/Observed_grid_1980.csv"),h=T,sep=";")
FishMed_grid <- rast(FishMed_grid, type="xyz", crs="EPSG:4326")
pus %>% st_geometry %>% st_centroid %>% st_transform(crs=st_crs(4326)) %>% st_coordinates() %>% terra::extract(FishMed_grid, .) -> species_presence
pus <- cbind(pus, species_presence)

feature_names <- names(FishMed_grid)

# Descriptive statistics
# Add species richness
pus %>%
  st_drop_geometry %>%
  select(feature_names) %>%
  rowSums %>%
  mutate(pus, richness=.) -> pus
plot(pus["richness"])

# Calculate area (number of cells) of local spatial distribution for each species
pus %>%
  st_drop_geometry %>%
  select(feature_names) %>%
  colSums(na.rm=T) ->
local_area
hist(local_area)
hist(local_area/nrow(pus)) # few species spread over all the pilot area, but most have a restricted distribution 

# Calculate area (number of cells) of spatial distribution in the whole Mediterranean for each species
global(FishMed_grid, "sum", na.rm=T) %>% pull(sum) -> med_area
hist(med_area)
hist(med_area / ncell(FishMed_grid))

# Save
save(pus, pus_centroid,
     feature_names,
     local_area, med_area,
     file="Planning_units_SoS.RData")

