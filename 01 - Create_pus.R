# 01 - Create planning units

rm(list=ls())

library(terra)
library(sf)
library(tidyverse)
library(units)

# Define extent of planning region (domain) from long-lat coordinates in WGS84 CRS
region <- st_sfc(
    st_polygon(
        list(rbind(c(-5.5,29), c(36.5,29), c(36.5,45), c(-5.5,45), c(-5.5,29)))
    ),
    crs=st_crs(4326))

# Transform to Albers Equal Area CRS - to get planning units of the same area (otherwise southern PUs will be larger in long-lat CRS)
region <- st_transform(region, crs=st_crs("ESRI:102013"))
# EPSG:19986 Lambers equal area recommended by the EU

# Make planning units
st_make_grid(region,cellsize=10000) %>% st_sf() -> pus
rm(region)

# Read ETOPO1 depth data
etopo1 <- rast("C:/Users/marco/Il mio Drive/Maps/ETOPO1-ice-Mediterranean.tif")
# Extract depth from ETOPO1
pus %>% st_centroid() %>% st_transform(crs=st_crs(4326)) %>% st_coordinates() -> pu_coord_4326
terra::extract(etopo1, pu_coord_4326) -> depth
# Add depth to pus
pus$depth <- depth$`ETOPO1-ice-Mediterranean`
# remove land cells
pus %>% filter(depth < 0) -> pus_temp 
st_write(pus_temp, dsn="pus_temp.shp")
# In QGis, delete cells out of domain
pus <- st_read("pus_temp.shp")
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
pus %>% st_geometry %>% st_centroid %>% st_transform(crs=st_crs(4326)) %>% st_coordinates() %>% terra::extract(FishMed_grid, .) -> whole_Med_presence
pus <- cbind(pus, whole_Med_presence)

# Save
save(pus,pus_centroid,file="Planning_units.RData")
