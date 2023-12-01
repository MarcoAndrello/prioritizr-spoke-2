# Simulate connectivity between PUs
rm(list=ls())

library(tidyverse)
library(sf)
library(RNetCDF)

load("Planning_units_SoS.RData")

pus %>%
  st_geometry %>%
  st_centroid %>%
  st_transform(st_crs(4326)) %>%
  st_coordinates %>%
  replicate(100, ., simplify = FALSE) %>%
  do.call("rbind", .) -> release_points

write.table(release_points, file=paste0(getwd(),"/Release_points_SoS.txt"),
            row.names=F, col.names=F, quote=F)

# Run larval dispersal in Ichtyop with MEDSEA
curr <- open.nc(paste0(getwd(),"/ichthyop/AllMed_currents_2021_June.nc"),write=T)
print.nc(curr)
dim.rename.nc(curr, "lon", "longitude")
dim.rename.nc(curr, "lat", "latitude")
close.nc(curr)
rm(curr)
