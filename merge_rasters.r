library(sf)
library(terra)

F_IN_1 <- "./Data/Aria_Guo/summaries/z11/elevation_percentage_first_returns_above_1pnt37.tif"
F_IN_2 <- "./Data/Aria_Guo/summaries/z12/elevation_percentage_first_returns_above_1pnt37.tif"
F_OUT <- "./Data/Aria_Guo/combined/elevation_percentage_first_returns_above_1pnt37.tif"

raster_1 <- rast(F_IN_1)
raster_2 <- rast(F_IN_2)
raster_2 <- project(raster_2, raster_1)
rlist <- list(raster_1, raster_2)
rsrc <- sprc(rlist)

merge(rsrc, filename=F_OUT, overwrite=TRUE)

# for (i in 1:3){
#     F_SHAPE <- paste(c("./Data/Areas_of_Interest/AoI", i, ".shp"), sep="", collapse="")
#     F_OUT <- paste(c("./Data/Species_classification_2019/Species_classification_2019_aoi", i, ".tif"), sep="", collapse="")

#     AoI <- vect(F_SHAPE)
#     species <- rast(F_RASTER)
#     species_aoi <- crop(species, AoI, snap = "near", mask=TRUE, filename = F_OUT, overwrite=TRUE)

#     plot(species_aoi)
# }
