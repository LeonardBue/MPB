library(sf)
library(terra)

F_RASTER <- "./Data/Species_classification_2019/Species_classification_2019_Alberta.tif"


for (i in 1:3){
    F_SHAPE <- paste(c("./Data/Areas_of_Interest/AoI", i, ".shp"), sep="", collapse="")
    F_OUT <- paste(c("./Data/Species_classification_2019/Species_classification_2019_aoi", i, ".tif"), sep="", collapse="")

    AoI <- vect(F_SHAPE)
    species <- rast(F_RASTER)
    species_aoi <- crop(species, AoI, snap = "near", mask=TRUE, filename = F_OUT, overwrite=TRUE)

    plot(species_aoi)
}


# plot(species)

# plot(species_aoi)

# # read two shape files, transform into another crs and write a new file
# subregions <- st_read('./Data/Aria_Guo/summaries/Natural_Regions_Subregions_of_Alberta.shp')
# alberta <- st_read('./Data/alberta_outline/alberta_outline.shp')

# alberta_pro <- st_transform(alberta, st_crs(subregions))
# plot(st_geometry(alberta_pro))
# plot(st_geometry(subregions), add=TRUE)

# st_write(alberta_pro, './Data/alberta_outline/alberta_pro.shp')