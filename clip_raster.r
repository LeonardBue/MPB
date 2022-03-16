library(sf)
library(terra)


SHAPEFILE <- "./Data/Areas_of_Interest/Areas_of_Interest.shp"
RASTERFILE <- "./Data/Species_classification_2019/Species_classification_2019_Alberta.tif"
OUPUTFILE <- "./Data/Species_classification_2019/Species_classification_2019_aoi.tif"

AoIs <- vect(SHAPEFILE)
species <- rast(RASTERFILE)

species_aoi <- crop(species, AoIs, snap = "near", mask=TRUE, filename = OUPUTFILE, overwrite=TRUE)

plot(AoIs)
plot(species_aoi)

# # read two shape files, transform into another crs and write a new file
# subregions <- st_read('./Data/Aria_Guo/summaries/Natural_Regions_Subregions_of_Alberta.shp')
# alberta <- st_read('./Data/alberta_outline/alberta_outline.shp')

# alberta_pro <- st_transform(alberta, st_crs(subregions))
# plot(st_geometry(alberta_pro))
# plot(st_geometry(subregions), add=TRUE)

# st_write(alberta_pro, './Data/alberta_outline/alberta_pro.shp')