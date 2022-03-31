library(sf)
library(terra)
library(ggplot2)

F_RASTER <- vector(mode = "list", length = 3)
for (i in 1:3) {
    F_RASTER[[i]] <- paste(c("./Data/Species_classification_2019/Species_classification_2019_aoi", i, ".tif"), sep="", collapse="")

    # species_aoi <- rast(F_RASTER
}
fx <- lapply(F_RASTER, rast)

summaries <- lapply(1:3, function(i) {summary(fx[[i]])})
var <- lapply(1:3, function(i) {global(fx[[i]], sd, na.rm=TRUE)})
tab_species <- lapply(1:3, function(i) {table(fx[[i]])})

hist(fx[[1]])
df_s <- df(fx)
