library(sf)

subregions <- st_read('./Data/Aria_Guo/summaries/Natural_Regions_Subregions_of_Alberta.shp')
alberta <- st_read('./Data/alberta_outline/alberta_outline.shp')


alberta_pro <- st_transform(alberta, st_crs(subregions))
plot(st_geometry(alberta_pro))
plot(st_geometry(subregions), add=TRUE)

st_write(alberta_pro, './Data/alberta_outline/alberta_pro.shp')