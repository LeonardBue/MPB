# -----------------------------------------------------------------------------------------------------------------
# Reading Arcgis REST servers
# -----------------------------------------------------------------------------------------------------------------



#####
# Pre-Processing
# Load road layers from Alberta government into memory #
# Save to disk as .shp
#####

from geopandas import GeoDataFrame
from arcgis.features import FeatureLayer
from arcgis.features import GeoAccessor, GeoSeriesAccessor

URL1 = 'https://geospatial.alberta.ca/titan/rest/services/utility/access/MapServer/14' # paved roads at 1:20,000
URL2 = 'https://geospatial.alberta.ca/titan/rest/services/utility/access/MapServer/23' # gravel roads at 1:20,000
URL3 = 'https://geospatial.alberta.ca/titan/rest/services/utility/access/MapServer/26' # other roads at 1:20,000


def downloadFromRest(url, inputCrs, outputCrs, outputPath):
    layer = FeatureLayer(url)
    features = layer.query() # Get Features
    # FeatureLayer --> Pandas --> GeoPandas
    features_gdf = GeoDataFrame(features.sdf[["FEATURE_CODE","SHAPE"]], geometry="SHAPE")
    features_gdf = features_gdf.set_crs(inputCrs) # Set CRS
    features_gdf = features_gdf.to_crs(outputCrs) # Reproject
    features_gdf.to_file(outputPath)  # Write to file

downloadFromRest(URL3,"EPSG:3400","EPSG:26911","./Data/Road_Network/_otherRoads.shp")