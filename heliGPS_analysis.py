# %% 
## MPB heli-GPS Analysis

from cgi import test
import affine
import fiona
import geopandas as gpd
import georasters as gr
from locale import normalize
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterstats as rs
import rasterio
from rasterio.plot import plotting_extent
import rioxarray as rxr


F_GDB = './Data/MPB_AERIAL_SURVEY.gdb'
F_SPECIES = './Data/Species_classification_2019/Species_classification_2019_aoi'
F_AOI = './Data/Areas_of_Interest/Areas_of_Interest.shp'
F_GRAVEL_ROADS = './Data/Road_Network/_gravelRoads.shp'
F_PAVED_ROADS = './Data/Road_Network/_pavedRoads.shp'
F_BUFFER = './Data/MPB_buffer/heliGPS_buffer.shp'
F_MPB_BUFFERED = './Data/MPB_buffer/heliGPS_species_df'

CRS = 'EPSG:3400'

# %%

# geopandas included map, filtered to just this hemisphere
# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
# worlsd = world.to_crs(CRS)
# westhem = world[(world['continent'] == 'North America') | 
#                 (world['continent'] == 'South America')]
# base = westhem.plot(color='white', edgecolor='black', figsize=(11, 11))

layers = fiona.listlayers(F_GDB)
heliGPS_layers = sorted([lyr for lyr in layers if lyr.endswith('x')])
polygons = sorted([lyr for lyr in layers if lyr.endswith('p')])
heliGPS_last_3_years = heliGPS_layers[-3:]

heliGPS = None
for l in heliGPS_last_3_years:
    if heliGPS is None:
        heliGPS = gpd.read_file(F_GDB, layer=l)
    else:
        survey = gpd.read_file(F_GDB, layer=l)
        heliGPS = pd.concat([heliGPS, survey])
        
heliGPS = heliGPS.to_crs(CRS)
heliGPS = heliGPS.reset_index()
heliGPS

# imported data should be of geomerty type point. check with: heliGPS.geom_type.head()

# interactively visualize the data
# heliGPS.explore()

# %%
# First Order Statistics with the Heli-GPS Data
heliGPS.info()
heliGPS.describe(include = 'all')

print(heliGPS.dmg_desc.value_counts(ascending=True))
print(heliGPS.survyear.value_counts(ascending=True))
print(heliGPS.att_stage.value_counts(ascending=True))
print(heliGPS.zone.value_counts(ascending=True))
print(heliGPS.loc[heliGPS['att_stage'] == 'Green'].survyear.value_counts(ascending=True))

# %% 
# These data points do not have an attack stage (att_stage) assigned to them, thus they will be excluded from further anlysis.
heliGPS.loc[heliGPS['att_stage'] == ''].explore()
heliGPS.drop(heliGPS[heliGPS.att_stage  == ''].index, inplace=True)
# All remaining points are associated either to attack stage "Green" or "Red".


# Visualize some statistics of the data

# Distribution of affected trees, grouped by att_stage
ax1 = heliGPS.plot.hist(column=['num_trees'], by='att_stage', figsize=(10, 8))

# boxplot of the spacial distribution for red and green attack trees
fig, axes = plt.subplots(1, 2, figsize=(10, 8))
for i, c in enumerate(['longitude', 'latitude']):
    heliGPS.boxplot(column=c, by='att_stage', ax=axes[i])
fig.suptitle('Spacial Distribution of Red and Green Attack Trees', fontsize=16)

# %% 
# Create buffer around heli-GPS points
heliGPS_poly = heliGPS.copy()
heliGPS_poly["geometry"] = heliGPS.geometry.buffer(, cap_style = 3)
heliGPS_poly.head()
heliGPS_poly.to_file(F_BUFFER)

# %% 
# Read tree species data
species = []
gdf_species = None
for i in [1]:
    species = rxr.open_rasterio(F_SPECIES + str(i) + '.tif', masked=True).squeeze()
    species.rio.reproject(CRS)
#     species = gr.from_file(F_SPECIES + str(i) + '.tif')
#     gdf = species.to_geopandas()
#     if gdf_species is None:
#         gdf_species = gdf
#     else:
#         gdf_species = pd.concat([gdf_species, gdf])

# gdf_species.to_crs(CRS)
# # cleanup: set values of 0, corresponding to no trees, to NAN
# # gdf_species_no_zeros = gdf_species.where(gdf_species == 0, np.nan)
# ax = gdf_species.value.plot.hist(bins=max(gdf_species.value.unique()), figsize=(8,8))
print('finished loading species')
# %% 
# extract pixel values for each polygon of the bufferes heli-GPS data
# bounds = gdf_species.total_bounds
# north, west = bounds[0:2]
# xsize = bounds[2] - bounds[0]
# ysize = bounds[3] - bounds[1]
# species = rxr.open_rasterio(F_SPECIES + str(i) + '.tif', masked=True).squeeze().rio.reproject(CRS)
affine=species.rio.transform()

heliGPS_species = rs.zonal_stats(F_BUFFER,
                                    species.values, # gdf_species_no_zeros.values,
                                    nodata=0,
                                    # affine=rasterio.transform.from_origin(north, west, xsize, ysize),
                                    affine=affine,
                                    geojson_out=True,
                                    copy_properties=True,
                                    stats=['count', 'min', 'mean', 'max', 'median', 'majority'])
                                
heliGPS_species_df = gpd.GeoDataFrame.from_features(heliGPS_species)

if heliGPS_species_df.crs is None:
    heliGPS_species_df.set_crs(CRS)
else:
    heliGPS_species_df.to_crs(CRS)

heliGPS_species_df.head()

# save file to GeoJSON and shapefile
heliGPS_species_df.to_file(F_MPB_BUFFERED + '.shp')
heliGPS_species_df.to_file(F_MPB_BUFFERED + '.shp')

# %%
fig, ax = plt.subplots(figsize=(10, 10))

species.plot(ax=ax)
heliGPS.plot(ax=ax, marker='o', markersize=2, color='red')
# ax.set_axis_off()
plt.show()

# %%
