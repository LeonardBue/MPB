# %% 
## MPB heli-GPS Analysis

import affine
import fiona
import geopandas as gpd
import georasters as gr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import re
import rioxarray as rxr
from rasterio.plot import show
from osgeo import gdal
from rasterstats import zonal_stats
from shapely.geometry import box, Point


F_GDB = './Data/MPB_AERIAL_SURVEY.gdb'
F_SPECIES = './Data/Species_classification_2019/Species_classification_2019_aoi'
F_AOI = './Data/Areas_of_Interest/AoI'
F_GRAVEL_ROADS = './Data/Road_Network/_gravelRoads.shp'
F_PAVED_ROADS = './Data/Road_Network/_pavedRoads.shp'
F_BUFFER = './Data/MPB_buffer/heliGPS_buffer.shp'
F_MPB_BUFFERED = './Data/MPB_buffer/heliGPS_species'
F_CMAP_SPECIES = './Data/Species_classification_2019/updated_species_list_alberta.csv'

CRS = 'EPSG:3400'

POINTBUFFER = 2

# %% 
# load heli-GPS data

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
heliGPS = heliGPS.reset_index(drop=False)
heliGPS # imported data should be of geometry type point. check with: heliGPS.geom_type.head()

# These data points do not have an attack stage (att_stage) assigned to them, thus they will be excluded from further anlysis.
heliGPS.loc[heliGPS['att_stage'] == ''].explore()
heliGPS.drop(heliGPS[heliGPS.att_stage  == ''].index, inplace=True)
# All remaining points are associated either to attack stage "Green" or "Red".

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
# Visualize some statistics of the data

# Distribution of affected trees, grouped by att_stage
ax1 = heliGPS.plot.hist(column=['num_trees'], by='att_stage', figsize=(10, 8))

# boxplot of the spacial distribution for red and green attack trees
fig, axes = plt.subplots(1, 2, figsize=(10, 8))
for i, c in enumerate(['longitude', 'latitude']):
    heliGPS.boxplot(column=c, by='att_stage', ax=axes[i])
fig.suptitle('Spacial Distribution of Red and Green Attack Trees', fontsize=16)


# %%
# Extract species of trees for points in area of interest
# add a column for the area of interest; defaults to 0 --> not in AoI
heliGPS['AoI'] = np.zeros(max(heliGPS.count()), dtype=int)
heliGPS['species'] = np.zeros(max(heliGPS.count()))

# iterate through AoI number
for i in [1, 2, 3]:
    with rasterio.open(F_SPECIES + str(i) + '.tif') as species:
        bounds = species.bounds # get bounding box raster
        boundsGdf = gpd.GeoDataFrame({"geometry":[box(*bounds).buffer(-60)]}, crs=species.crs) # to GeoDataFrame

    points_within = heliGPS.to_crs(species.crs)
    points_within = heliGPS.overlay(boundsGdf, how='intersection') # Keep points within raster bounds -  https://geopandas.org/en/stable/docs/user_guide/set_operations.html
    points_within.set_index('index', inplace=True)

    heliGPS.loc[points_within.index.astype('uint64'), 'AoI'] = i # set AoI identifier
    # points_within = points_within.to_crs(CRS)
    points_within['geometry'] = points_within.geometry.buffer(distance=60, resolution=1, cap_style = 3) # Do some buffering
    zs_species = zonal_stats(points_within, F_SPECIES + str(i) + '.tif', categorical=True, all_touched=True)
        # https://pythonhosted.org/rasterstats/_modules/rasterstats/main.html#gen_zonal_stats
    # print(zs_species)
    heliGPS.loc[points_within.index.astype('uint64'), 'species'] = zs_species # set AoI identifier
    points_within.to_file(F_MPB_BUFFERED + str(i) + '.shp')

heliGPS.to_file(F_MPB_BUFFERED + '.shp')
# %%
# # visualize tree species data with heliGPS points and corresponding buffers
for i in [3]:
    fig, ax = plt.subplots(figsize=(10, 10))
    species_array = rxr.open_rasterio(F_SPECIES + str(i) + '.tif')
    species_array = species_array.rio.reproject(CRS)
    species_array.plot(cmap = 'cividis')
    boundsGdf.to_crs(CRS)
    boundsGdf.plot(ax=ax, alpha=0)
    # with rasterio.open(F_SPECIES + str(i) + '.tif') as species:
    #     # rasterio.plot.show(species.read(1), ax=ax, cmap = 'cividis')
    #     ax.imshow(species.read(1), cmap = 'cividis', interpolation ='nearest', extent=ax.get_window_extent)
    
    # heliGPS.plot(ax=ax, marker='o', markersize=2, color='red')
    points_within.to_crs(CRS)
    points_within.plot(ax=ax, color='red', alpha=0.5)
    ax.set_xlim([438141-500, 438141+500])
    ax.set_ylim([5814418-500, 5814418+500])
    # ax.set_axis_off()
    # plt.show()





# %%
# create pandas df just with species for selected points
df_species= pd.DataFrame(heliGPS[heliGPS['AoI']!=0]['species'].to_list())
# df_species['sum'] = df_species.sum(axis=1)
df_species.describe()
# %%
# visualize species data for sampled points
fig, ax = plt.subplots(1,3, figsize=(10,10))
for i in range(3):
    df_species = pd.DataFrame(heliGPS[heliGPS['AoI']==i+1]['species'].to_list())
    ax[i].bar(df_species.columns, df_species.count()/sum(df_species.count()))