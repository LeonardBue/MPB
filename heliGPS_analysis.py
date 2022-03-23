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
from shapely.geometry import box


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
heliGPS = heliGPS.reset_index()
heliGPS # imported data should be of geometry type point. check with: heliGPS.geom_type.head()

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
# load data and extract species of trees for points in area of interest
# add a column for the area of interest; defaults to 0 --> not in AoI
heliGPS['AoI'] = np.zeros(max(heliGPS.count()), dtype=int)
heliGPS['species'] = np.zeros(max(heliGPS.count()))

# iterate through AoI number
for i in [1, 2, 3]:
    with rasterio.open(F_SPECIES + str(i) + '.tif') as species:
        bounds = species.bounds # get bounding box raster
        boundsGdf = gpd.GeoDataFrame({"id":1,"geometry":[box(*bounds)]}, crs=species.crs) # to GeoDataFrame

    points_within = heliGPS.to_crs(species.crs)
    points_within = heliGPS.overlay(boundsGdf, how='intersection') # Keep points within raster bounds -  https://geopandas.org/en/stable/docs/user_guide/set_operations.html
    heliGPS.loc[points_within.index, 'AoI'] = i # set AoI identifier

    points_within = points_within.to_crs(CRS)
    points_within['geometry'] = points_within.geometry.buffer(distance=60, resolution=1, cap_style = 3) # Do some buffering
    zs_species = zonal_stats(points_within, F_SPECIES + str(i) + '.tif', categorical=True)
    # print(zs_species)
    # add species to heli-GPS dataset
    heliGPS.loc[points_within.index, 'species'] = zs_species # set AoI identifier

# heliGPS.to_file(F_MPB_BUFFERED + '.geojson')
# heliGPS.to_file(F_MPB_BUFFERED + '.shp')

# # visualize tree species data with heliGPS points and corresponding buffers
# for i in [1, 2, 3]:
#     fig, ax = plt.subplots(figsize=(10, 10))
#     with rasterio.open(F_SPECIES + str(i) + '.tif') as species:
#         species_data = species.read().reproject(species, CRS)
#         rasterio.plot.show(species.read(), ax=ax, cmap = 'cividis')
#         # heliGPS.plot(ax=ax, marker='o', markersize=2, color='red')
#         # points_within.plot(ax=ax, color='red', alpha=0.5)
#         # ax.set_xlim([254000, 256000])
#         # ax.set_ylim([6.016e6, 6.018e6])
#         # ax.set_axis_off()
#         plt.show()


# %%
def test_reproject_epsg():
    with rasterio.Env():
        with rasterio.open('tests/data/RGB.byte.tif') as src:
            source = src.read(1)

        dst_crs = {'init': 'EPSG:3857'}
        out = np.empty(src.shape, dtype=np.uint8)
        rasterio.warp.reproject(
            source,
            out,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=DST_TRANSFORM,
            dst_crs=dst_crs,
            resampling=Resampling.nearest)

# %%
def pixel_from_coords(gdal_data, data_array, pos):
    """
    maps coordinates from georeferenced space to pixel column row values
    GDALDataset.GetGeoTransform() returns affine transformation. more information is available here:
    https://gdal.org/tutorials/geotransforms_tut.html
    GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
    GT(1) w-e pixel resolution / pixel width.
    GT(2) row rotation (typically zero).
    GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
    GT(4) column rotation (typically zero).
    GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).
    """

    gt = gdal_data.GetGeoTransform()
    col = int((pos[0] - gt[0]) / gt[1])
    row = int((pos[1] - gt[4]) / -gt[5])
    data_array
    return col, row, data_array[row][col]

def retrieve_pixel_value(pos, gdal_data):
    """Return floating-point value that corresponds to given point."""
    x, y = pos[0], pos[1]
    forward_transform = affine.Affine.from_gdal(*gdal_data.GetGeoTransform())
    reverse_transform = ~forward_transform
    col, row = reverse_transform * (x, y)
    col, row = int(col + 0.5), int(row + 0.5)
    pixel_coord = col, row

    data_array = np.array(gdal_data.GetRasterBand(1).ReadAsArray())
    return data_array[pixel_coord[0]][pixel_coord[1]]


# %%
