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
from osgeo import gdal
import pandas as pd
import rasterstats as rs
import rasterio
from rasterio.plot import plotting_extent
import rioxarray as rxr


F_GDB = './Data/MPB_AERIAL_SURVEY.gdb'
F_SPECIES = './Data/Species_classification_2019/Species_classification_2019_aoi'
F_AOI = './Data/Areas_of_Interest/AoI'
F_GRAVEL_ROADS = './Data/Road_Network/_gravelRoads.shp'
F_PAVED_ROADS = './Data/Road_Network/_pavedRoads.shp'
F_BUFFER = './Data/MPB_buffer/heliGPS_buffer.shp'
F_MPB_BUFFERED = './Data/MPB_buffer/heliGPS_species_gdf'
F_CMAP_SPECIES = './Data/Species_classification_2019/updated_species_list_alberta.csv'

CRS = 'EPSG:4269'

POINTBUFFER = 2

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
# add a column for the area of interest; defaults to 0 --> not in AoI
heliGPS['AoI'] = np.zeros(max(heliGPS.count()))

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
# resolution+1 corresponds to a square
# https://shapely.readthedocs.io/en/latest/manual.html#object.buffer
heliGPS_poly["geometry"] = heliGPS.geometry.buffer(1, resolution=1, cap_style = 3) 
heliGPS_poly.head()
heliGPS_poly.to_file(F_BUFFER)

# %% 
# Read tree species data
species = []
heliGPS_species_gdf = None

for i in [1]:
    species = rxr.open_rasterio(F_SPECIES + str(i) + '.tif', masked=True).squeeze()
    species = species.rio.reproject(CRS)
#     species = gr.from_file(F_SPECIES + str(i) + '.tif')

    # cleanup: set values of nan, to 0 corresponding to no trees
    species_clean = species.fillna(0)
    # ax = species.value.plot.hist(bins=max(species.value.unique()), figsize=(8,8))
    print(f'Loaded species data for AoI{i}')

    # extract pixel values for each polygon of the bufferes heli-GPS data
#     affine=species.rio.transform()

#     heliGPS_species = rs.zonal_stats(F_BUFFER,``
#                                         species_clean.values, # gdf_species_no_zeros.values,
#                                         nodata=0,
#                                         # affine=rasterio.transform.from_origin(north, west, xsize, ysize),
#                                         affine=affine,
#                                         geojson_out=True,
#                                         copy_properties=True,
#                                         categorical=True, 
# )                               

    # get points within AoI and update AoI column value in heliGPS data
    aoi = gpd.read_file(F_AOI + str(i) + '.shp')
    points_within = gpd.sjoin(heliGPS, aoi, op='within')
    heliGPS.loc[points_within.index, 'AoI'] = i
    
    gdal_species = gdal.Open(F_SPECIES + str(i) + '.tif')
    gdal.Warp(gdal_species, gdal_species,dstSRS = CRS)
    gdal_band = gdal_species.GetRasterBand(1)
    no_data_val = gdal_band.GetNoDataValue() # seems to be 255
    data_array = gdal_species.ReadAsArray().astype(np.float)
    data_array
    # replacing missing values by nan
    if np.any(data_array == no_data_val):
        data_array[data_array == no_data_val] = np.nan

    plt.figure(figsize = (8, 8))
    img = plt.imshow(data_array, cmap = 'cividis') 

    prj = gdal_species.GetProjection()

    coords = points_within[['latitude', 'longitude']].to_numpy()

    row, col = pixel_from_coords(gdal_species, data_array, pos = coords[1]) 



#     if heliGPS_species_gdf is None:
#         heliGPS_species_gdf = gpd.GeoDataFrame.from_features(heliGPS_species)
#     else:
#         heliGPS_species_gdf = pd.concat([heliGPS_species_gdf, gpd.GeoDataFrame.from_features(heliGPS_species)])
    
#     print(f'Calculated correlated buffers and species data in AoI{i}') 
        
# if heliGPS_species_gdf.crs is None:
#     heliGPS_species_gdf.set_crs(CRS, inplace=True) 
# else:
#     heliGPS_species_gdf.to_crs(CRS, inplace=True)

# cmap_species = pd.read_csv(F_CMAP_SPECIES)

# for col in heliGPS_species_gdf.columns:
#     if col in [float(val) for val in cmap_species['Value Code']]:
#         heliGPS_species_gdf.rename(columns = {col : cmap_species[cmap_species['Value Code'] == col]['NFI Code'].item()}, inplace = True)
        
# heliGPS_species_gdf.head()

# print('Writing files.')
# save file to GeoJSON and shapefile
# heliGPS_species_gdf.to_file(F_MPB_BUFFERED + '.geojson')
# heliGPS_species_gdf.to_file(F_MPB_BUFFERED + '.shp')

# %%
# visualize tree species data with heliGPS points and corresponding buffers
fig, ax = plt.subplots(figsize=(10, 10))
species.plot(ax=ax, cmap = 'cividis')
heliGPS.plot(ax=ax, marker='o', markersize=2, color='red')
heliGPS_poly.plot(ax=ax, color='red', alpha=0.5)
# ax.set_xlim([254000, 256000])
# ax.set_ylim([6.016e6, 6.018e6])
# ax.set_axis_off()
plt.show()

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
