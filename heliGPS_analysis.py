# %% 
# MPB heli-GPS Analysis

import fiona
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import rioxarray as rxr
from osgeo import gdal
from rasterstats import zonal_stats
from rasterio.plot import show
from scipy import stats
from shapely.geometry import box, Point
from sympy import rotations


F_GDB = './Data/MPB_AERIAL_SURVEY.gdb'
F_SPECIES = './Data/Species_classification_2019/Species_classification_2019_aoi'
F_AOI = './Data/Areas_of_Interest/AoI'
F_GRAVEL_ROADS = './Data/Road_Network/_gravelRoads.shp'
F_PAVED_ROADS = './Data/Road_Network/_pavedRoads.shp'
F_ROADS = './Data/Road_Network/bufferedRoads.shp'
F_BUFFER = './Data/MPB_buffer/heliGPS_buffer.shp'
F_MPB_BUFFERED = './Data/MPB_buffer/heliGPS_species'
F_CMAP_SPECIES = './Data/Species_classification_2019/updated_species_list_alberta.csv'


CRS = 'EPSG:3400'

# %% 
# load heli-GPS data

layers = fiona.listlayers(F_GDB)
heliGPS_layers = sorted([lyr for lyr in layers if lyr.endswith('x')])
polygons = sorted([lyr for lyr in layers if lyr.endswith('p')])
relevant_years = [19]
relevant_layers = [lyr for lyr in heliGPS_layers if any([str(yr) in lyr for yr in relevant_years])]

heliGPS = None
for l in relevant_layers:
    if heliGPS is None:
        heliGPS = gpd.read_file(F_GDB, layer=l)
    else:
        survey = gpd.read_file(F_GDB, layer=l)
        heliGPS = pd.concat([heliGPS, survey])
        
heliGPS = heliGPS.to_crs(CRS)
heliGPS.drop(heliGPS[(heliGPS[('att_stage').casefold()].isin(['', ' ']))].index, inplace=True)
heliGPS = heliGPS.reset_index(drop=False)
heliGPS.columns = heliGPS.columns.str.lower()
heliGPS # imported data should be of geometry type point. check with: heliGPS.geom_type.head()

# These data points do not have an attack stage (att_stage) assigned to them, thus they will be excluded from further anlysis.
# heliGPS.loc[heliGPS['att_stage'] == ''].explore()

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
# print(heliGPS.zone.value_counts(ascending=True))
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
heliGPS['aoi'] = np.zeros(max(heliGPS.count()), dtype=int)
heliGPS['species'] = np.zeros(max(heliGPS.count()))

# read tree species names
df_names = pd.read_csv(F_CMAP_SPECIES)
cmap_species = dict(zip(df_names['Value Code'], df_names['NFI Code']))

# iterate through AoI number
for i in [1]:#, 2, 3]:
    with rasterio.open(F_SPECIES + str(i) + '.tif') as species:
        bounds = species.bounds # get bounding box raster
        boundsGdf = gpd.GeoDataFrame({"geometry":[box(*bounds).buffer(-60)]}, crs=species.crs) # to GeoDataFrame

    points_within = heliGPS.to_crs(species.crs)
    points_within = heliGPS.overlay(boundsGdf, how='intersection') # Keep points within raster bounds -  https://geopandas.org/en/stable/docs/user_guide/set_operations.html
    points_within.set_index('index', inplace=True)

    heliGPS.loc[points_within.index.astype('uint64'), 'aoi'] = i # set AoI identifier
    # points_within = points_within.to_crs(CRS)
    points_within['geometry'] = points_within.geometry.buffer(distance=60, resolution=1, cap_style = 3) # Do some buffering
    try:
        zs_species = zonal_stats(points_within, F_SPECIES + str(i) + '.tif', categorical=True, 
                                    all_touched=True, category_map = cmap_species)
            # https://pythonhosted.org/rasterstats/_modules/rasterstats/main.html#gen_zonal_stats
        # print(zs_species)
        heliGPS.loc[points_within.index.astype('uint64'), 'species'] = zs_species # set AoI identifier
        points_within.to_file(F_MPB_BUFFERED + str(i) + '.shp')
    except Exception as e: 
        print(e)

# heliGPS.to_file(F_MPB_BUFFERED + '.shp')
# %%
# # visualize tree species data with heliGPS points and corresponding buffers
df_names = pd.read_csv(F_CMAP_SPECIES)
cmap_labels = dict(zip(df_names['Value Code'], df_names['Common Species Name']))

for i in [1]:
    fig, ax = plt.subplots(figsize=(10, 10))
    species_array = rxr.open_rasterio(F_SPECIES + str(i) + '.tif', masked=True).squeeze()
    species_array = species_array.rio.reproject('EPSG:4269')
    species_array.plot(cmap = 'cividis', alpha=1, 
                        cbar_kwargs={'ticks': list(cmap_labels.keys()), 'spacing': 'proportional', 
                        'label': 'Species Classification 2019 in Alberta', 'shrink': 0.6})
    # boundsGdf.to_crs('EPSG:4269', inplace+True)
    # boundsGdf.plot(ax=ax, alpha=0)
    # with rasterio.open(F_SPECIES + str(i) + '.tif') as species:
    #     # rasterio.plot.show(species.read(1), ax=ax, cmap = 'cividis')
    #     ax.imshow(species.read(1), cmap = 'cividis', interpolation ='nearest', extent=ax.get_window_extent)
    
    heliGPS.to_crs('EPSG:4269', inplace=True)
    # heliGPS.plot(ax=ax, marker='o', markersize=2, color='k')
    points_within.to_crs('EPSG:4269', inplace=True)
    points_within.plot(ax=ax, color='red', alpha=0.5)
    ax.set_xlim([-118.66, -118.6])
    ax.set_ylim([54.25, 54.28])
    ax.set_title('Buffered GPS Locations and Species Raster Data', fontsize=14)

    # ax.set_axis_off()
    # fig.suptitle('Buffered GPS Locations and Species Raster Data', fontsize=14)
    plt.tight_layout()    
    graphics_file = './graphics/bufferedPointsExample'
    # plt.savefig(graphics_file + '.pdf')
    # plt.savefig(graphics_file + '.png')


# %%
# create pandas df just with species for selected points
df_species= pd.DataFrame(heliGPS[heliGPS['aoi']!=0]['species'].to_list())
df_species = df_species.reindex(sorted(df_species.columns), axis=1)
# df_species['sum'] = df_species.sum(axis=1)
df_species.describe()


# %%
# visualize species data for sampled points
color_map = plt.get_cmap("cividis")


# data preparation
df_names = pd.read_csv(F_CMAP_SPECIES)
cmap_labels = dict(zip(df_names['NFI Code'], df_names['Common Species Name']))
labels = []
for l in df_species.columns.to_list():
    labels.append(cmap_labels[l])


rel_tree_species = {}
species_aoi = {}
for i in range(3):
    species_aoi[i] = pd.DataFrame(heliGPS[heliGPS['aoi']==i+1]['species'].to_list())
    species_aoi[i] = species_aoi[i].reindex(sorted(species_aoi[i].columns), axis=1)
    rel_tree_species[i] = species_aoi[i]/(25*len(species_aoi[i]))   
    # ax[i].bar(species_aoi[i].columns, species_aoi[i].count())#/sum(df_species.count()))    
# rel_tree_species = np.sum(rel_tree_species, axis=0)


# plotting
# ax[-1].bar(labels,
# rel_tree_species)
bar_width = 0.2
bar_position = {}
fig, axes = plt.subplots(2,1, figsize=(10,10))
for i in range(len(species_aoi.keys())):
    bar_position[i] = np.arange(len(species_aoi[i].columns)) + i *bar_width
    axes[0].bar(bar_position[i], species_aoi[i].count(), width=bar_width, 
                color=color_map.colors[int(i*255/(len(species_aoi.keys())-1))], 
                label=f'Area of Interest {i+1}')
    axes[0].set_ylabel('Absolute number of pixels')
    axes[1].bar(bar_position[i], rel_tree_species[i].sum(), 
                width=bar_width, color=color_map.colors[int(i*255/(len(species_aoi.keys())-1))], label=f'Area of Interest {i+1}')

# style
for ax in axes.flat:
    ax.set_xticks(bar_position[1])
    ax.set_xticklabels(labels)
    ax.set_xlabel('Most dominant species')
    ax.legend()
    ax.label_outer()
    ax.grid(axis='y')
axes[1].set_ylabel('Relative number of pixels')
# axes[1].set_ylim([0,1])
plt.xticks()
fig.suptitle('Most Dominant Tree Species in areas of interest in 2011', fontsize=14)
plt.tight_layout()

graphics_file = './graphics/freq_species_11'
# plt.savefig(graphics_file + '.pdf')
# plt.savefig(graphics_file + '.png')


#TODO titles, and add second line relative to beetle infestation, separate for cumulative
#TODO repeat for 2011, check if issues for earlier years


# %%
# first order statistics of species
bands={}
for i in [1, 2, 3]:
    with rasterio.open(F_SPECIES + str(i) + '.tif') as sp:
        bands['aoi'+str(i)] = sp.read(1)[0:1665, 0:1665].reshape((1, -1)).squeeze()
        print(bands['aoi'+str(i)].shape)
df = pd.DataFrame(bands)

for i in range(len(bands.keys())):
    print(df['aoi'+str(i+1)].value_counts()/len(bands['aoi1']))
    print('std', np.std(df['aoi'+str(i+1)]))


# %%
# load road network 

paved_rds = gpd.read_file(F_PAVED_ROADS)
paved_rds['type'] = 'paved'
paved_rds.to_crs(CRS, inplace = True)
gravel_rds = gpd.read_file(F_GRAVEL_ROADS)
gravel_rds['type'] = 'gravel'
gravel_rds.to_crs(CRS, inplace = True)
# %%
# buffer roads to flyable range
try:
    roads = pd.concat([paved_rds, gravel_rds])
except Exception as e: 
    print(e)
    roads = paved_rds
roads['aoi'] = np.zeros(max(roads.count()), dtype=int)
roads = roads.reset_index(drop=False)

for i in [1, 2, 3]:
    with rasterio.open(F_SPECIES + str(i) + '.tif') as species:
        bounds = species.bounds # get bounding box raster
        boundsGdf = gpd.GeoDataFrame({"geometry":[box(*bounds)]}, crs=species.crs) # to GeoDataFrame

    rds_within = roads.to_crs(species.crs)
    rds_within = roads.overlay(boundsGdf, how='intersection') # Keep roads within raster bounds
    rds_within.set_index('index', inplace=True)

    roads.loc[rds_within.index.astype('uint64'), 'aoi'] = i # set AoI identifier
    # rds_within = rds_within.to_crs(CRS)

roads.drop(roads[(roads['aoi'].isin([0]))].index, inplace=True)
roads['geometry'] = roads.geometry.buffer(distance=250, cap_style = 1) # Do some buffering
roads.to_file(F_ROADS + '.shp')
# %%
# calculate LPP percentage for every gps point
# get all ponits with percentage above threshold and within buffered roads
# for every area rank the resulting points by number of trees
# group by year 
# define nuymber of flight areas and flight ares themselves