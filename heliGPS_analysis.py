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
from sympy import are_similar, rotations


F_GDB = './Data/MPB_AERIAL_SURVEY.gdb'
F_SPECIES = './Data/Species_classification_2019/Species_classification_2019_aoi'
F_AOI = './Data/Areas_of_Interest/AoI'
F_FLIGHT_AREA = './Data/Areas_of_Interest/FA_250_reltop50.shp'
F_GRAVEL_ROADS = './Data/Road_Network/_gravelRoads.shp'
F_PAVED_ROADS = './Data/Road_Network/_pavedRoads.shp'
F_ROADS = './Data/Road_Network/bufferedRoads_250.shp'
F_BUFFER = './Data/MPB_buffer/heliGPS_buffer.shp'
F_MPB_BUFFERED = './Data/MPB_buffer/heliGPS_buffered_'
F_CMAP_SPECIES = './Data/Species_classification_2019/updated_species_list_alberta.csv'


CRS = 'EPSG:3400'

# %% 
# load heli-GPS data

layers = fiona.listlayers(F_GDB)
heliGPS_layers = sorted([lyr for lyr in layers if lyr.endswith('x')])
polygons = sorted([lyr for lyr in layers if lyr.endswith('p')])
relevant_years = [11, 19, 20, 21]
relevant_layers = [lyr for lyr in heliGPS_layers if any([str(yr) in lyr for yr in relevant_years])]

# load relevant layers as specified and concatenate them to one gdb
heliGPS = None
for l in relevant_layers:
    if heliGPS is None:
        heliGPS = gpd.read_file(F_GDB, layer=l)
    else:
        survey = gpd.read_file(F_GDB, layer=l)
        heliGPS = pd.concat([heliGPS, survey])
        
# project to desired crs and clean up gdb from data without attack stages and non-MPB damages
# All remaining points are associated either to attack stage "Green" or "Red".
heliGPS = heliGPS.to_crs(CRS)
heliGPS.drop(heliGPS[(heliGPS[('att_stage').casefold()].isin(['', ' ']))].index, inplace=True)
heliGPS.drop(heliGPS[(heliGPS[('dmg_desc').casefold()] != 'Mountain pine beetle')].index, inplace=True)
# reset index to have no duplicates in index
heliGPS = heliGPS.reset_index(drop=True)
heliGPS['id'] = heliGPS.index
heliGPS.columns = heliGPS.columns.str.lower()
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
# print(heliGPS.zone.value_counts(ascending=True))
print(heliGPS.loc[heliGPS['att_stage'] == 'Green'].survyear.value_counts(ascending=True))


# %% 
# Visualize some statistics of the data

# Distribution of affected trees, grouped by att_stage
ax1 = heliGPS.plot.hist(column=['num_trees'], by='att_stage', figsize=(10, 8))

# boxplot of the spatial distribution for red and green attack trees
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
for i in [1, 2, 3]:
    _raster_file = F_SPECIES + str(i) + '.tif'
    with rasterio.open(F_SPECIES + str(i) + '.tif') as species:
        bounds = species.bounds # get bounding box raster
        boundsGdf = gpd.GeoDataFrame({"geometry":[box(*bounds).buffer(-60)]}, crs=species.crs) # to GeoDataFrame

    points_within = heliGPS.to_crs(species.crs)
    points_within = heliGPS.overlay(boundsGdf, how='intersection') # Keep points within raster bounds -  https://geopandas.org/en/stable/docs/user_guide/set_operations.html
    points_within.set_index('id', inplace=True)

    heliGPS.loc[points_within.index.astype('uint64'), 'aoi'] = i # set AoI identifier
    # points_within = points_within.to_crs(CRS)
    points_within['geometry'] = points_within.geometry.buffer(distance=60, resolution=1, cap_style = 3) # Do some buffering
    try:
        zs_species = zonal_stats(points_within, F_SPECIES + str(i) + '.tif', categorical=True, 
                                    all_touched=True, category_map = cmap_species)
            # https://pythonhosted.org/rasterstats/_modules/rasterstats/main.html#gen_zonal_stats
        # print(zs_species)
        heliGPS.loc[points_within.index.astype('uint64'), 'species'] = zs_species # set AoI identifier
        # points_within.to_file(F_MPB_BUFFERED + str(i) + '.shp')
    except Exception as e: 
        print(e)

# for yr in relevant_years:
#     heliGPS.loc[heliGPS.survyear == 2000+yr].to_file(f'{F_MPB_BUFFERED}{yr}.shp', )

# %%
# visualize tree species data with heliGPS points and corresponding buffers
df_names = pd.read_csv(F_CMAP_SPECIES)
cmap_labels = dict(zip(df_names['Value Code'], df_names['Common Species Name']))

for i in [1,2,3]:
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
    # leg = ax.get_legend()
    # leg.legendHandles[2].set_color(color_map.colors[int(i*255/(len(species_aoi.keys())-1))])
axes[1].set_ylabel('Relative number of pixels')



# axes[1].set_ylim([0,1])
plt.xticks()
fig.suptitle('Most Dominant Tree Species in areas of interest in 2011', fontsize=14)
plt.tight_layout()

graphics_file = './graphics/freq_species_11'
# plt.savefig(graphics_file + '.pdf')
# plt.savefig(graphics_file + '.png')


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
# get roads within area of interest
try:
    roads = pd.concat([paved_rds, gravel_rds])
except Exception as e: 
    print(e)
    roads = paved_rds
roads['aoi'] = np.zeros(max(roads.count()), dtype=int)
roads = roads.reset_index(drop=True)
roads['id'] = roads.index

for i in [1, 2, 3]:
    aoi = gpd.read_file(F_AOI + str(i) + '.shp')
    aoi.to_crs(CRS, inplace = True)
    rds_within = roads.clip(aoi)

    roads.loc[rds_within['id'].astype('uint64'), 'aoi'] = i # set AoI identifier

roads.drop(roads[(roads['aoi'].isin([0]))].index, inplace=True)
# roads['geometry'] = roads.geometry.buffer(distance=250, cap_style = 1) # Do some buffering
roads.plot()

# %%
# roads that contain heliGPS points, buffered to range of flight
roads.to_crs(CRS, inplace = True)
buffered_points = heliGPS.to_crs(CRS, inplace = True)
buffered_points = heliGPS.geometry.buffer(distance=250, cap_style = 1)

# roads.loc[rds_within['id'].astype('uint64'), 'aoi'] = i # set AoI identifier
flight_area = roads.clip(buffered_points, True)
flight_area.drop(flight_area[flight_area.geom_type == 'MultiLineString'].index, inplace=True)
flight_area['geometry'] = flight_area.geometry.buffer(distance=250, cap_style = 1)

_combined_area = flight_area.geometry.unary_union
flight_area = gpd.GeoDataFrame([polygon for polygon in _combined_area.geoms])
flight_area.rename(columns={0: 'geometry'}, inplace = True)
flight_area.set_geometry('geometry', inplace=True)
flight_area.set_crs(CRS, inplace = True)
# flight_area.to_file(F_ROADS)
# %%
# calculate LPP percentage for every gps point from species
heliGPS['lpp'] = 0.
heliGPS.loc[heliGPS.species != 0.0, 'lpp'] = heliGPS[heliGPS.species != 0.0].species.apply(lambda x: x.get('PINU.CON')).div(25)
heliGPS['lpp'] = heliGPS['lpp'].fillna(0.)

# associate points with flight areas
flight_area['aoi'] = 0
flight_area['lpp'] = 0
flight_area['red']= 0
flight_area['grey']= 0
flight_area['fallen']= 0
flight_area['num_trees'] = 0.
flight_area['survyear'] = 'None'
flight_area['id'] = flight_area.index

for i in [1, 2, 3]:    
    for p in flight_area.index:
        red, grey, fallen = 0., 0., 0.
        yrs = {2011: 0., 2019: 0., 2020: 0., 2021: 0.}
        _polygon = gpd.GeoDataFrame(index=[p], geometry=[flight_area.geometry[p]], crs=CRS)
        _pts = gpd.overlay(_polygon, heliGPS[heliGPS.aoi == i], how='intersection', keep_geom_type=False)
        if len(_pts) > 0:
            for idx in _pts.index:
                _pt = _pts.iloc[[idx]]
                if _pt.survyear.item() == 2011:
                    fallen += _pt.num_trees.item()
                elif _pt.survyear.item() == 2019:
                    grey += _pt.num_trees.item()
                elif (_pt.survyear.item() == 2020 or _pt.survyear.item() == 2021):
                    red += _pt.num_trees.item()
                yrs[_pt.survyear.item()] += 1
            flight_area.at[p, 'red'] = red
            flight_area.at[p, 'grey'] = grey
            flight_area.at[p, 'fallen'] = fallen
            flight_area.at[p, 'lpp'] = _pts.lpp.mean()
            flight_area.at[p, 'num_trees'] = _pts.num_trees.sum()
            flight_area.at[p, 'aoi'] = _pts.aoi.mode()
            flight_area.at[p, 'survyear'] = yrs
flight_area.num_trees = flight_area[['red', 'grey', 'fallen']].sum(axis=1)
flight_area.drop(flight_area[flight_area.aoi == 0].index, inplace=True) 
# flight_area.drop(flight_area[flight_area[['red', 'grey', 'fallen']].sum(axis=1)<=3].index, inplace=True)
# flight_area.drop(flight_area[flight_area.num_trees < 9].index, inplace=True) # only top 50% (9) or top 25% (18)
# flight_area.drop(flight_area[flight_area.lpp != 0.].index, inplace=True)
# flight_area.to_file(F_FLIGHT_AREA)

# %%
# Comparison of the different flight areas
_fa = pd.DataFrame(flight_area.survyear.to_list())
flight_area[_fa.columns] = _fa.values
# find threshold for upper 50 % and 25 %
_cols = [2011, 2019, 2020, 2021]
_cols_n = ['2011n', '2019n', '2020n', '2021n']
_flight_area = (flight_area[_cols] - flight_area[_cols].mean())/flight_area[_cols].std()
_flight_area['combined']=_flight_area.sum(axis=1)
_flight_area.describe() # -6.772848e-01 corresponds to top 50 % and 3.136384e-01 to top 25 %

flight_area[_cols_n] = _flight_area.drop(columns='combined') # add columns with normalized tree counts
flight_area.rename(columns={2011: '2011', 2019: '2019', 2020: '2020', 2021: '2021'}, inplace=True) 
# (flight_area.loc[(flight_area[_cols_n].sum(axis=1) > -6.772848e-01)]).to_file(F_FLIGHT_AREA)

# %%
