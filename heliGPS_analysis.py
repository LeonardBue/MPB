# %% 
## MPB heli-GPS Analysis

from cgi import test
from locale import normalize
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt

# import arcpy
# aprx = arcpy.mp.ArcGISProject(r"C:\Projects\YosemiteNP\Yosemite.aprx")
# insertLyr = arcpy.mp.LayerFile(r"C:\Projects\YosemiteNP\LayerFiles\Ranger Stations.lyrx")
# m = aprx.listMaps("Yosemite National Park")[0]
# refLyr = m.listLayers("Points of Interest")[0]
# m.insertLayer(refLyr, insertLyr, "BEFORE")
# aprx.saveACopy(r"C:\Projects\YosemiteNP\Yosemite_updated.aprx")


# geopandas included map, filtered to just this hemisphere
GDB_FILENAME = './Data/MPB_AERIAL_SURVEY.gdb'
CRS = 4269

print(GDB_FILENAME)

# geopandas included map, filtered to just this hemisphere
# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
# worlsd = world.to_crs(CRS)
# westhem = world[(world['continent'] == 'North America') | 
#                 (world['continent'] == 'South America')]
# base = westhem.plot(color='white', edgecolor='black', figsize=(11, 11))

layers = fiona.listlayers(GDB_FILENAME)
heliGPS_layers = sorted([lyr for lyr in layers if lyr.endswith('x')])
polygons = sorted([lyr for lyr in layers if lyr.endswith('p')])
heliGPS_last_3_years = heliGPS_layers[-3:]

heliGPS = None
for l in heliGPS_last_3_years:
    if heliGPS is None:
        heliGPS = gpd.read_file(GDB_FILENAME, layer=l)
    else:
        survey = gpd.read_file(GDB_FILENAME, layer=l)
        heliGPS = pd.concat([heliGPS, survey])
        
heliGPS = heliGPS.to_crs(CRS)
heliGPS = heliGPS.reset_index()
heliGPS

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

