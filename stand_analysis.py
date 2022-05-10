# %% 
# Stand Analysis

import statistics
import fiona
import geopandas as gpd
import matplotlib.pyplot as plt
import math
import multiprocessing as mp
import numpy as np
import time
import sys

from functools import partial
from multiprocessing.pool import ThreadPool as Pool
from util.raster_analysis import sample_raster


F_AOI = './Data/Areas_of_Interest/AoI'
F_CMAP_SPECIES = './Data/Species_classification_2019/updated_species_list_alberta.csv'
F_GDB = './Data/MPB_AERIAL_SURVEY.gdb'
F_SPECIES = './Data/Species_classification_2019/Species_classification_2019_aoi'
F_LIDAR = './Data/Aria_Guo/summaries/'
F_ELEV_FIRST_RETURN = 'elevation_percentage_first_returns_above_1pnt37.tif'
F_ELEV_P99 = 'elevation_elev_p99.tif'
F_ELEV_CV = 'elevation_elev_cv.tif'
F_CLUSTER_CLEAN ='./Data/Aria_Guo/block_8_clusters_clean1.tif'


CRS = 'EPSG:3400'

n_cpus = int(mp.cpu_count()/2)
# %% 
# load heli-GPS data
yrs = [19] # range(11, 22) # [11, 19, 20, 21]
from util.read_gdb import read_heliGPS

MPB_locations = read_heliGPS(F_GDB, yrs, CRS)

# %%
# the quick and dirty (and maybe to memory hungry)
batchsize = 128

# MPB_locations['first_return'] = np.nan
buffer = False
bufferargs = {'distance': 60, 'resolution': 1, 'cap_style': 3}
feature_stats = 'mean'

test_df = MPB_locations.iloc[:1000, :].copy()

# tic = time.process_time()
# batches = np.array_split(range(len(test_df)), n_cpus) #math.ceil(len(test_df)/batchsize))
# for batch in batches:
#     print(f'Batch begins at index {batch[0]} and ends at {batch[-1]}')
    
#     test_df.loc[batch,'first_return'] = sample_raster(test_df.iloc[batch,:], 
#         raster_path=F_LIDAR+'z11/'+F_ELEV_FIRST_RETURN, value_name='first_return', 
#         buffer=buffer, bufferargs=bufferargs, stats=feature_stats)
#     test_df.loc[batch,'first_return'] = sample_raster(test_df.iloc[batch,:], 
#         raster_path=F_LIDAR+'z12/'+F_ELEV_FIRST_RETURN, value_name='first_return', 
#         buffer=buffer, bufferargs=bufferargs, stats=feature_stats)
# toc=time.process_time()
# print(f'Elapsed time with single core processing: {toc-tic} seconds.')

tic = time.process_time()
with Pool(processes=n_cpus) as pool:
    batches = np.array_split(range(len(test_df)), n_cpus) #math.ceil(len(test_df)/batchsize))
    input_data = []
    for batch in batches:
        input_data.append(test_df.iloc[batch,:])
        print(f'Batch begins at index {batch[0]} and ends at {batch[-1]}')
        
    elev_first_return = pool.map_async(partial(sample_raster, raster_path=F_LIDAR+'z11/'+F_ELEV_FIRST_RETURN,
                value_name='elev_first_return', buffer=buffer, bufferargs=bufferargs, stats=feature_stats), input_data).get()
    elev_first_return = pool.map_async(partial(sample_raster, raster_path=F_LIDAR+'z12/'+F_ELEV_FIRST_RETURN, 
                value_name='elev_first_return', buffer=buffer, bufferargs=bufferargs, stats=feature_stats), input_data).get()
    test_df.loc[:,'elev_first_return'] = np.hstack(np.asarray(elev_first_return, dtype=object))

    elev_p99 = pool.map_async(partial(sample_raster, raster_path=F_LIDAR+'z11/'+F_ELEV_P99, 
                value_name='elev_p99', buffer=buffer, bufferargs=bufferargs, stats=feature_stats), input_data).get()
    elev_p99 = pool.map_async(partial(sample_raster, raster_path=F_LIDAR+'z12/'+F_ELEV_P99, 
                value_name='elev_p99', buffer=buffer, bufferargs=bufferargs, stats=feature_stats), input_data).get()
    test_df.loc[:,'elev_p99'] = np.hstack(np.asarray(elev_p99, dtype=object))
    
    elev_cv = pool.map_async(partial(sample_raster, raster_path=F_LIDAR+'z11/'+F_ELEV_CV, 
                value_name='elev_cv', buffer=buffer, bufferargs=bufferargs, stats=feature_stats), input_data).get()
    elev_cv = pool.map_async(partial(sample_raster, raster_path=F_LIDAR+'z12/'+F_ELEV_CV, 
                value_name='elev_cv', buffer=buffer, bufferargs=bufferargs, stats=feature_stats), input_data).get()
    test_df.loc[:,'elev_cv'] = np.hstack(np.asarray(elev_cv, dtype=object))

    cluster_clean = pool.map_async(partial(sample_raster, raster_path=F_CLUSTER_CLEAN, 
                value_name='cluster_clean', buffer=buffer, bufferargs=bufferargs, stats='mean'), input_data).get()
    test_df.loc[:,'cluster_clean'] = np.hstack(np.asarray(cluster_clean, dtype=object))
    
pool.close()
pool.join()
toc=time.process_time()
print(f'Elapsed time with multi core processing on {n_cpus} cores: {toc-tic} seconds.')

# %%
