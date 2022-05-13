import rasterio as rio
from rasterstats import zonal_stats
import pandas as pd
import sys

def buffer_points(gdf, distance, **kwargs):
    """
    Applies buffer of distance to points or polygons and returns the geometry.
    See http://shapely.readthedocs.io/en/latest/manual.html#object.buffer for details about buffering
    """
    bp_geometry = gdf.geometry.buffer(distance, **kwargs)
    return bp_geometry

def sample_raster(gdf, raster_path, value_name ='value_name', buffer = False, bufferargs = {}, **kwargs):
    """
    Returns raster data given for point coordinates.
    Instead of points, polygons can be used as well.
    For details refer to:
    https://pythonhosted.org/rasterstats/rasterstats.html?highlight=zonal_stats#rasterstats.gen_zonal_stats
    """
    
    if any(gdf.geometry.type != 'Point'):
        print('Some geometries were no Points, and will be dropped')
        gdf.drop(gdf[gdf.geometry.type != 'Point'].index, inplace = True)
    with rio.open(raster_path) as r:
        points = gdf.to_crs(r.crs)

    if buffer:
        points.loc[:, 'geometry'] = buffer_points(points, **bufferargs)
    # else:
        # kwargs['stats'] = 'mean' # because returning all other values for one pixel doesn't make sense.
    try:
        values = zonal_stats(points, raster_path, **kwargs) # categorical=False, all_touched=True)
        df_values = pd.DataFrame(values)
        df_values.set_index(points.index, inplace = True)
        # if kwargs.get('categorical') == True:
        #     if len(values[0].keys()) == 1:
        #         gdf.loc[points[df_values.notnull().any(axis=1)].index.astype('uint64'), value_name] = df_values[df_values.notnull().all(axis=1)].keys
        #     else:
        #         gdf.loc[points[df_values.notnull().any(axis=1)].index.astype('uint64'), value_name] = df_values[df_values.notnull().any(axis=1)].to_dict(orient='records')
                
        # else:
        if len(values[0].keys()) == 1:
            gdf.loc[points[df_values.notnull().any(axis=1)].index.astype('uint64'), value_name] = df_values[df_values.notnull().any(axis=1)].values
            
        else:
            gdf.loc[points[df_values.notnull().any(axis=1)].index.astype('uint64'), value_name] = df_values[df_values.notnull().any(axis=1)].to_dict(orient='records') # set specified value
            
    except Exception as e: 
        print(e)
    
    return gdf[value_name].values