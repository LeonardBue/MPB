import fiona
import geopandas as gpd
import pandas as pd

def read_heliGPS(filepath, years, crs='EPSG:3400'):
    layers = fiona.listlayers(filepath)

    heliGPS_layers = sorted([lyr for lyr in layers if lyr.endswith('x')])
    # polygons = sorted([lyr for lyr in layers if lyr.endswith('p')])
    relevant_layers = [lyr for lyr in heliGPS_layers if any([str(yr) in lyr for yr in years])]

    # load relevant layers as specified and concatenate them to one gdb
    heliGPS = None
    for l in relevant_layers:
        if heliGPS is None:
            heliGPS = gpd.read_file(filepath, layer=l)
        else:
            survey = gpd.read_file(filepath, layer=l)
            heliGPS = pd.concat([heliGPS, survey])
            
    # project to desired crs and clean up gdb from data without attack stages and non-MPB damages
    # All remaining points are associated either to attack stage "Green" or "Red".
    heliGPS = heliGPS.to_crs(crs)
    heliGPS.drop(heliGPS[(heliGPS[('att_stage').casefold()].isin(['', ' ']))].index, inplace=True)
    heliGPS.drop(heliGPS[(heliGPS[('dmg_desc').casefold()] != 'Mountain pine beetle')].index, inplace=True)
    # reset index to have no duplicates in index
    heliGPS = heliGPS.reset_index(drop=True)
    heliGPS['id'] = heliGPS.index
    heliGPS.columns = heliGPS.columns.str.lower()
    return heliGPS