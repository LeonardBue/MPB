from osgeo import gdal


# def clip_raster(shp_path, full_raster_path, clipped_raster_path):

# Clip raster using GDAL
# https://gis.stackexchange.com/questions/45053/gdalwarp-cutline-along-with-shapefile

rasPath = "D:\\TEMP\\_output\\_tif\\_Manning_Invasives\\Struc_Spec_1cm_near_fill.tif"

for i in range(1,4):
    poly =  f'E:\\Sync\\_Documents\\_Letter_invasives\\_Data\\_quadrats\\{i}.shp'
    ds = gdal.Warp(f"E:\\Sync\\_Documents\\_Letter_invasives\\_Data\\_quadrats\\clippedQd{i}.tif", rasPath, cropToCutline=True, cutlineDSName=poly, format="GTiff")
    ds = None # Close object