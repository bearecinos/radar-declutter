import arcpy
import numpy as np

_SetupRun = False
_CellSize = 30.0

def Setup():
    """Performs most of the initialisation needed by arcpy."""
    global _raster, _aspectRaster, _slopeRaster, _SetupRun, _CellSize
    try:
        _raster = arcpy.Raster("default.gdb/srtm30_dem_utm45n_crop")
        _slopeRaster = arcpy.Raster("default.gdb/Slope_srtm")
        _aspectRaster = arcpy.Raster("default.gdb/Aspect_srtm")
    except RuntimeError as e:
        if "ERROR 000732" in e.message:
            print "Error in Generate.py, could not load one of the input rasters."
            print "Please check folder default.gdb exists."
            print "Must contain srtm30_dem_utm45n_crop, Slope_srtm, Aspect_srtm."
            print "This can be checked by looking at folder in arcMap."
            return -1
    _CellSize = = float(arcpy.GetRasterProperties_management(_raster,"CELLSIZEX").getOutput(0))
    _SetupRun = True
    return 0

def storeRasters(): 
    """Save full rasters used by arcpy as numpy arrays to avoid arcpy use in future"""
    if not _SetupRun:
        if Setup():
            return -1
    height = arcpy.RasterToNumPyArray(_raster)
    aspect = arcpy.RasterToNumPyArray(_aspectRaster)
    slope = arcpy.RasterToNumPyArray(_slopeRaster)
    ex = _raster.extent
    left = ex.XMin
    low = ex.YMin
    line = str(left)+","+str(low)+","+str(_CellSize)
    np.savez_compressed("maps.npz",heightmap=height,aspect=aspect,slope=slope)
    with open("info.txt","w") as f:
        f.write(line)
    return 0
    

