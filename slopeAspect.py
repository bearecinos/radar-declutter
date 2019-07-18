import arcpy
import numpy as np
from projection import *
import os
_NODATA = np.nan

# Warning: previous versions were split into several parts due to getting memory error
# For the same reason, resampling is done by arcpy
def makeRasters(source, sampleLat, sampleLon, cellSize = None): 
    arcpy.env.workspace = os.getcwd()
    try:
        raster = arcpy.Raster(source)
    except RuntimeError as e:
        if "ERROR 000732" in e.message:
            print "input raster not recognised: "+source
            return -1
    
    inSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) # assume square
    if cellSize is None:
        cellSize = inSize
    elif cellSize % inSize != 0:
        print "Please use a cell size which is muliple of the original: "+str(inSize)
        return -1
    else:
        print "Resampling to a cell size of {0}".format(cellSize)
        arcpy.Resample_management(raster, "resampled", str(cellSize), "CUBIC")
        raster = arcpy.Raster("resampled")
        
    e = raster.extent
    # need utm coordinates to relate to path data
    if "WGS_1984_UTM_Zone" not in e.spatialReference.name:
        project(raster,determineSystem(sampleLat,sampleLon))
        raster = arcpy.Raster("projected")
        e = raster.extent
    
    x0,y0,x1,y1 = e.XMin, e.YMin, e.XMax, e.YMax
    corner = e.lowerLeft
    arcpy.env.snapRaster = raster

    arcpy.CheckOutExtension("3D")
    slope = arcpy.sa.Slope(raster)
    aspect = arcpy.sa.Aspect(raster)
    arcpy.CheckInExtension("3D")


    raster = arcpy.RasterToNumPyArray(raster,corner,nodata_to_value=_NODATA)
    slope = arcpy.RasterToNumPyArray(slope,corner,nodata_to_value=_NODATA)
    aspect = arcpy.RasterToNumPyArray(aspect,corner,nodata_to_value=_NODATA)

    line = str(x0)+","+str(y0)+","+str(cellSize)
    np.savez_compressed("maps.npz",heightmap=raster,aspect=aspect,slope=slope)
    with open("info.txt","w") as f:
        f.write(line)
    return 0
