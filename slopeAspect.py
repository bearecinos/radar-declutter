import arcpy
import numpy as np
from projection import *
import os
_NODATA = -100000

def makeRasters(source, sampleLat, sampleLon, cellSize = None): # currently can't change stepsize
    arcpy.env.workspace = os.getcwd()
    try:
        raster = arcpy.Raster(source)
    except RuntimeError as e:
        if "ERROR 000732" in e.message:
            print "input raster not recognised: "+source
            return -1
    e = raster.extent
    # need utm coordinates to relate to path data
    if "WGS_1984_UTM_Zone" not in e.spatialReference.name:
        project(raster,determineSystem(sampleLat,sampleLon))
        raster = arcpy.Raster("projected")
        e = raster.extent
    
    inSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) # assume square
    if cellSize is None:
        cellSize = inSize
    elif cellSize % inSize != 0:
        print "Please use a cell size which is muliple of the original: "+str(inSize)
        return -1
    step = int(cellSize/inSize)
    
    x0,y0,x1,y1 = e.XMin, e.YMin, e.XMax, e.YMax
    corner = e.lowerLeft
    arcpy.env.snapRaster = raster

    arcpy.CheckOutExtension("3D")
##    slope = arcpy.sa.Slope(raster)
##    aspect = arcpy.sa.Aspect(raster)
    arcpy.sa.Slope(raster).save("slope")
    arcpy.sa.Aspect(raster).save("aspect")
    arcpy.CheckInExtension("3D")

##    raster = arcpy.RasterToNumPyArray(raster,corner,nodata_to_value=_NODATA)[::step,::step]
##    slope = arcpy.RasterToNumPyArray(slope,corner,nodata_to_value=_NODATA)[::step,::step]
##    aspect = arcpy.RasterToNumPyArray(aspect,corner,nodata_to_value=_NODATA)[::step,::step]

    line = str(x0)+","+str(y0)+","+str(cellSize)
    #np.savez_compressed("maps.npz",heightmap=raster,aspect=aspect,slope=slope)
    with open("info.txt","w") as f:
        f.write(line)
    return 0

def moveToNumpy(source):
    arcpy.env.workspace = os.getcwd()
    corner = arcpy.Raster(source).extent.lowerLeft
    raster = arcpy.RasterToNumPyArray(source,corner,nodata_to_value=_NODATA)
    slope = arcpy.RasterToNumPyArray("slope",corner,nodata_to_value=_NODATA)
    aspect = arcpy.RasterToNumPyArray("aspect",corner,nodata_to_value=_NODATA)
    np.savez_compressed("maps.npz",heightmap=raster,aspect=aspect,slope=slope)

    arcpy.Delete_management("slope")
    arcpy.Delete_management("aspect")
