import arcpy
import numpy as np
from projection import *
import os
_NODATA = np.nan

# divided into several different function calls to try to reduce strain on memory
# each returns name of any raster needed for next stage
def loadRaster(source):
    arcpy.env.workspace = os.getcwd()
    try:
        return arcpy.Raster(source)
    except RuntimeError as e:
        if "ERROR 000732" in e.message:
            print "input raster not recognised: "+source
            return -1

def resample(source,cellSize):
    raster = loadRaster(source)
    
    inSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) # assume square
    if cellSize % inSize != 0:
        print "Please use a cell size which is muliple of the original: "+str(inSize)
        return -1
    elif cellSize != inSize:
        print "Resampling to a cell size of {0}".format(cellSize)
        arcpy.Resample_management(raster, "resampled", str(cellSize), "CUBIC")
        return "resampled"
    else:
        print "Already that size."
        return source

def coordinateSystem(source,sampleLat,sampleLon):
    raster = loadRaster(source)
    
    e = raster.extent
    # need utm coordinates to relate to path data
    if "WGS_1984_UTM_Zone" not in e.spatialReference.name:
        project(raster,determineSystem(sampleLat,sampleLon))
        print "Projecting to UTM coordinate system."
        return "projected"
    print "Already in UTM coordinates"
    return source

def makeSlope(source):
    raster = loadRaster(source)
    e = raster.extent
    x0,y0,x1,y1 = e.XMin, e.YMin, e.XMax, e.YMax
    corner = e.lowerLeft
    arcpy.env.snapRaster = raster
    arcpy.CheckOutExtension("3D")
    slope = arcpy.sa.Slope(raster)
    arcpy.CheckInExtension("3D")

    slope = arcpy.RasterToNumPyArray(slope,corner,nodata_to_value=_NODATA)
    save = {"slope":slope}
    if os.path.exists("maps.npz"):
        f = np.load("maps.npz")
        for name in f:
            if name != "slope":
                save[name] = f[name]
    np.savez_compressed("maps.npz",**save)  
        
def makeAspect(source):
    raster = loadRaster(source)
    e = raster.extent
    x0,y0,x1,y1 = e.XMin, e.YMin, e.XMax, e.YMax
    corner = e.lowerLeft
    arcpy.env.snapRaster = raster
    arcpy.CheckOutExtension("3D")
    aspect = arcpy.sa.Aspect(raster)
    arcpy.CheckInExtension("3D")

    aspect = arcpy.RasterToNumPyArray(aspect,corner,nodata_to_value=_NODATA)
    save = {"aspect":aspect}
    if os.path.exists("maps.npz"):
        f = np.load("maps.npz")
        for name in f:
            if name != "aspect":
                save[name] = f[name]
    np.savez_compressed("maps.npz",**save) 

# set heightmap cells to NaN if NaN in slope (or aspect)
def storeHeightmap():
    raster = loadRaster(source)
    e = raster.extent
    x0,y0,x1,y1 = e.XMin, e.YMin, e.XMax, e.YMax
    corner = e.lowerLeft
    heightmap = arcpy.RasterToNumPyArray(heightmap,corner,nodata_to_value=_NODATA)

    save = {"heightmap":heightmap}
    if os.path.exists("maps.npz"):
        f = np.load("maps.npz")
        for name in f:
            if name != "heightmap":
                save[name] = f[name]
    if "aspect" in save:
        heightmap[np.isnan(save["aspect"])] = np.nan
    elif "slope" in save:
        heightmap[np.isnan(save["slope"])] = np.nan
    np.savez_compressed("maps.npz",**save) 

# Warning: can fail due to memory error, if so, try separate parts above
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
    elif cellSize != inSize:
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
    # avoid cells which couldn't have slope/aspect calculated for being used as data
    raster[np.isnan(slope)] = _NODATA
    
    line = str(x0)+","+str(y0)+","+str(cellSize)
    np.savez_compressed("maps.npz",heightmap=raster,aspect=aspect,slope=slope)
    with open("info.txt","w") as f:
        f.write(line)
    return 0
