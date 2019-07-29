import arcpy
import numpy as np
from projection import *
import os
import h5py
__all__ = ["justAspect","justSlope","justHeightmap","makeAll","rastersToNumpy"]

_NODATA = np.nan

def loadRaster(source):
    """Takes the name of a raster and returns the arcpy object."""
    if type(source) == arcpy.Raster:
        return source
    arcpy.env.workspace = os.getcwd()
    try:
        return arcpy.Raster(source)
    except RuntimeError as e:
        if "ERROR 000732" in e.message:
            print "input raster not recognised: "+source
            return -1

def resample(raster,cellSize,outDir=None):
    """Returns the original raster or a resampled version if required. Size to sample
    at must be a multiple of the original cell size. Always saved as 'resampled'"""
    if cellSize is None:
        return raster
    if outDir is None:
        out = "resampled"
    else:
        out = outDir+"/resampled"
    # assume square cells    
    inSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0))
    if cellSize % inSize != 0:
        print "Please use a cell size which is muliple of the original: "+str(inSize)
        return -1
    elif cellSize != inSize:
        print "Resampling to a cell size of {0}".format(cellSize)
        arcpy.Resample_management(raster, out, str(cellSize), "CUBIC")
        return arcpy.Raster(out)
    else:
        return raster

def coordinateSystem(raster,sampleLat,sampleLon,outDir=None):
    """Takes the name of a raster and returns the same, unless the coordinate system
    needs changing, in which case the name of the projected raster is returned.
    Always saved as 'projected'"""
    e = raster.extent
    if outDir is None:
        out = "projected"
    else:
        out = outDir+"/projected"
    # need utm coordinates to relate to path data
    if -79.5 <= sampleLat and sampleLat <= 83.5:
        if "WGS_1984_UTM_Zone" not in e.spatialReference.name:
            print "Projecting to UTM coordinate system."
            return project(raster,determineSystem(sampleLat,sampleLon),out)
    else:
        if "UPS" not in e.spatialReference.name:
            print "Projecting to UPS coordinate system."
            return project(raster,determineSystem(sampleLat,sampleLon),out)
    return raster

def _makeSlope(raster):
    """Creates and returns the slope raster for a given raster object."""
    arcpy.CheckOutExtension("3D")
    slope = arcpy.sa.Slope(raster) 
    arcpy.CheckInExtension("3D")
    return slope

def _makeAspect(raster):
    """Creates and returns the aspect raster for a given raster object."""
    arcpy.CheckOutExtension("3D")
    aspect = arcpy.sa.Aspect(raster) 
    arcpy.CheckInExtension("3D")
    return aspect

# These 3 methods exist in case there are issues with the amount of
# memory needed to create and store all 3 at once using storeAll().
def justAspect(source,sampleLat,sampleLon,cellSize=None,outDir=None):
    arcpy.env.overwriteOutput = True
    raster = coordinateSystem(resample(loadRaster(source),cellSize,outDir),sampleLat,sampleLon,outDir)
    arcpy.env.snapRater = heightmap
    corner = raster.extent.lowerLeft
    cellSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) 
    aspect = arcpy.RasterToNumPyArray(_makeAspect(raster),corner,nodata_to_value=_NODATA)
    
    with h5py.File("maps.hdf5","a") as f:
        f["aspect"] = aspect
        f["meta"] = np.array([corner.X,corner.Y,cellSize])

def justSlope(source,sampleLat,sampleLon,cellSize=None,outDir=None):
    arcpy.env.overwriteOutput = True
    raster = coordinateSystem(resample(loadRaster(source),cellSize,outDir),sampleLat,sampleLon,outDir)
    arcpy.env.snapRater = heightmap
    corner = raster.extent.lowerLeft
    cellSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) 
    slope = arcpy.RasterToNumPyArray(_makeSlope(raster),corner,nodata_to_value=_NODATA)
    
    with h5py.File("maps.hdf5","a") as f:
        f["slope"] = slope
        f["meta"] = np.array([corner.X,corner.Y,cellSize])

def justHeightmap(source,sampleLat,sampleLon,cellSize=None,outDir=None):
    arcpy.env.overwriteOutput = True
    raster = coordinateSystem(resample(loadRaster(source),cellSize,outDir),sampleLat,sampleLon,outDir)
    corner = raster.extent.lowerLeft
    cellSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) 
    heightmap = arcpy.RasterToNumPyArray(heightmap,corner,nodata_to_value=_NODATA)


    with h5py.File("maps.hdf5","a") as f:
        if "aspect" in f:
            heightmap[np.isnan(f["aspect"][()])] = np.nan
        elif "slope" in f:
            heightmap[np.isnan(f["slope"][()])] = np.nan
        else:
            print "Note: Should run after storing aspect/slope."
            print "This ensures a point with missing data in any 1 array has missing data according to all 3."
        f["heightmap"] = heightmap
        f["meta"] = np.array([corner.X,corner.Y,cellSize])


def makeAll(source,sampleLat,sampleLon,cellSize = None,outDir = None):
    """Takes the name of a raster and generates numpy arrays for that raster as well as its
    slope and aspect. Latitude and longitude values are needed to correctly project the raster
    to another coordinate system if not already in UTM zone coordinates.
    
    Parameters:
        source string : name of geodatabase to work off of.
        sampleLat float : the latitude of a point in the same UTM zone.
        sampleLon float : the longitude of a point in the same UTM zone.
        cellSize float (optional) : the cell size used to produce the arrays. Default option
            keeps the original size. Must be a multiple of the original size.
        outDir string (optional) : the directory to place projected or resampled rasters in
            if created. By default, will be placed in current directory."""
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = os.getcwd()
    # projected and resampled if required
    raster = coordinateSystem(resample(loadRaster(source),cellSize,outDir),sampleLat,sampleLon,outDir)
    arcpy.env.snapRater = raster
    corner = raster.extent.lowerLeft
    cellSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) 
    aspect = arcpy.RasterToNumPyArray(_makeAspect(raster),corner,nodata_to_value=_NODATA)
    slope = arcpy.RasterToNumPyArray(_makeSlope(raster),corner,nodata_to_value=_NODATA)
    heightmap = arcpy.RasterToNumPyArray(raster,corner,nodata_to_value=_NODATA)
    # may move to stateless - after viewshed called
    heightmap[np.isnan(slope) | np.isnan(aspect)] = _NODATA

    with h5py.File("maps.hdf5","w") as f:
        f["heightmap"] = heightmap
        f["aspect"] = aspect
        f["slope"] = slope
        f["meta"] = np.array([corner.X,corner.Y,cellSize])
    return 0

def rastersToHdf(heightmap,aspect=None,slope=None):
    """Converts all rasters to numpy arrays, generating them if they don't already exist.
    Note: This expects the strings for the names of the rasters and will not do any resampling etc.
    Sometimes get errors trying to run tools on projected/resampled rasters or accessing environment.
    If so, may have to run through steps manually or call other functions e.g. 'justSlope'."""
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = os.getcwd()
    heightmap = loadRaster(heightmap)
    arcpy.env.snapRater = heightmap
    corner = heightmap.extent.lowerLeft
    cellSize = float(arcpy.GetRasterProperties_management(heightmap,"CELLSIZEX").getOutput(0)) 

    if aspect is None:
        aspect = arcpy.RasterToNumPyArray(_makeAspect(heightmap),corner,nodata_to_value=_NODATA)
    else:
        aspect = arcpy.RasterToNumpyArray(loadRaster(aspect),corner,nodata_to_value=_NODATA)
    if slope is None:
        slope = arcpy.RasterToNumPyArray(_makeSlope(heightmap),corner,nodata_to_value=_NODATA)
    else:
        slope = arcpy.RasterToNumpyArray(loadRaster(slope),corner,nodata_to_value=_NODATA)
    heightmap = arcpy.RasterToNumPyArray(heightmap,corner,nodata_to_value=_NODATA)
    heightmap[np.isnan(aspect)] = _NODATA

    with h5py.File("maps.hdf5","w") as f:
        f["heightmap"] = heightmap
        f["aspect"] = aspect
        f["slope"] = slope
        f["meta"] = np.array([corner.X,corner.Y,cellSize])
    return 0
