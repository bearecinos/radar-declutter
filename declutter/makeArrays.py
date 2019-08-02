"""Takes a raster which can be viewed in arcMap and produces a slope and aspect
array. The original raster and two new arrays are all saved in numpy format in
'maps.hdf5'."""
import arcpy
import numpy as np
from projection import *
import os
import h5py
from errors import RasterError
__all__ = ["justAspect","justSlope","justHeightmap","makeAll","rastersToNumpy"]

_NODATA = np.nan

def loadRaster(source):
    """Takes the name of a raster and returns the arcpy object.
    Returns custom RasterError if an issue occurs when loading the raster.
    If the input is already a raster, this is returned directly."""
    if type(source) == arcpy.Raster:
        return source
    arcpy.env.workspace = os.getcwd()
    try:
        return arcpy.Raster(source)
    except RuntimeError as e:
        if "ERROR 000732" in e.message:
            raise RasterError("input raster not recognised: "+source)

def resample(raster,cellSize,outDir=None):
    """Returns the original raster or a resampled version if required. Size to sample
    at must be a multiple of the original cell size. Always saved as 'resampled'.

    Parameters
    raster Raster : Original raster to resample.
    cellSize float : New cell size in metres.
    outDir string : Directory to save resampled raster in. Not used if raster already correct size.

    Returns
    Original raster if already correct size. Resampled raster if not correct size. Raises RasterError
    if cellSize is not multiple of original size."""
    if cellSize is None:
        return raster
    if outDir is None:
        out = "resampled"
    else:
        out = outDir+"/resampled"
    # assume square cells    
    inSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0))
    if cellSize % inSize != 0:
        raise RasterError("Please use a cell size which is muliple of the original: "+str(inSize))
    elif cellSize != inSize:
        print "Resampling to a cell size of {0}".format(cellSize)
        arcpy.Resample_management(raster, out, str(cellSize), "CUBIC")
        return arcpy.Raster(out)
    else:
        return raster

def coordinateSystem(raster,sampleLat,sampleLon,outDir=None):
    """Takes the name of a raster and returns the same, unless the coordinate system
    needs changing, in which case the name of the projected raster is returned.
    Always saved as 'projected'. Projection is to a UTM zone or UPS.

    Parameters
    raster Raster : Original raster to project.
    sampleLat float : Latitude of a point in the zone to project to.
    sampleLon float : Longitude of a point in the zone to project to.
    outDir string : Directory to save resampled raster in. Not used if raster doesn't need projecting.

    Returns
    Original raster if already in correct coordinate system, and the projected raster if not."""
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


def justAspect(source,sampleLat,sampleLon,cellSize=None,outDir=None):
    """Stores a numpy array of the aspect of the input raster.

    Parameters
    source string/Raster : name/instance of raster to generate array for.
    sampleLat float : Latitude of a point in the zone to project to.
    sampleLon float : Longitude of a point in the zone to project to.
    cellSize float (optional) : the cell size used to produce the arrays. Default option
        keeps the original size. Must be a multiple of the original size.
    outDir string (optional) : the directory to place projected or resampled rasters in
        if created. By default, will be placed in current directory.

    Returns
    0 if successful, otherwise -1."""
    arcpy.env.overwriteOutput = True
    try:
        raster = coordinateSystem(resample(loadRaster(source),cellSize,outDir),sampleLat,sampleLon,outDir)
    except RasterError as e:
        print "Error loading raster : "+e.message
        return -1    
    arcpy.env.snapRater = heightmap
    corner = raster.extent.lowerLeft
    cellSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) 
    aspect = arcpy.RasterToNumPyArray(_makeAspect(raster),corner,nodata_to_value=_NODATA)

    try:
        with h5py.File("maps.hdf5","a") as f:
            f["aspect"] = aspect
            f["meta"] = np.array([corner.X,corner.Y,cellSize])
    except IOError as e:
        print "Error appending to 'maps.hdf5' : "+e.message
        return -1
    return 0
    

def justSlope(source,sampleLat,sampleLon,cellSize=None,outDir=None):
    """Stores a numpy array of the slope of the input raster.

    Parameters
    source string/Raster : name/instance of raster to generate array for.
    sampleLat float : Latitude of a point in the zone to project to.
    sampleLon float : Longitude of a point in the zone to project to.
    cellSize float (optional) : the cell size used to produce the arrays. Default option
        keeps the original size. Must be a multiple of the original size.
    outDir string (optional) : the directory to place projected or resampled rasters in
        if created. By default, will be placed in current directory.

    Returns
    0 if successful, otherwise -1."""
    arcpy.env.overwriteOutput = True
    try:
        raster = coordinateSystem(resample(loadRaster(source),cellSize,outDir),sampleLat,sampleLon,outDir)
    except RasterError as e:
        print "Error loading raster : "+e.message
        return -1
    
    arcpy.env.snapRater = heightmap
    corner = raster.extent.lowerLeft
    cellSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) 
    slope = arcpy.RasterToNumPyArray(_makeSlope(raster),corner,nodata_to_value=_NODATA)

    try:
        with h5py.File("maps.hdf5","a") as f:
            f["slope"] = slope
            f["meta"] = np.array([corner.X,corner.Y,cellSize])
    except IOError as e:
        print "Error appending to 'maps.hdf5' : "+e.message
        return -1
    return 0

def justHeightmap(source,sampleLat,sampleLon,cellSize=None,outDir=None):
    """Stores the input raster in numpy format.

    Parameters
    source string/Raster : name/instance of raster to generate array for.
    sampleLat float : Latitude of a point in the zone to project to.
    sampleLon float : Longitude of a point in the zone to project to.
    cellSize float (optional) : the cell size used to produce the arrays. Default option
        keeps the original size. Must be a multiple of the original size.
    outDir string (optional) : the directory to place projected or resampled rasters in
        if created. By default, will be placed in current directory.

    Returns
    0 if successful, otherwise -1."""
    arcpy.env.overwriteOutput = True
    try:
        raster = coordinateSystem(resample(loadRaster(source),cellSize,outDir),sampleLat,sampleLon,outDir)
    except RasterError as e:
        print "Error loading raster : "+e.message
        return -1
    
    corner = raster.extent.lowerLeft
    cellSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) 
    heightmap = arcpy.RasterToNumPyArray(heightmap,corner,nodata_to_value=_NODATA)

    try:
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
    except IOError as e:
        print "Error appending to 'maps.hdf5' : "+e.message
        return -1
    return 0


def makeAll(source,sampleLat,sampleLon,cellSize = None,outDir = None):
    """Takes a raster and generates numpy arrays for that raster as well as its
    slope and aspect. Latitude and longitude values are needed to project the raster
    to an appropriate coordinate system, either a UTM zone or UPS.
    
    Parameters
    source string/Raster : name/instance of raster to generate arrays for.
    sampleLat float : Latitude of a point in the zone to project to.
    sampleLon float : Longitude of a point in the zone to project to.
    cellSize float (optional) : the cell size used to produce the arrays. Default option
        keeps the original size. Must be a multiple of the original size.
    outDir string (optional) : the directory to place projected or resampled rasters in
        if created. By default, will be placed in current directory.

    Returns
    0 if successful, otherwise -1."""
    if outDir is not None and not os.path.exists(outDir):
        os.makedirs(outDir)

    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = os.getcwd()
    # projected and resampled if required
    try:
        raster = coordinateSystem(resample(loadRaster(source),cellSize,outDir),sampleLat,sampleLon,outDir)
        arcpy.env.snapRater = raster
        corner = raster.extent.lowerLeft
        cellSize = float(arcpy.GetRasterProperties_management(raster,"CELLSIZEX").getOutput(0)) 
        aspect = arcpy.RasterToNumPyArray(_makeAspect(raster),corner,nodata_to_value=_NODATA)
        slope = arcpy.RasterToNumPyArray(_makeSlope(raster),corner,nodata_to_value=_NODATA)
        heightmap = arcpy.RasterToNumPyArray(raster,corner,nodata_to_value=_NODATA)
        # may move to pointData - after viewshed called
        heightmap[np.isnan(slope) | np.isnan(aspect)] = _NODATA

        with h5py.File("maps.hdf5","w") as f:
            f["heightmap"] = heightmap
            f["aspect"] = aspect
            f["slope"] = slope
            f["meta"] = np.array([corner.X,corner.Y,cellSize])
    except (RasterError,IOError) as e:
        print "Error making arrays for 'maps.hdf5' : "+e.message
        return -1    
    return 0

def rastersToHdf(heightmap,aspect=None,slope=None):
    """Converts all rasters to numpy arrays, generating them if they don't already exist.
    Note: This method will not project or resample the rasters given.

    Parameters
    heightmap string/Raster : name/instance of raster to generate arrays for.
    aspect string/Raster (optional) : name/instance of aspect raster. Will be generated from 'heightmap' by default.
    slope string/Raster (optional) : name/instance of slope raster. Will be generated from 'heightmap' by default.

    Returns
    0 if successful, -1 otherwise."""
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = os.getcwd()
    try:
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
    except RasterError as e:
        print "Error loading raster : " +e.message
        return -1
    
    heightmap = arcpy.RasterToNumPyArray(heightmap,corner,nodata_to_value=_NODATA)
    heightmap[np.isnan(aspect)] = _NODATA

    try:
        with h5py.File("maps.hdf5","w") as f:
            f["heightmap"] = heightmap
            f["aspect"] = aspect
            f["slope"] = slope
            f["meta"] = np.array([corner.X,corner.Y,cellSize])
    except IOError:
        print "Error: Could not write arrays to 'maps.hdf5'"
        return -1
    return 0
