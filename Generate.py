import arcpy
import math
import numpy as np
import os
from time import clock

_pointx, _pointy = 469900.0, 3095000.0
_above_ground, _isOffset = 100.0, True
_PATH, _antennaDir = "point1", None

def setPoint(x,y,name,above=100.0,isOffset=True, antennaDir=None):
    """ Preparation for generating the data for a particular point.
    Parameters:
    x float : x-coordinate of point.
    y float : y-coordiante of point.
    name string : name of folder to save data in.
    above float (optional) : Either actual altitude of point or elevation above the ground. Default = 100.
    isOffset boolean (optional) : Whether the given 'above' value is an offset from the ground or not. Default = True.
    antennaDir float (optional) : Bearing radar is facing in degrees. By default, not used."""
    #"""Enter the coordinates of the point and filename. Optionally set the height above the ground."""
    global _pointx,_pointy,_above_ground,_PATH, _isOffset
    _pointx = x
    _pointy = y
    _above_ground, _isOffset = above, isOffset
    _PATH, _antennaDir = name, antennaDir


def getPoint():
    """Returns the currently set point to generate data for."""
    return {"x":_pointx,"y":_pointy,"above_ground":_above_ground, "isOffset":_isOffset,"antennaDir":_antennaDir,"Filename":_PATH}

_cropWidth, _cropHeight = 210, 210
_RANGE = 3000
_SetupRun = False
_CellSize = 30.0

def Setup():
    """Performs most of the initialisation needed by arcpy."""
    global _raster, _aspectRaster, _slopeRaster, _SetupRun
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("WGS 1984 UTM Zone 45N") #  for use defining coordinates

    _raster = arcpy.Raster("default.gdb/srtm30_dem_utm45n_crop")
    _slopeRaster = arcpy.Raster("default.gdb/Slope_srtm")
    _aspectRaster = arcpy.Raster("default.gdb/Aspect_srtm")
    
    arcpy.env.snapRaster = _raster
    os.chdir("C:\\Users\\dboyle\\OneDrive\\generation") # sometimes changes

    arcpy.CheckOutExtension("3D")
    arcpy.CheckOutExtension("Spatial")
    _SetupRun = True

def _crop(name="visibletemp"):
    """Calculates the corner of the smaller raster to use to generate data.
    Parameters:
    name string (optional) : Name of raster to use to calculate the new corner from. Choice shouldn't have any effect on results."""
    global _cropLeft, _cropLow
    rasterData = arcpy.Raster(name).extent
    low = rasterData.YMin
    left = rasterData.XMin

    pointCoords = [int((_pointx-left)/30),int((_pointy-low)/30)]

    _cropLeft = left + (pointCoords[0]-_cropWidth/2)*30
    _cropLow = low + (pointCoords[1]-_cropHeight/2)*30

# creates 'viewshed' for radar location
def _createPoint():
    """Adds a POINT Featureclass and adds a point at the specified coordinates."""
    arcpy.CreateFeatureclass_management("\\","viewshed","POINT")

    cur = arcpy.InsertCursor("viewshed.shp") # allow write access
    p = arcpy.Point(_pointx,_pointy) # x and y coordinates
    row = cur.newRow()
    row.setValue("Shape",p) # set only editable field of row
    cur.insertRow(row)

# creates npArray of compass direction to all visible points
def _getDirections():
    """Calculates the bearing from the radar to each visible point on the raster."""
    global _directions
    _directions = np.full((_cropHeight,_cropWidth),-1,"float32")
    y,x = np.indices([_cropHeight,_cropWidth])
    _directions[_vis==1] = 180.0/math.pi * np.arctan2(_pointx-_cropLeft-30.0*x[_vis==1], _pointy-_cropLow-(_cropHeight-y[_vis==1]-1)*30)+180

def _mcos(x):
    return math.cos(math.radians(x))
def _msin(x):
    return math.sin(math.radians(x))

# generates map of incidence angle for all visible surfaces
def _makeIncidence():
    """Generates the incidence angle from the radar for each point of the raster. Final array is in radians."""
    incidence = np.full_like(_vis,-1.0,"float64")
    m = (_vis==1)
    cosTheta = (_elevation-_heightmap[m])/_trueDist[m]
    sinTheta = _distances[m]/_trueDist[m]
    cosAng = cosTheta * np.cos(math.pi/180.0*_slope[m]) - sinTheta * np.sin(np.pi/180.0*_slope[m]) * np.cos(math.pi/180.0*(_directions[m]-_aspect[m]))
    
    incidence[m] = np.arccos(cosAng) # if errors - caused by FP errors in cosAng
    # wrap in arcpy.sa.Float(...) to avoid being converted to integers
    inc = arcpy.sa.Float(arcpy.NumPyArrayToRaster(incidence,arcpy.Point(_cropLeft,_cropLow),30,30,-1.0))
    inc.save("incidence")

# generates map of 3d distance to visible surfaces from 2d distance + height difference
def _makeDistance():
    """Generates the 3D distance from the radar to each point of the raster."""
    global _trueDist
    _trueDist = np.full_like(_vis,-1,"float32") # use float64 instead??
    _trueDist[_vis==1] = (_distances[_vis==1]**2+(_heightmap[_vis==1]-_elevation)**2)**0.5
    trued = arcpy.NumPyArrayToRaster(_trueDist,arcpy.Point(_cropLeft,_cropLow),30,30,-1)
    trued.save("distance")

def _makeDist2D():
    """Calculates the xy distance to each visible point on the raster, ignoring differences in height."""
    distances = np.full_like(_vis,-1,"float32")

    y,x = np.indices([_cropHeight,_cropWidth])
    distances[_vis==1] = ((_pointx - _cropLeft - _CellSize*x[_vis==1])**2 + (_pointy - _cropLow - _CellSize*y[_vis==1])**2)**0.5
    
    return distances

def _makeAntenna():
    """Generate theta, the pitch, and phi, bearing relative to antenna direction, of each point on raster. Result is in radians."""
    theta = np.full_like(_vis,-1,"float32")
    phi = np.full_like(_vis,-1,"float32")

    m = (_vis==1)
    # Both use radians
    theta[m] = np.arcsin((_heightmap[m]-_elevation)/_trueDist[m])
    phi[m] = ((_directions[m]-_antennaDir)%360)*math.pi/180.0
    
    theta = arcpy.NumPyArrayToRaster(theta,arcpy.Point(_cropLeft,_cropLow),30,30,-1)
    theta.save("antennaTheta")
    phi = arcpy.NumPyArrayToRaster(phi,arcpy.Point(_cropLeft,_cropLow),30,30,-1)
    phi.save("antennaPhi")

def generateMaps():
    """Produces rasters for which points are visible, their distance and incidence angles to the radar, and optionally antenna orientation data."""
    global _vis,_distances,_slope,_aspect,_heightmap, _elevation, _isOffset, _above_ground
    if not _SetupRun:
        Setup()

    os.makedirs(_PATH)
    arcpy.env.workspace = os.getcwd() + "\\" + _PATH
    
    _createPoint()
    
    # elevation can be given absolutely or relative to ground
    if _isOffset:
        visibletemp = arcpy.Viewshed2_3d(_raster,"viewshed.shp","visibletemp","",
                                 "OBSERVERS","","viewtable","","","",str(_above_ground)+" Meters",
                                 "","",str(_RANGE)+" Meters","3D","","","","","ALL_SIGHTLINES")
    else:
        visibletemp = arcpy.Viewshed2_3d(_raster,"viewshed.shp","visibletemp","",
                                 "OBSERVERS","","viewtable","","",str(_above_ground)+" Meters","",
                                 "","",str(_RANGE)+" Meters","3D","","","","","ALL_SIGHTLINES")
    
    _crop("visibletemp")
    corner = arcpy.Point(_cropLeft,_cropLow)

    arcpy.Clip_management("visibletemp",
        str(_cropLeft)+" "+str(_cropLow)+" "+str(_cropLeft+30*_cropWidth)+" "+str(_cropLow+30*_cropHeight),"visible")
    _vis = arcpy.RasterToNumPyArray("visible",corner,_cropWidth,_cropHeight, -1)

    _heightmap = arcpy.RasterToNumPyArray(_raster,corner,_cropWidth,_cropHeight, -1).astype("float64")

    pointCoords = [int((_pointx-_cropLeft)/30),int((_pointy-_cropLow)/30)]
    if _isOffset:
        _elevation = _heightmap[_cropHeight-pointCoords[1],pointCoords[0]]+_above_ground
    else:
        _elevation = _above_ground
        _above_ground = _elevation - _heightmap[_cropHeight-pointCoords[1],pointCoords[0]]
        _isOffset = True

    _distances = _makeDist2D()
    _getDirections()

    _slope = arcpy.RasterToNumPyArray(_slopeRaster,corner,_cropWidth,_cropHeight, -1)
    _aspect = arcpy.RasterToNumPyArray(_aspectRaster,corner,_cropWidth,_cropHeight, -2)

    _makeDistance()
    _makeIncidence()


    if _antennaDir is not None:
        _makeAntenna()
    

    # stores coordinates, z being against reference and elevation being above ground
    with open(_PATH+"\\x_y_z_elevation","w") as f:
        f.write(str(_pointx)+","+str(_pointy)+","+str(_elevation)+","+str(_above_ground))

    # unnecessary objects deleted at end
    arcpy.Delete_management("viewtable")
    arcpy.Delete_management("visibletemp")

def finish():
    """Check the arcpy extensions used back in (done by default when program ends)."""
    arcpy.CheckInExtension("3D")
    arcpy.CheckInExtension("Spatial")
