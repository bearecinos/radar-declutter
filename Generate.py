import arcpy
import math
import numpy as np
import os
from time import clock

_pointx, _pointy = 469900.0, 3095000.0
_above_ground, _isOffset = 100.0, True
_PATH, _antennaDir = "point1", None

## FOR USE CHECKING TIMINGS
t1 = 0
t2 = 0
t3 = 0
t4 = 0
current = 0

def setPoint(x,y,name,above=100,isOffset=True, antennaDir=None):
    """Enter the coordinates of the point and filename. Optionally set the height above the ground."""
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
    global _raster, _aspectRaster, _slopeRaster, _onesRaster, _SetupRun
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("WGS 1984 UTM Zone 45N") #  for use defining coordinates

    _raster = arcpy.Raster("default.gdb/srtm30_dem_utm45n_crop")
    _slopeRaster = arcpy.Raster("default.gdb/Slope_srtm")
    _aspectRaster = arcpy.Raster("default.gdb/Aspect_srtm")
    # No longer needed _onesRaster = arcpy.Raster("default.gdb/srtm30_dem_utm45n_crop_ones")
    
    arcpy.env.snapRaster = _raster
    os.chdir("C:\\Users\\dboyle\\OneDrive\\generation") # sometimes changes

    arcpy.CheckOutExtension("3D")
    arcpy.CheckOutExtension("Spatial")
    _SetupRun = True

def _crop(name="visibletemp"):
    """Calculates the corner of the smaller raster to use to generate data."""
    global _cropLeft, _cropLow
    rasterData = arcpy.Raster(name).extent
    low = rasterData.YMin
    left = rasterData.XMin

    pointCoords = [int((_pointx-left)/30),int((_pointy-low)/30)]

    _cropLeft = left + (pointCoords[0]-_cropWidth/2)*30
    _cropLow = low + (pointCoords[1]-_cropHeight/2)*30

# creates 'viewshed' for radar location
def _createPoint():
    """Adds a POINT Featureclass and adds a point at the selected coordinates."""
    arcpy.CreateFeatureclass_management("\\","viewshed","POINT")

    cur = arcpy.InsertCursor("viewshed.shp") # allow write access
    p = arcpy.Point(_pointx,_pointy) # x and y coordinates
    row = cur.newRow()
    row.setValue("Shape",p) # set only editable field of row
    cur.insertRow(row)

# creates npArray of compass direction to all visible points
def _getDirections():
    """Calculates the bearing from the radar to each point on the raster."""
    global _directions
    _directions = np.full((_cropHeight,_cropWidth),-1,"float32")
    # not yet tested this ndindex works, can remove old version once checked once
    for x,y in np.ndindex(_cropWidth,_cropHeight):
        if _vis[y][x] != -1:
            _directions[y][x] = math.degrees(math.atan2(_pointx-_cropLeft-30.0*x, _pointy-_cropLow-(_cropHeight-y-1)*30))+180

def _mcos(x):
    return math.cos(math.radians(x))
def _msin(x):
    return math.sin(math.radians(x))

def getAngle(vis,height,trueDist,dist,direction,slope,aspect):
    """Calcuates the angle of incidence given single values."""
    if vis == -1:
        return -1
    else:
        cosTheta = float(_elevation-height)/trueDist # sure about sign?
        sinTheta = dist/trueDist
        cosAng = cosTheta * _mcos(slope) - sinTheta * _msin(slope) * _mcos(direction-aspect)
        return  math.degrees(math.acos(cosAng)) # sometimes invalid, assuming FP error

def getDistance(vis,distance,height):
    if vis == -1:
        return -1.0
    return (distance**2+(height-_elevation)**2)**0.5

# generates map of incidence angle for all visible surfaces
def _makeIncidence():
    """Generates the incidence angle for each point of the raster."""
    vectorized = np.vectorize(getAngle)
    incidence = vectorized(_vis,_heightmap,_trueDist,_distances,_directions,_slope,_aspect)
    inc = arcpy.NumPyArrayToRaster(incidence,arcpy.Point(_cropLeft,_cropLow),30,30,-1)
    inc.save("incidence")

# generates map of 3d distance to visible surfaces from 2d distance + height difference
def _makeDistance():
    """Generates the 3D distance from the radar to each point of the raster."""
    global _trueDist
    _trueDist = np.full_like(_vis,-1,"float32") # use float64 instead??
    _trueDist[_vis==1] = (_distances[_vis==1]**2+(_heightmap[_vis==1]-_elevation)**2)**0.5
    #vectorized = np.vectorize(getDistance)
    #_trueDist = vectorized(_vis,_distances,_heightmap)
    trued = arcpy.NumPyArrayToRaster(_trueDist,arcpy.Point(_cropLeft,_cropLow),30,30,-1)
    trued.save("distance")

def _makeDist2D():
    distances = np.full_like(_vis,-1,"float32")
    for x,y in np.ndindex(_cropWidth,_cropHeight):
        if _vis[y][x] != -1:  
            distances[y][x] = ((_pointx - _cropLeft - _CellSize*x)**2 + (_pointy - _cropLow - _CellSize*y)**2)**0.5
    return distances

def _makeAntenna():
    """Find theta from horizon and phi from antenna direction (as bearing)"""
    theta = np.full_like(_vis,-1,"float32")
    phi = np.full_like(_vis,-1,"float32")
    # theta[m] = _toDegrees*np.arcsin((_heightmap[m]-_elevation)/_trueDist[m])
    # phi[m] = (_directions[m]-_antennaDir)%360
    for i,j in np.ndindex(_cropWidth,_cropHeight):
        if _vis[j][i] == 1:
            theta[j][i] = math.degrees(math.asin((_heightmap[j][i]-_elevation)/_trueDist[j][i]))
            phi[j][i] = (_directions[j][i]-_antennaDir)%360
    theta = arcpy.NumPyArrayToRaster(theta,arcpy.Point(_cropLeft,_cropLow),30,30,-1)
    theta.save("antennaTheta")
    phi = arcpy.NumPyArrayToRaster(phi,arcpy.Point(_cropLeft,_cropLow),30,30,-1)
    phi.save("antennaPhi")

def generateMaps():
    """Produces the visibility, distance and incidence rasters for the set point."""
    global _vis,_distances,_slope,_aspect,_heightmap, _elevation, _isOffset, _above_ground, t1, t2, t3, t4, current
    if not _SetupRun:
        Setup()

    os.makedirs(_PATH)
    arcpy.env.workspace = os.getcwd() + "\\" + _PATH
    
    _createPoint()
    current = clock()
    
    # elevation can be given absolutely or relative to ground
    if _isOffset:
        visibletemp = arcpy.Viewshed2_3d(_raster,"viewshed.shp","visibletemp","",
                                 "OBSERVERS","","viewtable","","","",str(_above_ground)+" Meters",
                                 "","",str(_RANGE)+" Meters","3D","","","","","ALL_SIGHTLINES")
    else:
        visibletemp = arcpy.Viewshed2_3d(_raster,"viewshed.shp","visibletemp","",
                                 "OBSERVERS","","viewtable","","",str(_above_ground)+" Meters","",
                                 "","",str(_RANGE)+" Meters","3D","","","","","ALL_SIGHTLINES")
        
    t1 += clock()-current
    current = clock()
    
    _crop("visibletemp")
    corner = arcpy.Point(_cropLeft,_cropLow)

    arcpy.Clip_management("visibletemp",
        str(_cropLeft)+" "+str(_cropLow)+" "+str(_cropLeft+30*_cropWidth)+" "+str(_cropLow+30*_cropHeight),"visible")
    _vis = arcpy.RasterToNumPyArray("visible",corner,_cropWidth,_cropHeight, -1)

    _heightmap = arcpy.RasterToNumPyArray(_raster,corner,_cropWidth,_cropHeight, -1)

    t2 += clock() - current
    current = clock()


    pointCoords = [int((_pointx-_cropLeft)/30),int((_pointy-_cropLow)/30)]
    if _isOffset:
        _elevation = _heightmap[_cropHeight-pointCoords[1],pointCoords[0]]+_above_ground
    else:
        _elevation = _above_ground
        _above_ground = _elevation - _heightmap[_cropHeight-pointCoords[1],pointCoords[0]]
        _isOffset = True

    # ArcMap approach much slower and appears to create artefacts
    # distances = arcpy.sa.CostDistance("viewshed.shp",_onesRaster)
    #_distances = arcpy.RasterToNumPyArray(distances,corner,_cropWidth,_cropHeight, -1)
    _distances = _makeDist2D()
    

    t3 += clock() - current
    current = clock()

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

    t4 += clock() - current

def finish():
    """Check the chosen extensions back in (done by default when program ends."""
    arcpy.CheckInExtension("3D")
    arcpy.CheckInExtension("Spatial")
