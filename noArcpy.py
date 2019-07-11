# Expect to find 2 files
# maps.npz with heightmap, aspect and scope in
# info.txt with data such as corner coordinates and cell size
# IN FUTURE: could store other data in info.txt such as for modelling
import math
import numpy as np
import os
import viewshed
from time import clock

_pointx, _pointy = 469900.0, 3095000.0
_above_ground, _isOffset = 100.0, True
_PATH, _antennaDir = "point2", None

_cropWidth, _cropHeight = 210, 210
_RANGE = 3000
_SetupRun = False
_CellSize = 30.0

_save = {"visible":None,"distance":None,"incidence":None,"antennaTheta":None,"antennaPhi":None}

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

def Setup(maps="maps.npz",info="info.txt"):
    """Loads the full numpy arrays to be cropped for each point."""
    global _fullHeightmap, _fullSlope, _fullAspect , _SetupRun, low, left, _CellSize

    arrays = np.load(maps)
    _fullHeightmap = arrays["heightmap"].astype(float)
    _fullSlope = arrays["slope"].astype(float)
    _fullAspect = arrays["aspect"].astype(float)

    with open(info,"r") as f:
        data = f.read().split(",")
        left = float(data[0])
        low = float(data[1])
        _CellSize = float(data[2])
    _SetupRun = True

# creates npArray of compass direction to all visible points
def _getDirections():
    """Calculates the bearing from the radar to each visible point on the raster."""
    global _directions
    _directions = np.full((_cropHeight,_cropWidth),-1,float)
    y,x = np.indices([_cropHeight,_cropWidth])
    _directions[_vis==1] = 180.0/math.pi * np.arctan2(_pointx-_cropLeft-30.0*x[_vis==1], _pointy-_cropLow-(_cropHeight-y[_vis==1]-1)*30)+180

def _mcos(x):
    return math.cos(math.radians(x))
def _msin(x):
    return math.sin(math.radians(x))

# generates map of incidence angle for all visible surfaces
def _makeIncidence():
    global _save
    """Generates the incidence angle from the radar for each point of the raster. Final array is in radians."""
    incidence = np.full_like(_vis,-1.0,float)
    m = (_vis==1)
    cosTheta = (_elevation-_heightmap[m])/_trueDist[m]
    sinTheta = _distances[m]/_trueDist[m]
    cosAng = cosTheta * np.cos(math.pi/180.0*_slope[m]) - sinTheta * np.sin(np.pi/180.0*_slope[m]) * np.cos(math.pi/180.0*(_directions[m]-_aspect[m]))
    
    incidence[m] = np.arccos(cosAng) # if errors - caused by FP errors in cosAng
    _save["incidence"] = incidence


# generates map of 3d distance to visible surfaces from 2d distance + height difference
def _makeDistance():
    """Generates the 3D distance from the radar to each point of the raster."""
    global _trueDist, _save
    _trueDist = np.full_like(_vis,-1,float)
    _trueDist[_vis==1] = (_distances[_vis==1]**2+(_heightmap[_vis==1]-_elevation)**2)**0.5

    _save["distance"] = _trueDist


def _makeDist2D():
    """Calculates the xy distance to each visible point on the raster, ignoring differences in height."""
    distances = np.full_like(_vis,-1,float)

    y,x = np.indices([_cropHeight,_cropWidth])
    distances[_vis==1] = ((_pointx - _cropLeft - _CellSize*x[_vis==1])**2 + (_pointy - _cropLow - _CellSize*y[_vis==1])**2)**0.5
    
    return distances

def _makeAntenna():
    global _save
    """Generate theta, the pitch, and phi, bearing relative to antenna direction, of each point on raster. Result is in radians."""
    theta = np.full_like(_vis,-1,float)
    phi = np.full_like(_vis,-1,float)

    m = (_vis==1)
    # Both use radians
    theta[m] = np.arcsin((_heightmap[m]-_elevation)/_trueDist[m])
    phi[m] = ((_directions[m]-_antennaDir)%360)*math.pi/180.0

    _save["antennaTheta"] = theta
    _save["antennaPhi"] = phi

def crop(ar):
    """Crops numpy array to just grid of size cropWidth x cropHeight, centered around point."""
    # cropping - numpy arrays formatted so that printing gives map view i.e.
    # [[(0,1),  (1,1)],
    #  [(0,0),  (1,0)]]
    # so point (x,y) actually at (height-1-y,x)
    # yCoord = int(_pointy-ex.YMin)/30
    # xCoord = int(_pointx-ex.XMin)/30
    # array[height-coord-width/2:height-coord+width/2,coord-width/2:coord+width/2]
    pointCoords = [int(_pointx-left)/30,int(_pointy-low)/30]
    height = ar.shape[0]
    return ar[height-pointCoords[1]-_cropHeight/2:height-pointCoords[1]+_cropHeight/2,
                                pointCoords[0]-_cropWidth/2:pointCoords[0]+_cropWidth/2]
total = 0
part = 0
def generateMaps():
    """Produces rasters for which points are visible, their distance and incidence angles to the radar, and optionally antenna orientation data."""
    global _vis,_distances,_slope,_aspect,_heightmap, _elevation, _isOffset, _above_ground,  _cropLeft, _cropLow, total, part
    t = clock()
    if not _SetupRun:
        Setup()
        
    os.makedirs(_PATH)

    # rounds down the corner
    pointCoords = [int(_pointx-left)/30,int(_pointy-low)/30]
    _cropLeft = left + (pointCoords[0]-_cropWidth/2)*30
    _cropLow = low + (pointCoords[1]-_cropHeight/2)*30

    _heightmap = crop(_fullHeightmap)
    _slope = crop(_fullSlope)
    _aspect = crop(_fullAspect)

    c = clock()
    _vis = viewshed.viewshed(_heightmap,(_pointx-_cropLeft)/30.0, (_pointy-_cropLow)/30.0,_above_ground,_isOffset)
    part += clock() - c
    
    pointCoords = [int(_pointx-_cropLeft)/30,int(_pointy-_cropLow)/30]
    _save["visible"] = _vis
    if _isOffset:
        _elevation = _heightmap[_cropHeight-pointCoords[1],pointCoords[0]]+_above_ground
    else:
        _elevation = _above_ground
        _above_ground = _elevation - _heightmap[_cropHeight-pointCoords[1],pointCoords[0]]
        _isOffset = True

    _distances = _makeDist2D()
    _getDirections()

    _makeDistance()
    _makeIncidence()

    if _antennaDir is not None:
        _makeAntenna()
        np.savez_compressed(_PATH+"/arrays.npz",visible=_save["visible"],distance=_save["distance"],
                            incidence=_save["incidence"],antennaTheta=_save["antennaTheta"],antennaPhi=_save["antennaPhi"]) 
    else:
        np.savez_compressed(_PATH+"/arrays.npz",visible=_save["visible"],distance=_save["distance"],
                            incidence=_save["incidence"])   

    # stores coordinates, z being against reference and elevation being above ground
    with open(_PATH+"\\x_y_z_elevation","w") as f:
        f.write(str(_pointx)+","+str(_pointy)+","+str(_elevation)+","+str(_above_ground))
    total += clock()-t

