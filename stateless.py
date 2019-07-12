# Expect to find 2 files
# maps.npz with heightmap, aspect and scope in
# info.txt with data such as corner coordinates and cell size
# IN FUTURE: could store other data in info.txt such as for modelling
import math
import numpy as np
import os
import viewshed
from time import clock

#_pointx, _pointy = 469900.0, 3095000.0
#_above_ground, _isOffset = 100.0, True
#_PATH, _antennaDir = "point2", None

_cropWidth, _cropHeight = 210, 210
_RANGE = 3000
_SetupRun = False
_CellSize = 30.0


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
def _getDirections(vis,pointx,pointy,cropLeft,cropLow):
    """Calculates the bearing from the radar to each visible point on the raster."""
    directions = np.full((_cropHeight,_cropWidth),-1,float)
    y,x = np.indices([_cropHeight,_cropWidth])
    directions[vis==1] = 180.0/math.pi * np.arctan2(pointx-cropLeft-30.0*x[vis==1], pointy-cropLow-(_cropHeight-y[vis==1]-1)*30)+180
    return directions

def _mcos(x):
    return math.cos(math.radians(x))
def _msin(x):
    return math.sin(math.radians(x))

# generates map of incidence angle for all visible surfaces
def _makeIncidence(vis,elevation,heightmap,trueDist,distances,slope,aspect,directions):
    """Generates the incidence angle from the radar for each point of the raster. Final array is in radians."""
    incidence = np.full_like(vis,-1.0,float)
    m = (vis==1)
    cosTheta = (elevation-heightmap[m])/trueDist[m]
    sinTheta = distances[m]/trueDist[m]
    cosAng = cosTheta * np.cos(math.pi/180.0*slope[m]) - sinTheta * np.sin(np.pi/180.0*slope[m]) * np.cos(math.pi/180.0*(directions[m]-aspect[m]))
    
    incidence[m] = np.arccos(cosAng) # if errors - caused by FP errors in cosAng
    return incidence


# generates map of 3d distance to visible surfaces from 2d distance + height difference
def _makeDistance(vis,distances,heightmap,elevation):
    """Generates the 3D distance from the radar to each point of the raster."""
    trueDist = np.full_like(vis,-1,float)
    trueDist[vis==1] = (distances[vis==1]**2+(heightmap[vis==1]-elevation)**2)**0.5
    return trueDist
    


def _makeDist2D(vis,pointx,pointy,cropLeft,cropLow):
    """Calculates the xy distance to each visible point on the raster, ignoring differences in height."""
    distances = np.full_like(vis,-1,float)

    y,x = np.indices([_cropHeight,_cropWidth])
    distances[vis==1] = ((pointx - cropLeft - _CellSize*x[vis==1])**2 + (pointy - cropLow - _CellSize*y[vis==1])**2)**0.5
    
    return distances

def _makeAntenna(vis,heightmap,elevation,trueDist,directions,antennaDir):
    """Generate theta, the pitch, and phi, bearing relative to antenna direction, of each point on raster. Result is in radians."""
    theta = np.full_like(vis,-1,float)
    phi = np.full_like(vis,-1,float)

    m = (vis==1)
    # Both use radians
    theta[m] = np.arcsin((heightmap[m]-elevation)/trueDist[m])
    phi[m] = ((directions[m]-antennaDir)%360)*math.pi/180.0

    return theta,phi

def generateMaps(pointx,pointy,path,above_ground=100.0,isOffset=True,antennaDir=None):
    """Produces rasters for which points are visible, their distance and incidence angles to the radar, and optionally antenna orientation data."""
    if not _SetupRun:
        print "called"
        Setup()
        
    os.makedirs(path)
    
    # rounds down the corner
    pointCoords = [int(pointx-left)/30,int(pointy-low)/30]
    cropLeft = left + (pointCoords[0]-_cropWidth/2)*30
    cropLow = low + (pointCoords[1]-_cropHeight/2)*30

    height = _fullHeightmap.shape[0]

    heightmap = _fullHeightmap[height-pointCoords[1]-_cropHeight/2:height-pointCoords[1]+_cropHeight/2,
                                pointCoords[0]-_cropWidth/2:pointCoords[0]+_cropWidth/2]
    slope = _fullSlope[height-pointCoords[1]-_cropHeight/2:height-pointCoords[1]+_cropHeight/2,
                                pointCoords[0]-_cropWidth/2:pointCoords[0]+_cropWidth/2]
    aspect = _fullAspect[height-pointCoords[1]-_cropHeight/2:height-pointCoords[1]+_cropHeight/2,
                                pointCoords[0]-_cropWidth/2:pointCoords[0]+_cropWidth/2]

    vis = viewshed.viewshed(heightmap,(pointx-cropLeft)/30.0, (pointy-cropLow)/30.0,above_ground,isOffset)
    
    
    pointCoords = [int(pointx-cropLeft)/30,int(pointy-cropLow)/30]
    if isOffset:
        elevation = heightmap[_cropHeight-pointCoords[1],pointCoords[0]]+above_ground
    else:
        elevation = above_ground
        above_ground = elevation - heightmap[_cropHeight-pointCoords[1],pointCoords[0]]
        isOffset = True

    distances = _makeDist2D(vis,pointx,pointy,cropLeft,cropLow)
    directions = _getDirections(vis,pointx,pointy,cropLeft,cropLow)

    trueDist = _makeDistance(vis,distances,heightmap,elevation)
    
    incidence = _makeIncidence(vis,elevation,heightmap,trueDist,distances,slope,aspect,directions)

    if antennaDir is not None:
        theta, phi = _makeAntenna(vis,heightmap,elevation,trueDist,directions,antennaDir)
        np.savez_compressed(path+"/arrays.npz",visible=vis,distance=trueDist,
                            incidence=incidence,antennaTheta=theta,antennaPhi=phi) 
    else:
        np.savez_compressed(path+"/arrays.npz",visible=vis,distance=trueDist,
                            incidence=incidence)   

    # stores coordinates, z being against reference and elevation being above ground
    with open(path+"/x_y_z_elevation","w") as f:
        f.write(str(pointx)+","+str(pointy)+","+str(elevation)+","+str(above_ground))
    return 0
    
if __name__=="__main__":
    Setup()
    print generateMaps(469900.0, 3095000.0,"point1")
