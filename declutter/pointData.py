# Expect to find maps.hdf5 with surface data e.g. height/slope/aspect
import math
import numpy as np
import os
import viewshed
from time import clock
import h5py

__all__ = ["Setup", "generateMaps"]

#_pointx, _pointy = 469900.0, 3095000.0
#_above_ground, _isOffset = 100.0, True
#_PATH, _antennaDir = "point2", None

#_cropWidth, _cropHeight = 210, 210
_RANGE = 3000.0
_SetupRun = False
_NODATA = np.nan


def Setup(): 
    """Loads the full numpy arrays to be cropped for each point."""
    global _fullHeightmap, _fullSlope, _fullAspect , _SetupRun, low, left, _CellSize, _cropSize

    with h5py.File('maps.hdf5',"r") as f:
        _fullHeightmap = f["heightmap"][()]
        _fullSlope = f["slope"][()]
        _fullAspect = f["aspect"][()]
        left,low,_CellSize = f["meta"][()]
   
    _cropSize = int(2.0*_RANGE/_CellSize) + 10 # padding to avoid rounding cropping out cells
    _SetupRun = True
    return 0

# creates npArray of compass direction to all visible points
def _getDirections(vis,pointx,pointy,cropLeft,cropLow):
    """Calculates the bearing from the radar to each visible point on the raster."""
    directions = np.full_like(vis,_NODATA,float)
    y,x = np.indices(vis.shape) 
    directions[vis==1] = 180.0/math.pi * np.arctan2(pointx-cropLeft-_CellSize*x[vis==1], pointy-cropLow-(_cropSize-y[vis==1]-1)*_CellSize)+180
    return directions

# generates map of incidence angle for all visible surfaces
def _makeIncidence(vis,elevation,heightmap,trueDist,distances,slope,aspect,directions):
    """Generates the incidence angle from the radar for each point of the raster. Final array is in radians."""
    incidence = np.full_like(vis,_NODATA,float)
    m = (vis==1)
    cosTheta = (elevation-heightmap[m])/trueDist[m]
    sinTheta = distances[m]/trueDist[m]
    cosAng = cosTheta * np.cos(math.pi/180.0*slope[m]) - sinTheta * np.sin(np.pi/180.0*slope[m]) * np.cos(math.pi/180.0*(directions[m]-aspect[m]))
    incidence[m] = np.arccos(cosAng) # if errors - caused by FP errors in cosAng or NaN values
    return incidence


# generates map of 3d distance to visible surfaces from 2d distance + height difference
def _makeDistance(vis,distances,heightmap,elevation):
    """Generates the 3D distance from the radar to each point of the raster."""
    trueDist = np.full_like(vis,_NODATA,float)
    trueDist[vis==1] = (distances[vis==1]**2+(heightmap[vis==1]-elevation)**2)**0.5
    return trueDist

def _makeDist2D(vis,pointx,pointy,cropLeft,cropLow):
    """Calculates the xy distance to each visible point on the raster, ignoring differences in height."""
    distances = np.full_like(vis,_NODATA,float)

    y,x = np.indices(vis.shape) 
    distances[vis==1] = ((pointx - cropLeft - _CellSize*x[vis==1])**2 + (pointy - cropLow - _CellSize*y[vis==1])**2)**0.5
    
    return distances

def _makeAntenna(vis,heightmap,elevation,trueDist,directions,distances,antennaDir):
    # finds theta and phi for spherical system where direction antenna points is z-axis
    # expecting directivity to only depend on theta
    # assume antenna level, pointing in direction antennaDir
    """ Generate theta and phi for a spherical coordinate system with the antenna along the z-axis. Result is in radians."""
    m = (vis==1)
    worldTheta = np.arccos((heightmap[m]-elevation)/trueDist[m]) # angle from vertical down to point
    worldPhi = ((directions[m]-antennaDir)%360)*math.pi/180.0 # angle antenna would turn to point closest

    theta = np.full_like(vis,_NODATA,float)
    phi = np.full_like(vis,_NODATA,float)
    theta[m] = np.arccos(np.sin(worldTheta)*np.cos(worldPhi))
    phi[m] = np.arctan2(distances[m]*np.sin(worldPhi),heightmap[m]-elevation)

    return theta,phi

def generateMaps(pointx,pointy,above_ground=100.0,isOffset=True,antennaDir=None):
    """Produces rasters for which points are visible, their distance and incidence angles to the radar, and optionally antenna orientation data.
    Parameters:
    pointx float : x-coordinate of point.
    pointy float : y-coordinate of point.
    path string : name of folder to save data in.
    above_ground float (optional) : Either actual altitude of radar or elevation above ground. Default = 100.0.
    isOffset boolean (optional) : Indicates the given 'above_ground' is relative to the ground. Default = True.
    antennaDir float (optional) : The direction the radar was facing in degrees. By default, not used.    
    """
    
    if not _SetupRun:
        if Setup():
            return -1
    
    # rounds down the corner
    pointCoords = [int((pointx-left)/_CellSize),int((pointy-low)/_CellSize)]
    
    height = _fullHeightmap.shape[0]
    # stop at edges, end up with point no longer in center - affects coordinates of point in grid.
    # in y direction lower index is greater y position
    upper, lower = max(0,height-pointCoords[1]-_cropSize/2), min(height,height-pointCoords[1]+_cropSize/2)
    l, r = max(0,pointCoords[0]-_cropSize/2), min(height,pointCoords[0]+_cropSize/2)

    cropLeft = left + l*_CellSize
    cropLow = low + (height-1-lower)*_CellSize
    
    heightmap = _fullHeightmap[upper:lower,l:r]
    slope = _fullSlope[upper:lower,l:r]
    aspect = _fullAspect[upper:lower,l:r]
    cropHeight,cropWidth = heightmap.shape
    
    pointCoords = np.array([(pointx-cropLeft)/_CellSize,(pointy-cropLow)/_CellSize])
    groundHeight = viewshed.quadHeight(heightmap,np.array([pointCoords[0]]),np.array([cropHeight-1.0-pointCoords[1]]))[0]

    theta,phi = None, None
    
    if np.isnan(groundHeight): # Assumes use of NaN, have to use different checks to regular values
        # not above mapped ground so don't generate anything, as if nothing visible
        # TODO: replace returning above_ground as elevation with NaN then check for in other methods
        if antennaDir is not None:
            theta, phi = np.empty(0), np.empty(0)
        return np.full(heightmap.shape,0,int), np.empty(0), np.empty(0), theta, phi, above_ground
    
    if isOffset:
        elevation = groundHeight+above_ground
    else:
        elevation = above_ground
    
    ys, xs = np.indices(heightmap.shape,float)
    ys = cropHeight - 1 - ys
    xs = cropLeft + xs*_CellSize
    ys = cropLow + ys*_CellSize
    mask = (((xs-pointx)**2 + (ys-pointy)**2) <= _RANGE**2) & ~np.isnan(heightmap) ## leave noData points

    distances = _makeDist2D(mask,pointx,pointy,cropLeft,cropLow)
    directions = _getDirections(mask,pointx,pointy,cropLeft,cropLow)
    
    trueDist = _makeDistance(mask,distances,heightmap,elevation)
    
    incidence = _makeIncidence(mask,elevation,heightmap,trueDist,distances,slope,aspect,directions)
    mask = incidence <= np.pi/2 # font facing cells only, also ignores NaNs

    
    vis = viewshed.viewshed(heightmap,(pointx-cropLeft)/_CellSize, (pointy-cropLow)/_CellSize,mask,
                            elevation,False,gridsize=_CellSize) # return bool array
    
    
    if antennaDir is not None:
        theta, phi = _makeAntenna(vis,heightmap,elevation,trueDist,directions,distances,antennaDir)
        theta, phi = theta[vis], phi[vis]
    trueDist, incidence = trueDist[vis], incidence[vis]
    return vis,trueDist,incidence,theta,phi,elevation


def store(path,dist,incidence,x,y,elevation,vis=None,antennaDir=None,theta=None,phi=None):
    '''Expects all arrays except vis to have already been reduced by generateMaps()
    i.e. 1D array of only valid points.'''
    if not path[-4:] == ".hdf5":
        path = path+".hdf5"
    with h5py.File(path,"w") as f:
        if vis is not None:
            f["visible"] = vis
        f["distance"] = dist
        f["incidence"] = incidence
        if antennaDir is not None:
            f["antennaTheta"] = theta
            f["antennaPhi"] = phi
        f["meta"] = np.array([x,y,elevation,antennaDir])
    return 0
    
if __name__=="__main__" and False:
    if Setup():
        print "Failed"
    else:
        #import matplotlib.pyplot as plt
        t = clock()
        generateMaps(469900.0, 3095000.0,"timing/point0")
        # For new data:
        # x,y =  436137.27838013606, 8757751.77326004
        # elevation = 526.0
        print clock()  - t
        input()
