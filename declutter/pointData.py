'''Uses maps.hdf5 to generate the data needed to model the radar response
from a single point.'''
import math
import numpy as np
import os
import viewshed
from time import clock
import h5py
from modelling import parameters
from version import version

__all__ = ["Setup", "generateMaps"]

_SetupRun = False
_NODATA = np.nan

def Setup(): 
    """Loads the full numpy arrays to be cropped for each point.
`   Raises IOError if maps.hdf5 cannot be loaded."""
    global _fullHeightmap, _fullSlope, _fullAspect , _SetupRun
    global low, left, _CellSize, _cropSize, _RANGE
    _RANGE = parameters.env.getMaxDist()
    
    with h5py.File('maps.hdf5',"r") as f:
        _fullHeightmap = f["heightmap"][()]
        _fullSlope = f["slope"][()]
        _fullAspect = f["aspect"][()]
        left,low,_CellSize = f["meta"][()]

    # padding to avoid rounding cropping out cells
    _cropSize = int(2.0 * _RANGE / _CellSize) + 10
    _SetupRun = True
    return 0


def _getDirections(vis,pointx,pointy,cropLeft,cropLow):
    """Calculates the bearing from the radar to each visible point
    on the raster."""
    directions = np.full_like(vis,_NODATA,float)
    y, x = np.indices(vis.shape)
    directions[vis==1] = (180.0/math.pi *
                          np.arctan2(pointx-cropLeft-_CellSize*x[vis==1],
                                     pointy-cropLow-y[vis==1]*_CellSize)+180)
    return directions


def _makeIncidence(vis,elevation,heightmap,trueDist,distances,
                   slope,aspect,directions):
    """Generates the incidence angle from the radar for each
    point of the raster. Final array is in radians."""
    incidence = np.full_like(vis,_NODATA,float)
    m = (vis==1)
    cosTheta = (elevation-heightmap[m])/trueDist[m]
    sinTheta = distances[m]/trueDist[m]
    cosAng = (cosTheta * np.cos(math.pi/180.0*slope[m]) -
              sinTheta * np.sin(np.pi/180.0*slope[m]) *
              np.cos(math.pi/180.0*(directions[m]-aspect[m])))
    incidence[m] = np.arccos(cosAng) # potential NaN due to FP errors
    return incidence


def _makeDistance(vis,distances,heightmap,elevation):
    """Generates the 3D distance from the radar to each point of the raster."""
    trueDist = np.full_like(vis,_NODATA,float)
    trueDist[vis==1] = np.hypot(distances[vis==1],heightmap[vis==1]-elevation) 
    return trueDist

def _makeDist2D(vis,pointx,pointy,cropLeft,cropLow):
    """Calculates the xy distance to each visible point on the raster,
    ignoring differences in height."""
    distances = np.full_like(vis,_NODATA,float)
    y,x = np.indices(vis.shape) 
    distances[vis==1] = (((pointx - cropLeft - _CellSize*x[vis==1])**2 +
                          (pointy - cropLow - _CellSize*y[vis==1])**2)**0.5)
    return distances

def _makeAntenna(vis,heightmap,elevation,trueDist,directions,
                 distances,antennaDir):
    """Generate theta and phi for a spherical coordinate system with
    the antenna along the z-axis. Result is in radians."""
    m = (vis==1)
    # angle from vertical down to point
    worldTheta = np.arccos((heightmap[m]-elevation)/trueDist[m])
    # bearing difference converted to radians
    worldPhi = ((directions[m]-antennaDir)%360)*math.pi/180.0
    
    theta = np.full_like(vis,_NODATA,float)
    phi = np.full_like(vis,_NODATA,float)
    theta[m] = np.arccos(np.sin(worldTheta)*np.cos(worldPhi))
    phi[m] = np.arctan2(distances[m]*np.sin(worldPhi),heightmap[m]-elevation)

    return theta,phi

def generateMaps(pointx,pointy,elevation=100.0,antennaDir=None):
    """Produces rasters for which points are visible, their distance and
    incidence angles to the radar, and optionally antenna orientation data.
    
    Parameters
    ----------
    pointx - float : x-coordinate of point.
    pointy - float : y-coordinate of point.
    elevation - float (optional) : Elevation of the radar, default = 100.0.
    antennaDir - float (optional) : The direction the radar was facing in
        degrees. By default, not used.    

    Returns
    -------
    vis - 2D float array : An array of which points on the surface around
        the radar are visible or not.
    corner - [float, float] : Coordinates of the corner of the vis array.
    trueDist - float array : the distance in metres to every visible point.
    incidence - float array : the surface incidence angle (in radians) of every
        visible point.
    theta, phi - float arrays : the direction to each visible point in spherical
        coordinates, with the ends of the antenna being the poles. These will
        be None if antennaDir was not given.
    elevation - the elevation of the radar.

    Only vis is a 2D array, all others arrays are 1D and contain information
    for just the visible points.

    Note: If the elevation of the ground is undefined directly below the radar,
        vis is set to all 0, empty arrays are returned for the others, and
        above_ground is returned as the elevation, regardless of isOffset.

    Raises IOError if 'maps.hdf5' cannot be loaded.
    """
    if not _SetupRun:
        Setup()
    
    # rounds down the corner
    pointCoords = [int((pointx-left)/_CellSize),int((pointy-low)/_CellSize)]
    
    height, width = _fullHeightmap.shape
    lower = max(0,pointCoords[1]-_cropSize/2)
    upper = min(height,pointCoords[1]+_cropSize/2)
    l = max(0,pointCoords[0]-_cropSize/2)
    r = min(width,pointCoords[0]+_cropSize/2)

    cropLeft = left + l*_CellSize
    cropLow = low + lower*_CellSize
    
    heightmap = _fullHeightmap[lower:upper,l:r]
    slope = _fullSlope[lower:upper,l:r]
    aspect = _fullAspect[lower:upper,l:r]
    cropHeight,cropWidth = heightmap.shape
    
    pointCoords = np.array([(pointx-cropLeft)/_CellSize,
                            (pointy-cropLow)/_CellSize])
    groundHeight = viewshed.quadHeight(heightmap,np.array([pointCoords[0]]),
                                       np.array([pointCoords[1]]))[0]
    
    theta,phi = None, None
    
    if np.isnan(groundHeight): # Assumes use of NaN for noData, can't use ==.
        # not above mapped ground so don't generate anything,
        # as if nothing visible
        if antennaDir is not None:
            theta, phi = np.empty(0), np.empty(0)
        return (np.full(heightmap.shape,0,int),(cropLeft,cropLow), np.empty(0),
                np.empty(0), theta, phi, elevation)
    
    ys, xs = np.indices(heightmap.shape,float)
    xs = cropLeft + xs*_CellSize
    ys = cropLow + ys*_CellSize
    # avoid generating data where not needed, leave noData points
    mask = (((xs-pointx)**2 + (ys-pointy)**2) <= _RANGE**2) & ~np.isnan(slope)
    
    distances = _makeDist2D(mask,pointx,pointy,cropLeft,cropLow)
    directions = _getDirections(mask,pointx,pointy,cropLeft,cropLow)
    
    trueDist = _makeDistance(mask,distances,heightmap,elevation)
    
    incidence = _makeIncidence(mask,elevation,heightmap,trueDist,distances,
                               slope,aspect,directions)
    mask = incidence <= np.pi/2 # font facing cells only, also ignores NaNs

    # actually visible cells
    vis = viewshed.viewshed(heightmap,(pointx-cropLeft)/_CellSize,
                            (pointy-cropLow)/_CellSize,mask,
                            elevation,gridsize=_CellSize) # returns bool array
    
    if antennaDir is not None:
        theta, phi = _makeAntenna(vis,heightmap,elevation,trueDist,directions,
                                  distances,antennaDir)
        theta, phi = theta[vis], phi[vis]
    trueDist, incidence = trueDist[vis], incidence[vis]
    return vis,(cropLeft,cropLow),trueDist,incidence,theta,phi,elevation


def store(path,dist,incidence,x,y,elevation,vis=None,visCorner = None,
          antennaDir=None,theta=None,phi=None):
    """Saves the data generated for each point in a .hdf5 file.

    Parameters
    ----------
    path - string : Name of the file to save the data in. If it does not
        end in .hdf5 already, this will be appended to the file name.
    dist - float array : the distance in metres to every visible point.
    incidence - float array : the surface incidence angle (in radians) of
        every visible point.
    x,y - floats : The coordinates of the point the data was generated for.
    elevation - float : The elevation of the point the data was generated for.
    vis - 2D int array (optional) : An array of which points on the surface
        around the radar were visible or not. By default this is None,
        in which case it is not saved.
    antennaDir - float (optional) : The bearing of the radar antenna. By
        default this is None so is not saved, in which case theta and phi
        are also not saved.
    theta, phi - float arrays : the direction to each visible point in spherical
        coordinates, with the ends of the antenna being the poles. By default
        these are None. These should be None iff antennaDir is None.

    Returns
    -------
    0 if successful, otherwise -1.
    """
    if not path[-4:] == ".hdf5":
        path = path+".hdf5"
    try:
        with h5py.File(path,"w") as f:
            if vis is not None:
                f.create_dataset("visible",compression="gzip",data = vis)
                f["corner"] = visCorner
            f.create_dataset("distance",compression="gzip",data = dist)
            f.create_dataset("incidence",compression="gzip",data = incidence)
            if antennaDir is not None:
                f.create_dataset("antennaTheta",compression="gzip",data = theta)
                f.create_dataset("antennaPhi",compression="gzip",data = phi)
            f["meta"] = np.array([x,y,elevation,antennaDir])
            f["version"] = version
    except IOError as e:
        print "Could not write to hdf5 file : "+e.message
        return -1
    return 0
    
