'''This module combines the methods used in path.py, pointData.py and
models.py to display a radargram for data without storing the files
for all points.
This is useful where a single output is needed and the data for each
point would be deleted immediately after.'''
import path
import numpy as np
import os
import models
from progress import progress
import multiprocessing as mp
import pointData
import viewshed
import matplotlib.pyplot as plt

def processData(filename,crop=[0,0],outName=None,style=None,adjusted=False,save=True):
    """Takes a gps path and displays a radargram for that path.

    Parameters
    filename string  : The name of the file to generate a radargram for.
    crop [int, int] (optional) : crop = [A,B] ignores the first A and last B points of the path.
    outName string (optional) : The name of the file to save the radargram in. By default, this is
        taken from 'filename' unless 'save' is set to False.
    style string (optional) : The format of the input file, either 'gpx', 'dst' or 'xyz'. By default,
        the loadData method determines the format from the file extension and assumes gpx if the extension
        is not recognised.
    adjusted bool (optional) : Default False. If True, slices of the radargram are adjusted by the elvation
        of the radar, meaning the response from the surface appears similar to the profile of the ground
        beneath the radar path rather than being flat.
    save bool (optional) : Default True. If True, the radargram output is saved automatically.

    Returns
    The radargram output if successful, otherwise -1.

    """
    try:
        xs, ys, zs = path.loadData(filename, crop, style)
    except IOError:
        print "Could not load data from file : "+filename
        if style is not None:
            print "Is "+style+" the correct format?"
        return -1
    if len(xs) == 0:
        return -1
    if outName is None and save:
        outName = filename[:-4]+".png"
    return _genPath(xs,ys,zs,outName,False,adjusted)

def _makeDirections(xs,ys):
    """Calculates the direction the antenna is facing at each point based on the adjacent points.

    Parameters
    xs float array : The x-coordinates of the path.
    ys float array : The y-coordinates of the path.

    Returns
    direction float array : The bearing of the antenna for each point."""
    direction = np.full_like(xs,0,float) # degrees
    direction[0] = 180.0/np.pi*(np.arctan2(xs[0]-xs[1], ys[0]-ys[1]))+180.0
    direction[-1] = 180.0/np.pi*(np.arctan2(xs[-2]-xs[-1], ys[-2]-ys[-1]))+180.0
    m = np.full_like(xs,True,bool)
    m[[0,-1]] = False
    m2 = np.roll(m,-1)
    m3 = np.roll(m,1)
    direction[m] = 180.0/np.pi*np.arctan2(xs[m2]-xs[m3], ys[m2]-ys[m3])+180.0
    return direction

def _genPath(xs,ys,zs,name,isOffset=True,adjusted=False):
    global pool
    """Displays the radargram for a path.

    Parameters
    xs float array : Array of x coordinates of path.
    ys float array : Array of y coordinates of path.
    zs float array : Array of altitude/elevation above ground along path.
    name string : Name of file to save radargram as. Not saved if name is None.
    isOffset bool (optional) : Whether the given z coordinates are altitude or relative to the ground. Default is relative.
    adjusted bool (optional) : Default False. If True, slices of the radargram are adjusted by the elvation
        of the radar, meaning the response from the surface appears similar to the profile of the ground
        beneath the radar path rather than being flat.

    Returns
    returnData 2D float array : The radargram output.
    Returns -1 if unsuccessful."""
    direction = _makeDirections(xs,ys)
    n = len(xs)

    _MAXDIST = models._MAXDIST
    _SPACE_GRANULARITY = models._SPACE_GRANULARITY
    _steps = models._steps
    reflectionModels = models.models
    titles = models.titles
    
    returnData = np.full((len(reflectionModels),n,_steps),0,float) # 3D - many plots
    plt.rcParams['axes.formatter.limits'] = [-4,4]
    plt.figure(figsize=(12,8))
    heights = []

    p = mp.Pool(mp.cpu_count())
    data = [(x,y,z,i,isOffset,angle,reflectionModels) for x,y,i,z,angle in
                            zip(xs,ys,np.arange(n),zs,direction)]
    try:
        for i, h, ars in progress(p.imap_unordered(_worker,data),n):
            returnData[:,i] = ars
            heights.append(h)
    except IOError as e:
        p.close()
        print "\nError reading 'maps.hdf5' :\n" + e.message
        return -1
    p.close()

    returnData = returnData[:,np.any(returnData[0] != 0,1)] # remove invalid rows
    n = returnData.shape[1]
    
    highest,lowest = 0,0
    if adjusted:
        highest = max(heights)
        lowest = min(heights)
        draw = np.full((len(reflectionModels),n,_steps + int((highest-lowest)/_SPACE_GRANULARITY)),0)
        for i in range(n):
            start = int((highest-heights[i])/_SPACE_GRANULARITY)
            draw[:,i,start:start+_steps] = returnData[:,i]
        returnData = draw
    
    cells = int(np.ceil(np.sqrt(len(reflectionModels))))*110
    ys = np.linspace(0,(_MAXDIST+highest-lowest)*2.0/3e8,_steps+(highest-lowest)/_SPACE_GRANULARITY)
    for j in range(len(reflectionModels)):
        plt.subplot(cells+j+1)
        plt.ylim((_MAXDIST+highest-lowest)*2.0/3e8,0)
        draw = np.swapaxes(returnData[j],0,1)
        plt.contourf(np.arange(n), ys, draw, 100,norm=models.MidNorm(np.mean(draw)), cmap="Greys")
        plt.title(titles[j])
        plt.colorbar()
    if name is not None:
        plt.savefig(name)
    plt.show()
    return returnData

def _worker(args): # x,y,z,i,offset,angle,models array
    # returns i, height, each response array.
    # can raise IOError when reading maps.hdf5 in generateMaps
    # or in pointData for pointX.hdf5
    pointx,pointy,pointz,i,isOffset,angle,reflectionModels = args

    _,dist,incidence,theta,phi,elevation = pointData.generateMaps(pointx,pointy,pointz,
                                                                  isOffset,angle)
    ars = np.full((len(reflectionModels),models._steps),0,float)
    for j in range(len(reflectionModels)):
        ars[j] = models.processSlice(dist,incidence,theta,phi,reflectionModels[j])

    return i, elevation, ars
