'''This module combines the methods used in path.py, pointData.py and
models.py to display a radargram for data without storing the files
for all points.
This is useful where a single output is needed and the data for each
point would be deleted immediately after.'''
import path
import numpy as np
import os
import radar
from progress import progress
import multiprocessing as mp
import pointData
import viewshed
import matplotlib.pyplot as plt
import align
from modelling import parameters

def processData(filename,crop=[0,0],outName=None,style=None,adjusted=False,save=True, parallel=True):
    """Takes a gps path and displays a radargram for that path.

    Parameters
    filename string  : The name of the file to generate a radargram for.
    crop [int, int] (optional) : crop = [A,B] ignores the first A and last B points of the path.
    outName string (optional) : The name of the file to save the radargram in. By default, this is
        taken from 'filename' unless 'save' is set to False.
    style string (optional) : The format of the input file, either 'gpx', 'dst' or 'xyz'. By default,
        the loadData method determines the format from the file extension and assumes gpx if the extension
        is not recognised.
    adjusted bool (optional) : Shift the data for each point to align the response from the surface
        directly beneath the radar with the top of the plot.
    save bool (optional) : Default True. If True, the radargram output is saved automatically.

    Returns
    The radargram output if successful, otherwise -1.

    """
    print parameters.env
    try:
        xs, ys, zs = path.loadData(filename, crop, style)
    except IOError:
        print "Could not load data from file : "+filename
        if style is not None:
            print "Is "+style+" the correct format?"
        return -1
    if len(xs) < 2: # not enough points for a path
        return -1
    if outName is None and save:
        outName = filename[:-4]+".png"
    return _genPath(xs,ys,zs,outName,False,adjusted, parallel)

def _makeDirections(xs,ys):
    """Calculates the direction the antenna is facing at each point based on the adjacent points.

    Parameters
    xs float array : The x-coordinates of the path.
    ys float array : The y-coordinates of the path.

    Returns
    direction float array : The bearing of the antenna for each point.

    Fails if path is less than 2 points due to index errors."""
    direction = np.full_like(xs,0,float) # degrees
    # endpoints are edge cases handled separately
    direction[0] = 180.0/np.pi*(np.arctan2(xs[0]-xs[1], ys[0]-ys[1]))+180.0
    direction[-1] = 180.0/np.pi*(np.arctan2(xs[-2]-xs[-1], ys[-2]-ys[-1]))+180.0
    
    m = np.full_like(xs,True,bool)
    m[[0,-1]] = False
    
    m2 = np.roll(m,-1) # point behind
    m3 = np.roll(m,1)  # point ahead
    # central approximation
    direction[m] = 180.0/np.pi*np.arctan2(xs[m2]-xs[m3], ys[m2]-ys[m3])+180.0
    return direction

def _genPath(xs,ys,zs,name,isOffset=True,adjusted=False, parallel=True):
    global pool
    """Displays the radargram for a path.

    Parameters
    xs float array : Array of x coordinates of path.
    ys float array : Array of y coordinates of path.
    zs float array : Array of altitude/elevation above ground along path.
    name string : Name of file to save radargram as. Not saved if name is None.
    isOffset bool (optional) : Whether the given z coordinates are altitude or relative to the ground. Default is relative.
    adjusted bool (optional) : Shift the data for each point to align the response from the surface
        directly beneath the radar with the top of the plot.

    Returns
    returnData 2D float array : The radargram output.
    Returns -1 if unsuccessful."""
    direction = _makeDirections(xs,ys)
    n = len(xs)

    # env holds timestep/range to sample over for radargram and granularity of samples
    env = parameters.env
    reflectionModels = radar.models
    titles = radar.titles
    
    returnData = np.full((len(reflectionModels),n,env.getSteps()),0,float) # 3D - many subplots
    plt.rcParams['axes.formatter.limits'] = [-4,4] # use standard form
    plt.figure(figsize=parameters.figsize)

    p = mp.Pool(mp.cpu_count())
    # arguments needed by processors as global state not shared
    data = [(x,y,z,i,isOffset,angle,reflectionModels, env) for x,y,i,z,angle in
                            zip(xs,ys,np.arange(n),zs,direction)]
    try: # calculate output across multiple processors
        if parallel:
            for i, ars in progress(p.imap_unordered(_worker,data),n):
                returnData[:,i] = ars
        else: # non-parallel option
            for j in range(n):
                i,ars = _worker(data[j])
                returnData[:,i] = ars
    except IOError as e: # most likely couldn't find maps.hdf5 in current directory
        p.close()
        print "\nError reading 'maps.hdf5' :\n" + e.message
        return -1
    p.close()

    returnData = returnData[:,np.any(returnData[0] != 0,1)] # remove invalid rows
    n = returnData.shape[1]
    
    if adjusted: # align first response of each point at top of plot
        returnData = align.minAlign(returnData, env.getDx())

    # allows pyplot to show subplots in square grid
    cells = int(np.ceil(np.sqrt(len(reflectionModels))))*110
    ys = np.linspace(0, env.getMaxTime(), env.getSteps())
    for j in range(len(reflectionModels)):
        plt.subplot(cells+j+1) # subplot index
        plt.ylim(env.getMaxTime(),0) # t=0 at top of plot
        draw = np.swapaxes(returnData[j],0,1)
        # colors adjusted so that mean value is 50% grey
        plt.contourf(np.arange(n), ys, draw, 100,norm=radar.MidNorm(np.mean(draw)), cmap="Greys")
        plt.title(titles[j])
        plt.colorbar()
    if name is not None:
        plt.savefig(name)
    plt.show()
    return returnData

def _worker(args): 
    # can raise IOError when reading maps.hdf5 in generateMaps
    # or in pointData for pointX.hdf5
    pointx,pointy,pointz,i,isOffset,angle,reflectionModels, env = args
    
    radar.env = env
    
    _,_,dist,incidence,theta,phi,elevation = pointData.generateMaps(pointx,pointy,pointz,
                                                                  isOffset,angle)
    
    ars = np.full((len(reflectionModels), env.getSteps()),0,float)
    for j in range(len(reflectionModels)):
        ars[j] = radar.processSlice(dist,incidence,theta,phi,reflectionModels[j])

    return i, ars
