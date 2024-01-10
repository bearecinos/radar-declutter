'''This module combines the methods used in path.py, pointData.py and
radar.py to display a radargram for data without storing the files
for all points.
This is useful where a single output is needed and the data for each
point would be deleted immediately after.'''
from declutter import path
import numpy as np
import os
from declutter import radar
from declutter.progress import progress
import multiprocessing as mp
from declutter import pointData
from declutter import viewshed
import matplotlib.pyplot as plt
from declutter import align
from declutter.modelling import parameters
from declutter.modelling.defaults import default


def processData(filename, crop=[0, 0], outName=None, style=None,
                offset=0, adjusted=False, save=True, parallel=True):
    """Takes a gps path and displays a radargram for that path.

    Parameters
    ----------
    filename - string  : The name of the file to generate a radargram for.
    crop - [int, int] (optional) : crop = [A,B] ignores the first A and last B
        points of the path.
    outName - string (optional) : The name of the file to save the radargram
        in. By default, this is taken from 'filename' unless 'save' is set
        to False.
    style - string (optional) : The format of the input file, either 'gpx',
        'dst' or 'xyz'. By default, the loadData method determines the format
        from the file extension and assumes gpx if the extension is not
        recognised.
    offset - float (optional) : Height to correct for where gps data taken
        from helicopter, not on radar. Default 0.
    adjusted - bool (optional) : Shift the data for each point to align the
        response from the surface directly beneath the radar with the
        top of the plot.
    save - bool (optional) : Default True. If True, the radargram output
        is saved automatically.

    Returns
    -------
    The radargram output if successful, otherwise -1.

    """
    print(parameters.env)
    try:
        xs, ys, zs = path.loadData(filename, crop, outName, style)
        xs -= offset
    except IOError:
        print("Could not load data from file : "+filename)
        if style is not None:
            print("Is "+style+" the correct format?")
        return -1
    if len(xs) < 2:  # not enough points for a path
        return -1
    if outName is None and save:
        outName = filename[:-4]+".png"
    return _genPath(xs, ys, zs, outName, adjusted, parallel)


def _genPath(xs, ys, zs, name, adjusted=False, parallel=True):
    global pool
    """
    Displays the radargram for a path.

    Parameters
    ----------
    xs - float array : Array of x coordinates of path.
    ys - float array : Array of y coordinates of path.
    zs - float array : Array of altitude/elevation above ground along path.
    name - string : Name of file to save radargram as. The plot is not
        saved automatically if None.
    adjusted - bool (optional) : Shift the data for each point to align the
        response from the surface directly beneath the radar with the top of
        the plot.

    Returns
    -------
    returnData - 2D float array : The radargram output.
    Returns -1 if unsuccessful.
    """

    direction = path._makeDirections(xs, ys)
    n = len(xs)

    # env holds timestep/range to sample over
    env = parameters.env

    intensityModel = default.getIntensity()
    wave = default.getWave()
    directional = default.getDirectivity()

    returnData = np.full((n, env.getSteps()), 0, float)
    plt.rcParams['axes.formatter.limits'] = [-4, 4]  # use standard form
    plt.figure(figsize=env.figsize)

    p = mp.Pool(mp.cpu_count())
    # arguments needed by processors as global state not shared
    data = [(x, y, z, i, angle, wave, intensityModel, directional, env) for
            x, y, i, z, angle in zip(xs, ys, np.arange(n), zs, direction)]
    try:  # calculate output across multiple processors
        if parallel:
            for i, ar in progress(p.imap_unordered(_worker, data), n):
                returnData[i] = ar
        else:  # non-parallel option
            for j in range(n):
                i, ar = _worker(data[j])
                returnData[i] = ar
    except IOError as e:  # likely couldn't find maps.hdf5 in current directory
        p.close()
        print("\nError reading 'maps.hdf5' :\n" + e.message)
        return -1
    p.close()

    n = returnData.shape[0]

    if adjusted:  # align first response of each point at top of plot
        returnData = align.minAlign(returnData, env.getDx())

    ys = np.linspace(0, env.getMaxTime(), env.getSteps())
    plt.ylim(env.getMaxTime(), 0)  # t=0 at top of plot
    draw = np.swapaxes(returnData, 0, 1)
    # colors adjusted so that mean value is 50% grey
    plt.contourf(np.arange(n), ys, draw, 100,
                 norm=radar.MidNorm(np.mean(draw)), cmap="Greys")

    plt.colorbar()
    if name is not None:
        plt.savefig(name)
    plt.show()
    return returnData


def _worker(args):
    # can raise IOError when reading maps.hdf5 in generateMaps
    # or in pointData for pointX.hdf5
    pointx, pointy, pointz, i, angle, wave, intensityModel, directional, env \
            = args

    # New process so env may have been reset, updates to passed parameters.
    # This means pointData and radar both get the correct range/steps etc.
    parameters.setEnv(env)

    _, _, dist, incidence, theta, phi, _ = pointData.generateMaps(pointx,
                                                                  pointy,
                                                                  pointz,
                                                                  angle)
    ar = np.full((env.getSteps()), 0, float)
    ar = radar.processSlice(dist, incidence, theta, phi, intensityModel, wave,
                            directional=directional)

    return i, ar
