import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py
from declutter import radar
from declutter import modelling
from modelling import analysisFilters
import os
'''Methods for identifying the cause of peaks in a radargram. Also a method
for finding the incidence angle to the glacier wall at each point.'''


def flyBy(dirname, above=False, stepsize=3):
    '''Attempts to identify the glacier wall and returns the incidence angle
    to the wall for each point. The points are grouped into nearby sections
    and the start indices of each group is also returned.

    Parameters
    ----------
    dirname - string : Name of the directory containing pointX.hdf5 files for
        the path.
    above - bool (optional) : By default a surface must be above the radar to
        be considered part of the wall. Setting this means it must be at
        least 20m above.'''
    global fig, ax
    fig = plt.figure(figsize=(10, 7))
    ax = fig.gca(projection='3d')
    plt.subplots_adjust(0, 0, 1, 1)
    xs, ys, zs = [], [], []  # radar points
    surfxs, surfys, surfzs, aspects = [], [], [], []  # wall points
    pointCols = []
    grid, left, low, cellSize = loadMap()
    with h5py.File("maps.hdf5", "r") as f:
        slope = f["slope"][()]
        aspect = f["aspect"][()]

    height = grid.shape[0]
    files = os.listdir(dirname)
    files.sort(key=lambda x: int(x[5:-5]))
    cols = plt.cm.jet

    i = 0.0
    for f in files:
        i += 1.0
        visible, visCorner, distance, angle, _, _, pointx, pointy, pointz, \
            antDir = loadPoint(dirname+"/"+f)
        if len(distance) == 0:
            continue  # invalid point
        xs.append(pointx)
        ys.append(pointy)
        zs.append(pointz)

        # maps cropped to area around radar point
        s = matchMap(slope, left, low, visCorner[0], visCorner[1],
                     visible.shape[1], visible.shape[0], cellSize)
        a = matchMap(aspect, left, low, visCorner[0], visCorner[1],
                     visible.shape[1], visible.shape[0], cellSize)
        g = matchMap(grid, left, low, visCorner[0], visCorner[1],
                     visible.shape[1], visible.shape[0], cellSize)

        dx = int((visCorner[0]-left)/cellSize)
        dy = int((visCorner[1] - low)/cellSize)

        # enforce how far above the point the wall must be
        if above:
            offset = 20
        else:
            offset = 0

        # Use several filters to determine where the glacier wall is.
        # Note - first filter will create issues if radar goes within
        # 400m of wall.
        m = analysisFilters.compose(
            [analysisFilters.minDist(distance, 400),  # > 400m away
             # 5 degrees to horizon
             analysisFilters.horizon(g, visible, distance, pointz, 5),
             analysisFilters.steepness(s, 45),  # > 45 degree gradient
             # above point
             analysisFilters.height(g, pointz, offset)], visible)
        idx = np.argmin(distance + (~m)*5000)  # nearest valid point

        # y then x indices into 2D grid
        inds = [a[idx] for a in np.where(visible)]

        # x,y coordinates of point on wall
        inds = [inds[0]*cellSize+visCorner[1], inds[1]*cellSize+visCorner[0]]

        z = grid[int((inds[0]-low)/cellSize), int((inds[1]-left)/cellSize)]
        aspects.append(aspect[int((inds[0]-low)/cellSize),
                              int((inds[1]-left)/cellSize)])
        c = cols(i/len(files))

        surfxs.append(inds[1])
        surfys.append(inds[0])
        surfzs.append(z)
        pointCols.append(c)

    surfxs, surfys, surfzs = (np.array(surfxs), np.array(surfys),
                              np.array(surfzs))
    xs, ys, zs = np.array(xs), np.array(ys), np.array(zs)
    pointCols, aspects = np.array(pointCols), np.array(aspects)
    groups = np.arange(len(pointCols))
    # ax.scatter(surfxs,surfys,surfzs,color="r",linewidths=2)

    # Distance to plot surface for around path
    extend = 800

    # plot surface
    height, width = grid.shape
    xcoords = (np.amin(xs)-left)/cellSize, (np.amax(xs)-left)/cellSize
    xcoords = (max(0, int(xcoords[0]-extend/cellSize)),
               min(width, int(xcoords[1]+extend/cellSize)))
    ycoords = (np.amin(ys)-low)/cellSize, (np.amax(ys)-low)/cellSize
    ycoords = (max(0, int(ycoords[0]-extend/cellSize)),
               min(height, int(ycoords[1]+extend/cellSize)))
    grid = grid[ycoords[0]:ycoords[1], xcoords[0]:xcoords[1]]
    Y, X = np.indices(grid.shape)
    Y += ycoords[0]
    X += xcoords[0]

    size = grid.size
    if size > 200**2:
        factor = int(np.ceil(size**0.5/200.0))
        grid = grid[::factor, ::factor]
        X, Y = X[::factor, ::factor], Y[::factor, ::factor]
    X = left + X*cellSize
    Y = low + Y*cellSize

    my_col = plt.cm.terrain((grid-np.nanmin(grid)) /
                            (np.nanmax(grid)-np.nanmin(grid)))
    ax.plot_surface(X, Y, grid, facecolors=my_col, linewidth=0,
                    antialiased=False, zorder=3)
    # end

    # determine groups
    distances = []  # distance between surface points
    for i in range(len(xs)-1):
        distances.append(((surfxs[i]-surfxs[i+1])**2 +
                          (surfys[i]-surfys[i+1])**2 +
                          (surfzs[i]-surfzs[i+1])**2)**0.5)
    distances = np.array(distances)

    lengths = np.hypot(np.hypot(xs-surfxs, ys-surfys), zs-surfzs)

    for i in range(len(xs)-1):
        if distances[i] < 800:  # group surfaces within 800m of each other
            groups[i+1] = groups[i]
    pointCols = pointCols[groups]

    # update surface for each radar point to single point per group
    for g in np.unique(groups):
        i = np.argmin(lengths[groups == g])
        j = np.where(groups == g)[0][0]
        k = np.where(groups == g)[0][-1]+1
        surfxs[j:k] = surfxs[j+i]
        surfys[j:k] = surfys[j+i]
        surfzs[j:k] = surfzs[j+i]

        # approximation of surface normal as average vector of each normal
        aspects[j:k] = (np.arctan2(-np.sum(np.sin(aspects[j:k]*np.pi/180)),
                                   -np.sum(np.cos(aspects[j:k]*np.pi/180))) *
                        180.0/np.pi + 180.0)
    # end

    # Plot path
    ax.scatter(xs[::stepsize], ys[::stepsize], zs[::stepsize], color="k",
               linewidths=2, zorder=10)
    # plot lines to walls
    for i in range(0, len(xs), stepsize):
        ax.plot([xs[i], surfxs[i]], [ys[i], surfys[i]], [zs[i], surfzs[i]],
                color=pointCols[i], zorder=5)

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # finding incidence
    aspects *= np.pi / 180.0  # convert to radians
    incidence = np.array([(xs[i]-surfxs[i])*np.sin(aspects[i]) +
                          (ys[i]-surfys[i])*np.cos(aspects[i])
                          for i in range(len(xs))])
    incidence /= np.hypot([(xs[i]-surfxs[i]) for i in range(len(xs))],
                          [(ys[i]-surfys[i]) for i in range(len(xs))])
    incidence = np.arccos(incidence)*180.0/np.pi

    plt.show()
    plt.close()
    return incidence, np.unique(groups)


def matchMap(grid, left, low, cropLeft, cropLow, w, h, cellSize):
    '''Crops the whole array to just the extent of the visible array
    for a single point.'''
    l = int((cropLeft-left)/cellSize)
    r = int(l+w)
    lower = int((cropLow-low)/cellSize)
    upper = int(lower+h)
    return grid[max(0, lower):upper, max(0, l):r]


def markSurfaces(filename, start, end, alpha=0.2):
    '''Indicate which surfaces cause a response within the given interval
    for the pointX.hdf file given.
    Note: Time is taken as distance/1.5e8, ignoring any shift caused by
    the wave used for convolution (see radar.showWave()). Therefore intervals
    less than about 5e-7 may give unexpected results.

    Parameters
    ----------
    filename - string : The pointX.hdf5 file to plot the response of.
    start - float : The start of the time interval to consider.
    end - float : The end of the time interval to consider.
    alpha - float (optional) : Transparency of the surface plot from 0 to 1.'''
    global fig, ax
    load = loadPoint(filename)
    if load == -1:  # Can't use file, see loadPoint()
        return -1
    visible, visCorner, distance, angle, theta, phi, pointx, pointy, pointz, \
        antDir = load

    grid, left, low, cellSize = loadMap()

    height = grid.shape[0]

    l = int((visCorner[0]-left)/cellSize)
    r = int(l+visible.shape[1])
    lower = int((visCorner[1]-low)/cellSize)
    upper = int(lower + visible.shape[0])

    grid = grid[lower:upper, l:r]  # crop to extent of visible array

    m = (radar._time(distance) > start) & (radar._time(distance) < end)
    # default intensity model used
    intensityModel = modelling.defaults.default.getIntensity()
    # grid to store response intensities in
    intensity = np.full(grid.shape, 0.0, float)
    which = tuple([a[m] for a in np.where(visible)])
    intensity[which] = intensityModel(angle[m])

    height = visible.shape[0]
    ys, xs = np.indices(visible.shape)

    xs = visCorner[0] + xs*cellSize
    ys = visCorner[1] + ys*cellSize
    m = intensity != 0
    print "{0} points in this interval".format(np.sum(m))
    print "Total intensity: {0}".format(np.sum(intensity))

    fig = plt.figure(figsize=(10, 7))
    ax = fig.gca(projection='3d')
    plt.subplots_adjust(0, 0, 1, 1)

    drawMesh([pointx], [pointy], end*1.8e8, True, alpha)
    drawPoint(pointx, pointy, pointz, antDir)
    p = ax.scatter(xs[m], ys[m], grid[m], c=intensity[m], cmap="jet")
    plt.colorbar(p, shrink=0.7)
    plt.show()
    plt.close()
    return


def drawPoint(x, y, z, antDir=None, length=150):
    ax.scatter([x], [y], [z], color="k", linewidths=2, zorder=20)
    if antDir is not None:
        ax.plot([x, x+antDir[0]*length], [y, y+antDir[1]*length], [z, z],
                linewidth=2, color="k")


def drawMesh(x, y, distance, back=False, alpha=0.2):
    '''Plots a wireframe mesh of the surrounding surface.

    Parameters
    ----------
    x,y - floats : Coordinates of point.
    distance - float : Range around the point to draw up to.
    back - bool (optional) : Whether this should create a new
        figure or use the existing ax and fig. Also whether the
        plot should be displayed immediately or not.
    alpha - float (optional) : Transparency of the wireframe plot.'''
    global ax, fig
    grid, left, low, cellSize = loadMap()
    height, width = grid.shape

    minCoords = (np.amin(x)-left)/cellSize, (np.amin(y)-low)/cellSize
    maxCoords = (np.amax(x)-left)/cellSize, (np.amax(y)-low)/cellSize
    dist = distance/cellSize
    grid = grid[int(minCoords[1]-dist):int(maxCoords[1]+dist),
                int(minCoords[0]-dist):int(maxCoords[0]+dist)]

    ys, xs = np.indices(grid.shape)

    xs = left + cellSize*(xs+int(minCoords[0]-dist))
    ys = low + cellSize*(ys+int(minCoords[1]-dist))

    if not back:  # create new figure
        fig = plt.figure(figsize=(10, 7))
        ax = fig.gca(projection='3d')
        plt.subplots_adjust(0, 0, 1, 1)
    ax.plot_wireframe(xs, ys, grid, rcount=30, ccount=30, zorder=1,
                      alpha=alpha)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    if not back:
        plt.show()


def loadMap():
    with h5py.File("maps.hdf5", "r") as f:
        grid = f["heightmap"][()]
        left, low, cellSize = f["meta"][:3]
    return grid, left, low, cellSize


def loadPoint(filename):
    antDir = None
    with h5py.File(filename, "r") as f:
        if "visible" not in f:
            print "Only possible for points which stored visible array."
            return -1
        visible = f["visible"][()]
        if "corner" not in f:
            print "Please recreate point files, format has changed."
            return -1
        visCorner = f["corner"][()]
        distance = f["distance"][()]
        angle = f["incidence"][()]
        theta, phi = None, None
        if "antennaTheta" in f:
            theta = f["antennaTheta"][()]
            phi = f["antennaPhi"][()]
        pointx, pointy, pointz, antDir = f["meta"][:4]
    if antDir is not None:
        antDir = (np.sin(antDir*np.pi/180.0), np.cos(antDir*np.pi/180.0))
    return (visible, visCorner, distance, angle, theta, phi, pointx,
            pointy, pointz, antDir)
