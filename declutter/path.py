import pointData
from progress import progress
import numpy as np
import xml.etree.cElementTree as ET
import utm
import pyproj
import multiprocessing as mp
import h5py
import os
#https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element


def workerCall(args):
    # args = (x,y,z,isOffset,angle,pathName)
    vis,visCorner,dist,incidence,theta,phi,elevation = pointData.generateMaps(*args[:-2])

    if args[5]:
        return pointData.store(args[6],dist,incidence,args[0],args[1],
                           elevation, vis,visCorner, args[4], theta, phi)
    else:
        return pointData.store(args[6], dist,incidence,args[0],args[1],
                           elevation, None,None, args[4], theta, phi)
    
def _genPath(xs,ys,zs,name,isOffset=True,save_visible=True,parallel=True):
    global pool
    """Generates path data for the specified points, including antenna orientation data.
    Parameters:
    xs - float array : array of x coordinates of path.
    ys - float array : array of y coordinates of path.
    zs - float array : array of altitude/elevation above ground along path.
    isOffset - bool (optional) : whether then given z coordinates are altitude or relative to the ground. Default is relative."""
    os.makedirs(name)
    direction = np.full_like(xs,0,float) # degrees
    direction[0] = 180.0/np.pi*(np.arctan2(xs[0]-xs[1], ys[0]-ys[1]))+180.0
    direction[-1] = 180.0/np.pi*(np.arctan2(xs[-2]-xs[-1], ys[-2]-ys[-1]))+180.0
    m = np.full_like(xs,True,bool)
    m[[0,-1]] = False
    m2 = np.roll(m,-1)
    m3 = np.roll(m,1)
    direction[m] = 180.0/np.pi*np.arctan2(xs[m2]-xs[m3], ys[m2]-ys[m3])+180.0
    steps = len(xs)
    pool = mp.Pool(mp.cpu_count())
    data = [(x,y,z,isOffset,angle,save_visible,name+"/point"+str(i)) for
             x,y,z,angle,i in zip(xs,ys,zs,direction,np.arange(steps))]
    fail = False
    try:
        if parallel:
            for r in progress(pool.imap_unordered(workerCall,data),steps):
                pass
        else:
            for i in progress(range(steps)):
                result = workerCall(data[i])
    except IOError:
        pool.close()
        print "\nError reading 'maps.hdf5' :\n" + e.message
        return -1
    pool.close()
    return 0     

def processData(filename,crop=[0,0],outName=None,style=None, save_visible=True, parallel = True):
    """Generates path data for the points in the given file.

    Parameters
    filename - string : The file to generate data for.
    crop - [int,int] (optional) : crop=[A,B] ignores the first A and last B points of the input file.
    outName - string (optional) : The directory to save the generated data in. By default, this is taken
        from the input filename.
    style - string (optional) : One of 'gpx', 'dst' or 'xyz', indicating the format of 'filename'. By default,
        this is determined by the file extension and assumed to be 'gpx' if that is unclear.
    save_visible - bool (optional) : Whether or not the array of which points are visible and which aren't should
        be stored (significant proportion of the storage where only a few points are visible. Default is True.

    Returns
    0 if successful, otherwise -1.
    """
    try:
        xs, ys, zs = loadData(filename, crop, style)
    except IOError:
        print "Could not load data from file : "+filename
        if style is not None:
            print "Is "+style+" the correct format?"
        return -1
    if len(xs) == 0:
        return -1
    if outName is None:
        outName = filename[:-4]
    return _genPath(xs,ys,zs,outName,False,save_visible, parallel)   

def loadData(filename, crop = [0,0], style = None):
    """Takes a file of points and returns three arrays of x-coordinate, y-coordinate, and elevation.

    Parameters
    filename - string : The file to generate data for.
    crop - [int,int] (optional) : crop=[A,B] ignores the first A and last B points of the input file.
    style - string (optional) : One of 'gpx', 'dst' or 'xyz', indicating the format of 'filename'. By default,
        this is determined by the file extension and assumed to be 'gpx' if that is unclear.

    Returns
    xs, ys - The x and y coordinates of the input points. If the input file had longitude and latitude coordinates,
        these have been mapped to a UTM zone or UPS for coordinates near the poles.
    zs - Elevation values taken directly from the file.
    """
    if style is None:
        if filename[-3:] == "dst":
            style = "dst"
        elif filename[-3:] == "xyz":
            style = "xyz"
        else:
            style = "gpx"
    if style == "gpx":
        lons,lats,zs,outName = _loadGpx(filename,crop)
        xs,ys = gpsToXY(lons,lats)
    elif style == "dst":
        lons,lats,zs = _loadDst(filename,crop)
        xs,ys = gpsToXY(lons,lats)
    elif style == "xyz":
        try:
            data = np.loadtxt(filename)
        except ValueError as e:
            print "Could not read as xyz file: " + e.message
            raise IOError("Could not read file: "+filename)
        n = len(data)
        xs,ys,zs = data[crop[0]:n-crop[1]].swapaxes(0,1)
    else:
        print "Format not recognised. should be 'gpx', 'dst' or 'xyz'"
        xs, ys, zs = np.array([[],[],[]])
    return xs, ys, zs
    
    
def _loadGpx(filename,crop=[0,0],outName=None):
    root = ET.parse(filename).getroot()
    prefix = ""
    if root.get("targetNamespace") is not None:
        prefix = "{"+root.get("targetNamespace")+"}"
    
    if outName is None:
        try:
            n = root.iter(prefix+"name").next()
            outName = n.text
        except StopIteration:
            outName = filename
        if outName[-4:] in [".gpx", ".xml"]:
            outName = outName[:-4] # removes file extension from name
    lats,lons,zs = [],[],[]
    for pt in root.iter(prefix+"trkpt"): # recursively searches for points in document order
        lats.append(pt.get("lat"))
        lons.append(pt.get("lon"))
        zs.append(pt.find(prefix+"ele").text)
    n = len(lats)
    lats = np.array(lats[crop[0]:n-crop[1]],float)
    lons = np.array(lons[crop[0]:n-crop[1]],float)
    zs = np.array(zs[crop[0]:n-crop[1]],float)
    return lons,lats,zs,outName

def _loadDst(filename,crop=[0,0],noData=0.0):
    try:
        data = np.loadtxt(filename)
    except ValueError as e:
        print "Could not read as dst file: " + e.message
        raise IOError("Could not read file: "+filename)
    n = len(data)
    data = data[crop[0]:n-crop[1]].swapaxes(0,1)
    zs = data[2]
    m = zs != noData
    zs = zs[m]
    lats = data[0][m]
    lons = data[1][m]
    return lons,lats,zs

_northProj = pyproj.Proj("+proj=ups")
_southProj = pyproj.Proj("+proj=ups +south")
_gpsProj = pyproj.Proj("+init=EPSG:4326") # wgs 84

def gpsToXY(lons,lats):
    """Converts an array of longitude and latitude coordinates to x,y coordinates
    for the relevant zone, either a UTM zone or UPS."""
    if -79.5 <= lats[0] and lats[0] <= 83.5:
        xs,ys,zoneNum,zoneLet = utm.from_latlon(lats,lons)
    elif -79.5 > lats[0]:
        print "Need to use UPS coordinates."
        print "Not tested if elevation values still accurate on UPS map."
        xs,ys = pyproj.transform(_gpsProj,_southProj,lons,lats)
    else:
        print "Need to use UPS coordinates."
        print "Not tested if elevation values still accurate on UPS map."
        xs,ys = pyproj.transform(_gpsProj,_northProj,lons,lats)
    return xs,ys

def showPath(filename,crop=[0,0],style=None):
    """Plots a given path in 3D.

    Parameters
    filename - string : The file to generate data for.
    crop - [int,int] (optional) : crop=[A,B] ignores the first A and last B points of the input file.
    style - string (optional) : One of 'gpx', 'dst' or 'xyz', indicating the format of 'filename'. By default,
        this is determined by the file extension and assumed to be 'gpx' if that is unclear.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    try:
        xs, ys, zs = loadData(filename, crop, style)
    except IOError:
        print "Could not load data from file : "+filename
        return -1
    if len(xs) == 0:
        return -1
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs,ys,zs)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.show()
    return 0

def showAboveGround(filename,crop=[0,0],style=None):
    """Uses 'maps.hdf5' to show the radar elevation relative to the ground directly beneath the radar along the path.

    Parameters
    filename - string : The file to generate data for.
    crop - [int,int] (optional) : crop=[A,B] ignores the first A and last B points of the input file.
    style - string (optional) : One of 'gpx', 'dst' or 'xyz', indicating the format of 'filename'. By default,
        this is determined by the file extension and assumed to be 'gpx' if that is unclear.
    """
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import viewshed
    
    try:
        xs, ys, zs = loadData(filename, crop, style)
    except IOError:
        print "Could not load data from file : "+filename
        return -1
    if len(xs) == 0:
        return -1

    with h5py.File("maps.hdf5","r") as f:
        grid = f["heightmap"][()]
        left,low,cellSize = f["meta"][()]
    height,width = grid.shape

    xs = (xs-left)/cellSize
    ys = height-(ys-low)/cellSize
    groundHeights = viewshed.quadHeight(grid,xs,ys)
    plt.plot(zs,label="radar")
    plt.plot(groundHeights,label="ground")
    plt.legend()
    plt.show()
    return 0
    

def showOnSurface(filename,crop=[0,0],extend=10,style=None):
    """Plots the radar path in 3D, using 'maps.hdf5' to plot the surrounding terrain for reference.

    Parameters
    filename - string : The file to generate data for.
    crop - [int,int] (optional) : crop=[A,B] ignores the first A and last B points of the input file.
    style - string (optional) : One of 'gpx', 'dst' or 'xyz', indicating the format of 'filename'. By default,
        this is determined by the file extension and assumed to be 'gpx' if that is unclear."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    try:
        xs, ys, zs = loadData(filename, crop, style)
    except IOError:
        print "Could not load data from file : "+filename
        return -1
    if len(xs) == 0:
        return -1

    with h5py.File("maps.hdf5","r") as f:
        heightmap = f["heightmap"][()]
        left,low,cellSize = f["meta"][()]
    height,width = heightmap.shape
    
    xcoords = (np.amin(xs)-left)/cellSize, (np.amax(xs)-left)/cellSize
    xcoords = max(0,int(xcoords[0]-extend)), min(width,int(xcoords[1]+extend))
    ycoords = height - 1 - (np.amax(ys)-low)/cellSize, height - 1 - (np.amin(ys)-low)/cellSize
    ycoords = max(0,int(ycoords[0]-extend)), min(height,int(ycoords[1]+extend))
    Y,X = np.indices(heightmap.shape)
    Y = Y[::-1]
    heightmap = heightmap[ycoords[0]:ycoords[1],xcoords[0]:xcoords[1]]
    Y = Y[ycoords[0]:ycoords[1],xcoords[0]:xcoords[1]]
    X = X[ycoords[0]:ycoords[1],xcoords[0]:xcoords[1]]

    # scale resolution to display in reasonable time
    size = heightmap.size
    if size > 200**2:
        factor = int(np.ceil(size**0.5/200.0))
        heightmap = heightmap[::factor,::factor]
        X,Y = X[::factor,::factor], Y[::factor,::factor]
    X = left + X*cellSize
    Y = low + Y*cellSize
    
    my_col = cm.terrain((heightmap-np.nanmin(heightmap))/(np.nanmax(heightmap)-np.nanmin(heightmap)))
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X,Y,heightmap,facecolors=my_col,linewidth=0,antialiased=False)
    ax.plot(xs,ys,zs)
    plt.show()
    return 0

# Assumes NaN is undefined value, otherwise need to change equality test
# (Can't use == for NaNs, need np.isnan(...) instead)
def checkValid(filename,crop = [0,0],style="gpx"):
    """Indicate if any of the map is undefined within 3km (current fixed range) of a point on the path.
    Also highlights if the map is undefined directly beneath any points on the path.

    Parameters
    filename - string : The file to generate data for.
    crop - [int,int] (optional) : crop=[A,B] ignores the first A and last B points of the input file.
    style - string (optional) : One of 'gpx', 'dst' or 'xyz', indicating the format of 'filename'. By default,
        this is determined by the file extension and assumed to be 'gpx' if that is unclear."""
    try:
        xs, ys, zs = loadData(filename, crop, style)
    except IOError:
        print "Could not load data from file : "+filename
        return -1
    if len(xs) == 0:
        return -1
    
    with h5py.File("maps.hdf5","r") as f:
        heightmap = f["heightmap"][()]
        left,low,cellSize = f["meta"][()]
    height,width = heightmap.shape
    Y,X = np.indices(heightmap.shape)
    Y = Y[::-1]
    X = left + X*cellSize
    Y = low + Y*cellSize
    right = left + width * cellSize
    high = low + height * cellSize
    xBounds = [left+3000,right-3000]
    yBounds = [low+3000,high-3000]
    
    ar = np.full(len(xs),True)
    
    notFullRange = 0
    undefinedGround = 0
    for i in range(len(xs)):
        x,y  = xs[i], ys[i]
        # points within 3km
        m = ((X-x)**2+(Y-y)**2)<3000**2
        if np.any(np.isnan(heightmap[m])) or x<xBounds[0] or x>xBounds[1] or y<yBounds[0] or y>yBounds[1]:
            notFullRange += 1
        # points used for interpolating ground height
        m = ((X-x)**2 <= cellSize**2) & ((Y-y)**2 <= cellSize**2)
        if np.any(np.isnan(heightmap[m])) or x < left or x > right or y < low or y > high:
            undefinedGround += 1
            ar[i] = False
    print "{0} of {1} points have part of range undefined.".format(notFullRange,len(xs))
    print "{0} of {1} points are above an undefined part of the map.".format(undefinedGround,len(xs))
    return ar
    
