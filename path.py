import stateless
from progress import progress
import numpy as np
import xml.etree.cElementTree as ET
import utm
import multiprocessing as mp
#https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element


def workerCall(args):
    return stateless.generateMaps(*args)
    
def genPath(xs,ys,zs,name,isOffset=True):
    global pool
    """Generates path data for the specified points, including antenna orientation data.
    Parameters:
    xs float array : array of x coordinates of path.
    ys float array : array of y coordinates of path.
    zs float array : array of altitude/elevation above ground along path.
    isOffset boolean (optional) : whether then given z coordinates are altitude or relative to the ground. Default is relative."""
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
    data = [(x,y,name+"/point"+str(i),z,isOffset,angle) for x,y,i,z,angle in
                            zip(xs,ys,np.arange(steps),zs,direction)]
    fail = False
    for r in progress(pool.imap_unordered(workerCall,data),steps):
        if r == -1:
            fail = True
            break
    pool.close()
    if fail:
        print "Failed to generate one or more points using stateless.generateMaps"
        return -1
    return 0     
    
# only display from crop[0] to end-crop[1] e.g. [250,220]
def loadGpx(filename,crop=[0,0],outName=None):
    try:
        root = ET.parse(filename).getroot()
    except IOError:
        print "Error in path.py, could not load file: "+filename
        print "Check filename correct and data exists in gpx format."
        return -1
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
    # Convert Lon/lat to x,y for coordinate system
    # Assuming z should be unchanged, check by 'showOnSurface()'
    xs,ys,zoneNum,zoneLet = utm.from_latlon(lats,lons)
    return genPath(xs,ys,zs,outName,False)

def showPath(filename,crop=[0,0]):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    root = ET.parse(filename).getroot()
    prefix = ""
    if root.get("targetNamespace") is not None:
        prefix = "{"+root.get("targetNamespace")+"}"
    lats,lons,zs = [],[],[]
    for pt in root.iter(prefix+"trkpt"): # recursively searches for points in document order
        lats.append(pt.get("lat"))
        lons.append(pt.get("lon"))
        zs.append(pt.find(prefix+"ele").text)
    n = len(lats)
    lats = np.array(lats[crop[0]:n-crop[1]],float)
    lons = np.array(lons[crop[0]:n-crop[1]],float)
    zs = np.array(zs[crop[0]:n-crop[1]],float)
    # Convert Lon/lat to x,y for coordinate system
    xs,ys,zoneNum,zoneLet = utm.from_latlon(lats,lons)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs,ys,zs)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.show()
    
def showOnSurface(filename,crop=[0,0],extend=10):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    root = ET.parse(filename).getroot()
    prefix = ""
    if root.get("targetNamespace") is not None:
        prefix = "{"+root.get("targetNamespace")+"}"
    lats,lons,zs = [],[],[]
    for pt in root.iter(prefix+"trkpt"): # recursively searches for points in document order
        lats.append(pt.get("lat"))
        lons.append(pt.get("lon"))
        zs.append(pt.find(prefix+"ele").text)
    n = len(lats)
    lats = np.array(lats[crop[0]:n-crop[1]],float)
    lons = np.array(lons[crop[0]:n-crop[1]],float)
    zs = np.array(zs[crop[0]:n-crop[1]],float)
    xs,ys,zoneNum,zoneLet = utm.from_latlon(lats,lons)

    heightmap = np.load("maps.npz")["heightmap"]
    height,width = heightmap.shape
    with open("info.txt","r") as f:
        line = f.read().split(",")
    left = float(line[0])
    low = float(line[1])
    cellSize = float(line[2])
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

    # replaced noData value with nan so no longer needed
    #heightmap[heightmap < -50000] = np.nan
    
    my_col = cm.jet((heightmap-np.nanmin(heightmap))/(np.nanmax(heightmap)-np.nanmin(heightmap)))
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X,Y,heightmap,facecolors=my_col,linewidth=0,antialiased=False)
    ax.plot(xs,ys,zs)
    plt.show()

# Assumes NaN is undefined value, otherwise need to change equality test
# (Can't use == for NaNs, need np.isnan(...) instead)
def checkValid(filename,crop = [0,0]):
    """Indicate if any of the map is undefined within 3km (current fixed range) of a point on the path.
    Also highlights if the map is undefined directly beneath any points on the path."""
    root = ET.parse(filename).getroot()
    prefix = ""
    if root.get("targetNamespace") is not None:
        prefix = "{"+root.get("targetNamespace")+"}"
    lats,lons,zs = [],[],[]
    for pt in root.iter(prefix+"trkpt"): # recursively searches for points in document order
        lats.append(pt.get("lat"))
        lons.append(pt.get("lon"))
        zs.append(pt.find(prefix+"ele").text)
    n = len(lats)
    lats = np.array(lats[crop[0]:n-crop[1]],float)
    lons = np.array(lons[crop[0]:n-crop[1]],float)
    zs = np.array(zs[crop[0]:n-crop[1]],float)
    xs,ys,zoneNum,zoneLet = utm.from_latlon(lats,lons)

    heightmap = np.load("maps.npz")["heightmap"]
    height,width = heightmap.shape
    with open("info.txt","r") as f:
        line = f.read().split(",")
    left = float(line[0])
    low = float(line[1])
    cellSize = float(line[2])
    Y,X = np.indices(heightmap.shape)
    Y = Y[::-1]
    X = left + X*cellSize
    Y = low + Y*cellSize
    right = left + width * cellSize
    high = low + height * cellSize
    xBounds = [left+3000,right-3000]
    yBounds = [low+3000,high-3000]
    
    notFullRange = 0
    undefinedGround = 0
    for x,y in zip(xs,ys):
        # points within 3km
        m = ((X-x)**2+(Y-y)**2)<3000**2
        if np.any(np.isnan(heightmap[m])) or x<xBounds[0] or x>xBounds[1] or y<yBounds[0] or y>yBounds[1]:
            notFullRange += 1
        # points used for interpolating ground height
        m = ((X-x)**2 <= 1) & ((Y-y)**2 <= 1)
        if np.any(np.isnan(heightmap[m])) or x < left or x > right or y < low or y > high:
            undefinedGround += 1
    print "{0} of {1} points have part of range undefined.".format(notFullRange,n)
    print "{0} of {1} points are above an undefined part of the map.".format(undefinedGround,n)    
    
