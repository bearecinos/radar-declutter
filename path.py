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
    
# format: 1st line is "True"/"False" for isOffset
# each other line is x,y,z
def loadFromFile(filename,outName=None):
    """Reads a file of path coordinates and generates the required data for those points.
    Parameters:
    filename string : path of file to read coordinates from.
    outName string (optional) : name of folder to store generated data in. Defaults to same path as filename if not given.
    File format:
    First line - 'True' or 'False' to indicate whether z-coordinates are altitude or offset from surface.
    Each next line - 'x,y,z' comma separated coordinates of next point along path."""
    if outName is None:
        outName = filename
    try:
        with open(filename,"r") as f:
            isOffset = eval(f.readline())
            data = np.genfromtxt(filename,delimiter=",",skip_header=1).swapaxes(0,1)
    except IOError:
        print "Error in path.py, could not load file: "+filename
        print "Check filename correct and file in specified format."
        return -1
    return genPath(data[0],data[1],data[2],outName,isOffset)

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
    # How to convert z?
    xs,ys,zoneNum,zoneLet = utm.from_latlon(lats,lons)
    #print "Warning : Currently treating elevation as relative to surface."
    # CHANGE 'True' as coordinates not relative - only workaround to make rest work
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
    # How to convert z?
    xs,ys,zoneNum,zoneLet = utm.from_latlon(lats,lons)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs,ys,zs)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.show()
    
