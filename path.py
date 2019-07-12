import stateless
import progressbar
from numpy import genfromtxt
import xml.etree.cElementTree as ET
import utm
import multiprocessing as mp
#https://docs.python.org/3/library/xml.etree.elementtree.html#xml.etree.ElementTree.Element

def workerCall(args):
    x,y,name,z,isOffset,angle = args
    stateless.generateMaps(x,y,name,z,isOffset,angle)
    
def genPath(xs,ys,zs,name,isOffset=True):
    """Generates path data for the specified points, including antenna orientation data.
    Parameters:
    xs float array : array of x coordinates of path.
    ys float array : array of y coordinates of path.
    zs float array : array of altitude/elevation above ground along path.
    isOffset boolean (optional) : whether then given z coordinates are altitude or relative to the ground. Default is relative."""
    direction = np.full_like(xs,0,float) # degrees
    direction[0] = 180.0/np.pi*(np.arctan2(xs[0]-xs[1], ys[0]-ys[1]))+180.0
    direction[-1] = 180.0/np.pi*(np.arctan2(xs[-2]-xs[-1], ys[-2]-ys[-1]))+180.0
    m = np.full_like(xs,1,int)
    m[[0,-1]] = 0
    m2 = np.roll(m,-1)
    m3 = np.roll(m,1)
    direction[m] = 180.0/np.pi*(xs[m2]-xs[m3], ys[m2]-ys[m3]))+180.0

    pool = mp.Pool(mp.cpu_count())
    data = [(x,y,name+"/point"+str(i),z,isOffset,angle) for x,y,i,z,angle in
                            zip(xs,ys,np.arange(_steps),zs,direction)]
    fail = False
    for r in progressbar(pool.imap_unordered(workerCall,data),max_value=_steps):
        if r == -1:
            fail = True
    pool.close()
    if fail:
        print "Failed, setup didn't run correctly. Likely issue with state copying between processes"
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
    with open(filename,"r") as f:
        isOffset = eval(f.readline())
    data = genfromtxt(filename,delimiter=",",skip_header=1).swapaxes(0,1)
    return isOffset
    genPath(data[0],data[1],data[2],outName,isOffset)

def loadGpx(filename,outName=None):
    root = ET.parse(filename).getroot()
    if outName is None:
        n = root.find("name")
        if n is None:
            outName = filename[:-4] # removes .xml ending
        else:
            outName = n.text
    lats,lons,zs = [],[],[]
    track = root.find("trk")
    for seg in track.findall("trkseg"):
	for pt in seg.findall("trkpt"):
	    lats.append(pt.get("lat"))
	    lons.append(pt.get("lon"))
	    zs.append(pt.find("ele").text)
    lats = np.array(lats,float)
    lons = np.array(lons,float)
    zs = np.array(zs,float)
    # Convert Lon/lat to x,y for coordinate system
    xs,ys,zoneNum,zoneLet = utm.from_latlon(lats,lons) # says zone 45R?
    # x and y coordinates still appear correct though
    gePath(xs,ys,zs,outName,False)
    
