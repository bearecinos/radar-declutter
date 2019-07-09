import Generate
import progressbar
from numpy import genfromtxt

# alternative for z coordinate vs height above ground
def genPath(xs,ys,zs,name,isOffset=True):
    """Generates path data for the specified points, including antenna orientation data.
    Parameters:
    xs float array : array of x coordinates of path.
    ys float array : array of y coordinates of path.
    zs float array : array of altitude/elevation above ground along path.
    isOffset boolean (optional) : whether then given z coordinates are altitude or relative to the ground. Default is relative."""
    for i in progressbar.progressbar(range(len(xs))):
        if i == 0:
            direction = math.degrees(math.atan2(xs[i]-xs[i+1], ys[i]-ys[i+1]))+180
        elif i == len(xs)-1:
            direction = math.degrees(math.atan2(xs[i-1]-xs[i], ys[i-1]-ys[i]))+180
        else:
            direction = math.degrees(math.atan2(xs[i-1]-xs[i+1], ys[i-1]-ys[i+1]))+180
        Generate.setPoint(xs[i],ys[i],name+"\\point"+str(i),zs[i],isOffset,direction)
        Generate.generateMaps()
    Generate.finish()

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
