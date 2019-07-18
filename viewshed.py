import numpy as np
#import matplotlib.pyplot as plt
#import scipy.signal as signal
#from mpl_toolkits.mplot3d import Axes3D

# y,x = np.indices([50,50])
# fig = plt.figure()
# ax = fig.add_subplot(1, 3, 1, projection='3d')
# ax.plot_wireframe(x,y,res)

# To simplify initially, assume grid is n by n
# and have view point in grid coordinates e.g. 105.4,105.3 (decimal for accuracy)
# also assume outer bounds of grid never of interest so automatically not considered
# i.e. only do visible for grid[1:-1,1:-1]

### UPDATE
#
# Can no longer assume given point even within grid, have to handle edges
#
###

_NODATA = -100000

# should ignore points outside grid, and any point interpolating a NODATA value
# should also be ignored
def quadHeight(grid,x,y): # bi-quadratic interpolation of height
    res = np.full_like(x,_NODATA,float) # can mask points out of grid
    h,w = grid.shape
    m = (x >= 0) & (y >= 0) & (x < w-1) & (y < h-1) # valid points mask

    # later use nan in grid
    x,y = x[m], y[m]
        
    cornerx = x.astype(int)
    cornery = y.astype(int)
    # output NODATA if any cell to interpolate between is already NODATA
    # use np.nan in the future
    w = ((grid[cornery,cornerx] != _NODATA) & (grid[cornery+1,cornerx] != _NODATA) &
        (grid[cornery,cornerx+1] != _NODATA) & (grid[cornery+1,cornerx+1] != _NODATA))

    x,y = x[w], y[w]
    cornerx,cornery = cornerx[w], cornery[w]
    
    dx = x-cornerx
    tx = dx**2+(1-dx)**2
    dy = y-cornery
    ty = dy**2+(1-dy)**2
    a = (dx**2*grid[cornery,cornerx+1] + (1-dx)**2*grid[cornery,cornerx])/tx
    b = (dx**2*grid[cornery+1,cornerx+1] + (1-dx)**2*grid[cornery+1,cornerx])/tx
    res[tuple(mask[w] for mask in np.where(m))] = (dy**2*b + (1-dy)**2*a)/ty
    return res

def visible(grid,startx,starty,x,y,elevation=100,isOffset=True,stepSize=1.0):
    vis = np.full_like(x,1,"float64")
    if len(x) == 0:
        return vis # else get error when attempting to call np.amax on 0 length array
    endh = quadHeight(grid,x,y)
    if isOffset:
        starth = quadHeight(grid,np.array([startx]),np.array([starty]))+elevation
    else:
        starth = elevation
    dx = x-startx
    dy = y-starty
    norm = (dx**2+dy**2)**0.5

    # ignores divide by 0 cases in coming loop
    norm[norm==0] = 0.1
    
    dx /= norm
    dy /= norm
    dh = (endh-starth)/norm
    # more fine grain steps?

    thisx = startx
    thisy = starty
    thish = np.full_like(norm,starth)
    for i in range(1,int(np.amax(norm)/stepSize)):
        m = (norm >= i+1)&(vis==1) # equivalent to a non-inclusive upper bound of int(norm)
        # should this be the case? or should it just be >= i?
        thisx += dx*stepSize
        thisy += dy*stepSize
        thish += dh*stepSize
        h = quadHeight(grid,thisx[m],thisy[m]) 
        w = h >= thish[m]
        vis[tuple(a[w] for a in np.where(m))] = _NODATA
    return vis

# y passed in world coordinates so grid coordinates are height-1-y
def viewshed(grid,pointx,pointy,elevation=100,isOffset=True,maxRange=3000,gridsize=30.0,stepSize=None):
    gheight, gwidth = grid.shape
    
    pointy = gheight - pointy - 1
    
    view = np.full_like(grid,_NODATA,int)
    if pointx < 0 or pointy < 0 or pointx > gwidth - 1 or pointy > gheight - 1:
        return view # point off grid, can't tell what is visible or not
    maxSquaredCells = (maxRange/float(gridsize))**2

    ys, xs = np.indices(grid.shape,float)

    squareDist = (xs-pointx)**2+(ys-pointy)**2
    inRange = (squareDist < maxSquaredCells) & (grid != _NODATA) # ignore noData points

    if stepSize is None:
        stepSize = gridsize
    stepSize /= gridsize # step in grid coordinates
    
    view[inRange] = visible(grid,pointx,pointy,xs[inRange],ys[inRange],elevation,isOffset,stepSize)
    
    return view
  

