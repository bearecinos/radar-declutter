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

def quadHeight(grid,x,y): # bi-quadratic interpolation of height
    cornerx = x.astype(int)
    cornery = y.astype(int)
    dx = x-cornerx
    tx = dx**2+(1-dx)**2
    dy = y-cornery
    ty = dy**2+(1-dy)**2
    a = (dx**2*grid[cornery,cornerx+1] + (1-dx)**2*grid[cornery,cornerx])/tx
    b = (dx**2*grid[cornery+1,cornerx+1] + (1-dx)**2*grid[cornery+1,cornerx])/tx
    return (dy**2*b + (1-dy)**2*a)/ty

def visible(grid,startx,starty,x,y,elevation=100,isOffset=True):
    vis = np.full_like(x,1,"float64")
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
    for i in range(1,int(np.amax(norm))):
        m = (norm >= i+1)&(vis==1) # equivalent to a non-inclusive upper bound of int(norm)
        # should this be the case? or should it just be >= i?
        thisx += dx
        thisy += dy
        thish += dh
        h = quadHeight(grid,thisx[m],thisy[m]) 
        w = h >= thish[m]
        vis[[a[w] for a in np.where(m)]] = -1
    return vis

# y passed in world coordinates so grid coordinates are height-1-y
def viewshed(grid,pointx,pointy,elevation=100,isOffset=True,maxRange=3000,gridsize=30):
    gheight, gwidth = grid.shape
    
    pointy = gheight - pointy - 1
    
    view = np.full_like(grid,-1,int)
    maxSquaredCells = (maxRange/gridsize)**2

    ys, xs = np.indices(grid.shape,float)

    squareDist = (xs-pointx)**2+(ys-pointy)**2
    inRange = squareDist < maxSquaredCells
    
    view[inRange] = visible(grid,pointx,pointy,xs[inRange],ys[inRange],elevation,isOffset)
    
    return view
  

