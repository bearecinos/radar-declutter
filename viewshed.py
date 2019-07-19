import numpy as np


_NODATA = np.nan

# should ignore points outside grid, and any point interpolating a NODATA value
# should also be ignored
# expects points in grid coordinates
def quadHeight(grid,x,y): # bi-quadratic interpolation of height
    res = np.full_like(x,_NODATA,float) 
    h,w = grid.shape
    m = (x >= 0) & (y >= 0) & (x < w-1) & (y < h-1) # valid points mask
    x,y = x[m], y[m]
        
    cornerx = x.astype(int)
    cornery = y.astype(int)
    # output NODATA if any cell to interpolate between is already NODATA
    # assumes use of NaN
    w = ~(np.isnan(grid[cornery,cornerx]) | np.isnan(grid[cornery,cornerx+1]) |
          np.isnan(grid[cornery+1,cornerx]) | np.isnan(grid[cornery+1,cornerx+1]))
##    w = ((grid[cornery,cornerx] != _NODATA) & (grid[cornery+1,cornerx] != _NODATA) &
##        (grid[cornery,cornerx+1] != _NODATA) & (grid[cornery+1,cornerx+1] != _NODATA))

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
    vis = np.full_like(x,1,int)
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
        vis[tuple(a[w] for a in np.where(m))] = 0
    return vis

# y passed in world coordinates so grid coordinates are height-1-y
# assumes point actually inside grid (checked by call to quadheight to get groundHeight first in stateless.py)
def viewshed(grid,pointx,pointy,mask,elevation=100,isOffset=True,maxRange=3000,gridsize=30.0,stepSize=None):
    gheight, gwidth = grid.shape
    
    pointy = gheight - pointy - 1
    
    view = np.full_like(grid,0,int)

    ys, xs = np.indices(grid.shape,float)
    
    if stepSize is None:
        stepSize = gridsize
    stepSize /= gridsize # step in grid coordinates
    
    view[mask] = visible(grid,pointx,pointy,xs[mask],ys[mask],elevation,isOffset,stepSize)
    
    return view
  

