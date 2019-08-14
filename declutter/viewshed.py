"""Determines which points on a raster are visible from a given viewpoint."""
import numpy as np
import math

_NODATA = np.nan

def quadHeight(grid,x,y):
    """Calculates the interpolated height of the surface below the given point.
    Quadratic interpolation in both axes.

    Parameters
    grid - 2D float array : A heightmap of the surface.
    x,y - float arrays :Indices of points on the grid. i.e. 1,2 refers to grid[2,1].
    """
    res = np.full_like(x,_NODATA,float) 
    h,w = grid.shape
    m = (x >= 0) & (y >= 0) & (x < w-1) & (y < h-1) # valid points mask
    x,y = x[m], y[m]
        
    cornerx = x.astype(int)
    cornery = y.astype(int)

    # any NaN values are carried through and picked up in calling function    
    dx = x-cornerx
    tx = dx*dx+(1-dx)*(1-dx)
    dy = y-cornery
    ty = dy*dy+(1-dy)*(1-dy)
    a = (dx*dx*grid[cornery,cornerx+1] + (1-dx)*(1-dx)*grid[cornery,cornerx])/tx
    b = (dx*dx*grid[cornery+1,cornerx+1] + (1-dx)*(1-dx)*grid[cornery+1,cornerx])/tx
    res[m] = (dy*dy*b + (1-dy)*(1-dy)*a)/ty
    return res

def _visible(grid,startx,starty,x,y,elevation=100,isOffset=True,stepSize=1.0):
    vis = np.full_like(x,True,bool)
    if len(x) == 0:
        return vis # else get error when attempting to call np.amax on 0 length array
    endh = quadHeight(grid,x,y)
    vis[np.isnan(endh)] = False # won't be detected later as this will be NaN
    if isOffset: # have to wrap as an array and unpack again
        starth = quadHeight(grid,np.array([startx]),np.array([starty]))[0]+elevation
    else:
        starth = elevation
    dx = x-startx
    dy = y-starty
    norm = np.hypot(dx,dy)

    # avoids divide by 0 cases in coming loop
    # anything less than stepSize is automatically visible
    norm[norm==0] = stepSize/2.0
    
    dx *= stepSize/norm
    dy *= stepSize/norm
    dh = (endh-starth)*stepSize/norm

    thisx = startx
    thisy = starty
    thish = starth
    for i in range(1,int(np.amax(norm)/stepSize)):
        m = (norm >= (i+1)*stepSize)&(vis==1) 
        thisx += dx
        thisy += dy
        thish += dh
        h = quadHeight(grid,thisx[m],thisy[m])
        # if height unknow get NaN so also mark as not visible
        w = (h >= thish[m]) | np.isnan(h)
        vis[tuple(a[w] for a in np.where(m))] = False
    return vis


# assumes point actually inside grid (checked by call to quadheight to get groundHeight first in stateless.py)
def viewshed(grid,pointx,pointy,mask,elevation=100,isOffset=True,gridsize=30.0,stepSize=None):
    """Calculates which points on a grid are visible from a given viewpoint.
    If the ray from the viewpoint to a surface passes over a part of the grid where the height
    is unknown (NaN), that surface is treated as not being visible. This should only affect
    the borders of the grid. This also applies to points outside the grid.

    Parameters
    grid - 2D float array : A heightmap of the surface.
    pointx, pointy - floats : Coordinates of the viewpoint in cell size units from the lower left corner.
        i.e. [1,1] refers to grid[-2,1].
    mask - 2D bool array : A mask of the points to consider. Everything else is assumed to be not visible.
        e.g. because those points are known to be out of range or have back facing surfaces.
    elevation - float (optional) : The elevation of the viewpoint. By default this is 100m.
    isOffset - bool (true) : Whether the elevation is relative to the surface beneath the radar or not. Default True.
    gridsize - float (optional) : The resolution of the grid. Default is 30m. Only important if stepSize is set.
    stepSize - float (optional) : The smallest resolution steps to take along the ray to each surface point to determine
        if it is visible. This is in world coordinates i.e. stepSize = 30 and gridSize = 30 means the smallest step
        is 1 cell at a time.

    Returns
    view - 2D bool array : A mask where only points with value True are visible.
    """
    gheight, gwidth = grid.shape
    
    view = np.full_like(grid,False,bool)

    ys, xs = np.indices(grid.shape,float)

    if stepSize is None:
        stepSize = gridsize
    stepSize /= gridsize # step in grid coordinates

    m = np.where(mask)
    
    # step size decreases in multiples of 8 to reduce overall work.
    # factor chosen by experiment, checking total number of points
    # tested for each factor.
    scaleUp = int(math.log(gheight/(2.0*stepSize),8))
    for s in range(scaleUp,-1,-1):
        result = _visible(grid,pointx,pointy,xs[m],ys[m],elevation,isOffset,stepSize*8.0**s).astype(bool)
        view[m] = result # clears points known to not be visible
        m = tuple(a[result] for a in m) # smaller mask for next loop

    return view

