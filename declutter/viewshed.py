import numpy as np
import math

_NODATA = np.nan

def quadHeight(grid,x,y): # bi-quadratic interpolation of height
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

def visible(grid,startx,starty,x,y,elevation=100,isOffset=True,stepSize=1.0):
    vis = np.full_like(x,True,bool)
    if len(x) == 0:
        return vis # else get error when attempting to call np.amax on 0 length array
    endh = quadHeight(grid,x,y)
    vis[np.isnan(endh)] = False # won't be detected later as thish will be NaN
    if isOffset:
        starth = quadHeight(grid,np.array([startx]),np.array([starty]))+elevation
    else:
        starth = elevation
    dx = x-startx
    dy = y-starty
    norm = np.hypot(dx,dy)
    #norm = (dx**2+dy**2)**0.5

    # avoids divide by 0 cases in coming loop
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

cost = 0
# y passed in world coordinates so grid coordinates are height-1-y
# assumes point actually inside grid (checked by call to quadheight to get groundHeight first in stateless.py)
def viewshed(grid,pointx,pointy,mask,elevation=100,isOffset=True,gridsize=30.0,stepSize=None):
    global cost
    gheight, gwidth = grid.shape
    
    pointy = gheight - pointy - 1
    
    view = np.full_like(grid,False,bool)

    ys, xs = np.indices(grid.shape,float)

    if stepSize is None:
        stepSize = gridsize
    stepSize /= gridsize # step in grid coordinates

    m = np.where(mask)
    
    # step size decreases in multiples of 8 to reduce overall work
    scaleUp = int(math.log(gheight/(2.0*stepSize),8))
    for s in range(scaleUp,-1,-1):
        result = visible(grid,pointx,pointy,xs[m],ys[m],elevation,isOffset,stepSize*8.0**s).astype(bool)
        view[m] = result # clears points known to not be visible
        m = tuple(a[result] for a in m) # smaller mask for next loop

    return view


if __name__=="__main__" and False:
    import timeit
    grid = np.full((1000,1000),0,float)
    mask = np.full((1000,1000),0,bool)
    px = 499.9
    py = 500.1
    t = timeit.Timer(lambda : viewshed(grid,px,py,mask,100,False,gridsize=2.0))
    print t.timeit(500)

