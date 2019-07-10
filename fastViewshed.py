import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.signal as signal
from mpl_toolkits.mplot3d import Axes3D

# y,x = np.indices([50,50])
# fig = plt.figure()
# ax = fig.add_subplot(1, 3, 1, projection='3d')
# ax.plot_wireframe(x,y,res)


# To simplify initially, assume grid is n by n
# and have view point in grid coordinates e.g. 105.4,105.3 (decimal for accuracy)
# also assume outer bounds of grid never of interest so automatically not considered
# i.e. only do visible for grid[1:-1,1:-1]

#grid = np.full((210,210),0.0)

def bandpass_filter(data, lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return signal.lfilter(b, a, data)

##grid = np.random.rand(220,220)*50.0
##grid = bandpass_filter(grid,0,0.05,1).swapaxes(0,1)
##grid = bandpass_filter(grid,0,0.05,1)[10:,10:]

def quadHeight(x,y):
    cornerx = x.astype(int)
    cornery = y.astype(int)
    dx = x-cornerx
    tx = dx**2+(1-dx)**2
    dy = y-cornery
    ty = dy**2+(1-dy)**2
    a = (dx**2*grid[cornery,cornerx+1] + (1-dx)**2*grid[cornery,cornerx])/tx
    b = (dx**2*grid[cornery+1,cornerx+1] + (1-dx)**2*grid[cornery+1,cornerx])/tx
    return (dy**2*b + (1-dy)**2*a)/ty

def visible(x,y,elevation=100,isOffset=True): # need to flip y???
    vis = np.full_like(x,1,"float64")
    endh = quadHeight(x,y)
    if isOffset:
        starth = quadHeight(np.array([px]),np.array([py]))+elevation
    else:
        starth = elevation
    dx = x-px
    dy = y-py
    norm = (dx**2+dy**2)**0.5

    # ignores divide by 0 cases in coming loop
    norm[norm==0] = 0.1
    
    dx /= norm
    dy /= norm

    for i in range(1,int(np.amax(norm))):
        m = norm >= i+1 # equivalent to a non-inclusive upper bound of int(norm)
        thisx = px+i*dx[m]
        thisy = px+i*dy[m]
        h = quadHeight(thisx,thisy)
        w = h >= (i/norm[m])*endh[m]+(1.0-i/norm[m])*starth
        vis[[a[w] for a in np.where(m)]] = -1
    return vis

gridsize = 30
maxRange = 3000
px,py = 105.0, 105.0

def viewshed(inGrid,px,py,elevation=100,isOffset=True,maxrange=3000,size=30):
    global grid, maxRange, gridsize
    maxRange = maxrange
    gridsize = size
    grid = inGrid
    gheight, gwidth = grid.shape
    view = np.full_like(grid,-1,int)
    maxSquaredCells = (maxRange/gridsize)**2

    ys, xs = np.indices(grid.shape)
    #ys = gheight - ys

    squareDist = (xs-px)**2+(gheight-ys-py)**2
    inRange = squareDist < maxSquaredCells

    view[inRange] = visible(xs[inRange],ys[inRange],elevation,isOffset)
    
    return view
    
##oth = np.full_like(grid,False)
##for x in range(1,209):
##    for y in range(1,209):
##        oth[y][x] = visible(x,y)
##plt.subplot("121")
##plt.contourf(grid,100)
##plt.subplot("122")
##plt.imshow(oth)
##plt.show()
