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

def getHeight(x,y):
    cornerx = int(x)
    cornery = int(y)
    num = 0
    den = 0
    for cx in range(2):
        for cy in range(2):
            d = ((cornerx+cx-x)**2+(cornery+cy-y)**2)**0.5
            if d == 0:
                return grid[cornerx+cx,cornery+cy]
            num += grid[cornery+cy,cornerx+cx]/d
            den += 1/d
    return num/den

def linHeight(x,y):
    cornerx = int(x)
    cornery = int(y)
    dx = x-cornerx
    dy = y-cornery
    a = dx*grid[cornery][cornerx+1] + (1-dx)*grid[cornery][cornerx]
    b = dx*grid[cornery+1][cornerx+1] + (1-dx)*grid[cornery+1][cornerx]
    return dy*b + (1-dy)*a

def quadHeight(x,y):
    cornerx = int(x)
    cornery = int(y)
    dx = x-cornerx
    tx = dx**2+(1-dx)**2
    dy = y-cornery
    ty = dy**2+(1-dy)**2
    a = (dx**2*grid[cornery][cornerx+1] + (1-dx)**2*grid[cornery][cornerx])/tx
    b = (dx**2*grid[cornery+1][cornerx+1] + (1-dx)**2*grid[cornery+1][cornerx])/tx
    return (dy**2*b + (1-dy)**2*a)/ty

def visible(x,y,elevation=100,isOffset=True):
    endh = quadHeight(x,y)
    if isOffset:
        starth = quadHeight(px,py)+elevation
    else:
        starth = elevation
    dx = x-px
    dy = y-py
    norm = (dx**2+dy**2)**0.5
    if norm == 0:
        return True
    dx /= norm
    dy /= norm
    for i in range(1,int(norm)):
        thisx = px+i*dx
        thisy = px+i*dy
        h = quadHeight(thisx,thisy)
        if h >= (i/norm)*endh+(1.0-i/norm)*starth:
            return -1
    return 1

gridsize = 30
maxRange = 3000
px,py = 105.0, 105.0
def viewshed(inGrid,px,py,elevation=100,isOffset=True,maxrange=3000,size=30):
    global grid, maxRange, gridsize
    maxRange = maxrange
    gridsize = size
    grid = inGrid
    gheight, gwidth = grid.shape
    xmin = max(0,int(px-maxRange/gridsize))
    ymin = max(0,int(py-maxRange/gridsize))
    xmax = min(int(px+maxRange/gridsize+1),gwidth)
    ymax = min(int(py+maxRange/gridsize+1),gwidth)
    view = np.full_like(grid,-1,int)
    maxSquaredCells = (maxRange/gridsize)**2

    #ys, xs = np.indices(grid.shape)
    #ys = gheight - ys

    #squareDist = (xs-px)**2+(ys-py)**2
    #inRange = sqaureDist < maxSquaredCells
    
    for x in range(xmin,xmax):
        for y in range(ymin,ymax):
            if (x-px)**2+(gheight-y-py)**2 < maxSquaredCells:
                view[y][x] =  visible(x,y,elevation,isOffset) 
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
