import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py
from declutter import models
from declutter import modelling
import os

###########################
# Takes point data (requires vis generated) and tracks bearing of
# wall. Currently this means nearest surface within 50m height of radar
# and horizon angle between 30 to 150 degrees to radar with smallest
# horizontal incidence close by (bearing incidence).
# May have to change to only consider fixed areas.
# 3D plot or just line.
# Try using below method to show which point it is first. Check reasonable
def passingAngle(fileName):
    pass


############################
# Interactive way of showing response from surfaces around
# points on path and stepping through path.
# May need graphics library rather than just pyplot.
def flyBy(pathName, steps = 10, timeInterval = [0,1.5e-5]):
    pass

def pathOnSurface(dirname):
    global fig, ax
    #fig=plt.figure(figsize=(10,7))
    #ax=fig.gca(projection='3d')
    #plt.subplots_adjust(0,0,1,1)
    xs,ys = [],[]

    grid, left, low, cellSize = loadMap()
    height = grid.shape[0]
    files = os.listdir(dirname)
    i = 0.0
    for f in files:
        i += 1.0
        visible,visCorner,distance,angle,_,_,pointx,pointy,pointz,antDir = loadPoint(dirname+"/"+f)
        xs.append(pointx)
        ys.append(pointy)
        #drawPoint(pointx,pointy,pointz,antDir,50)

        # mark nearest surface for each point
        idx = np.argmin(distance)
        inds = [a[idx] for a in np.where(visible)]
        print inds
        inds[1] = visible.shape[0]-1-inds[1]
        inds = [inds[0]*cellSize+visCorner[0],inds[1]*cellSize+visCorner[1]]

        print (pointx-visCorner[0])/cellSize
        print np.amin(np.hypot(np.where(visible)[0]-605,np.where(visible)[1]-605))
        print i
        if i == 6:
            plt.close()
            plt.imshow(visible)
            plt.figure()
            a = np.full_like(visible,0.0,float)
            a[visible] = distance
            a[a>5000] = np.nan
            print np.sum(a < 0)
            a[a < 0] = np.nan
            plt.imshow(a)
            print "here"
            
            #plt.colorbar()
            plt.show()
            return
            print inds
            print pointx
        else:
            continue
        
        z = grid[int(height-1-(inds[1]-low)/cellSize),int((inds[0]-left)/cellSize)]
        c = i/len(files)
        #ax.scatter([inds[0]],[inds[1]],[z],color=[c,0,1-c],linewidths=2,zorder=20)
        ax.plot([pointx,inds[0]],[pointy,inds[1]],[pointz,z],color="r")
        
        
    ##### how far extra around bounds of path
    #extend = 1000
    extend = 100

    drawMesh(xs,ys,extend,True,0.5) 
        
    plt.show()
    plt.close()

def markSurfaces(filename,start,end):
    global fig,ax
    load = loadPoint(filename)
    visible,visCorner,distance,angle,theta,phi,pointx,pointy,pointz,antDir = load

    grid, left, low, cellSize = loadMap()

    height = grid.shape[0]

    l = int((visCorner[0]-left)/cellSize)
    r = int(l+visible.shape[1])
    lower = int(height - (visCorner[1]-low)/cellSize)
    upper = int(lower - visible.shape[0])
        
    grid = grid[upper:lower,l:r]
  
    if antDir is None:
        if theta is not None:
            height = grid.shape[0]
            pointCoords = [int((pointx-visCorner[0])/cellSize),int((pointy-visCorner[1])/cellSize)]
            
            for y in range(height-1-pointCoords[1]):
                if visible[y,pointCoords[0]]:
                    idx = (np.where(visible)[1] == pointCoords[0]) & (np.where(visible)[0] == y)
                    z = np.cos(theta[idx])[()]
                    y = np.sin(theta[idx])*np.sin(phi[idx])[()]
                    antDir = (-y,z)
                    break
                
    m = (models._time(distance) > start) & (models._time(distance) < end)
    intensityModel = modelling.backscatter.raySpecular
    intensity = np.full(grid.shape,0.0,float)
    which = tuple([a[m] for a in np.where(visible)])
    intensity[which] = intensityModel(angle[m])

    height = visible.shape[0]
    ys,xs = np.indices(visible.shape)

    ys = height - 1 - ys
    xs = visCorner[0] + xs*cellSize
    ys = visCorner[1]+ ys*cellSize
    m = intensity != 0
    print "{0} points in this interval".format(np.sum(m))
    print "Total intensity: {0}".format(np.sum(intensity))

    fig=plt.figure(figsize=(10,7))
    ax=fig.gca(projection='3d')
    plt.subplots_adjust(0,0,1,1)
    
    drawMesh([pointx],[pointy],end*1.8e8,True)
    drawPoint(pointx,pointy,pointz,antDir)
    p = ax.scatter(xs[m],ys[m],grid[m],c=intensity[m],cmap="jet")
    plt.colorbar(p,shrink=0.7)
    plt.show()
    plt.close()
    return

def drawPoint(x,y,z,antDir=None,length = 150):
    ax.scatter([x],[y],[z],color="k",linewidths=2,zorder=20)
    if antDir is not None:
        ax.plot([x,x+antDir[0]*length],[y,y+antDir[1]*length],[z,z],linewidth=2,color="k")

def drawMesh(x,y,distance,back=False,alpha=0.2):
    global ax,fig
    grid, left, low, cellSize = loadMap()
    height,width = grid.shape

    ys,xs = np.indices(grid.shape)
    ys = height-1-ys
    xs = left + cellSize*xs
    ys = low + cellSize*ys

    minCoords = (np.amin(x)-left)/cellSize, height-1-(np.amin(y)-low)/cellSize
    maxCoords = (np.amax(x)-left)/cellSize, height-1-(np.amax(y)-low)/cellSize
    #pointCoords = (x-left)/cellSize, height-1-(y-low)/cellSize
    dist = distance/cellSize
    grid = grid[int(maxCoords[1]-dist):int(minCoords[1]+dist),
                int(minCoords[0]-dist):int(maxCoords[0]+dist)]
    xs = xs[int(maxCoords[1]-dist):int(minCoords[1]+dist),
                int(minCoords[0]-dist):int(maxCoords[0]+dist)]
    ys = ys[int(maxCoords[1]-dist):int(minCoords[1]+dist),
                int(minCoords[0]-dist):int(maxCoords[0]+dist)]

    if not back:
        fig=plt.figure(figsize=(10,7))
        ax=fig.gca(projection='3d')
    ax.plot_wireframe(xs,ys,grid,rcount=30,ccount=30,zorder=1,alpha=alpha)
    ax.set_xlabel("x"),
    ax.set_ylabel("y")


##    ax.scatter([x],[y],[z],color="k",linewidths=2,zorder=20)
##    if antDir is not None:
##        ax.plot([x,x+antDir[0]*150],[y,y+antDir[1]*150],[z,z],linewidth=2,color="k")
    if not back:
        plt.show()

def loadMap():
    with h5py.File("maps.hdf5","r") as f:
        grid = f["heightmap"][()]
        left,low,cellSize = f["meta"][:3]
    return grid, left, low, cellSize

def loadPoint(filename):
    antDir = None
    with h5py.File(filename,"r") as f:
        if "visible" not in f:
            print "Only possible for points which stored visible array."
            return -1
        visible = f["visible"][()]
        if "corner" not in f:
            print "Please recreate point files, format has changed."
            return -1
        visCorner = f["corner"][()]
        distance = f["distance"][()]
        angle = f["incidence"][()]
        theta,phi = None, None
        if "antennaTheta" in f:
            theta = f["antennaTheta"][()]
            phi = f["antennaPhi"][()]
        pointx,pointy,pointz,antDir = f["meta"][:4]
    if antDir is not None:
        antDir = (np.sin(antDir*np.pi/180.0), np.cos(antDir*np.pi/180.0))
    return visible,visCorner,distance,angle,theta,phi,pointx,pointy,pointz,antDir

