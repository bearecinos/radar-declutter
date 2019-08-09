import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py
from declutter import radar
from declutter import modelling
import os

area = 10
steps = 5
def smoothness(filename):
    grid, left, low, cellSize = loadMap()
    visible,visCorner,_,_,_,_,_,_,_,_ = loadPoint(filename)
    
    
    results = np.full((visible.shape[0]/steps,visible.shape[1]/steps),0.0,float)
    height = grid.shape[0]
    l = int((visCorner[0]-left)/cellSize)
    r = int(l+visible.shape[1])
    lower = int((visCorner[1] - low)/cellSize)
    upper = int(lower + visible.shape[0])
    # order reversed now
##    lower = int(height - (visCorner[1]-low)/cellSize)
##    upper = int(lower - visible.shape[0])
    grid = grid[lower:upper,l:r]
    for x in range(0,visible.shape[1]-area,steps):
        for y in range(0,visible.shape[0]-area,steps):
                ybox,xbox = np.indices((area,area))
                ybox = ybox.reshape(-1)
                xbox = xbox.reshape(-1)
                gbox = grid[y:y+area,x:x+area].reshape(-1,)
                if np.sum(np.isnan(gbox))>0:
                        results[y/steps,x/steps] = np.nan
                        continue
                A = np.vstack([xbox,ybox,np.ones(len(xbox))]).T
                residuals = np.linalg.lstsq(A,gbox)[1]
                results[y/steps,x/steps] = 1.0 - residuals/np.sum((gbox-np.mean(gbox))**2)

    
    return results 

def spread(result,visible,cutoff=0.95):
    out = np.full_like(visible,0,int)
    y,x = np.indices(visible.shape)
    out[y,x] = result[y/5,x/5] > cutoff
    return out

# returns direction of glacier in terms of theta = arctan(dy/dx)
# grid should already be cropped
# can ignore points where line would be within 30 degrees of approximated direction?
# error of at least 20 degrees
# doesn't help remove points on edges of glacier
def glacierDir(grid,groundHeight):
    m = (grid < groundHeight + 30) & (grid > groundHeight-30)
    ys,xs = np.indices(m.shape)[:,m]
    grad = np.polyfit(xs,ys,1)[0] # 1st order coefficient
    return np.arctan(-grad) # account for world y coordinates being in other direction

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
def flyBy(dirname, above=False, stepsize=3):
    global fig, ax
    fig=plt.figure(figsize=(10,7))
    ax=fig.gca(projection='3d')
    plt.subplots_adjust(0,0,1,1)
    xs,ys,zs = [], [],[]
    surfxs, surfys, surfzs, aspects = [],[],[],[]
    pointCols = []
    grid, left, low, cellSize = loadMap()
    with h5py.File("maps.hdf5","r") as f:
        slope = f["slope"][()]
        aspect = f["aspect"][()]

    height = grid.shape[0]
    files = os.listdir(dirname)
    files.sort(key=lambda x : int(x[5:-5]))

    cols = plt.cm.jet
    
    i = 0.0
    for f in files:
        i += 1.0
        visible,visCorner,distance,angle,_,_,pointx,pointy,pointz,antDir = loadPoint(dirname+"/"+f)
        if len(distance) == 0:
            continue # invalid point
        xs.append(pointx)
        ys.append(pointy)
        zs.append(pointz)
        #drawPoint(pointx,pointy,pointz,antDir,50)

        s = matchMap(slope,left,low,visCorner[0],visCorner[1],visible.shape[1],visible.shape[0],cellSize)
        a = matchMap(aspect,left,low,visCorner[0],visCorner[1],visible.shape[1],visible.shape[0],cellSize)
        g = matchMap(grid,left,low,visCorner[0],visCorner[1],visible.shape[1],visible.shape[0],cellSize)
        
        dx = int((visCorner[0]-left)/cellSize)
        dy = int((visCorner[1] - low)/cellSize)
        
        
        # way of detecting wall (may be very data dependent)
        # closes point which is : visible, angle less than 5 degrees to horizon from radar,
        #                         at least 400m away, at least 20m above radar
        if above:
            offset = 20
        else:
            offset = 0
        idx = np.argmin(distance + (abs((g[visible]-pointz)/distance) > np.sin(5*np.pi/180.0))*5000
                        + (distance < 400)*5000 + (g[visible] < pointz+offset)*5000)
        
        # y then x
        inds = [a[idx] for a in np.where(visible)]

        #inds[0] = visible.shape[0]-1-inds[0] grid reversed now so [0,0] is corner
        inds = [inds[0]*cellSize+visCorner[1],inds[1]*cellSize+visCorner[0]]
        
        #z = grid[int(height-1-(inds[0]-low)/cellSize),int((inds[1]-left)/cellSize)]
        z = grid[int((inds[0]-low)/cellSize),int((inds[1]-left)/cellSize)]
        aspects.append(aspect[int((inds[0]-low)/cellSize),int((inds[1]-left)/cellSize)])
        c = cols(i/len(files))

        surfxs.append(inds[1])
        surfys.append(inds[0])
        surfzs.append(z)
        pointCols.append(c)

    surfxs, surfys, surfzs = np.array(surfxs), np.array(surfys), np.array(surfzs)
    xs, ys, zs = np.array(xs), np.array(ys), np.array(zs)
    pointCols, aspects = np.array(pointCols), np.array(aspects)
    groups = np.arange(len(pointCols))
    #ax.scatter(surfxs,surfys,surfzs,color="r",linewidths=2)
        
    ##### how far extra around bounds of path
    #extend = 1000
    extend = 700


    # PLOTTING SURFACE
    height,width = grid.shape
    xcoords = (np.amin(xs)-left)/cellSize, (np.amax(xs)-left)/cellSize
    xcoords = max(0,int(xcoords[0]-extend/cellSize)), min(width,int(xcoords[1]+extend/cellSize))
    ycoords = (np.amin(ys)-low)/cellSize, (np.amax(ys)-low)/cellSize
    ycoords = max(0,int(ycoords[0]-extend/cellSize)), min(height,int(ycoords[1]+extend/cellSize))
    grid = grid[ycoords[0]:ycoords[1],xcoords[0]:xcoords[1]]
    Y,X = np.indices(grid.shape)
    Y += ycoords[0]
    X += xcoords[0]

    size = grid.size
    if size > 200**2:
        factor = int(np.ceil(size**0.5/200.0))
        grid = grid[::factor,::factor]
        X,Y = X[::factor,::factor], Y[::factor,::factor]
    X = left + X*cellSize
    Y = low + Y*cellSize
    
    my_col = plt.cm.terrain((grid-np.nanmin(grid))/(np.nanmax(grid)-np.nanmin(grid)))
    ax.plot_surface(X,Y,grid,facecolors=my_col,linewidth=0,antialiased=False,zorder=3)
 
    ## / PLOTTING SURFACE
    
    # Groups
    distances = []
    for i in range(len(xs)-1):
        distances.append(((surfxs[i]-surfxs[i+1])**2+(surfys[i]-surfys[i+1])**2+(surfzs[i]-surfzs[i+1])**2)**0.5)
    distances = np.array(distances)

    lengths = np.hypot(np.hypot(xs-surfxs,ys-surfys),zs-surfzs)

    for i in range(len(xs)-1):
        if distances[i] < 800: # usure about choice or if better way to check (e.g. opposite direction to surface)
            groups[i+1] = groups[i]
    pointCols = pointCols[groups]
    
    for g in np.unique(groups):
        i = np.argmin(lengths[groups == g])
        j = np.where(groups == g)[0][0]
        k = np.where(groups == g)[0][-1]+1
        surfxs[j:k] = surfxs[j+i]
        surfys[j:k] = surfys[j+i]
        surfzs[j:k] = surfzs[j+i]
        aspects[j:k] = aspects[j+i]
    # /Groups


    # surfxs, surfys, surfzs gives surface points corresponding to xs, ys, zs
    ## Points of path
    ax.scatter(xs[::stepsize],ys[::stepsize],zs[::stepsize],color="k",linewidths=2,zorder=10)
    # plotting lines to surfaces
    for i in range(0,len(xs),stepsize):
        ax.plot([xs[i],surfxs[i]],[ys[i],surfys[i]],[zs[i],surfzs[i]],color=pointCols[i],zorder=5)

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    

    # finding incidence
    aspects *= np.pi / 180.0 # convert to radians
    incidence = np.array([(xs[i]-surfxs[i])*np.sin(aspects[i])+(ys[i]-surfys[i])*np.cos(aspects[i]) for i in range(len(xs))])
    incidence /= np.hypot([(xs[i]-surfxs[i]) for i in range(len(xs))],[(ys[i]-surfys[i]) for i in range(len(xs))])

    incidence = np.arccos(incidence)*180.0/np.pi
    # / calculating incidence

    plt.show()
    plt.close()

    return incidence, np.unique(groups)
#    return distances, lengths, groups

    

def matchMap(grid,left,low,cropLeft,cropLow,w,h,cellSize):
    l = int((cropLeft-left)/cellSize)
    r = int(l+w)
    lower = int((cropLow-low)/cellSize)
    upper = int(lower+h)
##    lower = int(height - (cropLow-low)/cellSize)
##    upper = int(lower - h)
##    return grid[upper:lower,l:r]
    return grid[lower:upper,l:r]


def pathOnSurface(dirname,twoD = False):
    global fig, ax
    if not twoD:
        fig=plt.figure(figsize=(10,7))
        ax=fig.gca(projection='3d')
        plt.subplots_adjust(0,0,1,1)
    xs,ys,zs = [], [],[]
    surfxs, surfys, surfzs, aspects = [],[],[],[]
    grid, left, low, cellSize = loadMap()
    with h5py.File("maps.hdf5","r") as f:
        slope = f["slope"][()]
        aspect = f["aspect"][()]

    height = grid.shape[0]
    files = os.listdir(dirname)
    files.sort(key=lambda x : int(x[5:-5]))

    cols = plt.cm.jet
    
    i = 0.0
    for f in files:
        i += 1.0
        visible,visCorner,distance,angle,_,_,pointx,pointy,pointz,antDir = loadPoint(dirname+"/"+f)
        if len(distance) == 0:
            continue # invalid point
        xs.append(pointx)
        ys.append(pointy)
        zs.append(pointz)
        #drawPoint(pointx,pointy,pointz,antDir,50)

        s = matchMap(slope,left,low,visCorner[0],visCorner[1],visible.shape[1],visible.shape[0],cellSize)
        a = matchMap(aspect,left,low,visCorner[0],visCorner[1],visible.shape[1],visible.shape[0],cellSize)
        g = matchMap(grid,left,low,visCorner[0],visCorner[1],visible.shape[1],visible.shape[0],cellSize)
        
        dx = int((visCorner[0]-left)/cellSize)
        dy = int((visCorner[1] - low)/cellSize)
        
        
        # way of detecting wall (may be very data dependent)
        # closes point which is : visible, angle less than 5 degrees to horizon from radar,
        #                         at least 400m away, at least 20m above radar
        idx = np.argmin(distance + (abs((g[visible]-pointz)/distance) > np.sin(5*np.pi/180.0))*5000
                        + (distance < 400)*5000 + (g[visible] < pointz+20)*5000)

        # y then x
        inds = [a[idx] for a in np.where(visible)]

        #inds[0] = visible.shape[0]-1-inds[0] grid reversed now so [0,0] is corner
        inds = [inds[0]*cellSize+visCorner[1],inds[1]*cellSize+visCorner[0]]
        
        #z = grid[int(height-1-(inds[0]-low)/cellSize),int((inds[1]-left)/cellSize)]
        z = grid[int((inds[0]-low)/cellSize),int((inds[1]-left)/cellSize)]
        aspects.append(aspect[int((inds[0]-low)/cellSize),int((inds[1]-left)/cellSize)])
        c = cols(i/len(files))

        if twoD:
            plt.plot([pointx,inds[1]],[pointy,inds[0]],color=c)
        else:
            ax.plot([pointx,inds[1]],[pointy,inds[0]],[pointz,z],color=c)


        surfxs.append(inds[1])
        surfys.append(inds[0])
        surfzs.append(z)

    if twoD:
        plt.scatter(xs,ys,color="k",linewidths=2)
    else:
        ax.scatter(xs,ys,zs,color="k",linewidths=2)
    #ax.scatter(surfxs,surfys,surfzs,color="r",linewidths=2)
        
    ##### how far extra around bounds of path
    #extend = 1000
    extend = 700

    if not twoD:
        drawMesh(xs,ys,extend,True,0.5) 
        
    plt.show()
    plt.close()

    incidence = np.array([(xs[i]-surfxs[i])*np.sin(aspects[i])+(ys[i]-surfys[i])*np.cos(aspects[i]) for i in range(len(xs))])
    incidence /= np.hypot([(xs[i]-surfxs[i]) for i in range(len(xs))],[(ys[i]-surfys[i]) for i in range(len(xs))])

    incidence = np.arccos(incidence)*180.0/np.pi
    return incidence

def markSurfaces(filename,start,end):
    global fig,ax
    load = loadPoint(filename)
    visible,visCorner,distance,angle,theta,phi,pointx,pointy,pointz,antDir = load

    grid, left, low, cellSize = loadMap()

    height = grid.shape[0]

    l = int((visCorner[0]-left)/cellSize)
    r = int(l+visible.shape[1])
##    lower = int(height - (visCorner[1]-low)/cellSize)
##    upper = int(lower - visible.shape[0])
    lower = int((visCorner[1]-low)/cellSize)
    upper = int(lower + visible.shape[0])
        
##    grid = grid[upper:lower,l:r]
    grid = grid[lower:upper,l:r]
  
    if antDir is None:
        if theta is not None:
            height = grid.shape[0]
            pointCoords = [int((pointx-visCorner[0])/cellSize),int((pointy-visCorner[1])/cellSize)]
            
##            for y in range(height-1-pointCoords[1]): - now really want to step backwards from top
            for y in range(height-1,pointCoords[1],-1):
                if visible[y,pointCoords[0]]:
                    idx = (np.where(visible)[1] == pointCoords[0]) & (np.where(visible)[0] == y)
                    z = np.cos(theta[idx])[()]
                    y = np.sin(theta[idx])*np.sin(phi[idx])[()]
                    antDir = (-y,z)
                    break
                
    m = (radar._time(distance) > start) & (radar._time(distance) < end)
    intensityModel = modelling.backscatter.raySpecular
    intensity = np.full(grid.shape,0.0,float)
    which = tuple([a[m] for a in np.where(visible)])
    intensity[which] = intensityModel(angle[m])

    height = visible.shape[0]
    ys,xs = np.indices(visible.shape)

    #ys = height - 1 - ys ## map reversed so not needed
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

    # map reversed so now corner is [0,0]
##    minCoords = (np.amin(x)-left)/cellSize, height-1-(np.amin(y)-low)/cellSize
##    maxCoords = (np.amax(x)-left)/cellSize, height-1-(np.amax(y)-low)/cellSize
    minCoords = (np.amin(x)-left)/cellSize, (np.amin(y)-low)/cellSize
    maxCoords = (np.amax(x)-left)/cellSize, (np.amax(y)-low)/cellSize
    #pointCoords = (x-left)/cellSize, height-1-(y-low)/cellSize
    dist = distance/cellSize
    grid = grid[int(minCoords[1]-dist):int(maxCoords[1]+dist),
                int(minCoords[0]-dist):int(maxCoords[0]+dist)]

    ys,xs = np.indices(grid.shape)
##    ys = grid.shape[0]-1-ys # now grid reversed so not needed
    
    xs = left + cellSize*(xs+int(minCoords[0]-dist))
   ## ys = low + cellSize*(ys+height-int(minCoords[1]+dist))
    ys = low + cellSize*(ys+int(minCoords[1]-dist))


    if not back:
        fig=plt.figure(figsize=(10,7))
        ax=fig.gca(projection='3d')
    ax.plot_wireframe(xs,ys,grid,rcount=30,ccount=30,zorder=1,alpha=alpha)
    ax.set_xlabel("x")
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

