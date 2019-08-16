"""Provides methods for reducing the size of the maps.hdf5 file."""
import h5py
from modelling import parameters
import path
import numpy as np
from version import version

def pathCrop(pathName, crop = [0,0]):
    '''Removes areas of the arrays which are not needed to process the given
    path. crop = [A,B] also avoids considering the first A and last B points
    of the path. Leaves 2 extra cell-widths padding to avoid boundary/rounding
    issues later.'''
    d = parameters.env.getMaxDist()
    xs,ys,_ = path.loadData(pathName,crop)
    with h5py.File("maps.hdf5","r") as f:
        xmin,ymin,cellsize = f["meta"][()]
        d += 2*cellsize
        print "Cropping to a range of {0}m from path.".format(d)
        hmap = f["heightmap"][()]
        height,width = hmap.shape

        xBounds = [max(0,int((np.amin(xs)-d-xmin)/cellsize)),
                   min(width,int((np.amax(xs)+d-xmin)/cellsize)+1)]
        
        yBounds = [max(0,int((np.amin(ys)-d-ymin)/cellsize)),
                   min(height,int((np.amax(ys)+d-ymin)/cellsize)+1)]
        slope = f["slope"][yBounds[0]:yBounds[1],xBounds[0]:xBounds[1]]
        aspect = f["aspect"][yBounds[0]:yBounds[1],xBounds[0]:xBounds[1]]
        
    with h5py.File("maps.hdf5","w") as f:
        f.create_dataset("heightmap",compression="gzip",data = hmap[yBounds[0]:yBounds[1],xBounds[0]:xBounds[1]])
        f.create_dataset("slope",compression="gzip",data = slope)
        f.create_dataset("aspect",compression="gzip",data = aspect)
        f["meta"] = np.array([xmin+xBounds[0]*cellsize, ymin+yBounds[0]*cellsize,cellsize])
        f["version"] = version
    return 0


def resize(cellsize): 
    """Changes the size of the cells in an existing numpy array. The lower left
    corner is unchanged i.e. array[0,0].
    
    Parameters
    ----------
    cellsize - float : The size of each new cell in metres. This must be
        a multiple of the original cell size.

    Returns
    -------
    0 if successful. -1 otherwise (due to incorrect cell size)."""
    with h5py.File("maps.hdf5","r") as f:
        line = f["meta"][()] # min-x coord, min-y coord, cellSize
        originalSize = line[-1]
        if cellsize % originalSize != 0: # not a multiple
            print "Please use a multiple of the original size: "+str(originalSize)
            return -1
        print "Changing from cell size of "+str(originalSize)+" to "+str(cellsize)
        factor = int(cellsize/originalSize)
        
        heightmap = f["heightmap"][::factor,::factor]
        aspect = f["aspect"][::factor,::factor]
        slope = f["slope"][::factor,::factor]

        # update metaData with new cell size
        line[-1] = cellsize
        
    with h5py.File("maps.hdf5","w") as f:
        f.create_dataset("heightmap", compression="gzip", data = heightmap)
        f.create_dataset("slope", compression="gzip", data = slope)
        f.create_dataset("aspect", compression="gzip", data = aspect)
        f["meta"] = line
        f["version"] = version
    return 0
    

    
