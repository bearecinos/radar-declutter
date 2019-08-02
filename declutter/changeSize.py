"""Provides a method to resample the arrays stored in maps.hdf5."""
import h5py

def resize(cellsize): 
    """Changes the size of the cells in an existing numpy array. The lower left
    corner is unchanged i.e. array[-1,0].
    
    Parameters
    cellsize float : The size of each new cell in metres. This must be
        a multiple of the original cell size.

    Returns
    0 if successful. -1 otherwise (due to incorrect cell size)."""
    with h5py.File("maps.hdf5","r+") as f:
        line = f["meta"][()]
        originalSize = line[-1]
        if cellsize % originalSize != 0:
            print "Please use a multiple of the original size: "+str(originalSize)
            return -1
        print "Changing from cell size of "+str(originalSize)+" to "+str(cellsize)
        factor = int(cellsize/originalSize)
        f["heightmap"] = f["heightmap"][::-factor,::factor][::-1]
        f["aspect"] = f["aspect"][::-factor,::factor][::-1]
        f["slope"] = f["slope"][::-factor,::factor][::-1]

        line[-1] = cellsize
        f["meta"] = line
    return 0
    

    
