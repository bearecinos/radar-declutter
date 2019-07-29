import h5py

def resize(cellsize): # keeps bottom right corner - a[-1:0] unchanged
    """Changes size of an existing numpy array."""
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
    

    
