import numpy as np


def resize(cellsize): # keeps bottom right corner - a[-1:0] unchanged
    """Changes size of an existing numpy array."""
    with open("info.txt","r") as f:
        line = f.read().split(",")
        originalSize = float(line[-1])
    if cellsize % originalSize != 0:
        print "Please use a multiple of the original size: "+str(originalSize)
        return -1
    print "Changing from cell size of "+str(originalSize)+" to "+str(cellsize)
    factor = int(cellsize/originalSize)
    f = np.load("maps.npz")
    heightmap = f["heightmap"][::-factor,::factor][::-1]
    aspect = f["aspect"][::-factor,::factor][::-1]
    slope = f["slope"][::-factor,::factor][::-1]

    line[-1] = str(cellsize)
    with open("info.txt","w") as f:
        f.write(",".join(line))
    np.savez_compressed("maps.npz",heightmap=heightmap,
                        aspect=aspect,slope=slope)
    

    
