import sys
import path
import models
import makeArrays
def main(args=None):
    if args is None:
        args = sys.argv[1:]
    print args
    if len(args) < 2:
        print "Wrong number of arguments. Options: "
        return showHelp()
    # load rasters into numpy format
    if args[0] == "load":
        filename = args[1]
        lat,lon = float(args[2]),float(args[3])
        out,cellSize = None,None
        if len(args) == 4: # no options
            pass 
        elif len(args) == 6: # 1 option
            if args[4] == "-c":
                cellSize = float(args[5])
            elif args[4] == "-d":
                out = args[5]
            else:
                return badOptions()
        elif len(args) == 8: # 2 options
            if args[4] == "-c":
                cellSize = float(args[5])
                if args[6] == "-d":
                    out = args[7]
                else:
                    return badOptions() 
            elif args[4] == "-d":
                out = args[5]
                if args[6] == "-c":
                    cellSize = float(args[7])
                else:
                    return badOptions() 
            else:
                return badOptions()
        else:
            print "Wrong number of arguments. Options: "
            return showHelp()
        makeArrays.makeAll(filename,lat,lon,cellSize,out)
    # process radar path and display/return result
    elif args[0] == "model":
        filename =  args[1]
        if len(args) == 4:
            if args[2] == "-d":
                out = args[3]
            else:
                return badOptions()
        elif len(args) > 4:
            print "Too many optional arguments. Options:"
            return showHelp()
        else:
            out = filename[:-4]+"Out"
        path.loadGpx(filename,[0,0],out)
        return models.compare(out)
    else:
        print "Command not recognised. Options: "
        return showHelp()

def badOptions():
    print "Optional arguments not recognised. Options: "
    showHelp()

def showHelp():
    print "load - generates numpy arrays from a given raster."
    print "    arguments: name latitude longitude (-c cellSize) (-d outDir) "
    print "model - simluates the response for a given path (from a gpx).:"
    print "    arguments: filename (-d outDir)"

if __name__=="__main__":
    main()
