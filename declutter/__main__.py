"""Provides a way of running the module by command line arguments, using the
three methods defined below. The -h option provides information on the whole
interface and also on each of the methods.
"""
import argparse

def load(args):
    """Converts a raster to numpy format in maps.hdf5, projecting and resampling
    if needed."""
    print "Loading..."
    import makeArrays
    return makeArrays.makeAll(args.filename,args.latitude,args.longitude,args.cellSize,args.out)

def model(args):
    """Produces a radargram from a gps path. The intermediate data generated can
    be saved to make displaying the result faster if run several times, say to
    test variations on the model."""
    import fullModel, path, models
    if (args.out is not None or args.view) and not args.files:
        print """Cannot set directory for intermediate files (-o/--out) unless also setting option
            to store these files (-f/--filles)."""
        return -1
    if not args.files: # no intermediate files
        return fullModel.processData(args.filename,[args.start,args.end],args.save,args.type,save=args.no)
    else: # intermediate files
        if args.out is None:
            if "." in args.filename[-4:]:
                args.out = args.filename[:-4]
            else:
                args.out = args.filename
        if path.processData(args.filename,[args.start,args.end],args.out,args.type, args.view):
            print "Could not generate point data files."
            return -1
        if args.no and args.save is None:
            if "." in args.filename[-4:]:
                args.save = args.filename[:-4]+".png"
            else:
                args.save = args.filename+".png"
        return models.compare(args.out,save=args.save)

def display(args):
    """Produces a radargram from existing intermediate data."""
    import models
    if args.no and args.save is None:
        print "set"
        args.save = args.directory + ".png"
    return models.compare(args.directory,args.adjusted,save=args.save)

if __name__ == "__main__":   
    parser = argparse.ArgumentParser(prog="python -m declutter")
    subparsers = parser.add_subparsers(description="""Select one of the commands listed below.
                        Use '<command> -h' to display help for that particular command.
                        For more detailed instructions, see the gitlab page
                        https://gitlab.data.bas.ac.uk/dboyle/radar-declutter/wikis/Contents""")

    loadParse = subparsers.add_parser('load', description="""
                        Takes a raster and generates the slope and aspect of the surface, saving all three
                        these along with the original raster in 'maps.hdf5'. The raster will be projected
                        if a coordinate system other than UTM is used and the raster can also be changed to
                        use a cell size which is a multiple of the original.""",
                        help='Loads a raster into numpy arrays.')
    loadParse.add_argument("filename",help="Name of the raster file to load.")
    loadParse.add_argument("latitude",type=float,help="Latitude of the region, used to determine UTM zone for projection.")
    loadParse.add_argument("longitude",type=float,help="Longitude of the region, used to determine UTM zone for projection.")
    loadParse.add_argument("-c","--cellSize",type=float,help="Allows the raster to be stored with a larger cell size than the original.")
    loadParse.add_argument("-o","--out",help="Specifies a directory to store the raster in if it is resampled or projected.")
    loadParse.set_defaults(func=load)

    modelParse = subparsers.add_parser('model', description="""
                        Assuming 'maps.hdf5' exists in the current directory, this takes a series of points and generates
                        an estimate of the radargram expected.""",
                        help='Displays the radargram for a path.')
    modelParse.add_argument("filename",help="Name of the path file to load.")
    modelParse.add_argument("-t","--type",choices=["gpx","dst","xyz"],help="""Specifies the file type, one of 'gpx', 'dst' or 'xyz'.
                            This is identified by the file extension if not given.""")
    modelParse.add_argument("--start",type = int,default=0, help = "How many points to ignore from the start of the path.")
    modelParse.add_argument("--end",type = int,default=0, help = "How many points to ignore from the end of the path.")
    modelParse.add_argument("-f","--files",action="store_true", help = "Store the intermediate data generated for each point.")
    modelParse.add_argument("-o","--out", help = "Specifies the directory to store the intermediate point data in.") #  only valid if -f selected
    exclusive = modelParse.add_mutually_exclusive_group()
    modelParse.add_argument("-v","--view",action="store_true", help = "Location of visible points stored too, helpful for analysis.") # only valid if -f
    exclusive.add_argument("-n","--no",action="store_false", help = "Do not save the radargram produced, but still displays it.")
    exclusive.add_argument("-s","--save", help = "The name to save the radargram as.")
    modelParse.set_defaults(func=model)

    displayParse = subparsers.add_parser('display',description="""
                        Assuming 'maps.hdf5' exists in the current directory, this takes a directory of point data and generates
                        an estimate of the radargram expected.""",
                        help='Displays a radargram from existing point data.')
    displayParse.add_argument('directory',help="Name of the directory containing the point data.")
    displayParse.add_argument('-a','--adjusted',action="store_true", help="Aligns points of equal elevation on the radargram rather than time.")
    exclusive = displayParse.add_mutually_exclusive_group()
    exclusive.add_argument("-n","--no",action="store_false", help = "Do not save the radargram produced, but still displays it.")
    exclusive.add_argument("-s","--save", help = "The name to save the radargram as.")
    displayParse.set_defaults(func=display)

    a = parser.parse_args()
    a.func(a)
    
