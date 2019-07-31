import argparse

def load(args):
    print "Loading..."
    import makeArrays
    makeArrays.makeAll(args.filename,args.latitude,args.longitude,args.cellSize,args.out)
    return 0

def model(args):
    import fullModel, path, models
    if args.out is not None and not args.files:
        print """Cannot set directory for intermediate files (-o/--out) unless also setting option
            to store these files (-f/--filles)."""
        return -1
    if not args.files: # no intermediate files
        fullModel.processData(args.filename,[args.start,args.end],args.save,args.type,save=args.no)
    else: # intermediate files
        if args.out is None:
            if "." in args.filename[-4:]:
                args.out = args.filename[:-4]
            else:
                args.out = args.filename
        path.processData(args.filename,[args.start,args.end],args.out,args.type)
        if args.no and args.save is None:
            if "." in args.filename[-4:]:
                args.save = args.filename[:-4]+".png"
            else:
                args.save = args.filename+".png"
        models.compare(args.out,save=args.save)
    return 0

def display(args):
    import models
    if args.no and args.save is None:
        print "set"
        args.save = args.directory + ".png"
    models.compare(args.directory,args.adjusted,save=args.save)
    return 0

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
    modelParse.add_argument("--start",type = int,default=0, help = "Can specify how many points to ignore from the start of the path.")
    modelParse.add_argument("--end",type = int,default=0, help = "Can specify how many points to ignore from the end of the path.")
    modelParse.add_argument("-f","--files",action="store_true", help = "Store the intermediate data generated for each point.")
    modelParse.add_argument("-o","--out", help = "Specifies the directory to store the intermediate point data in.") #  only valid if -f selected
    exclusive = modelParse.add_mutually_exclusive_group()
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
    
