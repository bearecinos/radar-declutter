"""Provides a way of running the module by command line arguments, using the
three methods defined below. The -h option provides information on the whole
interface and also on each of the methods.
"""
import argparse

def load(args):
    """Converts a raster to numpy format in maps.hdf5, projecting and resampling
    if needed.
    Returns 0 if successful, else -1."""
    print "Loading..."
    import makeArrays
    return makeArrays.makeAll(args.filename,args.latitude,args.longitude,args.cellSize,args.out)

def model(args):
    """Produces a radargram from a gps path. The intermediate data generated can
    be saved to make displaying the result faster if run several times, say to
    test variations on the model.
    Returns 0 if successful, else -1."""
    import fullModel, path, models
    from modelling import parameters
    
    if (args.out is not None or args.view) and not args.files: # arguments don't make sense
        print """Must also set (-f/--files) to enable storing files."""
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
        if args.no and args.save is None: # set name to save radargram as
            if "." in args.filename[-4:]:
                args.save = args.filename[:-4]+".png"
            else:
                args.save = args.filename+".png"
        #models.loadParameters()
        parameters.loadParameters()
        return models.compare(args.out,save=args.save)

def display(args):
    """Produces a radargram from existing intermediate data.
    Returns 0 if successful, else -1."""
    import models
    from modelling import parameters
    #models.loadParameters()
    parameters.loadParameters() 
    if args.no and args.save is None:
        args.save = args.directory + ".png"
    return models.compare(args.directory,args.adjusted,save=args.save)

def setParams(args):
    """Overwrites any existing parameters stored which will revert to their default values if not fixed by arguments passed."""
    import numpy as np
    import os
    if args.show  and (args.maxdist is not None or args.maxtime is not None or
                                  args.steps is not None or args.dx is not None or args.dt is not None):
        print "Cannot show and set at same time. Do not use --show with other options."
        return -1
    elif args.show:
        d = np.load(os.path.dirname(__file__)+"/modelling/config.npy", allow_pickle=True).item()
        for key, val in d.items():
            if val is not None:
                print key + " = "+str(val)
        return 0
    if (args.steps is not None and (args.maxdist is not None or args.maxtime is not None) and
        (args.dx is not None or args.dt is not None)):
        if args.maxdist != args.dx * 1.5e8:
            print "Contradiction, cannot set range, granularity and total number of steps."
            return -1
    
    d = {"steps":args.steps, "maxDist":args.maxdist, "maxTime":args.maxtime, "dx":args.dx, "dt":args.dt}
    np.save(os.path.dirname(__file__)+"/modelling/config.npy", d)
    return 0

# Defines a parser for reading command line arguments
if __name__ == "__main__":
    # base command for any part of package
    parser = argparse.ArgumentParser(prog="python -m declutter")

    # will contain separate parsers for 3 different subcommands
    subparsers = parser.add_subparsers(description="""Select one of the commands listed below.
                        Use '<command> -h' to display help for that particular command.
                        For more detailed instructions, see the gitlab page
                        https://gitlab.data.bas.ac.uk/dboyle/radar-declutter/wikis/Contents""")

    # loading a raster into numpy format and storing in maps.hdf5
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

    # processing a path file to create a radargram
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

    # processing already created point data to create a radargram
    displayParse = subparsers.add_parser('display',description="""
                        Assuming 'maps.hdf5' exists in the current directory, this takes a directory of point data and generates
                        an estimate of the radargram expected.""",
                        help='Displays a radargram from existing point data.')
    displayParse.add_argument('directory',help="Name of the directory containing the point data.")
    displayParse.add_argument('-a','--adjusted',action="store_true", help="Aligns the surface below the radar with the top of the plot.")
    exclusive = displayParse.add_mutually_exclusive_group()
    exclusive.add_argument("-n","--no",action="store_false", help = "Do not save the radargram produced, but still displays it.")
    exclusive.add_argument("-s","--save", help = "The name to save the radargram as.")
    displayParse.set_defaults(func=display)

    paramParse = subparsers.add_parser('config',description="""
                        Sets/displays default options for the radargrams produced by command line calls to 'models'.""",
                        help='Sets the range/granularity of radargrams.')
    
    
    paramParse.add_argument('--show',action="store_true",help="Display any existing parameters set.")
    paramParse.add_argument('-n','--steps',type=int,help="Number of time samples in the radargram.")

    exclusive = paramParse.add_mutually_exclusive_group()
    exclusive.add_argument('--maxdist',type=float, help="Furthest surfaces to consider in the radargram.")
    exclusive.add_argument('--maxtime',type=float, help="Recording duration of the radargram.")

    exclusive = paramParse.add_mutually_exclusive_group()    
    exclusive.add_argument("--dx",type=float, help = "Spatial resolution of the radargram.")
    exclusive.add_argument("--dt",type=float, help = "Time resolution of the radargram.")
    paramParse.set_defaults(func=setParams)

    # parse sys.argv (command line input) and call appropriate method
    a = parser.parse_args()
    a.func(a)
    
