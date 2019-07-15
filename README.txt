Overview of each file/folder:
default.gdb - contains the heightmap, slope and aspect maps for the whole
	      area at 30m resolution, to be used by arcpy.
generation.mxd - arcmap profile used for examining maps.
maps.npz - same data but stored more densely in the numpy comressed format.
info.txt - structure information needed to interpret the numpy arrays.
other folders - 'info' used by arcpy. 'images' has a range of plots stored
		such as wiggle plots and comparisons of different methods
Generate.py - Call setPoint(...) then generateMaps() to create the
	      visibility/distance/incidence numpy arrays needed to simulate
	      the radar for a point. Relies on default.gdb
stateless.py - As above but faster as it doesn't rely on the arcmap rasters.
	     Instead takes its data from maps.npz and info.txt. Single call
	     as arguments passed rather than global to enable parallel solving.
viewshed.py - used by both to identify visible point, about 10x faster
	      than relying on arcpy to do so. The code is the main
	      bottleneck of the project, taking up about 90% of execution
	      time.
test.py - setPath(...) and setDetailed() prepare simple paths which 
          generateAll() then creates data for. 
path.py - similar but meant for processing actual data e.g. samplePath.txt
models.py - compare(...) displays the response expected for a given path
	    using a particular model or set of models. Can also do wiggle
	    plots and vary the wave modelled from the antenna.
edges.py - started looking into detecting edges in images to consider
 	   automating bed indentification. So far not useful.