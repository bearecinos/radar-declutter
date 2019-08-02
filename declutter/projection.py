"""Enables projection of an arcMap raster into a different coordinate system
using arcpy."""
import arcpy
import numpy as np
import utm
import pyproj
import os
from errors import RasterError

def determineSystem(lat,lon):
    """Returns a string arcpy recognises as a spatial reference.
    This wil be for the UTM zone of the provided coordinates, except if
    they are near the poles, in which case the UPS system is used."""
    if lat > 83.5:
        return north
    if lat < -79.5:
        return south
    _,_,num,_ = utm.from_latlon(lat,lon)
    if lat > 0:
        let = "N"
    else:
        let = "S"
    return prefix+str(num)+let

north = 'Projected Coordinate Systems/Polar/UPS North'
south = 'Projected Coordinate Systems/Polar/UPS South'
prefix = 'Projected Coordinate Systems/UTM/WGS 1984/Northern Hemisphere/WGS 1984 UTM Zone '

# utm limits are -80 and +84, have 0.5 degree overlap with polar so convert if <-79.5 or >83.5
# (solves issue of being right on boundary with first gps point and crossing later)

# all mappings should be based on same gps path so should be consistent.
# still need utm to determine the zone
# possible only arcMap will know used coordinate system so need reference gps point

def project(source, systemName, saveAs = "projected"):
    """Converts the input raster to the given coordinate system. Note that this
    operation is attempted even if the input raster is already in that coordinate
    system.

    Parameters:
    source - string/Raster : The name/instance of the raster to project.
    systemName - string/SpatialReference : The name/instance of the spatial
        reference to project to. A string will be passed to arcpy to find
        the corresponding spatial reference object.
    saveAs - string (optional) : The name to save the projected raster as.
        This will be 'projected' by default.
        
    Returns
    The projected arcpy raster.
    Raises RasterError if the projection cannot be performed. This is usually
    because the wrong coordinate system has been provided.
    """
    if type(source) == str:
        source = arcpy.Raster(source)
    sr = systemName
    if type(systemName) == str:
        sr = arcpy.SpatialReference(systemName)
    arcpy.env.workspace = os.getcwd()
    arcpy.env.overwriteOutput = True
    try:
        arcpy.ProjectRaster_management(source,saveAs,sr,"CUBIC")
    except arcpy.ExecuteError as e:
        if "attempted on an empty geometry" in e.message:
            raise RasterError("arcpy error while projecting, check input latitude and longitude.")
    return arcpy.Raster(saveAs)

