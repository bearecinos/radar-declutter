"""Enables projection of an arcMap raster into a different coordinate system
using arcpy.
Although this projects to UPS coordinates where the point provided
is outside the range of UTM zones, this has not been tested.
In particular, elevation values may no longer be treated correctly.
"""
import arcpy
import numpy as np
import utm
import pyproj
import os
from declutter.errors import RasterError


def determineSystem(lat, lon):
    """Returns a string arcpy recognises as a spatial reference.
    This wil be for the UTM zone of the provided coordinates, or a
    UPS system if outside the bounds of UTM zones."""
    if lat > 83.5:
        return north
    if lat < -79.5:
        return south
    _, _, num, _ = utm.from_latlon(lat, lon)
    if lat > 0:
        let = "N"
    else:
        let = "S"
    return prefix+str(num)+let

north = 'Projected Coordinate Systems/Polar/UPS North'
south = 'Projected Coordinate Systems/Polar/UPS South'
prefix = ('Projected Coordinate Systems/UTM/WGS 1984/Northern Hemisphere/'
          'WGS 1984 UTM Zone ')


def project(source, systemName, saveAs="projected"):
    """Converts the input raster to the given coordinate system. Note that this
    operation is attempted even if the input raster is already in that
    coordinate system.

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
        arcpy.ProjectRaster_management(source, saveAs, sr, "CUBIC")
    except arcpy.ExecuteError as e:
        if "attempted on an empty geometry" in e.message:
            raise RasterError("arcpy error while projecting, check input "
                              "latitude and longitude.")
    return arcpy.Raster(saveAs)
