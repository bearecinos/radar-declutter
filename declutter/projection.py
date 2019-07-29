import arcpy
import numpy as np
import utm
import pyproj
import os

def determineSystem(lat,lon):
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
    if type(source) == str:
        source = arcpy.Raster(source)
    sr = systemName
    if type(systemName) == str:
        sr = arcpy.SpatialReference(systemName)
    arcpy.env.workspace = os.getcwd()
    arcpy.env.overwriteOutput = True
    arcpy.ProjectRaster_management(source,saveAs,sr,"CUBIC")
    return arcpy.Raster(saveAs)

