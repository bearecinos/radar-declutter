import arcpy
import numpy as np
import utm
import os

def determineSystem(lat,lon):
    _,_,num,_ = utm.from_latlon(lat,lon)
    if lat > 0:
        let = "N"
    else:
        let = "S"
    return prefix+str(num)+let

default = 'Projected Coordinate Systems/UTM/WGS 1984/Northern Hemisphere/WGS 1984 UTM Zone 33N'
prefix = 'Projected Coordinate Systems/UTM/WGS 1984/Northern Hemisphere/WGS 1984 UTM Zone '
def project(source, systemName = default, saveAs = "projected"):
    if type(source) == str:
        source = arcpy.Raster(source)
    sr = systemName
    if type(systemName) == str:
        sr = arcpy.SpatialReference(systemName)
    arcpy.env.workspace = os.getcwd()
    arcpy.env.overwriteOutput = True
    arcpy.ProjectRaster_management(source,saveAs,sr,"CUBIC")
    return arcpy.Raster(saveAs)
