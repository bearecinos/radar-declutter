"""A custom error for if arcpy fails. Allows well known problems to be
identified despite often not helpful error codes/messages from arcpy."""
class Error(Exception):
    pass

class RasterError(Error):
    pass
