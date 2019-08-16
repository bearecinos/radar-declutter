"""A custom error class which allows well known problems to be
identified despite often not helpful error codes/messages from arcpy."""


class DeclutterError(Exception):
    '''Parent for all package defined errors. Allows more to be defined in
    future and be caught be single "excpet DeclutterError" statement.'''
    pass


class RasterError(DeclutterError):
    '''Raised in serveral cases where arcpy fails to run.'''
    pass
