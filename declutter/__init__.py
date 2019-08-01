"""Package for modelling the radar response of a surface along a path.
__init__ just declares which modules visible with import * .
lazy_import makes it seem like arcpy is available but really just wrapper
function which is imported if needed. Hence avoid import errors if not
available and don't have to wait for import if never used."""

import warnings
# ignore warnings when comparing arrays containing NaNs.
# Always get False from ar > x where NaNs present
# which is desired behaviour.
warnings.filterwarnings("ignore")

__all__ = ["changeSize","models","path","pointData","makeArrays"]  

try:
    import lazy_import
    arcpy = lazy_import.lazy_module("arcpy")
except ImportError:
    print "lazy-import module not available, attempting to import arcpy."
    try:
        import arcpy
    except ImportError:
        print "arcpy not available, cannot use makeArrays.py."
        __all__.remove("makeArrays")

  
