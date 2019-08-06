"""Package for modelling the radar response of a surface along a path.
Warnings are disabled as getting False from testing less/greater than on a
NaN is the desired result.
lazy_import makes it seem like arcpy is available but is just wrapper
function which imports if needed. This avoids import errors if not
available and the slow setup time for arcpy if it is never needed."""

import warnings
# ignore warnings when comparing arrays containing NaNs.
# Always get False from NaN > x which is desired behaviour.
warnings.filterwarnings("ignore")

# modules which should be immediately visible
__all__ = ["changeSize","models","path","pointData","makeArrays"]  

# runs whenever any part of package is run.
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

  
