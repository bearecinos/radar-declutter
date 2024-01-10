"""Package for modelling the radar response of a surface along a path.
Warnings are disabled as getting False from testing less/greater than on a
NaN is the desired result.
lazy_import makes it seem like arcpy is available but is just wrapper
function which imports if needed. This avoids import errors if not
available and the slow setup time for arcpy if it is never needed."""

import warnings
# ignore warnings when comparing arrays containing NaNs.
# Always get False from NaN > x which is the desired behaviour.
# | np.isnan(x) is included where the result should be true.
warnings.filterwarnings("ignore")

from declutter.modelling import *
from declutter import changeSize
from declutter import radar
from declutter import path
from declutter import pointData
from declutter import makeArrays
from declutter import analysis
from declutter import compare
from declutter import fullModel

# modules to make visible with import *
__all__ = ["changeSize", "radar", "path", "pointData", "makeArrays",
           "modelling", "analysis", "compare", "fullModel"]

# runs whenever any part of package is run.
import arcpy
# try:
#     import lazy_import
#     arcpy = lazy_import.lazy_module("arcpy")
# except ImportError:
#     print("lazy-import module not available, attempting to import arcpy.")
#     try:
#         import arcpy
#     except ImportError:
#         print("arcpy not available, cannot use makeArrays.py.")
#         __all__.remove("makeArrays")
