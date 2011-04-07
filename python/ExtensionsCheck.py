# Marine InVEST: Python Extensions Check
# Author: Gregg Verutes
# 02/16/11

# import basic extensions
import sys, string, os
import arcgisscripting
# create the geoprocessor object
gp = arcgisscripting.create()

# messages
msgPythonYes = "Python version "+str(sys.version_info[0])+"."+str(sys.version_info[1])+" is installed on your machine."
msgPythonNo = "- Python was NOT found.  Please install."
msgNumPyYes = "- NumPy has been installed.  However, make sure to install the latest version of NumPy."
msgNumPyNo = "- NumPy was NOT found.  Please install."
msgSciPyYes = "- SciPy has been successfully installed."
msgSciPyNo = "- SciPy was NOT found.  Please install."
msgWin32ComYes = "- PythonWin (allows Python to open Windows applications) has been successfully installed."
msgWin32ComNo = "- PythonWin (allows Python to open Windows applications) was NOT found.  Please install."
msgSpatialYes = "- Spatial Analyst is available."
msgSpatialNo = "- Spatial Analyst is NOT available.  Please check ESRI license agreement."

gp.AddMessage("\n")

try:
    gp.AddMessage(msgPythonYes)
except:
    gp.AddError(msgPythonNo)
    
gp.AddMessage("\n")
gp.AddMessage("Checking Python Libraries/Extensions...")

try:
    import numpy
    gp.AddMessage(msgNumPyYes)
except:
    gp.AddError(msgNumPyNo)
try:
    import scipy
    gp.AddMessage(msgSciPyYes) 
except:
    gp.AddError(msgSciPyNo)
try:
    import win32com.client
    gp.AddMessage(msgWin32ComYes) 
except:
    gp.AddError(msgWin32ComNo)

gp.AddMessage("\n")
gp.AddMessage("Checking ArcGIS Extensions...")

try:
    gp.CheckOutExtension("spatial")
    gp.AddMessage(msgSpatialYes)
except:
    gp.AddError(msgSpatialNo)

gp.AddMessage("\n")

del gp