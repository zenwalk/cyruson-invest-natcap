# Marine InVEST: Coastal Protection Tier 1 Python and GIS Extensions Check
# Author: Gregg Verutes
# 10/20/11

# import basic extensions
import sys, string, os
import arcgisscripting
# create the geoprocessor object
gp = arcgisscripting.create()

# error messages
msgArguments = "Problem with arguments."
msgPythonYes = "Python version "+str(sys.version_info[0])+"."+str(sys.version_info[1])+" is installed on your machine."
msgPythonNo = "- Python was NOT found.  Please install."
msgNumPyYes = "- NumPy extension has been installed.\n  However, make sure to install the latest version of NumPy."
msgNumPyNo = "- NumPy extension was NOT found.  Please install."
msgSciPyYes = "- SciPy extension has been successfully installed."
msgSciPyNo = "- SciPy extension was NOT found.  Please install."
msgWin32ComYes = "- PythonWin extension has been successfully installed."
msgWin32ComNo = "- PythonWin extension was NOT found.  Please install."
msgMatplotlibYes = "- Matplotlib extension has been successfully installed."
msgMatplotlibNo = "- Matplotlib extension was NOT found.  Please install."
msgSpatialYes = "- Spatial Analyst extension for ArcGIS is available.\n  However, make sure this extension is activated.\n  In ArcMap, select 'Customize' >> 'Extensions' >> check box next to 'Spatial Analyst'"
msgSpatialNo = "- Spatial Analyst extension for ArcGIS is NOT available.  Please check your ESRI license agreement."

try:
    gp.AddMessage(msgPythonYes)
except:
    gp.AddError(msgPythonNo)
    
gp.AddMessage("\n")
gp.AddMessage("Checking Python extensions needed for the selected model(s)...")

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

try:
	import matplotlib
	gp.AddMessage(msgMatplotlibYes) 
except:
    gp.AddError(msgMatplotlibNo)


gp.AddMessage("\n")
gp.AddMessage("Checking ArcGIS extensions needed for the selected model(s)...")

try:
	gp.CheckOutExtension("spatial")
	gp.AddMessage(msgSpatialYes)
except:
    gp.AddError(msgSpatialNo)

gp.AddMessage("\n")

del gp