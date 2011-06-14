# Marine InVEST: Python and GIS Extensions Check
# Author: Gregg Verutes
# 05/09/11

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
    # get parameters
    AestheticQuality = gp.GetParameterAsText(0)
    CP_Tier0 = gp.GetParameterAsText(1)
    FinfishAquaculture = gp.GetParameterAsText(2)
    HabitatRiskAssessment = gp.GetParameterAsText(3)
    OverlapAnalysis = gp.GetParameterAsText(4)
    WaveEnergy = gp.GetParameterAsText(5)
except:
    raise Exception, msgArguments + gp.GetMessages(2)

try:
    gp.AddMessage(msgPythonYes)
except:
    gp.AddError(msgPythonNo)
    
gp.AddMessage("\n")
gp.AddMessage("Checking Python extensions needed for the selected model(s)...")

try:
    if AestheticQuality == "true" or CP_Tier0 == "true" or FinfishAquaculture  == "true" or HabitatRiskAssessment  == "true" or OverlapAnalysis == "true" or WaveEnergy  == "true":
        import numpy
        gp.AddMessage(msgNumPyYes)
except:
    gp.AddError(msgNumPyNo)
    
try:
    if AestheticQuality == "true" or CP_Tier0 == "true" or WaveEnergy  == "true":
        import scipy
        gp.AddMessage(msgSciPyYes) 
except:
    gp.AddError(msgSciPyNo)
    
try:
    if CP_Tier0 == "true" or FinfishAquaculture  == "true" or HabitatRiskAssessment  == "true" or OverlapAnalysis == "true" or WaveEnergy  == "true":
        import win32com.client
        gp.AddMessage(msgWin32ComYes) 
except:
    gp.AddError(msgWin32ComNo)

try:
    if HabitatRiskAssessment  == "true":
        import matplotlib
        gp.AddMessage(msgMatplotlibYes) 
except:
    gp.AddError(msgMatplotlibNo)


gp.AddMessage("\n")
gp.AddMessage("Checking ArcGIS extensions needed for the selected model(s)...")

try:
    if AestheticQuality == "true" or CP_Tier0 == "true" or HabitatRiskAssessment  == "true" or OverlapAnalysis == "true" or WaveEnergy  == "true":
        gp.CheckOutExtension("spatial")
        gp.AddMessage(msgSpatialYes)
except:
    gp.AddError(msgSpatialNo)

gp.AddMessage("\n")

del gp