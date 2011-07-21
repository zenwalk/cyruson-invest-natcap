# Marine InVEST: Blue Carbon Model
# Authors: Gregg Verutes, Joey Bernhardt, Nasser Olwero
# 06/24/11

# import modules
import sys, string, os, datetime, shlex
import arcgisscripting
from math import *

# create the geoprocessor object
gp = arcgisscripting.create()
# set output handling
gp.OverwriteOutput = 1
# check out any necessary extensions
gp.CheckOutExtension("management")
gp.CheckOutExtension("analysis")
gp.CheckOutExtension("conversion")

# error messages
msgArguments = "Problem with arguments."

# import modules
try:
    import numpy as np
except:
    gp.AddError(msgNumPyNo)
    raise Exception
    
try:
    from win32com.client import Dispatch
except:
    gp.AddError(msgWin32ComNo)
    raise Exception

try:
try:
    # get parameters
    parameters = []
    now = datetime.datetime.now()
    parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
    gp.workspace = gp.GetParameterAsText(0)
    parameters.append("Workspace: "+ gp.workspace)
    GriddedSeascape = gp.GetParameterAsText(1)
    parameters.append("Gridded Seascape Output: "+ GriddedSeascape)
    Input_Directory = gp.GetParameterAsText(2)
    parameters.append("Input Layers Directory: "+ Input_Directory)
    ## xxx
except:
    raise Exception, msgArguments + gp.GetMessages(2)


try:
    thefolders=["intermediate","Output"]
    for folder in thefolders:
        if not gp.exists(gp.workspace+folder):
            gp.CreateFolder_management(gp.workspace, folder)
        
except:
    raise Exception, "Error creating folders"

# local variables 
outputws = gp.GetParameterAsText(0) + os.sep + "Output" + os.sep
interws = gp.GetParameterAsText(0) + os.sep + "intermediate" + os.sep

try:
    thefolders=["maps", "html_plots"]
    for folder in thefolders:
        if not gp.exists(outputws+folder):
            gp.CreateFolder_management(outputws, folder)
        
except:
    raise Exception, "Error creating folders"

# intermediate
GS_HQ = interws + "GS_HQ.shp"
GS_rst = interws + "GS_rst"

# output
##xxx = maps + "xxx"


def AddField(FileName, WghtFieldName, Type, Precision, Scale):
    fields = gp.ListFields(FileName, WghtFieldName)
    field_found = fields.Next()
    if field_found:
        gp.DeleteField_management(FileName, WghtFieldName)
    gp.AddField_management(FileName, WghtFieldName, Type, Precision, Scale, "", "", "NON_NULLABLE", "NON_REQUIRED", "")
    return FileName

def checkProjections(thedata):
    dataDesc = gp.describe(thedata)
    spatreflc = dataDesc.SpatialReference
    if spatreflc.Type <> 'Projected':
        gp.AddError(thedata+" does not appear to be projected.  It is assumed to be in meters.")
    if spatreflc.LinearUnitName <> 'Meter':
        gp.AddError("This model assumes that "+thedata+" is projected in meters for area calculations.  You may get erroneous results.")
        raise Exception

try:
    gp.AddMessage("\nChecking and preparing inputs...")
    # copy and populate input FC attribute table
    gp.CopyFeatures_management(GriddedSeascape, GS_HQ, "", "0", "0", "0")

    # grab cellsize info
    cur = gp.UpdateCursor(GS_HQ)
    row = cur.Next()
    cellsize = int(row.GetValue("CELL_SIZE"))
    del cur
    del row

    gp.FeatureToRaster_conversion(GS_HQ, "VALUE", GS_rst, str(cellsize))
except:
    gp.AddError(msgCheckInputs)
    raise Exception

try:
    # build vat
    sDesc = gp.describe(GS_rst)
    sInRasterDS = sDesc.catalogpath
    sExpression = "BuildVat " + sInRasterDS
    gp.MultiOutputMapAlgebra_sa(sExpression)
except:
    gp.AddError(msgBuildVAT)
    raise Exception		


# create parameter file
parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
parafile = open(outputws+"parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
parafile.writelines("BLUE CARBON MODEL PARAMETERS\n")
parafile.writelines("____________________________\n\n")
     
for para in parameters:
    parafile.writelines(para+"\n")
    parafile.writelines("\n")
parafile.close()

# delete superfluous intermediate data
##gp.workspace = interws
##del1 = []
##del2 = []
##deletelist = del1 + del2
##for data in deletelist:
##    if gp.exists(data):
##        gp.delete_management(data)
##del gp
