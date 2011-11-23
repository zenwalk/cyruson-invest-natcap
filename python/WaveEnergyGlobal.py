# Marine InVEST: Global Wave Energy Model
# Authors: CK Kim, Gregg Verutes, Mike Papenfus
# Coded for ArcGIS 9.3 and 10
# 10/12/11

# import modules
import sys, string, os, datetime
import shlex
from math import *
import arcgisscripting

# create the geoprocessor object
gp = arcgisscripting.create()
# set output handling
gp.OverwriteOutput = 1
# check out extensions
gp.CheckOutExtension("spatial")
gp.CheckOutExtension("management")
gp.CheckOutExtension("analysis")

# error messages
msgArguments = "\nProblem with arguments."
msgNoFeatures = "\nNo wave point features exist within the area of interest (AOI) specified."
msgMachineTables = "\nError reading in machine parameters."
msgWEcalc = "\nError during wave energy calculations."
msgValuation = "\nError during economic valuation."
msgClipWaveData = "\nError clipping wave data."
msgCalcDepth = "\nError calculating depth of wave points for wave power calculations.  Check that DEM covers the extent of the selected wave points."

# import modules
try:
    import numpy as np
except:
    gp.AddError(msgNumPyNo)
    raise Exception

try:
    from scipy.interpolate import interp2d
    from scipy import stats
except:
    gp.AddError(msgSciPyNo)
    raise Exception
    
try:
    from win32com.client import Dispatch
except:
    gp.AddError(msgWin32ComNo)
    raise Exception

try:
    # get parameters
    parameters = []
    now = datetime.datetime.now()
    parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
    gp.workspace = gp.GetParameterAsText(0)
    parameters.append("Workspace: "+ gp.workspace)
    gp.scratchWorkspace = gp.GetParameterAsText(0)
    parameters.append("Scratch Workspace: "+ gp.scratchWorkspace)
    WaveData = gp.GetParameterAsText(1)
    parameters.append("Projected Wave Data: "+ WaveData)
    EconParameter = gp.GetParameterAsText(2)
    parameters.append("Machine Economic Parameters Table: "+ EconParameter)
    NumUnits = gp.GetParameterAsText(3)
    parameters.append("Number of Machine Units: "+ str(NumUnits))
    if NumUnits:
        NumUnits = int(gp.GetParameterAsText(3))
except:
    raise Exception, msgArguments + gp.GetMessages(2)

##################################################################################################################

try:
    thefolders=["Output"]
    for folder in thefolders:
        if not gp.exists(gp.workspace+folder):
            gp.CreateFolder_management(gp.workspace, folder)
except:
    raise Exception, "Error creating folders"

outputws = gp.workspace + os.sep + "Output" + os.sep
WaveData_prj = outputws + "WaveData_"+str(NumUnits)+".shp"


##############################################
###### COMMON FUNCTION AND CHECK INPUTS ######
##############################################

def AddField(FileName, FieldName, Type, Precision, Scale):
    fields = gp.ListFields(FileName, FieldName)
    field_found = fields.Next()
    if field_found:
        gp.DeleteField_management(FileName, FieldName)
    gp.AddField_management(FileName, FieldName, Type, Precision, Scale, "", "", "NON_NULLABLE", "NON_REQUIRED", "")
    return FileName

g = 9.81 # meter / second square
d = 1028 # water density: kilogram / cubic meter
alfa = 0.86 # wave period parameter

##############################################################
####################### ECONOMIC ANALYSIS ####################
##############################################################

gp.AddMessage("\nPerforming economic valuation..." )
xlApp = Dispatch("Excel.Application")
xlApp.Visible = 0
xlApp.DisplayAlerts=0
xlBook1 = xlApp.Workbooks.Open(EconParameter[:-(1+len(EconParameter.split("\\")[-1]))])
econpath = EconParameter.split("\\")
econsheet = econpath[-1]
xlSheet1 = xlBook1.Worksheets(econsheet[:-1])

capRate = xlSheet1.Cells(2,2).Value
cc = xlSheet1.Cells(3,2).Value
cml = xlSheet1.Cells(4,2).Value
cul = xlSheet1.Cells(5,2).Value
col = xlSheet1.Cells(6,2).Value
omc = xlSheet1.Cells(7,2).Value
p = xlSheet1.Cells(8,2).Value
r = xlSheet1.Cells(9,2).Value
smlpm = xlSheet1.Cells(10,2).Value

year = 25.0
T = np.linspace(0.0, year-1.0, year)
rho = 1.0/(1.0+r)

# computes net present value of wave energy facility
def NPVWE(annualRevenue, annualCost):
    NPV = []
    for i in range(len(T)):
        NPV.append(rho**i * (annualRevenue[i] - annualCost[i]))
    return sum(NPV)        

# copy wave data file
gp.CopyFeatures_management(WaveData, WaveData_prj, "", "0", "0", "0")

# add various economic fields...
WaveData_prj = AddField(WaveData_prj, "UNITS", "SHORT", "0", "0")
WaveData_prj = AddField(WaveData_prj, "CAPWE_ALL", "DOUBLE", "8", "2")
WaveData_prj = AddField(WaveData_prj, "NPV_25Y", "DOUBLE", "8", "2")
gp.CalculateField_management(WaveData_prj, "UNITS", NumUnits, "VB")

# calculate net present value for each wave farm site            
cur = gp.UpdateCursor(WaveData_prj, "", "", "")
row = cur.Next()
while row:
    capturedWE = np.ones(len(T))*int(row.CAPWE_MWHY)*1000.0  # kWh here
    capturedWE[0] = 0.0

    # initial costs
    lenml = 3.0*abs(row.DEPTH_M)
    q = row.UNITS
    installCost = q*capRate*cc
    mooringCost = smlpm * lenml * cml * q
    transCost = (row.NEAR_DIST*cul/1000.0)
    IC = installCost + mooringCost + transCost

    # annual revenue and costs
    annualRevenue = p*q*capturedWE
    annualCost = omc*capturedWE*q
    annualCost[0] = IC
    
    row.SetValue("UNITS", q)
    row.SetValue("CAPWE_ALL", (row.CAPWE_MWHY * q))
    row.SetValue("NPV_25Y", (NPVWE(annualRevenue, annualCost))) # do NOT convert NPV to value "in thousands"
    cur.UpdateRow(row)
    row = cur.Next()
    
xlApp.ActiveWorkbook.Close(SaveChanges=0)
xlApp.Quit()
del xlApp
del cur, row
           