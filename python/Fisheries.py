# Marine InVEST: Fisheries GIS Tool
# Author: Gregg Verutes
# 10/06/11

## DISSOLVE ZONES SO MULTIPLE FEATURES GET PUT INTO ONE ##
## DISSOLVING TAKES TOO LONG, CLIP FIRST? NEED TO DISSOLVE AT ALL? ##
## 'ZONESFIELD' NEEDS TO BE NUMRIC NOW.  WANT TO ALLOW FOR STRINGS TOO --- DICTIONARY? ##

# import modules
import sys, string, os, datetime, csv
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
gp.CheckOutExtension("spatial")

# error messages
msgArguments = "\nProblem with arguments."

# import modules
try:
    import numpy as np
except:
    gp.AddError(msgNumPyNo)
    raise Exception

try:
    # get parameters
    parameters = []
    now = datetime.datetime.now()
    parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
    gp.workspace = gp.GetParameterAsText(0)
    parameters.append("Workspace: "+ gp.workspace)
    InputZones = gp.GetParameterAsText(1)
    parameters.append("Zones: "+ InputZones)
    ZonesField = gp.GetParameterAsText(2)
    parameters.append("Zones Identifier Field: "+ ZonesField)
    FishingArea = gp.GetParameterAsText(3)
    parameters.append("Fishing Area: "+ FishingArea)    
    DEM = gp.GetParameterAsText(4)
    parameters.append("Digital Elevation Model (DEM): "+ DEM)
    Depth = gp.GetParameterAsText(5)
    if Depth:
        Depth = int(gp.GetParameterAsText(5))
    parameters.append("Depth Threshold for Fishing (meters): "+ str(Depth))
    Mangrove = gp.GetParameterAsText(6)
    parameters.append("Mangrove: "+ Mangrove)
    Coral = gp.GetParameterAsText(7)
    parameters.append("Coral: "+ Coral)
    Seagrass = gp.GetParameterAsText(8)
    parameters.append("Seagrass: "+ Seagrass)        
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
outputws = gp.workspace + os.sep + "Output" + os.sep
interws = gp.workspace + os.sep + "intermediate" + os.sep

DEM_rc = interws + "DEM_rc"
removeFA = interws + "removeFA.shp"
unionFC = interws + "unionFC.shp"
eraseFA = interws + "eraseFA.shp"

FishingArea_diss = interws + "FishingArea_diss.shp"
unionFishingAreaZones = interws + "unionFishingAreaZones.shp"
Mangrove_diss = interws + "Mangrove_diss.shp"
unionMangroveZones = interws + "unionMangroveZones.shp"
Coral_diss = interws + "Coral_diss.shp"
unionCoralZones = interws + "unionCoralZones.shp"
Seagrass_diss = interws + "Seagrass_diss.shp"
unionSeagrassZones = interws + "unionSeagrassZones.shp"

Zones = outputws + "Zones.shp"
FishingAreaZones = outputws + "FishingAreaZones.shp"
MangroveZones = outputws + "MangroveZones.shp"
CoralZones = outputws + "CoralZones.shp"
SeagrassZones = outputws + "SeagrassZones.shp"

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

def ckProjection(data):
    dataDesc = gp.describe(data)
    spatreflc = dataDesc.SpatialReference
    if spatreflc.Type <> 'Projected':
        gp.AddError(data +" does not appear to be projected.  It is assumed to be in meters.")
        raise Exception
    if spatreflc.LinearUnitName <> 'Meter':
        gp.AddError("This model assumes that "+data+" is projected in meters for area calculations.  You may get erroneous results.")
        raise Exception


###########################################
######## CHECK INPUTS & DATA PREP #########
###########################################

gp.AddMessage("\nChecking and preparing the inputs...")
ckProjection(InputZones)
ckProjection(FishingArea)
ckProjection(DEM)
gp.CopyFeatures_management(InputZones, Zones, "", "0", "0", "0")

gp.AddMessage("\nProcessing DEM and modifying fishing areas...")

gp.Reclassify_sa(DEM, "Value", "-1000000 "+str(-Depth)+" 1;0 100000000 2", DEM_rc, "NODATA")
gp.RasterToPolygon_conversion(DEM_rc, removeFA, "SIMPLIFY", "VALUE")

# add marker to eraseFC and unite the two FCs
fields = gp.ListFields(FishingArea, "ERASE")
field_found = fields.Next()
if field_found:
    gp.DeleteField_management(FishingArea, "ERASE")

removeFA = AddField(removeFA, "ERASE", "SHORT", "0", "0")
gp.CalculateField_management(removeFA, "ERASE", "1", "VB")
UnionExpr = FishingArea+" 1; "+removeFA+" 2"        
gp.Union_analysis(UnionExpr, unionFC)

# select features where "ERASE = 0"
gp.Select_analysis(unionFC, eraseFA, "\"ERASE\" = 0")

# area of input 'Zones'
Zones = AddField(Zones, "AREA", "FLOAT", "0", "0")
gp.CalculateField_management(Zones, "AREA", "!shape.area@squarekilometers!", "PYTHON", "")

# create 'AreaStatsArray' to store calcs
NumCols = 3
if Mangrove:
    NumCols += 1
if Coral:
    NumCols += 1
if Seagrass:
    NumCols += 1
AreaStatsList = np.zeros(gp.GetCount_management(Zones)*NumCols, dtype=np.float64)
AreaStatsArray = np.reshape(AreaStatsList, (gp.GetCount_management(Zones),NumCols))

cur = gp.UpdateCursor(Zones)
row = cur.Next()
count = 0
while row:
    AreaStatsArray[count][0] = row.GetValue(ZonesField) ## has to be float now ##
    AreaStatsArray[count][1] = row.GetValue("AREA") 
    count += 1
    row = cur.next()
del row, cur

gp.AddMessage("\nCalculating area of fishing and habitat within each zone...")  
gp.Extent = Zones

# add fields, dissolve, union, select, calculate area...
Zones = AddField(Zones, "FISH_AREA", "SHORT", "0", "0")
gp.CalculateField_management(Zones, "FISH_AREA", "1", "VB")
gp.Dissolve_management(eraseFA, FishingArea_diss, "")
FishingArea_diss = AddField(FishingArea_diss, "ZONES", "SHORT", "0", "0")
gp.CalculateField_management(FishingArea_diss, "ZONES", "1", "VB")
UnionFishingAreaExpr =  FishingArea_diss+" 1; "+Zones+" 2"        
gp.Union_analysis(UnionFishingAreaExpr, unionFishingAreaZones)
gp.Select_analysis(unionFishingAreaZones, FishingAreaZones, "\"ZONES\" = 1 AND \"FISH_AREA\" = 1")
FishingAreaZones = AddField(FishingAreaZones, "AREA", "FLOAT", "0", "0")
gp.CalculateField_management(FishingAreaZones, "AREA", "!shape.area@squarekilometers!", "PYTHON", "")
cur = gp.UpdateCursor(FishingAreaZones)
row = cur.Next()
count = 0
while row:
    AreaStatsArray[count][2] = row.GetValue("AREA") 
    count += 1
    row = cur.next()
del row, cur

if Mangrove:
    Zones = AddField(Zones, "MANGRV", "SHORT", "0", "0")
    gp.CalculateField_management(Zones, "MANGRV", "1", "VB")
    gp.Dissolve_management(Mangrove, Mangrove_diss, "")
    Mangrove_diss = AddField(Mangrove_diss, "ZONES", "SHORT", "0", "0")
    gp.CalculateField_management(Mangrove_diss, "ZONES", "1", "VB")
    UnionMangroveExpr =  Mangrove_diss+" 1; "+Zones+" 2"        
    gp.Union_analysis(UnionMangroveExpr, unionMangroveZones)
    gp.Select_analysis(unionMangroveZones, MangroveZones, "\"ZONES\" = 1 AND \"MANGRV\" = 1")
    MangroveZones = AddField(MangroveZones, "AREA", "FLOAT", "0", "0")
    gp.CalculateField_management(MangroveZones, "AREA", "!shape.area@squarekilometers!", "PYTHON", "")
    cur = gp.UpdateCursor(MangroveZones)
    row = cur.Next()
    count = 0
    while row:
        AreaStatsArray[count][3] = row.GetValue("AREA") 
        count += 1
        row = cur.next()
    del row, cur
    
if Coral:
    Zones = AddField(Zones, "CORAL", "SHORT", "0", "0")
    gp.CalculateField_management(Zones, "CORAL", "1", "VB")
    gp.Dissolve_management(Coral, Coral_diss, "")
    Coral_diss = AddField(Coral_diss, "ZONES", "SHORT", "0", "0")
    gp.CalculateField_management(Coral_diss, "ZONES", "1", "VB")
    UnionCoralExpr =  Coral_diss+" 1; "+Zones+" 2"        
    gp.Union_analysis(UnionCoralExpr, unionCoralZones)
    gp.Select_analysis(unionCoralZones, CoralZones, "\"ZONES\" = 1 AND \"CORAL\" = 1")
    CoralZones = AddField(CoralZones, "AREA", "FLOAT", "0", "0")
    gp.CalculateField_management(CoralZones, "AREA", "!shape.area@squarekilometers!", "PYTHON", "")
    cur = gp.UpdateCursor(CoralZones)
    row = cur.Next()
    count = 0
    while row:
        AreaStatsArray[count][4] = row.GetValue("AREA") 
        count += 1
        row = cur.next()
    del row, cur
    
if Seagrass:
    Zones = AddField(Zones, "SEAGRASS", "SHORT", "0", "0")
    gp.CalculateField_management(Zones, "SEAGRASS", "1", "VB")
    gp.Dissolve_management(Seagrass, Seagrass_diss, "")
    Seagrass_diss = AddField(Seagrass_diss, "ZONES", "SHORT", "0", "0")
    gp.CalculateField_management(Seagrass_diss, "ZONES", "1", "VB")
    UnionSeagrassExpr =  Seagrass_diss+" 1; "+Zones+" 2"        
    gp.Union_analysis(UnionSeagrassExpr, unionSeagrassZones)
    gp.Select_analysis(unionSeagrassZones, SeagrassZones, "\"ZONES\" = 1 AND \"SEAGRASS\" = 1")
    SeagrassAreaZones = AddField(SeagrassZones, "AREA", "FLOAT", "0", "0")
    gp.CalculateField_management(SeagrassZones, "AREA", "!shape.area@squarekilometers!", "PYTHON", "")
    cur = gp.UpdateCursor(SeagrassZones)
    row = cur.Next()
    count = 0
    while row:
        AreaStatsArray[count][5] = row.GetValue("AREA") 
        count += 1
        row = cur.next()
    del row, cur


gp.AddMessage("\nCreating output CSV...")  
# write to CSV
AreaStatsCSV  = open(outputws+"AreaStats.csv", "wb")
writer = csv.writer(AreaStatsCSV, delimiter=',', quoting=csv.QUOTE_NONE)
count = -1
while count < gp.GetCount_management(Zones):
    if count == -1:
        writer.writerow(['ZONE_ID', 'ZONE_AREA', 'FISHING_AREA', 'MANGROVE_AREA', 'CORAL_AREA', 'SEAGRASS_AREA'])
    else:
        writer.writerow(AreaStatsArray[count])
    count += 1
AreaStatsCSV.close()

# create parameter file
parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
parafile = open(gp.workspace+"\\Output\\parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
parafile.writelines("FISHERIES TIER 1 MODEL\n")
parafile.writelines("______________________\n\n")
for para in parameters:
    parafile.writelines(para+"\n")
    parafile.writelines("\n")
parafile.close()

##    # delete superfluous intermediate data
##    gp.workspace = interws
##    del1 = []
##    del2 = []
##    deletelist = del1 + del2
##    for data in deletelist:
##        if gp.exists(data):
##            gp.Delete_management(data)
##    del gp

##except Exception, ErrorDesc:
    ##gp.AddMessage(gp.GetMessages())