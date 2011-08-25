# Marine InVEST: Remove Habitat
# Author: Gregg Verutes
# 08/23/11

# import modules
import sys, string, os, datetime
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
msgArguments = "\nProblem with arguments."

# get parameters
parameters = []
RiskThreshList = []
HabRemoveList = []
now = datetime.datetime.now()
parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
gp.workspace = gp.GetParameterAsText(0)
parameters.append("Workspace: "+ gp.workspace)
HRA_Workspace = gp.GetParameterAsText(1)
parameters.append("Gridded Seascape Output: "+ HRA_Workspace)
Hab_Directory = gp.GetParameterAsText(2)
parameters.append("Habitat Data Directory: "+ Hab_Directory)
H1_RiskThresh = gp.GetParameterAsText(3)
parameters.append("Habitat #1 Risk Threshold: "+ H1_RiskThresh)
if float(H1_RiskThresh) > 0.0:
    RiskThreshList.append(float(gp.GetParameterAsText(3)))
    HabRemoveList.append(1)
H2_RiskThresh = gp.GetParameterAsText(4)
parameters.append("Habitat #2 Risk Threshold: "+ H2_RiskThresh)
if float(H2_RiskThresh) > 0.0:
    RiskThreshList.append(float(gp.GetParameterAsText(4)))
    HabRemoveList.append(2)  
H3_RiskThresh = gp.GetParameterAsText(5)
parameters.append("Habitat #3 Risk Threshold: "+ H3_RiskThresh)
if float(H3_RiskThresh) > 0.0:
    RiskThreshList.append(float(gp.GetParameterAsText(5)))
    HabRemoveList.append(3)   
H4_RiskThresh = gp.GetParameterAsText(6)
parameters.append("Habitat #4 Risk Threshold: "+ H4_RiskThresh)
if float(H4_RiskThresh) > 0.0:
    RiskThreshList.append(float(gp.GetParameterAsText(6)))
    HabRemoveList.append(4)  
H5_RiskThresh = gp.GetParameterAsText(7)
parameters.append("Habitat #5 Risk Threshold: "+ H5_RiskThresh)
if float(H5_RiskThresh) > 0.0:
    RiskThreshList.append(float(gp.GetParameterAsText(7)))
    HabRemoveList.append(5) 
H6_RiskThresh = gp.GetParameterAsText(8)
parameters.append("Habitat #6 Risk Threshold: "+ H6_RiskThresh)
if float(H6_RiskThresh) > 0.0:
    RiskThreshList.append(float(gp.GetParameterAsText(8)))
    HabRemoveList.append(6)  
H7_RiskThresh = gp.GetParameterAsText(9)
parameters.append("Habitat #7 Risk Threshold: "+ H7_RiskThresh)
if float(H7_RiskThresh) > 0.0:
    RiskThreshList.append(float(gp.GetParameterAsText(9)))
    HabRemoveList.append(7)
H8_RiskThresh = gp.GetParameterAsText(10)
parameters.append("Habitat #8 Risk Threshold: "+ H8_RiskThresh)
if float(H8_RiskThresh) > 0.0:
    RiskThreshList.append(float(gp.GetParameterAsText(10)))
    HabRemoveList.append(8)

try:
    thefolders=["intermediate","Output"]
    for folder in thefolders:
        if not gp.exists(gp.workspace+folder):
            gp.CreateFolder_management(gp.workspace, folder)
        
except:
    raise Exception, "Error creating folders"

# local variables
interws = gp.GetParameterAsText(0) + os.sep + "intermediate" + os.sep
outputws = gp.GetParameterAsText(0) + os.sep + "Output" + os.sep

# fetch calculator run variables
GS_HQ_area = HRA_Workspace + os.sep + "intermediate" + os.sep + "GS_HQ_area.shp"

# add "ERASE" field and populate with "0"
fields = gp.ListFields(GS_HQ_area, "ERASE")
field_found = fields.Next()
if field_found:
    gp.DeleteField_management(GS_HQ_area, "ERASE")

# variables
GS_H1_lyr = interws + "GS_H1_lyr.lyr"
GS_H2_lyr = interws + "GS_H2_lyr.lyr"
GS_H3_lyr = interws + "GS_H3_lyr.lyr"
GS_H4_lyr = interws + "GS_H4_lyr.lyr"
GS_H5_lyr = interws + "GS_H5_lyr.lyr"
GS_H6_lyr = interws + "GS_H6_lyr.lyr"
GS_H7_lyr = interws + "GS_H7_lyr.lyr"
GS_H8_lyr = interws + "GS_H8_lyr.lyr"

# various checks
def AddField(FileName, WghtFieldName, Type, Precision, Scale):
    fields = gp.ListFields(FileName, WghtFieldName)
    field_found = fields.Next()
    if field_found:
        gp.DeleteField_management(FileName, WghtFieldName)
    gp.AddField_management(FileName, WghtFieldName, Type, Precision, Scale, "", "", "NON_NULLABLE", "NON_REQUIRED", "")
    return FileName

def checkInteger(thedata):
    if thedata.find("0") == -1 and thedata.find("1") == -1 and thedata.find("2") == -1 and thedata.find("3") == -1 and thedata.find("4") == -1 and thedata.find("5") == -1 and thedata.find("6") == -1 and thedata.find("7") == -1 and thedata.find("8") == -1 and thedata.find("9") == -1:
        gp.AddError(thedata +" must contain an underscore followed by an integer ID at the end of it's name (e.g. filename_1.shp). This is necessary to properly link it with the input table.")
        raise Exception

def checkGeometry(thedata, Type, Message):
    if gp.Describe(thedata).ShapeType <> Type:
        raise Exception, "\nInvalid input: "+thedata+"\n"+Message+" must be of geometry type "+Type+"."

# compare "HabRemoveList" with "GS_HQ_area" to confirm that "CUMRISK_H" field exists
for i in range(0,len(HabRemoveList)):
    fields = gp.ListFields(GS_HQ_area, "CUMRISK_H"+str(HabRemoveList[i]))
    field_found = fields.Next()
    if not field_found:
        gp.AddWarning("Unable to remove habitat #"+str(int(HabRemoveList[i]))+" because there is no associated risk score from the Habitat Risk Assessment model run.  Please confirm this habitat layer exists.")
        raise Exception

#  prepare habitat input data 
gp.workspace = Hab_Directory
fcList = gp.ListFeatureClasses("*", "all")
fc = fcList.Next()
HabLyrList = []
HabIDList = []
HabCount = 0
while fc:
    # match SS ID with naming convention (_ID)
    checkInteger(fc)
    checkGeometry(fc, "Polygon", "Input feature class")
    HabLyrList.append(fc)
    fc2 = fc[::-1]
    j = fc2.find('_')
    indexS = len(fc)-j-1
    indexE = fc.find(".")
    fc_ID = fc[indexS+1:indexE]
    HabIDList.append(int(fc_ID))
    fc = fcList.Next()
    HabCount = HabCount + 1
del fc
HabZip = zip(HabIDList, HabLyrList)
HabZip.sort()
HabIDList, HabLyrList = zip(*HabZip)

gp.Extent = GS_HQ_area # set extent to "GS_HQ_area"

inputFieldList = []
fields = gp.ListFields(GS_HQ_area, "*")
fc_field = fields.Next()
while fc_field:
    inputFieldList.append(str(fc_field.name))
    fc_field = fields.Next()
del fc_field

# create expression so that feature layer only contains 'CUMRISK' fields
FieldShowExpr = inputFieldList[0]+" "+inputFieldList[0]+" HIDDEN NONE"
for i in range(1,len(inputFieldList)):
    if inputFieldList[i].find("CUMRISK") == 0:
        FieldShowExpr = FieldShowExpr+";"+inputFieldList[i]+" "+inputFieldList[i]+" VISIBLE NONE"
    else:
        FieldShowExpr = FieldShowExpr+";"+inputFieldList[i]+" "+inputFieldList[i]+" HIDDEN NONE"

NewHabLyrList = []
NewHabRemoveList = []
for i in range(0,len(HabRemoveList)):
    HabLyrListIndex = int(HabRemoveList[i])-1
    # make feature layers from "GS_HQ_area" for cells that exceed each habitat's risk threshold
    gp.MakeFeatureLayer_management(GS_HQ_area, "GS_H"+str(HabRemoveList[i])+"_lyr.lyr", "\"CUMRISK_H"+str(HabRemoveList[i])+"\" >= "+str(RiskThreshList[i]), interws, FieldShowExpr)
    # check that feature layer has at least one entry
    if gp.GetCount_management("GS_H"+str(HabRemoveList[i])+"_lyr.lyr") == 0:
        gp.AddMessage("No cells found that exceed specified risk threshold for habitat #"+str(HabRemoveList[i])+" ("+HabLyrList[HabLyrListIndex]+").")
    else:
        gp.AddMessage("Erasing habitat #"+str(HabRemoveList[i])+" ("+HabLyrList[HabLyrListIndex]+") in cells where cumulative risk score exceeds "+str(RiskThreshList[i]))
        gp.AddField_management("GS_H"+str(HabRemoveList[i])+"_lyr.lyr", "ERASE", "SHORT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        gp.CalculateField_management("GS_H"+str(HabRemoveList[i])+"_lyr.lyr", "ERASE", "1", "VB")
        # create new lists based on whether new feature layers have entries
        NewHabLyrList.append(HabLyrList[HabLyrListIndex])
        NewHabRemoveList.append(HabRemoveList[i])

# erase habitat (using "Union" and "Select" to avoid "Erase" tool)
for i in range(0,len(NewHabRemoveList)):
    HabVariable = Hab_Directory+"\\"+NewHabLyrList[i]
    UnionExpr = HabVariable+" 1; "+"GS_H"+str(NewHabRemoveList[i])+"_lyr.lyr"+" 2"        
    gp.Union_analysis(UnionExpr, interws+"GS_H"+str(NewHabRemoveList[i])+"_Intersect.shp")
    gp.Select_analysis(interws+"GS_H"+str(NewHabRemoveList[i])+"_Intersect.shp", outputws+NewHabLyrList[i], "\"ERASE\" = 0")