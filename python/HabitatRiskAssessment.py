# Marine InVEST: Habitat Risk Assessment Model
# Authors: Joey Bernhardt, Katie Arkema, Gregg Verutes, Jeremy Davies, Martin Lacayo
# 12/01/11

# import modules
import sys, string, os, datetime, shlex, csv
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
msgCheckInputs = "\nError checking and preparing inputs."
msgBuildVAT = "\nError building VAT for zonal statistics raster.  Make sure there are no spaces in your path names.  Try opening a new ArcMap session and re-run the HRA model."
msgPrepHSLayers ="\nError preparing habitat and stressor input layers."
msgCheckHSLayers = "\nError checking habitat and stressor input layers."
msgBuffRastHULayers = "\nError buffering and rasterizing stressor and habitat input layers."
msgOverlap = "\nError calculating spatial overlap."
msgGetHabStressRatings = "\nError obtaining habitat and stressor ratings from table."
msgCalcRiskMap = "\nError calculating risk scores within analysis grid (GS)."
msgPlotArrays = "\nError preparing exposure and consequence matrix values for plotting."
msgMapOutputs = "\nError generating map outputs."
msgPlotHTMLOutputs = "\nError generating plot and HTML outputs."
msgNumPyNo = "NumPy extension is required to run the Habitat Risk Assessment Model.  Please consult the Marine InVEST FAQ document for instructions on how to install."
msgWin32ComNo = "PythonWin extension is required to run the Habitat Risk Assessment Model.  Please consult the Marine InVEST FAQ document for instructions on how to install."
msgSciPyNo = "SciPy extension is required to run the Habitat Risk Assessment Model.  Please consult the Marine InVEST FAQ for instructions on how to install."
msgMatplotlibNo = "Matplotlib extension (version 1.0 or newer) is required to run the Habitat Risk Assessment Model.  Please consult the Marine InVEST FAQ document for instructions on how to install."

# import modules
try:
    import numpy as np
except:
    gp.AddError(msgNumPyNo)
    raise Exception

try:
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
    try:
        # get parameters
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
        gp.workspace = gp.GetParameterAsText(0)
        parameters.append("Workspace: "+ gp.workspace)
        GriddedSeascape = gp.GetParameterAsText(1)
        parameters.append("Gridded Seascape Output: "+ GriddedSeascape)
        Hab_Directory = gp.GetParameterAsText(2)
        parameters.append("Habitat Data Directory: "+ Hab_Directory)
        Stress_Directory = gp.GetParameterAsText(3)
        parameters.append("Stressor Data Directory: "+ Stress_Directory)
        HabStressRate_Table = gp.GetParameterAsText(4)
        parameters.append("Stressor Data Directory: "+ HabStressRate_Table)
        PlotBoolean = gp.GetParameterAsText(5)
        parameters.append("Create HTML output with risk plots: "+ PlotBoolean)
        RiskBoolean = gp.GetParameterAsText(6)
        parameters.append("Generate Habitat Risk Maps from Inputs: "+ RiskBoolean)
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
    interws = gp.GetParameterAsText(0) + os.sep + "intermediate" + os.sep
    outputws = gp.GetParameterAsText(0) + os.sep + "Output" + os.sep
    maps = outputws + os.sep + "maps" + os.sep
    html_plots = outputws + os.sep + "html_plots" + os.sep

    try:
        if PlotBoolean == "true":
            thefolders=["maps", "html_plots"]
        else:    
            thefolders=["maps"]
        for folder in thefolders:
            if not gp.exists(outputws+folder):
                gp.CreateFolder_management(outputws, folder)
    except:
        raise Exception, "Error creating folders"

    # intermediate
    GS_HQ = interws + "GS_HQ.shp"
    GS_HQ_lyr = interws + "GS_HQ_lyr.lyr"
    GS_rst = interws + "gs_rst"
    GS_HQ_risk = interws + "GS_HQ_risk.shp"
    GS_HQ_predom = interws + "GS_HQ_predom.shp"
    GS_HQ_area = interws + "GS_HQ_area.shp"
    GS_HQ_intersect = interws + "GS_HQ_intersect.shp"

    # output
    ecosys_risk = maps + "ecosys_risk"
    recov_potent = maps + "recov_potent"

    risk_plots = html_plots + "plots_risk.png"
    ecosysRisk_plots = html_plots + "plot_ecosys_risk.png"
    outputHTML = html_plots + "output.html"


    def AddField(FileName, WghtFieldName, Type, Precision, Scale):
        fields = gp.ListFields(FileName, WghtFieldName)
        field_found = fields.Next()
        if field_found:
            gp.DeleteField_management(FileName, WghtFieldName)
        gp.AddField_management(FileName, WghtFieldName, Type, Precision, Scale, "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        return FileName

    def checkPathLength(path):
        pathLength = len(path)
        if pathLength > 52:
            gp.AddError("In order to properly build a value attribute table (VAT) for the gridded seascape, the path length of your workspace must be shorten by at least "+str(pathLength-52)+" character(s).")
            raise Exception
    
    def checkProjections(thedata):
        dataDesc = gp.describe(thedata)
        spatreflc = dataDesc.SpatialReference
        if spatreflc.Type <> 'Projected':
            gp.AddError(thedata+" does not appear to be projected.  It is assumed to be in meters.")
        if spatreflc.LinearUnitName <> 'Meter':
            gp.AddError("This model assumes that "+thedata+" is projected in meters for area calculations.  You may get erroneous results.")
            raise Exception

    def checkInteger(thedata):
        if thedata.find("0") == -1 and thedata.find("1") == -1 and thedata.find("2") == -1 and thedata.find("3") == -1 and thedata.find("4") == -1 and thedata.find("5") == -1 and thedata.find("6") == -1 and thedata.find("7") == -1 and thedata.find("8") == -1 and thedata.find("9") == -1:
            gp.AddError(thedata +" must contain an underscore followed by an integer ID at the end of it's name (e.g. filename_1.shp). This is necessary to properly link it with the input table.")
            raise Exception

    # percentiles list (33%, 66%)
    def getPercentiles(list):
        PctList = []
        PctList.append(stats.scoreatpercentile(list, 1.0/3.0))
        PctList.append(stats.scoreatpercentile(list, 2.0/3.0))
        return PctList

    try:
        gp.AddMessage("\nChecking and preparing inputs...")
        # check to make sure workspace path length doesn't exceed limit for building VAT
        checkPathLength(gp.GetParameterAsText(0))
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


    ####################################################            
    ####### PREPARE HABITAT AND STRESSOR LAYERS ########
    ####################################################
    try:
        gp.workspace = Hab_Directory
        fcList = gp.ListFeatureClasses("*", "all")
        fc = fcList.Next()
        HabLyrList = []
        HabIDList = []
        HabCount = 0
        while fc:
            # match SS ID with naming convention (_ID)
            checkInteger(fc)
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


        gp.workspace = Stress_Directory
        fcList = gp.ListFeatureClasses("*", "all")
        fc = fcList.Next()
        StressLyrList = []
        StressLyrBuffList = []
        StressIDList = []
        StressCount = 0
        while fc:
            # match SS ID with naming convention (_ID)
            checkInteger(fc)
            StressLyrList.append(fc)
            fc0 = fc.replace(".", "")
            fc1 = fc0.replace("_", "")
            StressLyrBuffList.append(fc1[:8]+"_buff.shp")
            fc2 = fc[::-1]
            j = fc2.find('_')
            indexS = len(fc)-j-1
            indexE = fc.find(".")
            fc_ID = fc[indexS+1:indexE]
            StressIDList.append(int(fc_ID))
            fc = fcList.Next()
            StressCount = StressCount + 1
        del fc

        StressZip = zip(StressIDList, StressLyrList, StressLyrBuffList)
        StressZip.sort()
        StressIDList, StressLyrList, StressLyrBuffList = zip(*StressZip)
    except:
        gp.AddError(msgPrepHSLayers)
        raise Exception


    ####################################################            
    ###### CK CONSISTENCY WITH HAB/STRESS LAYERS #######
    ####################################################
    try:

        # get data from CSV
        rawList = []
        csvReader = csv.reader(open(HabStressRate_Table, 'rb'), delimiter=',', quotechar='|');
        for row in csvReader:
            rawList.append(row)

        # get habitat and stressor counts
        HabInputCount = int(rawList[14][0])
        StressInputCount = int(rawList[21][0])

        # put values in temporary lists
        HabVar = []
        StressVar = []
        HabStressVar = []
        CritWeights = []
        for i in range(15,len(rawList)):
            if i > 15 and i < 15 + HabCount+1:
                HabVar.append(rawList[i][2:])
            elif i > 15 + HabInputCount+3 and i < 15 + HabInputCount+3 + StressInputCount+1:
                StressVar.append(rawList[i][2:])
            elif i > 15 + HabInputCount+3 + StressInputCount+3 and i < 15 + HabInputCount+3 + StressInputCount+3 + (HabInputCount*StressInputCount)+1:
                HabStressVar.append(rawList[i][4:])
            elif i > 15 + HabInputCount+3 + StressInputCount+3 + (HabInputCount*StressInputCount)+2 and i < 15 + HabInputCount+3 + StressInputCount+3 + (HabInputCount*StressInputCount)+12:
                CritWeights.append(int(rawList[i][0]))

        # adjust inter-criteria weights
        CritWeights = [0.5 if i == 0 else i for i in CritWeights]; CritWeights = [1.5 if i == 2 else i for i in CritWeights]
        ExpCritWeights = []
        ExpIndex = [1,2,0,3]
        for i in range(0,4)
            ExpCritWeights.append(CritWeights[ExpIndex[i]])
        ConsCritWeights = []
        ConsIndex = [7,8,9,10,4,5,6]            
        for i in range(0,7)
            ConsCritWeights.append(CritWeights[ConsIndex[i]])

        # populate buffer distance list
        StressBuffDistList = []
        for i in range(0,StressInputCount):
            StressBuffDistList.append(float(StressVar[i][5]))

        if HabInputCount <> HabCount:
            gp.AddError("There is an inconsistency between the number of habitat layers in the specified directory and the input spreadsheet.")
            raise Exception
        if StressInputCount <> StressCount:
            gp.AddError("There is an inconsistency between the number of stressor layers in the specified directory and the input spreadsheet.")
            raise Exception
    except:
        gp.AddError(msgCheckHSLayers)
        raise Exception


    ########################################################            
    ###### RASTERIZE AND BUFFER HAB AND STRESS LAYERS ######
    ########################################################
    try:
        gp.workspace = interws
        gp.Extent = GS_rst
        gp.MakeFeatureLayer_management(GS_HQ, GS_HQ_lyr, "", "", "")
        HabNoDataList = []
        del_hab = []
        for i in range(0,len(HabLyrList)):
            HabVariable = Hab_Directory+"\\"+HabLyrList[i]
            checkProjections(HabVariable)
            HabVariable = AddField(HabVariable, "VID", "SHORT", "0", "0")        
            gp.CalculateField_management(HabVariable, "VID", 1, "VB")
            gp.MakeFeatureLayer_management(HabVariable, HabLyrList[i][:-4]+".lyr", "", interws, "")
            SelectHab = gp.SelectLayerByLocation_management(HabLyrList[i][:-4]+".lyr", "INTERSECT", GS_HQ_lyr, "", "NEW_SELECTION")
            if gp.GetCount_management(SelectHab) == 0:
                HabNoDataList.append("yes")
            else:
                gp.FeatureToRaster_conversion(HabVariable, "VID", "hab_"+str(i+1), "50")
                HabNoDataList.append("no")
                del_hab.append("hab_"+str(i+1))
            gp.SelectLayerByAttribute_management(HabLyrList[i][:-4]+".lyr", "CLEAR_SELECTION", "")

        StressNoDataList = []
        del_stress = []
        for i in range(0,len(StressLyrList)):
            StressVariable = Stress_Directory+"\\"+StressLyrList[i]
            checkProjections(StressVariable)
            StressVariable = AddField(StressVariable, "VID", "SHORT", "0", "0")        
            gp.CalculateField_management(StressVariable, "VID", 1, "VB")
            if StressBuffDistList[i] > 0:
                gp.Buffer_analysis(StressVariable, StressLyrBuffList[i][:-4]+"_s"+str(i+1)+StressLyrBuffList[i][-4:], str(StressBuffDistList[i]) + " Meters", "FULL", "ROUND", "NONE", "")
                gp.MakeFeatureLayer_management(StressLyrBuffList[i][:-4]+"_s"+str(i+1)+StressLyrBuffList[i][-4:], StressLyrList[i][:-4]+".lyr", "", interws, "")
            else:
                gp.MakeFeatureLayer_management(StressVariable, StressLyrList[i][:-4]+".lyr", "", interws, "")
                
            SelectStress = gp.SelectLayerByLocation_management(StressLyrList[i][:-4]+".lyr", "INTERSECT", GS_HQ_lyr, "", "NEW_SELECTION")
            if gp.GetCount_management(SelectStress) == 0:                
                StressNoDataList.append("yes")
            else:
                StressNoDataList.append("no")
                del_stress.append("stress_"+str(i+1))
                if StressBuffDistList[i] > 0:
                    gp.FeatureToRaster_conversion(StressLyrBuffList[i][:-4]+"_s"+str(i+1)+StressLyrBuffList[i][-4:], "VID", "stress_"+str(i+1), "50")
                    gp.CopyFeatures_management(StressLyrBuffList[i][:-4]+"_s"+str(i+1)+StressLyrBuffList[i][-4:], maps+"s"+str(i+1)+"_"+StressLyrList[i][:-5]+"buff.shp", "", "0", "0", "0")
                else:
                    gp.FeatureToRaster_conversion(StressVariable, "VID", "stress_"+str(i+1), "50")
                    gp.CopyFeatures_management(StressVariable, maps+"s"+str(i+1)+"_"+StressLyrList[i][:-6]+".shp", "", "0", "0", "0")
            gp.SelectLayerByAttribute_management(StressLyrList[i][:-4]+".lyr", "CLEAR_SELECTION", "")
    except:
        gp.AddError(msgBuffRastHULayers)
        raise Exception        


    ############################################# 
    ############ CALCULATE OVERLAP ##############
    #############################################
    try:
        gp.AddMessage("\nCalculating spatial overlap...")

        def difference(a, b): # show whats in list b which isn't in list a
            return list(set(b).difference(set(a)))

        # determine which hab and stress rasters weren't in GS AOI        
        potHabList = range(1,HabCount+1)
        potStressList = range(1,StressCount+1)        
        rasterHabList = []
        rasterStressList = []
        gp.workspace = interws
        rasters = gp.ListRasters("hab_*", "GRID")
        rasters.Reset()
        for i in range(0,HabCount):
            raster = rasters.Next()
            if raster != None:
                rasterHabList.append(raster[4:])
        del rasters
        rasters = gp.ListRasters("stress_*", "GRID")
        rasters.Reset()
        for i in range(0,StressCount):
            raster = rasters.Next()
            if raster != None:
                rasterStressList.append(raster[7:])                  
        del rasters
        rasterHabList = [int(s) for s in rasterHabList]        
        rasterStressList = [int(s) for s in rasterStressList]
        diffHabList = difference(rasterHabList, potHabList)
        diffStressList = difference(rasterStressList, potStressList)

        # combine hab and stress rasters that overlap
        OverlapList = []
        OverlapNoDataList = []
        for i in range(0,len(HabLyrList)):
            for j in range(0,len(StressLyrList)):
                if i+1 not in diffHabList and j+1 not in diffStressList:
                    CmbExpr = "hab_"+str(i+1)+";"+"stress_"+str(j+1)
                    gp.Combine_sa(CmbExpr, "H"+str(i+1)+"S"+str(j+1))
                    OverlapList.append("H"+str(i+1)+"S"+str(j+1))
                    gp.BuildRasterAttributeTable_management("H"+str(i+1)+"S"+str(j+1), "Overwrite")                   
                    if gp.GetCount("H"+str(i+1)+"S"+str(j+1)) == 0:
                        OverlapNoDataList.append("yes")
                    else:
                        cur = gp.UpdateCursor("H"+str(i+1)+"S"+str(j+1))
                        row = cur.Next()
                        rstValue = row.GetValue("VALUE")
                        del cur, row
                        if rstValue < 0:
                            OverlapNoDataList.append("yes")
                        else:
                            OverlapNoDataList.append("no")
                else:
                    OverlapList.append("H"+str(i+1)+"S"+str(j+1))
                    OverlapNoDataList.append("yes")

        for i in range(0,len(HabLyrList)):
            GS_HQ = AddField(GS_HQ, "H"+str(i+1)+"_A", "DOUBLE", "8", "2")
        for i in range(0,len(OverlapList)):
            GS_HQ = AddField(GS_HQ, OverlapList[i]+"_A", "DOUBLE", "8", "2")
            GS_HQ = AddField(GS_HQ, OverlapList[i]+"_PCT", "DOUBLE", "8", "2")

        gp.snapRaster = GS_rst
        gp.cellsize = "MINOF"
      
        for i in range(0,len(HabLyrList)):
            if HabNoDataList[i] == "no":
                gp.ZonalStatisticsAsTable_sa(GS_rst, "VALUE", "hab_"+str(i+1), "zs_H"+str(i+1)+".dbf", "DATA")
                gp.AddJoin_management(GS_HQ_lyr, "VALUE", "zs_H"+str(i+1)+".dbf", "VALUE", "KEEP_COMMON")
                gp.CalculateField_management(GS_HQ_lyr, "GS_HQ.H"+str(i+1)+"_A", "[zs_H"+str(i+1)+".AREA]", "VB", "")
                gp.RemoveJoin_management(GS_HQ_lyr, "zs_H"+str(i+1))

        for i in range(0,len(OverlapList)):
            if OverlapNoDataList[i] == "no":
                gp.ZonalStatisticsAsTable_sa(GS_rst, "VALUE", OverlapList[i], "zs_"+OverlapList[i]+".dbf", "DATA")
                gp.AddJoin_management(GS_HQ_lyr, "VALUE", "zs_"+OverlapList[i]+".dbf", "VALUE", "KEEP_COMMON")
                gp.CalculateField_management(GS_HQ_lyr, "GS_HQ."+OverlapList[i]+"_A", "[zs_"+OverlapList[i]+".AREA]", "VB", "")
                gp.RemoveJoin_management(GS_HQ_lyr, "zs_"+OverlapList[i])

        gp.FeatureClassToFeatureClass_conversion(GS_HQ_lyr, gp.workspace, "GS_HQ_area.shp", "")
        GS_HQ_area = AddField(GS_HQ_area, "PREDOM_HAB", "SHORT", "", "")

        cur = gp.UpdateCursor(GS_HQ_area)          
        row = cur.Next()
        while row:
            # calculate percent overlap for all permutations
            for i in range(0,len(HabLyrList)):
                for j in range(0,len(StressLyrList)):
                    HabArea = row.GetValue("H"+str(i+1)+"_A")
                    if HabArea <> 0:
                        HabStrOverlap = row.GetValue("H"+str(i+1)+"S"+str(j+1)+"_A")
                        PctOverlap = float(HabStrOverlap/HabArea)*100.0
                        row.SetValue("H"+str(i+1)+"S"+str(j+1)+"_PCT", PctOverlap)
                    else:
                        row.SetValue("H"+str(i+1)+"S"+str(j+1)+"_PCT", 0)
                        
            # determine predominate habitat
            PredomIDList = []
            for i in range(0,len(HabLyrList)):
                PredomIDList.append(row.GetValue("H"+str(i+1)+"_A"))

            if np.sum(PredomIDList) == 0:
                row.PREDOM_HAB = 0
            else:
                PredomIDMax = max(PredomIDList)
                row.PREDOM_HAB = (PredomIDList.index(PredomIDMax))+1
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur               
    except:
        gp.AddError(msgOverlap)
        raise Exception

    # delete some fields to avoid 250 max
    DelExpr = ""
    for i in range(0,len(HabLyrList)):
        for j in range(0,len(StressLyrList)):
            DelExpr = DelExpr + "H"+str(i+1)+"S"+str(j+1)+"_A;"
    DelExpr = DelExpr[:-1]
    gp.workspace = interws
    gp.DeleteField_management(GS_HQ_area, DelExpr)

    ##########################################################   
    ############ GRAB RATINGS FROM EXCEL TABLE  ##############
    ##########################################################
    try:
        gp.AddMessage("\nObtaining ratings for risk scoring and plotting...")

        # determine total number of permutations for habitat and stressor layers
        TotalHSCombo = int(HabCount*StressCount)
        ExposureList = np.zeros(TotalHSCombo*5, dtype=np.float64)
        ExposureArray = np.reshape(ExposureList, (TotalHSCombo,5))
        ExpQualityList = np.zeros(TotalHSCombo*5, dtype=np.float64)
        ExpQualityArray = np.reshape(ExpQualityList, (TotalHSCombo,5)) # 1-3 range
        ExpQualityList2 = np.zeros(TotalHSCombo*5, dtype=np.float64)
        ExpQualityArray2 = np.reshape(ExpQualityList2, (TotalHSCombo,5)) # 1-3 range

        ConsequenceList = np.zeros(TotalHSCombo*8, dtype=np.float64)
        ConsequenceArray = np.reshape(ConsequenceList, (TotalHSCombo,8))
        ConsQualityList = np.zeros(TotalHSCombo*8, dtype=np.float64)
        ConsQualityArray = np.reshape(ConsQualityList, (TotalHSCombo,8)) # 1-3 range
        ConsQualityList2 = np.zeros(TotalHSCombo*8, dtype=np.float64)
        ConsQualityArray2 = np.reshape(ConsQualityList2, (TotalHSCombo,8)) # 1-3 range

        rowCount = 0
        for i in range(0,HabCount):
            for j in range(0,StressCount):
                HabDataQual = float(HabVar[i][0]) # DQ of spatial overlap for habitat
                StressDataQual = float(StressVar[j][0]) # DQ of spatial overlap for stressor
                HabDataQual2 = float(HabVar[i][0])+1 # DQ of spatial overlap for habitat
                StressDataQual2 = float(StressVar[j][0])+1 # DQ of spatial overlap for stressor

                # ConsequenceArray: index 7 is blank for now (7 = weighted average (Exp) or sum (DQ))                 
                ConsequenceArray[rowCount][0] = int(HabVar[i][1]) # natural mortality rate
                ConsQualityArray[rowCount][0] = int(HabVar[i][2]) # DQ of natural mortality rate
                ConsQualityArray2[rowCount][0] = int(HabVar[i][2])+1 # DQ of natural mortality rate
                ConsequenceArray[rowCount][1] = int(HabVar[i][3]) # recruitment pattern
                ConsQualityArray[rowCount][1] = int(HabVar[i][4]) # DQ of recruitment pattern
                ConsQualityArray2[rowCount][1] = int(HabVar[i][4])+1 # DQ of recruitment pattern
                ConsequenceArray[rowCount][2] = int(HabVar[i][5]) # connectivity
                ConsQualityArray[rowCount][2] = int(HabVar[i][6]) # DQ of connectivity
                ConsQualityArray2[rowCount][2] = int(HabVar[i][6])+1 # DQ of connectivity
                ConsequenceArray[rowCount][3] = int(HabVar[i][7]) # age at maturity or recovery time
                ConsQualityArray[rowCount][3] = int(HabVar[i][8]) # DQ of age at maturity or recovery time
                ConsQualityArray2[rowCount][3] = int(HabVar[i][8])+1 # DQ of age at maturity or recovery time
                ConsequenceArray[rowCount][4] = int(HabStressVar[rowCount][0]) # change in area
                ConsQualityArray[rowCount][4] = int(HabStressVar[rowCount][1]) # DQ of change in area
                ConsQualityArray2[rowCount][4] = int(HabStressVar[rowCount][1])+1 # DQ of change in area
                ConsequenceArray[rowCount][5] = int(HabStressVar[rowCount][2]) # change in structure
                ConsQualityArray[rowCount][5] = int(HabStressVar[rowCount][3]) # DQ of change in structure
                ConsQualityArray2[rowCount][5] = int(HabStressVar[rowCount][3])+1 # DQ of change in structure
                ConsequenceArray[rowCount][6] = int(HabStressVar[rowCount][4]) # frequency of natural disturbance
                ConsQualityArray[rowCount][6] = int(HabStressVar[rowCount][5]) # DQ of frequency of natural disturbance
                ConsQualityArray2[rowCount][6] = int(HabStressVar[rowCount][5])+1 # DQ of frequency of natural disturbance
                
                # ExposureArray: index 2 and 4 are blank for now (2 = spatial overlap rank; 4 = weighted average (Exp) or sum (DQ))               
                ExposureArray[rowCount][0] = int(StressVar[j][1]) # intensity
                ExpQualityArray[rowCount][0] = int(StressVar[j][2]) # DQ of intensity
                ExpQualityArray2[rowCount][0] = int(StressVar[j][2])+1 # DQ of intensity                
                ExposureArray[rowCount][1] = int(StressVar[j][3]) # management
                ExpQualityArray[rowCount][1] = int(StressVar[j][4]) # DQ of management
                ExpQualityArray2[rowCount][1] = int(StressVar[j][4])+1 # DQ of management
                ExpQualityArray[rowCount][2] = ((StressDataQual+HabDataQual)/2.0)# DQ average of spatial overlap
                ExpQualityArray2[rowCount][2] = ((StressDataQual2+HabDataQual2)/2.0)# DQ average of spatial overlap
                ExposureArray[rowCount][3] = int(HabStressVar[rowCount][6]) # overlap time
                ExpQualityArray[rowCount][3] = int(HabStressVar[rowCount][7]) # DQ of overlap time
                ExpQualityArray2[rowCount][3] = int(HabStressVar[rowCount][7])+1 # DQ of overlap time
                rowCount += 1

        ## FIX RECOVERY SCORING ##
        
        # hard-code the spatial overlap rank breaks
        UB_1 = 10
        UB_2 = 30

        # delete non-recovery specific columns from consequence array and redo calcs (except one for sums/avgs)
        ConsNumArray = np.where(ConsQualityArray == 0.0, 0.0,  ConsequenceArray/ConsQualityArray)
        ConsDenomArray = np.where(ConsNumArray == 0.0, 0.0,  1.0/ConsQualityArray)

        # create RecoveryArrays from [4-6] of ConsequenceArrays
        RecoveryArray = np.delete(ConsequenceArray, [4,5,6], axis=1)
        RecovQualityArray = np.delete(ConsQualityArray, [4,5,6], axis=1)
        RecovNumArray = np.delete(ConsNumArray, [4,5,6], axis=1)
        RecovDenomArray = np.delete(ConsDenomArray, [4,5,6], axis=1)
        RecovNumArray = np.where(RecovQualityArray == 0.0, 0.0,  RecoveryArray/RecovQualityArray)
        RecovDenomArray = np.where(RecovNumArray == 0.0, 0.0,  1.0/RecovQualityArray)

        for i in range(0,TotalHSCombo):
            RecovNumArray[i][4] = np.sum(RecovNumArray, axis=1)[i]
            RecovDenomArray[i][4] = np.sum(RecovDenomArray, axis=1)[i]
            if RecovDenomArray[i][4] == 0.0:
                RecoveryArray[i][4] = 0.0
            else:
                RecoveryArray[i][4] = RecovNumArray[i][4]/RecovDenomArray[i][4]

        ExpNumArray = np.where(ExpQualityArray == 0.0, 0.0,  ExposureArray/ExpQualityArray)
        ExpDenomArray = np.where(ExpNumArray == 0.0, 0.0,  1.0/ExpQualityArray)

        # sum up rows; divide rows and place weighted average value in original
        for i in range(0,TotalHSCombo):  
            ConsNumArray[i][7] = np.sum(ConsNumArray[i][:-1])
            ConsDenomArray[i][7] = np.sum(ConsDenomArray[i][:-1])
            if ConsDenomArray[i][7] == 0.0:
                ConsequenceArray[i][7] = 0.0
            else:
                ConsequenceArray[i][7] = ConsNumArray[i][7]/ConsDenomArray[i][7]

    except:
        gp.AddError(msgGetHabStressRatings)
        raise Exception


    ###############################################################    
    ############ OVERLAP RANKING AND RISK SCORING  ##############
    ###############################################################
    try:
        # add fields for risk calculations
        for i in range(0,StressCount):
            GS_HQ_area = AddField(GS_HQ_area, "OLP_RNK_S"+str(i+1), "SHORT", "", "")  ###
        for j in range(0,len(OverlapList)):
            GS_HQ_area = AddField(GS_HQ_area, "RISK_"+OverlapList[j], "DOUBLE", "8", "2") ###
        for k in range(0,HabCount):
            GS_HQ_area = AddField(GS_HQ_area, "CUMRISK_H"+str(k+1), "DOUBLE", "", "")  ###            
        GS_HQ_area = AddField(GS_HQ_area, "RECOV_HAB", "DOUBLE", "", "")
        GS_HQ_area = AddField(GS_HQ_area, "ECOS_RISK", "DOUBLE", "", "")

        cur = gp.UpdateCursor(GS_HQ_area)          
        row = cur.Next()
        while row:
            PredomHab = row.GetValue("PREDOM_HAB")
            if PredomHab <> 0:
                offset = 0
                OverlapCell = "no"
                for i in range(0,HabCount):
                    if row.GetValue("H"+str(i+1)+"_A") > 0.0:
                        for j in range(0,StressCount):
                            OverlapPct = row.GetValue(OverlapList[i+j+offset]+"_PCT")
                            if OverlapPct < UB_1 and OverlapPct <> 0.0: 
                                OverlapRank = 1.0
                            elif OverlapPct < UB_2 and OverlapPct >= UB_1:
                                OverlapRank = 2.0
                            elif OverlapPct >= UB_2:
                                OverlapRank = 3.0
                            else:
                                OverlapRank = 0.0
                            row.SetValue("OLP_RNK_S"+str(j+1), OverlapRank)
                            if OverlapRank > 0:
                                OverlapCell = "yes"
                                # calculate exposure and consequence averages                    
                                ExposureArray[i+j+offset][2] = OverlapRank
                                if ExpQualityArray[i+j+offset][2] == 0.0:
                                    ExpNumArray[i+j+offset][2] = 0.0
                                else:
                                    ExpNumArray[i+j+offset][2] = ExposureArray[i+j+offset][2]/ExpQualityArray[i+j+offset][2]
                                     
                                if ExpNumArray[i+j+offset][2] == 0.0:
                                    ExpDenomArray[i+j+offset][2] = 0.0
                                else: 
                                    ExpDenomArray[i+j+offset][2] = 1.0/ExpQualityArray[i+j+offset][2]
                                    
                                # sum up rows; divide rows and place weighted average value in original
                                ExpNumArray[i+j+offset][4] = np.sum(ExpNumArray[i+j+offset][:-1])
                                ExpDenomArray[i+j+offset][4] = np.sum(ExpDenomArray[i+j+offset][:-1])
                                if ExpDenomArray[i+j+offset][4] == 0.0:
                                    ExposureArray[i+j+offset][4] = 0.0
                                else:
                                    ExposureArray[i+j+offset][4] = ((ExpNumArray[i+j+offset][4])/(ExpDenomArray[i+j+offset][4]))

                                # get specific habitat-stressor risk score    
                                RiskScore = np.sqrt(((ExposureArray[i+j+offset][4]-1)**2) + ((ConsequenceArray[i+j+offset][7]-1)**2))
                                row.SetValue("RISK_"+OverlapList[i+j+offset], RiskScore)
                                              
                    offset = offset + (StressCount-1)

                if OverlapCell == "yes":
                    HabRiskList = []
                    for k in range(0,len(OverlapList)):
                        HabRiskList.append(row.GetValue("RISK_"+OverlapList[k]))
                
                    HabCounter = 0
                    while HabCounter < HabCount:
                        HabRiskTally = 0.0
                        for m in range((HabCounter*StressCount),((HabCounter*StressCount)+StressCount)):
                            HabRiskTally = HabRiskTally + HabRiskList[m]
                        row.SetValue("CUMRISK_H"+str(HabCounter+1), HabRiskTally)
                        HabCounter += 1
                        
                # set recovery score for each habitat that is present
                row.SetValue("RECOV_HAB", RecoveryArray[(PredomHab*StressCount)-StressCount][4])

            # sum up all habitat risk scores
            EcoRiskScore = 0.0
            for k in range(0,HabCount):
                EcoRiskScore = EcoRiskScore + row.GetValue("CUMRISK_H"+str(k+1))
            row.SetValue("ECOS_RISK", EcoRiskScore)
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur
        
    except:
        gp.AddError(msgCalcRiskMap)
        raise Exception

    try:
        # delete spatial overlap column from arrays and redo calcs
        ExposureArray = np.delete(ExposureArray, 2, axis=1)
        ExpQualityArray = np.delete(ExpQualityArray, 2, axis=1)
        ExpNumArray = np.delete(ExpNumArray, 2, axis=1)
        ExpDenomArray = np.delete(ExpDenomArray, 2, axis=1)
        ExpNumArray = np.where(ExpQualityArray == 0.0, 0.0,  ExposureArray/ExpQualityArray)
        ExpDenomArray = np.where(ExpNumArray == 0.0, 0.0,  1.0/ExpQualityArray)
 
        for i in range(0,TotalHSCombo):
            ExpNumArray[i][3] = np.sum(ExpNumArray, axis=1)[i]
            ExpDenomArray[i][3] = np.sum(ExpDenomArray, axis=1)[i]
            if ExpDenomArray[i][3] == 0.0:
                ExposureArray[i][3] = 0.0
            else:
                ExposureArray[i][3] = ExpNumArray[i][3]/ExpDenomArray[i][3]
            
            # Average Data Quality Arrays
            ExpQualitySum = np.sum(ExpQualityArray, axis=1)[i]
            ExpQualityCountList = list(ExpQualityArray[i])
            ExpQualityZeroCount = ExpQualityCountList.count(0)
            if ExpQualityZeroCount == 4:  # 4 and not 3 b/c last column is zero
                ExpQualityArray[i][3] = 0.0
            else:
                ExpQualityArray[i][3] = (ExpQualitySum / (4.0 - ExpQualityZeroCount))

            ConsQualitySum = np.sum(ConsQualityArray, axis=1)[i]
            ConsQualityCountList = list(ConsQualityArray[i])
            ConsQualityZeroCount = ConsQualityCountList.count(0)
            if ConsQualityZeroCount == 8: # 8 and not 7 b/c last column is zero
                ConsQualityArray[i][7] = 0.0
            else:
                ConsQualityArray[i][7] = (ConsQualitySum / (8.0 - ConsQualityZeroCount))

            # Average Data Quality Plotting Arrays
            ExpQualitySum2 = np.sum(ExpQualityArray2, axis=1)[i]
            ExpQualityArray2[i][4] = ExpQualitySum2/4.0

            ConsQualitySum2 = np.sum(ConsQualityArray2, axis=1)[i]
            ConsQualityArray2[i][7] = ConsQualitySum2/7.0

    except:
        gp.AddError(msgPlotArrays)
        raise Exception

    try:
        # create outputs
        gp.workspace = maps
        # create raster outputs of risk (where habitat exists)
        gp.Select_analysis(GS_HQ_area, GS_HQ_risk, "\"ECOS_RISK\" > 0")
        gp.Select_analysis(GS_HQ_area, GS_HQ_predom, "\"PREDOM_HAB\" > 0")

        # create raster output of each individual habitat's risk score from all stressors
        del3 = []
        for k in range(0,HabCount):
            if HabNoDataList[k] == "no":
                gp.FeatureToRaster_conversion(GS_HQ_area, "CUMRISK_H"+str(k+1), interws+"cr_h"+str(k+1), cellsize)       
                SetNullExp = "setnull("+interws+"cr_h"+str(k+1)+" <= 0 , "+interws+"cr_h"+str(k+1)+")"
                gp.SingleOutputMapAlgebra_sa(SetNullExp, "cum_risk_h"+str(k+1))
                del3.append("cr_h"+str(k+1))

        # create raster output for ecosystem risk and habitat recovery      
        gp.FeatureToRaster_conversion(GS_HQ_risk, "ECOS_RISK", ecosys_risk, cellsize)
        gp.FeatureToRaster_conversion(GS_HQ_predom, "RECOV_HAB", recov_potent, cellsize)

    except:
        gp.AddError(msgMapOutputs)
        raise Exception


    ##########################################################   
    ############## MATPLOT LIBRARY FUNCTIONS  ################
    ##########################################################
    try:
        if PlotBoolean == "true":
            try:
                from matplotlib import *
                from pylab import plt
            except:
                gp.AddError(msgMatplotlibNo)
                raise Exception

            gp.workspace = html_plots
            
            ExposureList = []
            ConsequenceList = []
            DataQualityList = []

            CumExposureList = np.zeros(HabCount, dtype=np.float64)
            CumConsequenceList = np.zeros(HabCount, dtype=np.float64)
            CumDataQualityList = np.zeros(HabCount, dtype=np.float64)

            for i in range(0,TotalHSCombo):
                if ExposureArray[i][3] == 0:
                    ExposureList.append(ExposureArray[i][3]+1)
                else:
                    ExposureList.append(ExposureArray[i][3])

                if ConsequenceArray[i][7] == 0:
                    ConsequenceList.append(ConsequenceArray[i][7]+1)
                else:
                    ConsequenceList.append(ConsequenceArray[i][7])

                DataQualityList.append(np.ceil((ExpQualityArray2[i][4]+ConsQualityArray2[i][7])/2.0))          

            # create risk plots
            CountX = 0
            CountY = 0
            plt.figure(1)
            if HabCount < 5:
                for i in range(0,HabCount): #1,2,3,4 = 2,2
                    if i > 0 and i < 2:
                        CountY = CountY + 1
                    elif i == 2:
                        CountX = CountX + 1
                        CountY = 0
                    elif i == 3:
                        CountY = CountY + 1
                        
                    plt.subplot2grid((2,2), (CountX,CountY))
                    
                    if HabCount <= 2:
                        plt.xlabel('Exposure')
                    elif HabCount > 2 and i > 1:
                        plt.xlabel('Exposure')
                    if i == 0:
                        plt.ylabel('Consequence')
                    elif HabCount > 2 and i == 2:
                        plt.ylabel('Consequence')
                    
                    plt.title(HabLyrList[i])
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=5, fc='0.15')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=4.75, fc='0.25')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=4.5, fc='0.25')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=4.25, fc='0.35')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=4, fc='0.35')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=3.75, fc='0.45')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=3.5, fc='0.45')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=3.25, fc='0.55')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=3, fc='0.55')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=2.75, fc='0.65')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=2.5, fc='0.65')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=2.25, fc='0.75')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=2, fc='0.75')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=1.75, fc='0.85')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=1.5, fc='0.85')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=1.25, fc='0.95')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=1, fc='0.95')
                    plt.gca().add_patch(cir)

                    for j in range((i*StressCount),(i*StressCount)+StressCount):
                        CumExposureList[i] = CumExposureList[i] + ExposureList[j]
                        CumConsequenceList[i] = CumConsequenceList[i] + ConsequenceList[j]
                        CumDataQualityList[i] = CumDataQualityList[i] + DataQualityList[j]
                        
                        if DataQualityList[j] == 1:
                            plt.plot(ExposureList[j],ConsequenceList[j], 'k^', markerfacecolor='black', markersize=8) # no score
                            plt.annotate('S'+str(j%StressCount+1),xy=(ExposureList[j]-.15,ConsequenceList[j]-.15),xytext=(ExposureList[j]+.05,ConsequenceList[j]))
                        elif DataQualityList[j] == 2:
                            plt.plot(ExposureList[j],ConsequenceList[j], 'gd', markerfacecolor='green', markersize=8) # best data
                            plt.annotate('S'+str(j%StressCount+1),xy=(ExposureList[j]-.15,ConsequenceList[j]-.15),xytext=(ExposureList[j]+.05,ConsequenceList[j]))
                        elif DataQualityList[j] == 3:
                            plt.plot(ExposureList[j],ConsequenceList[j], 'yo', markerfacecolor='yellow', markersize=8) # adequate data
                            plt.annotate('S'+str(j%StressCount+1),xy=(ExposureList[j]-.15,ConsequenceList[j]-.15),xytext=(ExposureList[j]+.05,ConsequenceList[j]))
                        else:
                            plt.plot(ExposureList[j],ConsequenceList[j], 'rs', markerfacecolor='red', markersize=8) # limited data
                            plt.annotate('S'+str(j%StressCount+1),xy=(ExposureList[j]-.15,ConsequenceList[j]-.15),xytext=(ExposureList[j]+.05,ConsequenceList[j]))

                    plt.plot(range(5), [1,1,1,1,1], 'k-')
                    plt.plot([1,1,1,1,1], range(5), 'k-')
                    plt.xlim(0.5,3.5)
                    plt.ylim(0.5,3.5)
                    plt.grid()

            else:
                for i in range(0,HabCount): #5,6,7,8 = 3,3
                    if i > 0 and i < 3:
                        CountY = CountY + 1
                    elif i == 3:
                        CountX = CountX + 1
                        CountY = 0
                    elif i > 3 and i < 6:
                        CountY = CountY + 1
                    elif i == 6:
                        CountX = CountX + 1
                        CountY = 0
                    elif i == 7:
                        CountY = CountY + 1

                    ax = plt.subplot2grid((3,3), (CountX,CountY))
                    plt.subplot2grid((3,3), (CountX,CountY))
                    if HabCount <= 3:
                        plt.xlabel('Exposure')
                    elif HabCount > 3 and HabCount < 7 and i > 2 and i < 6:
                        plt.xlabel('Exposure')
                    elif HabCount > 6 and i > 5:
                        plt.xlabel('Exposure')

                    if i == 0 or i == 3 or i == 6:
                        plt.ylabel('Consequence')

                    plt.title(HabLyrList[i])
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=5, fc='0.15')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=4.75, fc='0.25')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=4.5, fc='0.25')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=4.25, fc='0.35')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=4, fc='0.35')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=3.75, fc='0.45')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=3.5, fc='0.45')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=3.25, fc='0.55')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=3, fc='0.55')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=2.75, fc='0.65')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=2.5, fc='0.65')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=2.25, fc='0.75')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=2, fc='0.75')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=1.75, fc='0.85')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=1.5, fc='0.85')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), radius=1.25, fc='0.95')
                    plt.gca().add_patch(cir)
                    cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=1, fc='0.95')
                    plt.gca().add_patch(cir)

                    for j in range((i*StressCount),(i*StressCount)+StressCount):
                        CumExposureList[i] = CumExposureList[i] + ExposureList[j]
                        CumConsequenceList[i] = CumConsequenceList[i] + ConsequenceList[j]
                        CumDataQualityList[i] = CumDataQualityList[i] + DataQualityList[j]
                        
                        if DataQualityList[j] == 1:
                            plt.plot(ExposureList[j],ConsequenceList[j], 'k^', markerfacecolor='black', markersize=8) # no score
                            plt.annotate('S'+str(j%StressCount+1),xy=(ExposureList[j]-.15,ConsequenceList[j]-.15),xytext=(ExposureList[j]+.05,ConsequenceList[j]))
                        elif DataQualityList[j] == 2:
                            plt.plot(ExposureList[j],ConsequenceList[j], 'gd', markerfacecolor='green', markersize=8) # best data
                            plt.annotate('S'+str(j%StressCount+1),xy=(ExposureList[j]-.15,ConsequenceList[j]-.15),xytext=(ExposureList[j]+.05,ConsequenceList[j]))
                        elif DataQualityList[j] == 3:
                            plt.plot(ExposureList[j],ConsequenceList[j], 'yo', markerfacecolor='yellow', markersize=8) # adequate data
                            plt.annotate('S'+str(j%StressCount+1),xy=(ExposureList[j]-.15,ConsequenceList[j]-.15),xytext=(ExposureList[j]+.05,ConsequenceList[j]))
                        else:
                            plt.plot(ExposureList[j],ConsequenceList[j], 'rs', markerfacecolor='red', markersize=8) # limited data
                            plt.annotate('S'+str(j%StressCount+1),xy=(ExposureList[j]-.15,ConsequenceList[j]-.15),xytext=(ExposureList[j]+.05,ConsequenceList[j]))

                    plt.plot(range(5), [1,1,1,1,1], 'k-')
                    plt.plot([1,1,1,1,1], range(5), 'k-')       
                    plt.xlim(0.5,3.5)
                    plt.ylim(0.5,3.5)
                    xtl = ax.get_xticklabels()
                    ytl = ax.get_yticklabels()
                    xtl[1].set_visible(False)
                    xtl[2].set_visible(False)       
                    xtl[3].set_visible(False)
                    ytl[1].set_visible(False)
                    ytl[3].set_visible(False)  
                    plt.grid()

            plt.savefig(risk_plots, dpi=(1280/8))
            plt.clf()

            # create cumulative risk plot
            plt.figure(2)
            plt.subplot
            plt.xlabel('Exposure (Cumulative)')
            plt.ylabel('Consequence (Cumulative)')

            # divide cumulative quality total by number of habitat
            for i in range(0,len(CumDataQualityList)):
                CumDataQualityList[i] = np.ceil(CumDataQualityList[i]/HabCount)

            # find the max of both lists and set axes to be that length + 1
            if max(CumConsequenceList) >= max(CumExposureList):
                maxListValue = max(CumConsequenceList)
            else:
                maxListValue = max(CumExposureList)
                
            cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=(maxListValue+4), fc='0.15')
            plt.gca().add_patch(cir)
            cir = plt.Circle((0,0), radius=(maxListValue+4)*0.9, fc='0.25')
            plt.gca().add_patch(cir)
            cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=(maxListValue+4)*0.8, fc='0.35')
            plt.gca().add_patch(cir)
            cir = plt.Circle((0,0), radius=(maxListValue+4)*0.7, fc='0.45')
            plt.gca().add_patch(cir)
            cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=(maxListValue+4)*0.6, fc='0.55')
            plt.gca().add_patch(cir)
            cir = plt.Circle((0,0), radius=(maxListValue+4)*0.5, fc='0.65')
            plt.gca().add_patch(cir)
            cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=(maxListValue+4)*0.4, fc='0.75')
            plt.gca().add_patch(cir)
            cir = plt.Circle((0,0), radius=(maxListValue+4)*0.3, fc='0.85')
            plt.gca().add_patch(cir)
            cir = plt.Circle((0,0), edgecolor='.25', linestyle ='dashed', radius=(maxListValue+4)*0.2, fc='0.95')
            plt.gca().add_patch(cir)
            for i in range(0,HabCount):
                if CumDataQualityList[i] == 1:
                    plt.plot(CumExposureList[i],CumConsequenceList[i], 'k^', markerfacecolor='black', markersize=8) # no score
                    plt.annotate('H'+str(i+1),xy=(CumExposureList[i]-.15,CumConsequenceList[i]-.15),xytext=(CumExposureList[i]+.05,CumConsequenceList[i]))
                elif CumDataQualityList[i] == 2:
                    plt.plot(CumExposureList[i],CumConsequenceList[i], 'gd', markerfacecolor='green', markersize=8) # best data
                    plt.annotate('H'+str(i+1),xy=(CumExposureList[i]-.15,CumConsequenceList[i]-.15),xytext=(CumExposureList[i]+.05,CumConsequenceList[i]))
                elif CumDataQualityList[i] == 3:
                    plt.plot(CumExposureList[i],CumConsequenceList[i], 'yo', markerfacecolor='yellow', markersize=8) # adequate data
                    plt.annotate('H'+str(i+1),xy=(CumExposureList[i]-.15,CumConsequenceList[i]-.15),xytext=(CumExposureList[i]+.05,CumConsequenceList[i]))
                else:
                    plt.plot(CumExposureList[i],CumConsequenceList[i], 'rs', markerfacecolor='red', markersize=8) # limited data
                    plt.annotate('H'+str(i+1),xy=(CumExposureList[i]-.15,CumConsequenceList[i]-.15),xytext=(CumExposureList[i]+.05,CumConsequenceList[i]))

            OnesList = np.ones(maxListValue+4, dtype=np.float64)
            plt.plot(range(maxListValue+4), OnesList, 'k-')
            plt.plot(OnesList, range(maxListValue+4), 'k-')
            plt.xlim(0,maxListValue+1)
            plt.ylim(0,maxListValue+1)
            plt.grid()
            plt.savefig(ecosysRisk_plots, dpi=(1280/8))
            plt.clf()

            # create html file
            htmlfile = open(outputHTML, "w")
            htmlfile.write("<html>\n")
            htmlfile.write("<title>Marine InVEST</title>\n")
            htmlfile.write("<CENTER><H1>Visualizing the InVEST Habitat Risk Assessment Model</H1></CENTER>\n")
            htmlfile.write("<br><b>A note on data quality and uncertainty:</b>  Ecological risk assessment is an integrative process, \
                            which requires a substantial amount of data on many attributes of human and ecological systems. \
                            It is likely that some aspects of the risk assessment will be supported high quality data and others \
                            aspects will be subject to limited data availability and high uncertainty. To increase the transparency \
                            of the model results, we color-code the results in the figures below according to the average quality of \
                            the data that was used to generate each score. Points in the figures below represented by red squares are \
                            derived from low quality data and should be interpreted with caution, while those represented by green diamonds\
                            are derived from high quality data and may be interpreted with high confidence. We hope that by displaying data\
                            quality, users will be aware of some sources of uncertainty in the risk assessment, and will therefore be cautious \
                            when using results derived from low quality data. In addition, this information can be used to guide research and \
                            monitoring effects to improve data quality and availability.<br>\n")
            htmlfile.write("<br><HR><H2>Cumulative Ecosystem Risk Plot</H2>\n")
            htmlfile.write("<table border=\"0\"><tr><td>")
            htmlfile.write("<img src=\"plot_ecosys_risk.png\" width=\"960\" height=\"720\">")
            htmlfile.write("</td><td>")
            htmlfile.write("This figure shows the cumulative risk for each habitat in the study region. This figure can be used to determine \
                            which habitats are at highest risk from human activities, and if this risk is mostly due to high cumulative exposure\
                            (exogenous factors which can be mitigated by management) or high cumulative consequence (endogenous factors which are\
                            less responsive to human intervention).  Cumulative exposure scores are derived by summing the exposure scores for each\
                            stressor. They represent the total exposure of each habitat to all stressors in study region. For example, the cumulative\
                            exposure for eelgrass in a study region with shellfish aquaculture and destructive fishing is the sum of the shellfish \
                            aquaculture exposure score for eelgrass and the destructive fishing exposure score for eelgrass. Cumulative consequence \
                            scores are derived in the same way. Habitats with high cumulative exposure and high cumulative consequence are at the \
                            highest risk from human activities.<p>")
            htmlfile.write("<img src=\"file:\\"+os.path.dirname(sys.argv[0])+"\\HQM_plotLegend.png\" width=\"152\" height=\"88\"><p>\n")
            for i in range(0,len(HabLyrList)):
                htmlfile.write("<big><u>H"+str(i+1)+"</u>: "+str(HabLyrList[i])+"</big><br>\n")
            htmlfile.write("</td></tr></table>\n")
            htmlfile.write("<br><HR><H2>Risk Plots for Each Habitat</H2>\n")
            htmlfile.write("<table border=\"0\"><tr><td>")
            htmlfile.write("<img src=\"plots_risk.png\" width=\"960\" height=\"720\">")
            htmlfile.write("</td><td>")
            htmlfile.write("These figures show the exposure and consequence scores for each stressor and habitat combination in the study region. \
                            Stressors that have high exposure scores and high consequence scores pose the greatest risk to habitats. Reducing risk \
                            through management is likely to be more effective in situations where high risk is driven by high exposure, not high consequence.<p>")
            htmlfile.write("<img src=\"file:\\"+os.path.dirname(sys.argv[0])+"\\HQM_plotLegend.png\" width=\"152\" height=\"88\"><p>\n")
            for j in range(0,len(StressLyrList)):
                htmlfile.write("<big><u>S"+str(j+1)+"</u>: "+str(StressLyrList[j])+"</big><br>\n")
            htmlfile.write("</td></td></table>\n")
            htmlfile.write("</html>")
            htmlfile.close()
    except:
        gp.AddError(msgPlotHTMLOutputs)
        raise Exception


    ##########################################################   
    ######## GENERATE HABITAT MAPS OF RISK HOTSPOTS  #########
    ##########################################################

    if RiskBoolean == "true":
        gp.AddMessage("\nGenerating habitat maps of risk hotspots...")
        
        # copy 'GS_HQ_area' and erase superfluous attributes before intersection
        gp.CopyFeatures_management(GS_HQ_area, GS_HQ_intersect, "", "0", "0", "0")
        keepFieldList = ["FID", "Shape", "CELL_SIZE"]
        eraseFieldList = []
        for i in range(0,HabCount):
            keepFieldList.append("CUMRISK_H"+str(i+1))
            for j in range(0,StressCount):
                keepFieldList.append("RISK_H"+str(i+1)+"S"+str(j+1))
        
        fields = gp.ListFields(GS_HQ_intersect, "*")
        fc_field = fields.Next()
        while fc_field:
            if fc_field.name not in keepFieldList:
                eraseFieldList.append(fc_field.name)
            fc_field = fields.Next()
        del fc_field
 
        EraseFieldExpr = eraseFieldList[0]
        for i in range(1,len(eraseFieldList)):
            EraseFieldExpr = EraseFieldExpr+";"+str(eraseFieldList[i])
        gp.DeleteField_management(GS_HQ_intersect, EraseFieldExpr)

        # intersect 'GS_HQ_area' with each habitat input and genrate risk hotspots
        for i in range(0,len(HabLyrList)):
            HabVariable = Hab_Directory+"\\"+HabLyrList[i]
            IntersectExpr = HabVariable+" 1; "+GS_HQ_intersect+" 2"
            gp.Intersect_analysis(IntersectExpr, maps+"h"+str(i+1)+"_"+HabLyrList[i][:-5]+"Risk.shp", "NO_FID", "", "INPUT")
            for j in range(0,StressCount):
                gp.AddField_management(maps+"h"+str(i+1)+"_"+HabLyrList[i][:-5]+"Risk.shp", "S"+str(j+1)+"RISKNUM", "SHORT", "0", "0", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
            gp.AddField_management(maps+"h"+str(i+1)+"_"+HabLyrList[i][:-5]+"Risk.shp", "CRISK_NUM", "SHORT", "0", "0", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
            gp.AddField_management(maps+"h"+str(i+1)+"_"+HabLyrList[i][:-5]+"Risk.shp", "RISK_NUM", "SHORT", "0", "0", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
            gp.AddField_management(maps+"h"+str(i+1)+"_"+HabLyrList[i][:-5]+"Risk.shp", "RISK_QUAL", "TEXT", "10", "0", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
            
            cur = gp.UpdateCursor(maps+"h"+str(i+1)+"_"+HabLyrList[i][:-5]+"Risk.shp")          
            row = cur.Next()
            while row:
                # individual stressor risk logic
                for j in range(0,StressCount):
                    if row.GetValue("RISK_H"+str(i+1)+"S"+str(j+1)) < (np.sqrt(8.0)*(1.0/3.0)):
                        row.SetValue("S"+str(j+1)+"RISKNUM", 1)
                    elif row.GetValue("RISK_H"+str(i+1)+"S"+str(j+1)) >= (np.sqrt(8.0)*(1.0/3.0)) and row.GetValue("RISK_H"+str(i+1)+"S"+str(j+1)) < (np.sqrt(8.0)*(2.0/3.0)):
                        row.SetValue("S"+str(j+1)+"RISKNUM", 2)
                    else:
                        row.SetValue("S"+str(j+1)+"RISKNUM", 3)        
                    
                # cumulative risk logic
                if row.GetValue("CUMRISK_H"+str(i+1)) < (np.sqrt(8.0*StressCount)*(1.0/3.0)):
                    row.SetValue("CRISK_NUM", 1)
                elif row.GetValue("CUMRISK_H"+str(i+1)) >= (np.sqrt(8.0*StressCount)*(1.0/3.0)) and row.GetValue("CUMRISK_H"+str(i+1)) < (np.sqrt(8.0*StressCount)*(2.0/3.0)):
                    row.SetValue("CRISK_NUM", 2)
                else:
                    row.SetValue("CRISK_NUM", 3)

                # determine risk hotspots ratings
                HotSpotRatingList = []
                for k in range(0,StressCount):
                    HotSpotRatingList.append(row.GetValue("S"+str(k+1)+"RISKNUM"))
                HotSpotRatingList.append(row.GetValue("CRISK_NUM"))    
                row.SetValue("RISK_NUM", max(HotSpotRatingList))
                if row.GetValue("RISK_NUM") == 1:
                    row.SetValue("RISK_QUAL", "Low")
                elif row.GetValue("RISK_NUM") == 2:
                    row.SetValue("RISK_QUAL", "Medium")
                else:
                    row.SetValue("RISK_QUAL", "High")            
                cur.UpdateRow(row)
                row = cur.next()
            del row
            del cur
            
    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(outputws+"parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("HABITAT RISK ASSESSMENT MODEL PARAMETERS\n")
    parafile.writelines("___________________________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

    # delete superfluous intermediate data
    gp.workspace = interws
    del1 = [GS_HQ, GS_HQ_lyr, GS_rst, GS_HQ_risk, GS_HQ_predom]
    del2 = []
    for i in range(0,len(HabLyrList)):
        for j in range(0,len(StressLyrList)):
            del2.append("H"+str(i+1)+"S"+str(j+1))
    deletelist = del1 + del2 + del3 + del_hab + del_stress
    for data in deletelist:
        if gp.exists(data):
            gp.delete_management(data)
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())