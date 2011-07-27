# Marine InVEST: Habitat Risk Assessment Model
# Authors: Gregg Verutes, Joey Bernhardt, Katie Arkema, Jeremy Davies 
# 05/09/11

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
msgCheckInputs = "\nError checking and preparing inputs."
msgBuildVAT = "\nError building VAT for zonal statistics raster.  Try opening a new ArcMap session and re-run the HRA model."
msgPrepHSLayers ="\nError preparing habitat and stressor input layers."
msgCheckHSLayers = "\nError checking habitat and stressor input layers."
msgBuffRastHULayers = "\nError buffering and rasterizing stressor and habitat input layers."
msgOverlapPredomHab = "\nError calculating spatial overlap and predominant habitat."
msgGetHabStressRatings = "\nError obtaining habitat and stressor ratings from table."
msgCalcRiskMap = "\nError calculating risk scores within analysis grid (GS)."
msgPlotArrays = "\nError preparing exposure and consequence matrix values for plotting."
msgMapOutputs = "\nError generating map outputs."
msgPlotHTMLOutputs = "\nError generating plot and HTML outputs."
msgNumPyNo = "NumPy extension is required to run the Habitat Risk Assessment Model.  Please consult the Marine InVEST FAQ document for instructions on how to install."
msgWin32ComNo = "PythonWin extension is required to run the Habitat Risk Assessment Model.  Please consult the Marine InVEST FAQ document for instructions on how to install."
msgMatplotlibNo = "Matplotlib extension (version 1.0 or newer) is required to run the Habitat Risk Assessment Model.  Please consult the Marine InVEST FAQ document for instructions on how to install."

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
        Hab_Directory = gp.GetParameterAsText(2)
        parameters.append("Habitat Data Directory: "+ Hab_Directory)
        Stress_Directory = gp.GetParameterAsText(3)
        parameters.append("Stressor Data Directory: "+ Stress_Directory)
        HabStressRate_Table = gp.GetParameterAsText(4)
        parameters.append("Stressor Data Directory: "+ HabStressRate_Table)
        PlotBoolean = gp.GetParameterAsText(5)
        parameters.append("Create HTML output with risk plots?: "+ PlotBoolean)
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
    maps = outputws + os.sep + "maps" + os.sep
    html_plots = outputws + os.sep + "html_plots" + os.sep
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
    GS_HQ_lyr = interws + "GS_HQ_lyr.lyr"
    GS_rst = interws + "gs_rst"
    GS_HQ_risk = interws + "GS_HQ_risk.shp"
    GS_HQ_predom = interws + "GS_HQ_predom.shp"
    GS_HQ_area = interws + "GS_HQ_area.shp"

    # output
    ecosys_risk = maps + "ecosys_risk"
    predom_hab = maps + "predom_hab"
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
        xlApp = Dispatch("Excel.Application")
        xlApp.Visible=0
        xlApp.DisplayAlerts=0
        xlApp.Workbooks.Open(HabStressRate_Table)
        cell = xlApp.Worksheets("(1)")
        HabInputCount = int(cell.Range("g11").Value)
        StressInputCount = int(cell.Range("b13").Value)

        gp.AddMessage("... Habitat Count = "+str(HabInputCount))
        gp.AddMessage("... Stressor Count = "+str(StressInputCount))

        StressBuffDistList = []
        counter = 2
        while counter < StressInputCount+2:
            StressBuffDistList.append(int(cell.Range("e"+str(counter)).Value))
            counter += 1
            
        xlApp.ActiveWorkbook.Close(SaveChanges=0)
        xlApp.Quit()

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

        HabNoDataList = []
        for i in range(0,len(HabLyrList)):
            HabVariable = Hab_Directory+"\\"+HabLyrList[i]
            checkProjections(HabVariable)
            HabVariable = AddField(HabVariable, "VID", "SHORT", "0", "0")        
            gp.CalculateField_management(HabVariable, "VID", "[FID]+1", "VB")
            gp.FeatureToRaster_conversion(HabVariable, "VID", "Hab_"+str(i+1), "50")
            try:
                if gp.GetCount("Hab_"+str(i+1)) == 0:
                    HabNoDataList.append("yes")
                else:
                    HabNoDataList.append("no")
            except:
                gp.BuildRasterAttributeTable_management("Hab_"+str(i+1), "Overwrite")
                if gp.GetCount("Hab_"+str(i+1)) == 0:
                    HabNoDataList.append("yes")
                else:
                    HabNoDataList.append("no")
            
        for i in range(0,len(StressLyrList)):
            StressVariable = Stress_Directory+"\\"+StressLyrList[i]
            checkProjections(StressVariable)
            StressVariable = AddField(StressVariable, "VID", "SHORT", "0", "0")        
            gp.CalculateField_management(StressVariable, "VID", "[FID]+1", "VB")
            if StressBuffDistList[i] > 0:
                gp.Buffer_analysis(StressVariable, StressLyrBuffList[i], str(StressBuffDistList[i]) + " Meters", "FULL", "ROUND", "NONE", "")        
                gp.FeatureToRaster_conversion(StressLyrBuffList[i], "VID", "Stress_"+str(i+1), "50")
            else:
                gp.FeatureToRaster_conversion(StressVariable, "VID", "Stress_"+str(i+1), "50")
            
    except:
        gp.AddError(msgBuffRastHULayers)
        raise Exception        

    #############################################################    
    ############ CALCULATE OVERLAP AND PREDOM HAB  ##############
    #############################################################
    try:
        gp.AddMessage("\nCalculating spatial overlap and predominant habitat...")
        OverlapList = []
        OverlapNoDataList = []
        for i in range(0,len(HabLyrList)):
            for j in range(0,len(StressLyrList)):
                CmbExpr = "Hab_"+str(i+1)+";"+"Stress_"+str(j+1)
                gp.Combine_sa(CmbExpr, "H"+str(i+1)+"S"+str(j+1))
                OverlapList.append("H"+str(i+1)+"S"+str(j+1))
                try:
                    if gp.GetCount("H"+str(i+1)+"S"+str(j+1)) == 0:
                        OverlapNoDataList.append("yes")
                    else:
                        OverlapNoDataList.append("no")
                except:
                    gp.BuildRasterAttributeTable_management("H"+str(i+1)+"S"+str(j+1), "Overwrite")
                    if gp.GetCount("H"+str(i+1)+"S"+str(j+1)) == 0:
                        OverlapNoDataList.append("yes")
                    else:
                        OverlapNoDataList.append("no")

        for i in range(0,len(HabLyrList)):
            GS_HQ = AddField(GS_HQ, "H"+str(i+1)+"_A", "DOUBLE", "8", "2")
        for i in range(0,len(OverlapList)):
            GS_HQ = AddField(GS_HQ, OverlapList[i]+"_A", "DOUBLE", "8", "2")
            GS_HQ = AddField(GS_HQ, OverlapList[i]+"_PCT", "DOUBLE", "8", "2")

        gp.MakeFeatureLayer_management(GS_HQ, GS_HQ_lyr, "", "", "")
        gp.snapRaster = GS_rst
        gp.cellsize = "MINOF"

        for i in range(0,len(HabLyrList)):
            if HabNoDataList[i] == "no":
                gp.ZonalStatisticsAsTable_sa(GS_rst, "VALUE", "Hab_"+str(i+1), "zs_H"+str(i+1)+".dbf", "DATA")
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
        gp.AddError(msgOverlapPredomHab)
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
        ExpQualityArray2 = np.reshape(ExpQualityList2, (TotalHSCombo,5)) # 1-4 range
        ConsequenceList = np.zeros(TotalHSCombo*8, dtype=np.float64)
        ConsequenceArray = np.reshape(ConsequenceList, (TotalHSCombo,8))
        ConsQualityList = np.zeros(TotalHSCombo*8, dtype=np.float64)
        ConsQualityArray = np.reshape(ConsQualityList, (TotalHSCombo,8)) # 1-3 range
        ConsQualityList2 = np.zeros(TotalHSCombo*8, dtype=np.float64)
        ConsQualityArray2 = np.reshape(ConsQualityList2, (TotalHSCombo,8)) # 1-3 range

        x2App = Dispatch("Excel.Application")
        x2App.Visible=0
        x2App.DisplayAlerts=0
        x2App.Workbooks.Open(HabStressRate_Table)
        cell = x2App.Worksheets("Calc (exp)")
        rowCount = 0
        for i in range(3,3+(HabCount*10),10):
            for j in range(0,StressCount):
                # index 2 and 3 are blank for now (2 = spatial overlap; 4 = weighted average (Exp) or sum (DQ))
                ExposureArray[rowCount][0] = cell.Range("q"+str(i+j)).Value # intensity 
                ExposureArray[rowCount][1] = cell.Range("r"+str(i+j)).Value # management
                ExpQualityArray[rowCount][0] = cell.Range("s"+str(i+j)).Value # DQ of intensity 
                ExpQualityArray[rowCount][1] = cell.Range("t"+str(i+j)).Value # DQ of management
                StressDataQual = cell.Range("u"+str(i+j)).Value # DQ of spatial overlap for stressor
                HabDataQual = cell.Range("v"+str(i+j)).Value # DQ of spatial overlap for habitat
                ExpQualityArray[rowCount][2] = ((StressDataQual+HabDataQual)/2.0)# DQ average of spatial overlap
                ExpQualityArray2[rowCount][0] = cell.Range("i"+str(i+j)).Value # DQ of intensity 
                ExpQualityArray2[rowCount][1] = cell.Range("j"+str(i+j)).Value # DQ of management
                StressDataQual2 = cell.Range("e"+str(i+j)).Value # DQ of spatial overlap for stressor
                HabDataQual2 = cell.Range("f"+str(i+j)).Value # DQ of spatial overlap for habitat
                ExpQualityArray2[rowCount][2] = ((StressDataQual2+HabDataQual2)/2.0)# DQ average of spatial overlap
                rowCount += 1
            
        cell2 = x2App.Worksheets("Calc (con)")
        rowCount = 0
        for i in range(3,3+(HabCount*10),10):
            for j in range(0,StressCount):
                # index 7 is blank (7 = weighted average (Cons) or sum (DQ))
                ConsequenceArray[rowCount][0] = cell2.Range("x"+str(i+j)).Value # natural mortality rate
                ConsequenceArray[rowCount][1] = cell2.Range("y"+str(i+j)).Value # recruitment pattern
                ConsequenceArray[rowCount][2] = cell2.Range("z"+str(i+j)).Value # connectivity
                ConsequenceArray[rowCount][3] = cell2.Range("aa"+str(i+j)).Value # age at maturity or recovery time
                ConsequenceArray[rowCount][4] = cell2.Range("ab"+str(i+j)).Value # change in area
                ConsequenceArray[rowCount][5] = cell2.Range("ac"+str(i+j)).Value # change in structure
                ConsequenceArray[rowCount][6] = cell2.Range("ad"+str(i+j)).Value # frequency of natural disturbance
                ExposureArray[rowCount][3] = cell2.Range("ae"+str(i+j)).Value # overlap time
                ConsQualityArray[rowCount][0] = cell2.Range("af"+str(i+j)).Value # DQ of natural mortality rate
                ConsQualityArray[rowCount][1] = cell2.Range("ag"+str(i+j)).Value # DQ of recruitment pattern
                ConsQualityArray[rowCount][2] = cell2.Range("ah"+str(i+j)).Value # DQ of connectivity
                ConsQualityArray[rowCount][3] = cell2.Range("ai"+str(i+j)).Value # DQ of age at maturity or recovery time
                ConsQualityArray[rowCount][4] = cell2.Range("aj"+str(i+j)).Value # DQ of change in area
                ConsQualityArray[rowCount][5] = cell2.Range("ak"+str(i+j)).Value # DQ of change in structure
                ConsQualityArray[rowCount][6] = cell2.Range("al"+str(i+j)).Value # DQ of frequency of natural disturbance
                ExpQualityArray[rowCount][3] = cell2.Range("am"+str(i+j)).Value # DQ of overlap time
                ConsQualityArray2[rowCount][0] = cell2.Range("m"+str(i+j)).Value # DQ of natural mortality rate
                ConsQualityArray2[rowCount][1] = cell2.Range("n"+str(i+j)).Value # DQ of recruitment pattern
                ConsQualityArray2[rowCount][2] = cell2.Range("o"+str(i+j)).Value # DQ of connectivity
                ConsQualityArray2[rowCount][3] = cell2.Range("p"+str(i+j)).Value # DQ of age at maturity or recovery time
                ConsQualityArray2[rowCount][4] = cell2.Range("q"+str(i+j)).Value # DQ of change in area
                ConsQualityArray2[rowCount][5] = cell2.Range("r"+str(i+j)).Value # DQ of change in structure
                ConsQualityArray2[rowCount][6] = cell2.Range("s"+str(i+j)).Value # DQ of frequency of natural disturbance
                ExpQualityArray2[rowCount][3] = cell2.Range("t"+str(i+j)).Value # DQ of overlap time
                rowCount += 1

        # grab ranges for spatial overlap ratings
        cell3 = x2App.Worksheets("Rating legends")
        UB_1 = cell3.Range("d43").Value
        UB_2 = cell3.Range("d44").Value 
            
        x2App.ActiveWorkbook.Close(SaveChanges=0)
        x2App.Quit()

        # delete non-recovery specific columns from consequence array and redo calcs (except one for sums/avgs)
        ConsNumArray = np.where(ConsQualityArray == 0.0, 0.0,  ConsequenceArray/ConsQualityArray)
        ConsDenomArray = np.where(ConsNumArray == 0.0, 0.0,  1.0/ConsQualityArray)
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
            
            # Average Data Quality Arrays
            RecovQualitySum = np.sum(RecovQualityArray, axis=1)[i]
            RecovQualityCountList = list(RecovQualityArray[i])
            RecovQualityZeroCount = RecovQualityCountList.count(0)
            if RecovQualityZeroCount == 3:
                RecovQualityArray[i][4] = 0.0
            else:
                RecovQualityArray[i][4] = (RecovQualitySum / (3.0 - RecovQualityZeroCount))

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
            GS_HQ_area = AddField(GS_HQ_area, "OLP_RNK_S"+str(i+1), "SHORT", "", "")
        for j in range(0,len(OverlapList)):
            GS_HQ_area = AddField(GS_HQ_area, "RISK_"+OverlapList[j], "DOUBLE", "8", "2")
        for k in range(0,HabCount):
            GS_HQ_area = AddField(GS_HQ_area, "CUMRISK_H"+str(k+1), "DOUBLE", "", "")            
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

        for k in range(0,HabCount):   
            gp.FeatureToRaster_conversion(GS_HQ_area, "CUMRISK_H"+str(k+1), "cum_risk_h"+str(k+1), cellsize)

        gp.FeatureToRaster_conversion(GS_HQ_risk, "ECOS_RISK", ecosys_risk, cellsize)
        gp.FeatureToRaster_conversion(GS_HQ_predom, "PREDOM_HAB", predom_hab, cellsize)
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
                ExposureList.append(ExposureArray[i][3])
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
    deletelist = del1 + del2
    for data in deletelist:
        if gp.exists(data):
            gp.delete_management(data)
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())