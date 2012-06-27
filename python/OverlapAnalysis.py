# Marine InVEST: Overlap Analysis
# Authors: Gregg Verutes, Mike Papenfus, Jodie Toft
# 08/17/11

# import modules
import sys, string, os, datetime, csv, shlex
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
msgCheckInputs = "\nError checking and preparing inputs."
msgBuildVAT = "\nError building VAT for zonal statistics raster.  Try opening a new ArcMap session and re-run the Overlap Analysis model."
msgPrepHULayers = "\nError preparing human uses input layers."
msgCheckHULayers = "\nError checking human uses input layers."
msgBuffRastHULayers = "\nError buffering and rasterizing human uses input layers."
msgDistanceDecay = "\nError calculating distances for decay function."
msgImportScoring = "\nError calculating importance scoring."
msgCreateOutputs = "\nError creating outputs."
msgNumPyNo = "NumPy extension is required to run the Human Uses Overlay Analysis.  Please consult the Marine InVEST FAQ for instructions on how to install."
msgWin32ComNo = "PythonWin extension is required to run the Human Uses Overlay Analysis.  Please consult the Marine InVEST FAQ for instructions on how to install."

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
        AnalysisType = gp.GetParameterAsText(1)
        parameters.append("Type of Analysis Zones: "+ AnalysisType)
        AnalysisZones = gp.GetParameterAsText(2)
        parameters.append("Analysis Zones Layer: "+ AnalysisZones)
        HU_Directory = gp.GetParameterAsText(3)
        parameters.append("Overlap Analysis Data Directory: "+ HU_Directory)
        HULayersTable = gp.GetParameterAsText(4)
        parameters.append("Overlap Analysis Layers Table: "+ HULayersTable)
        WghtFieldName = gp.GetParameterAsText(5)
        parameters.append("Importance Score Field Name: "+ WghtFieldName)
        PopPlaces = gp.GetParameterAsText(6)
        parameters.append("Points Layer Indicating Location of Human Use Hubs: "+ PopPlaces)
        DecayRate = gp.GetParameterAsText(7)
        parameters.append("Distance Decay Rate: "+ str(DecayRate))
        if DecayRate:
            DecayRate = float(gp.GetParameterAsText(7))
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

    # intermediate
    AnalysisZonesHU = interws + "AnalysisZonesHU.shp"
    AnalysisZonesHU_Pts = interws + "AnalysisZonesHU_Pts.shp"
    AnalysisZonesHU_lyr = interws + "AnalysisZonesHU_lyr.lyr"
    AZHU_rst = interws + "azhu_rst"
    AnalysisZonesHU_area = interws + "AnalysisZonesHU_area.shp"
    hu_freq_all = interws + "hu_freq_all"

    # output
    hu_freq = outputws + "hu_freq"
    hu_impscore = outputws + "hu_impscore"
    HU_calcs_csv = outputws + "HU_calcs.csv"


    ##############################################
    ###### COMMON FUNCTION AND CHECK INPUTS ######
    ##############################################

    def checkProjections(thedata):
        dataDesc = gp.describe(thedata)
        spatreflc = dataDesc.SpatialReference
        if spatreflc.Type <> 'Projected':
            gp.AddError(thedata+" does not appear to be projected.  It is assumed to be in meters.")
        if spatreflc.LinearUnitName <> 'Meter':
            gp.AddError("This model assumes that "+thedata+" is projected in meters for area calculations.  You may get erroneous results.")
            raise Exception

    def checkPathLength(path):
        pathLength = len(path)
        if pathLength > 50:
            gp.AddError("In order to properly build a value attribute table (VAT) for the gridded seascape, the path length of your workspace must be shorten by at least "+str(pathLength-50)+" character(s).")
            raise Exception

    def grabProjection(data):
        dataDesc = gp.describe(data)
        sr = dataDesc.SpatialReference
        gp.OutputCoordinateSystem = sr
        strSR = str(gp.OutputCoordinateSystem)
        return strSR

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

    def FieldExists(FileName, FieldName):
        if gp.Exists(FileName):
            fieldcount = 0
            fieldList = gp.ListFields(FileName, FieldName, "All")
            field = fieldList.Next()
            while field <> None:
                if field.Name == FieldName:
                    fieldcount = fieldcount + 1
                field = fieldList.Next()
            del fieldList
            del field
            if fieldcount <> 0:
                pass
            else:
                gp.AddError("A field called "+FieldName+" does not exist in "+FileName+".")
                raise Exception 
        else:
            gp.AddError(FileName+"does not exist in the human uses layers directory.")
            raise Exception

    try:
        gp.AddMessage("\nPreparing and checking human uses input data...")
        # check to make sure workspace path length doesn't exceed limit for building VAT
        checkPathLength(gp.GetParameterAsText(0))
        # copy and populate input FC attribute table
        gp.CopyFeatures_management(AnalysisZones, AnalysisZonesHU, "", "0", "0", "0")

        if AnalysisType == "Gridded Seascape (GS)":
            # grab cellsize info
            cur = gp.UpdateCursor(AnalysisZonesHU)
            row = cur.Next()
            cellsize = int(row.GetValue("CELL_SIZE"))
            del cur
            del row
        else:
            checkProjections(AnalysisZonesHU)
            AnalysisZonesHU = AddField(AnalysisZonesHU, "VALUE", "LONG", "0", "0")
            gp.CalculateField_management(AnalysisZonesHU, "VALUE", "[FID]+1", "VB")
            cellsize = 250

        # make feature to check that inputs overlap
        gp.MakeFeatureLayer_management(AnalysisZonesHU, AnalysisZonesHU_lyr, "", "", "")
            
        if PopPlaces:
            checkProjections(PopPlaces) 
        # grab projection
        projection = grabProjection(AnalysisZonesHU)
        gp.FeatureToRaster_conversion(AnalysisZonesHU, "VALUE", AZHU_rst, str(cellsize))
    except:
        gp.AddError(msgCheckInputs)
        raise Exception

    try:
        # build vat
        sDesc = gp.describe(AZHU_rst)
        sInRasterDS = sDesc.catalogpath
        sExpression = "BuildVat " + sInRasterDS
        gp.MultiOutputMapAlgebra_sa(sExpression)
    except:
        gp.AddError(msgBuildVAT)
        raise Exception


    ############################################           
    ######## PREPARE HUMAN USES LAYERS #########
    ############################################
    try:
        xlApp = Dispatch("Excel.Application")
        xlApp.Visible = 0
        xlApp.DisplayAlerts = 0
        xlBook = xlApp.Workbooks.Open(HULayersTable[:-(1+len(HULayersTable.split("\\")[-1]))])
        HUPath = HULayersTable.split("\\")
        HUSheet = HUPath[-1]
        xlSheet1 = xlBook.Worksheets(HUSheet[:-1])
        HUInputCount = int(xlSheet1.Range("b21").Value) # total number of human uses layers
        HUSSIDList = []
        InterHUImpList = []
        HUBuffDistList = []
        counter = 0
        row = 2
        col = 1
        while row < 20:
            if xlSheet1.Cells(row, col).Value <> None and xlSheet1.Cells(row, col).Value <> '':
                HUSSIDList.append(int(xlSheet1.Cells(row,2).Value))
                if xlSheet1.Cells(row,3).Value == None or xlSheet1.Cells(row,3).Value == '':
                    InterHUImpList.append(1.0)
                else:
                    InterHUImpList.append(float(xlSheet1.Cells(row,3).Value))
                if xlSheet1.Cells(row,4).Value == None or xlSheet1.Cells(row,4).Value == '':
                    HUBuffDistList.append(0)
                else:
                    HUBuffDistList.append(int(xlSheet1.Cells(row,4).Value))
                counter += 1
            row += 1
        xlApp.ActiveWorkbook.Close(SaveChanges=0)
        xlApp.Quit()

        HU1Zip = zip(HUSSIDList, InterHUImpList, HUBuffDistList)
        HU1Zip.sort()
        HUSSIDList, InterHUImpList, HUBuffDistList = zip(*HU1Zip)

        gp.workspace = HU_Directory
        fcList = gp.ListFeatureClasses("*", "all")
        fc = fcList.Next()
        HULyrList = []
        HULyrAbbrevList = []
        HULyrBuffList = []
        HUIDList = []
        HUCount = 0
        while fc:
            # match SS ID with naming convention (_ID)
            checkInteger(fc)
            HULyrList.append(fc)
            fc0 = fc.replace(".", "")
            fc1 = fc0.replace("_", "")
            HULyrAbbrevList.append(fc1[:6]+str(HUCount + 1))
            HULyrBuffList.append(fc1[:8]+str(HUCount + 1)+"_buff.shp")
            fc2 = fc[::-1]
            j = fc2.find('_')
            indexS = len(fc)-j-1
            indexE = fc.find(".")
            fc_ID = fc[indexS+1:indexE]
            HUIDList.append(int(fc_ID))
            fc = fcList.Next()
            HUCount = HUCount + 1
        del fc

        HU2Zip = zip(HUIDList, HULyrList, HULyrAbbrevList, HULyrBuffList)
        HU2Zip.sort()
        HUIDList, HULyrList, HULyrAbbrevList, HULyrBuffList = zip(*HU2Zip)

        # process lists from table
        InterHUImpMax = max(InterHUImpList)
        InterHUNormList = []

        for i in range(0,len(InterHUImpList)):
            InterHUNormList.append(InterHUImpList[i]/InterHUImpMax) # Sj
    except:
        gp.AddError(msgPrepHULayers)
        raise Exception


    ################################################
    ###### CHECKS FOR HUMAN USES LAYER INPUTS ######
    ################################################
    try:
        # check that number of human uses layers is 18 or less
        if len(HULyrList) > 18:
            gp.AddError("\nThis model permits a maximum of 18 human uses layers.\nThe human uses layer directory you specified contains "+str(len(HULyrList))+" layers.")
            raise Exception

        if len(HULyrList) == 0:
            gp.AddError("\nThe human uses layer directory you specified contains no layers.")
            raise Exception

        if HUInputCount <> HUCount:
            gp.AddError("There is an inconsistency between the number of human uses layers in the specified directory and the input spreadsheet.")
            raise Exception
    except:
        gp.AddError(msgCheckHULayers)
        raise Exception


    ################################################
    #### BUFFER AND RASTERIZE HUMAN USES LAYERS ####
    ################################################
    try:    
        gp.workspace = interws
        gp.Extent = AZHU_rst
     
        for i in range(0,len(HULyrList)):
            HUVariable = HU_Directory+"\\"+HULyrList[i]
            # check each human uses layer's projection
            checkProjections(HUVariable)
            gp.MakeFeatureLayer_management(HUVariable, HULyrAbbrevList[i]+".lyr", "", "", "")
            SelectHULayer = gp.SelectLayerByLocation(AnalysisZonesHU_lyr, "INTERSECT", HULyrAbbrevList[i]+".lyr", "", "NEW_SELECTION")
            if gp.GetCount(SelectHULayer) > 0:
                pass
            else:
                gp.AddError(HULyrList[i]+" does not overlap the analysis grid created by the Grid the Seascape tool (GS).")
                raise Exception 
            
            if WghtFieldName:
                # check that "WghtFieldName" field exists
                FieldExists(HU_Directory+"\\"+HULyrList[i], WghtFieldName)
                if HUBuffDistList[i] > 0:
                    gp.Buffer_analysis(HUVariable, HULyrBuffList[i], str(HUBuffDistList[i]) + " Meters", "FULL", "ROUND", "NONE", "")        
                    gp.FeatureToRaster_conversion(HULyrBuffList[i], WghtFieldName, HULyrAbbrevList[i], "50")
                else:
                    gp.FeatureToRaster_conversion(HUVariable, WghtFieldName, HULyrAbbrevList[i], "50")

            else:
                HUVariable = AddField(HUVariable, "VID", "SHORT", "0", "0")        
                gp.CalculateField_management(HUVariable, "VID", 1, "VB")
                if HUBuffDistList[i] > 0:
                    gp.Buffer_analysis(HUVariable, HULyrBuffList[i], str(HUBuffDistList[i]) + " Meters", "FULL", "ROUND", "NONE", "")        
                    gp.FeatureToRaster_conversion(HULyrBuffList[i], "VID", HULyrAbbrevList[i], "50")
                else:
                    gp.FeatureToRaster_conversion(HUVariable, "VID", HULyrAbbrevList[i], "50")

        # clear selections from checks
        gp.SelectLayerByAttribute(AnalysisZonesHU_lyr, "CLEAR_SELECTION", "")

    except:
        gp.AddError(msgBuffRastHULayers)
        raise Exception    


    ################################
    ## DISTANCE DECAY FUNCTION #####
    ################################
    try:
        if PopPlaces:
            gp.AddMessage("\nMeasuring distances for decay function...")
            if AnalysisType == "Gridded Seascape (GS)":
                gp.RasterToPoint_conversion(AZHU_rst, AnalysisZonesHU_Pts, "VALUE")
            else:
                gp.FeatureToPoint_management(AnalysisZonesHU, AnalysisZonesHU_Pts, "CENTROID") # requires ArcInfo license
                
            AnalysisZonesHU = AddField(AnalysisZonesHU, "DIST_G2P", "DOUBLE", "0", "0")
            
            # convertes point data to x, y coordinates
            def GrabPts(xy, xyShp):
                cur = gp.UpdateCursor(xyShp)
                row = cur.Next()
                while row:
                    feat = row.Shape
                    midpoint1 = feat.Centroid
                    midList1 = shlex.split(midpoint1)
                    midList1 = [float(s) for s in midList1]
                    xy.append(midList1)
                    cur.UpdateRow(row)
                    row = cur.Next()
                del cur, row
                xy = np.array(xy)
                return xy

            # calculate distance between two numpy arrays
            def CalcDist(xy_1, xy_2):
                mindist=np.zeros(len(xy_1))
                minid=np.zeros(len(xy_1))
                for i,xy in enumerate(xy_1):
                    dists=np.sqrt(np.sum((xy-xy_2)**2,axis=1))
                    mindist[i],minid[i]=dists.min(),dists.argmin()
                return mindist, minid

            xy1 = []
            xy2 = []
            xy1 = GrabPts(xy1, AnalysisZonesHU_Pts) # projected analysis zones points
            xy2 = GrabPts(xy2, PopPlaces) # projected populated places
            G2P_Dist, G2P_ID = CalcDist(xy1,xy2) # analysis zones to pop places dist and ID

            cur = gp.UpdateCursor(AnalysisZonesHU, "", "", "*")
            row = cur.Next()
            j = 0
            while row:
                row.SetValue("DIST_G2P", G2P_Dist[j])
                cur.UpdateRow(row)
                row = cur.Next()
                j = j + 1
            del cur
            del row
    except:
        gp.AddError(msgDistanceDecay)
        raise Exception


    ##############################################     
    ######## CALCULATE IMPORTANCE SCORING ########
    ##############################################
    try:
        gp.AddMessage("\nConducting overlap analysis...")

        for i in range(0,HUInputCount):
            FieldNameArea = HULyrAbbrevList[i]+"_A"
            FieldNameMean = HULyrAbbrevList[i]+"_M"
            FieldNameIS = HULyrAbbrevList[i]+"_I"
            AnalysisZonesHU = AddField(AnalysisZonesHU, FieldNameArea, "DOUBLE", "8", "2")
            AnalysisZonesHU = AddField(AnalysisZonesHU, FieldNameMean, "DOUBLE", "8", "2")
            AnalysisZonesHU = AddField(AnalysisZonesHU, FieldNameIS, "DOUBLE", "8", "2")

        gp.snapRaster = AZHU_rst
        gp.cellsize = "MINOF"

        for i in range(0,len(HULyrAbbrevList)):
            gp.ZonalStatisticsAsTable_sa(AZHU_rst, "VALUE", HULyrAbbrevList[i], "zs_"+HULyrAbbrevList[i]+".dbf", "DATA")
            gp.AddJoin_management(AnalysisZonesHU_lyr, "VALUE", "zs_"+HULyrAbbrevList[i]+".dbf", "VALUE", "KEEP_COMMON")
            gp.CalculateField_management(AnalysisZonesHU_lyr, "AnalysisZonesHU."+HULyrAbbrevList[i]+"_A", "[zs_"+HULyrAbbrevList[i]+".AREA]", "VB", "")
            gp.CalculateField_management(AnalysisZonesHU_lyr, "AnalysisZonesHU."+HULyrAbbrevList[i]+"_M", "[zs_"+HULyrAbbrevList[i]+".MEAN]", "VB", "")
            gp.RemoveJoin_management(AnalysisZonesHU_lyr, "zs_"+HULyrAbbrevList[i])

        gp.FeatureClassToFeatureClass_conversion(AnalysisZonesHU_lyr, gp.workspace, "AnalysisZonesHU_area.shp", "")

        # max area of each cell
        Area_SqM = (cellsize*cellsize)
        # create list for max RIS across cells for each layer
        HUMaxIS = np.zeros(HUInputCount, dtype=np.float64)

        # grab max of importance score
        cur = gp.UpdateCursor(AnalysisZonesHU_area)
        row = cur.Next()
        while row:
            for i in range(0,len(HULyrAbbrevList)):
                if row.GetValue(HULyrAbbrevList[i]+"_M") > HUMaxIS[i]:
                    HUMaxIS[i] = row.GetValue(HULyrAbbrevList[i]+"_M") 
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur

        cur = gp.UpdateCursor(AnalysisZonesHU_area)
        row = cur.Next()
        while row:
            if AnalysisType == "Gridded Seascape (GS)":
                for i in range(0,HUInputCount):
                    FieldNameArea = row.GetValue(HULyrAbbrevList[i]+"_A")
                    FieldNameMean = row.GetValue(HULyrAbbrevList[i]+"_M")
                    row.SetValue(HULyrAbbrevList[i]+"_I", (FieldNameArea*FieldNameMean)/(Area_SqM*HUMaxIS[i]))
            else:
                for i in range(0,HUInputCount):
                    FieldNameMean = row.GetValue(HULyrAbbrevList[i]+"_M")
                    row.SetValue(HULyrAbbrevList[i]+"_I", (FieldNameMean/HUMaxIS[i]))
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur

        AnalysisZonesHU_area = AddField(AnalysisZonesHU_area, "FREQ", "SHORT", "0", "0")
        AnalysisZonesHU_area = AddField(AnalysisZonesHU_area, "IMP_SCORE", "DOUBLE", "8", "2")

        StatsList = []
        StatsArray = []
        cur = gp.UpdateCursor(AnalysisZonesHU_area)
        row = cur.Next()
        while row:
            count = 0
            ImpScoreExpr = 0.0
            for i in range(0,len(HULyrAbbrevList)):
                ImpValue = row.GetValue(HULyrAbbrevList[i]+"_I")
                if ImpValue > 1.0:
                    ImpValue = 1.0
                # for FREQ calc
                if ImpValue > 0:
                    count = count + 1
                # for IMP_SCORE calc        
                ImpScoreExpr = ImpScoreExpr + (ImpValue * InterHUNormList[i])

            row.SetValue("FREQ", count)
            if PopPlaces:
                MinDistance = row.GetValue("DIST_G2P")
                row.SetValue("IMP_SCORE", float(ImpScoreExpr*np.exp((-1)*float(DecayRate)*(MinDistance/1000.0))/HUInputCount))
            else:
                row.SetValue("IMP_SCORE", (ImpScoreExpr/HUInputCount))

            StatsList.append(row.GetValue("VALUE"))
            StatsList.append(row.GetValue("FREQ"))
            if WghtFieldName:
                StatsList.append(row.GetValue("IMP_SCORE"))
            StatsArray.append(StatsList)
            StatsList = []
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur
    except:
        gp.AddError(msgImportScoring)
        raise Exception

    try:
        # create raster outputs
        gp.FeatureToRaster_conversion(AnalysisZonesHU_area, "FREQ", hu_freq_all, cellsize)
        gp.Reclassify_sa(hu_freq_all, "VALUE", "1 1;0 NODATA", hu_freq, "DATA")
        if WghtFieldName:
            gp.FeatureToRaster_conversion(AnalysisZonesHU_area, "IMP_SCORE", hu_impscore, cellsize)

        calcs = open(HU_calcs_csv, "wb")
        writer = csv.writer(calcs, delimiter=',', quoting=csv.QUOTE_NONE)
        writer.writerow(['ZONE_ID', 'ACTIVITY_COUNT', 'IMPORTANCE_SCORE'])
        for i in range(0,len(StatsArray)):
            writer.writerow(StatsArray[i])
        calcs.close()
    except:
        gp.AddError(msgCreateOutputs)
        raise Exception    

    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(outputws+"parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("HUMAN USES OVERLAY ANALYSIS MODEL PARAMETERS\n")
    parafile.writelines("____________________________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

    # delete superfluous intermediate data
    del1 = [AnalysisZonesHU, AnalysisZonesHU_Pts, AnalysisZonesHU_lyr, AZHU_rst, hu_freq_all]
    deletelist = del1
    for data in deletelist:
        if gp.exists(data):
            gp.delete_management(data)

    gp.AddMessage("")
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())