# Marine InVEST: Wave Energy Model
# Authors: CK Kim, Gregg Verutes, Apollo Xi, Mike Papenfus
# Coded for ArcGIS 9.3 and 10
# 09/02/11

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
    try:
        # get parameters
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
        gp.workspace = gp.GetParameterAsText(0)
        parameters.append("Workspace: "+ gp.workspace)
        gp.scratchWorkspace = gp.GetParameterAsText(0)
        parameters.append("Scratch Workspace: "+ gp.scratchWorkspace)
        WaveDataFolder = gp.GetParameterAsText(1)
        parameters.append("Path to Folder with Wave Base Data: "+ WaveDataFolder)
        AnalysisArea = gp.GetParameterAsText(2)
        parameters.append("Wave Base Data"+ AnalysisArea)
        AOI = gp.GetParameterAsText(3)
        parameters.append("Area of Interest (AOI): "+ AOI)
        MachinePerform = gp.GetParameterAsText(4)
        parameters.append("Machine Performance Table: "+ MachinePerform)
        MachineParameter = gp.GetParameterAsText(5)
        parameters.append("Machine Parameters Table: "+ MachineParameter)
        DEM = gp.GetParameterAsText(6)
        parameters.append("Digital Elevation Model (DEM): "+ DEM)
        EconBoolean = gp.GetParameterAsText(7)
        parameters.append("Compute Economic Valuation? "+ EconBoolean)
        EconParameter = gp.GetParameterAsText(8)
        parameters.append("Machine Economic Parameters Table: "+ EconParameter)
        LandGridSheet = gp.GetParameterAsText(9)
        parameters.append("Landing and Grid Points: "+ LandGridSheet)
        NumUnits = gp.GetParameterAsText(10)
        parameters.append("Number of Machine Units: "+ str(NumUnits))
        if NumUnits:
            NumUnits = int(gp.GetParameterAsText(10))
    except:
        raise Exception, msgArguments + gp.GetMessages(2)

##################################################################################################################

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

    def checkDatum(thedata):
        desc = gp.describe(thedata)
        SR = desc.SpatialReference
        if SR.Type == "Geographic":
            strDatum = SR.DatumName         
        else:
            gp.OutputCoordinateSystem = SR
            strSR = str(gp.OutputCoordinateSystem)
            gp.OutputCoordinateSystem = ""
            n1 = strSR.find("DATUM[\'")
            n2 = strSR.find("\'",n1+7)
            strDatum = strSR[n1+7:n2]
        if strDatum == "D_WGS_1984":
            pass
        else:
            gp.AddError(thedata+" is not a valid input.\nThe model requires data inputs and a projection with the \"WGS84\" datum.\nSee InVEST FAQ document for how to reproject datasets.")
            raise Exception

    def checkProjections(thedata):
        dataDesc = gp.describe(thedata)
        spatreflc = dataDesc.SpatialReference
        if spatreflc.Type <> 'Projected':
            gp.AddError(thedata+" does not appear to be projected.  It is assumed to be in meters.")
        if spatreflc.LinearUnitName <> 'Meter':
            gp.AddError("This model assumes that "+thedata+" is projected in meters for area calculations.  You may get erroneous results.")

    def grabProjection(data):
        dataDesc = gp.describe(data)
        sr = dataDesc.SpatialReference
        gp.OutputCoordinateSystem = sr
        strSR = str(gp.OutputCoordinateSystem)
        return strSR

    # percentiles list (25, 50, 75, 90)
    def getPercentiles(list):
        PctList = []
        PctList.append(stats.scoreatpercentile(list, 25))
        PctList.append(stats.scoreatpercentile(list, 50))
        PctList.append(stats.scoreatpercentile(list, 75))
        PctList.append(stats.scoreatpercentile(list, 90))
        return PctList

    # checks of AOI and DEM
    if AOI:
        checkDatum(AOI)
        projection = grabProjection(AOI)
    checkDatum(DEM)

    # do not run econ analysis if an AOI is not specified
    if AOI == "" and EconBoolean == "true":
        EconBoolean = "false"
        gp.AddWarning("\nCannot conduct economic valuation without specifying an AOI.")
    
    # if conducting economic valuation, check that all econ data exists
    if EconBoolean == "true":
        inputs = [LandGridSheet, EconParameter, AOI]
        for x in inputs:
            if not gp.Exists(x):
                gp.AddError("\nOne or more of the required economic valuation input parameters was not defined.")
                raise Exception
        if NumUnits == "":
                gp.AddError("\nThe number of machine units is required to conduct economic valuation.")
                raise Exception

    # check and create folders
    try:
        thefolders=["intermediate", "Output"]
        for folder in thefolders:
            if not gp.exists(gp.workspace + folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        raise Exception, "Error creating folders"
    
    # local variables 
    try:        
        # intermediate and output directories
        outputws = gp.workspace + os.sep + "Output" + os.sep
        interws = gp.workspace + os.sep + "intermediate" + os.sep
        
        g = 9.81 # meter / second square
        d = 1028 # water density: kilogram / cubic meter
        alfa = 0.86 # wave period parameter

        # local analysis
        LandPtsTxt = interws + "LandPts.txt"
        LandPts_WGS84 = interws + "LandPts_WGS84.shp"
        GridPtTxt = interws + "GridPt.txt"
        GridPt_WGS84 = interws + "GridPt_WGS84.shp"
        WaveDataWCLyr = interws + "WaveDataWC.lyr"
        WaveDataECLyr = interws + "WaveDataEC.lyr"
        AOILyr = interws + "AOILyr.lyr"
        PointsAOICount = interws + "PointsAOICount.shp"
        DEM_prj = interws + "DEM_prj"
        WaveDataWC = WaveDataFolder + os.sep + "NAmerica_WestCoast_4m.shp"
        WaveDataEC = WaveDataFolder + os.sep + "NAmerica_EastCoast_4m.shp"
        WCNA_barrier = WaveDataFolder + os.sep + "WCNA_barrier.shp"
        WCNA_extract = WaveDataFolder + os.sep + "WCNA_extract.shp"
        ECNA_barrier = WaveDataFolder + os.sep + "ECNA_barrier.shp"
        ECNA_extract = WaveDataFolder + os.sep + "ECNA_extract.shp"
        WaveData_clip = interws + "WaveData_clip.shp"
        WaveData_clipZ = interws + "WaveData_clipZ.shp"
        WaveData_prj = interws + "WEM_InputOutput_Pts.shp"

        # global analysis       
        GlobalBox_EastHemi = WaveDataFolder + os.sep + "GlobalBox_EastHemi.shp"
        GlobalBox_WestHemi = WaveDataFolder + os.sep + "GlobalBox_WestHemi.shp"
        WaveDataEHLyr = interws + "WaveDataEHLyr.lyr"
        WaveDataWHLyr = interws + "WaveDataWHLyr.lyr"
        WaveDataEH = WaveDataFolder + os.sep + "Global_EastHemi_30m.shp"
        WaveDataWH = WaveDataFolder + os.sep + "Global_WestHemi_30m.shp"
        WH_barrier = WaveDataFolder + os.sep + "WH_barrier.shp"
        WH_extract = WaveDataFolder + os.sep + "WH_extract.shp"
        EH_barrier = WaveDataFolder + os.sep + "EH_barrier.shp"
        EH_extract = WaveDataFolder + os.sep + "EH_extract.shp"
        WaveData_clipEH = interws + "WaveData_clipEH.shp"
        WaveData_clipWH = interws + "WaveData_clipWH.shp"
        WaveData_clipZ_EH = interws + "WaveData_clipZ_EH.shp"
        WaveData_clipZ_WH = interws + "WaveData_clipZ_WH.shp"

        # variables for both global and regional/local analysis
        splineWP = interws + "splineWP"
        splineCWE = interws + "splineCWE"
        neighbIntp = interws + "neighbIntp"
        splineWPext = interws + "splineWPext"    
        splineCWEext = interws + "splineCWEext"
        splineWPint = interws + "splineWPint"
        splineCWEint = interws + "splineCWEint"
        splineNPV = interws + "splineNPV"
        neighbIntp2 = interws + "neighbIntp2"
        splineNPVext = interws + "splineNPVext"

        # variables for output        
        outputWP = outputws + "wp_kw"
        outputWP_rc = outputws + "wp_rc"
        outputCWE = outputws + "capwe_mwh"
        outputCWE_rc = outputws + "capwe_rc"
        outputNPV = outputws + "npv_usd"
        outputNPV_rc = outputws + "npv_rc"
        GridPt_prj = outputws + "GridPt_prj.shp"
        LandPts_prj = outputws + "LandPts_prj.shp"
    
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages())
        raise Exception

##############################################################################################################################

    ######################################################
    ############# WAVE POWER AND CWE FUNCTIONS ###########
    ######################################################
    
    # func1 calculates the wave numbers
    def func1(sigma, h):
        kestimated = (sigma**2)/(g*(sqrt(tanh((sigma**2)*h/g))))
        kprevious = 0.0000001
        count = 0
        while (abs(kestimated-kprevious) > 0.000005) and (count < 1000):
            count += 1
            kh = kestimated*h
            kcalculated = (sigma**2)/(tanh(kh)*g)
            kprevious = kestimated
            kestimated = kcalculated
        k = kcalculated
        return k

    # func2 calculates wave group velocity
    def func2(k2, h2):
        if 2*k2*h2 <= 300:
            cg = 0.5*(1+2*k2*h2/sinh(2*k2*h2))*sqrt(g*tanh(k2*h2)/k2)
        else:
            cg = 0.5*sqrt(g*tanh(k2*h2)/k2)
        return cg

    # func3 calculates wave energy in kilowatts
    def func3(hs, c):
        p = (d*g/16)*(hs**2)*c/1000
        return p
    
    # calculate WE_kWM
    def WPcalc(WaveData_clipZ):
        cur = gp.UpdateCursor(WaveData_clipZ, "", "", "HSAVG_M; TPAVG_S; DEPTH_M; WE_kWM")
        row = cur.Next()
        while row:
            if row.DEPTH_M > -9999 and row.DEPTH_M < 0:
                tem = 2.0*pi/(row.TPAVG_S*float(alfa))
                kest = func1(tem, abs(row.DEPTH_M))   # function 1
                cest = func2(kest, abs(row.DEPTH_M))  # function 2
                pest = func3(row.HSAVG_M, cest)      # function 3
                row.SetValue("WE_kWM", pest)
            else:
                row.SetValue("WE_kWM", 0)
            cur.UpdateRow(row)
            row = cur.Next()
        del cur, row
        return WaveData_clipZ

    # calculate captured wave energy    
    def CWEcalc(WaveData_clipZ, SeastateTxt, index_Tp, index_Hs, x_array, y_array, z_array, PointCount):
        count = 0
        # run through wave data and populate array with point index values: I, J
        PtStrList = []
        PtTxtList = []
        FIDList = []
        CapWEList = [0]*PointCount
        cur = gp.UpdateCursor(WaveData_clipZ, "", "", "I; J")
        row = cur.Next()
        while row:
            valueI = int(row.GetValue("I"))
            valueJ = int(row.GetValue("J"))
            PtStrList.append("I = " +str(row.GetValue("I"))+" AND J = "+str(row.GetValue("J")))
            PtTxtString = 'I,'
            if valueI < 10:
                PtTxtString = PtTxtString + "  "+str(valueI)
            elif valueI >= 10 and valueI < 100:
                PtTxtString = PtTxtString + " "+str(valueI)
            else:
                PtTxtString = PtTxtString + str(valueI)
            PtTxtString = PtTxtString + ',J,'   
            if valueJ < 10:
                PtTxtString = PtTxtString + "  "+str(valueJ)
            elif valueJ >= 10 and valueJ < 100:
                PtTxtString = PtTxtString + " "+str(valueJ)
            else:
                PtTxtString = PtTxtString + str(valueJ)
            PtTxtList.append(PtTxtString+"\n")
            FIDList.append(int(row.GetValue("FID")))
            cur.UpdateRow(row)
            row = cur.Next()
        del cur    
        del row

        # reads in seastate data from WW3
        SSTables = open(SeastateTxt,"r")
        # read in x and y lists only once
        listID = SSTables.readline()
        arrayX = SSTables.readline()
        arrayX = [float(s) for s in arrayX.split(",")] 
        arrayY = SSTables.readline()
        arrayY = [float(s) for s in arrayY.split(",")]
        del listID
        SSTables.close()

        counter = 25
        for line in open(SeastateTxt,"r"):
            if line in PtTxtList:
                indexed = PtTxtList.index(line)
                counter = 0
                arrayZ = []
            if counter > 2 and counter < 24:
                temp_line = line.split(",")
                arrayZ.append(temp_line)
                
           # increase counter by 1
            counter += 1
            
            if counter == 24:
                arrayZ = np.array(arrayZ, dtype='f')

                # set parameter max limits for each seastate table
                for row in range(0,21):
                    for col in range(0,21):
                        if col >= index_Tp:
                            arrayZ[row][col] = 0.0
                        if row >= index_Hs:
                            arrayZ[row][col] = 0.0

                # divide Z data by 5 to get yearly average
                arrayZ = np.divide(arrayZ, 5.0)

                # interpolation and calculate cap wave energy
                ip = interp2d(x_array, y_array, z_array, kind = 'cubic', copy = True, bounds_error = False, fill_value = 0.0)
                z_Range_intp = ip(arrayX, arrayY)
                z_Array_Intp = np.array(z_Range_intp)
                capwave_Array = arrayZ*z_Array_Intp
                capwave_Array = np.where(capwave_Array < 0, 0, capwave_Array)

                # assign CWE value to 'CapWEList' 
                CapWEList[indexed] = (capwave_Array.sum()/1000)
                counter = 25              

                # benchmark
                count += 1
                if int(PointCount*0.25) == count:
                    gp.AddMessage("...25% completed")
                elif int(PointCount*0.50) == count:
                    gp.AddMessage("...50% completed")
                elif int(PointCount*0.75) == count:
                    gp.AddMessage("...75% completed")

        # sort lists by FID
        WEPtsZip = zip(FIDList, CapWEList, PtStrList, PtTxtList)
        WEPtsZip.sort()
        FIDList, CapWEList, PtStrList, PtTxtList = zip(*WEPtsZip)

        cur = gp.UpdateCursor(WaveData_clipZ)
        row = cur.Next()
        m = 0
        while row:
            row.SetValue("CAPWE_MWHY", CapWEList[m])
            cur.UpdateRow(row)
            row = cur.next()
            m += 1
        del cur, row
        del arrayZ
        return WaveData_clipZ

##############################################################################################################################

    ##############################################
    ################ PREPARE DATA ################
    ##############################################

    # set extent to AOI and projection to that of AOI
    if AOI:
        gp.Extent = AOI
        gp.OutputCoordinateSystem = projection
    else:
        gp.Extent = "MAXOF"

    # clip the wave data to drawing boundary or global hemispheres, grab depth values for each point, add fields
    try:
        gp.AddMessage ("\nPreparing input data..." )
        if AnalysisArea == "West Coast of North America and Hawaii":
            if AOI:
                gp.MakeFeatureLayer_management(WaveDataWC, WaveDataWCLyr, "", "", "")
                gp.MakeFeatureLayer_management(AOI, AOILyr, "", "", "")
                PointsSelectWC = gp.SelectLayerByLocation_management(WaveDataWCLyr, "INTERSECT", AOILyr, "", "NEW_SELECTION")
                gp.FeatureClassToFeatureClass_conversion(PointsSelectWC, interws, "PointsAOICount.shp", "")
                if gp.GetCount_management(PointsAOICount) == 0:
                    gp.AddError(msgNoFeatures)
                    raise Exception
                else:
                    gp.Clip_analysis(WaveDataWC, AOI, WaveData_clip)
            else:
                gp.Clip_analysis(WaveDataWC, WCNA_extract, WaveData_clip)
            SeastateTxt = WaveDataFolder + os.sep + "NAmerica_WestCoast_4m.txt"   

        if AnalysisArea == "East Coast of North America and Puerto Rico":
            if AOI:
                gp.MakeFeatureLayer_management(WaveDataEC, WaveDataECLyr, "", "", "")
                gp.MakeFeatureLayer_management(AOI, AOILyr, "", "", "")
                PointsSelectEC = gp.SelectLayerByLocation_management(WaveDataECLyr, "INTERSECT", AOILyr, "", "NEW_SELECTION")
                gp.FeatureClassToFeatureClass_conversion(PointsSelectEC, interws, "PointsAOICount.shp", "")
                if gp.GetCount_management(PointsAOICount) == 0:
                    gp.AddError(msgNoFeatures)
                    raise Exception
                else:
                    gp.Clip_analysis(WaveDataEC, AOI, WaveData_clip)
            else:
                gp.Clip_analysis(WaveDataEC, ECNA_extract, WaveData_clip)
            SeastateTxt = WaveDataFolder + os.sep + "NAmerica_EastCoast_4m.txt"

        if AnalysisArea == "Global (Eastern Hemisphere)":
            if AOI:
                gp.MakeFeatureLayer_management(WaveDataEH, WaveDataEHLyr, "", "", "")
                gp.MakeFeatureLayer_management(AOI, AOILyr, "", "", "")
                PointsSelectEH = gp.SelectLayerByLocation_management(WaveDataEHLyr, "INTERSECT", AOILyr, "", "NEW_SELECTION")
                gp.FeatureClassToFeatureClass_conversion(PointsSelectEH, interws, "PointsAOICount.shp", "")
                if gp.GetCount_management(PointsAOICount) == 0:
                    gp.AddError(msgNoFeatures)
                    raise Exception
                else:
                    gp.Clip_analysis(WaveDataEH, AOI, WaveData_clip)
            else:
                gp.Clip_analysis(WaveDataEH, GlobalBox_EastHemi, WaveData_clip)
            SeastateTxt = WaveDataFolder + os.sep + "Global_EastHemi_30m.txt"

        if AnalysisArea == "Global (Western Hemisphere)":
            if AOI:
                gp.MakeFeatureLayer_management(WaveDataWH, WaveDataWHLyr, "", "", "")
                gp.MakeFeatureLayer_management(AOI, AOILyr, "", "", "")
                PointsSelectWH = gp.SelectLayerByLocation_management(WaveDataWHLyr, "INTERSECT", AOILyr, "", "NEW_SELECTION")
                gp.FeatureClassToFeatureClass_conversion(PointsSelectWH, interws, "PointsAOICount.shp", "")
                if gp.GetCount_management(PointsAOICount) == 0:
                    gp.AddError(msgNoFeatures)
                    raise Exception
                else:
                    gp.Clip_analysis(WaveDataWH, AOI, WaveData_clip)
            else:
                gp.Clip_analysis(WaveDataWH, GlobalBox_WestHemi, WaveData_clip)
            SeastateTxt = WaveDataFolder + os.sep + "Global_WestHemi_30m.txt"
    except:
        raise Exception, msgClipWaveData


    try:
        # get number of points in clipped wave data
        PointCount = gp.GetCount_management(WaveData_clip)
        # project DEM, if not projected
        if AOI:
            gp.CopyRaster_management(DEM, DEM_prj, "", "", "", "", "")
            # add depth values to points
            gp.ExtractValuesToPoints_sa(WaveData_clip, DEM_prj, WaveData_clipZ, "INTERPOLATE")
        else:
            gp.ExtractValuesToPoints_sa(WaveData_clip, DEM, WaveData_clipZ, "INTERPOLATE") 
        # add depth field
        WaveData_clipZ = AddField(WaveData_clipZ, "DEPTH_M", "DOUBLE", "0", "0")
        gp.CalculateField_management(WaveData_clipZ, "DEPTH_M", "[RASTERVALU]", "VB")
        gp.DeleteField_management(WaveData_clipZ, "RASTERVALU")
        # add WE fields
        WaveData_clipZ = AddField(WaveData_clipZ, "WE_kWM", "DOUBLE", "8", "2")  
        WaveData_clipZ = AddField(WaveData_clipZ, "CAPWE_MWHY", "DOUBLE", "8", "2")
            
    except:
        raise Exception, msgCalcDepth


    #################################################################################
    ################ READ IN PERFORMANCE AND PARAMETER INFO #########################
    #################################################################################

    try:
        # read in machine performance array
        xlApp = Dispatch("Excel.Application")
        xlApp.Visible = 0
        xlApp.DisplayAlerts=0
        xlBook = xlApp.Workbooks.Open(MachinePerform[:-(1+len(MachinePerform.split("\\")[-1]))])
        machinepath = MachinePerform.split("\\")
        machinesheet = machinepath[-1]
        xlSheet = xlBook.Worksheets(machinesheet[:-1])
        
        row = 1
        col = 1
        bottom = row
        while xlSheet.Cells(bottom+1, col).Value not in [None, '']:
            bottom += 1
        right = col
        while xlSheet.Cells(row, right+1).Value not in [None, '']:
            right += 1
            
        x_Range = xlSheet.Range(xlSheet.Cells(1,2), xlSheet.Cells(1,right)).Value
        y_Range = xlSheet.Range(xlSheet.Cells(2,1), xlSheet.Cells(bottom,1)).Value
        z_Range = xlSheet.Range(xlSheet.Cells(2,2), xlSheet.Cells(bottom,right)).Value
        
        x_array_temp = list(x_Range)
        x_array = []
        for i in x_array_temp[0]:
            x_array.append(i)
        
        y_array_temp = list(y_Range)
        y_array = []
        for i in y_array_temp:
            y_array_temp2 = list(i)
            y_array.append(y_array_temp2[0])
            
        z_array_temp = list(z_Range)
        z_array = []
        for i in z_array_temp:
            z_array_temp2 = list(i)
            z_array.append(z_array_temp2)

        # read in machine parameter limits
        x2Book = xlApp.Workbooks.Open(MachineParameter[:-(1+len(MachineParameter.split("\\")[-1]))])
        parameter = MachineParameter.split("\\")
        parametersheet = parameter[-1]
        x2Sheet = x2Book.Worksheets(parametersheet[:-1])
           
        Hs_max = x2Sheet.Range("b3").Value
        Tp_max = x2Sheet.Range("b4").Value

        # close the Excel applications
        xlApp.ActiveWorkbook.Close(SaveChanges=0)
        xlApp.Quit()

        # WW3 set matrix
        Tp = [0.25, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 59.75]
        Hs = [0.13, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 54.88]
        
        for i in range(0,len(Tp)):
            if Tp[i] > Tp_max:
                index_Tp = Tp.index(Tp[i])
                break
        for j in range(0,len(Hs)):
            if Hs[j] > Hs_max:
                index_Hs = Hs.index(Hs[j])
                break
    except:
        raise Exception, msgMachineTables


    ###################################################################
    ##################### WE EXPRESSIONS #############################
    ###################################################################
    try:
        if AnalysisArea == "West Coast of North America and Hawaii" or AnalysisArea == "East Coast of North America and Puerto Rico":
            gp.AddMessage("\nPerforming wave energy calculations...")
            WaveData_clipZ = WPcalc(WaveData_clipZ)
            WaveData_clipZ = CWEcalc(WaveData_clipZ, SeastateTxt, index_Tp, index_Hs, x_array, y_array, z_array, PointCount)

            # logic whether to project wave data points
            if AOI:
                gp.AddMessage("...projecting wave data within AOI")
                gp.Project_management(WaveData_clipZ, WaveData_prj, projection, "")
                gp.AddMessage("...generating wave power and captured wave energy outputs\n")
                gp.Spline_sa(WaveData_prj, "WE_kWM", splineWP)
                gp.Spline_sa(WaveData_prj, "CAPWE_MWHY", splineCWE)
                gp.NaturalNeighbor_sa(WaveData_prj, "CAPWE_MWHY", neighbIntp)
                # output: wave power
                gp.ExtractByMask_sa(splineWP, neighbIntp, splineWPext)
                gp.Int_sa(splineWPext, splineWPint)
                gp.SetNull_sa(splineWPint, splineWPint, outputWP, '"VALUE" < 1')
                # output: captured wave energy
                gp.ExtractByMask_sa(splineCWE, neighbIntp, splineCWEext)
                gp.Int_sa(splineCWEext, splineCWEint)
                gp.SetNull_sa(splineCWEint, splineCWEint, outputCWE, '"VALUE" < 1')
                
            else:
                gp.AddMessage("...generating wave power and captured wave energy outputs\n")
                if AnalysisArea == "West Coast of North America and Hawaii":
                    gp.SplineWithBarriers_sa(WaveData_clipZ, "WE_kWM", WCNA_barrier, "0.05", splineWP, "0")
                    gp.SplineWithBarriers_sa(WaveData_clipZ, "CAPWE_MWHY", WCNA_barrier, "0.05", splineCWE, "0")
                    gp.ExtractByMask_sa(splineWP, WCNA_extract, splineWPext)
                    gp.Int_sa(splineWPext, splineWPint)
                    gp.SetNull_sa(splineWPint, splineWPint, outputWP, '"VALUE" < 1')                
                    gp.ExtractByMask_sa(splineCWE, WCNA_extract, splineCWEext)
                    gp.Int_sa(splineCWEext, splineCWEint)
                    gp.SetNull_sa(splineCWEint, splineCWEint, outputCWE, '"VALUE" < 1') 
                if AnalysisArea == "East Coast of North America and Puerto Rico":
                    gp.SplineWithBarriers_sa(WaveData_clipZ, "WE_kWM", ECNA_barrier, "0.05", splineWP, "0")
                    gp.SplineWithBarriers_sa(WaveData_clipZ, "CAPWE_MWHY", ECNA_barrier, "0.05", splineCWE, "0")
                    gp.ExtractByMask_sa(splineWP, ECNA_extract, splineWPext)
                    gp.Int_sa(splineWPext, splineWPint)
                    gp.SetNull_sa(splineWPint, splineWPint, outputWP, '"VALUE" < 1')                
                    gp.ExtractByMask_sa(splineCWE, ECNA_extract, splineCWEext)
                    gp.Int_sa(splineCWEext, splineCWEint)
                    gp.SetNull_sa(splineCWEint, splineCWEint, outputCWE, '"VALUE" < 1') 

        #######################################################################################################

        if AnalysisArea == "Global (Eastern Hemisphere)" or AnalysisArea == "Global (Western Hemisphere)":
            gp.AddMessage("\nPerforming wave energy calculations...")
            WaveData_clipZ = WPcalc(WaveData_clipZ)
            WaveData_clipZ = CWEcalc(WaveData_clipZ, SeastateTxt, index_Tp, index_Hs, x_array, y_array, z_array, PointCount)

            if AOI:
                gp.AddMessage("...projecting wave data within AOI")
                gp.Project_management(WaveData_clipZ, WaveData_prj, projection, "")
                gp.AddMessage("...generating wave power and captured wave energy outputs\n")
                gp.Spline_sa(WaveData_prj, "WE_kWM", splineWP)
                gp.Spline_sa(WaveData_prj, "CAPWE_MWHY", splineCWE)
                gp.NaturalNeighbor_sa(WaveData_prj, "CAPWE_MWHY", neighbIntp)
                # output: wave power
                gp.ExtractByMask_sa(splineWP, neighbIntp, splineWPext)
                gp.Int_sa(splineWPext, splineWPint)
                gp.SetNull_sa(splineWPint, splineWPint, outputWP, '"VALUE" < 1')
                # output: captured wave energy
                gp.ExtractByMask_sa(splineCWE, neighbIntp, splineCWEext)
                gp.Int_sa(splineCWEext, splineCWEint)
                gp.SetNull_sa(splineCWEint, splineCWEint, outputCWE, '"VALUE" < 1')                
                
            else:
                gp.AddMessage("...generating wave power and captured wave energy outputs\n")
                if AnalysisArea == "Global (Eastern Hemisphere)":
                    gp.SplineWithBarriers_sa(WaveData_clipZ, "WE_kWM", EH_barrier, "0.05", splineWP, "0")
                    gp.SplineWithBarriers_sa(WaveData_clipZ, "CAPWE_MWHY", EH_barrier, "0.05", splineCWE, "0")
                    gp.ExtractByMask_sa(splineWP, EH_extract, splineWPext)
                else:
                    gp.SplineWithBarriers_sa(WaveData_clipZ, "WE_kWM", WH_barrier, "0.05", splineWP, "0")
                    gp.SplineWithBarriers_sa(WaveData_clipZ, "CAPWE_MWHY", WH_barrier, "0.05", splineCWE, "0")
                    gp.ExtractByMask_sa(splineWP, WH_extract, splineWPext)
                gp.Int_sa(splineWPext, splineWPint)
                gp.SetNull_sa(splineWPint, splineWPint, outputWP, '"VALUE" < 1')
                if AnalysisArea == "Global (Eastern Hemisphere)":
                    gp.ExtractByMask_sa(splineCWE, EH_extract, splineCWEext)
                else:
                    gp.ExtractByMask_sa(splineCWE, WH_extract, splineCWEext)
                gp.Int_sa(splineCWEext, splineCWEint)
                gp.SetNull_sa(splineCWEint, splineCWEint, outputCWE, '"VALUE" < 1')

        # reclass WP & CWE rasters
        WPList = []
        cur = gp.UpdateCursor(outputWP)
        row = cur.Next()
        while row:
            CellCount = int(row.GetValue("COUNT"))
            for i in range(CellCount):
                WPList.append(row.GetValue("VALUE"))
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur

        CWEList = []
        cur = gp.UpdateCursor(outputCWE)
        row = cur.Next()
        while row:
            CellCount = int(row.GetValue("COUNT"))
            for i in range(CellCount):
                CWEList.append(row.GetValue("VALUE"))
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur            
        
        WPPctList = getPercentiles(WPList)
        CWEPctList = getPercentiles(CWEList)
        WPExpr = "1 "+str(int(WPPctList[0]))+" 1;"+str(int(WPPctList[0]))+" "+str(int(WPPctList[1]))+" 2;"+str(int(WPPctList[1]))+" "\
                 +str(int(WPPctList[2]))+" 3;"+str(int(WPPctList[2]))+" "+str(int(WPPctList[3]))+" 4;"+str(int(WPPctList[3]))+" "\
                 +str(int(max(WPList)))+" 5"
        CWEExpr = "1 "+str(int(CWEPctList[0]))+" 1;"+str(int(CWEPctList[0]))+" "+str(int(CWEPctList[1]))+" 2;"+str(int(CWEPctList[1]))+" "\
                 +str(int(CWEPctList[2]))+" 3;"+str(int(CWEPctList[2]))+" "+str(int(CWEPctList[3]))+" 4;"+str(int(CWEPctList[3]))+" "\
                 +str(int(max(CWEList)))+" 5"

        # reclassify wave power outputs
        gp.Reclassify_sa(outputWP, "VALUE", WPExpr, outputWP_rc, "DATA")
        outputWP_rc = AddField(outputWP_rc, "VAL_RANGE", "TEXT", "75", "")
        cur = gp.UpdateCursor(outputWP_rc)
        row = cur.Next()
        while row:
            WPValue = row.GetValue("VALUE")
            if WPValue == 1:
                row.SetValue("VAL_RANGE", "1 - "+str(WPPctList[0])+" kilowatts per square meter (kW/m)")
            elif WPValue == 2:
                row.SetValue("VAL_RANGE", str(WPPctList[0])+" - "+str(WPPctList[1])+" kW/m")
            elif WPValue == 3:
                row.SetValue("VAL_RANGE", str(WPPctList[1])+" - "+str(WPPctList[2])+" kW/m")
            elif WPValue == 4:
                row.SetValue("VAL_RANGE", str(WPPctList[2])+" - "+str(WPPctList[3])+" kW/m")
            elif WPValue == 5:
                row.SetValue("VAL_RANGE", "Greater than "+str(WPPctList[3])+" kW/m")
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur
        
        # reclassify captured wave energy outputs
        gp.Reclassify_sa(outputCWE, "VALUE", CWEExpr, outputCWE_rc, "DATA")
        outputCWE_rc = AddField(outputCWE_rc, "VAL_RANGE", "TEXT", "75", "")
        cur = gp.UpdateCursor(outputCWE_rc)
        row = cur.Next()
        while row:
            CWEValue = row.GetValue("VALUE")
            if CWEValue == 1:
                row.SetValue("VAL_RANGE", "1 - "+str(CWEPctList[0])+" megawatt hours per year (MWh/yr)")
            elif CWEValue == 2:
                row.SetValue("VAL_RANGE", str(CWEPctList[0])+" - "+str(CWEPctList[1])+" (MWh/yr)")
            elif CWEValue == 3:
                row.SetValue("VAL_RANGE", str(CWEPctList[1])+" - "+str(CWEPctList[2])+" (MWh/yr)")
            elif CWEValue == 4:
                row.SetValue("VAL_RANGE", str(CWEPctList[2])+" - "+str(CWEPctList[3])+" (MWh/yr)")
            elif CWEValue == 5:
                row.SetValue("VAL_RANGE", "Greater than "+str(CWEPctList[3])+" (MWh/yr)")
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur
    except:
        raise Exception, msgWEcalc
    
    
    ##############################################################
    ####################### ECONOMIC ANALYSIS ####################
    ##############################################################
    try:
        if EconBoolean == "true":
            gp.AddMessage("Performing economic valuation..." )
            # create a text file and write point coordinates to it
            land = open(LandPtsTxt,'a')
            grid = open(GridPtTxt,'a')
            thestring = "Point\n"
            land.writelines(thestring)
            grid.writelines(thestring)
            counter = 2

            xlApp = Dispatch("Excel.Application")
            xlApp.Visible = 0
            xlApp.DisplayAlerts=0
            xlBook1 = xlApp.Workbooks.Open(LandGridSheet[:-(1+len(LandGridSheet.split("\\")[-1]))])
            WECpath = LandGridSheet.split("\\")
            WECsheet = WECpath[-1]
            xlSheet = xlBook1.Worksheets(WECsheet[:-1])

            row = 1
            col = 2
            bottom = row
            while xlSheet.Cells(bottom+1, col).Value not in [None, '']:
                bottom += 1
            while counter <= bottom:
                xl = xlSheet.Cells(counter,2).Value
                yl = xlSheet.Cells(counter,3).Value
                xg = xlSheet.Cells(counter,5).Value
                yg = xlSheet.Cells(counter,6).Value
                landstring = str(counter-1)+" "+ str(yl)+" "+str(xl)+"\n"
                land.writelines(landstring)
                gridstring = str(counter-1)+" "+ str(yg)+" "+str(xg)+"\n"
                grid.writelines(gridstring)
                counter = counter + 1
            thestring = "END"
            land.writelines(thestring)
            grid.writelines(thestring)
            land.close()
            grid.close()

            # Process: Create Features From Text File...
            gp.CreateFeaturesFromTextFile_samples(LandPtsTxt, ".", LandPts_WGS84, "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];IsHighPrecision")
            gp.CreateFeaturesFromTextFile_samples(GridPtTxt, ".", GridPt_WGS84, "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];IsHighPrecision")

            # check that there is at least one grid and one landing point input
            if gp.GetCount_management(GridPt_WGS84) < 1:
                gp.AddError("Model takes only one grid point location as input.")
                raise Exception
            else:
                gp.Project_management(GridPt_WGS84, GridPt_prj, projection, "")
                
            if gp.GetCount_management(LandPts_WGS84) < 1:
                gp.AddError("Model must have at least one landing point location as input.")
                raise Exception
            else:
                gp.Project_management(LandPts_WGS84, LandPts_prj, projection, "")

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
            xy3 = []
            xy1 = GrabPts(xy1, WaveData_prj) # projected wave data
            xy2 = GrabPts(xy2, LandPts_prj) # projected land
            xy3 = GrabPts(xy3, GridPt_prj) # projected grid

            W2L_Dist, W2L_ID = CalcDist(xy1,xy2) # wave to land dist and ID
            L2G_Dist, L2G_ID = CalcDist(xy2,xy3) # land to grid dist and ID

            # add various economic fields...
            WaveData_prj = AddField(WaveData_prj, "W2L_MDIST", "DOUBLE", "8", "2")  
            WaveData_prj = AddField(WaveData_prj, "LAND_ID", "SHORT", "0", "0")
            WaveData_prj = AddField(WaveData_prj, "L2G_MDIST", "DOUBLE", "8", "2")
            WaveData_prj = AddField(WaveData_prj, "UNITS", "SHORT", "0", "0")
            WaveData_prj = AddField(WaveData_prj, "CAPWE_ALL", "DOUBLE", "8", "2")
            WaveData_prj = AddField(WaveData_prj, "NPV_25Y", "DOUBLE", "8", "2")
            gp.CalculateField_management(WaveData_prj, "UNITS", NumUnits, "VB")

            # populate projected wave data
            cur = gp.UpdateCursor(WaveData_prj, "", "", "W2L_MDIST; LAND_ID; L2G_MDIST")
            row = cur.Next()
            j = 0
            while row:
                row.SetValue("W2L_MDIST", W2L_Dist[j])
                LANDID = int(W2L_ID[j])
                row.SetValue("LAND_ID", LANDID)
                row.SetValue("L2G_MDIST", L2G_Dist[LANDID])
                cur.UpdateRow(row)
                row = cur.Next()
                j = j + 1
            del cur, row

            # calculate net present value for each wave farm site            
            cur = gp.UpdateCursor(WaveData_prj, "", "", "DEPTH_M; CAPWE_MWHY; UNITS; CAPWE_ALL; NPV_25Y; W2L_MDIST; LAND_ID; L2G_MDIST")
            row = cur.Next()
            while row:
                capturedWE = np.ones(len(T))*int(row.CAPWE_MWHY)*1000.0  # kWh here
                capturedWE[0] = 0.0
        
                # initial costs
                lenml = 3.0*abs(row.DEPTH_M)
                q = row.UNITS
                installCost = q*capRate*cc
                mooringCost = smlpm * lenml * cml * q
                transCost = ((row.W2L_MDIST)*cul/1000.0) + ((row.L2G_MDIST)*col/1000.0)
                IC = installCost + mooringCost + transCost
        
                # annual revenue and costs
                annualRevenue = p*q*capturedWE
                annualCost = omc*capturedWE*q
                annualCost[0] = IC
                
                row.SetValue("UNITS", q)
                row.SetValue("CAPWE_ALL", (row.CAPWE_MWHY * q))
                row.SetValue("NPV_25Y", (NPVWE(annualRevenue, annualCost)/1000.0)) # convert NPV to value "in thousands"
                cur.UpdateRow(row)
                row = cur.Next()
                
            xlApp.ActiveWorkbook.Close(SaveChanges=0)
            xlApp.Quit()
            del xlApp
            del cur, row
        
            #convert point NPV to raster surface
            gp.AddMessage("...generating valuation outputs\n")
            gp.Spline_sa(WaveData_prj, "NPV_25Y", splineNPV)
            gp.NaturalNeighbor_sa(WaveData_prj, "NPV_25Y", neighbIntp2)
            gp.ExtractByMask_sa(splineNPV, neighbIntp2, splineNPVext)
            gp.Int_sa(splineNPVext, outputNPV)

            # reclass NPV raster
            gp.BuildRasterAttributeTable_management(outputNPV, "Overwrite")
            outputNPV = AddField(outputNPV, "NPV_BREAKS", "SHORT", "", "")            
            NPVList = []
            cur = gp.UpdateCursor(outputNPV)
            row = cur.Next()
            while row:
                CellCount = int(row.GetValue("COUNT"))
                if row.GetValue("VALUE") > 0:
                    for i in range(CellCount):
                        NPVList.append(row.GetValue("VALUE"))
                cur.UpdateRow(row)
                row = cur.next()
            del row
            del cur
            
            NPVPctList = getPercentiles(NPVList)
            cur = gp.UpdateCursor(outputNPV)
            row = cur.Next()
            while row:
                NPVValue = row.GetValue("VALUE")
                if NPVValue < 1:
                    row.SetValue("NPV_BREAKS", 0)
                elif NPVValue < NPVPctList[0] and NPVValue >= 1:
                    row.SetValue("NPV_BREAKS", 1)
                elif NPVValue >= NPVPctList[0] and NPVValue < NPVPctList[1]:
                    row.SetValue("NPV_BREAKS", 2)
                elif NPVValue >= NPVPctList[1] and NPVValue <= NPVPctList[2]:
                    row.SetValue("NPV_BREAKS", 3)
                elif NPVValue >= NPVPctList[2] and NPVValue <= NPVPctList[3]:
                    row.SetValue("NPV_BREAKS", 4)
                elif NPVValue > NPVPctList[3]:
                    row.SetValue("NPV_BREAKS", 5)
                cur.UpdateRow(row)
                row = cur.next()
            del row
            del cur

            # reclass NPV output (5 classes)
            NPVExpr = "1 1;2 2;3 3;4 4;5 5"
            gp.Reclassify_sa(outputNPV, "NPV_BREAKS", NPVExpr, outputNPV_rc, "NODATA")
            outputNPV_rc = AddField(outputNPV_rc, "VAL_RANGE", "TEXT", "75", "")
            cur = gp.UpdateCursor(outputNPV_rc)
            row = cur.Next()
            while row:
                NPVValue = row.GetValue("VALUE")
                if NPVValue == 1:
                    row.SetValue("VAL_RANGE", "1 - "+str(NPVPctList[0])+" thousands of US dollars (US$)")
                elif NPVValue == 2:
                    row.SetValue("VAL_RANGE", str(NPVPctList[0])+" - "+str(NPVPctList[1])+" thousands of US$")
                elif NPVValue == 3:
                    row.SetValue("VAL_RANGE", str(NPVPctList[1])+" - "+str(NPVPctList[2])+" thousands of US$")
                elif NPVValue == 4:
                    row.SetValue("VAL_RANGE", str(NPVPctList[2])+" - "+str(NPVPctList[3])+" thousands of US$")
                elif NPVValue == 5:
                    row.SetValue("VAL_RANGE", "Greater than "+str(NPVPctList[3])+" thousands of US$")
                cur.UpdateRow(row)
                row = cur.next()
            del row
            del cur
            
    except:
        raise Exception, msgValuation
    
    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(gp.workspace+"\\Output\\parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("WAVE ENERGY MODEL PARAMETERS\n")
    parafile.writelines("____________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

    # delete superfluous intermediate data
    del1 = [splineWP, splineCWE, neighbIntp, splineWPext, splineCWEext, splineWPint, splineCWEint, splineNPV, neighbIntp2, splineNPVext]  
    del2 = [LandPts_WGS84, GridPt_WGS84, WaveDataWCLyr, WaveDataECLyr, WaveDataEHLyr, WaveDataWHLyr, AOILyr, WaveData_clip, PointsAOICount]
    del3 = [WaveData_clipZ, DEM_prj]

    if AOI:
        deletelist = del1 + del2 + del3
    else:
        deletelist = del1 + del2
    for data in deletelist:
        if gp.exists(data):
            gp.delete_management(data)
    gp.AddMessage("")       
    del gp
except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())