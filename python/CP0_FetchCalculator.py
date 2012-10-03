# Marine InVEST: Fetch Calculator Tool
# Authors: Gregg Verutes, Greg Guannel, Jeremy Davies 
# Coded for ArcGIS 9.3 and 10
# 09/27/12

from CPf_FetchTools import fetchGeoprocessor
import sys, string, os, datetime, array, time, datetime
from math import *
import arcgisscripting
import fpformat, operator

# create the geoprocessor object
gp = arcgisscripting.create()

# set output handling
gp.OverwriteOutput = 1
# check out extensions
gp.CheckOutExtension("spatial")
gp.CheckOutExtension("management")
gp.CheckOutExtension("analysis")
gp.CheckOutExtension("conversion")

# error messages
msgArguments = "\nProblem with arguments."
msgCheckDatum = "\nError checking the datum of inputs."
msgDataPrep = "\nError preparing the data inputs."
msgFetchCalc = "\nError calculating fetch distances."
msgFetchThreshDist = "\nError storing fetch threshold distance in a text file."
msgNumPyNo = "NumPy extension is required to run the Fetch Calculator tool.  Please consult the Marine InVEST FAQ document for instructions on how to install."
msgWin32ComNo = "PythonWin extension is required to run the Fetch Calculator tool.  Please consult the Marine InVEST FAQ document for instructions on how to install."

# import modules
try:
    import numpy
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
        gp.CreateFolder_management(gp.workspace, "scratch")
        gp.scratchWorkspace = gp.workspace + os.sep + "scratch" 
        parameters.append("Scratch Workspace: "+ gp.scratchWorkspace)
        landPoly = gp.GetParameterAsText(1)
        parameters.append("Land Polygon: "+ landPoly)
        landLine = gp.GetParameterAsText(2)
        parameters.append("Land Polyline: "+ landLine)
        areaFilter = gp.GetParameterAsText(3)
        if areaFilter:
            areaFilter = int(gp.GetParameterAsText(3))
        parameters.append("Land Area Filter (kilometers squared): "+ str(areaFilter))
        AOI = gp.GetParameterAsText(4)
        parameters.append("Area of Interest (AOI): "+ AOI)
        cellsize = int(gp.GetParameterAsText(5))
        parameters.append("Cell Size (meters): "+ str(cellsize))
        fetchDist = int(gp.GetParameterAsText(6))
        parameters.append("Fetch Distance Threshold (meters): "+ str(fetchDist))
    except:
        raise Exception, msgArguments + gp.GetMessages(2)
        
    # check and create folders
    try:
        thefolders=["intermediate", "Output"]
        for folder in thefolders:
            if not gp.exists(gp.workspace + folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        raise Exception, "Error creating folders"

    # intermediate and output directories
    outputws = gp.workspace + os.sep + "Output" + os.sep
    interws = gp.workspace + os.sep + "intermediate" + os.sep

    # prep variables
    landLine_clip = interws + "landLine_clip.shp"
    coast_select = interws + "coast_select.shp"
    coast_select_prj = interws + "coast_select_prj.shp"
    coast_rst = interws + "coast_rst"
    landPoly_clip = interws + "landPoly_clip.shp"
    coastpoly_rst = interws + "coastpoly_rst"
    aoi_prj = interws + "aoi_prj.shp"
    aoi_rst = interws + "aoi_rst"
    aoi_lyr = interws + "aoi_lyr.lyr"
    landsea = interws + "landsea"
    landsea_rst = interws + "landsea_rst"
    costsurf = interws + "costsurf"
    costsurf_e = interws + "costsurf_e"

    # output variables
    coastPoly_prj = outputws + "coastPoly_prj.shp"
    fetch_cmb = outputws + "fetch_cmb"
    fetch_threshDist = outputws + "fetch_threshDist.txt"


    ##############################################
    ###### COMMON FUNCTION AND CHECK INPUTS ######
    ##############################################
    # check that scratch workspace exists
    if gp.ScratchWorkspace == "":
        gp.AddError("A scratch workspace must be defined in you Environment Settings for this script to function properly.")
        raise Exception

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
            gp.AddError(thedata+" is not a valid input.\nThe model requires data inputs and a projection with the \"WGS84\" datum.\nPlease review InVEST FAQ guide on how to transform a layer's datum.")
            raise Exception

    def ckProjection(data):
        dataDesc = gp.describe(data)
        spatreflc = dataDesc.SpatialReference
        if spatreflc.Type <> 'Projected':
            gp.AddError(data +" does not appear to be projected.  It is assumed to be in meters.")
            raise Exception
        if spatreflc.LinearUnitName <> 'Meter':
            gp.AddError("This model assumes that "+data+" is projected in meters for area calculations.  You may get erroneous results.")
            raise Exception

    def grabProjection(data):
        dataDesc = gp.describe(data)
        sr = dataDesc.SpatialReference
        gp.OutputCoordinateSystem = sr
        strSR = str(gp.OutputCoordinateSystem)
        return strSR

    def AddField(FileName, FieldName, Type, Precision, Scale):
        fields = gp.ListFields(FileName, FieldName)
        field_found = fields.Next()
        if field_found:
            gp.DeleteField_management(FileName, FieldName)
        gp.AddField_management(FileName, FieldName, Type, Precision, Scale, "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        return FileName

    def fetchLength(inFloat, outFloat):
        # calculates the fetch from top to bottom
        # read the header information
        finHeaderName = inFloat[:-3] + "hdr"
        fin = open(finHeaderName, 'r')
        header = fin.readlines()
        fin.close()
        # store header fields we need for later
        ncols = int(header[0].split()[1])
        nrows = int(header[1].split()[1])
        xllcorner = float(header[2].split()[1])
        yllcorner = float(header[3].split()[1])
        fgp.cellsize = float(header[4].split()[1])
        # apply shift since final output seems to be shifted a bit
        xllcorner = xllcorner
        yllcorner = yllcorner
        header[2] = "xllcorner     %s\n" % str(xllcorner)
        header[3] = "yllcorner     %s\n" % str(yllcorner)
        # copy header into new float header
        foutHeaderName = outFloat[:-3] + "hdr"
        fout = open(foutHeaderName, 'w')
        fout.writelines(header)
        fout.close()
        # open data files for processing
        fin = open(inFloat, 'rb')
        fout = open(outFloat, 'wb')
        fetch = array.array("f", [(-1 * fgp.cellsize) for i in range(ncols)])
        # walk down the DEM from the top row
        for row in xrange(nrows):
            # get the next (first) elevation row from the DEM 
            land = array.array("f")
            land.fromfile(fin, ncols)
            # calculate the fetch for this row
            for col in xrange(ncols):
                if land[col] > 0:
                    # found land, reset fetch length to zero
                    fetch[col] = 0
                else:
                    # found water, increment fetch length by cellsize
                    if fetch[col] >= 0:
                        # bounded fetch, positive fetch length
                        fetch[col] = fetch[col] + fgp.cellsize
                    else:
                        # unbounded fetch (so far), negative fetch length
                        fetch[col] = fetch[col] - fgp.cellsize
            # write the fetch for this row out to the new float binary
            fetch.write(fout)
        # close up the files        
        fin.close()
        fout.close()

    def calcSingle(inGrid, windDir, windDirInt, fetchpath):
        # calculate fetch along one wind direction
        windDir = (-1 * windDir) % 360.0
        input = inGrid
        output = fgp.tempGrid()
        fgp.deleteGrid(output)
        pivot = fgp.rotateGrid(input, output, windDir)
        # converting raster to float binary
        input = output
        output = fgp.tempFloat()
        fgp.deleteFloat(output)
        fgp.rasterToFloat(input, output)
        fgp.deleteGrid(input)
        # calculate the fetch length
        input = output
        output = fgp.tempFloat()
        fgp.deleteFloat(output)
        fetchLength(input, output)
        fgp.deleteFloat(input)
        # converting fetch binary to raster
        input = output
        output = fgp.tempGrid()
        fgp.deleteGrid(output)
        fgp.floatToRaster(input, output)
        fgp.deleteFloat(input)
        # rotating raster into north alignment (reverse the rotation performed above)
        input = output
        output = fgp.tempGrid()
        fgp.deleteGrid(output)
        angle = (-1 * windDir) % 360.0
        fgp.rotateGrid(input, output, angle, pivot)
        fgp.deleteGrid(input)
        # clipping final grid
        input = output
        output = fgp.tempGrid()
        clipGrid = inGrid
        fgp.deleteGrid(output)
        fgp.clipToGrid(input, output, clipGrid)
        fgp.deleteGrid(input)
        # setting fetch length = 0.0 to NODATA
        input = output
        output = fgp.tempGrid()
        fgp.deleteGrid(output)
        gp.SetNull_sa(input, input, output, "Value = 0")
        fgp.deleteGrid(input)
        # convert grid to integer for allocation process
        input = output
        output = fgp.tempGrid()
        fgp.deleteGrid(output)
        gp.Int_sa(input, output)
        fgp.deleteGrid(input)
        # allocate a few cells length buffer around results to make sure we intersect the shoreline
        input = output
        RasterName = "dir" + str(windDirInt).zfill(3)
        outRaster = os.path.join(fetchpath, RasterName)
        tempFinal = fgp.tempGrid()
        fgp.deleteGrid(tempFinal)
        distance = str(4.0 * fgp.cellsize)
        Output_distance_raster = ""
        Output_direction_raster = ""
        gp.EucAllocation_sa(input, tempFinal, distance, "", str(fgp.cellsize), "VALUE", Output_distance_raster, Output_direction_raster)
        fgp.deleteGrid(input)
        # copy the final temp file into the finished directory and rename
        tempMask = fgp.tempGrid()
        fgp.deleteGrid(tempMask)
        gp.Con_sa(inGrid, "1", tempMask, "", "Value = 0")
        tempClip = fgp.tempGrid()
        fgp.deleteGrid(tempClip)
        gp.ExtractByMask_sa(tempFinal, tempMask, tempClip)
        gp.SingleOutputMapAlgebra_sa("INT ( " + tempClip + " + 0.5)", outRaster, "")
        # clean up the remaining temp files
        fgp.deleteGrid(tempFinal)
        fgp.deleteGrid(tempMask)
        fgp.deleteGrid(tempClip)

    def expandExtent(inGrid):
        dsc = gp.Describe(inGrid)          
        extentList = dsc.Extent.split(" ")
        XMin1 = float(extentList[0])
        YMin1 = float(extentList[1])
        XMax1 = float(extentList[2])
        YMax1 = float(extentList[3])
        XDiff = XMax1 - XMin1
        YDiff = YMax1 - YMin1
        gridCenterX = (XMin1 + (XDiff / 2))
        gridCenterY = (YMin1 + (YDiff / 2))
        hypDiff = pow(((XDiff * XDiff) + (YDiff * YDiff)), 0.5)
        hypDiff = (hypDiff / 2)
        XMin2 = gridCenterX - hypDiff
        YMin2 = gridCenterY - hypDiff
        XMax2 = gridCenterX + hypDiff
        YMax2 = gridCenterY + hypDiff 
        extentString = str(XMin2) + " " + str(YMin2) + " " + str(XMax2) + " " + str(YMax2)
        gp.Extent = extentString

        
    #########################################################################################
    #########################################################################################

    try:
        gp.AddMessage("\nChecking inputs and preparing data...")
        # check the datum of certain inputs
        theinputs=[landPoly, landLine, AOI]
        for input in theinputs:
            if input:
                checkDatum(input)
    except:
        gp.AddError(msgCheckDatum)
        raise Exception

    try:
        # write fetch threshold distance to a text file
        txtfile = open(fetch_threshDist, "w")
        txtfile.write(str(fetchDist))
        txtfile.close()
    except:
        gp.AddError(msgFetchThreshDist)
        raise Exception

    try:
        # prepare the data
        ckProjection(AOI)
        projection = grabProjection(AOI)
        gp.Clip_analysis(landLine, AOI, landLine_clip, "")
        if areaFilter:
            gp.Select_analysis(landLine_clip, coast_select, "\"area\" >= "+str(areaFilter))
        else:
            coast_select = landLine_clip
        gp.Project_management(coast_select, coast_select_prj, projection)
        coast_select_prj = AddField(coast_select_prj, "ID", "SHORT", "0", "")
        gp.CalculateField_management(coast_select_prj, "ID", 1, "VB")
        gp.FeatureToRaster_conversion(coast_select_prj, "ID", coast_rst, str(cellsize))
        gp.Clip_analysis(landPoly, AOI, landPoly_clip, "")
        gp.Project_management(landPoly_clip, coastPoly_prj, projection)
        coastPoly_prj = AddField(coastPoly_prj, "land", "SHORT", "0", "")
        gp.CalculateField_management(coastPoly_prj, "land", 2, "VB")
        gp.Project_management(AOI, aoi_prj, projection)
        AOI = AddField(aoi_prj, "land", "SHORT", "0", "")
        gp.CalculateField_management(aoi_prj, "land", 1, "VB")
        gp.FeatureToRaster_conversion(aoi_prj, "land", aoi_rst, str(cellsize))
        gp.FeatureToRaster_conversion(coastPoly_prj, "land", coastpoly_rst, str(cellsize))
        # set extent and cellsize
        gp.Extent = aoi_rst
        gp.CellSize = cellsize
        MergeExpr = "Merge("+coastpoly_rst+","+aoi_rst+")"
        gp.SingleOutputMapAlgebra_sa(MergeExpr, landsea)
        gp.Reclassify_sa(landsea, "VALUE", "1 0;2 1", landsea_rst, "DATA")
        gp.Reclassify_sa(landsea_rst, "VALUE", "1 NODATA;0 1", costsurf, "DATA")
        gp.Expand_sa(costsurf, costsurf_e, "1", "1")
    except:
        gp.AddError(msgDataPrep)
        raise Exception

    # list for the 16 directions (used for fetch, wind and veg calculations)
    inp = [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 202.5, 225, 247.5, 270, 292.5, 315, 337.5]

    try:
        gp.AddMessage("\nCalculating fetch distances...")
        gp = Dispatch("esriGeoprocessing.GPDispatch.1")
        fgp = fetchGeoprocessor(gp)
        dirName = "fetch"
        fetchpath = os.path.join(interws, dirName)
        gp.CreateFolder_management(interws, dirName)
        expandExtent(landsea_rst)
        dirList = ["dir000", "dir022", "dir045", "dir067", "dir090", \
                   "dir112", "dir135", "dir157", "dir180", "dir202", \
                   "dir225", "dir247", "dir270", "dir292", "dir315", "dir337"]
        for IndWindDir in range(0,len(inp)):
            windDir = float(inp[IndWindDir])
            windDirInt = int(windDir)
            calcSingle(landsea_rst, windDir, windDirInt, fetchpath)

        for i in range(0,len(dirList)):
            gp.Reclassify_sa(fetchpath+"\\"+dirList[i], "Value", "-1000000 0 65000;0 250 0; 250 1000 625;1000 2000 1500;\
                                                                  2000 3000 2500;3000 4000 3500;4000 5000 4500;5000 6000 5500;\
                                                                  6000 7000 6500;7000 8000 7500;8000 9000 8500;9000 10000 9500;\
                                                                  10000 11000 10500;11000 12000 11500;12000 13000 12500;\
                                                                  13000 14000 13500;14000 15000 14500;15000 16000 15500;\
                                                                  16000 17000 16500;17000 18000 17500;18000 19000 18500;\
                                                                  19000 20000 19500;20000 10000000 50000", fetchpath+"\\"+dirList[i]+"_rc", "DATA")
            
            gp.Expand_sa(fetchpath+"\\"+dirList[i]+"_rc", fetchpath+"\\"+dirList[i]+"_e", "2", "0; \
                                                          625;1500;2500;3500;4500;5500;6500;7500;8500;9500;10500;\
                                                          11500;12500;13500;14500;15500;16500;17500;18500;19500;50000")
        del gp
        
    except:
        gp.AddError(msgFetchCalc)
        raise Exception

    try:
        # recreate geoprocessor
        gp = arcgisscripting.create()
        # combine various index rasters
        gp.AddMessage("\nCombining fetch results and determining coastal exposure...")    
        gp.Mask = coast_rst
        gp.snapRaster = coast_rst
        gp.workspace = interws
        
        CmbExpr = fetchpath+"\\"+dirList[0]+"_e"
        for i in range(1,len(dirList)):
            CmbExpr = CmbExpr+";"+fetchpath+"\\"+dirList[i]+"_e"

        gp.Combine_sa(CmbExpr, fetch_cmb)
        fetch_cmb = AddField(fetch_cmb, "FFILT", "SHORT", "", "")

        dirCapList = ["DIR000_E", "DIR022_E", "DIR045_E", "DIR067_E", "DIR090_E", \
                      "DIR112_E", "DIR135_E", "DIR157_E", "DIR180_E", "DIR202_E", \
                      "DIR225_E", "DIR247_E", "DIR270_E", "DIR292_E", "DIR315_E", "DIR337_E"]

        cur = gp.UpdateCursor(fetch_cmb)
        row = cur.Next()
        while row:
            CellFetchList = []
            for i in range(0,len(dirCapList)):
                CellFetchList.append(row.GetValue(dirCapList[i]))
                
            count = 0
            for i in range(0,len(dirCapList)):
                if int(row.GetValue(dirCapList[i])) >= fetchDist:
                    count = count + 1
              
            if count == 0 or count == 1:
                row.SetValue("FFILT", 0)
            else:
                row.SetValue("FFILT", 1)
                
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur

    except:
        gp.AddError(msgCombineIndex)
        raise Exception

    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(outputws+"\\parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("FETCH CALCULATOR PARAMETERS\n")
    parafile.writelines("___________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())