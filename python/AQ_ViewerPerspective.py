# Marine InVEST: Viewer Perspective (Aesthetic Quality)
# Author: Gregg Verutes
# Coded for ArcGIS 9.3 and 10
# 10/05/11

# import modules
import sys, string, os, datetime, shlex
import arcgisscripting
from math import *

# create the geoprocessor object
gp = arcgisscripting.create()
# set output handling
gp.OverwriteOutput = 1
# check out any necessary licenses
gp.CheckOutExtension("spatial")
gp.CheckOutExtension("management")
gp.CheckOutExtension("analysis")
gp.CheckOutExtension("conversion")

# error messages
msgArguments = "\nProblem with arguments."
msgDataPrep = "\nError in preparing data."
msgVShed = "\nError conducting viewshed and population analysis."
msgNumPyNo = "NumPy extension is required to run the Aesthetic Quality model.  Please consult the Marine InVEST FAQ for instructions on how to install."

# import modules
try:
    import numpy as np
except:
    gp.AddError(msgNumPyNo)
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
        AOI = gp.GetParameterAsText(1)
        parameters.append("Area of Interest (AOI): "+ AOI)
        cellsize = gp.GetParameterAsText(2)
        if cellsize:
            cellsize = int(gp.GetParameterAsText(2))
        parameters.append("Cell Size (meters): "+ str(cellsize))
        DEM = gp.GetParameterAsText(3)
        parameters.append("Digital Elevation Model (DEM): "+ DEM)
        ObserverPts = gp.GetParameterAsText(4)
        parameters.append("Viewer (Point) Locations: "+ ObserverPts)
        visualPolys = gp.GetParameterAsText(5)
        parameters.append("Areas (Polygon) Impacting Aesthetic Quality: "+ visualPolys)
        CurWgtField = gp.GetParameterAsText(6)
        parameters.append("Current Weight Field for Areas Impacting Visual Quality: "+ CurWgtField)
        FutWgtField = gp.GetParameterAsText(7)
        parameters.append("Future Weight Field for Areas Impacting Visual Quality: "+ FutWgtField)
        RefractCoeff = float(gp.GetParameterAsText(8))
        parameters.append("Refractivity Coefficient: "+ str(RefractCoeff))
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

    DEM_1_rc = interws + "DEM_1_rc"
    DEM_2poly = interws + "DEM_2poly.shp"
    AOI_lyr = interws + "AOI_buff_lyr.lyr"
    DEM_2poly_lyr = interws + "DEM_2poly_lyr.lyr"
    DEM_land = interws + "DEM_land"
    DEM_sea = interws + "DEM_sea"
    DEM_sea_rc = interws + "DEM_sea_rc"
    DEM_vs = interws + "DEM_vs"

    AOI_geo = interws + "AOI_geo.shp"
    vp_cur = interws + "vp_cur"
    vp_fut = interws + "vp_fut"
    
    vshedcmb_pts = outputws + "vshedcmb_pts"
    ObserverPts_geo = outputws + "ObserverPts_geo.shp"
    vshedStatsHTML = outputws + "viewshedStats_"+now.strftime("%Y-%m-%d-%H-%M")+".html"


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

    def getDatum(thedata):
        desc = gp.describe(thedata)
        SR = desc.SpatialReference
        if SR.Type == "Geographic":
            strDatum = SR.DatumName         
        else:
            gp.OutputCoordinateSystem = SR
            strSR = str(gp.OutputCoordinateSystem)
            gp.OutputCoordinateSystem = ""
            n1 = strSR.find("GEOGCS")
            n2 = strSR.find("PROJECTION")
            strDatum = strSR[n1:n2-1]
        return strDatum

    def ckProjection(data):
        dataDesc = gp.describe(data)
        spatreflc = dataDesc.SpatialReference
        if spatreflc.Type <> 'Projected':
            gp.AddError(data +" does not appear to be projected.  It is assumed to be in meters.")
            raise Exception
        if spatreflc.LinearUnitName <> 'Meter':
            gp.AddError("This model assumes that "+data+" is projected in meters for area calculations.  You may get erroneous results.")
            raise Exception

    def compareProjections(NegPoints, DEM):
        if gp.describe(NegPoints).SpatialReference.name <> gp.describe(DEM).SpatialReference.name:
            gp.AddError("Projection Error: "+NegPoints+" is in a different projection from the DEM input data.  The two inputs must be the same projection to conduct viewshed analysis.")
            raise Exception

    def grabProjection(data):
        dataDesc = gp.describe(data)
        sr = dataDesc.SpatialReference
        gp.OutputCoordinateSystem = sr
        strSR = str(gp.OutputCoordinateSystem)
        return strSR

    def checkCellSize(thedata):
         desc=gp.Describe(thedata)
         CellWidth = desc.MeanCellWidth
         CellHeight = desc.MeanCellHeight
         return int((CellHeight+CellWidth)/2.0)

    def findNearest(array,value):
        idx=(np.abs(array-value)).argmin()
        return array[idx]


    ##############################################################
    ######## CHECK INPUTS, SET ENVIRONMENTS, & DATA PREP #########
    ##############################################################
    try:
        gp.AddMessage("\nChecking the inputs...")  
        # call various checking functions
        if gp.GetCount_management(ObserverPts) > 10 and cellsize < 250:
            gp.AddMessage("\nMany observer points and/or low cell size may take a long time to complete...")  
        ckProjection(ObserverPts)
        ckProjection(DEM)
        compareProjections(ObserverPts, DEM)
        if visualPolys:
            ckProjection(visualPolys)
        cellsizeDEM = checkCellSize(DEM)
        if cellsize:
            if int(cellsize) < int(cellsizeDEM):
                gp.AddError("The cell size input is too small.\nThe model requires the cell size to be equal to or larger than DEM input's cell size.")
                raise Exception
        
        # set environments
        gp.Extent = AOI
        gp.Mask = AOI

        # ID and grab poly info to lists
        if visualPolys and CurWgtField <> '':
            visualPolys = AddField(visualPolys, "POLY_ID", "LONG", "0", "0")
            gp.CalculateField_management(visualPolys, "POLY_ID", "[FID]+1", "VB")
            PolyList = []
            CurWghtList = []
            FutWghtList = []  
            cur = gp.UpdateCursor(visualPolys)
            row = cur.Next()   
            while row:
                PolyList.append(row.GetValue("POLY_ID"))
                CurWghtList.append(row.GetValue(CurWgtField))
                if FutWgtField <> '':
                    FutWghtList.append(row.GetValue(FutWgtField))
                cur.UpdateRow(row)
                row = cur.next()
            del row
            del cur

            # normalize the values
            # convert lists to numpy array
            CurWghtArray = np.array(CurWghtList)
            FutWghtArray = np.array(FutWghtList)
            rescaleNeg = 1.0
            rescalePos = 1.0
            if FutWgtField <> '':
                if max(CurWghtList) >= max(FutWghtList):
                    maxWght = max(CurWghtList)
                else:
                    maxWght = max(FutWghtList)
                    
                if min(CurWghtList) <= min(FutWghtList):
                    minWght = min(CurWghtList)
                else:
                    minWght = min(FutWghtList)

                if minWght < 0 and maxWght > 0:
                    if np.abs(minWght) < np.abs(maxWght):
                        rescaleNeg = np.abs(minWght)/np.abs(maxWght)
                    else:
                        rescalePos = np.abs(maxWght)/np.abs(minWght)
                CurWghtArray = np.where(CurWghtArray > 0.0, ((CurWghtArray/maxWght)*rescalePos), CurWghtArray)
                CurWghtArray = np.where(CurWghtArray < 0.0, ((CurWghtArray/(minWght*-1.0))*rescaleNeg), CurWghtArray)
                FutWghtArray = np.where(FutWghtArray > 0.0, ((FutWghtArray/maxWght)*rescalePos), FutWghtArray)            
                FutWghtArray = np.where(FutWghtArray < 0.0, ((FutWghtArray/(minWght*-1.0))*rescaleNeg), FutWghtArray) 
        
            else:
                if min(CurWghtList) < 0 and max(CurWghtList) > 0:
                    if np.abs(min(CurWghtList)) < np.abs(max(CurWghtList)):
                        rescaleNeg = np.abs((min(CurWghtList))/np.abs(max(CurWghtList)))
                    else:
                        rescalePos = np.abs((max(CurWghtList))/np.abs(min(CurWghtList)))
                CurWghtArray = np.where(CurWghtArray > 0.0, ((CurWghtArray/max(CurWghtList))*rescalePos), CurWghtArray)
                CurWghtArray = np.where(CurWghtArray < 0.0, ((CurWghtArray/(min(CurWghtList)*-1.0))*rescaleNeg), CurWghtArray)  

            # add normalized value fields to 'visualPolys'
            visualPolys = AddField(visualPolys, "CUR_NORM", "FLOAT", "0", "0")
            if FutWgtField <> '':
                visualPolys = AddField(visualPolys, "FUT_NORM", "FLOAT", "0", "0")
            # return normalized values to 'visualPolys'                                    
            cur = gp.UpdateCursor(visualPolys)
            row = cur.Next()
            id = 0
            while row:
                row.SetValue("CUR_NORM", CurWghtArray[id])
                if FutWgtField <> '':
                    row.SetValue("FUT_NORM", FutWghtArray[id])
                cur.UpdateRow(row)
                row = cur.next()
                id = id + 1
            del row
            del cur
       
            # convert 'visualPolys' to current and future rasters with weights as values
            gp.FeatureToRaster_conversion(visualPolys, "CUR_NORM", vp_cur, "50")
            if FutWgtField <> '':
                gp.FeatureToRaster_conversion(visualPolys, "FUT_NORM", vp_fut, "50")


        # change resolution if specified by user
        if cellsize:
            gp.CellSize = int(cellsize)
        else:
            gp.CellSize = int(cellsizeDEM)

        # check that AOI intersects DEM extent
        gp.Reclassify_sa(DEM, "Value", "-1000000 100000 1", DEM_1_rc, "DATA")
        gp.RasterToPolygon_conversion(DEM_1_rc, DEM_2poly, "SIMPLIFY", "Value")
        gp.MakeFeatureLayer_management(AOI, AOI_lyr, "", "", "")
        gp.MakeFeatureLayer_management(DEM_2poly, DEM_2poly_lyr, "", "", "")
        SelectAOI = gp.SelectLayerByLocation(AOI_lyr, "INTERSECT", DEM_2poly_lyr, "", "NEW_SELECTION")
        if gp.GetCount(SelectAOI) > 0:
            pass
        else:
            gp.AddError("The extent of the input area of interest does not intersect the DEM input.\nResize the AOI to fit within the DEM's extent.")
            raise Exception 

        gp.AddMessage("\nPreparing the DEM...")     
        LandExpr = "setnull("+DEM+" < 0, "+DEM+")"
        gp.SingleOutputMapAlgebra_sa(LandExpr, DEM_land, "")
        SeaExpr = "setnull("+DEM+" >= 0, "+DEM+")"
        gp.SingleOutputMapAlgebra_sa(SeaExpr, DEM_sea, "")
        gp.Reclassify_sa(DEM, "Value", "-1000000 0 0", DEM_sea_rc, "DATA")

        # merge DEM sea with DEM land; DEM sea has been flattened to 0
        MergeExpr = "Merge("+DEM_land+","+DEM_sea_rc+")"
        gp.SingleOutputMapAlgebra_sa(MergeExpr, DEM_vs)
    except:
        raise Exception, msgDataPrep


    ##############################################
    ############## VIEWSHED ANALYSIS #############
    ##############################################
    try:
        gp.AddMessage("\nConducting the viewshed analysis...")
        
        # ID and grab point info to lists
        ObserverPts = AddField(ObserverPts, "PTS_ID", "SHORT", "0", "0")        
        gp.CalculateField_management(ObserverPts, "PTS_ID", "[FID]+1", "VB")
        PtsList = []
        cur = gp.UpdateCursor(ObserverPts)
        row = cur.Next()   
        while row:
            PtsList.append(row.GetValue("PTS_ID"))
            cur.UpdateRow(row)
            row = cur.next()
        del row
        del cur

        # snap outputs to the (resampled) DEM resolution
        gp.snapRaster = DEM_vs
        
        # create array for zonal stats results (PTID, AREA NEG, AREA POS, SUM, EXISTENCE, LAT, LONG)
        NumPoints = int(len(PtsList))

        if visualPolys and CurWgtField <> '':
            AQStatsCurList = np.zeros(NumPoints*7, dtype=np.float64)
            AQStatsCurArray = np.reshape(AQStatsCurList, (NumPoints,7))
            if FutWgtField <> '':
                AQStatsFutList = np.zeros(NumPoints*7, dtype=np.float64)
                AQStatsFutArray = np.reshape(AQStatsFutList, (NumPoints,7))


        # start expression for combining individual vshed raster outputs
        SumExpr = outputws+"vshed"+str(PtsList[0])+"_pt"

        for i in range(0,len(PtsList)):
            # run viewshed on all points, one at a time
            gp.MakeFeatureLayer_management(ObserverPts, interws+"Observer_pt"+str(PtsList[i])+".lyr", "\"PTS_ID\" = "+str(PtsList[i]), "", "")
            gp.Viewshed_sa(DEM_vs, interws+"Observer_pt"+str(PtsList[i])+".lyr", outputws+"vshed"+str(PtsList[i])+"_pt", "1", "CURVED_EARTH", str(RefractCoeff))   
            if i <> 0:
                SumExpr = SumExpr + " + " + outputws+"vshed"+str(PtsList[i])+"_pt" # continuing appending to expression

            if visualPolys and CurWgtField <> '':                 
                gp.Reclassify_sa(outputws+"vshed"+str(PtsList[i])+"_pt", "VALUE", "0 NODATA;1 1", interws+"vshd_rc"+str(PtsList[i]), "DATA")
                gp.SingleOutputMapAlgebra_sa("Float("+interws+"vshd_rc"+str(PtsList[i])+")", interws+"vshd_flt"+str(PtsList[i]))

                # create 'current' outputs; multiply by vp_cur
                CurMultExpr = interws+"vshd_flt"+str(PtsList[i])+" * "+vp_cur
                gp.SingleOutputMapAlgebra_sa(CurMultExpr, interws+"vshd_cur"+str(PtsList[i]))
                gp.Int_sa(interws+"vshd_cur"+str(PtsList[i]), interws+"vshd_int"+str(PtsList[i]))
                gp.BuildRasterAttributeTable_management(interws+"vshd_int"+str(PtsList[i]), "Overwrite")
                # if not empty, mark array and perform zonal stats   
                if gp.GetCount(interws+"vshd_int"+str(PtsList[i])) > 0:
                    AQStatsCurArray[i][4] = 1
                    gp.SetNull_sa(interws+"vshd_cur"+str(PtsList[i]), interws+"vshd_cur"+str(PtsList[i]), outputws+"vshed"+str(PtsList[i])+"_cur", '"VALUE" = 0')
                    gp.Reclassify_sa(outputws+"vshed"+str(PtsList[i])+"_cur", "VALUE", "0 2 1;-2 0 -1", interws+"vshd_rcc"+str(PtsList[i]), "DATA")
                    gp.ZonalStatisticsAsTable_sa(interws+"vshd_rcc"+str(PtsList[i]), "VALUE", outputws+"vshed"+str(PtsList[i])+"_cur", interws+"zstats_cur"+str(PtsList[i]), "DATA")
                    gp.Delete_management(interws+"vshd_rcc"+str(PtsList[i])) 
                else:
                    gp.Delete_management(interws+"vshd_cur"+str(PtsList[i])) # if empty raster, delete output

                if FutWgtField <> '':
                    # create 'future' outputs; multiply by vp_fut
                    FutMultExpr = interws+"vshd_flt"+str(PtsList[i])+" * "+vp_fut
                    gp.SingleOutputMapAlgebra_sa(FutMultExpr, interws+"vshd_fut"+str(PtsList[i]))
                    gp.Int_sa(interws+"vshd_fut"+str(PtsList[i]), interws+"vshd_int"+str(PtsList[i]))
                    gp.BuildRasterAttributeTable_management(interws+"vshd_int"+str(PtsList[i]), "Overwrite")
                    # if not empty, mark array and perform zonal stats                      
                    if gp.GetCount(interws+"vshd_int"+str(PtsList[i])) > 0:
                        AQStatsFutArray[i][4] = 1
                        gp.SetNull_sa(interws+"vshd_fut"+str(PtsList[i]), interws+"vshd_fut"+str(PtsList[i]), outputws+"vshed"+str(PtsList[i])+"_fut", '"VALUE" = 0')
                        gp.Reclassify_sa(outputws+"vshed"+str(PtsList[i])+"_fut", "VALUE", "0 2 1;-2 0 -1", interws+"vshd_rcf"+str(PtsList[i]), "DATA")
                        gp.ZonalStatisticsAsTable_sa(interws+"vshd_rcf"+str(PtsList[i]), "VALUE", outputws+"vshed"+str(PtsList[i])+"_fut", interws+"zstats_fut"+str(PtsList[i]), "DATA")
                        gp.Delete_management(interws+"vshd_rcf"+str(PtsList[i])) 
                    else:
                        gp.Delete_management(interws+"vshd_fut"+str(PtsList[i])) # if empty raster, delete output

        # create expression for combining individual vshed raster outputs
        gp.SingleOutputMapAlgebra_sa(SumExpr, vshedcmb_pts)


        if visualPolys and CurWgtField <> '': 
            for i in range(0,len(PtsList)):
                # search through zonal stats table for population within current viewshed
                if AQStatsCurArray[i][4] == 1:
                    cur = gp.UpdateCursor(interws+"zstats_cur"+str(PtsList[i]))
                    row = cur.Next()
                    while row:
                        AQStatsCurArray[i][0] = PtsList[i]
                        if row.GetValue("VALUE") < 0:
                            AQStatsCurArray[i][1] = (row.GetValue("AREA")/1000)
                        else:
                            AQStatsCurArray[i][2] = (row.GetValue("AREA")/1000)
                        AQStatsCurArray[i][3] = AQStatsCurArray[i][3] + row.GetValue("SUM")
                        row = cur.next()
                    del row, cur

                if FutWgtField <> '':
                    if AQStatsFutArray[i][4] == 1:
                        cur = gp.UpdateCursor(interws+"zstats_fut"+str(PtsList[i]))
                        row = cur.Next()
                        while row:
                            AQStatsFutArray[i][0] = PtsList[i]
                            if row.GetValue("VALUE") < 0:
                                AQStatsFutArray[i][1] = (row.GetValue("AREA")/1000)
                            else:
                                AQStatsFutArray[i][2] = (row.GetValue("AREA")/1000)
                            AQStatsFutArray[i][3] = AQStatsFutArray[i][3] + row.GetValue("SUM")
                            row = cur.next()
                        del row, cur

            # grab datum and reprojected 'ObserverPts' in geographic/unprojected
            geo_projection = getDatum(ObserverPts)
            gp.Project_management(ObserverPts, ObserverPts_geo, geo_projection)
            # grab latitude value of AOI polygons's centroid
            cur = gp.UpdateCursor(ObserverPts_geo)
            row = cur.Next()
            PtIndex = 0
            while row:
                feat = row.Shape
                midpoint = feat.Centroid
                midList = shlex.split(midpoint)
                AQStatsCurArray[PtIndex][5] = float(midList[1])
                AQStatsCurArray[PtIndex][6] = float(midList[0])
                if FutWgtField <> '':
                    AQStatsFutArray[PtIndex][5] = float(midList[1])
                    AQStatsFutArray[PtIndex][6] = float(midList[0])
                PtIndex = PtIndex + 1
                row = cur.Next()
            del cur, row

            # get projection from AOI
            projection = grabProjection(AOI)
            # return projected AOI to geographic (unprojected)
            geo_projection = getDatum(AOI)
            gp.Project_management(AOI, AOI_geo, geo_projection)
            # grab latitude value of AOI polygons's centroid
            cur = gp.UpdateCursor(AOI_geo)
            row = cur.Next()
            feat = row.Shape
            midpoint = feat.Centroid
            midList = shlex.split(midpoint)
            midList = [float(s) for s in midList]
            del cur, row
            latAOI = midList[1]
            longAOI = midList[0]

            # create html file output
            htmlfile = open(vshedStatsHTML, "w")
            htmlfile.write("<html>\n")
            htmlfile.write("<title>Marine InVEST</title>")
            htmlfile.write("<center><H1>Aesthetic Quality Model</H1><H2>(Viewer Perspective)</H2></center><br>")          
            htmlfile.write("<br><H2>Viewshed Statistics (by Observer)</H2><table border=\"1\", cellpadding=\"8\">\
                            <tr><th colspan=1></th><th colspan=6>EXTENT<br>Area of Polygons within Viewshed</th><th colspan=3>INTENSITY<br>Cumulative Rating <br>(All Polygons within Viewshed)</th></tr>\
                            <tr><th colspan=1></th><th colspan=3>Negative Views (sq. km)</th><th colspan=3>Positive Views (sq. km)</th></th><th colspan=3></th></tr>\
                            <td align=\"center\"><b>Point ID</b></td><td align=\"center\"><b>Current</b></td><td align=\"center\"><b>Future</b></td>\
                            <td align=\"center\"><b>Change</b></td><td align=\"center\"><b>Current</b></td>\
                            <td align=\"center\"><b>Future</b></td><td align=\"center\"><b>Change</b></td>\
                            <td align=\"center\"><b>Current</b></td><td align=\"center\"><b>Future</b></td><td align=\"center\"><b>Change</b></td></tr>")
            for i in range(0,len(PtsList)):
                htmlfile.write("<tr><td align=\"center\">"+str(PtsList[i])+"</td>")
                if FutWgtField <> '':
                    # negative area
                    if (AQStatsFutArray[i][1] - AQStatsCurArray[i][1]) > 0.0:
                        htmlfile.write("<td align=\"center\">"+str(AQStatsCurArray[i][1])+"</td><td align=\"center\">"+str(AQStatsFutArray[i][1])+"</td>\
                                        <td align=\"center\" bgcolor=\"#FF0000\">"+str(AQStatsFutArray[i][1]-AQStatsCurArray[i][1])+"</td>")
                    elif (AQStatsFutArray[i][1] - AQStatsCurArray[i][1]) == 0.0:
                        htmlfile.write("<td align=\"center\">"+str(AQStatsCurArray[i][1])+"</td><td align=\"center\">"+str(AQStatsFutArray[i][1])+"</td>\
                                        <td align=\"center\">"+str(AQStatsFutArray[i][1]-AQStatsCurArray[i][1])+"</td>")                     
                    else:
                        htmlfile.write("<td align=\"center\">"+str(AQStatsCurArray[i][1])+"</td><td align=\"center\">"+str(AQStatsFutArray[i][1])+"</td>\
                                        <td align=\"center\" bgcolor=\"#00FF00\">"+str(AQStatsFutArray[i][1]-AQStatsCurArray[i][1])+"</td>")

                    # positive area                    
                    if (AQStatsFutArray[i][2] - AQStatsCurArray[i][2]) > 0.0:
                        htmlfile.write("<td align=\"center\">"+str(AQStatsCurArray[i][2])+"</td><td align=\"center\">"+str(AQStatsFutArray[i][2])+"</td>\
                                        <td align=\"center\" bgcolor=\"#00FF00\">"+str(AQStatsFutArray[i][2]-AQStatsCurArray[i][2])+"</td>")
                    elif (AQStatsFutArray[i][2] - AQStatsCurArray[i][2]) == 0.0:
                        htmlfile.write("<td align=\"center\">"+str(AQStatsCurArray[i][2])+"</td><td align=\"center\">"+str(AQStatsFutArray[i][2])+"</td>\
                                        <td align=\"center\">"+str(AQStatsFutArray[i][2]-AQStatsCurArray[i][2])+"</td>")                     
                    else:
                        htmlfile.write("<td align=\"center\">"+str(AQStatsCurArray[i][2])+"</td><td align=\"center\">"+str(AQStatsFutArray[i][2])+"</td>\
                                        <td align=\"center\" bgcolor=\"#FF0000\">"+str(AQStatsFutArray[i][2]-AQStatsCurArray[i][2])+"</td>")

                else:
                    htmlfile.write("<td align=\"center\">"+str(AQStatsCurArray[i][1])+"</td><td align=\"center\">N/A</td><td align=\"center\">N/A</td>")
                    htmlfile.write("<td align=\"center\">"+str(AQStatsCurArray[i][2])+"</td><td align=\"center\">N/A</td><td align=\"center\">N/A</td>")
                    
                if FutWgtField <> '':
                    # intensity
                    if (AQStatsFutArray[i][3] - AQStatsCurArray[i][3]) > 0.0:
                        htmlfile.write("<td align=\"center\">"+str(round(AQStatsCurArray[i][3],3))+"</td><td align=\"center\">"+str(round(AQStatsFutArray[i][3],3))+"</td>\
                                        <td align=\"center\" bgcolor=\"#00FF00\">"+str(round(AQStatsFutArray[i][3]-AQStatsCurArray[i][3],3))+"</td>")

                    elif (AQStatsFutArray[i][3] - AQStatsCurArray[i][3]) == 0.0:
                        htmlfile.write("<td align=\"center\">"+str(round(AQStatsCurArray[i][3],3))+"</td><td align=\"center\">"+str(round(AQStatsFutArray[i][3],3))+"</td>\
                                        <td align=\"center\">"+str(round(AQStatsFutArray[i][3]-AQStatsCurArray[i][3],3))+"</td>")
                    else:
                        htmlfile.write("<td align=\"center\">"+str(round(AQStatsCurArray[i][3],3))+"</td><td align=\"center\">"+str(round(AQStatsFutArray[i][3],3))+"</td>\
                                        <td align=\"center\" bgcolor=\"#FF0000\">"+str(round(AQStatsFutArray[i][3]-AQStatsCurArray[i][3],3))+"</td>")
                else:
                    htmlfile.write("<td align=\"center\">"+str(round(AQStatsCurArray[i][3],3))+"</td><td align=\"center\">N/A</td><td align=\"center\">N/A</td>")
                htmlfile.write("</tr>")
            htmlfile.write("</tr></table></html>")
            htmlfile.close()
    except:
        raise Exception, msgVShed

    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(gp.workspace+"\\Output\\parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("VIEWER PERSPECTIVE (AESTHETIC QUALITY MODEL) PARAMETERS\n")
    parafile.writelines("_______________________________________________________\n\n")
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

    # delete superfluous intermediate data
    gp.workspace = interws
    del1 = [DEM_1_rc, DEM_2poly, AOI_lyr, DEM_2poly_lyr, DEM_land, DEM_sea, DEM_sea_rc, AOI_geo]
    del2 = []
    for i in range(0,len(PtsList)):
        del2.append("Observer_pt"+str(PtsList[i])+".lyr")
        del2.append("vshd_int"+str(PtsList[i]))
        del2.append("vshd_flt"+str(PtsList[i]))
        del2.append("vshd_rc"+str(PtsList[i]))
        del2.append("vshd_cur"+str(PtsList[i]))
        del2.append("vshd_fut"+str(PtsList[i]))
    deletelist = del1 + del2
    for data in deletelist:
        if gp.exists(data):
            gp.Delete_management(data)
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())