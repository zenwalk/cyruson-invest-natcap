# Marine InVEST: Visual Impact from Objects (Aesthetic Quality)
# Authors: Gregg Verutes, Mike Papenfus
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
msgVShed = "\nError conducting viewshed analysis."
msgIntersect = "\nError calculating overlap between viewshed output and visual polygons."
msgPopStats = "\nError conducting population statistics."
msgNumPyNo = "NumPy extension is required to run the Aesthetic Quality model.  Please consult the Marine InVEST FAQ for instructions on how to install."
msgSciPyNo = "SciPy extension is required to run the Aesthetic Quality model.  Please consult the Marine InVEST FAQ for instructions on how to install."

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
        NegPoints = gp.GetParameterAsText(3)
        parameters.append("Point Features Impacting Aesthetic Quality: "+ NegPoints)
        DEM = gp.GetParameterAsText(4)
        parameters.append("Digital Elevation Model (DEM): "+ DEM)
        RefractCoeff = float(gp.GetParameterAsText(5))
        parameters.append("Refractivity Coefficient: "+ str(RefractCoeff))
        globalPop = gp.GetParameterAsText(6)
        parameters.append("Population Raster: "+ globalPop)
        visualPolys = gp.GetParameterAsText(7)
        parameters.append("Polygon Features for Overlap Analysis: "+ visualPolys)

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
    vshed = interws + "vshed"
    vshed_rc1 = interws + "vshed_rc1"
    vshed_rc2 = interws + "vshed_rc2"
    zstatsPop = interws + "zstatsPop.dbf"
    vshed_2poly = interws + "vshed_2poly.shp"
    vp_prj_lyr = interws + "vp_prj.lyr"
    vshed_vp_intrsct = interws + "vshed_vp_intrsct.shp"
    vshed_vp_intrsct_lyr = interws + "vshed_vp_intrsct.lyr"
    vp_prj = outputws + "vp_overlap.shp"

    pop_prj = outputws + "pop_prj"
    vshed_qual = outputws + "vshed_qual"
    PopHTML = outputws + "PopulationStats_"+now.strftime("%Y-%m-%d-%H-%M")+".html"


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

    def ckDatum(thedata):
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
            gp.AddError(thedata+" is not a valid input.\nThe model requires a projected AOI with \"WGS84\" datum.\nSee InVEST FAQ document for how to reproject datasets.")
            raise Exception

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

    def getQuartiles(list):
        QrtList = []
        QrtList.append(stats.scoreatpercentile(list, 25))
        QrtList.append(stats.scoreatpercentile(list, 50))
        QrtList.append(stats.scoreatpercentile(list, 75))
        QrtList.append(stats.scoreatpercentile(list, 100))
        return QrtList

    def findNearest(array,value):
        idx=(np.abs(array-value)).argmin()
        return array[idx]


    ##############################################################
    ######## CHECK INPUTS, SET ENVIRONMENTS, & DATA PREP #########
    ##############################################################
    
    try:
        gp.AddMessage("\nChecking the inputs...")  
        # call various checking functions
        ckDatum(AOI)
        ckProjection(AOI)
        ckProjection(NegPoints)
        ckProjection(DEM)
        compareProjections(NegPoints, DEM)
        if visualPolys:
            ckProjection(visualPolys)  

        # set environments
        gp.Extent = AOI
        gp.Mask = AOI

        cellsizeDEM = checkCellSize(DEM)
        if cellsize:
            if int(cellsize) < int(cellsizeDEM):
                gp.AddError("The cell size input is too small.\nThe model requires the cell size to be equal to or larger than DEM input's cell size.")
                raise Exception

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


    ###########################################################
    ################### VIEWSHED ANALYSIS #####################
    ###########################################################
    
    try:
        gp.AddMessage("\nConducting the viewshed analysis...")
        gp.Viewshed_sa(DEM_vs, NegPoints, vshed, "1", "CURVED_EARTH", RefractCoeff)
        # get number of points in 'NegPoints'
        NegPointsCount = gp.GetCount_management(NegPoints)
        # reclassify 'vshed' output for overlap analysis and population stats
        if visualPolys:
            gp.Reclassify_sa(vshed, "VALUE", "0 NODATA;1 "+str(NegPointsCount)+" 1", vshed_rc1, "DATA")
        gp.Reclassify_sa(vshed, "VALUE", "0 0;1 "+str(NegPointsCount)+" 1", vshed_rc2, "DATA")
        
        # populate vshed values in list
        gp.BuildRasterAttributeTable_management(vshed, "Overwrite")
        cur = gp.UpdateCursor(vshed)
        row = cur.Next()
        vshedList = []      
        while row:
             cellValue = int(row.GetValue("VALUE"))
             cellCount = int(row.GetValue("COUNT"))
             if row.VALUE > 0:
                 for i in range(0,cellCount):
                     vshedList.append(cellValue)
             cur.UpdateRow(row)
             row = cur.next()
        del row, cur

        # create a list for breaks (25, 50, 75, 100 Percentiles)
        gp.AddMessage("...classifying the results in quartiles") 
        QuartList = getQuartiles(vshedList)
        QuartExpr = "0 0;1 "+str(int(QuartList[0]))+" 1;"\
                    +str(int(QuartList[0]))+" "+str(int(QuartList[1]))+" 2;"\
                    +str(int(QuartList[1]))+" "+str(int(QuartList[2]))+" 3;"\
                    +str(int(QuartList[2]))+" "+str(int(QuartList[3]))+" 4"
        gp.Reclassify_sa(vshed, "VALUE", QuartExpr, vshed_qual, "DATA")
        gp.BuildRasterAttributeTable_management(vshed_qual, "Overwrite")
        
        vshed_qual = AddField(vshed_qual, "VIS_QUAL", "TEXT", "50", "0")
        vshed_qual = AddField(vshed_qual, "VAL_BREAKS", "TEXT", "50", "0")

        # run through "vshed_qual" raster and add qual labels
        cur = gp.UpdateCursor(vshed_qual, "", "", "VALUE; VIS_QUAL; VAL_BREAKS")
        row = cur.Next()
        while row:
            QualValue = int(row.GetValue("VALUE"))
            if QualValue == 0:
                row.SetValue("VIS_QUAL", "No Visual Impact")
                row.SetValue("VAL_BREAKS", "0 Features Visible")   
            elif QualValue == 1:
                row.SetValue("VIS_QUAL", "Low Visual Impact")
                row.SetValue("VAL_BREAKS", "0 to "+str(int(QuartList[0]))+" Features Visible")
            elif QualValue == 2:
                row.SetValue("VIS_QUAL", "Low-Medium Visual Impact")
                row.SetValue("VAL_BREAKS", str(int(QuartList[0]))+" to "+str(int(QuartList[1]))+" Features Visible")
            elif QualValue == 3:
                row.SetValue("VIS_QUAL", "Medium-High Visual Impact")
                row.SetValue("VAL_BREAKS", str(int(QuartList[1]))+" to "+str(int(QuartList[2]))+" Features Visible")
            else:
                row.SetValue("VIS_QUAL", "High Visual Impact")
                row.SetValue("VAL_BREAKS", str(int(QuartList[2]))+" to "+str(int(QuartList[3]))+" Features Visible")
            cur.UpdateRow(row)
            row = cur.Next()
        del row, cur
    except:
        raise Exception, msgVShed


    ##############################################
    ############## INTERSECT ANALYSIS ############
    ##############################################
    
    try:
        if visualPolys:
            gp.AddMessage("\nCalculating overlap between viewshed output and visual polygons...\n")
            gp.RasterToPolygon_conversion(vshed_rc1, vshed_2poly, "NO_SIMPLIFY", "VALUE")
            gp.Select_analysis(visualPolys, vp_prj, "") # copy 'visualPolys' into outputws
            
            # intersect clipped visual polys with reclassed viewshed poly
            expr = str(vshed_2poly)+" 1; "+str(vp_prj)+" 2"
            gp.Intersect_analysis(expr, vshed_vp_intrsct, "ALL", "", "INPUT")

            # add three fields
            vshed_vp_intrsct = AddField(vshed_vp_intrsct, "AreaVS", "DOUBLE", "0", "0")
            vp_prj = AddField(vp_prj, "AreaVP", "DOUBLE", "0", "0")
            vp_prj = AddField(vp_prj, "AreaVShed", "SHORT", "0", "0")

            # calculate two fields
            gp.CalculateField_management(vshed_vp_intrsct, "AreaVS", "!shape.area@acres!", "PYTHON", "")
            gp.CalculateField_management(vp_prj, "AreaVP", "!shape.area@acres!", "PYTHON", "")

            # make feature layer for viewshed and visual poly fc intersection
            gp.MakeFeatureLayer_management(vshed_vp_intrsct, vshed_vp_intrsct_lyr, "", "", "")

            # run through visual poly layer one by one and calculate area overlap
            cur = gp.UpdateCursor(vp_prj, "", "", "AreaVP; AreaVShed; FID")
            row = cur.Next()

            while row:
                vpArea = float(row.GetValue("AreaVP"))
                FID = row.GetValue("FID")
                gp.MakeFeatureLayer_management(vp_prj, vp_prj_lyr, "\"FID\" = "+str(FID), "", "")
                selectVshed = gp.SelectLayerByLocation_management(vshed_vp_intrsct_lyr, "INTERSECT", vp_prj_lyr, "", "NEW_SELECTION")

                cur2 = gp.UpdateCursor(selectVshed, "", "", "AreaVS")
                row2 = cur2.Next()
                AreaVSSum = 0.0
                while row2:
                    AreaVSSum = AreaVSSum + float(row2.GetValue("AreaVS"))
                    cur2.UpdateRow(row2)
                    row2 = cur2.Next()
                del cur2    
                del row2
                 
                PctOverlap = ceil((AreaVSSum/vpArea)*100)
                if PctOverlap > 100:
                    PctOverlap = 100
                if PctOverlap == 0:
                    row.SetValue("AreaVShed", 0)
                else:
                    row.SetValue("AreaVShed", PctOverlap)
                cur.UpdateRow(row)
                row = cur.Next()
            del cur    
            del row
    except:
        raise Exception, msgIntersect


    ###########################################################
    ################# POPULATION SUMMARY ######################
    ###########################################################
    try:
        if globalPop:
            gp.AddMessage("\nCalculating population statistics...\n")
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
            del row, cur
            latAOI = midList[1]
            longAOI = midList[0]

            # based on centroid of AOI, calculate latitude for approximate projection cell size
            latList = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0]
            cellSizeList = [927.99, 924.39, 920.81, 917.26, 913.74, 902.85, 892.22, 881.83, 871.69, 853.83, 836.68, 820.21, 804.38, 778.99, 755.17, 732.76, 711.64, 679.16, 649.52, 622.35, 597.37, 557.69, 522.96, 492.30, 465.03, 416.93, 377.84, 345.46, 318.19, 256.15, 214.36, 184.29, 161.62]
            latList = np.array(latList)
            latList_near = findNearest(latList,latAOI)
            latList_index = np.where(latList==latList_near)

            # project and clip global population raster
            gp.ProjectRaster_management(globalPop, pop_prj, projection, "BILINEAR", str(cellSizeList[latList_index[0]]), "", "", "")

            gp.CellSize = "" # set environment's cell size to default
            
            gp.BuildRasterAttributeTable_management(vshed_rc2, "Overwrite")
            gp.ZonalStatisticsAsTable_sa(vshed_rc2, "VALUE", pop_prj, zstatsPop, "DATA") # get population statistics for vshed

            # search through zonal stats table for population within viewshed
            cur = gp.UpdateCursor(zstatsPop)
            row = cur.Next()  
            while row:
                if row.GetValue("VALUE") == 0:
                    PopZero = int(row.GetValue("SUM"))
                else:
                    PopSome = int(row.GetValue("SUM"))
                row = cur.next()
            del row, cur

            # create html file output
            htmlfile = open(PopHTML, "w")
            htmlfile.write("<html>\n")
            htmlfile.write("<title>Marine InVEST</title>")
            htmlfile.write("<center><H1>Aesthetic Quality Model</H1><H2>(Visual Impact from Objects)</H2></center><br>")
            htmlfile.write("<br><HR><br><H2>Population Statistics</H2><table border=\"1\", cellpadding=\"5\"><tr>")
            htmlfile.write("<td align=\"center\"><b>Number of Features Visible</b></td><td align=\"center\"><b>Population (estimate)</b></td></tr>")
            htmlfile.write("<tr><td align=\"center\">None visible<br> (unaffected)</td><td align=\"center\">"+str(PopZero)+"</td>")
            htmlfile.write("<tr><td align=\"center\">1 or more<br>visible</td><td align=\"center\">"+str(PopSome)+"</td>")
            htmlfile.write("</table>")
            htmlfile.write("</html>")
            htmlfile.close()
    except:
        raise Exception, msgPopStats
    
    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(gp.workspace+"\\Output\\parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("VISUAL IMPACT FROM OBJECTS (AESTHETIC QUALITY MODEL) PARAMETERS\n")
    parafile.writelines("_______________________________________________________________\n\n")
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

    # delete superfluous intermediate data
    gp.workspace = interws
    del1 = [AOI_lyr, DEM_2poly, DEM_2poly_lyr, DEM_1_rc, DEM_land, DEM_sea, DEM_sea_rc]
    del2 = [vshed_rc1, vshed_2poly, vp_prj_lyr, vshed_vp_intrsct, vshed_vp_intrsct_lyr]
    del3 = [AOI_geo, vshed_rc2, zstatsPop]
    deletelist = del1 + del2 + del3
    for data in deletelist:
        if gp.exists(data):
            gp.delete_management(data)
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())