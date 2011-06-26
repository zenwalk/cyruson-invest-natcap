# Marine InVEST: Aesthetic Views Model (Viewshed Analysis)
# Authors: Gregg Verutes, Mike Papenfus
# Coded for ArcGIS 9.3 and 10
# 06/26/11

# import modules
import sys, string, os, datetime
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
msgIntersect = "\nError calculating overlap between viewshed output and visual polygons."
msgNumPyNo = "NumPy extension is required to run the Aesthetic Quality model.  Please consult the Marine InVEST FAQ for instructions on how to install."
msgSciPyNo = "SciPy extension is required to run the Aesthetic Quality model.  Please consult the Marine InVEST FAQ for instructions on how to install."

# import modules
try:
    import numpy
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
        gp.scratchWorkspace = gp.GetParameterAsText(0)
        parameters.append("Scratch Workspace: "+ gp.scratchWorkspace)
        AOI = gp.GetParameterAsText(1)
        parameters.append("Area of Interest (AOI): "+ AOI)
        cellsize = gp.GetParameterAsText(2)
        if cellsize:
            cellsize = int(gp.GetParameterAsText(2))
        parameters.append("Cell Size (meters): "+ str(cellsize))
        NegPointsCur = gp.GetParameterAsText(3)
        parameters.append("Current Point Features Contributing to Negative Aesthetic Quality: "+ NegPointsCur)
        NegPointsFut = gp.GetParameterAsText(4)
        parameters.append("Future Point Features Contributing to Negative Aesthetic Quality: "+ NegPointsFut)
        DEM = gp.GetParameterAsText(5)
        parameters.append("Digital Elevation Model (DEM): "+ DEM)
        RefractCoeff = float(gp.GetParameterAsText(6))
        parameters.append("Refractivity Coefficient: "+ str(RefractCoeff))
        globalPop = gp.GetParameterAsText(7)
        parameters.append("Population Raster: "+ globalPop)
        visualPolys = gp.GetParameterAsText(8)
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
    vshed_cur = interws + "vshed_cur"
    vshed_rcc = interws + "vshed_rcc"
    vshed_fut = interws + "vshed_fut"
    vshed_rcf = interws + "vshed_rcf"
    vshed_poly = interws + "vshed_poly.shp"    
    vshed_2poly = interws + "vshed_2poly.shp"
    vp_prj = interws + "vp_inter.shp"
    vp_prj_lyr = interws + "vp_prj_lyr.lyr"
    vp_prj_lyr2 = interws + "vp_prj_lyr2.lyr"
    vshed_vp_intrsct = interws + "vshed_vp_intrsct.shp"
    vshed_vp_intrsct_lyr = interws + "vshed_vp_intrsct.lyr"
    zstatsPop_cur = interws + "zstatsPop_cur.dbf"
    zstatsPop_fut = interws + "zstatsPop_fut.dbf"
    zstats_vp_cur = interws + "zstats_vp_cur.dbf"
    zstats_vp_fut = interws + "zstats_vp_fut.dbf"

    vp_overlap = outputws + "vp_overlap.shp"
    vshed_qualc = outputws + "vshed_qualc"
    vshed_qualf = outputws + "vshed_qualf"
    vshed_diff = outputws + "vshed_diff"

    PopHTML = outputws + "populationStats_"+now.strftime("%Y-%m-%d-%H-%M")+".html"


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

    def checkGeometry(thedata, Type, Message):
        if gp.Describe(thedata).ShapeType <> Type:
            raise Exception, "\nInvalid input: "+thedata+"\n"+Message+" must be of geometry type "+Type+"."

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


    ##############################################################
    ######## CHECK INPUTS, SET ENVIRONMENTS, & DATA PREP #########
    ##############################################################
    try:
        gp.AddMessage("\nChecking the inputs...")  
        # call various checking functions
        checkGeometry(AOI, "Polygon", "Area of Interest (AOI)")
        checkGeometry(NegPointsCur, "Point", "Features contributing to negative aesthetic quality")
        ckProjection(NegPointsCur)
        ckProjection(DEM) 
        compareProjections(NegPointsCur, DEM)
        if NegPointsFut:
            checkGeometry(NegPointsFut, "Point", "Features contributing to negative aesthetic quality")
            ckProjection(NegPointsFut) 
            compareProjections(NegPointsFut, DEM)
        
        cellsizeDEM = checkCellSize(DEM)
        if cellsize:
            if int(cellsize) < int(cellsizeDEM):
                gp.AddError("The cell size input is too small.\nThe model requires the cell size to be equal to or larger than DEM input's cell size.")
                raise Exception

        if visualPolys:
            checkGeometry(visualPolys, "Polygon", "Features for overlap analysis")
            ckProjection(visualPolys)

        if cellsize:
            gp.CellSize = int(cellsize)
        else:
            gp.CellSize = int(cellsizeDEM)

        # set environments
        gp.Extent = AOI
        gp.Mask = AOI

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


    ##########################################################
    ############## VIEWSHED & POPULATION ANALYSIS#############
    ##########################################################
    try:
        gp.AddMessage("\nConducting the viewshed analysis...")  
        gp.Viewshed_sa(DEM_vs, NegPointsCur, vshed_cur, "1", "CURVED_EARTH", str(RefractCoeff))
        if NegPointsFut:
            gp.Viewshed_sa(DEM_vs, NegPointsFut, vshed_fut, "1", "CURVED_EARTH", str(RefractCoeff))
            # subtract future from current
            gp.SingleOutputMapAlgebra_sa("( "+vshed_fut+" - "+vshed_cur+" )", vshed_diff, "")

        # get number of points in "NegPoints"
        # conduct zonal stats for population summary
        NegPointsCountCur = gp.GetCount_management(NegPointsCur)
        gp.Reclassify_sa(vshed_cur, "VALUE", "0 1;1 "+str(NegPointsCountCur)+" 2", vshed_rcc, "DATA")
        gp.BuildRasterAttributeTable_management(vshed_rcc, "Overwrite")
        gp.ZonalStatisticsAsTable_sa(vshed_rcc, "VALUE", globalPop, zstatsPop_cur, "DATA")
        if NegPointsFut:
            NegPointsCountFut = gp.GetCount_management(NegPointsFut)
            gp.Reclassify_sa(vshed_fut, "VALUE", "0 1;1 "+str(NegPointsCountFut)+" 2", vshed_rcf, "DATA")
            gp.BuildRasterAttributeTable_management(vshed_rcf, "Overwrite")
            gp.ZonalStatisticsAsTable_sa(vshed_rcf, "VALUE", globalPop, zstatsPop_fut, "DATA")

        # populate vshed_cur values in list
        gp.BuildRasterAttributeTable_management(vshed_cur, "Overwrite")
        cur = gp.UpdateCursor(vshed_cur)
        row = cur.Next()
        vshedList_cur = []      
        while row:
             cellValue = row.VALUE
             cellCount = row.COUNT
             if row.VALUE > 0:
                 for i in range(0,cellCount):
                     vshedList_cur.append(cellValue)
             cur.UpdateRow(row)
             row = cur.next()
        del row
        del cur

        # populate vshed_fut values in list
        if NegPointsFut:
            gp.BuildRasterAttributeTable_management(vshed_fut, "Overwrite")
            cur = gp.UpdateCursor(vshed_fut)
            row = cur.Next()
            vshedList_fut = []      
            while row:
                 cellValue = row.VALUE
                 cellCount = row.COUNT
                 if row.VALUE > 0:
                     for i in range(0,cellCount):
                         vshedList_fut.append(cellValue)
                 cur.UpdateRow(row)
                 row = cur.next()
            del row
            del cur        

        # if future points provided, determine which max value is higher
        maxDist = 1
        if NegPointsFut:
            if max(vshedList_fut) > max(vshedList_cur):
                maxDist = 2

        # create a list for breaks (25, 50, 75, 100 Percentiles)
        gp.AddMessage("...classifying the results in quartiles")
        if maxDist == 1:
            QuartList = getQuartiles(vshedList_cur)
        else:
            QuartList = getQuartiles(vshedList_fut)
        QuartExpr = "0 0;1 "+str(int(QuartList[0]))+" 1;"+str(int(QuartList[0]))+" "+str(int(QuartList[1]))+" 2;"+str(int(QuartList[1]))+" "+str(int(QuartList[2]))+" 3;"+str(int(QuartList[2]))+" "+str(int(QuartList[3]))+" 4"
        gp.Reclassify_sa(vshed_cur, "VALUE", QuartExpr, vshed_qualc, "DATA")
        gp.BuildRasterAttributeTable_management(vshed_qualc, "Overwrite")
        vshed_qualc = AddField(vshed_qualc, "VIS_QUAL", "TEXT", "50", "0")
        vshed_qualc = AddField(vshed_qualc, "VAL_BREAKS", "TEXT", "50", "0")
        # run through "vshed_qualc" raster and add qual labels
        cur = gp.UpdateCursor(vshed_qualc, "", "", "VALUE; VIS_QUAL; VAL_BREAKS")
        row = cur.Next()
        while row:
            QualValue = int(row.GetValue("VALUE"))
            if QualValue == 0:
                row.SetValue("VIS_QUAL", "0 UNAFFECTED (No Visual Impact)")
                row.SetValue("VAL_BREAKS", "0 Sites Visible")   
            if QualValue == 1:
                row.SetValue("VIS_QUAL", "1 HIGH (Low Visual Impact)")
                row.SetValue("VAL_BREAKS", "0 < Sites Visible <= "+str(int(QuartList[0])))
            if QualValue == 2:
                row.SetValue("VIS_QUAL", "2 MEDIUM (Moderate Visual Impact)")
                row.SetValue("VAL_BREAKS", str(int(QuartList[0]))+" < Sites Visible <= "+str(int(QuartList[1])))
            if QualValue == 3:
                row.SetValue("VIS_QUAL", "3 LOW (High Visual Impact)")
                row.SetValue("VAL_BREAKS", str(int(QuartList[1]))+" < Sites Visible <= "+str(int(QuartList[2])))
            if QualValue == 4:
                row.SetValue("VIS_QUAL", "4 VERY LOW/POOR (Very High Visual Impact)")
                row.SetValue("VAL_BREAKS", str(int(QuartList[2]))+" < Sites Visible <= "+str(int(QuartList[3])))
            cur.UpdateRow(row)
            row = cur.Next()
        del cur    
        del row      
        
        if NegPointsFut:
            gp.Reclassify_sa(vshed_fut, "VALUE", QuartExpr, vshed_qualf, "DATA")
            gp.BuildRasterAttributeTable_management(vshed_qualf, "Overwrite")
            vshed_qualf = AddField(vshed_qualf, "VIS_QUAL", "TEXT", "50", "0")
            vshed_qualf = AddField(vshed_qualf, "VAL_BREAKS", "TEXT", "50", "0")
            # run through "vshed_qualc" raster and add qual labels
            cur = gp.UpdateCursor(vshed_qualf, "", "", "VALUE; VIS_QUAL; VAL_BREAKS")
            row = cur.Next()
            while row:
                QualValue = int(row.GetValue("VALUE"))
                if QualValue == 0:
                    row.SetValue("VIS_QUAL", "0 UNAFFECTED (No Visual Impact)")
                    row.SetValue("VAL_BREAKS", "0 Sites Visible")   
                if QualValue == 1:
                    row.SetValue("VIS_QUAL", "1 HIGH (Low Visual Impact)")
                    row.SetValue("VAL_BREAKS", "0 < Sites Visible <= "+str(int(QuartList[0])))
                if QualValue == 2:
                    row.SetValue("VIS_QUAL", "2 MEDIUM (Moderate Visual Impact)")
                    row.SetValue("VAL_BREAKS", str(int(QuartList[0]))+" < Sites Visible <= "+str(int(QuartList[1])))
                if QualValue == 3:
                    row.SetValue("VIS_QUAL", "3 LOW (High Visual Impact)")
                    row.SetValue("VAL_BREAKS", str(int(QuartList[1]))+" < Sites Visible <= "+str(int(QuartList[2])))
                if QualValue == 4:
                    row.SetValue("VIS_QUAL", "4 VERY LOW/POOR (Very High Visual Impact)")
                    row.SetValue("VAL_BREAKS", str(int(QuartList[2]))+" < Sites Visible <= "+str(int(QuartList[3])))
                cur.UpdateRow(row)
                row = cur.Next()
            del cur    
            del row

        # search through zonal stats table for population within current viewshed
        cur = gp.UpdateCursor(zstatsPop_cur)
        row = cur.Next()
        PopListCur = []
        while row:
             PopListCur.append(int(row.SUM))
             cur.UpdateRow(row)
             row = cur.next()
        del row
        del cur

        if NegPointsFut:
            # search through zonal stats table for population within future viewshed
            cur = gp.UpdateCursor(zstatsPop_fut)
            row = cur.Next()
            PopListFut = []
            while row:
                 PopListFut.append(int(row.SUM))
                 cur.UpdateRow(row)
                 row = cur.next()
            del row
            del cur      

        # create html file output
        htmlfile = open(PopHTML, "w")
        htmlfile.write("<html>\n")
        htmlfile.write("<title>Marine InVEST</title>")
        htmlfile.write("<center><H1>Aesthetic Quality Model</H1></center><br>")
        htmlfile.write("This page contains population results from running the Marine InVEST Aesthetic Quality Model.")
        htmlfile.write("<br><HR><H2>Population Statistics</H2><table border=\"1\", cellpadding=\"5\"><tr>")
        if NegPointsFut:
            htmlfile.write("<td align=\"center\"><b>Number of Visible Sites</b></td><td align=\"center\"><b>Current Population</b></td><td align=\"center\"><b>Future Population</b></td></tr>")
            htmlfile.write("<tr><td align=\"center\">0<br> (unaffected)</td><td align=\"center\">"+str(PopListCur[0])+"</td><td align=\"center\">"+str(PopListFut[0])+"</td>")
            htmlfile.write("<tr><td align=\"center\">1 or more<br>sites visible</td><td align=\"center\">"+str(PopListCur[1])+"</td><td align=\"center\">"+str(PopListFut[1])+"</td>")
        else:
            htmlfile.write("<td align=\"center\"><b>Number of Visible Sites</b></td><td align=\"center\"><b>Current Population</b></td><td align=\"center\"><b>Future Population</b></td></tr>")
            htmlfile.write("<tr><td align=\"center\">0<br> (unaffected)</td><td align=\"center\">"+str(PopListCur[0])+"</td>")
            htmlfile.write("<tr><td align=\"center\">1 or more<br>sites visible</td><td align=\"center\">"+str(PopListCur[1])+"</td>")
        htmlfile.write("</table>")
        htmlfile.write("</html>")
        htmlfile.close()
    except:
        raise Exception, msgVShed


    ##############################################
    ############## INTERSECT ANALYSIS ############
    ##############################################
    try:
        if visualPolys:
            gp.AddMessage("\nCalculating overlap between viewshed output and visual polygons...\n")
            gp.RasterToPolygon_conversion(vshed_rcc, vshed_poly, "NO_SIMPLIFY", "Value")
            gp.Select_analysis(vshed_poly, vshed_2poly, "\"GRIDCODE\" = 2")
            gp.Select_analysis(visualPolys, vp_prj, "")
            HabVariable = AddField(vp_prj, "VALUE", "SHORT", "0", "0")  
            gp.CalculateField_management(vp_prj, "VALUE", "[FID]+1", "VB")            
    
            # intersect clipped visual polys with reclassed viewshed poly
            expr = str(vshed_2poly)+" 1; "+str(vp_prj)+" 2"
            gp.Intersect_analysis(expr, vshed_vp_intrsct, "ALL", "", "INPUT")

            # add three fields
            vshed_vp_intrsct = AddField(vshed_vp_intrsct, "AreaVS", "DOUBLE", "0", "0")
            vp_prj = AddField(vp_prj, "AreaVP", "DOUBLE", "0", "0")
            vp_prj = AddField(vp_prj, "AreaVShed", "SHORT", "0", "0")
            vp_prj = AddField(vp_prj, "RateCur", "DOUBLE", "0", "0")
            if NegPointsFut:
                vp_prj = AddField(vp_prj, "RateFut", "DOUBLE", "0", "0")
                vp_prj = AddField(vp_prj, "RateDiff", "DOUBLE", "0", "0")

            # calculate two fields
            gp.CalculateField_management(vshed_vp_intrsct, "AreaVS", "!shape.area@acres!", "PYTHON", "")
            gp.CalculateField_management(vp_prj, "AreaVP", "!shape.area@acres!", "PYTHON", "")

            # make feature layer for viewshed and visual poly fc intersection
            gp.MakeFeatureLayer_management(vshed_vp_intrsct, vshed_vp_intrsct_lyr, "", "", "")

            # run through visual poly layer one by one and calculate area overlap
            cur = gp.UpdateCursor(vp_prj, "", "", "AreaVP; AreaVShed; VALUE")
            row = cur.Next()

            while row:
                vpArea = float(row.GetValue("AreaVP"))
                Value = row.GetValue("VALUE")
                gp.MakeFeatureLayer_management(vp_prj, vp_prj_lyr, "\"VALUE\" = "+str(Value), "", "")
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

            gp.MakeFeatureLayer_management(vp_prj, vp_prj_lyr2, "", "", "")
            gp.ZonalStatisticsAsTable_sa(vp_prj, "VALUE", vshed_qualc, zstats_vp_cur, "DATA")
            gp.AddJoin_management(vp_prj_lyr2, "VALUE", zstats_vp_cur, "VALUE", "KEEP_COMMON")
            gp.CalculateField_management(vp_prj_lyr2, "vp_inter.RateCur", "[zstats_vp_cur.MEAN]", "VB", "")
            gp.RemoveJoin_management(vp_prj_lyr2, "zstats_vp_cur")
            
            if NegPointsFut:
                gp.ZonalStatisticsAsTable_sa(vp_prj, "VALUE", vshed_qualf, zstats_vp_fut, "DATA")
                gp.AddJoin_management(vp_prj_lyr2, "VALUE", zstats_vp_fut, "VALUE", "KEEP_COMMON")
                gp.CalculateField_management(vp_prj_lyr2, "vp_inter.RateFut", "[zstats_vp_fut.MEAN]", "VB", "")
                gp.RemoveJoin_management(vp_prj_lyr2, "zstats_vp_fut")

            gp.FeatureClassToFeatureClass_conversion(vp_prj_lyr2, outputws, "vp_overlap.shp", "")
            if NegPointsFut:
                gp.CalculateField_management(vp_overlap, "RateDiff", "!RateCur! - !RateFut!", "PYTHON") 
    except:
        raise Exception, msgIntersect

    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(gp.workspace+"\\Output\\parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("AESTHETIC QUALITY MODEL PARAMETERS\n")
    parafile.writelines("__________________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

##    # delete superfluous intermediate data
##    del1 = [AOI_lyr, DEM_2poly, DEM_2poly_lyr, DEM_1_rc, DEM_land, DEM_sea, DEM_sea_rc, vshed_rcc, vshed_rcf]
##    del2 = [vshed_poly, vshed_2poly, vp_prj_lyr, vp_prj_lyr2, vp_prj, vshed_vp_intrsct,vshed_vp_intrsct_lyr, zstatsPop_cur, zstatsPop_fut]
##    deletelist = del1 + del2
##    for data in deletelist:
##        if gp.exists(data):
##            gp.delete_management(data)
    del gp
    
except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())