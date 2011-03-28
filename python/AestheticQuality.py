# Marine InVEST: Aesthetic Views Model (Viewshed Analysis)
# Authors: Gregg Verutes, Mike Papenfus
# Coded for ArcGIS 9.3 and 10
# 02/16/11

# import modules
import sys, string, os, datetime
import arcgisscripting, numpy
from scipy import stats
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
        parameters.append("Cell Size (meters): "+ cellsize)
        NegPoints = gp.GetParameterAsText(3)
        parameters.append("Point Features Contributing to Negative Aesthetic Quality: "+ NegPoints)
        DEM = gp.GetParameterAsText(4)
        parameters.append("Digital Elevation Model (DEM): "+ DEM)
        RefractCoeff = gp.GetParameterAsText(5)
        parameters.append("Refractivity Coefficient: "+ RefractCoeff)
        ClassType = gp.GetParameterAsText(6)
        parameters.append("Viewshed Visual Quality Classification Type: "+ ClassType)
        globalPop = gp.GetParameterAsText(7)
        parameters.append("Global Population Raster: "+ globalPop)
        visualPolys = gp.GetParameterAsText(8)
        parameters.append("Polygon Features for Overlap Analysis: "+ visualPolys)
        projection = gp.GetParameterAsText(9)
        parameters.append("Projection for Overlap Analysis: "+ projection)  
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
    vshed = interws + "vshed_raw"
    vshed_rc = interws + "vshed_rc"
    vshed_2poly = interws + "vshed_2poly.shp"
    vp_prj_lyr = interws + "vp_prj.lyr"
    vshed_vp_intrsct = interws + "vshed_vp_intrsct.shp"
    vshed_vp_intrsct_lyr = interws + "vshed_vp_intrsct.lyr"
    zstatsPop = interws + "zstatsPop.dbf"

    vshed_qual = outputws + "vshed_qual"
    vp_prj = outputws + "vp_overlap.shp"
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

    def getNatBreaks( dataList, numClass ):
        dataList.sort()
        mat1 = []
        for i in range(0,len(dataList)+1):
            temp = []
            for j in range(0,numClass+1):
                temp.append(0)
            mat1.append(temp)
        mat2 = []
        for i in range(0,len(dataList)+1):
            temp = []
            for j in range(0,numClass+1):
                temp.append(0)
            mat2.append(temp)
        for i in range(1,numClass+1):
            mat1[1][i] = 1
            mat2[1][i] = 0
            for j in range(2,len(dataList)+1):
                mat2[j][i] = float(1e30000) # 'inf' for Python 2.6
        v = 0.0
        for l in range(2,len(dataList)+1):
            s1 = 0.0
            s2 = 0.0
            w = 0.0
            for m in range(1,l+1):
                i3 = l - m + 1
                val = float(dataList[i3-1])
                s2 += val * val
                s1 += val
                w += 1
                v = s2 - (s1 * s1) / w
                i4 = i3 - 1
                if i4 != 0:
                    for j in range(2,numClass+1):
                        if mat2[l][j] >= (v + mat2[i4][j - 1]):
                            mat1[l][j] = i3
                            mat2[l][j] = v + mat2[i4][j - 1]
            mat1[l][1] = 1
            mat2[l][1] = v
        k = len(dataList)
        kclass = []
        for i in range(0,numClass+1):
            kclass.append(0)
        kclass[numClass] = float(dataList[len(dataList) - 1])
        countNum = numClass
        while countNum >= 2:
            id = int((mat1[k][countNum]) - 2)
            kclass[countNum - 1] = dataList[id]
            k = int((mat1[k][countNum] - 1))
            countNum -= 1
        return kclass


    ##############################################################
    ######## CHECK INPUTS, SET ENVIRONMENTS, & DATA PREP #########
    ##############################################################
    try:
        gp.AddMessage("\nChecking the inputs...")  
        # call various checking functions
        checkGeometry(AOI, "Polygon", "Area of Interest (AOI)")
        checkGeometry(NegPoints, "Point", "Features contributing to negative aesthetic quality")
        ckProjection(NegPoints)
        ckProjection(DEM) 
        compareProjections(NegPoints, DEM)

        if projection:
            checkDatum(projection)

        cellsizeDEM = checkCellSize(DEM)
        if cellsize:
            if int(cellsize) < int(cellsizeDEM):
                gp.AddError("The cell size input is too small.\nThe model requires the cell size to be equal to or larger than DEM input's cell size.")
                raise Exception

        if visualPolys:
            checkGeometry(visualPolys, "Polygon", "Features for overlap analysis")
            checkDatum(visualPolys)
            ckProjection(visualPolys)
            if not gp.exists(projection):
                gp.AddError("A projection input is required when specifying polygon features for overlap analysis.")
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

        # set environments
        gp.Extent = AOI
        gp.Mask = AOI

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
        gp.Viewshed_sa(DEM_vs, NegPoints, vshed, "1", "CURVED_EARTH", RefractCoeff)
        gp.ZonalStatisticsAsTable_sa(vshed, "VALUE", globalPop, zstatsPop, "DATA")

        # populate vshed values in list
        cur = gp.UpdateCursor(vshed)
        row = cur.Next()
        vshedList = []      
        while row:
             cellValue = row.VALUE
             cellCount = row.COUNT
             if row.VALUE > 0:
                 for i in range(0,cellCount):
                     vshedList.append(cellValue)
             cur.UpdateRow(row)
             row = cur.next()
        del row
        del cur

        # get number of points in "NegPoints"
        NegPointsCount = gp.GetCount_management(NegPoints)
        gp.Reclassify_sa(vshed, "Value", "0 NODATA;1 "+str(NegPointsCount)+" 1", vshed_rc, "DATA")

        if ClassType == "Quartiles":
            # create a list for breaks (25, 50, 75, 100 Percentiles)
            gp.AddMessage("...classifying the results in quartiles") 
            QuartList = getQuartiles(vshedList)
            QuartExpr = "0 0;1 "+str(int(QuartList[0]))+" 1;"+str(int(QuartList[0]))+" "+str(int(QuartList[1]))+" 2;"+str(int(QuartList[1]))+" "+str(int(QuartList[2]))+" 3;"+str(int(QuartList[2]))+" "+str(int(QuartList[3]))+" 4"
            gp.Reclassify_sa(vshed, "Value", QuartExpr, vshed_qual, "DATA")
        if ClassType == "Natural Breaks":
            # create 4 Natural Breaks
            gp.AddMessage("...classifying the results using natural breaks") 
            NatBreak4 = getNatBreaks(vshedList,4)
            NBExpr = "0 0;1 "+str(int(NatBreak4[1]))+" 1;"+str(int(NatBreak4[1]))+" "+str(int(NatBreak4[2]))+" 2;"+str(int(NatBreak4[2]))+" "+str(int(NatBreak4[3]))+" 3;"+str(int(NatBreak4[3]))+" "+str(int(NatBreak4[4]))+" 4"
            gp.Reclassify_sa(vshed, "Value", NBExpr, vshed_qual, "DATA")
        
        gp.BuildRasterAttributeTable_management(vshed_qual, "Overwrite")
        vshed_qual = AddField(vshed_qual, "VIS_QUAL", "TEXT", "50", "0")
        vshed_qual = AddField(vshed_qual, "VAL_BREAKS", "TEXT", "50", "0")

        # run through "vshed_qual" raster and add qual labels
        cur = gp.UpdateCursor(vshed_qual, "", "", "VALUE; VIS_QUAL; VAL_BREAKS")
        row = cur.Next()
        while row:
            QualValue = int(row.GetValue("VALUE"))
            if QualValue == 0:
                row.SetValue("VIS_QUAL", "0 UNAFFECTED (No Visual Impact)")
                row.SetValue("VAL_BREAKS", "0 Sites Visible")   
            if QualValue == 1:
                row.SetValue("VIS_QUAL", "1 HIGH (Low Visual Impact)")
                if ClassType == "Quartiles":
                    row.SetValue("VAL_BREAKS", "0 < Sites Visible <= "+str(int(QuartList[0])))
                if ClassType == "Natural Breaks":
                    row.SetValue("VAL_BREAKS", "0 < Sites Visible <= "+str(int(NatBreak4[1])))
            if QualValue == 2:
                row.SetValue("VIS_QUAL", "2 MEDIUM (Moderate Visual Impact)")
                if ClassType == "Quartiles":
                    row.SetValue("VAL_BREAKS", str(int(QuartList[0]))+" < Sites Visible <= "+str(int(QuartList[1])))
                if ClassType == "Natural Breaks":
                    row.SetValue("VAL_BREAKS", str(int(NatBreak4[1]))+" < Sites Visible <= "+str(int(NatBreak4[2]))) 
            if QualValue == 3:
                row.SetValue("VIS_QUAL", "3 LOW (High Visual Impact)")
                if ClassType == "Quartiles":
                    row.SetValue("VAL_BREAKS", str(int(QuartList[1]))+" < Sites Visible <= "+str(int(QuartList[2])))
                if ClassType == "Natural Breaks":
                    row.SetValue("VAL_BREAKS", str(int(NatBreak4[2]))+" < Sites Visible <= "+str(int(NatBreak4[3]))) 
            if QualValue == 4:
                row.SetValue("VIS_QUAL", "4 VERY LOW/POOR (Very High Visual Impact)")
                if ClassType == "Quartiles":
                    row.SetValue("VAL_BREAKS", str(int(QuartList[2]))+" < Sites Visible <= "+str(int(QuartList[3])))
                if ClassType == "Natural Breaks":
                    row.SetValue("VAL_BREAKS", str(int(NatBreak4[3]))+" < Sites Visible <= "+str(int(NatBreak4[4])))  
            cur.UpdateRow(row)
            row = cur.Next()
        del cur    
        del row

        # search through zonal stats table for population within viewshed
        cur = gp.UpdateCursor(zstatsPop)
        row = cur.Next()
        NumSitesList = []
        PopList = []      
        while row:
             NumSitesList.append(int(row.VALUE))
             PopList.append(int(row.SUM))
             cur.UpdateRow(row)
             row = cur.next()
        del row
        del cur
        PopSum = numpy.sum(PopList[1::])

        # create html file output
        htmlfile = open(PopHTML, "w")
        htmlfile.write("<html>\n")
        htmlfile.write("<title>Marine InVEST</title>")
        htmlfile.write("<center><H1>Viewshed Model (Aesthetic Quality)</H1></center><br>")
        htmlfile.write("This page contains population results from running the Marine InVEST Viewshed model.")
        htmlfile.write("<br><HR><H2>Population Statistics</H2><table border=\"1\", cellpadding=\"5\"><tr>")
        htmlfile.write("<td align=\"center\"><b>Number of Visible Sites</b></td><td align=\"center\"><b>Population</b></td></tr>")
        htmlfile.write("<tr><td align=\"center\">"+str(NumSitesList[0])+"<br> (unaffected)</td><td align=\"center\">"+str(PopList[0])+"</td>")
        htmlfile.write("<tr><td align=\"center\">1 or more<br>sites visible</td><td align=\"center\">"+str(PopSum)+"</td>")
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
            gp.RasterToPolygon_conversion(vshed_rc, vshed_2poly, "NO_SIMPLIFY", "Value")
            gp.Project_management(visualPolys, vp_prj, projection)
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

    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(gp.workspace+"\\Output\\parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("VIEWSHED ANALYSIS MODEL PARAMETERS\n")
    parafile.writelines("__________________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()


    # delete superfluous intermediate data
    del1 = [AOI_lyr, DEM_2poly, DEM_2poly_lyr, DEM_1_rc, DEM_land, DEM_sea, DEM_sea_rc, vshed_rc, ]
    del2 = [vshed_2poly, vp_prj_lyr, vshed_vp_intrsct,vshed_vp_intrsct_lyr, zstatsPop]
    deletelist = del1 + del2
    for data in deletelist:
        if gp.exists(data):
            gp.delete_management(data)
    del gp
except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())