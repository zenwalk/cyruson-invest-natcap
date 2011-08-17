# Marine InVEST: Grid the Seascape Tool (GS)
# Author: Gregg Verutes
# 05/09/11

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
msgCreateGrid = "\nError creating grid for seascape.  Make sure your AOI input is a polygon shapefile, projected in meters."
msgGISConversion = "\nError utilizing GIS tools to convert grid to feature class."

try:
    try:
        # get parameters
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
        gp.workspace = gp.GetParameterAsText(0)
        parameters.append("Workspace: "+ gp.workspace)
        AOI = gp.GetParameterAsText(1)
        parameters.append("Area of Interest (AOI): "+ AOI)
        cellsize = gp.GetParameterAsText(2)
        parameters.append("Analysis Cell Size (meters): "+ cellsize)
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

    const_rst = interws + "const_rst"
    consec_rst = interws + "consec_rst"
    constant_asc = interws + "constant_asc.txt"
    consec_asc = interws + "consec_asc.txt"
    GS_cellsize = outputws + "GS_"+cellsize+".shp"

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

    def grabProjection(data):
        dataDesc = gp.describe(data)
        sr = dataDesc.SpatialReference
        gp.OutputCoordinateSystem = sr
        strSR = str(gp.OutputCoordinateSystem)
        return strSR

    try:
        gp.AddMessage("\nCreating analysis grid...")
        if cellsize < 250:
            gp.AddError("Analysis cell size cannot be smaller than 250 meters.")
            raise Exception

        # various checks and preps
        gp.Extent = AOI
        checkProjections(AOI)
        projection = grabProjection(AOI)
        gp.OutputCoordinateSystem = projection
        AOI = AddField(AOI, "VALUE", "SHORT", "0", "0")
        gp.CalculateField_management(AOI, "VALUE", "1", "VB")
        gp.FeatureToRaster_conversion(AOI, "VALUE", const_rst, cellsize)
        gp.RasterToASCII_conversion(const_rst, constant_asc)

        # read in constant ascii and create consecutive ascii  
        inputfile = open(constant_asc, 'r')
        outputfile = open(consec_asc, 'w')
        GridID = []
        count = 1
        for line in inputfile.readlines():
            if count <= 6:
                outputfile.write(line)
                count = count + 1
            GridID.append([])
            for i in line.split():
                GridID[-1].append(i)

        # get number of cols and rows        
        NumCols = int(GridID[0][1])
        NumRows = int(GridID[1][1])

        GridID = GridID[6:]
        GridID = [[(int(y)) for y in x] for x in GridID]
        inputfile.close()

        # populate unique IDs in place of "1"
        counter = 1
        for i in range(0,NumRows):
            for j in range(0,NumCols):
                outputfile.write(str(counter) + " ")
                counter = counter + 1
            outputfile.write('\n')
        outputfile.close()
    except:
        gp.AddError(msgCreateGrid)
        raise Exception

    try:
        # run ArcGIS conversion tools
        gp.ASCIIToRaster_conversion(consec_asc, consec_rst, "INTEGER")
        gp.RasterToPolygon_conversion(consec_rst, GS_cellsize, "NO_SIMPLIFY")
        gp.DefineProjection_management(GS_cellsize, projection)
        GS_cellsize = AddField(GS_cellsize, "VALUE", "LONG", "0", "0")
        GS_cellsize = AddField(GS_cellsize, "CELL_SIZE", "LONG", "0", "0")
        gp.CalculateField_management(GS_cellsize, "VALUE", "[FID]+1", "VB")
        gp.CalculateField_management(GS_cellsize, "CELL_SIZE", cellsize, "VB")
        gp.DeleteField_management(GS_cellsize, "ID;GRIDCODE")
    except:
        gp.AddError(msgGISConversion)
        raise Exception


    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(outputws+"parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("GRID THE SEASCAPE PARAMETERS\n")
    parafile.writelines("____________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

    # delete superfluous intermediate data
    del1 = [const_rst, consec_rst, constant_asc, consec_asc, interws+"constant_asc.prj"]
    deletelist = del1
    for data in deletelist:
        if gp.exists(data):
            gp.delete_management(data)
    gp.AddMessage("")
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())