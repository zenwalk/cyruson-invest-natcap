# Marine InVEST: Difference Tool
# Author: Gregg Verutes
# 09/29/11

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
msgDiffFunction = "\nError creating difference raster."

try:
    try:
        # get parameters
        currentRst = gp.GetParameterAsText(0)
        futureRst = gp.GetParameterAsText(1)
        Field = gp.GetParameterAsText(2)
        outputRst = gp.GetParameterAsText(3)
        delIntermediate = gp.GetParameterAsText(4)
    except:
        raise Exception, msgArguments + gp.GetMessages(2)

    gp.workspace = outputRst[:-(1+len(outputRst.split("\\")[-1]))]

    # intermediate variables
    cur_field = gp.workspace+"\\cur_field"
    fut_field = gp.workspace+"\\fut_field"
    cur_isNull = gp.workspace+"\\cur_isnull"
    fut_isNull = gp.workspace+"\\fut_isnull"
    cur_con = gp.workspace+"\\cur_con"
    fut_con = gp.workspace+"\\fut_con"
    diff_rst = gp.workspace+"\\diff_rst"

    gp.Extent = "MAXOF"

    try:
        if Field:
            # check that attribute table and specified table exists in both inputs
            fields = gp.ListFields(currentRst, Field)
            field_found = fields.Next()
            if not field_found:
                gp.AddError("Specified field: "+str(Field)+" does not exist within input raster "+str(currentRst)+".\nPlease make sure to build a raster attribute table for raster inputs or leave the 'Attribute Field' input blank.")
                raise Exception

            fields = gp.ListFields(futureRst, Field)
            field_found = fields.Next()
            if not field_found:
                gp.AddError("Specified field: "+str(Field)+" does not exist within input raster "+str(futureRst)+".\nPlease make sure to build a raster attribute table for raster inputs or leave the 'Attribute Field' input blank.")
                raise Exception

            gp.Lookup_sa(currentRst, Field, cur_field)
            gp.Lookup_sa(futureRst, Field, fut_field)
            
            gp.IsNull_sa(cur_field, cur_isNull)
            gp.Con_sa(cur_isNull, "0", cur_con, cur_field, "VALUE = 1")
            gp.IsNull_sa(fut_field, fut_isNull)
            gp.Con_sa(fut_isNull, "0", fut_con, fut_field, "VALUE = 1")

        else: # just subtract VALUE fields, don't need to build attribute table
            # convert no data to "0" without losing floating point accuracy
            gp.IsNull_sa(currentRst, cur_isNull)
            gp.Con_sa(cur_isNull, "0", cur_con, currentRst, "VALUE = 1")
            gp.IsNull_sa(futureRst, fut_isNull)
            gp.Con_sa(fut_isNull, "0", fut_con, futureRst, "VALUE = 1")
            
        # subtract current from future
        gp.Minus_sa(fut_con, cur_con, diff_rst)
        # remove "0" values
        SetNullExp = "setnull("+diff_rst+" == 0, "+diff_rst+")"
        gp.SingleOutputMapAlgebra_sa(SetNullExp, outputRst)
        
        # delete intermediate data
        if delIntermediate == 'true':
            if Field:
                gp.delete_management(cur_field)
                gp.delete_management(fut_field)
            gp.delete_management(cur_isNull)
            gp.delete_management(fut_isNull)
            gp.delete_management(cur_con)
            gp.delete_management(fut_con)
            gp.delete_management(diff_rst)

    except:
        gp.AddError(msgDiffFunction)
        raise Exception

    gp.AddMessage("")
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())