# Marine InVEST: Difference Tool
# Author: Gregg Verutes
# 08/25/11

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
msgEraseFunction = "\nError creating difference raster."

try:
    try:
        # get parameters
        currentRst = gp.GetParameterAsText(0)
        futureRst = gp.GetParameterAsText(1)
        outputRst = gp.GetParameterAsText(2)
    except:
        raise Exception, msgArguments + gp.GetMessages(2)

    gp.workspace = outputRst[:-(1+len(outputRst.split("\\")[-1]))]

    # intermediate variables
    cur_isNull = gp.workspace+"\\cur_isnull"
    fut_isNull = gp.workspace+"\\fut_isnull"
    cur_con = gp.workspace+"\\cur_con"
    fut_con = gp.workspace+"\\fut_con"

    gp.Extent = "MAXOF"

    try:
        # convert no data to "0" without losing floating point accuracy
        gp.IsNull_sa(currentRst, cur_isNull)
        gp.Con_sa(cur_isNull, "0", cur_con, currentRst, "VALUE = 1")
        gp.IsNull_sa(futureRst, fut_isNull)
        gp.Con_sa(fut_isNull, "0", fut_con, futureRst, "VALUE = 1")
        
        # subtract current from future
        gp.Minus_sa(fut_con, cur_con, outputRst)

        # delete intermediate data
        gp.delete_management(cur_isNull)
        gp.delete_management(fut_isNull)
        gp.delete_management(cur_con)
        gp.delete_management(fut_con)
    except:
        gp.AddError(msgEraseFunction)
        raise Exception

    gp.AddMessage("")
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())