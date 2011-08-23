# Marine InVEST: Difference Tool
# Author: Gregg Verutes
# 08/23/11

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
    gp.Extent = "MAXOF"

    try:
        # subtract current from future
        DiffExpr = futureRst+" - "+currentRst
        gp.SingleOutputMapAlgebra_sa(DiffExpr, outputRst)
    except:
        gp.AddError(msgEraseFunction)
        raise Exception

    gp.AddMessage("")
    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())