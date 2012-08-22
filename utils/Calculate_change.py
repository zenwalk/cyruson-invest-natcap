# ---------------------------------------------------------------------------
# Calculate_change.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
#
# Calculates absolute and/or percent change between scenario outputs
# Optionally splits the results into two rasters, one with positive
#   values and the other with negative, to make symbolizing easier
# ---------------------------------------------------------------------------

# Import system modules
import sys, string, os, arcgisscripting, time, datetime

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out any necessary licenses
gp.CheckOutExtension("spatial")

# Allow overwriting of output files
gp.overwriteoutput = 1

verbose = True
cleanup = True

try:
    
    try:

        # Make parameters array, and later write input parameter values to an output file
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))

        # Script arguments
        gp.AddMessage ("\nValidating arguments..." )
        
        # Directory where output folder and files will be created 
        gp.workspace = gp.GetParameterAsText(0)
        parameters.append("Workspace: " + gp.workspace)

        # Scenario 1
        scenario1 = gp.GetParameterAsText(1)
        parameters.append("Scenario 1: " + scenario1)

        # Scenario 2
        scenario2 = gp.GetParameterAsText(2)
        parameters.append("Scenario 2: " + scenario2)

        # Calculate change?
        do_change = gp.GetParameterAsText(3)
        if do_change =='true':
            do_change = True
            parameters.append("Calculate change: Yes")
        else:
            do_change = False
            parameters.append("Calculate change: No")

        # Calculate percent change?
        do_percent_change = gp.GetParameterAsText(4)
        if do_percent_change =='true':
            do_percent_change = True
            parameters.append("Calculate percent change: Yes")
        else:
            do_percent_change = False
            parameters.append("Calculate percent change: No")

        # Split the results into positive and negative rasters?
        do_split = gp.GetParameterAsText(5)
        if do_split =='true':
            if (do_change or do_percent_change):
                do_split = True
                parameters.append("Split results: Yes")
            else:
                gp.AddError("\nError: If Split is selected, Change and/or Percent Change must also be selected.\n")
                raise Exception
        else:
            do_split = False
            parameters.append("Split results: No")

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(6)
        parameters.append("Suffix: " + Suffix)

        if (Suffix == "") or (Suffix == string.whitespace) or (Suffix == "#"):
            Suffix = ""
        else:
            Suffix = "_" + Suffix
            
    except:
        gp.AddError("Error in input arguments: " + gp.GetMessages(2))
        raise Exception


    # Check and create output folders
    try:
        thefolders=["Output", "Intermediate"]
        for folder in thefolders:
            if not gp.Exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        gp.AddError("Error creating folders: " + gp.GetMessages())
        raise Exception


    # Output files
    try:
        # Base output/working directories
        postprocws = gp.workspace + os.sep + "Output" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        
        # Intermediate variables
        x100 = "100"
        
        # Output layers
        change = postprocws + "change_" + Suffix + ".tif"
        percent_change =  postprocws + "percent_change_" + Suffix + ".tif"
        pchange_lt_zero = postprocws + "percent_change_lt0_" + Suffix + ".tif"
        pchange_gte_zero = postprocws + "percent_change_gte0_" + Suffix + ".tif"
        change_lt_zero = postprocws + "change_lt0_" + Suffix + ".tif"
        change_gte_zero = postprocws + "change_gte0_" + Suffix + ".tif"
        
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages(2))
        raise Exception


##    # Add suffix to end of output filenames
##    # Constrain length of output raster filenames to 13 characters
##    try:
##        Outputnames = [change, percent_change, change_lt_zero, change_gte_zero, pchange_lt_zero, pchange_gte_zero]           
##        num = 0
##        for x in Outputnames:
##            y = x.split("\\")
##            z = y[-1]
##            lenz = len(z + Suffix)
##            if lenz > 13:
##                x2 = x.rstrip(z)
##                overlen = 13 - lenz 
##                newz = z[0:overlen]
##                x2 = x2 + newz
##                Outputnames[num] = x2 + Suffix
##            else:
##                Outputnames[num] = x + Suffix
##            num = num + 1
##
##        change = str(Outputnames[0])
##        percent_change = str(Outputnames[1])
##        change_lt_zero = str(Outputnames[2])
##        change_gte_zero = str(Outputnames[3])
##        pchange_lt_zero = str(Outputnames[4])
##        pchange_gte_zero = str(Outputnames[5])
##            
##    except:
##        gp.AddError ("Error validating output filenames: " + gp.GetMessages(2))
##        raise Exception


    # Absolute change between two scenarios
    try:

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
        
        if do_change:
            gp.AddMessage ("\nCalculating change...")
            gp.Minus_sa(scenario1, scenario2, change)
            gp.AddMessage("\tCreated change output file: \n\t" + str(change))
    except:
        gp.AddError ("Error calculating change: " + gp.GetMessages(2))
        raise Exception


    # Percent change between two scenarios
    try:
        if do_percent_change:
            gp.AddMessage ("\nCalculating percent change...")
            gp.SingleOutputMapAlgebra_sa("((" + scenario1 + " - " + scenario2 + ") / " + scenario2 + ") * 100", percent_change)
            gp.AddMessage("\tCreated percent change output file: \n\t" + str(percent_change))
    except:
        gp.AddError ("Error calculating percent change: " + gp.GetMessages(2))
        raise Exception


    # If requested, split percent change results into
    # positive and negative valued rasters
    try:
        if do_split:
            if do_change:
                gp.AddMessage ("\nSplitting change output...")
                gp.SingleOutputMapAlgebra_sa("CON(" + change + " < 0, " + change + ")", change_lt_zero)
                gp.SingleOutputMapAlgebra_sa("CON(" + change + " >= 0, " + change + ")", change_gte_zero)
                gp.AddMessage("\tCreated change less than zero output file: \n\t" + str(change_lt_zero))
                gp.AddMessage("\tCreated change greater than zero output file: \n\t" + str(change_gte_zero))
            if do_percent_change:
                gp.AddMessage ("\nSplitting percent change output...")
                gp.SingleOutputMapAlgebra_sa("CON(" + percent_change + " < 0, " + percent_change + ")", pchange_lt_zero)
                gp.SingleOutputMapAlgebra_sa("CON(" + percent_change + " >= 0, " + percent_change + ")", pchange_gte_zero)
                gp.AddMessage("\tCreated percent change less than zero output file: \n\t" + str(pchange_lt_zero))
                gp.AddMessage("\tCreated percent change greater than zero output file: \n\t" + str(pchange_gte_zero))
    except:
        gp.AddError ("Error splitting change rasters: " + gp.GetMessages(2))
        raise Exception


# Throws an error in 9.3.1 - no time to figure out why right now

    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Output\\Calculate_Change_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("CALCULATE CHANGE\n")
        parafile.writelines("________________\n\n")

        for para in parameters:
            parafile.writelines(para + "\n")
            parafile.writelines("\n")
        parafile.close()
    except:
        gp.AddError ("Error creating parameter file:" + gp.GetMessages(2))
        raise Exception


    # Clean up temporary files
    gp.AddMessage("\nCleaning up temporary files...\n")
    try:
        gp.Delete_management(interws)
    except:
        gp.AddError("Error cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception



except:
    gp.AddError ("Error running script")
    raise Exception
