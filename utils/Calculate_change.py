# ---------------------------------------------------------------------------
# Calculate_change.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
#
# Calculates absolute and/or percent change between scenario outputs
# Optionally splits the results into two rasters, one with positive
#   values and the other with negative, to make symbolizing easier
#
# NOTES, REMOVE:
# - pixel-based outputs: use old code, add mask for summarizing over whole area
# - sub-watershed outputs: use old code, shapefile for making table
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

        # Calculate change for whole area of interest?
        do_area = gp.GetParameterAsText(3)
        if do_area =='true':
            do_area = True
            parameters.append("Calculate change for whole area of interest: Yes")
        else:
            do_area = False
            parameters.append("Calculate change for whole area of interest: No")

        # If calculating change for whole area, need area shapefile
        area_mask = gp.GetParameterAsText(4)
        parameters.append("Area mask: " + area_mask)
        if do_area and (area_mask == "") or (area_mask == string.whitespace) or (area_mask == "#"):
            gp.AddError("\nError: If the area of interest is to be processed, a corresponding shapefile must be provided")
            raise Exception

        # Area of interest ID field - used in making output table
        area_id_field = gp.GetParameterAsText(5)
        parameters.append("Area of interest ID field: " + area_id_field)
        if do_area and (area_id_field == "") or (area_id_field == string.whitespace) or (area_id_field == "#"):
            gp.AddError("\nError: If the area of interest is to be processed, an ID field must be provided")
            raise Exception

        # Calculate output table per subwatershed?
        do_subwshed = gp.GetParameterAsText(6)
        if do_subwshed =='true':
            do_subwshed = True
            parameters.append("Calculate change for subwatersheds: Yes")
        else:
            do_subwshed = False
            parameters.append("Calculate change for subwatersheds: No")

        # If calculating change for subwatersheds, need shapefile
        subwsheds = gp.GetParameterAsText(7)
        parameters.append("Subwatersheds: " + area_mask)
        if do_subwshed and (subwsheds == "") or (subwsheds == string.whitespace) or (subwsheds == "#"):
            gp.AddError("\nError: If subwatersheds are to be processed, a subwatershed shapefile must be provided")
            raise Exception

        # Subwatershed ID field - used in making output table
        subwshed_id_field = gp.GetParameterAsText(8)
        parameters.append("Subwatershed ID field: " + subwshed_id_field)
        if do_subwshed and (subwshed_id_field == "") or (subwshed_id_field == string.whitespace) or (subwshed_id_field == "#"):
            gp.AddError("\nError: If subwatersheds are to be processed, a subwatershed ID field must be provided")
            raise Exception
            
        # Calculate change?
        do_change = gp.GetParameterAsText(9)
        if do_change =='true':
            do_change = True
            parameters.append("Calculate change: Yes")
        else:
            do_change = False
            parameters.append("Calculate change: No")

        # Calculate percent change?
        do_percent_change = gp.GetParameterAsText(10)
        if do_percent_change =='true':
            do_percent_change = True
            parameters.append("Calculate percent change: Yes")
        else:
            do_percent_change = False
            if not do_change:
                gp.AddError("\nError: Nothing to do. Please select \'Calculate change\' and/or \'Calculate percent change\'.")
                raise Exception
            parameters.append("Calculate percent change: No")

        # Split the results into positive and negative rasters?
        do_split = gp.GetParameterAsText(11)
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
        Suffix = gp.GetParameterAsText(12)
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
        postprocws = gp.workspace + os.sep + "Post_process" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        
        # Intermediate variables
        x100 = "100"
        change_zstat_table = interws + "change_zstat.dbf"
        pchange_zstat_table = interws + "pchange_zstat.dbf"
        
        # Output layers
        change = postprocws + "change" + Suffix + ".tif"
        percent_change =  postprocws + "percent_change" + Suffix + ".tif"
        pchange_lt_zero = postprocws + "percent_change_lt0" + Suffix + ".tif"
        pchange_gte_zero = postprocws + "percent_change_gte0" + Suffix + ".tif"
        change_lt_zero = postprocws + "change_lt0" + Suffix + ".tif"
        change_gte_zero = postprocws + "change_gte0" + Suffix + ".tif"


        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
        
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages(2))
        raise Exception



    # If requested, split percent change results into
    # positive and negative valued rasters
    def split(spraster):
        try:
            gp.AddMessage ("\nSplitting output...")
            gp.SingleOutputMapAlgebra_sa("CON(" + spraster + " < 0, " + spraster + ")", spraster_lt_zero)
            gp.SingleOutputMapAlgebra_sa("CON(" + spraster + " >= 0, " + spraster + ")", spraster_gte_zero)
            return(spraster_lt_zero, spraster_gte_zero)

        except:
            gp.AddError ("Error splitting change rasters: " + gp.GetMessages(2))
            raise Exception

    # Process inputs
    try:
        # Simple change (scenario1 - scenario2)
        if do_change:
            gp.AddMessage ("\nCalculating change...")
            gp.Minus_sa(scenario1, scenario2, change)
            gp.AddMessage("\tCreated change output file: \n\t" + str(change))

            if do_split:
                change_lt_zero, change_gte_zero = split(change)
                gp.AddMessage("\tCreated change less than zero output file: \n\t" + str(change_lt_zero))
                gp.AddMessage("\tCreated change greater than zero output file: \n\t" + str(change_gte_zero))

        # Percent change (((scenario1 - scenario2) / scenario2) * 100)
        if do_percent_change:
            gp.AddMessage("\nCalculating percent change...")
            gp.SingleOutputMapAlgebra_sa("((" + scenario1 + " - " + scenario2 + ") / " + scenario2 + ") * 100", percent_change)
            gp.AddMessage("\tCreated percent change output file: \n\t" + str(percent_change))

            if do_split:
                pchange_lt_zero, pchange_gte_zero = split(percent_change)
                gp.AddMessage("\tCreated percent change less than zero output file: \n\t" + str(pchange_lt_zero))
                gp.AddMessage("\tCreated percent change greater than zero output file: \n\t" + str(pchange_gte_zero))

        # Create table of change per subwatershed
        if do_subwsheds:
            gp.AddMessage("\nCreating change table for subwatersheds...")

            # Output table with change values per subwatershed
            change_subwshed_table_name = "change_subwatershed.dbf"
            change_subwshed_table = postprocws + change_subwshed_table_name            
            gp.CreateTable_management(postprocws, change_subwshed_table_name)

            # For subwatershed input, each subwatershed has the same value for each cell in it,
            # so do zonal stats and take the mean
                
            if do_change:
                gp.AddMessage("\n\tCalculating change per subwatershed...")
                gp.ZonalStatisticsAsTable_sa(subwsheds, subwshed_id_field, change, change_zstat_table, "DATA")
                gp.AddField(change_subwshed_table, "change", "double")
                
            if do_percent_change:
                gp.AddMessage("\n\tCalculating change per subwatershed...")
                gp.ZonalStatisticsAsTable_sa(subwsheds, subwshed_id_field, percent_change, pchange_zstat_table, "DATA")
                gp.AddField(change_subwshed_table, "pchange", "float")


        # Create table of change for whole area of interest
        if do_area:
            gp.AddMessage("\nCreating change table for whole area of interest...")

            # If subwatersheds are entered, need to sum them for whole area of interest
            if do_subwsheds:
                gp.AddMessage("\n\tSumming by subwatershed...")

            # If no subwatersheds, can just sum over whole area
            else:
                if do_change:
                    gp.ZonalStatisticsAsTable_sa(area_mask, subwshed_id_field, change, change_zstat_table, "DATA")
                    gp.AddField(change_subwshed_table, "change", "double")
                
                # Can't do zstat on pchange, need to sum over subwshed first, then get pchange
                if do_percent_change:
                    gp.AddField(change_subwshed_table, "pchange", "float")


           

    except:
        gp.AddError ("Error calculating change: " + gp.GetMessages(2))
        raise Exception


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


##    # Clean up temporary files
##    gp.AddMessage("\nCleaning up temporary files...\n")
##    try:
##        gp.Delete_management(interws)
##    except:
##        gp.AddError("Error cleaning up temporary files:  " + gp.GetMessages(2))
##        raise Exception


except:
    gp.AddError ("Error running script")
    raise Exception
