# ---------------------------------------------------------------------------
# Calculate_change.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
# 
# Last edit: 3/8/2010 
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

        # Shapefile of area(s) to calculate total change for (optional)
        zones = gp.GetParameterAsText(3)
        parameters.append("Zones: " + zones)

        # Zone field uniquely identifying each zone
        zones_id_field = gp.GetParameterAsText(4)
        parameters.append("Zones ID field: " + zones_id_field)

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(5)
        parameters.append("Suffix: " + Suffix)

        if (Suffix == "") or (Suffix == string.whitespace) or (Suffix == "#"):
            Suffix = ""
        else:
            Suffix = "_" + Suffix

        # Calculate change?
        do_change = gp.GetParameterAsText(6)
        if do_change =='true':
            do_change = True
            parameters.append("Calculate change: Yes")
        else:
            do_change = False
            parameters.append("Calculate change: No")

        # Calculate percent change?
        do_percent_change = gp.GetParameterAsText(7)
        if do_percent_change =='true':
            do_percent_change = True
            parameters.append("Calculate percent change: Yes")
        else:
            do_percent_change = False
            parameters.append("Calculate percent change: No")

        # Split the results into positive and negative rasters?
        do_split = gp.GetParameterAsText(8)
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

        # Calculate change over zones?
        do_zones = gp.GetParameterAsText(9)
        if do_zones =='true':
            if (zones != "") and (zones != string.whitespace) and (zones != "#"):
                if (zones_id_field != "") and (zones_id_field != string.whitespace) and (zones_id_field != "#"):
                    do_zones = True
                    parameters.append("Calculate change by zones: Yes")
                else:
                    gp.AddError("\nError: If change is to be calculated by zones, a Zone ID field must be specified.\n")
                    raise Exception
            else:
                gp.AddError("\nError: If change is to be calculated by zones, a Zones shapefile must be specified.\n")
                raise Exception
        else:
            do_zones = False
            parameters.append("Calculate change by zones: No")
            
    except:
        gp.AddError("Error in input arguments: " + gp.GetMessages(2))
        raise Exception


    # Check and create output folders
    try:
        thefolders=["Post_Process", "Intermediate"]
        for folder in thefolders:
            if not gp.Exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        gp.AddError("Error creating folders: " + gp.GetMessages())
        raise Exception


    # Output files
    try:
        # Base output/working directories
        postprocws = gp.workspace + os.sep + "Post_Process" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        
        # Intermediate variables
        x100 = "100"
        scen1_zstat_sum = interws + "scen1_zs_sum"
        scen2_zstat_sum = interws + "scen2_zs_sum"
        zone_sum_change_table = interws + "zone_sum_change.dbf"
        zone_sum_pchange_table = interws + "zone_sum_pchange.dbf"
        tmp_zone_output_table = interws + "tmp_zone_change.dbf"

        # Output field names
        change_sum_field = "change"
        pchange_sum_field = "pchange"
                
        # Output layers
        change = postprocws + "change"
        percent_change =  postprocws + "pchange"
        pchange_lt_zero = postprocws + "pch_lt0"
        pchange_gte_zero = postprocws + "pch_gte0"
        change_lt_zero = postprocws + "ch_lt0"
        change_gte_zero = postprocws + "ch_gte0"
        zone_sum_change = postprocws + "zone_ch"
        zone_sum_pchange = postprocws + "zone_pch"
        zone_output_table = postprocws + "zone_change" + Suffix + ".dbf"
        zone_output_table_name = "zone_change" + Suffix + ".dbf"

        
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [change, percent_change, change_lt_zero, change_gte_zero, pchange_lt_zero, pchange_gte_zero, zone_sum_change, zone_sum_pchange]           
        num = 0
        for x in Outputnames:
            y = x.split("\\")
            z = y[-1]
            lenz = len(z + Suffix)
            if lenz > 13:
                x2 = x.rstrip(z)
                overlen = 13 - lenz 
                newz = z[0:overlen]
                x2 = x2 + newz
                Outputnames[num] = x2 + Suffix
            else:
                Outputnames[num] = x + Suffix
            num = num + 1

        change = str(Outputnames[0])
        percent_change = str(Outputnames[1])
        change_lt_zero = str(Outputnames[2])
        change_gte_zero = str(Outputnames[3])
        pchange_lt_zero = str(Outputnames[4])
        pchange_gte_zero = str(Outputnames[5])
        zone_sum_change = str(Outputnames[6])
        zone_sum_pchange = str(Outputnames[7])
            
    except:
        gp.AddError ("Error validating output filenames: " + gp.GetMessages(2))
        raise Exception


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
            # %change not x100
##            gp.SingleOutputMapAlgebra_sa("((" + scenario1 + " - " + scenario2 + ") / " + scenario2 + ")", percent_change)
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

    # Table and raster with aggregated change results
    # Sum over whole area(s) of interest, calculate change from result
    # User enters a polygon layer to sum over - can be whole watershed, sub-watersheds, political districts, whatever
    try:
        if do_zones:
            gp.AddMessage("\nCalculating change by zones...")
            gp.ZonalStatistics_sa(zones, zones_id_field, scenario1, scen1_zstat_sum, "SUM")
            gp.ZonalStatistics_sa(zones, zones_id_field, scenario2, scen2_zstat_sum, "SUM")
            # Raster of absolute change by zones
            gp.Minus_sa(scen1_zstat_sum, scen2_zstat_sum, zone_sum_change)
            gp.AddMessage("\tCreated change by zones raster output file: \n\t" + str(zone_sum_change))
            # Raster of percent change by zones
            gp.SingleOutputMapAlgebra_sa("(" + zone_sum_change + " / " + scen2_zstat_sum + ") * 100", zone_sum_pchange)
            gp.AddMessage("\tCreated percent change by zones raster output file: \n\t" + str(zone_sum_pchange))

            # Turn them into a table
            gp.AddMessage("\nCreating output tables...")
            gp.ZonalStatisticsAsTable_sa(zones, zones_id_field, zone_sum_change, zone_sum_change_table, "DATA")
            gp.ZonalStatisticsAsTable_sa(zones, zones_id_field, zone_sum_pchange, zone_sum_pchange_table, "DATA")
            # Only one value per zone, so use the mean
            gp.AddField(zone_sum_change_table, change_sum_field, "FLOAT")
            gp.CalculateField_management(zone_sum_change_table, change_sum_field, "[MEAN]", "VB")
            gp.AddField(zone_sum_pchange_table, pchange_sum_field, "FLOAT")
            gp.CalculateField_management(zone_sum_pchange_table, pchange_sum_field, "[MEAN]", "VB")
            # Join change and pchange into one table
            gp.MakeTableView(zone_sum_change_table, "zone_sum_change_table_view")
            gp.MakeTableView(zone_sum_pchange_table, "zone_sum_pchange_table_view")
            gp.AddJoin("zone_sum_change_table_view", zones_id_field, "zone_sum_pchange_table_view", zones_id_field)
            gp.CopyRows("zone_sum_change_table_view", zone_output_table)
            # Of course, Arc has to pointlessly change my field name, so take it back
            gp.AddField(zone_output_table, pchange_sum_field, "FLOAT")
            gp.CalculateField_management(zone_output_table, pchange_sum_field, "[zone_sum_1]", "VB")
            # No easy way to clean up the unwanted fields, so do a workaround
            dropFields = ""
            fieldList = gp.ListFields(zone_output_table)
            field = fieldList.Next()
            keep_list = ["OID", zones_id_field, "AREA", change_sum_field, pchange_sum_field]
            while (field):
                if field.Name not in keep_list:
                    dropFields +=  field.Name + ";"
                field = fieldList.Next()
            dropFields = dropFields[:-1]
            gp.DeleteField_management(zone_output_table, "\"" + dropFields + "\"")
        
            gp.AddMessage("\tCreated change by zones table output file: \n\t" + str(zone_output_table))           
            
    except:
        gp.AddError ("Error calculating change by zones: " + gp.GetMessages(2))
        raise Exception


    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Post_Process\\PP_Calculate_Change_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("POST-PROCESSING - CALCULATE CHANGE\n")
        parafile.writelines("__________________________________\n\n")
        for para in parameters:
            parafile.writelines(para + "\n")
            parafile.writelines("\n")
        parafile.close()
    except:
        gp.AddError ("\nError creating parameter file:" + gp.GetMessages(2))
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
