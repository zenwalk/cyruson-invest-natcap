# ---------------------------------------------------------------------------
# Calculate_change.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
#
# Calculates absolute and/or percent change between scenario inputs
# Optionally splits the results into two rasters, one with positive
#   values and the other with negative, to make symbolizing more accurate
#
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

        # Is the input data raster?
        do_raster = gp.GetParameterAsText(1)
        if do_raster =='true':
            do_raster = True
            parameters.append("Calculate change for rasters: Yes")
        else:
            do_raster = False
            parameters.append("Calculate change for rasters: No")

        # Raster Scenario 1
        scenario1_ras = gp.GetParameterAsText(2)
        parameters.append("Raster Scenario 1: " + scenario1_ras)
        if do_raster and (scenario1_ras == "") or (scenario1_ras == string.whitespace) or (scenario1_ras == "#"):
            gp.AddError("\nError: If rasters are to be processed, Scenario 1 and 2 rasters must be provided")
            raise Exception

        # Raster Scenario 2
        scenario2_ras = gp.GetParameterAsText(3)
        parameters.append("Raster Scenario 2: " + scenario2_ras)
        if do_raster and (scenario2_ras == "") or (scenario2_ras == string.whitespace) or (scenario2_ras == "#"):
            gp.AddError("\nError: If rasters are to be processed, Scenario 1 and 2 rasters must be provided")
            raise Exception
        
        # Is the input data by subwatershed?
        do_subwshed = gp.GetParameterAsText(4)
        if do_subwshed =='true':
            do_subwshed = True
            parameters.append("Calculate change for subwatersheds: Yes")
        else:
            do_subwshed = False
            parameters.append("Calculate change for subwatersheds: No")

        # If calculating change for subwatersheds, need shapefile
        subwsheds = gp.GetParameterAsText(5)
        parameters.append("Subwatersheds: " + subwsheds)
        if do_subwshed and (subwsheds == "") or (subwsheds == string.whitespace) or (subwsheds == "#"):
            gp.AddError("\nError: If subwatersheds are to be processed, a subwatershed shapefile must be provided")
            raise Exception

        # Subwatershed ID field - used in making output table
        subwshed_id_field = gp.GetParameterAsText(6)
        parameters.append("Subwatershed ID field: " + subwshed_id_field)
        if do_subwshed and (subwshed_id_field == "") or (subwshed_id_field == string.whitespace) or (subwshed_id_field == "#"):
            gp.AddError("\nError: If subwatersheds are to be processed, a subwatershed ID field must be selected")
            raise Exception
        
        # Calculate change within an area of interest?
        do_area = gp.GetParameterAsText(7)
        if do_area =='true':
            do_area = True
            parameters.append("Calculate change for whole area of interest: Yes")
        else:
            do_area = False
            parameters.append("Calculate change for whole area of interest: No")

        # This tool does not aggregate subwatershed data to the watershed level
        # Use watershed tables instead
        if do_subwshed and do_area:
            gp.AddError("\nError: This tool does not aggregate subwatershed data within an area of interest.")
            gp.AddError("\tTo process watershed data, use the \'watershed output table\' option.")
            raise Exception

        # If calculating change for whole area, need area shapefile
        area_mask = gp.GetParameterAsText(8)
        parameters.append("Area mask: " + area_mask)
        if do_area and (area_mask == "") or (area_mask == string.whitespace) or (area_mask == "#"):
            gp.AddError("\nError: If the area of interest is to be processed, a corresponding shapefile must be provided")
            raise Exception

        # Area of interest ID field - used in making output table
        area_id_field = gp.GetParameterAsText(9)
        parameters.append("Area of interest ID field: " + area_id_field)
        if do_area and (area_id_field == "") or (area_id_field == string.whitespace) or (area_id_field == "#"):
            gp.AddError("\nError: If the area of interest is to be processed, an ID field must be selected")
            raise Exception

        # Is the input data watershed output tables?
        do_table = gp.GetParameterAsText(10)
        if do_table =='true':
            do_table = True
            parameters.append("Calculate change for watershed tables: Yes")
        else:
            do_table = False
            parameters.append("Calculate change for watershed tables: No")

        # Table Scenario 1
        scenario1_table = gp.GetParameterAsText(11)
        parameters.append("Scenario 1 table: " + scenario1_table)
        if do_table and (scenario1_table == "") or (scenario1_table == string.whitespace) or (scenario1_table == "#"):
            gp.AddError("\nError: If watershed tables are to be processed, Scenario 1 and 2 tables must be provided")
            raise Exception

        # Table Scenario 2
        scenario2_table = gp.GetParameterAsText(12)
        parameters.append("Scenario 2 table: " + scenario2_table)
        if do_table and (scenario2_table == "") or (scenario2_table == string.whitespace) or (scenario2_table == "#"):
            gp.AddError("\nError: If watershed tables are to be processed, Scenario 1 and 2 tables must be provided")
            raise Exception

        # Watershed ID field
        wshed_table_id_field = gp.GetParameterAsText(13)
        parameters.append("Watershed table ID field: " + wshed_table_id_field)
        if do_table and (wshed_table_id_field == "") or (wshed_table_id_field == string.whitespace) or (wshed_table_id_field == "#"):
            gp.AddError("\nError: If watershed tables are to be processed, a watershed ID field must be selected")
            raise Exception

        # Field to calculate change for
        wshed_table_data_field = gp.GetParameterAsText(14)
        parameters.append("Watershed table data field: " + wshed_table_data_field)
        if do_table and (wshed_table_data_field == "") or (wshed_table_data_field == string.whitespace) or (wshed_table_data_field == "#"):
            gp.AddError("\nError: If watershed tables are to be processed, a data field must be selected")
            raise Exception        
        
        # Calculate change?
        do_change = gp.GetParameterAsText(15)
        if do_change =='true':
            do_change = True
            parameters.append("Calculate change: Yes")
        else:
            do_change = False
            parameters.append("Calculate change: No")

        # Calculate percent change?
        do_percent_change = gp.GetParameterAsText(16)
        if do_percent_change == 'true':
            do_percent_change = True
            parameters.append("Calculate percent change: Yes")
        else:
            do_percent_change = False
            if not do_change:
                gp.AddError("\nError: Nothing to do. Please select \'Calculate change\' and/or \'Calculate percent change\'.")
                raise Exception
            parameters.append("Calculate percent change: No")

        # Split the results into positive and negative rasters?
        do_split = gp.GetParameterAsText(17)
        if do_split == 'true':
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
        Suffix = gp.GetParameterAsText(18)
        parameters.append("Suffix: " + Suffix)

        if (Suffix == "") or (Suffix == string.whitespace) or (Suffix == "#"):
            Suffix = ""
        else:
            Suffix = "_" + Suffix
            
    except:
        gp.AddError("\nError in input arguments: " + gp.GetMessages(2))
        raise Exception



    # Check and create output folders
    try:
        thefolders=["Post_process", "Intermediate"]
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
        
    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Set the Geoprocessing environment
    try:
        install_info = gp.GetInstallInfo("desktop")
        
        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
    except:
        gp.AddError( "\nError setting geoprocessing environment: " + gp.GetMessages(2))
        raise Exception


    # If requested, split percent change results into
    # positive and negative valued rasters
    def split(spraster, in_lt_zero, in_gte_zero):
        try:
        
            gp.AddMessage ("\n\tSplitting output...")
            gp.SingleOutputMapAlgebra_sa("CON(" + spraster + " < 0, " + spraster + ")", in_lt_zero)
            gp.SingleOutputMapAlgebra_sa("CON(" + spraster + " >= 0, " + spraster + ")", in_gte_zero)
            return(in_lt_zero, in_gte_zero)

        except:
            gp.AddError ("\nError splitting change rasters: " + gp.GetMessages(2))
            raise Exception

    # Process inputs
    try:
        # Simple change (scenario1_ras - scenario2_ras)
        if do_change and (do_raster or do_subwshed):
            try:
                change = postprocws + "change" + Suffix + ".tif"
                        
                gp.AddMessage ("\nCalculating change...")
                gp.Minus_sa(scenario1_ras, scenario2_ras, change)
                gp.AddMessage("\n\tCreated change output file: \n\t" + str(change))

                if do_split:
                    change_lt_zero = postprocws + "change_lt0" + Suffix + ".tif"
                    change_gte_zero = postprocws + "change_gte0" + Suffix + ".tif"
        
                    change_lt_zero, change_gte_zero = split(change, change_lt_zero, change_gte_zero)
                    gp.AddMessage("\t\tCreated change less than zero output file: \n\t" + str(change_lt_zero))
                    gp.AddMessage("\t\tCreated change greater than zero output file: \n\t" + str(change_gte_zero))

            except:
                gp.AddError ("\nError calculating change: " + gp.GetMessages(2))
                raise Exception

        # Percent change (((scenario1_ras - scenario2_ras) / scenario2_ras) * 100)
        if do_percent_change and (do_raster or do_subwshed):
            try:
                percent_change =  postprocws + "percent_change" + Suffix + ".tif"
                scenario2_ras_no0 = interws + "scen2ras_no0"
                
                gp.AddMessage("\nCalculating percent change...")
                # Can't divide by zero so add a very small number to scenario2_ras where there are zeros
                gp.SingleOutputMapAlgebra_sa("CON(" + scenario2_ras + " == 0, .0000001, " + scenario2_ras + ")", scenario2_ras_no0)
                gp.SingleOutputMapAlgebra_sa("((" + scenario1_ras + " - " + scenario2_ras_no0 + ") / " + scenario2_ras_no0 + ") * 100", percent_change)
                gp.AddMessage("\n\tCreated percent change output file: \n\t" + str(percent_change))

                if do_split:
                    pchange_lt_zero = postprocws + "percent_change_lt0" + Suffix + ".tif"
                    pchange_gte_zero = postprocws + "percent_change_gte0" + Suffix + ".tif"
                    
                    pchange_lt_zero, pchange_gte_zero = split(percent_change, pchange_lt_zero, pchange_gte_zero)
                    gp.AddMessage("\t\tCreated percent change less than zero output file: \n\t" + str(pchange_lt_zero))
                    gp.AddMessage("\t\tCreated percent change greater than zero output file: \n\t" + str(pchange_gte_zero))
            except:
                gp.AddError ("\nError calculating percent change: " + gp.GetMessages(2))
                raise Exception

        # Table of change per subwatershed
        if do_subwshed:
            try:
                gp.AddMessage("\nCreating change table for subwatersheds...")

                # Output table with change values per subwatershed
                change_subwshed_table_name = "change_subwatershed" + Suffix + ".dbf"
                change_subwshed_table = postprocws + change_subwshed_table_name            
                gp.CreateTable_management(postprocws, change_subwshed_table_name)
                gp.AddField(change_subwshed_table, "subws_id", "long")
                gp.DeleteField_management(change_subwshed_table, "Field1")


                # Zonal stats field name has changed in Arc 10
                if (install_info["Version"] == "10.0" or install_info["Version"] == "10.1"):
                    zstat_id_field = subwshed_id_field
                else:
                    zstat_id_field = "VALUE"
                    
                if do_change:
                    gp.AddMessage("\n\tCalculating change per subwatershed...")

                    change_zstat_table = interws + "change_zstat_sws.dbf"
                    
                    # Use per-pixel change calculated above.  Each subwatershed will have the same value
                    # for each cell in it, so do zonal stats and take the mean
                    gp.ZonalStatisticsAsTable_sa(subwsheds, subwshed_id_field, change, change_zstat_table, "DATA")
                    change_subws_zstat_rows = gp.SearchCursor(change_zstat_table)
                    gp.AddField(change_subwshed_table, "change", "double")
                    # Order zonal stats table by subwastershed id
                    change_zstat_rows = gp.SearchCursor(change_zstat_table, "", "", "", zstat_id_field + " A")
                    change_zstat_row = change_zstat_rows.Reset
                    change_zstat_row = change_zstat_rows.Next()

                    # Add change values to output table
                    change_subws_rows = gp.InsertCursor(change_subwshed_table)

                    while(change_zstat_row):
                        # Get mean value for this subwatershed
                        change_mean = float(change_zstat_row.getValue("MEAN"))
                        new_row = change_subws_rows.NewRow()
                        new_row.setValue("subws_id", change_zstat_row.getValue(zstat_id_field))
                        new_row.setValue("change", change_mean)
                        change_subws_rows.InsertRow(new_row)

                        change_zstat_row = change_zstat_rows.Next()

                    del new_row, change_subws_rows

                    
                if do_percent_change:
                    gp.AddMessage("\n\tCalculating percent change per subwatershed...")

                    pchange_zstat_table = interws + "pchange_zstat_sws.dbf"
                    
                    gp.ZonalStatisticsAsTable_sa(subwsheds, subwshed_id_field, percent_change, pchange_zstat_table, "DATA")
                    gp.AddField(change_subwshed_table, "pchange", "float")
                    
                    # Order zonal stats table ascending by subwastershed id
                    pchange_zstat_rows = gp.SearchCursor(pchange_zstat_table, "", "", "", zstat_id_field + " A")
                    pchange_zstat_row = pchange_zstat_rows.Reset
                    pchange_zstat_row = pchange_zstat_rows.Next()

                    # Add percent change values to output table

                    # If change was already added, update existing entries
                    if do_change:
                        change_subws_rows = gp.UpdateCursor(change_subwshed_table)
                        change_subws_row = change_subws_rows.Reset
                        change_subws_row = change_subws_rows.Next()
                        
                    # If not, add new entries
                    else:
                        change_subws_rows = gp.InsertCursor(change_subwshed_table)
                        

                    while(pchange_zstat_row):
                        # Get mean value for this subwatershed
                        # All pixel values within a subwatershed are the same, so take the mean
                        pchange_mean = float(pchange_zstat_row.getValue("MEAN"))

                        if not do_change:
                            # Add new entry
                            new_row = change_subws_rows.NewRow()
                            new_row.setValue("subws_id", pchange_zstat_row.getValue(zstat_id_field))
                            new_row.setValue("pchange", pchange_mean)
                            change_subws_rows.InsertRow(new_row)

                        else:
                            # Find existing entry

                            while(int(change_subws_row.getValue("subws_id")) <> int(pchange_zstat_row.getValue(zstat_id_field))):
                                change_subws_row = change_subws_rows.Next()

                            # Update entry with percent change value
                            change_subws_row.setValue("pchange", pchange_mean)
                            change_subws_rows.UpdateRow(change_subws_row)

                            change_subws_rows = gp.UpdateCursor(change_subwshed_table)
                            change_subws_row = change_subws_rows.Reset
                            change_subws_row = change_subws_rows.Next()

                        pchange_zstat_row = pchange_zstat_rows.Next()

                    if not do_change:
                        del change_subws_rows, new_row
                    else:
                        del change_subws_rows, change_subws_row

                gp.AddMessage("\n\tCreated change table for subwatersheds: \n\t" + change_subwshed_table)
                
            except:
                gp.AddError ("\nError creating change table for subwatersheds: " + gp.GetMessages(2))
                raise Exception
            
        # Create table of change for whole area of interest
        if do_area:
            try:
                
                gp.AddMessage("\nCreating change table within the area of interest...")

                # Create output table
                change_area_table_name = "change_area" + Suffix + ".dbf"
                change_area_table = postprocws + change_area_table_name            
                gp.CreateTable_management(postprocws, change_area_table_name)
                gp.AddField(change_area_table, "area_id", "long")
                gp.DeleteField_management(change_area_table, "Field1")

                # Zonal stats field name has changed in Arc 10
                if (install_info["Version"] == "10.0" or install_info["Version"] == "10.1"):
                    zstat_id_field = area_id_field
                else:
                    zstat_id_field = "VALUE"
             

                # For pixel data, can just calculate over whole area's pixels                 
                if do_change:

                    gp.AddMessage("\n\tCalculating change...")

                    gp.AddField(change_area_table, "change", "double")

                    change_zstat_table_area = interws + "change_zstat_table_area.dbf"
                    
                    gp.ZonalStatisticsAsTable_sa(area_mask, area_id_field, change, change_zstat_table_area, "DATA")
                    
                    # Order zonal stats table ascending by area id
                    change_zstat_area_rows = gp.SearchCursor(change_zstat_table_area, "", "", "", zstat_id_field + " A")
                    change_zstat_area_row = change_zstat_area_rows.Reset
                    change_zstat_area_row = change_zstat_area_rows.Next()

                    # Add change values to output table
                    change_area_rows = gp.InsertCursor(change_area_table)
                    
                    while(change_zstat_area_row):
                        # Get total change value for the area
                        change_sum = float(change_zstat_area_row.getValue("SUM"))
                        new_row = change_area_rows.NewRow()
                        new_row.setValue("area_id", change_zstat_area_row.getValue(zstat_id_field))
                        new_row.setValue("change", change_sum)
                        change_area_rows.InsertRow(new_row)

                        change_zstat_area_row = change_zstat_area_rows.Next()                        

                    del change_zstat_area_row, change_zstat_area_rows, change_area_rows, new_row
                
                # For percent change, need to sum over each scenario first, then calculate pchange
                if do_percent_change:

                    gp.AddMessage("\n\tCalculating percent change...")

                    s1_zstat = interws + "s1_zstat"
                    s2_zstat = interws + "s2_zstat"
                    s2_zstat_no0 = interws + "s2_zstat_no0"
                    pchange_area = interws + "pchange_area" + Suffix + ".tif"
                    pchange_zstat_table_area = interws + "pchange_zstat_table_area.dbf"

                    gp.AddField(change_area_table, "pchange", "float")

                    # Calculate percent change over whole area
                    gp.ZonalStatistics_sa(area_mask, area_id_field, scenario1_ras, s1_zstat, "SUM", "DATA")
                    gp.ZonalStatistics_sa(area_mask, area_id_field, scenario2_ras, s2_zstat, "SUM", "DATA")
                    # Can't divide by zero so add a very small number to scenario2_ras where there are zeros
                    gp.SingleOutputMapAlgebra_sa("CON(" + s2_zstat + " == 0, .0000001, " + s2_zstat + ")", s2_zstat_no0)
                    gp.SingleOutputMapAlgebra_sa("((" + s1_zstat + " - " + s2_zstat_no0 + ") / " + s2_zstat_no0 + ") * 100", pchange_area)

                    # Add percent change value to output table
                    
                    gp.ZonalStatisticsAsTable_sa(area_mask, area_id_field, pchange_area, pchange_zstat_table_area, "DATA")
                    pchange_zstat_table_area_rows = gp.SearchCursor(pchange_zstat_table_area)
                    # Order zonal stats table ascending by area id
                    pchange_zstat_area_rows = gp.SearchCursor(pchange_zstat_table_area, "", "", "", zstat_id_field + " A")
                    pchange_zstat_area_row = pchange_zstat_area_rows.Reset
                    pchange_zstat_area_row = pchange_zstat_area_rows.Next()

                    # Add percent change values to output table

                    # If change was already added, update existing entries
                    
                    if do_change:
                        change_area_rows = gp.UpdateCursor(change_area_table)
                        change_area_row = change_area_rows.Reset
                        change_area_row = change_area_rows.Next()
                        
                    # If not, add new entries
                    else:
                        change_area_rows = gp.InsertCursor(change_area_table)                        

                    while(pchange_zstat_area_row):
                        # Use mean for percent change across area
                        pchange_mean = float(pchange_zstat_area_row.getValue("MEAN"))

                        if not do_change:
                            # Add new entry
                            new_row = change_area_rows.NewRow()
                            new_row.setValue("area_id", pchange_zstat_area_row.getValue(zstat_id_field))
                            new_row.setValue("pchange", pchange_mean)
                            change_area_rows.InsertRow(new_row)
                            del new_row

                        else:
                            while(int(change_area_row.getValue("area_id")) <> int(pchange_zstat_area_row.getValue(zstat_id_field))):
                                change_area_row = change_area_rows.Next()

                            # Update entry with percent change value
                            change_area_row.setValue("pchange", pchange_mean)
                            change_area_rows.UpdateRow(change_area_row)
                            
                            change_area_rows = gp.UpdateCursor(change_area_table)
                            change_area_row = change_area_rows.Reset
                            change_area_row = change_area_rows.Next()

                        pchange_zstat_area_row = pchange_zstat_area_rows.Next()

#                            del change_area_row, change_area_rows, pchange_zstat_area_row, pchange_zstat_area_rows

                gp.AddMessage("\n\tCreated area of interest change table:\n\t" + str(change_area_table))
                    
            except:
                gp.AddError ("\nError creating change table within the area of interest: " + gp.GetMessages(2))
                raise Exception

        # If watershed table output is entered, process it
        try:
            if do_table:
                gp.AddMessage("\nCreating watershed change table...")

                # Create output table
                change_wshed_table_name = "change_watershed" + Suffix + ".dbf"
                change_wshed_table = postprocws + change_wshed_table_name            
                gp.CreateTable_management(postprocws, change_wshed_table_name)
                gp.AddField(change_wshed_table, wshed_table_id_field, "long")
                gp.DeleteField_management(change_wshed_table, "Field1")

                if do_change:
                    gp.AddField(change_wshed_table, "change", "double")
                if do_percent_change:
                    gp.AddField(change_wshed_table, "pchange", "double")

                change_wshed_table_rows = gp.InsertCursor(change_wshed_table)
                
                # Loop through input tables, calculating change for each watershed
                scen1_table_rows = gp.SearchCursor(scenario1_table, "", "", "", wshed_table_id_field + " A")
                scen1_table_row = scen1_table_rows.Reset
                scen1_table_row = scen1_table_rows.Next()

                while scen1_table_row:
                    # Loop through second scenario's table
                    scen2_table_rows = gp.SearchCursor(scenario2_table, "", "", "", wshed_table_id_field + " A")
                    scen2_table_row = scen2_table_rows.Reset
                    scen2_table_row = scen2_table_rows.Next()

                    while scen2_table_row:
                        # Find matching watershed ID in both input tables
                        if (int(scen2_table_row.getValue(wshed_table_id_field)) == int(scen1_table_row.getValue(wshed_table_id_field))):
                            new_row = change_wshed_table_rows.NewRow()
                            # Set watershed ID in output table
                            new_row.setValue(wshed_table_id_field, int(scen2_table_row.getValue(wshed_table_id_field)))

                            if do_change:
                                # Calculate change and add to output table
                                wshed_change = float(scen1_table_row.getValue(wshed_table_data_field)) - float(scen2_table_row.getValue(wshed_table_data_field))
                                new_row.setValue("change", float(wshed_change))
                            if do_percent_change:
                                # Calculate percent change and add to output table
                                scen1_data = float(scen1_table_row.getValue(wshed_table_data_field))
                                scen2_data = float(scen2_table_row.getValue(wshed_table_data_field))
                                # Can't divide by zero, so add a very small number
                                if scen2_data == 0:
                                    scen2_data = .0000001
                                wshed_pchange = ((scen1_data - scen2_data) / scen2_data) * 100
                                new_row.setValue("pchange", float(wshed_pchange))

                            change_wshed_table_rows.InsertRow(new_row)
                            
                        scen2_table_row = scen2_table_rows.Next()
                        
                    scen1_table_row = scen1_table_rows.Next()
                    
                gp.AddMessage("\n\tCreated watershed change table: " + change_wshed_table)

        except:
                gp.AddError ("\nError creating watershed change table: " + gp.GetMessages(2))
                raise Exception    

    except:
        gp.AddError ("Error calculating change: " + gp.GetMessages(2))
        raise Exception


    # Write input parameters to an output file for user reference
    try:
        gp.AddMessage("\nCreating parameter log file...")
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(postprocws + "Calculate_Change_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
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
