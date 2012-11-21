# ---------------------------------------------------------------------------
# Sum_To_Threshold.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
# 
# Last edit: 3/26/2010 
#
# Slices an input grid into small pieces, then sums them
#   one at a time until a threshold is reached.
# Outputs a grid of the selected cells
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

        # Score grid
        score_grid = gp.GetParameterAsText(1)
        parameters.append("Score grid: " + score_grid)

        # Cost grid
        cost_grid = gp.GetParameterAsText(2)
        parameters.append("Cost grid: " + cost_grid)

        # New activity grid
        na_grid = gp.GetParameterAsText(3)
        parameters.append("New activity grid: " + na_grid)

        # Original land use grid
        lu_grid = gp.GetParameterAsText(4)
        parameters.append("Land use grid: " + lu_grid)

        # Watershed mask
        wshed_mask = gp.GetParameterAsText(5)
        parameters.append("Watershed mask: " + wshed_mask)

        # Threshold
        threshold = float(gp.GetParameterAsText(6))
        parameters.append("Threshold: " + str(threshold))

        # Number of slices
        num_slices = int(gp.GetParameterAsText(7))
        parameters.append("Number of slices: " + str(num_slices))

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(8)
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
        sliced_grid = interws + "slices"
        slice_sum_table = interws + "slice_sum.dbf"
        slice_table_keep = interws + "slice_keep_" + Suffix + ".dbf"
        slice_join = interws + "slice_join"
        slice_keep = interws + "slice_keep"
        slice_extract = interws + "slice_ext"
        slice_ext_table = interws + "slice_ext_table.dbf"
        slice_grid_keep = interws + "sl_gr_keep"
        na_selected_null = interws + "nasel_null"
        
        # Output layers
        score_selected = postprocws + "sc_sel"
        na_selected = postprocws + "na_sel"
        final_landuse = postprocws + "new_lu"
        
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [score_selected, na_selected, final_landuse]           
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

        score_selected = str(Outputnames[0])
        na_selected = str(Outputnames[1])
        final_landuse = str(Outputnames[2])
            
    except:
        gp.AddError ("Error validating output filenames: " + gp.GetMessages(2))
        raise Exception


    try:

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws

        # Set mask to watershed
        gp.Mask = wshed_mask

        # Time this process for user reference
        start_time = time.clock()

        # Slice input grid into input number of slices
        gp.AddMessage("\nSlicing grid...")
        gp.Slice_sa(score_grid, sliced_grid, num_slices, "EQUAL_INTERVAL")

        # Make updates to a copy of the slice file
        gp.AddMessage("copy raster")
        gp.CopyRaster_management(sliced_grid, slice_grid_keep)

        # Add a field denoting whether or not the slice is selected as
        # part of the sum - set initially to zero for not selected
        gp.AddMessage("add field")
        gp.AddField_management(slice_grid_keep, "keep", "SHORT")
        gp.AddMessage("calculate field")
        gp.CalculateField_management(slice_grid_keep, "keep", "0")

        # Largest slice values have highest Value
        slrows = gp.SearchCursor(sliced_grid, "", "", "", "VALUE D")
        slrow = slrows.Next()

        psum = 0
        count = 0

        gp.AddMessage("\nLooping...")

        while (slrow):
            ##        Extract slice by highest (current) Value
##            gp.AddMessage("extract by attr")
            count += 1
            gp.AddMessage("count = " + str(count))
            gp.AddMessage("slice value = " + str(slrow.getValue("VALUE")))
            where_clause = "\"VALUE\" = " + str(slrow.getValue("VALUE"))
            gp.ExtractByAttributes_sa(sliced_grid, where_clause, slice_extract)

    ##        zonal sum of cost grid with extracted slice
            gp.AddMessage("zonal stats")
            gp.ZonalStatisticsAsTable_sa(slice_extract, "VALUE", cost_grid, slice_ext_table, "DATA")
            gp.AddMessage("search cursor zrows")
            zrows = gp.SearchCursor(slice_ext_table)
            zrow = zrows.Next()
        

            # Find matching slice in keep file
            # This is a workaround because UpdateCursor fails when using sort
            gp.AddMessage("update cursor urows")
            urows = gp.UpdateCursor(slice_grid_keep)
            urow = urows.Reset
            urow = urows.Next()

            gp.AddMessage("while")
            while (int(slrow.getValue("VALUE")) <> int(urow.getValue("VALUE"))):
                urow = urows.Next()

            gp.AddMessage("after while")

            # Skip over any of the largest point values that are larger than the threshold
##            gp.AddMessage("zrow value = " + str(zrow.getValue("SUM")))
            if psum == 0 and zrow.getValue("SUM") > threshold:
                gp.AddMessage("skipping value " + str(zrow.getValue("SUM")))
                slrow = slrows.Next()
                del zrow, zrows
            # Hit the threshold exactly
            elif (psum + zrow.getValue("SUM") == threshold):
                psum += zrow.getValue("SUM")
                urow.setValue("keep", "1")
                urows.UpdateRow(urow)
                gp.AddMessage("hit the threshold exactly")
                break
            # Went over the threshold, so don't add current value
            elif (psum + zrow.getValue("SUM") > threshold):
                gp.AddMessage("went over threshold")
                # test - looking for smaller values to fill in threshold gap
                slrow = slrows.Next()
                del zrow, zrows
##                break
            # Still under the threshold, add current value and loop again
            else:
                
                psum += zrow.getValue("SUM")
                gp.AddMessage("add value, psum = " + str(psum))
                urow.setValue("keep", "1")
                urows.UpdateRow(urow)
                
                slrow = slrows.Next()
                del zrow, zrows

        gp.AddMessage("del slrow")
        del slrow, slrows

        gp.AddMessage("\nThreshold = " + str(threshold) + "\nFinal sum = " + str(psum) + "\n")
        
        parameters.append("Output - final sum: " + str(psum))

    except:
        gp.AddError ("Error summing values: " + gp.GetMessages(2))
        raise Exception

    try:

        gp.AddMessage("\nCreating output files...")

##        # Join keep table to slice file
##        gp.MakeRasterLayer_management(sliced_grid, "slice_lyr", "#", "#")
##        gp.MakeTableView(slice_table_keep, "slice_keep_view")
##        gp.AddMessage("add join")
##        gp.AddJoin_management("slice_lyr", "VALUE", "slice_keep_view", "VALUE")
##        gp.CopyRaster_management("slice_lyr", slice_join)

        # Remap slice raster with keep values
        gp.AddMessage("lookup")
        gp.Lookup_sa(slice_grid_keep, "KEEP", slice_keep)
        gp.AddMessage("SOMA")
        gp.SingleOutputMapAlgebra_sa("CON( " + slice_keep + " == 1, " + score_grid + " )", score_selected)
        gp.AddMessage("\n\tCreated selected score grid output file: \n\t" + str(score_selected))

        # Remap original landuse raster with selected new activity cells
#        gp.SingleOutputMapAlgebra_sa("CON( " + slice_keep + " == 1, " + na_grid + ", -999 )", na_selected)
        gp.SingleOutputMapAlgebra_sa("CON( " + slice_keep + " == 1, " + na_grid + " )", na_selected_null)
        gp.SingleOutputMapAlgebra_sa("CON(IsNull(" + na_selected_null + "), 0, " + na_selected_null + " )", na_selected)
        gp.AddMessage("\n\tCreated selected new activity grid output file: \n\t" + str(na_selected))
        
#        gp.SingleOutputMapAlgebra_sa("CON( " + slice_keep + " == 1, " + na_selected + ", " + lu_grid + ")", final_landuse)
        gp.SingleOutputMapAlgebra_sa("CON( " + na_selected + " > 0, " + na_selected + ", " + lu_grid + ")", final_landuse)
        gp.AddMessage("\n\tCreated new land use output file: \n\t" + str(final_landuse))

        end_time = time.clock()
        # Write total time elapsed to reference output file
        tot_time_min = (end_time - start_time) / 60.0
        parameters.append("Elapsed time (minutes): " + str(tot_time_min))
        
    except:
        gp.AddError ("Error creating output files: " + gp.GetMessages(2))
        raise Exception




    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Post_Process\\PP_Sum_To_Threshold_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("POST-PROCESSING - SUM TO THRESHOLD\n")
        parafile.writelines("__________________________________\n\n")

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


    # REMOVE THIS
    gp.AddMessage("\nIGNORE THE FOLLOWING ERROR\n")
    raise Exception


except:
    gp.AddError ("Error running script")
    raise Exception
