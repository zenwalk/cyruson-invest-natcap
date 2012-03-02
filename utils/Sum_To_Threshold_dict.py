# ----------------------------------------------------------------------------------------
# Sum_To_Threshold.py
#
# By Stacie Wolny
# for the Natural Capital Project
# 
# Last edit: 10/10/2011 
#
# Chooses which grid cells on a landscape to select for water fund activities
# based on a score layer defining which cells are most important to target.
# Scores are ranked cell by cell, with a rank of 1 assigned to the highest
# (largest) score.  The ranked cells are looped through one by one, starting 
# with rank 1, adding up the activity cost associated with that cell, until the
# budget threshold is reached.
#
# Inputs:
#
#   Workspace: Folder where temporary and final output files are written
#   Score grid: Raster containing the scores/rankings assigned to each
#               cell on the landscape for a particular activity.
#               will be ordered largest to smallest, with the largest values
#               selected first, in descending order, for budget allocation.
#   Cost grid: Raster containing the cost per grid cell of the activity.
#               Cells without a cost should be assigned NoData.  The cost for each
#               Score-ordered cell will be added up until the budget Threshold is reached.
#   New activity grid: Raster containing activity land use code(s) for each
#               cell where the specified activity can occur.  Non-activity cells should
#               have the Value NoData.  Once the budget has been allocated based on Score, 
#               a raster will be output with the corresponding selected New Activity cells.
#   Original land use grid: Raster containing land use codes for each cell
#               on the current landscape.  The New Activity cells that are selected
#               in the budget process will be overlaid on this land use grid
#               to produce a scenario land use map that can be input into InVEST.
#   Watershed mask: Shapefile with a polygon defining the watershed where the budget
#               is to be allocated
#   Threshold: Number defining the amount of money that is to be allocated for
#               the specified activity in the Watershed entered above.
#
# Outputs:
#
#   score_sel: Scores corresponding to the cells selected for budget allocation
#   rank: Rank of the selected cells (from 1 -> N, where 1 has highest score)
#   na_sel: New activity cells selected for budget allocation
#   new_lu: Original land use grid, with selected new activity cells overlaid
#   
# ----------------------------------------------------------------------------------------

# Import system modules
import sys, string, os, arcgisscripting, time, datetime, operator
from operator import itemgetter

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

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(7)
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
        time_append = now.strftime("%Y%m%d%H%M")
        intermediate_name = "Intermediate" + time_append
        thefolders=["Post_Process", intermediate_name]
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
        interws = gp.workspace + os.sep + "Intermediate" + time_append + os.sep
        
        # Intermediate variables
        
        score_grid_clip = interws + "score_clip"
        cost_grid_clip = interws + "cost_clip"
        score_points = interws + "score_points" + Suffix + time_append + ".shp"
        cost_points = interws + "cost_points" + Suffix + time_append + ".shp"
        score_cost_keep = interws + "score_cost_keep" + Suffix + time_append + ".shp"
        score_cost_keep_sel = interws + "score_cost_keep_sel.shp"
        score_cost_keep_sel_update = interws + "score_cost_keep_sel_upd.shp"
        score_cost_keep_sfl = interws + "score_cost_keep_sfl.shp"
        score_cost_join = interws + "score_cost_join" + Suffix + time_append + ".shp"
        score_cost_join_fl = interws + "score_cost_join_fl.shp"
        score_cost_join_fl2 = interws + "score_cost_join_fl2.shp"
        score_cost_select = interws + "score_cost_sel.shp"
        summed_points = interws + "summed_points.shp"
        score_points_sort = interws + "score_points_sort.shp"
        na_selected_null = interws + "na_sel_null"

        
        # Output layers
        score_select_out = postprocws + "score_sel"
        rank_out = postprocws + "rank"
        na_selected = postprocws + "na_sel"
        final_landuse = postprocws + "new_lu"
        
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [score_select_out, rank_out, na_selected, final_landuse]           
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

        score_select_out = str(Outputnames[0])
        rank_out = str(Outputnames[1])
        na_selected = str(Outputnames[2])
        final_landuse = str(Outputnames[3])
            
    except:
        gp.AddError ("Error validating output filenames: " + gp.GetMessages(2))
        raise Exception


    try:

        gp.AddMessage("\nPreparing input grids...")

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
        gp.Mask = wshed_mask

        # Time this process for user reference
        start_time = time.clock()

        # Clip cost/score by mask
        gp.ExtractByMask_sa(cost_grid, wshed_mask, cost_grid_clip)
        gp.ExtractByMask_sa(score_grid, wshed_mask, score_grid_clip)

        # Can't easily work on a raster cell by cell, particularly
        # when multiple attribute fields are necessary, so turn into points
        gp.AddMessage("\nExtracting grids to points")
        # Make updates to a copy of the point grid file
        gp.RasterToPoint_conversion(cost_grid_clip, cost_points)
        gp.RasterToPoint_conversion(score_grid_clip, score_points)
        
        # Spatial Join merges fields with duplicate field names, so rename them
        gp.AddField_management(score_points, "sc_pointid", "LONG")
        gp.AddField_management(score_points, "sc_score", "DOUBLE")
        gp.CalculateField_management(score_points, "sc_pointid", "[POINTID]")
        gp.CalculateField_management(score_points, "sc_score", "[GRID_CODE]")
        gp.DeleteField_management(score_points, "POINTID;GRID_CODE")

        gp.AddField_management(cost_points, "co_pointid", "LONG")
        gp.AddField_management(cost_points, "co_cost", "DOUBLE")
        gp.CalculateField_management(cost_points, "co_pointid", "[POINTID]")
        gp.CalculateField_management(cost_points, "co_cost", "[GRID_CODE]")
        gp.DeleteField_management(cost_points, "POINTID;GRID_CODE")

 
        # Spatial join to merge cost grid with score grid
        gp.AddMessage("spatial join")
        gp.SpatialJoin_analysis(score_points, cost_points, score_cost_join, "JOIN_ONE_TO_ONE", "KEEP_COMMON", "", "INTERSECTS")
        
        # For those times when the cost/score layers don't quite line up        
##        gp.SpatialJoin_analysis(score_points, cost_points, score_cost_join, "JOIN_ONE_TO_ONE", "KEEP_COMMON", "", "CLOSEST")

    except:
        gp.AddError ("Error preparing input grids: " + gp.GetMessages(2))
        raise Exception


    try:
        # Sort points in descending order
        # Work around bug in 9.2 where UpdateCursor won't take an order argument
        
        scrows = gp.SearchCursor(score_cost_join, "", "", "", "sc_score D")
        scrow = scrows.Next()
        
        gp.AddMessage("Create ranking table...")
        
        table_name = "score_sort_tmp" + Suffix + time_append + ".dbf"
        score_sort_table = interws + table_name
        gp.CreateTable_management(interws, table_name, score_points)
        score_fields = ["sc_pointid", "sc_score"]
        gp.AddMessage("create ordered sort table")

        irows = gp.InsertCursor(score_sort_table)
        while scrow:
            irow = irows.newRow()
            for field in score_fields:
                irow.setValue(field, scrow.getValue(field))
            irows.insertRow(irow)
            scrow = scrows.Next()
        del scrow, scrows, irow, irows

        # Add rank field to sort table
        gp.AddField_management(score_sort_table, "rank_t", "LONG")

        gp.AddMessage("Assign score ranking...")
        
        # Add ranking value to each row in table (1 = highest (largest) score)
        usrows = gp.UpdateCursor(score_sort_table)

        usrow = usrows.Next()
        rank_val = 1

        while usrow:
            usrow.rank_t = rank_val
            usrows.UpdateRow(usrow)
            usrow = usrows.Next()
            rank_val += 1

        del usrow, usrows

        # Join rank table back to score/cost shapefile and transfer ranking over
        gp.AddMessage("copy to score_cost_keep")
        gp.MakeFeatureLayer_management(score_cost_join, score_cost_join_fl)
        gp.AddJoin_management(score_cost_join_fl, "sc_pointid", score_sort_table, "sc_pointid")
        gp.CopyFeatures_management(score_cost_join_fl, score_cost_keep)
        gp.AddMessage("add rank/cost/score fields")
        gp.AddField_management(score_cost_keep, "rank", "LONG")
        # CopyFeatures changes field names, so copy cost/point id/score to more useful field names too
        gp.AddField_management(score_cost_keep, "cost", "DOUBLE")
        gp.AddField_management(score_cost_keep, "score", "DOUBLE")

        # 9.3 changes CopyFields names yet again, somewhat for the better
        gp.AddMessage("calculate fields")
        gp.CalculateField_management(score_cost_keep, "rank", "[score_so_1]", "VB")
        gp.CalculateField_management(score_cost_keep, "cost", "[co_cost]", "VB")
        gp.CalculateField_management(score_cost_keep, "score", "[sc_score]", "VB")

        # Add a field denoting whether or not the point is selected as part of
        # the final budget allocation - set initially to zero for not selected
        gp.AddMessage("set keep to 0")
        gp.AddField_management(score_cost_keep, "keep", "SHORT")
        gp.CalculateField_management(score_cost_keep, "keep", "0")

    except:
        gp.AddError ("Error calculating ranking: " + gp.GetMessages(2))
        raise Exception

    
    try:

        # To save time, select out only the points with cost values
        # (those with a cost of 0 are not of interest)
        gp.AddMessage("select cost > 0")
        where_clause = "\"cost\" > 0 "
        gp.Select_analysis(score_cost_keep, score_cost_keep_sel, where_clause)

        # Go through score/cost/rank layer, starting with highest rank (1)
        # summing up cost until threshold is met 
        
        scrows = gp.SearchCursor(score_cost_keep_sel, "", "", "", "rank A")
        scrow = scrows.Next()

        # Working in a dictionary is much, much faster than cursors
        gp.AddMessage("populate dictionary")
        count = 0
        sc_dict = {}
        while scrow:
            sc_id = scrow.getValue("FID")
            sc_dict[sc_id] = [scrow.getValue("rank"), scrow.getValue("cost"), scrow.getValue("keep")]
            scrow = scrows.Next()
            count = count + 1

        # Sort by rank
        gp.AddMessage("sort values")
        values = sc_dict.values()
        values.sort()                             

        del scrow, scrows

        # Loop through point shapefile, summing one cost value at a time
        # until the budget threshold is reached

        psum = 0

        gp.AddMessage("\nSelecting cells for budget allocation...")


        for v in values:

            # Skip over any of the largest point values that are larger than the threshold
            if psum == 0 and v[1] > threshold:
                gp.AddMessage("Skipping cost value " + str(v[1]) + " - it's larger than the budget threshold")
            
            # Hit the threshold exactly
            elif (psum + v[1] == threshold):
                psum += v[1]
                
                # Set keep value to 1
                v[2] = 1

                gp.AddMessage("Done - hit the budget threshold exactly")
                break
            
            # Went over the threshold, so don't add current value
            elif (psum + v[1] > threshold):
                gp.AddMessage("Done - went over budget threshold")
                
                break

            # Still under the threshold, add current value and loop again
            else:
                psum += v[1]

                # Set keep value to 1
                v[2] = 1

        gp.AddMessage("\nBudget threshold = " + str(threshold) + "\nFinal allocated sum = " + str(psum) + "\n")
        parameters.append("Output - final allocated sum: " + str(psum))
            

    except:
        gp.AddError ("Error selecting cells for budget allocation: " + gp.GetMessages(2))
        raise Exception


    try:
        
        # Assign new keep values from dictionary to point shapefile

        gp.AddMessage("Create output layers...")

        # Sort values descending by 'keep', try to only process those that change to 1
        values.sort(key = itemgetter(2), reverse = True)

        lookupTblName = "rank_cost_keep_table"
        gp.CreateTable_management("in_memory", lookupTblName)
        gp.AddField_management("in_memory\\" + lookupTblName, "rank2", "LONG")
        gp.AddField_management("in_memory\\" + lookupTblName, "keep2", "LONG")
        insertRows = gp.insertcursor("in_memory\\" + lookupTblName)

        gp.AddMessage("for v in values")
        for v in values:
            insertRow = insertRows.newrow()
            insertRow.rank2 = v[0]
            insertRow.keep2 = v[2]
            insertRows.insertrow(insertRow)

        del insertRow
        del insertRows

        gp.AddMessage("copying to new table")
        newTable = interws + "rank_cost_keep_table2.dbf"
        gp.CopyRows_management("in_memory\\" + lookupTblName, newTable, "")
        gp.Delete_management("in_memory\\" + lookupTblName, "")

        gp.AddMessage("making new feature layer")
        gp.MakeFeatureLayer_management(score_cost_keep_sel, "sck_layer")
        gp.AddJoin_management("sck_layer", "rank", newTable, "rank2")
        gp.AddMessage("calculate field")
        gp.CopyFeatures_management("sck_layer", score_cost_keep_sel_update)
        gp.CalculateField_management(score_cost_keep_sel_update, "keep", "[rank_cos_1]")
        
        del sc_dict

        gp.AddMessage("Done with dictionary")
        
        # Save the selected points to a new file
        gp.AddMessage("select grid points to keep")
        select_exp = "\"keep\" = 1"
        gp.Select_analysis(score_cost_keep_sel_update, score_cost_select, select_exp)

        dsc = gp.Describe(score_grid)
        cell_size = str(dsc.MeanCellHeight)
        gp.Extent = dsc.Extent
        # Point to Raster doesn't seem to keep the point file's projection info
        # so set it explicitly
        gp.OutputCoordinateSystem = dsc.SpatialReference
        
        gp.AddMessage("point to raster outputs")
        gp.PointToRaster_conversion(score_cost_select, "score", score_select_out, "", "", cell_size)
        gp.AddMessage("\n\tCreated selected score output file: \n\t" + str(score_select_out))
        gp.PointToRaster_conversion(score_cost_keep, "rank", rank_out, "", "", cell_size)
        gp.AddMessage("\n\tCreated ranking output file: \n\t" + str(rank_out))

        # Remap original landuse raster with selected new activity cells
        gp.SingleOutputMapAlgebra_sa("CON( " + score_select_out + " > 0, " + na_grid + " )", na_selected_null)
        gp.SingleOutputMapAlgebra_sa("CON(IsNull(" + na_selected_null + "), 0, " + na_selected_null + " )", na_selected)
        gp.AddMessage("\n\tCreated selected new activity grid output file: \n\t" + str(na_selected))
        
        gp.SingleOutputMapAlgebra_sa("CON( " + na_selected + " > 0, " + na_selected + ", " + lu_grid + ")", final_landuse)
        gp.AddMessage("\n\tCreated new land use output file: \n\t" + str(final_landuse))

        end_time = time.clock()
        # Write total time elapsed to reference output file
        tot_time_min = (end_time - start_time) / 60.0
        parameters.append("Elapsed time (minutes): " + str(tot_time_min))

    except:
        gp.AddError ("Error creating output grids: " + gp.GetMessages(2))
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
##    gp.AddMessage("\nCleaning up temporary files...\n")
##    try:
##        gp.Delete_management(interws)
##    except:
##        gp.AddError("Error cleaning up temporary files:  " + gp.GetMessages(2))
##        raise Exception


except:
    gp.AddError ("\nError running script")
    raise Exception
