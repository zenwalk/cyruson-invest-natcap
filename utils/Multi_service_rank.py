# ---------------------------------------------------------------------------
# Multi_service_rank.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
#
# Ranks the landscape based on level of service provision,
#   combining ranks from multiple services to identify overall
#   areas of highest multi-service provision
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

        # Services - list of service outputs from an InVEST model
        services = gp.GetParameterAsText(1)
        parameters.append("Services: " + services)
        service_list = services.split(";")

        # Weights - list of floating point values,
        #  weight to multiply all values in corresponding services by
        weights = gp.GetParameterAsText(2)
        parameters.append("Weights: " + weights)
        weight_list = weights.split(";")

        # Make sure that number of weights matches number of services
        if len(service_list) <> len(weight_list):
            gp.AddError("\nError: The number of services entered must match the number of weights.")
            raise Exception

        # Optional percent value for slicing result (to see top X% of service provision)
        percent = gp.GetParameterAsText(3)
        parameters.append("Grouping percent" + str(percent))
        if ((percent == "") or (percent == string.whitespace) or (percent == "#")):
            do_percent = False
        else:
            do_percent = True

        # Group by X% of ranking value distribution?
        group_by_dist = gp.GetParameterAsText(4)
        parameters.append("Group by value distribution: " + str(group_by_dist))

        # Group by X% of area?
        group_by_area = gp.GetParameterAsText(5)
        parameters.append("Group by area: " + group_by_area)

        if do_percent and group_by_dist == 'false' and group_by_area == 'false':
            gp.AddError("\nError: If grouping is to be done, Group by Value and/or Group by Area must be selected.")
            raise Exception

        if not do_percent and (group_by_dist == 'true' or group_by_area == 'true') :
            gp.AddError("\nError: If grouping is to be done, a Grouping Percent must be provided.")
            raise Exception

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(6)
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

        postprocws = gp.workspace + os.sep + "Post_process" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
    except:
        gp.AddError("Error creating output folders: " + gp.GetMessages())
        raise Exception


    # Script-created files
    try:

        # Intermediate files
        service_nw = interws + "serv_nw"
        service_sum_tmp = interws + "serv_sum_tmp"
        weighted_service_slice = interws + "serv_slice"
        weighted_service_neg = interws + "serv_neg"
        weighted_service_oid = interws + "serv_oid"

        # Output files
        weighted_service_sum = postprocws + "combined_service_rank" + Suffix + ".tif"

        
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
    
    # Process inputs
    try:

        # Process all services entered to make the combined ranking raster
        
        service_num = 0
        
        for service in service_list:
            # Normalize the service by its maximum value
            service_norm = interws + "serv" + str(service_num) + "_norm"
            service_max = gp.GetRasterProperties_management(service, "MAXIMUM")
            gp.Divide_sa(service, service_max, service_norm)

            service_weight = float(weight_list[service_num])
            
            # Create string of service/weight pairs required for Arc's Weighted Sum function             
            overlay_input = service_norm + " VALUE " + str(service_weight) + ";"

            service_num += 1

        # Remove trailing semicolon from string of service/weight pairs
        overlay_input = overlay_input[:-1]
        
        # Use Weighted Sum to combine normalized services and weights
        gp.WeightedSum_sa(overlay_input, weighted_service_sum)

        gp.AddMessage("\nCreated combined services ranking output: \n\t" + weighted_service_sum)

        # If a percent value has been entered, slice up result to show
        #   divisions of X% of service provision

        if do_percent:

            try:

                group_dist_zstat_table = interws + "group_dist_zstat.dbf"
                group_area_zstat_table = interws + "group_area_zstat.dbf"
                group_dist_poly_dissolve = interws + "group_dist_poly_dissolve.shp"
                weighted_service_group_dist_poly = interws + "combined_service_group_by_dist.shp"
                weighted_service_group_area_poly = interws + "combined_service_group_by_area.shp"

                weighted_service_group_dist_ras = postprocws + "combined_service_group_by_distribution" + Suffix + ".tif"
                weighted_service_group_dist_poly_final = postprocws + "combined_service_group_by_distribution" + Suffix + ".shp"
                weighted_service_group_area_ras = postprocws + "combined_service_group_by_area" + Suffix + ".tif"
                weighted_service_group_area_poly_final = postprocws + "combined_service_group_by_area" + Suffix + ".shp"

                gp.AddMessage("\nGrouping results by " + str(percent) + "%...")
                num_slices = int(round(100 / int(percent)))

                # Need output slices to assign 1 to largest values, but by default 1 goes to
                #   lowest values.  So multiply by -1 and do slice to get desired rankings
                gp.Times_sa(weighted_service_sum, "-1.0", weighted_service_neg)
                
                ### If group by value is chosen, do a slice by equal interval
                if group_by_dist == 'true':

                    
                    gp.AddMessage("\n\tGrouping by value distribution...")
                    gp.Slice_sa(weighted_service_neg, weighted_service_group_dist_ras, num_slices, "EQUAL_INTERVAL")                
                    gp.AddMessage("\n\tCreated grouped raster output by value distribution: \n\t" + weighted_service_group_dist_ras)
                    gp.RasterToPolygon_conversion(weighted_service_group_dist_ras, weighted_service_group_dist_poly, "NO_SIMPLIFY")

                    # RasterToPolygon creates a lot of polygons with the same group value.
                    # Merge each group value into a single polygon
                    gp.Dissolve_management(weighted_service_group_dist_poly, weighted_service_group_dist_poly_final, "GRIDCODE")
                    
                    # Add min/max values per group to the polygon attribute table
                    gp.AddField_management(weighted_service_group_dist_poly_final, "group_id", "long")
                    gp.AddField_management(weighted_service_group_dist_poly_final, "group_min", "float")
                    gp.AddField_management(weighted_service_group_dist_poly_final, "group_max", "float")
                    # Set group id to GRIDCODE, which is created from the VALUE field in RasterToPolygon
                    gp.CalculateField_management(weighted_service_group_dist_poly_final, "group_id", "[GRIDCODE]", "VB")
                    gp.DeleteField_management(weighted_service_group_dist_poly_final, "GRIDCODE")
                    
                    # Get min/max values per group                
                    gp.ZonalStatisticsAsTable_sa(weighted_service_group_dist_ras, "VALUE", weighted_service_sum, group_dist_zstat_table, "DATA")

                    dist_poly_rows = gp.UpdateCursor(weighted_service_group_dist_poly_final, "", "", "", "group_id A")
                    dist_poly_row = dist_poly_rows.Reset
                    dist_poly_row = dist_poly_rows.Next()

                    # Loop through polygons
                    while dist_poly_row:

                        # and find matching zonal stats min/max values for each
                        dist_zstat_rows = gp.SearchCursor(group_dist_zstat_table, "", "", "", "VALUE A")
                        dist_zstat_row = dist_zstat_rows.Reset
                        dist_zstat_row = dist_zstat_rows.Next()

                        while (int(dist_zstat_row.getValue("VALUE")) <> int(dist_poly_row.getValue("group_id"))):
                            dist_zstat_row = dist_zstat_rows.Next()

                        dist_poly_row.setValue("group_min", dist_zstat_row.getValue("MIN"))
                        dist_poly_row.setValue("group_max", dist_zstat_row.getValue("MAX"))
                        dist_poly_rows.UpdateRow(dist_poly_row)
                        dist_poly_row = dist_poly_rows.Next()
                    
                    gp.AddMessage("\n\tCreated grouped polygon output by value distribution: \n\t" + weighted_service_group_dist_poly_final)                    

                if group_by_area == 'true':
                    
                    gp.AddMessage("\n\tGrouping by area...")
            
                    gp.Slice_sa(weighted_service_neg, weighted_service_group_area_ras, num_slices, "EQUAL_AREA")
                    gp.AddMessage("\n\tCreated grouped raster output by area: \n\t" + weighted_service_group_area_ras)
                    gp.RasterToPolygon_conversion(weighted_service_group_area_ras, weighted_service_group_area_poly, "NO_SIMPLIFY")

                    # RasterToPolygon creates a lot of polygons with the same group value.
                    # Merge each group value into a single polygon
                    gp.Dissolve_management(weighted_service_group_area_poly, weighted_service_group_area_poly_final, "GRIDCODE")
                    
                    # Add min/max values per group to the polygon attribute table
                    gp.AddField_management(weighted_service_group_area_poly_final, "group_id", "long")
                    gp.AddField_management(weighted_service_group_area_poly_final, "group_min", "float")
                    gp.AddField_management(weighted_service_group_area_poly_final, "group_max", "float")
                    # Set group id to GRIDCODE, which is created from the VALUE field in RasterToPolygon
                    gp.CalculateField_management(weighted_service_group_area_poly_final, "group_id", "[GRIDCODE]", "VB")
                    gp.DeleteField_management(weighted_service_group_area_poly_final, "GRIDCODE")
                    
                    # Get min/max values per group                
                    gp.ZonalStatisticsAsTable_sa(weighted_service_group_area_ras, "VALUE", weighted_service_sum, group_area_zstat_table, "DATA")

                    area_poly_rows = gp.UpdateCursor(weighted_service_group_area_poly_final, "", "", "", "group_id A")
                    area_poly_row = area_poly_rows.Reset
                    area_poly_row = area_poly_rows.Next()

                    # Loop through polygons
                    while area_poly_row:

                        # and find matching zonal stats min/max values for each
                        area_zstat_rows = gp.SearchCursor(group_area_zstat_table, "", "", "", "VALUE A")
                        area_zstat_row = area_zstat_rows.Reset
                        area_zstat_row = area_zstat_rows.Next()

                        while (int(area_zstat_row.getValue("VALUE")) <> int(area_poly_row.getValue("group_id"))):
                            area_zstat_row = area_zstat_rows.Next()
                        
                        area_poly_row.setValue("group_min", area_zstat_row.getValue("MIN"))
                        area_poly_row.setValue("group_max", area_zstat_row.getValue("MAX"))
                        area_poly_rows.UpdateRow(area_poly_row)

                        area_poly_row = area_poly_rows.Next()
                    
                    gp.AddMessage("\n\tCreated grouped polygon output by area: \n\t" + weighted_service_group_area_poly_final)

            except:
                gp.AddError ("Error grouping results: " + gp.GetMessages(2))
                raise Exception
            
    except:
        gp.AddError ("Error calculating change: " + gp.GetMessages(2))
        raise Exception


    # Write input parameters to an output file for user reference
    try:
        gp.AddMessage("\nCreating parameter log file...")
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(postprocws + "Multi_Service_Rank_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("MULTI SERVICE RANK\n")
        parafile.writelines("__________________\n\n")

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
