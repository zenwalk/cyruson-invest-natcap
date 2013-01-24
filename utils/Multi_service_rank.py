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

        # Service 1 - service output from an InVEST model
        service1 = gp.GetParameterAsText(1)
        parameters.append("Service 1: " + service1)

        # Weight 1 - floating point value, weight to multiply all values in Service 1 by
        weight1 = gp.GetParameterAsText(2)
        parameters.append("Weight 1: " + weight1)

        # At least one service/weight is required, so no need to check for existence
        # Normalize service
        service1_norm = interws + "serv1_norm"
        service_max = gp.GetRasterProperties_management(service1, "MAXIMUM")
        gp.Divide_sa(service1, service_max, service1_norm)
        # Create string of service/weight pairs required for Arc's Weighted Sum function             
        overlay_input = service1_norm + " VALUE " + str(weight1)

        # Service 2 - service output from an InVEST model
        service2 = gp.GetParameterAsText(3)
        parameters.append("Service 2: " + service2)

        # Weight 2 - floating point value, weight to multiply all values in Service 2 by
        weight2 = gp.GetParameterAsText(4)
        parameters.append("Weight 2: " + weight2)

        # If service 2 is entered, a corresponding weight must be entered too
        if ((service2 != "") and (service2 != string.whitespace) and (service2 != "#")) \
           and ((weight2 == "") or (weight2 == string.whitespace) or (weight2 == "#")):
            gp.AddError("\nError: If service 2 is to be processed, a corresponding weight must be provided")
            raise Exception
        # If both service and weight are entered, process them
        elif ((service2 != "") and (service2 != string.whitespace) and (service2 != "#")) \
           and ((weight2 != "") and (weight2 != string.whitespace) and (weight2 != "#")):
            # Normalize service
            service2_norm = interws + "serv2_norm"
            service_max = gp.GetRasterProperties_management(service2, "MAXIMUM")
            gp.Divide_sa(service2, service_max, service2_norm)
            # Create string of service/weight pairs required for Arc's Weighted Sum function             
            overlay_input = overlay_input + "; " + service2_norm + " VALUE " + str(weight2)

        # Service 3 - service output from an InVEST model
        service3 = gp.GetParameterAsText(5)
        parameters.append("Service 3: " + service3)

        # Weight 3 - floating point value, weight to multiply all values in Service 3 by
        weight3 = gp.GetParameterAsText(6)
        parameters.append("Weight 3: " + weight3)

        # If service 3 is entered, a corresponding weight must be entered too
        if ((service3 != "") and (service3 != string.whitespace) and (service3 != "#")) \
           and ((weight3 == "") or (weight3 == string.whitespace) or (weight3 == "#")):
            gp.AddError("\nError: If service 3 is to be processed, a corresponding weight must be provided")
            raise Exception
        # If both service and weight are entered, process them
        elif ((service3 != "") and (service3 != string.whitespace) and (service3 != "#")) \
           and ((weight3 != "") and (weight3 != string.whitespace) and (weight3 != "#")):
            # Normalize service
            service3_norm = interws + "serv3_norm"
            service_max = gp.GetRasterProperties_management(service3, "MAXIMUM")
            gp.Divide_sa(service3, service_max, service3_norm)
            # Create string of service/weight pairs required for Arc's Weighted Sum function             
            overlay_input = overlay_input + "; " + service3_norm + " VALUE " + str(weight3)

        # Service 4 - service output from an InVEST model
        service4 = gp.GetParameterAsText(7)
        parameters.append("Service 4: " + service4)

        # Weight 4 - floating point value, weight to multiply all values in Service 4 by
        weight4 = gp.GetParameterAsText(8)
        parameters.append("Weight 4: " + weight4)

        # If service 4 is entered, a corresponding weight must be entered too
        if ((service4 != "") and (service4 != string.whitespace) and (service4 != "#")) \
           and ((weight4 == "") or (weight4 == string.whitespace) or (weight4 == "#")):
            gp.AddError("\nError: If service 4 is to be processed, a corresponding weight must be provided")
            raise Exception
        # If both service and weight are entered, process them
        elif ((service4 != "") and (service4 != string.whitespace) and (service4 != "#")) \
           and ((weight4 != "") and (weight4 != string.whitespace) and (weight4 != "#")):
            # Normalize service
            service4_norm = interws + "serv4_norm"
            service_max = gp.GetRasterProperties_management(service4, "MAXIMUM")
            gp.Divide_sa(service4, service_max, service4_norm)
            # Create string of service/weight pairs required for Arc's Weighted Sum function             
            overlay_input = overlay_input + "; " + service4_norm + " VALUE " + str(weight4)

        # Service 5 - service output from an InVEST model
        service5 = gp.GetParameterAsText(9)
        parameters.append("Service 5: " + service5)

        # Weight 5 - floating point value, weight to multiply all values in Service 5 by
        weight5 = gp.GetParameterAsText(10)
        parameters.append("Weight 5: " + weight5)

        # If service 5 is entered, a corresponding weight must be entered too
        if ((service5 != "") and (service5 != string.whitespace) and (service5 != "#")) \
           and ((weight5 == "") or (weight5 == string.whitespace) or (weight5 == "#")):
            gp.AddError("\nError: If service 5 is to be processed, a corresponding weight must be provided")
            raise Exception
        # If both service and weight are entered, process them
        elif ((service5 != "") and (service5 != string.whitespace) and (service5 != "#")) \
           and ((weight5 != "") and (weight5 != string.whitespace) and (weight5 != "#")):
            # Normalize service
            service5_norm = interws + "serv5_norm"
            service_max = gp.GetRasterProperties_management(service5, "MAXIMUM")
            gp.Divide_sa(service5, service_max, service5_norm)
            # Create string of service/weight pairs required for Arc's Weighted Sum function             
            overlay_input = overlay_input + "; " + service5_norm + " VALUE " + str(weight5)

        # Service 6 - service output from an InVEST model
        service6 = gp.GetParameterAsText(11)
        parameters.append("Service 6: " + service6)

        # Weight 6 - floating point value, weight to multiply all values in Service 6 by
        weight6 = gp.GetParameterAsText(12)
        parameters.append("Weight 6: " + weight6)

        # If service 6 is entered, a corresponding weight must be entered too
        if ((service6 != "") and (service6 != string.whitespace) and (service6 != "#")) \
           and ((weight6 == "") or (weight6 == string.whitespace) or (weight6 == "#")):
            gp.AddError("\nError: If service 6 is to be processed, a corresponding weight must be provided")
            raise Exception
        # If both service and weight are entered, process them
        elif ((service6 != "") and (service6 != string.whitespace) and (service6 != "#")) \
           and ((weight6 != "") and (weight6 != string.whitespace) and (weight6 != "#")):
            # Normalize service
            service6_norm = interws + "serv6_norm"
            service_max = gp.GetRasterProperties_management(service6, "MAXIMUM")
            gp.Divide_sa(service6, service_max, service6_norm)
            # Create string of service/weight pairs required for Arc's Weighted Sum function             
            overlay_input = overlay_input + "; " + service6_norm + " VALUE " + str(weight6)

        # Service 7 - service output from an InVEST model
        service7 = gp.GetParameterAsText(13)
        parameters.append("Service 7: " + service7)

        # Weight 7 - floating point value, weight to multiply all values in Service 7 by
        weight7 = gp.GetParameterAsText(14)
        parameters.append("Weight 7: " + weight7)

        # If service 7 is entered, a corresponding weight must be entered too
        if ((service7 != "") and (service7 != string.whitespace) and (service7 != "#")) \
           and ((weight7 == "") or (weight7 == string.whitespace) or (weight7 == "#")):
            gp.AddError("\nError: If service 7 is to be processed, a corresponding weight must be provided")
            raise Exception
        # If both service and weight are entered, process them
        elif ((service7 != "") and (service7 != string.whitespace) and (service7 != "#")) \
           and ((weight7 != "") and (weight7 != string.whitespace) and (weight7 != "#")):
            # Normalize service
            service7_norm = interws + "serv7_norm"
            service_max = gp.GetRasterProperties_management(service7, "MAXIMUM")
            gp.Divide_sa(service7, service_max, service7_norm)
            # Create string of service/weight pairs required for Arc's Weighted Sum function             
            overlay_input = overlay_input + "; " + service7_norm + " VALUE " + str(weight7)
            
        # Optional percent value for slicing result (to see top X% of service provision)
        percent = gp.GetParameterAsText(15)
        parameters.append("Grouping percent" + str(percent))
        if ((percent == "") or (percent == string.whitespace) or (percent == "#")):
            do_percent = False
        else:
            do_percent = True

        # Group by X% of ranking value?
        group_by_value = gp.GetParameterAsText(16)
        parameters.append("Group by value: " + str(group_by_value))

        # Or group by X% of area?
        group_by_area = gp.GetParameterAsText(17)
        parameters.append("Group by area: " + group_by_area)

        if do_percent and group_by_value == 'false' and group_by_area == 'false':
            gp.AddError("\nError: If grouping is to be done, Group by Value and/or Group by Area must be selected.")
            raise Exception

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

    # Script-created files
    try:

        # Intermediate files
        service_norm = interws + "serv_norm"
        service_nw = interws + "serv_nw"
        service_sum_tmp = interws + "serv_sum_tmp"
        weighted_service_slice = interws + "serv_slice"
        weighted_service_neg = interws + "serv_neg"
        weighted_service_oid = interws + "serv_oid"

        # Output files
        weighted_service_sum = postprocws + "weighted_service_sum" + Suffix + ".tif"
        weighted_service_group_value_ras = postprocws + "weighted_service_group_by_value" + Suffix + ".tif"
        weighted_service_group_value_poly = postprocws + "weighted_service_group_by_value" + Suffix + ".shp"
        weighted_service_group_area_ras = postprocws + "weighted_service_group_by_area" + Suffix + ".tif"
        weighted_service_group_area_poly = postprocws + "weighted_service_group_by_area" + Suffix + ".shp"
        
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
        gp.WeightedSum_sa(overlay_input, weighted_service_sum)

        gp.AddMessage("\nCreated combined services ranking output: \n\t" + weighted_service_sum)

        # If a percent value has been entered, slice up result to show
        #   divisions of X% of service provision

        if do_percent:

            gp.AddMessage("\nGrouping results by " + str(percent) + "%...")
            num_slices = int(round(100 / int(percent)))

            # Need output slices to assign 1 to largest values, but by default 1 goes to
            #   lowest values.  So multiply by -1 and do slice to get desired rankings
            gp.Times_sa(weighted_service_sum, "-1.0", weighted_service_neg)
            
            ### If group by value is chosen, do a slice by equal interval
            if group_by_value == 'true':
                
                gp.AddMessage("\n\tGrouping by value...")
                gp.Slice_sa(weighted_service_neg, weighted_service_group_value_ras, num_slices, "EQUAL_INTERVAL")
                gp.AddMessage("\n\tCreated grouped raster output by value: \n\t" + weighted_service_group_value_ras)
                gp.RasterToPolygon_conversion(weighted_service_group_value_ras, weighted_service_group_value_poly, "NO_SIMPLIFY")
                gp.AddMessage("\n\tCreated grouped polygon output by value: \n\t" + weighted_service_group_value_poly)

            if group_by_area == 'true':
                
                gp.AddMessage("\n\tGrouping by area...")
                gp.Slice_sa(weighted_service_neg, weighted_service_group_area_ras, num_slices, "EQUAL_AREA")
                gp.AddMessage("\n\tCreated grouped raster output by area: \n\t" + weighted_service_group_area_ras)
                gp.RasterToPolygon_conversion(weighted_service_group_area_ras, weighted_service_group_area_poly, "NO_SIMPLIFY")
                gp.AddMessage("\n\tCreated grouped polygon output by area: \n\t" + weighted_service_group_area_poly)
            
##            # But there might be missing slice values (if no pixel values fall in that interval)
##            #   so re-assign result to the OID+1, which will give 1->whatever
##            # IS THIS NECESSARY, OR SHOULD PEOPLE KNOW THAT THERE ARE MISSING INTERVALS?
##            gp.Lookup_sa(weighted_service_slice, "Rowid", weighted_service_oid)
##            # OID is 0->x, so add 1 to get a ranking of 1->x
##            gp.Plus_sa(weighted_service_oid, "1", weighted_service_final)

            
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
