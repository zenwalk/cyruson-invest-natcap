# ---------------------------------------------------------------------------
# Sum_raster_by_overlapping_sshed.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
# 
# Last edit: 1/13/2012 
#
# Used for GEO BON tradeoffs
# Sums a raster (InVEST service output) by each polygon
# in a shapefile, adding them all up.  This gives the total service
# that each pixel provides, considering that it might be providing
# that service to several different overlapping servicesheds.
#
# NOTE: serviceshed shapefile must have an "sshed_id" field with unique ids
#
# Outputs a grid of the summed service 
#
# ---------------------------------------------------------------------------

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

        # Service raster
        service = gp.GetParameterAsText(1)
        parameters.append("Service: " + service)

        # Servicesheds to sum over
        ssheds = gp.GetParameterAsText(2)
        parameters.append("Servicesheds: " + ssheds)

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(3)
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
        
        service_clip = interws + "serv_clip"
        service_clip0 = interws + "serv_clip0"
        service_sum_tmp = interws + "ssum_tmp"
        sshed_sel = interws + "sshed_sel.shp"
        
        # Output layers
        service_sum = postprocws + "serv_sum"
        
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [service_sum]
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

        service_sum = str(Outputnames[0])

            
    except:
        gp.AddError ("Error validating output filenames: " + gp.GetMessages(2))
        raise Exception


    try:

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
        # Set working extent to service raster
        gp.Extent = service

        # Time this process for user reference
        start_time = time.clock()

        dsc = gp.describe(service)
        # Use cell size from service raster
        cell_size = str(dsc.MeanCellHeight)
        # Set extent of all intermediate outputs to full extent of input service raster
        extent = dsc.Extent

        # Create starter output sum file, containing all zeroes
        gp.CreateConstantRaster_sa(service_sum_tmp, 0, "FLOAT", cell_size, extent)

        # Loop through each serviceshed
        ss_rows = gp.SearchCursor(ssheds)
        ss_row = ss_rows.Reset
        ss_row = ss_rows.Next()

        while(ss_row):
            
            ss_id = str(int(ss_row.GetValue("sshed_id")))
            select_exp = "\"sshed_id\" = " + ss_id
            
            # Select one watershed and use as a mask for subsequent calculations
            gp.Select_analysis(ssheds, sshed_sel, select_exp)

            gp.AddMessage("\nProcessing serviceshed id " + str(ss_id) + "...")
           
            # Clip service raster by selected serviceshed
            gp.ExtractByMask_sa(service, sshed_sel, service_clip)

            # Set NoData to zero so rasters can be summed properly
            gp.AddMessage("SOMA")
            gp.SingleOutputMapAlgebra_sa("CON(isnull(" + service_clip + "), 0, " + service_clip + ")", service_clip0)

            # Add clipped service raster to running total
            gp.AddMessage("Plus")
            gp.Plus_sa(service_sum_tmp, service_clip0, service_sum)

            # Set tmp sum to total sum for processing the next round
            gp.AddMessage("CopyRaster")
            gp.CopyRaster_management(service_sum, service_sum_tmp)

            ss_row = ss_rows.Next()

        gp.AddMessage("\n\tCreated summed service output file: \n\t" + str(service_sum))

    except:
        gp.AddError ("Error summing by serviceshed: " + gp.GetMessages(2))
        raise Exception

    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Post_Process\\PP_Sum_raster_by_overlapping_sshed_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("POST-PROCESSING - SUM RASTER BY OVERLAPPING SERVICESHED\n")
        parafile.writelines("_______________________________________________________\n\n")

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


    # REMOVE THIS
##    gp.AddMessage("\nIGNORE THE FOLLOWING ERROR\n")
##    raise Exception


except:
    gp.AddError ("Error running script")
    raise Exception
