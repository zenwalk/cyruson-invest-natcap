# ---------------------------------------------------------------------------
# Prepare_DEM.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
#
# Performs several functions to help prepare the DEM to go into InVEST:
# - Fills sinks
# - Burns streams
# - Fills small areas of missing data 
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

        # Digital Elevation Model raster
        DEM = gp.GetParameterAsText(1)
        parameters.append("DEM: " + DEM)

        # Fill sinks?
        fill_sinks = gp.GetParameterAsText(2)
        if fill_sinks =='true':
            fill_sinks = True
            parameters.append("Fill sinks: Yes")
        else:
            fill_sinks = False
            parameters.append("Fill sinks: No")            

        # Burn streams?
        burn_streams = gp.GetParameterAsText(3)
        if burn_streams =='true':
            burn_streams = True
            parameters.append("Burn streams: Yes")
        else:
            burn_streams = False
            parameters.append("Burn streams: No")

        # Stream shapefile
        streams = gp.GetParameterAsText(4)
        parameters.append("Streams: " + streams)
        if burn_streams and ((streams == "") or (streams == string.whitespace) or (streams == "#")):
            gp.AddError("\nError: If streams are to be burned, a streams shapefile must be provided.\n")
            raise Exception

        # Depth to burn the streams
        stream_depth = gp.GetParameterAsText(5)
        parameters.append("Stream burn depth: " + str(stream_depth))
        if burn_streams and ((stream_depth == "") or (stream_depth == string.whitespace) or (stream_depth == "#")):
            gp.AddError("\nError: If streams are to be burned, a depth value must be provided.\n")
            raise Exception

        # Fill small holes (missing data)?
        fill_small_holes = gp.GetParameterAsText(6)
        if fill_small_holes =='true':
            fill_small_holes = True
            parameters.append("Fill small holes: Yes")
        else:
            fill_small_holes = False
            parameters.append("Fill small holes: No")

        # Fill large holes (missing data)?
        fill_large_holes = gp.GetParameterAsText(7)
        if fill_large_holes =='true':
            fill_large_holes = True
            parameters.append("Fill large holes: Yes")
        else:
            fill_large_holes = False
            parameters.append("Fill large holes: No")

        if not fill_large_holes and not fill_small_holes and not burn_streams and not fill_sinks:
            gp.AddError("\nError: Nothing to do - Please check boxes next to the DEM preparation steps to be done.\n")
            raise Exeption

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(8)
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
        thefolders=["Output", "Intermediate"]
        for folder in thefolders:
            if not gp.Exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        gp.AddError("\nError creating folders: " + gp.GetMessages())
        raise Exception


    # Output files
    try:
        # Base output/working directories
        postprocws = gp.workspace + os.sep + "Output" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        
        # Intermediate variables
        dem_fill_sinks = interws + "dem_fill_sinks.tif"
        dem_fill_sinks_copy = interws + "dem_fill_sinks_copy.tif"
        streams_copy = interws + "streams_copy.shp"
        streams_ras = interws + "streams_ras.tif"
        streams_ras0 = interws+ "streams_ras0.tif"
        dem_fill_holes = interws + "dem_fill_holes.tif"
        dem_burned = interws + "dem_burned.tif"
        dem_burned_filled = interws + "dem_burned_filled.tif"
        dem_5x5 = interws + "dem_fill_holes_5x5.tif"
        dem_11x11 = interws + "dem_fill_holes_11x11.tif"
        dem_points = interws + "dem_points.shp"
        dem_interp = interws + "dem_interp.tif"
        dem_lg_holes = interws + "dem_lg_holes.tif"
        flow_dir = interws + "flow_dir.tif"
        
        # Output layers
        dem_final = postprocws + "dem_prep" + Suffix + ".tif"       
        
    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception

    # Set environment and general variables
    try:
        dsc = gp.describe(DEM)
        cell_size = str(dsc.MeanCellHeight)
        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws                    
    except:
        gp.AddError ("\nError setting environment: " + gp.GetMessages(2))
        raise Exception


    # Fill small holes (missing data) in the DEM
    # Do two passes, 5x5 and 11x11
    try:
        if fill_small_holes:
            gp.AddMessage("\nFilling small holes...")
            # First, fill smallest holes with a 5x5 focal area
            gp.SingleOutputMapAlgebra_sa("CON(IsNull(" + DEM + "), FOCALMEAN(" + DEM + ", RECTANGLE, 5, 5), " + DEM + ")", dem_5x5)
            # Next, fill larger holes with an 11x11 focal area
            gp.SingleOutputMapAlgebra_sa("CON(IsNull(" + dem_5x5 + "), FOCALMEAN(" + dem_5x5 + ", RECTANGLE, 11, 11), " + dem_5x5 + ")", dem_11x11)            
    except:
        gp.AddError ("\nError filling small holes: " + gp.GetMessages(2))
        raise Exception

    # Fill large holes (missing data) in the DEM
    # Interpolate over the DEM and use interpolated values to fill holes
    try:
        if fill_large_holes:
            gp.AddMessage("\nFilling large holes (may take a long time)...")
            # Use DEM with small holes filled, if done, otherwise use input DEM
            if fill_small_holes:
                use_dem = dem_11x11
            else:
                use_dem = DEM
                
            # Turn DEM into points for interpolation
            gp.RasterToPoint_conversion(use_dem, dem_points)
            # Make sure interpolated output is same extent and cell size as DEM and snap to its alignment
            gp.extent = dsc.extent
            gp.SnapRaster = DEM
            # Interpolate over points using IDW interpolation
            # Use a power value of 3 to weight closer cells more heavily and do less smoothing
            gp.Idw_sa(dem_points, "GRID_CODE", dem_interp, cell_size, "3")
            
            # Fill DEM holes with interpolated data
            gp.SingleOutputMapAlgebra_sa("CON(IsNull(" + use_dem + ")," + dem_interp + ", " + use_dem + ")", dem_lg_holes)            
    except:
        gp.AddError ("\nError filling large holes: " + gp.GetMessages(2))
        raise Exception


    # Burn streams into the DEM
    try:
        if burn_streams:
            gp.AddMessage ("\nBurning streams...")
            stream_field = "str"
            gp.CopyFeatures_management(streams, streams_copy)
            # Create a field with single value to be used in creating raster streams
            gp.AddField_management(streams_copy, stream_field, "SHORT")
            gp.AddMessage("calculate field")
            gp.CalculateField_management(streams_copy, stream_field, stream_depth, "PYTHON")
            # Make sure stream layer is same extent and cell size as DEM and snap to its alignment
            gp.extent = dsc.extent
            gp.SnapRaster = DEM
            gp.FeatureToRaster_conversion(streams_copy, stream_field, streams_ras, cell_size)
            # Set NoDatas to zero
            gp.SingleOutputMapAlgebra_sa("CON(IsNull(" + streams_ras + "), 0, " + streams_ras + ")", streams_ras0)
            # Subtract stream raster from DEM
            if fill_large_holes:
                use_dem = dem_lg_holes
            elif fill_small_holes:
                use_dem = dem_11x11
            else:
                use_dem = DEM
            gp.Minus_sa(use_dem, streams_ras0, dem_burned)

    except:
        gp.AddError ("\nError burning streams: " + gp.GetMessages(2))
        raise Exception


    # Fill sinks
    # Multiple passes are sometimes necessary, so do up to 3
    #   stopping if/when flow directions are cardinal
    try:
        if fill_sinks:
            gp.AddMessage ("\nFilling sinks...")

            if burn_streams:
                use_dem = dem_burned
            elif fill_large_holes:
                use_dem = dem_lg_holes
            elif fill_small_holes:
                use_dem = dem_11x11
            else:
                use_dem = DEM

            fill_count = 1
            fill_max = 3

            while(fill_count <= fill_max):

                gp.AddMessage("\n\tPass " + str(fill_count) + " of " + str(fill_max))
                
                # If filling has already been done, use the filled raster
                if gp.exists(dem_fill_sinks):
                    gp.CopyRaster(dem_fill_sinks, dem_fill_sinks_copy)
                    use_dem2 = dem_fill_sinks_copy
                else:
                    use_dem2 = use_dem

                gp.Fill_sa(use_dem2, dem_fill_sinks)

                # Check result against flow direction
                # If there are non-cardinal values, do it again, up to 3 times
                gp.AddMessage("\n\tChecking for non-cardinal flow direction values...")
                gp.FlowDirection_sa(dem_fill_sinks, flow_dir, "NORMAL")
                flowdir_cardinals = [1, 2, 4, 8, 16, 32, 64, 128]

                gp.BuildRasterAttributeTable_management(flow_dir, "OVERWRITE")

                fd_rows = gp.SearchCursor(flow_dir)
                fd_row = fd_rows.Reset
                fd_row = fd_rows.Next()

                fill_again = False

                while(fd_row):
                    if flowdir_cardinals.count(fd_row.getValue("VALUE")) == 0:
                        gp.AddMessage("\n\tA non-cardinal flow direction value has been encountered, filling again... \n")
                        fill_again = True
                        
                    fd_row = fd_rows.Next()

                if fill_again:
                    fill_count += 1
                    if fill_count > fill_max:
                        gp.AddError("\n\tWarning: DEM has been filled 3 times and still has non-cardinal flow directions.")
                # Doesn't need to be filled again, break out of loop
                else:
                    break
        
    except:
        gp.AddError ("\nError filling sinks: " + gp.GetMessages(2))
        raise Exception

    # Final output
    try:
        gp.AddMessage("\nCreating final output...")

        if fill_sinks:
            use_dem = dem_fill_sinks
        elif burn_streams:
            use_dem = dem_burned
        elif fill_large_holes:
            use_dem = dem_lg_holes
        else:
            use_dem = dem_11x11

        gp.CopyRaster_management(use_dem, dem_final)
        gp.AddMessage("\nCreated prepared DEM output: " + dem_final + "\n")

    except:
        gp.AddError ("\nError creating final output: " + gp.GetMessages(2))
        raise Exception


    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Output\\Prepare_DEM_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("PREPARE DEM\n")
        parafile.writelines("___________\n\n")

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
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception

except:
    gp.AddError ("\nError running script")
    raise Exception
