# ---------------------------------------------------------------------------
# WP_2_Nutrient_Retention.py
# 
# Coded by :
# Stacie Wolny, Bradley Eichelberger
# Written by:
# Driss Ennaanay, Guillermo Mendoza, Marc Conte 
# for the Natural Capital Project
#
# Last edited: 4/25/2011
#
# Second script for Water Purification.  Calculates biophysical outputs corresponding
# to how pollutants flow with water through a watershed to a point of interest.
#
# ---------------------------------------------------------------------------

# Import system modules
import sys, string, os, arcgisscripting, time, datetime

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out any necessary licenses
gp.CheckOutExtension("spatial")

# Set output handling
gp.OverwriteOutput = 1

# Set to true to enable reporting
verbose = True
cleanup = True

try:
    # Script arguments
    gp.AddMessage ("\nValidating arguments..." )
    try:
        # Make parameters array, and later write input parameter values to an output file
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))

        # Directory where output folder and files will be created
        gp.workspace = sys.argv[1]
        parameters.append("Workspace: " + gp.workspace)

        # DEM raster
        DEM = sys.argv[2]
        parameters.append("DEM: " + DEM)

        # Water yield raster from step 1
        water_yield = sys.argv[3]
        parameters.append("Water yield: " + water_yield)

        # Land use raster
        Landuse = sys.argv[4]
        parameters.append("Landcover: " + Landuse)

        # Watersheds shapefile
        watershed = sys.argv[5]
        parameters.append("Watersheds: " + watershed)

        # Sub-watersheds shapefile
        sub_watersheds = sys.argv[6]
        parameters.append("Sub-watersheds: " + sub_watersheds)

        # Table containing model coefficients per landuse class
        Biophys_Table = sys.argv[7]
        parameters.append("Biophysical coefficient table: " + Biophys_Table)

        # Table containing water purification threshold values per point of interest
        threshold_table = sys.argv[8]
        parameters.append("Water purification threshold table: " + threshold_table)

        # Are we evaluating Nitrogen loading?
        nitrogen = sys.argv[9]
        if nitrogen =='true':
            nitrogen = True
            parameters.append("Nitrogen: Yes")
        else:
            nitrogen = False
            parameters.append("Nitrogen: No")

        # or Phosphorus loading?
        phosphorus = sys.argv[10]
        if phosphorus =='true':
            phosphorus = True
            parameters.append("Phosphorus: Yes")
        else:
            phosphorus = False
            parameters.append("Phosphorus: No")

        # Threshold flow accumulation integer - how many upstream cells must flow into a call
        # before it's considered part of a stream
        Threshold_Flow_Accumulation = sys.argv[11]
        parameters.append("Threshold flow accumulation: " + str(Threshold_Flow_Accumulation))

        # Output resolution
        resolution = sys.argv[12]
        parameters.append("Output resolution: " + str(resolution))

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = sys.argv[13]
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
        thefolders=["Output","Intermediate","Service"]
        for folder in thefolders:
            if not gp.Exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
        pixelfolder = "Pixel"
        pfparent = gp.workspace + os.sep + "Output"
        if not gp.exists(pfparent + folder):
            gp.CreateFolder_management(pfparent, pixelfolder)
    except:
        gp.AddError("\nError creating output folders: " + gp.GetMessages())
        raise Exception


    # Output files
    
    try:
        # Base output/working directories
        outputws = gp.workspace + os.sep + "Output" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        servicews = gp.workspace + os.sep + "Service" + os.sep
        pixelws = gp.workspace + os.sep + "Output" + os.sep + "Pixel" + os.sep

        # Intermediate variables
        cum_yield = interws + "cum_yield"
        frac_removed = interws + "frac_removed"
        mean_runoff_index = interws + "mean_rindex"
        loads_x1000 = interws + "loads_x1000"
        loads = interws + "loads"
        preALV = interws + "preALV"
        preALV2 = interws + "preALV2"
        dir2stream = interws + "dir2stream"
        runoff_index = interws + "runoff_idx"
        HSS = interws + "HSS"
        filt_eff_dec = interws + "filt_eff_d"
        frac_remaining = interws + "frac_remain"
        filtration_efficiency = interws + "filt_eff"
        watersheds_join = interws + "ws_join"
        watersheds_calib = interws + "ws_calib"
        flowdir_rs = interws + "flowdir_rs"
        flowdir_ext = interws + "flowdir_ext"
        frac_removed_rs = interws + "frac_rem_rs"
        frac_removed_ext = interws + "frac_rem_ext"
        loads_ext = interws + "loads_ext"
        daccum_grid_raster_tmp = interws + "dac_ras_tmp"
        daccum_grid = interws + "dac_grid"
        daccum_grid_ascii = interws + "dac_grid_ascii.asc"
        export_grid_raster_tmp = interws + "exp_ras_tmp"
        export_grid = interws + "exp_grid"
        export_grid_ascii = interws + "exp_grid_ascii.asc"
        flowdir_ascii = interws + "flowdir_ascii.asc"
        loads_ascii = interws + "loads_ascii.asc"
        frac_removed_ascii = interws + "frac_removed_ascii.asc"
        total_cum_yield = interws + "tot_cum_yld"
        ws_mask_poly = interws + "ws_mask_poly.shp"
        ws_mask_ras = interws + "ws_mask_ras"
        wshed_ras = interws + "wshed_ras"
        wshed_ras2 = interws + "wshed_ras2"
        wsheds_j = interws + "wsheds_j"
        wshed_nut_join = interws + "wshed_nutjoin"
        wshed_annual_load = interws + "wshed_annld"
        wshed_num_cells = interws + "ws_numcells"
        sws_total_retention_table = interws + "sws_tot_ret_nut.dbf"
        sws_export_table = interws + "sws_export_nut.dbf"
        sws_adj_load_table = interws + "sws_adj_load.dbf"
        ws_export_table = interws + "ws_exp.dbf"
        ws_retention_table = interws + "ws_ret.dbf"
        wsheds_totret = interws + "ws_totret"
        wsheds_export = interws + "ws_export"
        total_retention1 = interws + "tot_ret1"
        allowed_load_cell = interws + "allowed_load"
        watersheds_sjoin = interws + "wshed_sjoin.shp"
        
        # Input table field names
        # Biophysical Models table
        lucode_field = "lucode"
        if nitrogen and phosphorus:
            gp.AddError("\nError: Either Nitrogen or Phosphorus must be chosen in the parameter window, not both.  Exiting.\n")
            raise Exception
        elif nitrogen:
            loading_field = "load_n"
            veg_filt_field = "eff_n"
            wp_annload_field = "thresh_n"
        elif phosphorus:
            loading_field = "load_p"
            veg_filt_field = "eff_p"
            wp_annload_field = "thresh_p"
        else:
            gp.AddError("\nError: Either Nitrogen or Phosphorus must be chosen in the parameter window.  Exiting.\n")
            raise Exception
        
        # Threshold table
        wp_id_field = "ws_id"

        # ID field for watersheds/sub-watersheds inputs
        wshed_id_field = "ws_id"
        subwshed_id_field = "subws_id"

        # Output layers
        v_stream = pixelws + "v_stream"
        adjusted_load = pixelws + "adj_load"
        ws_out_table_name = "nutrient_watershed" + Suffix + ".dbf"
        sws_out_table_name = "nutrient_subwatershed" + Suffix + ".dbf"
        daccum_raster_fname = "n_retain1"
        export_raster_fname = "n_export"

        total_retention = pixelws + "n_retain"
        sws_adj_load_mean = outputws + "adjl_mn"
        sws_adj_load_sum = outputws + "adjl_sm"
        sws_export_mean = outputws + "nexp_mn"
        sws_export_sum = outputws + "nexp_sm"
        # Service outputs
        sws_total_retention_mean = servicews + "nret_mn"
        sws_total_retention_sum = servicews + "nret_sm"
        
    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [export_raster_fname, daccum_raster_fname, v_stream, sws_adj_load_mean, sws_adj_load_sum, sws_export_mean, sws_export_sum, sws_total_retention_mean, sws_total_retention_sum, total_retention, adjusted_load]
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

        export_raster_fname = str(Outputnames[0])
        daccum_raster_fname = str(Outputnames[1])
        v_stream = str(Outputnames[2])
        sws_adj_load_mean = str(Outputnames[3])
        sws_adj_load_sum = str(Outputnames[4])
        sws_export_mean = str(Outputnames[5])
        sws_export_sum = str(Outputnames[6])
        sws_total_retention_mean = str(Outputnames[7])
        sws_total_retention_sum = str(Outputnames[8])
        total_retention = str(Outputnames[9])
        adjusted_load = str(Outputnames[10])


    except:
        gp.AddError ("\nError validating output filenames: " + gp.GetMessages(2))
        raise Exception


    # Check input raster projections - they should all be the same
    try:
        gp.AddMessage("\nChecking input raster projections...")
        DEMDesc = gp.describe(DEM)
        DEMspatref = DEMDesc.SpatialReference

        rasters = [Landuse, water_yield]
        for x in rasters:   
            rasterDesc=gp.describe(x)
            spatreflc = rasterDesc.SpatialReference
            if spatreflc.Type <> 'Projected':
                gp.AddMessage(x + " does not appear to be projected.  It is assumed to be in meters")
            elif spatreflc.LinearUnitName <> 'Meter':
                gp.AddMessage("This model assumes that data in " + x + " is projected in meters.  You may get erroneous results")
            if str(DEMspatref.name) <> str(spatreflc.name):
                gp.AddError("\nError: " + x + " is not in the same coordinate system as the hydrology rasters.  " + x + " is projected in " + spatreflc.name + " while the hydrology layers are in " + DEMspatref.name + ".  Please project raster in the same projection as the hydrology layers.\n")
                raise Exception
    except:
        gp.AddError("\nError in validating projections: " + gp.GetMessages(2)) 
        raise Exception



    # Set the Geoprocessing environment 
    try:
        gp.cellSize = resolution
        gp.mask = watershed
        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
        install_info = gp.GetInstallInfo("desktop")

    except:
        gp.AddError( "\nError setting geoprocessing environment: " + gp.GetMessages(2))
        raise Exception
        


    # Preprocess DEM derivatives and check for hydrologically correct rasters
    try:

        gp.AddMessage("\nProcessing hydrological rasters...")

        DEMdesc = gp.describe(DEM)
        DEMpath = DEMdesc.CatalogPath
        DEMpath2 = DEMpath.split("\\")
        z = DEMpath2[-1]
        DEMws = DEMpath.rstrip(z)
        if not gp.exists(DEMws + "Hydro_layers"):
            gp.CreateFolder_management(DEMws, "Hydro_layers")
            gp.AddMessage("\nCreating hydrology layers...")
        Hydrows = DEMws + "Hydro_layers" + os.sep
        if not gp.exists(Hydrows + "Flow_dir"):
            # Create flow direction raster
            gp.FlowDirection_sa(DEM, Hydrows + "Flow_dir", "NORMAL", Hydrows + "Drop_ras")

        # Make sure that flow direction values are only the cardinal values
        # The tool will not work if there are non-cardinal values
        
        flowdir_cardinals = [1, 2, 4, 8, 16, 32, 64, 128]

        # Arc10 doesn't create an attribute table for flow dir by default
##        if (install_info["Version"] == "10.0"):
        gp.BuildRasterAttributeTable_management(Hydrows + "Flow_dir", "OVERWRITE")

        fd_rows = gp.SearchCursor(Hydrows + "Flow_dir")
        fd_row = fd_rows.Reset
        fd_row = fd_rows.Next()

        while(fd_row):
            if flowdir_cardinals.count(fd_row.getValue("VALUE")) == 0:
                gp.AddError("\nError: A non-cardinal flow direction value (" + str(fd_row.getValue("VALUE")) + ") has been encountered in " + Hydrows + "Flow_dir" + "\n")
                raise Exception
            fd_row = fd_rows.Next()
        
        if not gp.exists(Hydrows + "Slope"):
            # Create slope raster
            gp.Slope_sa(DEM, Hydrows + "Slope", "PERCENT_RISE", "1")
        if not gp.exists(Hydrows + "Flow_acc"):
            # Create flow accumulation raster
            gp.FlowAccumulation_sa(Hydrows + "Flow_dir", Hydrows + "Flow_acc", "", "FLOAT")

        if not gp.exists(Hydrows + "Sinks"):
            # Check the number of sinks in the DEM
            gp.Sink_sa(Hydrows + "Flow_dir", Hydrows + "Sinks")

        # Problems with GetRasterProperties in Arc 9.3.1 and 10, 
        # so do a workaround for determining number of sinks
        gp.BuildRasterAttributeTable_management(Hydrows + "Sinks")
        si_rows = gp.SearchCursor(Hydrows + "Sinks")
        si_row = si_rows.Reset
        si_row = si_rows.Next()
        sinkcount = 0
        while si_row:
            sinkcount = sinkcount + 1
            si_row = si_rows.Next()
        del si_row, si_rows

        if sinkcount == 0:
            gp.AddWarning("DEM has no sinks or areas of internal drainage.")
        elif sinkcount <> 0:
            gp.AddWarning("DEM has " + str(sinkcount) + " sinks.  The DEM may not be hydrologically correct.")

        # Hydrology layers
        slope = Hydrows + "Slope"
        flow_acc = Hydrows + "Flow_acc"
        flow_dir = Hydrows + "Flow_dir"

    except:
        gp.AddError("\nError processing hydrology layers" + gp.GetMessages(2))
        raise Exception



    # Verify required input table fields' existence and type

    def checkfields (fields, table):
        
        try:

            table_fields = gp.listfields(table, "*", "All")
            table_field = table_fields.next()

            foundfields = []
            while (table_field):
                foundfields.append(table_field.Name.upper())
                # Make sure that fields that need to be integer are
                if (table_field.Name.upper() == lucode_field.upper()) or (table_field.Name.upper() == loading_field.upper()) or (table_field.Name.upper() == veg_filt_field.upper()) or (table_field.Name.upper() == wp_id_field.upper()) :
                    if not "Integer" in table_field.Type:
                        gp.AddError("\nError: Field " + str(table_field.Name) + " in table " + str(table) + " must be of type Integer\n")
                        raise Exception
                table_field = table_fields.next()

            # Make sure that all required fields ('fields') are in the input table ('foundfields')
            for f in fields:
                if not f.upper() in foundfields:
                    gp.AddError("\nError: Required table field " + str(f) + " not found in input table " + str(table) + "\n")
                    raise Exception

        except:
            gp.AddError ("\nError verifying input table fields")
            raise Exception


    # Used for manual workaround of Workspace to New Mosaic
    def pixel_type(name):
        try:
            if (name == "U1"):
                return "1_bit"
            if (name == "U2"):
                return "2_bit"
            if (name == "U4"):
                return "4_bit"
            if (name == "U8"):
                return "8_bit_unsigned"
            if (name == "S8"):
                return "8_bit_signed"
            if (name == "U16"):
                return "16_bit_unsigned"
            if (name == "S16"):
                return "16_bit_signed"
            if (name =="U32"):
                return "32_bit_unsigned"
            if (name =="S32"):
                return "32_bit_signed"
            if (name =="F32"):
                return "32_bit_float"
            if (name == "D64"):
                return "64_bit"

        except:
            AddError("\nError determining sub-watershed pixel type...")
            raise Exception

    
    # Calculating runoff index
    # Ranking of geomorphological factors - water yield includes soil depth and available water content
    gp.AddMessage ("\nCalculating runoff index...")
    try:

        # Verify input table fields
        checkfields([lucode_field, loading_field, veg_filt_field], Biophys_Table)
        checkfields([wp_id_field, wp_annload_field], threshold_table)
        
        # Cumulative water yield to each cell
        gp.FlowAccumulation_sa(flow_dir, cum_yield, water_yield, "FLOAT")
        
        # Flow accumulation doesn't include the cell's contribution, so add it
        gp.Plus_sa(cum_yield, water_yield, total_cum_yield)
        gp.SingleOutputMapAlgebra_sa("LOG10( " + total_cum_yield + " )", runoff_index)

    except:
        gp.AddError ("\nError calculating runoff index: " + gp.GetMessages(2))
        raise Exception

    
    # Calculating Adjusted Loading Values
    # High values where geomorphology is conducive to runoff and pollutant loading is high
    gp.AddMessage ("\nCalculating adjusted loading values...")
    try:
        # Create loading grid from input table
        gp.ReclassByTable_sa(Landuse, Biophys_Table,  lucode_field, lucode_field, loading_field, loads_x1000, "DATA")
        # Divide by 1000, since input table values are x1000
        gp.SingleOutputMapAlgebra_sa("FLOAT(" + loads_x1000 + ") / 1000.0", loads)
        # Mean runoff index value for each watershed
        gp.SingleOutputMapAlgebra_sa("ZONALMEAN( " + watershed + ", " + runoff_index + " )", mean_runoff_index)
        # Hydrologic Sensitivity Score (HSS) - geomorphology only
        gp.Divide_sa(runoff_index, mean_runoff_index, HSS)
        # Adjusted Export Values - HSS and pollutant loading
        gp.SingleOutputMapAlgebra_sa("Float(" + loads + ") * " + HSS, preALV)
        # Convert kg/hectare to m^2, as provided by the input table loads field
        dsc = gp.describe(preALV)
        Cell_Size = str(dsc.MeanCellHeight)
        gp.SingleOutputMapAlgebra_sa(preALV + " * " + Cell_Size + " * " + Cell_Size, preALV2)
        gp.Divide_sa(preALV2, "10000.0", adjusted_load)

        # Aggregate by sub-watersheds

        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, adjusted_load, sws_adj_load_sum, "SUM", "DATA")

        area_ha_field = "AREA_HA"
        mean_ha_field = "MEAN_HA"

        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, adjusted_load, sws_adj_load_table, "DATA")
        gp.AddField(sws_adj_load_table, area_ha_field, "double")
        gp.AddField(sws_adj_load_table, mean_ha_field, "double")
        # Compute values per hectare
        gp.CalculateField_management(sws_adj_load_table, area_ha_field, "[AREA] / 10000", "VB")
        gp.CalculateField_management(sws_adj_load_table, mean_ha_field, "[SUM] / [" + area_ha_field + "]", "VB")

        gp.FeatureToRaster_conversion(sub_watersheds, subwshed_id_field, wshed_ras, Cell_Size)
        # Arc10 - attribute table not created by default
        gp.BuildRasterAttributeTable_management(wshed_ras, "OVERWRITE")
        gp.MakeRasterLayer_management(wshed_ras, "wsheds", "#", "#")

        # Zonal stats field name has changed in Arc 10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = subwshed_id_field
        else:
            zstat_id_field = "VALUE"

        # Map per-hectare values to sub-watersheds for raster output
        gp.AddJoin_management("wsheds", "VALUE", sws_adj_load_table, zstat_id_field)
        gp.CopyRaster_management("wsheds", wsheds_j)
        gp.Lookup_sa(wsheds_j, mean_ha_field, sws_adj_load_mean)

        gp.Delete_management("wsheds")


        gp.AddMessage("\n\tCreated adjusted load outputs:\n\t" + str(sws_adj_load_sum) + "\n\t" + str(sws_adj_load_mean))
        
    except:
        gp.AddError ("\nError calculating adjusted loading values: " + gp.GetMessages(2))
        raise Exception


    # Calculating filtration efficiency
    # Just a mapping from the model coefficient table's efficiency field
    gp.AddMessage ("\nCalculating fraction of nutrient removed...")
    try:
        gp.ReclassByTable_sa(Landuse, Biophys_Table, lucode_field, lucode_field, veg_filt_field, filtration_efficiency, "DATA")
        # Divide by 100 to change input table integer percent values to decimal
        gp.SingleOutputMapAlgebra_sa("FLOAT(" + filtration_efficiency + ") / 100.0", filt_eff_dec)
        # Fraction of pollutant not retained by a cell
        gp.SingleOutputMapAlgebra_sa("MAX(0, 1 - " + filt_eff_dec + ")", frac_remaining)
        # Fraction of pollutant retained by cell
        gp.Minus_sa("1.0", frac_remaining, frac_removed)

    except:
        gp.AddError ("\nError calculating fraction of nutrient removed: " + gp.GetMessages(2))
        raise Exception

#
# Calculate the amount of nutrient that is removed by each cell
# from what passes through the cell from upstream (daccum)
#
# Also keep track of how much of each cell's nutrient load makes it
# to the stream (export)
#
# Loop through sub-watersheds and process each separately,
# recombining them afterward to allow for arbitrarily large
# watersheds
#

    gp.AddMessage("\nCalculating nutrient removal...")
    try:

        subws_fields = gp.ListFields(sub_watersheds, subwshed_id_field)
        subws_field = subws_fields.Next()
        subws_type = subws_field.Type

        # Constrain watershed id type to integer, to facilitate making output
        # raster filenames without problematic characters
        if not "Integer" in subws_type:
            gp.AddError("\nError: Sub-watershed id field must be of type integer \n")
            raise Exception
        
        ws_rows = gp.SearchCursor(sub_watersheds)
        ws_row = ws_rows.Reset
        ws_row = ws_rows.Next()

        # Keep count of number of watersheds processed
        subws_count = 0

        # Strings to hold lists of sub-watershed filenames
        # so they can be mosaicked back together later
        export_mosaic_filenames = ""
        daccum_mosaic_filenames = ""

        # Loop through each sub-watershed in input watersheds file
        while (ws_row):

            subws_count += 1

            # Get watershed id
            subws_id = str(int(ws_row.GetValue(subwshed_id_field)))
            select_exp = "\"subws_id\" = " + subws_id

            # Use cell size from load raster
            dsc = gp.describe(adjusted_load)
            cell_size = str(dsc.MeanCellHeight)
            
            # Select one watershed and output to ws_mask_file for use by gp.Mask
            gp.Select_analysis(sub_watersheds, ws_mask_poly, select_exp)
            gp.FeatureToRaster_conversion(ws_mask_poly, subwshed_id_field, ws_mask_ras, cell_size)
            gp.Mask = ws_mask_ras

            gp.AddMessage("\nProcessing sub-watershed id " + str(subws_id) + "...")


            # Create sub-watershed intermediate files
            try:
                # Make new directory to hold sub-watershed outputs
                # Used to mosaic together later
                export_sub_wsheds_dir = "export_sub_wsheds" + Suffix  
                daccum_sub_wsheds_dir = "daccum_sub_wsheds" + Suffix

                thefolders=[export_sub_wsheds_dir, daccum_sub_wsheds_dir]
                for folder in thefolders:
                    if not gp.Exists(folder):
                        gp.CreateFolder_management(interws, folder)
                
                export_raster_subws = interws + export_sub_wsheds_dir + os.sep + "w_exp"
                daccum_grid_raster_subws = interws + daccum_sub_wsheds_dir + os.sep + "w_ret"

                # Add suffix to end of output filenames
                # Constrain length of output raster filenames to 13 characters
                Suffix_ws = "_" + str(subws_id) + Suffix
                Outputnames = [export_raster_subws, daccum_grid_raster_subws]
                num = 0
                for x in Outputnames:
                    y = x.split("\\")
                    z = y[-1]
                    lenz = len(z + Suffix_ws)
                    if lenz > 13:
                        x2 = x.rstrip(z)
                        overlen = 13 - lenz 
                        newz = z[0:overlen]
                        x2 = x2 + newz
                        Outputnames[num] = x2 + Suffix_ws
                    else:
                        Outputnames[num] = x + Suffix_ws
                    num = num + 1       
                export_raster_subws =  str(Outputnames[0])
                daccum_grid_raster_subws = str(Outputnames[1])

            except:
                gp.AddError ("\nError generating sub-watershed output filenames: " + gp.GetMessages(2))
                raise Exception           

            try:

                gp.AddMessage("\n\tPre-processing input/output grids...")

                # Count flow length to and from stream by changing flow direction
                # such that direction = 0 at the stream threshold value
                gp.LessThanEqual_sa(flow_acc, Threshold_Flow_Accumulation, v_stream)
                gp.Times_sa(flow_dir, v_stream, dir2stream)

                # Make sure all input rasters are the same resolution
                # as the loads raster and resample if necessary
                try:
                    gp.AddMessage("\n\tValidating input raster resolution...")
                    rasterDescL=gp.describe(adjusted_load)
                    loads_cell_width = rasterDescL.MeanCellWidth
                    
                    rasterDescFD=gp.describe(dir2stream)
                    if rasterDescFD.MeanCellWidth <> loads_cell_width:
                        gp.AddMessage("\n\tResampling flow direction raster...")
                        gp.Resample_management(dir2stream, flowdir_rs, loads_cell_width)
                    else:
                        gp.CopyRaster_management(dir2stream, flowdir_rs)
                    rasterDescFR=gp.describe(frac_removed)
                    if rasterDescFR.MeanCellWidth <> loads_cell_width:
                        gp.AddMessage("\n\tResampling removal efficiency raster...")
                        gp.Resample_management(frac_removed, frac_removed_rs, loads_cell_width)
                    else:
                        gp.CopyRaster_management(frac_removed, frac_removed_rs)
                except:
                    gp.AddError("\nError validating input raster resolution: " + gp.GetMessages(2)) 
                    raise Exception

                # Extract input rasters, using the input mask
                # so that all extents are the same

                # Pollutant loading at each cell
                gp.ExtractByMask_sa(adjusted_load, ws_mask_ras, loads_ext)
                # Flow direction 
                gp.ExtractByMask_sa(flowdir_rs, ws_mask_ras, flowdir_ext)
                gp.Delete_management(flowdir_rs)
                # Fraction of incoming pollutant removed by each cell
                gp.ExtractByMask_sa(frac_removed_rs, ws_mask_ras, frac_removed_ext)
                gp.Delete_management(frac_removed_rs)

                # Create accumulation and export grids from loads input raster, fill with zeroes
                gp.Times_sa(loads_ext, "0.0", daccum_grid)
                gp.Times_sa(loads_ext, "0.0", export_grid)

                # Convert flowdir, loads, frac_removed and daccum_grid to ASCII and read into arrays
                gp.RasterToASCII_conversion(flowdir_ext, flowdir_ascii)
                gp.RasterToASCII_conversion(daccum_grid, daccum_grid_ascii)
                gp.RasterToASCII_conversion(export_grid, export_grid_ascii)
                gp.RasterToASCII_conversion(loads_ext, loads_ascii)
                gp.RasterToASCII_conversion(frac_removed_ext, frac_removed_ascii)

                # Open ASCII files for reading/writing
                input_flowdir = open(flowdir_ascii, 'r')
                input_loads = open(loads_ascii, 'r')
                input_frac_removed = open(frac_removed_ascii, 'r')
                
                # New for Arc 10: read/write is now 'r+'
                if (install_info["Version"] == "10.0"):
                    daccum_recgrid = open(daccum_grid_ascii, 'r+')
                    export_recgrid = open(export_grid_ascii, 'r+')
                else:
                    daccum_recgrid = open(daccum_grid_ascii, 'rw')
                    export_recgrid = open(export_grid_ascii, 'rw')

                # Initialize input/output arrays
                flowdir_array = []
                loads_array = []
                frac_removed_array = []
                acc_load_removed_array = []
                export_array = []

                # Read ASCII rasters into arrays
                x = 0
                for line in input_flowdir.readlines():
                    x = x + 1
                    # First 6 lines are headers
                    if x == 1:
                        numcols_fd = int(line.split(" ")[-1])
                    if x == 2:
                        numrows_fd = int(line.split(" ")[-1])
                    if x == 3:
                        xcorner_fd = float(line.split(" ")[-1])
                    if x == 4:
                        ycorner_fd = float(line.split(" ")[-1])
                    if x == 5:
                        cellsize_fd = float(line.split(" ")[-1])
                    if x == 6:
                        nodataValue_fd = int(line.split(" ")[-1])
                    elif x > 6:
                        flowdir_array.append([])
                        for i in line.split():
                            flowdir_array[-1].append(int(i))

                x = 0
                for line in daccum_recgrid.readlines():
                    x = x + 1
                    # First 6 lines are headers
                    if x == 1:
                        numcols_out = int(line.split(" ")[-1])
                    if x == 2:
                        numrows_out = int(line.split(" ")[-1])
                    if x == 3:
                        xcorner_out = float(line.split(" ")[-1])
                    if x == 4:
                        ycorner_out = float(line.split(" ")[-1])
                    if x == 5:
                        cellsize_out = float(line.split(" ")[-1])
                    if x == 6:
                        nodataValue_out = int(line.split(" ")[-1])
                    elif x > 6:
                        acc_load_removed_array.append([])
                        for i in line.split():
                            acc_load_removed_array[-1].append(float(i))

                x = 0
                for line in export_recgrid.readlines():
                    x = x + 1
                    # First 6 lines are headers
                    if x == 1:
                        numcols_exp = int(line.split(" ")[-1])
                    if x == 2:
                        numrows_exp = int(line.split(" ")[-1])
                    if x == 3:
                        xcorner_exp = float(line.split(" ")[-1])
                    if x == 4:
                        ycorner_exp = float(line.split(" ")[-1])
                    if x == 5:
                        cellsize_exp = float(line.split(" ")[-1])
                    if x == 6:
                        nodataValue_exp = int(line.split(" ")[-1])
                    elif x > 6:
                        export_array.append([])
                        for i in line.split():
                            export_array[-1].append(float(i))
                            
                x = 0
                for line in input_frac_removed.readlines():
                    x = x + 1
                    # First 6 lines are headers
                    if x == 1:
                        numcols_fr = int(line.split(" ")[-1])
                    if x == 2:
                        numrows_fr = int(line.split(" ")[-1])
                    if x == 3:
                        xcorner_fr = float(line.split(" ")[-1])
                    if x == 4:
                        ycorner_fr = float(line.split(" ")[-1])
                    if x == 5:
                        cellsize_fr = float(line.split(" ")[-1])
                    if x == 6:
                        nodataValue_fr = int(line.split(" ")[-1])
                    elif x > 6:
                        frac_removed_array.append([])
                        for i in line.split():
                            frac_removed_array[-1].append(float(i))

                x = 0
                for line in input_loads.readlines():
                    x = x + 1
                    # First 6 lines are headers
                    if x == 1:
                        numcols_loads = int(line.split(" ")[-1])
                    if x == 2:
                        numrows_loads = int(line.split(" ")[-1])
                    if x == 3:
                        xcorner_loads = float(line.split(" ")[-1])
                    if x == 4:
                        ycorner_loads = float(line.split(" ")[-1])
                    if x == 5:
                        cellsize_loads = float(line.split(" ")[-1])
                    if x == 6:
                        nodataValue_loads = int(line.split(" ")[-1])
                    elif x > 6:
                        listLine = []
                        loads_array.append([])
                        for i in line.split():
                            loads_array[-1].append(float(i))

                input_flowdir.close()
                daccum_recgrid.close()
                export_recgrid.close()
                input_loads.close()
                input_frac_removed.close()
                
            except:
                gp.AddError("\nError pre-processing sub-watershed " + str(subws_id) + ": " + gp.GetMessages(2))
                gp.AddError("It's possible that the sub-watershed is too large.")
                gp.AddError("Skipping to next sub-watershed...\n")
                input_flowdir.close()
                daccum_recgrid.close()
                export_recgrid.close()
                input_loads.close()
                input_frac_removed.close()
                ws_row = ws_rows.Next()
                continue

            # Follow flow path of each cell's load,
            #   calculating nutrient removal along the way
            # Flow tracing based on logic from ArcScript AS15793 - Tracegrid
            #   by Christopher C. Duncan, duncan@geo.umass.edu         

            try:

                # Flow direction loop debugging
                #
                # If it looks like there's a circular flow path problem, uncomment the following three lines,
                # as well as two others later in the code beginning with 'outfile', also with the comment
                # 'Flow direction loop debugging'
                #
                # This will simply create a file where data values for each cell are written as it is processed.
                # Keep an eye on the last lines of this file as the script runs, and look for the same set of 
                # values repeating - this indicates a circular flow path.

                # Uncomment these three lines
                # outfile_name = gp.workspace + "\\Output\\wp_loop_debug_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt"
                # outfile = open(outfile_name, "w")
                # gp.AddMessage("\nFlow direction loop debugging is on, output file is " + str(outfile_name) + "\n")

                # Time sub-watershed process.  If it's taking too long, break out, it's probably in an infinite loop.
                start_time = time.clock()

                # Loop over all cells in the loads grid
                gp.AddMessage("\n\tCalculating nutrient removal per cell...")
                
                for y in range (0, numrows_loads):
                    for x in range (0, numcols_loads):
                        
                        more = 1

                        # Only process cells with loading data
                        # Note: Python list positions are denoted by [row, column], so
                        # the following matrix positions are [y,x]

                        # Skip cell if load = NoData or 0 or flow direction or fraction
                        # removed is NoData

                        if (loads_array[y][x] != nodataValue_loads) and (loads_array[y][x] != 0) and \
                           (flowdir_array[y][x] != nodataValue_fd) and (frac_removed_array[y][x] != nodataValue_fr):
                            
                            # Traverse flow path to end
                            first_cell = 1
                            while more:
                                # First cell in the flow path
                                if (first_cell == 1):
                                    cur_x = x
                                    cur_y = y
                                else:
                                    # Not first cell
                                    cur_x = next_x
                                    cur_y = next_y

                                # Flow direction loop debugging - to debug, uncomment the following line 
                                #   and 2 others beginning with 'outfile'
                                # outfile.writelines("\nloads: " + str(loads_array[cur_y][cur_x]) + "; flowdir: " + str(flowdir_array[cur_y][cur_x]) + "; frac_removed: " + str(frac_removed_array[cur_y][cur_x]))


                                # Assign x,y values for next cell in flow path based on flowdir value
                                # Flow direction values must only be the cardinal directions - 1, 2, 4, 8, 16, 32, 64, 128

                                if (flowdir_array[cur_y][cur_x] == 1 or flowdir_array[cur_y][cur_x] == 2 or flowdir_array[cur_y][cur_x] == 128):
                                    next_x = cur_x + 1
                                elif (flowdir_array[cur_y][cur_x] == 8 or flowdir_array[cur_y][cur_x] == 16 or flowdir_array[cur_y][cur_x] == 32):
                                    next_x = cur_x - 1
                                else:
                                    next_x = cur_x

                                if (flowdir_array[cur_y][cur_x] == 2 or flowdir_array[cur_y][cur_x] == 4 or flowdir_array[cur_y][cur_x] == 8):
                                    next_y = cur_y + 1
                                elif (flowdir_array[cur_y][cur_x] == 32 or flowdir_array[cur_y][cur_x] == 64 or flowdir_array[cur_y][cur_x] == 128):
                                    next_y = cur_y - 1
                                else:
                                    next_y = cur_y

                                # Process loading for next cell

                                # For first cell, no filtration, just pass it all along
                                if (first_cell == 1):
                                    load_passed = loads_array[cur_y][cur_x]
                                    load_retained = 0
                                    
                                # For the rest, each cell will remove a fraction of what comes into it from the cell upstream
                                # But only process if cell has frac_removed values
                                elif (frac_removed_array[cur_y][cur_x] != nodataValue_fr):
                                    load_retained = frac_removed_array[cur_y][cur_x] * load_passed
                                    load_passed = (1 - frac_removed_array[cur_y][cur_x]) * load_passed

                                    # Keep track of how much upstream pollutant is retained by the cell
                                    if (load_passed <= 0):
                                        load_passed == 0
                                        more = 0
                                    acc_load_removed_array[cur_y][cur_x] += load_retained
                                else:
                                    more = 0
                               
                                # Done if next cell goes off the edge of the loads map
                                if (next_x >= numcols_loads or next_y >= numrows_loads \
                                      or next_x < 0 or next_y < 0):
                                    more = 0
                                # or next cell in flow path has flow direction or frac_removed value NoData
                                elif (flowdir_array[next_y][next_x] == nodataValue_fd) or (frac_removed_array[next_y][next_x] == nodataValue_fr):
                                    more = 0
                                # or next cell has flow direction value of 0 (has hit a stream)
                                elif (flowdir_array[next_y][next_x] == 0):
                                    more = 0
                                    # Keep track of how much of each cell's nutrient makes it to
                                    # a stream, stored in the export array
                                    export_array[y][x] = load_passed

                                # or this sub-watershed has been processing for too long (45 minutes by default)
                                # - probably indicates that there's a circular flow path
                                # To check for circular flow paths, look for code lines beginning with 'outfile'
                                cur_time = time.clock()
                                
                                # 45 minutes = 2700 seconds (this can be changed)
                                too_long = 2700
                                if ((cur_time - start_time) >= too_long):
                                    gp.AddError("\nError: Sub-watershed " + str(subws_id) + " is taking too long (45 minutes).  This probably indicates that there's a circular flow path.\n")
                                    gp.AddError("\tSkipping to next sub-watershed...\n")
                                    # Free up memory for next sub-watershed
                                    del flowdir_array, loads_array, frac_removed_array, acc_load_removed_array, export_array
                                    # On to next row in sub-watershed file
                                    ws_row = ws_rows.Next()
                                    # Don't raise exception, skip rest of loop code and go on to next sub-watershed
                                    continue
                                
                                first_cell = 0 

            except:
                gp.AddError("\nError calculating nutrient removal per cell: " + gp.GetMessages(2))
                raise Exception


            # Process export and decay accumulation output files
            try:
                gp.AddMessage("\n\tProcessing nutrient removal output rasters...")

                # Flow direction loop debugging - to debug, uncomment the following line 
                #   and 2 others beginning with 'outfile'
                # outfile.close()

                # Write output array to ASCII file
                
                if (os.path.exists(daccum_grid_ascii) == True):
                    os.remove(daccum_grid_ascii)

                if (os.path.exists(export_grid_ascii) == True):
                    os.remove(export_grid_ascii)

                # Add headers to ASCII output file
                print >> open(daccum_grid_ascii, 'a'), "ncols " + " "*8 + str(numcols_out)
                print >> open(daccum_grid_ascii, 'a'), "nrows " + " "*8 + str(numrows_out)
                print >> open(daccum_grid_ascii, 'a'), "xllcorner" + " "*5 + str(xcorner_out)
                print >> open(daccum_grid_ascii, 'a'), "yllcorner" + " "*5 + str(ycorner_out)
                print >> open(daccum_grid_ascii, 'a'), "cellsize" + " "*6 + str(cellsize_out)
                print >> open(daccum_grid_ascii, 'a'), "NODATA_value" + " "*2 + str(nodataValue_out)

                print >> open(export_grid_ascii, 'a'), "ncols " + " "*8 + str(numcols_exp)
                print >> open(export_grid_ascii, 'a'), "nrows " + " "*8 + str(numrows_exp)
                print >> open(export_grid_ascii, 'a'), "xllcorner" + " "*5 + str(xcorner_exp)
                print >> open(export_grid_ascii, 'a'), "yllcorner" + " "*5 + str(ycorner_exp)
                print >> open(export_grid_ascii, 'a'), "cellsize" + " "*6 + str(cellsize_exp)
                print >> open(export_grid_ascii, 'a'), "NODATA_value" + " "*2 + str(nodataValue_exp)

                # Add output array, one row at a time
                nr = 0
                while nr < numrows_out:
                    x = 0
                    text = ""
                    for row in acc_load_removed_array[nr]:
                        text = text + " " + str(acc_load_removed_array[nr][x])
                        x = x + 1
                    print >> open(daccum_grid_ascii, 'a'), text[1:]
                    nr = nr + 1

                nr = 0
                while nr < numrows_exp:
                    x = 0
                    text = ""
                    for row in export_array[nr]:
                        text = text + " " + str(export_array[nr][x])
                        x = x + 1
                    print >> open(export_grid_ascii, 'a'), text[1:]
                    nr = nr + 1

                # Define output grid's projection as being the same as the input loads raster
                gp.ASCIIToRaster_conversion(daccum_grid_ascii, daccum_grid_raster_tmp, "FLOAT")
                gp.DefineProjection_management(daccum_grid_raster_tmp, daccum_grid)
                gp.ASCIIToRaster_conversion(export_grid_ascii, export_grid_raster_tmp, "FLOAT")
                gp.DefineProjection_management(export_grid_raster_tmp, export_grid)

                # Work around weird edge condition where decay accumulation is assigned
                # the value -9999, but it is not translated to NoData like other cells
                # with that value.  Happens rarely, and only along the raster edge, where 
                # layers don't necessarily line up exactly, even after extraction.

                gp.SingleOutputMapAlgebra_sa("SETNULL(" + daccum_grid_raster_tmp + " < 0, " + daccum_grid_raster_tmp + ")", daccum_grid_raster_subws)
                gp.SingleOutputMapAlgebra_sa("SETNULL(" + export_grid_raster_tmp + " < 0, " + export_grid_raster_tmp + ")", export_raster_subws)

                # Add these outputs to a string list of sub-watersheds
                # to be mosaicked back together later
                export_mosaic_filenames = export_mosaic_filenames + export_raster_subws + ";"
                daccum_mosaic_filenames = daccum_mosaic_filenames + daccum_grid_raster_subws + ";"

                # Free up memory for next sub-watershed
                del flowdir_array, loads_array, frac_removed_array, acc_load_removed_array, export_array

                # On to next row in sub-watershed file
                ws_row = ws_rows.Next()

            except:
                gp.AddError("\nError processing nutrient removal output raster: " + gp.GetMessages(2))
                raise Exception

    except:
        gp.AddError("\nError processing watershed id " + str(subws_id) + ": " + gp.GetMessages(2))
        raise Exception

     # Combine sub-watershed output into whole watershed grids for remaining processing
        
    try:

        # Reset mask to whole watershed
        gp.Mask = watershed

        export_raster = pixelws + export_raster_fname
        daccum_raster = interws + daccum_raster_fname

        if gp.Exists(export_raster):
            gp.Delete_management(export_raster)
        if gp.Exists(daccum_raster):
            gp.Delete_management(daccum_raster)
        
        if (subws_count > 1):
            # More than one watershed input, mosaic them together
            gp.AddMessage("\nCombining watershed outputs...")
            
            export_dir = interws + export_sub_wsheds_dir
            daccum_dir = interws + daccum_sub_wsheds_dir

            # WorkspaceToNewMosaic throws an error in ArcMap 9.3.1
            # so do a manual workaround 
            # Code based on the original ArcMap WorkspaceToNewMosaic.py

            # Remove trailing semi-colon from mosaic file lists
            export_mosaic_filenames = export_mosaic_filenames[:-1]
            daccum_mosaic_filenames = daccum_mosaic_filenames[:-1]

            # Split file list into an array and use one file to get
            # values for creating the mosaic
            export_mosaic_filenames_array = export_mosaic_filenames.split(';')
            tmpfile = export_mosaic_filenames_array[0]
            tmpfile_desc = gp.Describe(tmpfile)
            bandcount = tmpfile_desc.BandCount
            pixeltype = pixel_type(tmpfile_desc.PixelType)
            spatialref = tmpfile_desc.SpatialReference

            gp.CreateRasterDataset_management(pixelws, export_raster_fname, "#", pixeltype, spatialref, bandcount)
            gp.CreateRasterDataset_management(interws, daccum_raster_fname, "#", pixeltype, spatialref, bandcount)

            gp.Mosaic_management(export_mosaic_filenames, export_raster, "MAXIMUM")
            gp.Mosaic_management(daccum_mosaic_filenames, daccum_raster, "MAXIMUM")            
            
        else:
            # Only one watershed input, no mosaic
            gp.CopyRaster_management(export_raster_subws, export_raster)
            gp.CopyRaster_management(daccum_grid_raster_subws, daccum_raster)


    except:
        gp.AddError("\nError combining sub-watershed outputs: " + gp.GetMessages(2))
        raise Exception

    # Subtract amount of nutrient allowed by the annual load table input

    try:
    
        gp.AddMessage("\nCalculating nutrient allowance...")

        # Get sediment values from table, map to watershed
        # Copy to dbf to work around table OID requirement change in ArcMap 9.3
        dsc = gp.describe(daccum_raster)
        cell_size = str(dsc.MeanCellHeight)

        gp.CopyRows_management(threshold_table, "nut_val_tmp.dbf")
        gp.FeatureToRaster_conversion(watershed, wshed_id_field, wshed_ras2, cell_size)
        gp.BuildRasterAttributeTable_management(wshed_ras2, "OVERWRITE")
        gp.MakeRasterLayer_management(wshed_ras2, "wshed_tmp_nut")
        gp.AddJoin_management("wshed_tmp_nut", "Value", "nut_val_tmp.dbf", wshed_id_field)
        gp.CopyRaster_management("wshed_tmp_nut", wshed_nut_join)
        # Load allowed per cell = table wq_annload / number of cells in watershed
        gp.Lookup_sa(wshed_nut_join, wp_annload_field, wshed_annual_load)
        gp.Lookup_sa(wshed_nut_join, "COUNT", wshed_num_cells)
        gp.SingleOutputMapAlgebra_sa("FLOAT(" + wshed_annual_load + ") / FLOAT(" + wshed_num_cells + ")", allowed_load_cell)
        # Subtract allowed from total retention to get the service
        gp.Minus_sa(daccum_raster, allowed_load_cell, total_retention1)
        # Change negative values to zeroes
        gp.SingleOutputMapAlgebra_sa("CON(" + total_retention1 + " < 0, 0, " + total_retention1 + ")", total_retention)

        gp.Delete_management("nut_val_tmp.dbf")
        gp.Delete_management("wshed_tmp_nut")

    except:
        gp.AddError("\nError calculating nutrient allowance: " + gp.GetMessages(2))
        raise Exception

    # Aggregate by sub-watersheds
    try:

        # Raster and table outputs for retention and export
        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, total_retention, sws_total_retention_sum, "SUM", "DATA")
        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, export_raster, sws_export_sum, "SUM", "DATA")
        
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, total_retention, sws_total_retention_table, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, export_raster, sws_export_table, "DATA")

        gp.AddField(sws_total_retention_table, area_ha_field, "double")
        gp.AddField(sws_total_retention_table, mean_ha_field, "double")
        # Compute values per hectare
        gp.CalculateField_management(sws_total_retention_table, area_ha_field, "[AREA] / 10000", "VB")
        gp.CalculateField_management(sws_total_retention_table, mean_ha_field, "[SUM] / [" + area_ha_field + "]", "VB")

        gp.AddField(sws_export_table, area_ha_field, "double")
        gp.AddField(sws_export_table, mean_ha_field, "double")
        # Compute values per hectare
        gp.CalculateField_management(sws_export_table, area_ha_field, "[AREA] / 10000", "VB")
        gp.CalculateField_management(sws_export_table, mean_ha_field, "[SUM] / [" + area_ha_field + "]", "VB")

        gp.MakeRasterLayer_management(wshed_ras, "wsheds_tr", "#", "#")
        gp.MakeTableView(sws_total_retention_table, "totret_wq_table_view")
        gp.AddJoin_management("wsheds_tr", "VALUE", "totret_wq_table_view", zstat_id_field)
        gp.CopyRaster_management("wsheds_tr", wsheds_totret)
        gp.Lookup_sa(wsheds_totret, mean_ha_field, sws_total_retention_mean)

        gp.MakeRasterLayer_management(wshed_ras, "wsheds_exp", "#", "#")
        gp.MakeTableView(sws_export_table, "exp_table_view")
        gp.AddJoin_management("wsheds_exp", "VALUE", "exp_table_view", zstat_id_field)
        gp.CopyRaster_management("wsheds_exp", wsheds_export)
        gp.Lookup_sa(wsheds_export, mean_ha_field, sws_export_mean)

        if (install_info["Version"] == "10.0"):
            gp.Delete_management("totret_wq_table_view")
            gp.Delete_management("wsheds_tr")
            gp.Delete_management("exp_table_view")
            gp.Delete_management("wsheds_exp")
        
        gp.AddMessage("\nCreated sub-watershed outputs:\n\t" + str(sws_total_retention_sum) + "\n\t" + str(sws_total_retention_mean))
        gp.AddMessage("\n\t" + str(sws_export_mean) + "\n\t" + str(sws_export_sum))

    except:
        gp.AddError ("\nError aggregating by sub-watersheds: " + gp.GetMessages(2))
        raise Exception


    # Table to hold generated output values per watershed
    gp.AddMessage("\nCreating watershed nutrient output table...")
    try:
        # output table field names
        out_table_wshed_id_field = "ws_id"
        out_table_export_field = "nut_export"
        out_table_retention_field = "nut_retain"

        gp.CreateTable_management(outputws, ws_out_table_name)
        ws_out_table = outputws + ws_out_table_name

        gp.AddField(ws_out_table, out_table_wshed_id_field, "long")
        gp.AddField(ws_out_table, out_table_export_field, "double")
        gp.AddField(ws_out_table, out_table_retention_field, "double")
        # Remove Field1 - it's added by default and not used
        gp.DeleteField_management(ws_out_table, "Field1")

        out_table_rows = gp.InsertCursor(ws_out_table)

        # Zonal stats field name changed in Arc10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = wshed_id_field
        else:
            zstat_id_field = "Value"

        # Populate table with watershed export/retention values
        gp.ZonalStatisticsAsTable_sa(watershed, wshed_id_field, export_raster, ws_export_table, "DATA")
        wse_rows = gp.SearchCursor(ws_export_table, "", "", "", zstat_id_field + " A")
        wse_row = wse_rows.Reset
        wse_row = wse_rows.Next()

        gp.ZonalStatisticsAsTable_sa(watershed, wshed_id_field, total_retention, ws_retention_table, "DATA")
        wsr_rows = gp.SearchCursor(ws_retention_table, "", "", "", zstat_id_field + " A")
        wsr_row = wsr_rows.Reset
        wsr_row = wsr_rows.Next()

        while (wse_row):
            # Zonal stats field name changed in Arc 10
            if (install_info["Version"] == "10.0"):
                ws_id = wse_row.getValue(wshed_id_field)
            else:
                ws_id = wse_row.getValue("VALUE")
                
            ws_export = float(wse_row.getValue("SUM"))
            ws_retention = float(wsr_row.getValue("SUM"))
            new_row = out_table_rows.NewRow()
            new_row.setValue(out_table_wshed_id_field, ws_id)
            new_row.setValue(out_table_export_field, ws_export)
            new_row.setValue(out_table_retention_field, ws_retention)

            out_table_rows.InsertRow(new_row)
            wse_row = wse_rows.Next()
            wsr_row = wsr_rows.Next()

        gp.AddMessage("\n\tCreated watershed nutrient output table: \n\t" + str(ws_out_table))
        
    except:
        gp.AddError ("\nError creating watershed nutrient output table:" + gp.GetMessages(2))
        raise Exception

    # Table to hold generated output values per sub-watershed
    gp.AddMessage("\nCreating sub-watershed nutrient output table...")
    try:
        # output table field names
        out_table_wsid_field = "ws_id"
        out_table_subwshed_id_field = "subws_id"
        out_table_export_field = "nut_export"
        out_table_retention_field = "nut_retain"
        
        gp.CreateTable_management(outputws, sws_out_table_name)
        sws_out_table = outputws + sws_out_table_name

        gp.AddField(sws_out_table, out_table_subwshed_id_field, "long")
        gp.AddField(sws_out_table, out_table_wsid_field, "long")
        gp.AddField(sws_out_table, out_table_export_field, "double")
        gp.AddField(sws_out_table, out_table_retention_field, "double")
        # Remove Field1 - it's added by default and not used
        gp.DeleteField_management(sws_out_table, "Field1")

        sws_out_table_rows = gp.InsertCursor(sws_out_table)

        # Map sub-watersheds to their corresponding watersheds
        gp.SpatialJoin_analysis(sub_watersheds, watershed, watersheds_sjoin, "JOIN_ONE_TO_ONE", "KEEP_ALL", "#", "IS_WITHIN")

        # Zonal stats field name changed in Arc10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = subwshed_id_field
        else:
            zstat_id_field = "Value"
            
        # Populate output table with watershed export/retention values
        swse_rows = gp.SearchCursor(sws_export_table, "", "", "", zstat_id_field + " A")
        swse_row = swse_rows.Reset
        swse_row = swse_rows.Next()

        swsr_rows = gp.SearchCursor(sws_total_retention_table, "", "", "", zstat_id_field + " A")
        swsr_row = swsr_rows.Reset
        swsr_row = swsr_rows.Next()
        
        while (swse_row):

            sj_rows = gp.SearchCursor(watersheds_sjoin, "", "", "", subwshed_id_field + " A")
            sj_row = sj_rows.Reset
            sj_row = sj_rows.Next()

            while (int(sj_row.getValue(subwshed_id_field)) <> int(swse_row.getValue(zstat_id_field))):
                sj_row = sj_rows.Next()
        
            # Zonal stats field name changed in Arc 10
            if (install_info["Version"] == "10.0"):
                sws_id = swse_row.getValue(subwshed_id_field)
            else:
                sws_id = swse_row.getValue("VALUE")
                
            sws_export = float(swse_row.getValue("SUM"))
            sws_retention = float(swsr_row.getValue("SUM"))
            
            new_row = sws_out_table_rows.NewRow()
            new_row.setValue(out_table_subwshed_id_field, sws_id)
            new_row.setValue(out_table_wsid_field, int(sj_row.getValue(wshed_id_field)))
            new_row.setValue(out_table_export_field, sws_export)
            new_row.setValue(out_table_retention_field, sws_retention)

            sws_out_table_rows.InsertRow(new_row)
            swse_row = swse_rows.Next()
            swsr_row = swsr_rows.Next()
            
            del sj_row, sj_rows

        gp.AddMessage("\n\tCreated sub-watershed nutrient output table: \n\t" + str(sws_out_table))

    except:
        gp.AddError ("\nError creating sub-watershed nutrient output table:" + gp.GetMessages(2))
        raise Exception



    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Output\\WP_Nutrient_Retention_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("WATER PURIFICATION 2 - NUTRIENT RETENTION MODEL PARAMETERS\n")
        parafile.writelines("__________________________________________________________\n\n")

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
        del ws_row, ws_rows, new_row, out_table_rows, wse_row, wse_rows, fd_row, fd_rows
        del sws_out_table_rows, swse_row, swse_rows, swsr_row, swsr_rows
        gp.Delete_management(interws)
    except:
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception


except:
    gp.AddError ("\nError running script")
    raise Exception
