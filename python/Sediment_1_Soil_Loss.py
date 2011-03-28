# ---------------------------------------------------------------------------
# Sediment_1_Soil_Loss.py
# 
# Coded by :
# Stacie Wolny, Bradley Eichelberger
# Written by:
# Driss Ennaanay, Guillermo Mendoza, Marc Conte
# for the Natural Capital Project
#
# Last edited: 2/6/2011
#
# First script for Avoided Reservoir Sedimentation
# Calculates biophysical components of sedimentation using the
#    Universal Soil Loss Equation (USLE).  Traces sediment
#    movement downstream to determine sediment retention and
#    export.
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
    try:

        # Make parameters array, and later write input parameter values to an output file
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))

        gp.AddMessage ("\nValidating arguments..." )
        
        # Directory where output files will be written
        gp.workspace = sys.argv[1]
        parameters.append("Workspace: " + gp.workspace)
        
        # Digital elevation model raster
        DEM = sys.argv[2]
        parameters.append("DEM: " + DEM)
        
        # Soil/storm erosivity raster
        Erosivity = sys.argv[3]
        parameters.append("Erosivity: " + Erosivity)
        
        # Soil erodibility raster
        Erodibility = sys.argv[4]
        parameters.append("Erodibility: " + Erodibility)
        
        # Land use / land class raster
        Landuse = sys.argv[5]
        parameters.append("Landcover: " + Landuse)

        # Polygon of whole watershed
        watershed = sys.argv[6]
        parameters.append("Watershed: " + watershed)

        # Sub-watershed polygons
        sub_watersheds = sys.argv[7]
        parameters.append("Sub-watersheds: " + sub_watersheds)

        # Optional v_stream raster
        v_stream_in = sys.argv[8]
        parameters.append("v_stream: " + v_stream_in)
        
        # Table containing model coefficients per landuse class
        Biophys_Table = sys.argv[9]
        parameters.append("Biophysical coefficient table: " + Biophys_Table)

        # Table containing sediment valuation information - critical loading and/or dead volume
        value_table = sys.argv[10]
        parameters.append("Sediment valuation table: " + value_table)

        # Are we valuing dredging?
        value_dredging = sys.argv[11]

        # and/or water quality?
        value_wquality = sys.argv[12]
        
        if value_wquality == '#' and value_dredging == '#':
            gp.AddError("\nError: You must choose to value Dredging and/or Water Quality in the parameter window checkboxes.  Exiting.\n")
            raise Exception

        if value_wquality =='true':
            value_wquality = True
            parameters.append("Value water quality: Yes")
        else:
            value_wquality = False
            parameters.append("Value water quality: No")

        if value_dredging =='true':
            value_dredging = True
            parameters.append("Value dredging: Yes")
        else :
            value_dredging = False
            parameters.append("Value dredging: No")
        
        # Threshold flow accumulation integer
        # How many cells must flow into a cell before it's 
        # considered part of a stream
        threshold_flow_acc = sys.argv[13]
        parameters.append("Threshold flow accumulation: " + str(threshold_flow_acc))
        
        # Slope threshold for length-slope factor
        slope_threshold = sys.argv[14]
        parameters.append("Slope threshold: " + str(slope_threshold))
        
        # Resolution for output rasters
        resolution = sys.argv[15]
        parameters.append("Output resolution: " + str(resolution))
        
        # Suffix to add to end of output filenames, as <filename>_<suffix>
        Suffix = sys.argv[16]
        parameters.append("Suffix: " + Suffix)
        
        if (Suffix == "") or (Suffix == string.whitespace) or (Suffix == "#"):
            Suffix = ""
        else:
            Suffix = "_" + Suffix

    except:
        gp.AddError("\nError in specifying arguments: " + gp.GetMessages(2))
        raise Exception


    # Check and create output folders
    try:
        thefolders=["Output","Intermediate","Service"]
        for folder in thefolders:
            if not gp.exists(gp.workspace + folder):
                gp.CreateFolder_management(gp.workspace, folder)
        pixelfolder = "Pixel"
        pfparent = gp.workspace + os.sep + "Output"
        if not gp.exists(pfparent + folder):
            gp.CreateFolder_management(pfparent, pixelfolder)
    except:
        gp.AddError( "\nError creating folders: " + gp.GetMessages(2))
        raise Exception



    # Output files
    
    try:
        # Output directories
        outputws = gp.workspace + os.sep +  "Output" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        servicews = gp.workspace + os.sep + "Service" + os.sep
        pixelws = gp.workspace + os.sep + "Output" + os.sep + "Pixel" + os.sep

        # Intermediate files
        Lambda = interws + "lambda" 
        c_tmp = interws + "c_tmp"
        p_tmp = interws + "p_tmp"
        c_factor = interws + "c_factor"
        p_factor = interws + "p_factor"
        LS_le = interws + "ls_le"
        LS_gt = interws + "ls_gt"
        LS = interws + "LS"
        RKLS = interws  + "RKLS"
        dir2stream = interws + "dir2stream"
        degree_slope = interws + "deg_slope"
        tmp_ls1 = interws + "tmp_ls1"
        tmp_ls2 = interws + "tmp_ls2"
        tmp_pow1 = interws + "tmp_pow1"
        tmp_pow2 = interws + "tmp_pow2"
        retention_efficiency = interws + "filt_eff"
        frac_removed = interws + "frac_rem"
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
        ws_mask_poly = interws + "ws_mask_poly.shp"
        ws_mask_ras = interws + "ws_mask_ras"
        watersheds_ras = interws + "wshed_ras"
        watershed_loads_table = interws  + "ws_sedexp.dbf"
        sws_USLE_table = interws + "sws_USLE.dbf"
        sws_total_retention_table_wq = interws + "sws_tot_ret_wq.dbf"
        sws_total_retention_table_dr = interws + "sws_tot_ret_dr.dbf"
        sws_ups_retention_table = interws + "sws_ups_ret.dbf"
        sws_export_table = interws + "sws_exp.dbf"
        ws_export_table = interws + "ws_exp.dbf"
        wshed_ras = interws + "wshed_ras"
        wshed_ras2 = interws + "wshed_ras2"
        wshed_ras3 = interws + "wshed_ras3"
        wshed_wq_join = interws + "wshed_wqjoin"
        wshed_dr_join = interws + "wshed_drjoin"
        wsheds_j = interws + "wsheds_j"
        wsheds_totret_wq = interws + "ws_totret_wq"
        wsheds_totret_dr = interws + "ws_totret_dr"
        wsheds_upsret = interws + "wsheds_upsret"
        wsheds_export = interws + "wsheds_export"
        wshed_annual_load = interws + "wshed_annld"
        wshed_num_cells = interws + "ws_numcells"
        wshed_num_cells2 = interws + "ws_numcells2"
        allowed_load_cell_wq = interws + "allow_load_wq"
        allowed_load_cell_dr = interws + "allow_load_dr"
        wshed_dead_vol= interws + "ws_dead_vol"
        wshed_dr_time = interws + "ws_dr_time"
        total_retention_wq1 = interws + "sret_wq1"
        total_retention_dr1 = interws + "sret_dr1"
        watersheds_sjoin = interws + "wsheds_sjoin.shp"
        total_retention = interws + "s_retain1"


        # Input table field names
        # Biophysical Models table
        lucode_field = "lucode"
        c_field = "usle_c"
        p_field = "usle_p"
        retention_efficiency_field = "sedret_eff"

        # Sediment valuation table
        value_id_field = "ws_id"
        wq_ann_load_field = "wq_annload"
        dr_dead_vol_field = "dr_deadvol"
        dr_time_field = "dr_time"

        # ID field for watersheds/sub-watersheds inputs
        wshed_id_field = "ws_id"
        subwshed_id_field = "subws_id"

        # Output files
        v_stream = pixelws + "v_stream"
        USLE = pixelws + "USLE"
        total_retention_wq = pixelws + "sret_wq"
        total_retention_dr = pixelws + "sret_dr"
        # Output tables
        ws_out_table_name = "sediment_export_watershed" + Suffix + ".dbf"
        sws_out_table_name = "sediment_export_subwatershed" + Suffix + ".dbf"
        daccum_raster_fname = "ups_retain"
        export_raster_fname = "s_export"

        # Watershed/sub-watershed outputs
        sws_USLE_mean = outputws + "usle_sws_m"
        sws_USLE_sum = outputws + "usle_sws_s"
        sws_ups_retention_mean = outputws + "upret_sw_m"
        sws_ups_retention_sum = outputws + "upret_sw_s"
        sws_export_mean = outputws + "sexp_sw_m"
        sws_export_sum = outputws + "sexp_sw_s"
        
        # Service outputs
        sws_total_retention_mean_wq = servicews + "sret_sw_qm"
        sws_total_retention_sum_wq = servicews + "sret_sw_qs"
        sws_total_retention_mean_dr = servicews + "sret_sw_dm"
        sws_total_retention_sum_dr = servicews + "sret_sw_ds"

    except:
        gp.AddError( "\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [LS, RKLS, USLE, sws_USLE_mean, sws_USLE_sum, daccum_raster_fname, export_raster_fname, v_stream, total_retention, sws_ups_retention_mean, sws_ups_retention_sum, sws_export_mean, sws_export_sum, sws_total_retention_mean_wq, sws_total_retention_sum_wq, sws_total_retention_mean_dr, sws_total_retention_sum_dr]
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
        LS = str(Outputnames[0])
        RKLS = str(Outputnames[1])
        USLE = str(Outputnames[2])
        sws_USLE_mean = str(Outputnames[3])
        sws_USLE_sum = str(Outputnames[4])
        daccum_raster_fname = str(Outputnames[5])
        export_raster_fname = str(Outputnames[6])
        v_stream = str(Outputnames[7])
        total_retention = str(Outputnames[8])
        sws_ups_retention_mean = str(Outputnames[9])
        sws_ups_retention_sum = str(Outputnames[10])
        sws_export_mean = str(Outputnames[11])
        sws_export_sum = str(Outputnames[12])
        sws_total_retention_mean_wq = str(Outputnames[13])
        sws_total_retention_sum_wq = str(Outputnames[14])
        sws_total_retention_mean_dr= str(Outputnames[15])
        sws_total_retention_sum_dr = str(Outputnames[16])

    except:
        gp.AddError ("\nError validating output filenames: " + gp.GetMessages(2))
        raise Exception
    

    # Check input raster projections - they should all be the same
    try:
        gp.AddMessage("\nValidating input raster projections...")
        DEMDesc = gp.describe(DEM)
        DEMspatref = DEMDesc.SpatialReference
        
        if gp.Exists(v_stream_in):
            rasters = (Erosivity, Erodibility, Landuse, v_stream_in)
        else:
            rasters = (Erosivity, Erodibility, Landuse)
        for x in rasters:
            rasterDesc=gp.describe(x)
            spatreflc = rasterDesc.SpatialReference
            if spatreflc.Type <> 'Projected':
                gp.AddMessage(x + " does not appear to be projected.  It is assumed to be in meters.")
            elif spatreflc.LinearUnitName <> 'Meter':
                gp.AddMessage("This model assumes that data in " + x + " is projected in meters.  You may get erroneous results")
            if str(DEMspatref.name) <> str(spatreflc.name):
                gp.AddError("\nError: " + x + " is not in the same coordinate system as the hydrology rasters.  " + x + " is projected in " + spatreflc.name + " while the hydrology layers are in " + DEMspatref.name + ".  Please project raster in the same projection as the hydrology layers.\n")  
                raise Exception
    except:
        gp.AddError("\nError in validating input raster projections: " + gp.GetMessages(2)) 
        raise Exception

    
    # Set the Geoprocessing environment
    try:
        gp.cellSize = resolution
        gp.Mask = watershed
        install_info = gp.GetInstallInfo("desktop")
        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
    except:
        gp.AddError( "\nError configuring environment: " + gp.GetMessages(2))
        raise Exception



    # Preprocess DEM derivatives and check for hydrologically correct rasters...
    
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
        if (install_info["Version"] == "10.0"):
            gp.BuildRasterAttributeTable_management(Hydrows + "Flow_dir")

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

        # Hydrology rasters
        flow_dir = Hydrows + "Flow_dir"
        flow_acc = Hydrows + "Flow_acc"
        percent_slope = Hydrows + "Slope"
        drop_raster = Hydrows + "Drop_ras"
        
    except:
        gp.AddError("\nError processing hydrology layers: " + gp.GetMessages(2))
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
                if (table_field.Name.upper() == lucode_field.upper()) or (table_field.Name.upper() == p_field.upper()) or (table_field.Name.upper() == p_field.upper())or (table_field.Name.upper() == retention_efficiency_field.upper()) or (table_field.Name.upper() == value_id_field.upper()):
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
        
        
    # Calculate Length-Slope (LS) factor for the USLE
    gp.AddMessage ("\nCalculating Length-Slope (LS) factor..." )
    try:

        # Verify input table fields
        checkfields([lucode_field, c_field, p_field, retention_efficiency_field], Biophys_Table)

        if value_dredging:
            checkfields([value_id_field, dr_dead_vol_field, dr_time_field], value_table)
        if value_wquality:
            checkfields([value_id_field, wq_ann_load_field], value_table)
        
        # Define flow direction to streams, using
        # the input threshold flow accumulation value
        # If a v_stream layer has been entered, use it instead

        if gp.Exists(v_stream_in):
            gp.AddMessage("\nUsing input v_stream")
            gp.CopyRaster(v_stream_in, v_stream)
        else:
            gp.LessThanEqual_sa(flow_acc, threshold_flow_acc, v_stream)

        gp.Times_sa(flow_dir, v_stream, dir2stream)

        # Adjust LS equation value based on slope 
        if percent_slope >= 5:
            NN = 0.5
        elif percent_slope > 3.5 and percent_slope < 5:
            NN = 0.4
        elif percent_slope > 1 and percent_slope <= 3.5:
            NN = 0.3
        else:
            NN = 0.2

        # LS power and multiplication variables
        ls_power = 1.4
        ls_mult = 1.6
        
        # Use when slope <= the input slope threshold
        gp.Slope_sa(DEM, degree_slope, "DEGREE", "1")
        dsc = gp.describe(flow_acc)
        cell_size = str(dsc.MeanCellHeight)
        gp.SingleOutputMapAlgebra_sa(v_stream + " * " + flow_acc + " * " + cell_size + " / 22.13 ", tmp_ls1)
        gp.Power_sa(tmp_ls1, NN, tmp_pow1)
        gp.SingleOutputMapAlgebra_sa("SIN(" + degree_slope + " * 0.01745) / 0.09", tmp_ls2)
        gp.Power_sa(tmp_ls2, ls_power, tmp_pow2)
        gp.SingleOutputMapAlgebra_sa("(" + tmp_pow1 + " * " + tmp_pow2 + ") * " + str(ls_mult), LS_le)

        # Use when slope > the input slope threshold
        gp.SingleOutputMapAlgebra_sa("CON(" + dir2stream + " == 1 OR " + dir2stream + " == 4 OR " + dir2stream + " == 16 OR " + dir2stream + " == 64, " + cell_size + ", 1.4 * " + cell_size + ")", Lambda)
        gp.SingleOutputMapAlgebra_sa("0.08 * POW(" + Lambda + ", 0.35) * POW(" + percent_slope + ", 0.6)", LS_gt)

        # Assign LS values appropriately based on percent slope raster
        gp.SingleOutputMapAlgebra_sa("CON(" + percent_slope + " > " + slope_threshold + ", " + LS_gt + ", " + LS_le + ")", LS) 

    except:
        gp.AddError("\nError in calculating Length-Slope (LS) factor:  " + gp.GetMessages(2))
        raise Exception

    # Calculating RKLS, the biophysical component of USLE
    # Change units from tons per hectare to tons per cell = cell size ^2 / 10000
    gp.AddMessage ("\nCalculating RKLS...")
    try:
        gp.SingleOutputMapAlgebra_sa("(" + Erosivity + " * " + LS + " * " + Erodibility + ") * " + cell_size + " * " + cell_size + " / 10000", RKLS)
     
    except:
        gp.AddError("\nError calculating RKLS:  " + gp.GetMessages(2))
        raise Exception

        
    # Calculating Universal Soil Loss Equation (USLE)
    gp.AddMessage("\nCalculating USLE...") 
    try:
        # Get Practice/Cropping values from input table
        gp.ReclassByTable_sa(Landuse, Biophys_Table, lucode_field, lucode_field, c_field, c_tmp, "DATA")
        gp.ReclassByTable_sa(Landuse, Biophys_Table, lucode_field, lucode_field, p_field, p_tmp, "DATA")
        # Divide c and p values by 1000, since table values are x1000 to create integers,
        # which are needed by Reclass
        gp.SingleOutputMapAlgebra_sa("FLOAT(" + c_tmp + ") / 1000", c_factor)
        gp.SingleOutputMapAlgebra_sa("FLOAT(" + p_tmp + ") / 1000", p_factor)
        gp.SingleOutputMapAlgebra_sa(c_factor + " * " + p_factor + " * " + RKLS, USLE)

        # Aggregate by watersheds and sub-watersheds

        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, USLE, sws_USLE_sum, "SUM", "DATA")

        area_ha_field = "AREA_HA"
        mean_ha_field = "MEAN_HA"
        
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, USLE, sws_USLE_table, "DATA")
        gp.AddField(sws_USLE_table, area_ha_field, "double")
        gp.AddField(sws_USLE_table, mean_ha_field, "double")
        gp.CalculateField_management(sws_USLE_table, area_ha_field, "[AREA] / 10000", "VB")
        gp.CalculateField_management(sws_USLE_table, mean_ha_field, "[SUM] / [" + area_ha_field + "]", "VB")

        gp.FeatureToRaster_conversion(sub_watersheds, subwshed_id_field, wshed_ras, cell_size)
        gp.MakeRasterLayer_management(wshed_ras, "wsheds", "#", "#")
        gp.MakeTableView(sws_USLE_table, "usle_table_view")
        
        # Zonal stats field name has changed in Arc 10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = subwshed_id_field
        else:
            zstat_id_field = "VALUE"
            
        gp.AddJoin_management("wsheds", "VALUE", "usle_table_view", zstat_id_field)
        gp.CopyRaster_management("wsheds", wsheds_j)
        
        # Delete this now or Arc10 won't delete the Intermediate folder at the end
        if (install_info["Version"] == "10.0"):
            gp.Delete_management("wsheds")
            gp.Delete_management("usle_table_view")

        gp.Lookup_sa(wsheds_j, mean_ha_field, sws_USLE_mean)

        gp.AddMessage("\nCreated sub-watershed USLE outputs:\n\t" + str(sws_USLE_sum) + "\n\t" + str(sws_USLE_mean))

    except:
        gp.AddError("\nError calculating USLE: " + gp.GetMessages(2))
        raise Exception


#
# Calculate the amount of sediment that is retained by each cell
# from what passes through the cell from upstream (daccum)
#
# Also keep track of how much of each cell's sediment load makes it
# to the stream (export)
#
# Loop through sub-watersheds and process each separately,
# recombining them afterward to allow for arbitrarily large
# watersheds
#

    gp.AddMessage("\nCalculating sediment removal...")
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
            # TO DO: ADD SUPPORT FOR QUERYING GEODATABASES select_exp = "[subws_id]"
            subws_id = str(int(ws_row.GetValue(subwshed_id_field)))
            select_exp = "\"subws_id\" = " + subws_id

            # Use cell size from USLE raster
            dsc = gp.describe(USLE)
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
                    if not gp.exists(folder):
                        gp.CreateFolder_management(interws, folder)

                export_raster_subws = interws + export_sub_wsheds_dir + os.sep + "s_exp"
                daccum_grid_raster_subws = interws + daccum_sub_wsheds_dir + os.sep + "s_ret"

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

                # Create raster of fraction of sediment removed at each cell, from input table
                gp.ReclassByTable_sa(Landuse, Biophys_Table, lucode_field, lucode_field, retention_efficiency_field, retention_efficiency, "DATA")
                # Divide by 100 to change input table integer (as required by Reclass) percent values to decimal
                gp.SingleOutputMapAlgebra_sa("FLOAT(" + retention_efficiency + ") / 100.0", frac_removed)

                # Make sure all input rasters are the same resolution
                # as the loads (USLE) raster and resample if necessary
                try:
                    gp.AddMessage("\n\tValidating input raster resolution...")
                    rasterDescL=gp.describe(USLE)
                    loads_cell_width = rasterDescL.MeanCellWidth
                    
                    rasterDescFD=gp.describe(dir2stream)
                    if rasterDescFD.MeanCellWidth <> loads_cell_width:
                        gp.AddMessage("\n\tResampling flow direction...")
                        gp.Resample_management(dir2stream, flowdir_rs, loads_cell_width)
                    else:
                        flowdir_rs = dir2stream
                    rasterDescFR=gp.describe(frac_removed)
                    if rasterDescFR.MeanCellWidth <> loads_cell_width:
                        gp.AddMessage("\n\tResampling removal efficiency raster...")
                        gp.Resample_management(frac_removed, frac_removed_rs, loads_cell_width)
                    else:
                        frac_removed_rs = frac_removed
                
                except:
                    gp.AddError("\nError in validating input raster resolution: " + gp.GetMessages(2)) 
                    raise Exception

                # Extract input rasters, using the input mask
                # so that all extents are the same

                # Sediment loading at each cell
                gp.ExtractByMask_sa(USLE, ws_mask_poly, loads_ext)
                # Flow direction 
                gp.ExtractByMask_sa(flowdir_rs, ws_mask_poly, flowdir_ext)
                # Fraction of incoming sediment retained by each cell
                gp.ExtractByMask_sa(frac_removed_rs, ws_mask_poly, frac_removed_ext)

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
                continue
#                raise Exception

            # Follow flow path of each cell's load,
            #   calculating sediment retention along the way
            # Flow tracing based on logic from ArcScript AS15793 - Tracegrid
            #   by Christopher C. Duncan, duncan@geo.umass.edu         

            try:

                # Flow direction loop debugging
                #
                # If it looks like there's a circular flow direction problem, uncomment the following three lines,
                # as well as two others later in the code beginning with 'outfile', also with the comment
                # 'Flow direction loop debugging'
                #
                # This will simply create a file where data values for each cell are written as it is processed.
                # Keep an eye on the last lines of this file as the script runs, and look for the same set of 
                # values repeating - this indicates a circular flow path.

                # Uncomment these three lines
                # outfile_name = gp.workspace + "\\Output\\sed_loop_debug_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt"
                # outfile = open(outfile_name, "w")
                # gp.AddMessage("\nFlow direction loop debugging is on, output file is " + str(outfile_name) + "\n")


                # Time sub-watershed process.  If it's taking too long, break out, it's probably in an infinite loop.
                start_time = time.clock()

                # Loop over all cells in the loads grid
                gp.AddMessage("\n\tCalculating sediment retention per cell...")
                
                for y in range (0, numrows_loads):
                    for x in range (0, numcols_loads):

                        # Time sub-watershed process for each cell's trace down its flow path.
                        # If it's taking too long, break out, it's probably in an infinite loop.
                        #start_time = time.clock()
                        
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
                                if (first_cell == 1):
                                    cur_x = x
                                    cur_y = y
                                else:
                                    # Shift down one cell
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

                                    # Keep track of how much upstream sediment is retained by the cell
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
                                    # Keep track of how much of each cell's sediment makes it to
                                    # a stream, stored in the export array
                                    export_array[y][x] = load_passed

                                # or this sub-watershed has been processing for too long (45 minutes by default)
                                # - probably indicates that there's a circular flow path
                                # To check for circular flow paths, look for lines beginning with 'outfile'
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
                gp.AddError("\nError calculating sediment retention per cell: " + gp.GetMessages(2))
                raise Exception


            # Process export and retention output files
            try:
                gp.AddMessage("\n\tProcessing sediment retention output rasters...")

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
                gp.AddError("\nError processing sediment retention output rasters: " + gp.GetMessages(2))
                raise Exception

    except:
        gp.AddError("\nError processing sub-watershed id " + str(subws_id) + ": " + gp.GetMessages(2))
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


    # Total sediment retained, from upslope plus what originates in the cell itself
    gp.AddMessage("\nCalculating total sediment retention...")
    try:
        # Only include C factor here (not P), as we're evaluating the landscape, not farming practices
        gp.SingleOutputMapAlgebra_sa("( " + daccum_raster + " + ( " + RKLS + " - (" + RKLS + " * " + c_factor + " )) )", total_retention)
                                                                                                                                                                           
        # Subtract the amount of sediment allowed by the critical load and/or dead volume input
        if value_wquality:

            gp.AddMessage("\n\tCalculating sediment allowance for water quality...")

            # Get sediment values from table, map to watershed
            # Copy to dbf to work around table OID requirement change in ArcMap 9.3
            dsc = gp.describe(total_retention)
            cell_size = str(dsc.MeanCellHeight)

            gp.CopyRows_management(value_table, "sed_tmp_wq.dbf")
            gp.FeatureToRaster_conversion(watershed, wshed_id_field, wshed_ras2, cell_size)
            gp.MakeRasterLayer_management(wshed_ras2, "wshed_tmp_wq")
            gp.AddJoin_management("wshed_tmp_wq", "Value", "sed_tmp_wq.dbf", value_id_field)
            gp.CopyRaster_management("wshed_tmp_wq", wshed_wq_join)
            # Load allowed per cell = table wq_annload / number of cells in watershed
            gp.Lookup_sa(wshed_wq_join, wq_ann_load_field, wshed_annual_load)
            gp.Lookup_sa(wshed_wq_join, "COUNT", wshed_num_cells)
            gp.Divide_sa(wshed_annual_load, wshed_num_cells, allowed_load_cell_wq)
            gp.Minus_sa(total_retention, allowed_load_cell_wq, total_retention_wq1)
            # Change negative values to zeroes
            gp.SingleOutputMapAlgebra_sa("CON(" + total_retention_wq1 + " < 0, 0, " + total_retention_wq1 + ")", total_retention_wq)
                     
        if value_dredging:

            gp.AddMessage("\n\tCalculating sediment allowance for dredging...")

            # Get sediment values from table, map to watershed
            # Copy to dbf to work around table OID requirement change in ArcMap 9.3
            dsc = gp.describe(total_retention)
            cell_size = str(dsc.MeanCellHeight)

            gp.CopyRows_management(value_table, "sed_tmp_dr.dbf")
            gp.FeatureToRaster_conversion(watershed, wshed_id_field, wshed_ras3, cell_size)
            gp.MakeRasterLayer_management(wshed_ras3, "wshed_tmp_dr")
            gp.AddJoin_management("wshed_tmp_dr", "Value", "sed_tmp_dr.dbf", value_id_field)
            gp.CopyRaster_management("wshed_tmp_dr", wshed_dr_join)
            gp.Lookup_sa(wshed_dr_join, dr_dead_vol_field, wshed_dead_vol)
            gp.Lookup_sa(wshed_dr_join, dr_time_field, wshed_dr_time)
            gp.Lookup_sa(wshed_dr_join, "COUNT", wshed_num_cells2)
            gp.SingleOutputMapAlgebra_sa("(" + wshed_dead_vol + " * 1.26) / (" + wshed_num_cells2 + " * " + wshed_dr_time + ")", allowed_load_cell_dr)
            gp.Minus_sa(total_retention, allowed_load_cell_dr, total_retention_dr1)
            # Change negative values to zeroes
            gp.SingleOutputMapAlgebra_sa("CON(" + total_retention_dr1 + " < 0, 0, " + total_retention_dr1 + ")", total_retention_dr)
        
    except:
        gp.AddError ("\nError calculating total sediment retention: " + gp.GetMessages(2))
        raise Exception


    # Aggregate by sub-watersheds
    try:

        if (value_wquality):
            gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, total_retention_wq, sws_total_retention_sum_wq, "SUM", "DATA")
            gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, total_retention_wq, sws_total_retention_table_wq, "DATA")
            gp.AddField(sws_total_retention_table_wq, area_ha_field, "double")
            gp.AddField(sws_total_retention_table_wq, mean_ha_field, "double")
            gp.CalculateField_management(sws_total_retention_table_wq, area_ha_field, "[AREA] / 10000", "VB")
            gp.CalculateField_management(sws_total_retention_table_wq, mean_ha_field, "[SUM] / [" + area_ha_field + "]", "VB")
            gp.MakeRasterLayer_management(wshed_ras, "wsheds_tr_wq", "#", "#")
            gp.MakeTableView(sws_total_retention_table_wq, "totret_wq_table_view")
            gp.AddJoin_management("wsheds_tr_wq", "VALUE", "totret_wq_table_view", zstat_id_field)
            gp.CopyRaster_management("wsheds_tr_wq", wsheds_totret_wq)
            gp.Lookup_sa(wsheds_totret_wq, mean_ha_field, sws_total_retention_mean_wq)

        if (value_dredging):
            gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, total_retention_dr, sws_total_retention_sum_dr, "SUM", "DATA")
            gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, total_retention_dr, sws_total_retention_table_dr, "DATA")
            gp.AddField(sws_total_retention_table_dr, area_ha_field, "double")
            gp.AddField(sws_total_retention_table_dr, mean_ha_field, "double")
            gp.CalculateField_management(sws_total_retention_table_dr, area_ha_field, "[AREA] / 10000", "VB")
            gp.CalculateField_management(sws_total_retention_table_dr, mean_ha_field, "[SUM] / [" + area_ha_field + "]", "VB")
            gp.MakeRasterLayer_management(wshed_ras, "wsheds_tr_dr", "#", "#")
            gp.MakeTableView(sws_total_retention_table_dr, "totret_dr_table_view")
            gp.AddJoin_management("wsheds_tr_dr", "VALUE", "totret_dr_table_view", zstat_id_field)
            gp.CopyRaster_management("wsheds_tr_dr", wsheds_totret_dr)
            gp.Lookup_sa(wsheds_totret_dr, mean_ha_field, sws_total_retention_mean_dr)


        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, daccum_raster, sws_ups_retention_sum, "SUM", "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, daccum_raster, sws_ups_retention_table, "DATA")
        gp.AddField(sws_ups_retention_table, area_ha_field, "double")
        gp.AddField(sws_ups_retention_table, mean_ha_field, "double")
        gp.CalculateField_management(sws_ups_retention_table, area_ha_field, "[AREA] / 10000", "VB")
        gp.CalculateField_management(sws_ups_retention_table, mean_ha_field, "[SUM] / [" + area_ha_field + "]", "VB")
        gp.MakeRasterLayer_management(wshed_ras, "wsheds_ups", "#", "#")
        gp.MakeTableView(sws_ups_retention_table, "ups_table_view")
        gp.AddJoin_management("wsheds_ups", "VALUE", "ups_table_view", zstat_id_field)
        gp.CopyRaster_management("wsheds_ups", wsheds_upsret)
        gp.Lookup_sa(wsheds_upsret, mean_ha_field, sws_ups_retention_mean)
        
        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, export_raster, sws_export_sum, "SUM", "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, export_raster, sws_export_table, "DATA")
        gp.AddField(sws_export_table, area_ha_field, "double")
        gp.AddField(sws_export_table, mean_ha_field, "double")
        gp.CalculateField_management(sws_export_table, area_ha_field, "[AREA] / 10000", "VB")
        gp.CalculateField_management(sws_export_table, mean_ha_field, "[SUM] / [" + area_ha_field + "]", "VB")
        gp.MakeRasterLayer_management(wshed_ras, "wsheds_exp", "#", "#")
        gp.MakeTableView(sws_export_table, "exp_table_view")
        gp.AddJoin_management("wsheds_exp", "VALUE", "exp_table_view", zstat_id_field)
        gp.CopyRaster_management("wsheds_exp", wsheds_export)
        gp.Lookup_sa(wsheds_export, mean_ha_field, sws_export_mean)
        
        gp.AddMessage("\nCreated sub-watershed outputs:")
        if value_wquality:
            gp.AddMessage("\t" + str(sws_total_retention_sum_wq) + "\n\t" + str(sws_total_retention_mean_wq))
        if value_dredging:
            gp.AddMessage("\t" + str(sws_total_retention_sum_dr) + "\n\t" + str(sws_total_retention_mean_dr))
        gp.AddMessage("\t" + str(sws_ups_retention_sum) + "\n\t" + str(sws_ups_retention_mean))
        gp.AddMessage("\t" + str(sws_export_sum) + "\n\t" + str(sws_export_mean))

    except:
        gp.AddError ("\nError aggregating by sub-watersheds: " + gp.GetMessages(2))
        raise Exception

    

    # Table to hold generated output values per watershed
    gp.AddMessage("\nCreating watershed loads output table...")
    try:
        # output table field names
        out_table_wshed_id_field = "ws_id"
        out_table_load_field = "sed_load"
        
        gp.CreateTable_management(outputws, ws_out_table_name)
        ws_out_table = outputws + ws_out_table_name
        
        gp.AddField(ws_out_table, out_table_wshed_id_field, "long")
        gp.AddField(ws_out_table, out_table_load_field, "double")
        # Remove Field1 - it's added by default and not used
        gp.DeleteField_management(ws_out_table, "Field1")

        ws_out_table_rows = gp.InsertCursor(ws_out_table)

        gp.ZonalStatisticsAsTable_sa(watershed, wshed_id_field, export_raster, ws_export_table, "DATA")

        # Populate it with Watershed Loads values
        wse_rows = gp.SearchCursor(ws_export_table)
        wse_row = wse_rows.Reset
        wse_row = wse_rows.Next()

        while (wse_row):
            # Zonal stats field name has changed in Arc 10
            if (install_info["Version"] == "10.0"):
                ws_id = wse_row.getValue("ID")
            else:
                ws_id = wse_row.getValue("VALUE")
                
            ws_total_load = float(wse_row.getValue("SUM"))
            
            new_row = ws_out_table_rows.NewRow()
            new_row.setValue(out_table_wshed_id_field, ws_id)
            new_row.setValue(out_table_load_field, ws_total_load)

            ws_out_table_rows.InsertRow(new_row)
            wse_row = wse_rows.Next()

        gp.AddMessage("\tCreated watershed loads output table: \n\t" + str(ws_out_table))
        
    except:
        gp.AddError ("\nError creating watershed loads output table:" + gp.GetMessages(2))
        raise Exception

    # Table to hold generated output values per sub-watershed
    gp.AddMessage("\nCreating sub-watershed loads output table...")
    try:
        # output table field names
        out_table_wsid_field = "ws_id"
        out_table_subwshed_id_field = "subws_id"
        out_table_load_field = "sed_load"
        
        gp.CreateTable_management(outputws, sws_out_table_name)
        sws_out_table = outputws + sws_out_table_name
        
        gp.AddField(sws_out_table, out_table_subwshed_id_field, "long")
        gp.AddField(sws_out_table, out_table_wsid_field, "long")
        gp.AddField(sws_out_table, out_table_load_field, "double")
        # Remove Field1 - it's added by default and not used
        gp.DeleteField_management(sws_out_table, "Field1")

        sws_out_table_rows = gp.InsertCursor(sws_out_table)

        # Map sub-watersheds to their corresponding watersheds
        gp.SpatialJoin_analysis(sub_watersheds, watershed, watersheds_sjoin, "JOIN_ONE_TO_ONE", "KEEP_ALL", "#", "IS_WITHIN")

        # Populate output table with watershed load values
        swse_rows = gp.SearchCursor(sws_export_table, "", "", "", zstat_id_field + " A")
        swse_row = swse_rows.Reset
        swse_row = swse_rows.Next()

        while (swse_row):

            sj_rows = gp.SearchCursor(watersheds_sjoin, "", "", "", subwshed_id_field + " A")
            sj_row = sj_rows.Reset
            sj_row = sj_rows.Next()

            while (int(sj_row.getValue(subwshed_id_field)) <> int(swse_row.getValue(zstat_id_field))):
                sj_row = sj_rows.Next()
            
            # Zonal stats field name has changed in Arc 10
            if (install_info["Version"] == "10.0"):
                sws_id = swse_row.getValue("ID")
            else:
                sws_id = swse_row.getValue("VALUE")
                
            sws_total_load = float(swse_row.getValue("SUM"))
            
            new_row = sws_out_table_rows.NewRow()
            new_row.setValue(out_table_subwshed_id_field, sws_id)
            new_row.setValue(out_table_wsid_field, int(sj_row.getValue(wshed_id_field)))
            new_row.setValue(out_table_load_field, sws_total_load)

            sws_out_table_rows.InsertRow(new_row)
            swse_row = swse_rows.Next()
            del sj_row, sj_rows

        gp.AddMessage("\tCreated sub-watershed loads output table: \n\t" + str(sws_out_table))
        
    except:
        gp.AddError ("\nError creating output tables:" + gp.GetMessages(2))
        raise Exception


    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Output\\Sediment_1_Soil_Loss_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("SEDIMENT 1 - SOIL LOSS PARAMETERS\n")
        parafile.writelines("_________________________________\n\n")

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
        del ws_row, ws_rows, ws_out_table_rows, new_row, fd_row, fd_rows
        del wse_row, wse_rows, sws_out_table_rows, swse_row, swse_rows
        gp.Delete_management(interws)
    except:
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception


except:
    gp.AddError ("\nError running script")
    raise Exception


