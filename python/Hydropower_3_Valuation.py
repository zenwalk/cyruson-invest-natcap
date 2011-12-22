#-----------------------------------------------------------------
#
# Hydropower_Production_Reservoir.py
#
# Coded by :
# Stacie Wolny, Bradley Eichelberger
# Written by:
# Driss Ennaanay, Guillermo Mendoza, Marc Conte
# for the Natural Capital Project
#
# Last edit: 5/2/2011
#
# Calculates the value of a landscape for providing water yield
# for reservoir hydropower production.
#
#-----------------------------------------------------------------

# Import system modules
import sys, string, os, arcgisscripting, math, time, datetime

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out necessary licenses
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

        # Calibrated water yield (volume)
        cyield_vol = gp.GetParameterAsText(1)
        parameters.append("Calibrated water yield (volume): " + cyield_vol)

        # Water consumption (volume)
        consump_vol = gp.GetParameterAsText(2)
        parameters.append("Water consumption (volume): " + consump_vol)

        # Watersheds shapefile
        watersheds = gp.GetParameterAsText(3)
        parameters.append("Watersheds: " + watersheds)

        # Sub-watersheds shapefile
        sub_watersheds = gp.GetParameterAsText(4)
        parameters.append("Sub-watersheds: " + sub_watersheds)

        # Watershed scarcity output table from Scarcity script
        ws_scarcity_table = gp.GetParameterAsText(5)
        parameters.append("Watershed scarcity table: " + ws_scarcity_table)
        
        # Sub-watershed scarcity output table from Scarcity script
        sws_scarcity_table = gp.GetParameterAsText(6)
        parameters.append("Sub-watershed scarcity table: " + sws_scarcity_table)

        # Hydropower station table 
        station_table = gp.GetParameterAsText(7)
        parameters.append("Hydropower table: " + station_table)

        # Output resolution
        resolution = gp.GetParameterAsText(8)
        parameters.append("Output resolution: " + str(resolution))
        
        # Results suffix
        Suffix = gp.GetParameterAsText(9)
        parameters.append("Suffix: " + Suffix)

        if (Suffix == "") or (Suffix == string.whitespace) or (Suffix == "#"):
            Suffix = ""
        else:
            Suffix = "_" + Suffix 

    except:
        gp.AddError("\nError in specifying arguments: " + gp.GetMessages(2))
        raise Exception


    # Check and create folders

    try:
        thefolders=["Output","Intermediate","Service"]
        for folder in thefolders:
            if not gp.exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
                
    except:
        gp.AddError("\nError creating output folders: " + gp.GetMessages(2))
        raise Exception


    # Temporary files

    # Base output/working directories
    outputws = gp.workspace + os.sep + "Output" + os.sep
    interws = gp.workspace + os.sep + "Intermediate" + os.sep
    servicews = gp.workspace + os.sep + "Service" + os.sep

    # Intermediate variables
    hp_energy = interws + "energy"
    hp_wsupply = interws + "wsupply"
    hp_wconsump = interws + "wconsump"
    wsheds_tmp = interws + "hpsheds_tmp"
    wsheds_j = interws + "hp_sheds_j"
    wsheds_j2 = interws + "hp_sheds_j2"
    yld_m3 = interws + "yld_m3"
    yld2 = interws + "yld2"
    yld3 = interws + "yld3"
    consump = interws + "consump"
    consump0 = interws + "consump0"
    station_wsupply = interws + "wsupply_st"
    dam_timespan = interws + "dam_ts"
    ws_value1 = interws + "ws_val1"
    ws_energy1 = interws + "ws_en1"
    ws_energy2 = interws+ "ws_en2"
    stable_tmp = interws + "stable_tmpj1"
    wshed_ras = interws + "wshed_ras_hp"
    rsupply_vol = interws + "rsup_vol"
    rsupply_zstat = interws + "rsupply_zstat.dbf"
    cyield_zstat = interws + "cyield_zstat.dbf"
    consump_zstat = interws + "consump_zstat.dbf"
    rsupply_tot = interws + "rsupply_tot"
    ws_rsupply = interws + "ws_rsupply"
    subws_rsupply_fract = interws + "subws_fract"
    subws_energy_ann = interws + "sws_en_ann"
    subws_energy_zstat = interws + "subws_energy_zstat.dbf"
    subws_value_zstat = interws + "subws_value_zstat.dbf"
    watersheds_sjoin = interws + "wsheds_sjoin.shp"
    sub_watersheds_join = interws + "subws_join.shp"
    wsheds_sjoin_copy = interws + "wsheds_sjoin_copy.shp"
    ws_value = interws + "hp_val_ws"
    ws_energy = interws + "hp_en_ws"
    
    # Fields in hydropower station table
    station_id_field = "ws_id"
    kwval_field = "kw_price"
    height_field = "height"
    efficiency_field = "efficiency"
    fraction_field = "fraction"
    time_field = "time_span"
    discount_field = "discount"
    cost_field = "cost"
    # Water demand table
    lucode_field = "lucode"
    demand_field = "demand"
    # Output table, created in Scarcity script
    out_table_swid_field = "subws_id"

    # ID fields for watersheds/sub-watersheds inputs
    wshed_id_field = "ws_id"
    subwshed_id_field = "subws_id"
    
    # Output files

    subws_value = servicews + "hp_val"
    subws_energy = servicews + "hp_energy"
    hp_ws_table_name = "hydropower_value_watershed" + Suffix + ".dbf"
    hp_subws_table_name = "hydropower_value_subwatershed" + Suffix + ".dbf"


    # Add suffix to end of output filenames
    # Verify length of output rasters does not exceed 13 characters
    try:
        Outputnames = [subws_value, subws_energy]
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
        subws_value = str(Outputnames[0])
        subws_energy = str(Outputnames[1])
    except:
        gp.AddError ("\nError validating output filenames: " + gp.GetMessages(2))
        raise Exception


    # Set the Geoprocessing environment 

    try:
        gp.cellSize = resolution
        install_info = gp.GetInstallInfo("desktop")
        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws

    except:
        gp.AddError( "\nError configuring output resolution: " + gp.GetMessages(2))
        raise Exception
    

    # Verify required input table fields' existence and type

    def checkfields (fields, table):
        
        try:

            table_fields = gp.listfields(table, "*", "All")
            table_field = table_fields.next()

            foundfields = []
            while (table_field):
                foundfields.append(table_field.Name.upper())
                # Make sure that all fields that need to be integer are
                if (table_field.Name.upper() == lucode_field.upper()) or (table_field.Name.upper() == demand_field.upper()) or (table_field.Name.upper() == station_id_field.upper()):
                    if not "Integer" in table_field.Type:
                        gp.AddError("\nError: Field " + str(table_field.Name) + " in table " + str(table) + " must be of type Integer\n")
                        raise Exception
                table_field = table_fields.next()

            for f in fields:
                if not f.upper() in foundfields:
                    gp.AddError("\nError: Required table field " + str(f) + " not found in input table " + str(table) + "\n")
                    raise Exception

        except:
            gp.AddError ("\nError verifying input table fields" + gp.GetMessages(2))
            raise Exception
    

    # Create output tables - one for sub-basins, one for basins
    try:

        gp.AddMessage("\nCreating output tables...")

        # Table for whole watershed - add value fields to watershed scarcity table
        gp.CreateTable_management(servicews, hp_ws_table_name)
        hp_ws_table = servicews + hp_ws_table_name
        gp.CopyRows_management(ws_scarcity_table, hp_ws_table)

        # output table field names
        out_table_sid_field = "ws_id"
        out_table_value_field = "hp_value"
        out_table_energy_field = "hp_energy"
        out_table_rsupply_field = "rsupply_vl"

        gp.AddField(hp_ws_table, out_table_energy_field, "double")
        gp.AddField(hp_ws_table, out_table_value_field, "double")

        # Table for sub-watersheds - add value fields to sub-watershed scarcity table

        gp.CreateTable_management(servicews, hp_subws_table_name)
        hp_subws_table = servicews + hp_subws_table_name
        gp.CopyRows_management(sws_scarcity_table, hp_subws_table)

        gp.AddField(hp_subws_table, out_table_energy_field, "double")
        gp.AddField(hp_subws_table, out_table_value_field, "double")

    except:
        gp.AddError("\nError creating output tables: " + gp.GetMessages(2))
        raise Exception


    # Calculate Present Value for each hydropower station
    try:

        gp.AddMessage("\nCalculating present value...")

        # Verify input table fields
        checkfields([station_id_field, kwval_field, height_field, efficiency_field, fraction_field, time_field, discount_field, cost_field], station_table)

        # Total water supply and consumption for each watershed come from scarcity table
        y_rows = gp.SearchCursor(sws_scarcity_table)
        y_row = y_rows.Reset
        y_row = y_rows.Next()

        # Realized water supply volume per sub-watershed
        gp.Minus_sa(cyield_vol, consump_vol, rsupply_vol)

        wshed_rows = gp.SearchCursor(watersheds)
        wshed_row = wshed_rows.Reset
        wshed_row = wshed_rows.Next()

        # Map sub-watersheds to their corresponding watersheds using Spatial Join
        # Will sum volumes from sub-watersheds to create the Vin for corresponding watershed
        gp.MakeTableView_management(sws_scarcity_table, "scarcity_view")

        gp.SpatialJoin_analysis(sub_watersheds, watersheds, watersheds_sjoin, "JOIN_ONE_TO_ONE", "KEEP_ALL", "#", "IS_WITHIN")
        gp.MakeFeatureLayer_management(watersheds_sjoin, "wsheds_sjoin_layer")
        gp.AddJoin_management("wsheds_sjoin_layer", subwshed_id_field, "scarcity_view", out_table_swid_field)
        gp.CopyFeatures_management("wsheds_sjoin_layer", wsheds_sjoin_copy)

        if (install_info["Version"] == "10.0"):
            gp.Delete_management("scarcity_view")
            gp.Delete_management("wsheds_sjoin_layer")

        while (wshed_row):

            gp.AddMessage("\nCalculating value for station id " + str(wshed_row.getValue(wshed_id_field)) + "...")

            stable_rows = gp.SearchCursor(station_table)
            stable_row = stable_rows.Reset
            stable_row = stable_rows.Next()

            # Match hydropower station table watershed to a row in the watershed shapefile
            while (int(stable_row.getValue(station_id_field)) <> int(wshed_row.getValue(wshed_id_field))):
                stable_row = stable_rows.Next()

            # Calculate energy produced and value of energy produced at each 
            # hydropower station over the input time period

            station_efficiency = float(stable_row.getValue(efficiency_field))
            station_fraction = float(stable_row.getValue(fraction_field))
            station_height = float(stable_row.getValue(height_field))

            ws_sc_rows = gp.SearchCursor(ws_scarcity_table)
            ws_sc_row = ws_sc_rows.Reset
            ws_sc_row = ws_sc_rows.Next()

            # Match watershed in scarcity table to row in the watershed shapefile
            while (int(ws_sc_row.getValue(station_id_field)) <> int(wshed_row.getValue(wshed_id_field))):
                ws_sc_row = ws_sc_rows.Next()

            # Vin = realized water supply (yield - demand)
            Vin = float(ws_sc_row.getValue(out_table_rsupply_field))

            # Energy production 
            energy = station_efficiency * station_fraction * station_height * Vin * .00272

            # Output table
            o_rows = gp.UpdateCursor(hp_ws_table)
            o_row = o_rows.Reset
            o_row = o_rows.Next()

            # Match watershed in the output table with row in watershed shapefile
            while (int(o_row.getValue(station_id_field)) <> int(wshed_row.getValue(wshed_id_field))):
                o_row = o_rows.Next()
            
            station_kwval = float(stable_row.getValue(kwval_field))
            station_time = float(stable_row.getValue(time_field))
            station_discount = float(stable_row.getValue(discount_field))
            station_cost = float(stable_row.getValue(cost_field))

            # Calculate Present Value 

            dsum = 0

            for t in range (0, station_time):
                dsum += 1 / pow(1 + (station_discount / 100), t)

            NPV = ((station_kwval * energy) -  station_cost) * dsum
                        
            # Add new field values to output table
            o_row.setValue(out_table_energy_field, float(energy))
            o_row.setValue(out_table_value_field, float(NPV))
            o_rows.UpdateRow(o_row)
            wshed_row = wshed_rows.Next()

            del stable_row, stable_rows, ws_sc_row, ws_sc_rows
                        
    except:
        gp.AddError("\nError calculating present value: " + gp.GetMessages(2))
        raise Exception

    # Final values = energy and value over lifetime of dam

    try:

        gp.AddMessage("\nCreating watershed output...")
        
        dsc = gp.describe(cyield_vol)
        cell_size = str(dsc.MeanCellHeight)
        
        gp.FeatureToRaster_conversion(watersheds, wshed_id_field, wshed_ras, cell_size)
        gp.BuildRasterAttributeTable_management(wshed_ras)
        gp.MakeRasterLayer_management(wshed_ras, "wsheds", "#", "#")

        gp.CreateTable_management(interws, "test_hp")
        hp_test_table = interws + "test_hp"
        gp.CopyRows_management(hp_ws_table, hp_test_table)
        gp.MakeTableView(hp_test_table, "test_table_view")
        gp.AddJoin_management("wsheds", "VALUE", "test_table_view", out_table_sid_field)
        gp.CopyRaster_management("wsheds", wsheds_j)

        if (install_info["Version"] == "10.0"):
            gp.Delete_management("test_table_view")

        # Map watersheds to calculated energy/value to do subsequent raster math
        gp.Lookup_sa(wsheds_j, out_table_energy_field, ws_energy1)
        gp.Lookup_sa(wsheds_j, out_table_value_field, ws_value1)
        gp.Lookup_sa(wsheds_j, out_table_rsupply_field, ws_rsupply)

        if (install_info["Version"] == "10.0"):
            gp.Delete_management(wsheds_j)
            gp.Delete_management("wsheds")

        # Map dam life spans to corresponding watersheds
        gp.MakeRasterLayer_management(wshed_ras, wsheds_tmp, "#", "#")
        gp.CopyRows_management(station_table, stable_tmp)
        gp.AddJoin_management(wsheds_tmp, "VALUE", stable_tmp, station_id_field)
        gp.CopyRaster_management(wsheds_tmp, wsheds_j2)
        gp.Lookup_sa(wsheds_j2, time_field, dam_timespan)
        
        gp.Delete_management(stable_tmp)
        gp.Delete_management(wsheds_tmp)

        # Energy per watershed over lifetime of dam
        gp.Times_sa(ws_energy1, dam_timespan, ws_energy2)
        
        # Change negative values to zeroes
        gp.SingleOutputMapAlgebra_sa("CON(" + ws_energy2 + " < 0, 0, " + ws_energy2 + ")", ws_energy)
        gp.SingleOutputMapAlgebra_sa("CON(" +  ws_value1 + " < 0, 0, " + ws_value1 + ")", ws_value)

        gp.AddMessage("\n\tCreated watershed output table: \n\t" + str(hp_ws_table))
        
    except:
        gp.AddError("\nError creating watershed output: " + gp.GetMessages(2))
        raise Exception


    try:

        gp.AddMessage("\nCreating sub-watershed output...")

        # Fraction of total water supply provided by each sub-watershed
        gp.Divide_sa(rsupply_vol, ws_rsupply, subws_rsupply_fract)

        # Energy and value grids per sub-watershed
        gp.Times_sa(subws_rsupply_fract, ws_energy, subws_energy)
        gp.Times_sa(subws_rsupply_fract, ws_value, subws_value)
        # Also get annual energy for output table
        gp.Times_sa(subws_rsupply_fract, ws_energy1, subws_energy_ann)

        # Table of energy/value per sub-watershed, add to scarcity table

        # Zonal stats field name changed in Arc10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = subwshed_id_field
        else:
            zstat_id_field = "Value"

        # Aggregate energy/value by sub-watershed
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, subws_energy_ann, subws_energy_zstat, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, subws_value, subws_value_zstat, "DATA")

        se_rows = gp.SearchCursor(subws_energy_zstat, "", "", "", zstat_id_field + " A")
        se_row = se_rows.Reset
        se_row = se_rows.Next()

        sv_rows = gp.SearchCursor(subws_value_zstat, "", "", "", zstat_id_field + " A")
        sv_row = sv_rows.Reset
        sv_row = sv_rows.Next()

        # Set energy/value outputs in sub-watershed table
        while (se_row):

        # Find corresponding station in output (formerly scarcity) table

            swtable_rows = gp.UpdateCursor(hp_subws_table)
            swtable_row = swtable_rows.Reset
            swtable_row = swtable_rows.Next()

            # Match output table subwatershed to row in energy/value tables
            while (int(swtable_row.getValue(out_table_swid_field)) <> int(se_row.getValue(zstat_id_field))):
                swtable_row = swtable_rows.Next()

            # Can use mean because there's only one value per sub-watershed
            swtable_row.setValue(out_table_energy_field, float(se_row.getValue("MEAN")))
            swtable_row.setValue(out_table_value_field, float(sv_row.getValue("MEAN")))
            swtable_rows.UpdateRow(swtable_row)

            se_row = se_rows.Next()
            sv_row = sv_rows.Next()
            
        gp.AddMessage("\n\tCreated energy per sub-watershed output file: \n\t" + str(subws_energy))   
        gp.AddMessage("\n\tCreated value per sub-watershed output file: \n\t" + str(subws_value))
        gp.AddMessage("\n\tCreated sub-watershed output table:\n\t" + str(hp_subws_table))

    except:
        gp.AddError("\nError creating sub-watershed output: " + gp.GetMessages(2))
        raise Exception


    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(outputws + "\\Hydropower_Valuation_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("HYDROPOWER 3 - VALUATION MODEL PARAMETERS\n")
        parafile.writelines("_________________________________________\n\n")

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
        del wshed_row, wshed_rows, y_row, y_rows
        del o_row, o_rows, se_row, se_rows, sv_row, sv_rows, swtable_row, swtable_rows
        gp.Delete_management(interws)

    except:
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception
except:
    gp.AddError("\nError running script")
    raise Exception
