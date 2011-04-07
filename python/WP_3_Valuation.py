# ---------------------------------------------------------------------------
# WP_3_Valuation.py
# 
# Coded by :
# Stacie Wolny, Bradley Eichelberger
# Written by:
# Driss Ennaanay, Guillermo Mendoza, Marc Conte 
# for the Natural Capital Project
#
# Last edited: 2/6/2011
#
# Second script for Water Purification.
# Calculates the value of the landscape for keeping
#   nutrient pollution out of water 
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

    try:

        # Make parameters array, and later write input parameter values to an output file
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))

        # Script arguments
        gp.AddMessage ("\nValidating arguments..." )

        # Folder where intermediate and output files are written
        gp.workspace = sys.argv[1]
        parameters.append("Workspace: " + gp.workspace)

        # Shapefile of the watershed contributing to the point of interest
        watershed = sys.argv[2]
        parameters.append("Watershed: " + watershed)

        # Shapefile of sub-watersheds
        sub_watersheds = sys.argv[3]
        parameters.append("Sub-watersheds: " + sub_watersheds)

        # Nutrient retention (sum) output from biophysical script
        retention = sys.argv[4]
        parameters.append("Nutrient retention (sum): " + retention)

        # Watershed nutrient export table from biophysical script
        ws_export_table = sys.argv[5]
        parameters.append("Watershed nutrient export table: " + ws_export_table)

        # Sub-watershed nutrient export table from biophysical script
        sws_export_table = sys.argv[6]
        parameters.append("Sub-watershed nutrient export table: " + sws_export_table)

        # Table with model values for the watershed/point of interest
        wp_table = sys.argv[7]
        parameters.append("Water purification table: " + wp_table)

        # Resolution for output rasters
        resolution = sys.argv[8]
        parameters.append("Output resolution: " + str(resolution))

        # Suffix to add to end of output files, as <filename>_<suffix>
        Suffix = sys.argv[9]
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
        gp.AddError( "\nError creating folders: " + gp.GetMessages(2))
        raise Exception


    # Output files
    try:

        # Output/working directories
        outputws = gp.workspace + os.sep + "Output" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        servicews = gp.workspace + os.sep + "Service" + os.sep

        # Intermediate output
        watershed_join = interws + "ws_join"
        watershed_join2 = interws + "ws_join2"
        watershed_join3 = interws + "q_ws_join3"
        watershed_cost = interws + "ws_cost"
##        watershed_annual_load = interws + "ws_annld"
        watershed_max_cnl = interws + "ws_max_cnl"
        retention_value1 = interws + "ret_val1"
        wshed_ras = interws + "wshed_ras"
        watershed_pv = "q_wshed_pv"
        watersheds_sjoin = "wsheds_sjoin.shp"
        wsheds_sjoin_copy = "ws_sjoin_copy.shp"
        subws_val_zstat = "subws_val_zstat.dbf"
        subws_ret_zstat = "subws_ret_zstat.dbf"

        # Water Purification input table field names
        wp_id_field = "ws_id"
        calib_field = "calib"
        annual_load_field = "ann_load"
        cost_field = "cost"
        time_field = "time_span"
        discount_field = "discount"

        retention_value = servicews + "nut_val"
        nut_sws_value_table_name = "nutrient_value_subwatershed" + Suffix + ".dbf"
        nut_ws_value_table_name = "nutrient_value_watershed" + Suffix + ".dbf"

        # ID fields for watersheds/sub-watersheds
        wshed_id_field = "ws_id"
        subwshed_id_field = "subws_id"
        
    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Verify length of output rasters does not exceed 13 characters
    try:
        Outputnames = [retention_value]
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
        retention_value = str(Outputnames[0])
    except:
        gp.AddError ("\nError validating output filenames: " + gp.GetMessages(2))
        raise Exception


    # Set the Geoprocessing environment
    try:
        gp.cellSize = resolution
        gp.mask = watershed
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
                # Make sure that fields that need to be integer are
                if (table_field.Name.upper() == wp_id_field.upper()):
                    if not "Integer" in table_field.Type:
                        gp.AddError("\nError: Field " + str(table_field.Name) + " in table " + str(table) + " must be of type Integer\n")
                        raise Exception
                table_field = table_fields.next()

            for f in fields:
                if not f.upper() in foundfields:
                    gp.AddError("\nError: Required table field " + str(f) + " not found in input table " + str(table) + "\n")
                    raise Exception

        except:
            gp.AddError ("\nError verifying input table fields")
            raise Exception


    # Calculate present value
    
    try:

        # Verify input table fields
        checkfields([wp_id_field, calib_field, cost_field], wp_table)

        gp.AddMessage("\nCalculating present value...")
        
        # Get water quality values from table, map to watershed
        # Copy to dbf to work around table OID requirement change in ArcMap 9.3
        dsc = gp.describe(retention)
        cell_size = str(dsc.MeanCellHeight)
        
        gp.CopyRows_management(wp_table, "wp_tmp.dbf")
        gp.FeatureToRaster_conversion(watershed, wshed_id_field, wshed_ras, cell_size)
        gp.MakeRasterLayer_management(wshed_ras, "wsheds_tmp3")
        gp.AddJoin_management("wsheds_tmp3", "Value", "wp_tmp.dbf", wp_id_field)
        gp.CopyRaster_management("wsheds_tmp3", watershed_join)
##        gp.Lookup_sa(watershed_join, annual_load_field, watershed_annual_load)
        gp.Lookup_sa(watershed_join, cost_field, watershed_cost)
        # Delete now or Arc won't delete the Intermediate folder later
        gp.Delete_management("wsheds_tmp3")

        # Create temporary output table to hold discount values
        try:
            
            gp.CreateTable_management(interws, "wq_pv.dbf")
            pv_table = interws + "wq_pv.dbf"
            gp.CopyRows_management(wp_table, pv_table)
            
            # new output table field names
            pv_table_pv_field = "pv"
            gp.AddField(pv_table, pv_table_pv_field, "double")

            pv_rows = gp.UpdateCursor(pv_table)
            pv_row = pv_rows.Reset
            pv_row = pv_rows.Next()

        except:
            gp.AddError("\nError creating output table: " + gp.GetMessages(2))
            raise Exception
            
        ws_rows = gp.SearchCursor(watershed)
        ws_row = ws_rows.Reset
        ws_row = ws_rows.Next()

        while ws_row:

            wtable_rows = gp.SearchCursor(wp_table)
            wtable_row = wtable_rows.Reset
            wtable_row = wtable_rows.Next()

            while (int(wtable_row.getValue(wp_id_field)) <> int(ws_row.getValue(wshed_id_field))):
                wtable_row = wtable_rows.Next()

            gp.AddMessage("\n\tCalculating present value for watershed id " + str(wtable_row.getValue(wp_id_field)))

            # Present Value
            try:

                station_time = float(wtable_row.getValue(time_field))
                station_discount = float(wtable_row.getValue(discount_field))
                station_cost = float(wtable_row.getValue(cost_field))

                pv = 0

                for t in range (0, station_time):
                    pv += 1 / pow(1 + (station_discount / 100), t)
                            
            except:
                gp.AddError("\nError calculating Present Value: " + gp.GetMessages())
                raise Exception

            # Add new field values to discount table
            
            try:

                while (int(wtable_row.getValue(wp_id_field)) <> int(pv_row.getValue(wp_id_field))):
                    pv_row = pv_rows.Next()

                pv_row.setValue(pv_table_pv_field, float(pv))
                pv_rows.UpdateRow(pv_row)

                pv_row = pv_rows.Next()
                wtable_row = wtable_rows.Next()
            
            except:
                gp.AddError("\n\tError updating output table: " + gp.GetMessages(2))
                raise Exception

            ws_row = ws_rows.Next()

    except:
        gp.AddError("\nError calculating present value:" + gp.GetMessages(2))
        raise Exception

    try:
        gp.AddMessage("\nCalculating nutrient retention value...")

        # Map present values to watershed raster
        gp.CopyRows_management(pv_table, "pv_tmp_wq.dbf")
        gp.MakeRasterLayer_management(wshed_ras, "wshed_tmp_wq3")
        gp.AddJoin_management("wshed_tmp_wq3", "Value", "pv_tmp_wq.dbf", wp_id_field)
        gp.CopyRaster_management("wshed_tmp_wq3", watershed_join3)
        gp.Lookup_sa(watershed_join3, pv_table_pv_field, watershed_pv)
        # Delete now or Arc won't delete the Intermediate folder later
        gp.Delete_management("wshed_tmp_wq3")
        gp.Delete_management(watershed_join3)

        # Value of landscape for retaining pollution
        gp.SingleOutputMapAlgebra_sa(watershed_pv + " * " + watershed_cost + " * " + retention, retention_value1)

        # Can't have negative value, so change negatives to zero
        gp.SingleOutputMapAlgebra_sa("CON( " + retention_value1 + " < 0 , 0, " + retention_value1 + " )", retention_value)

        gp.AddMessage("\n\tCreated nutrient retention value output file: \n\t" + str(retention_value))

    except:
        gp.AddError("\nError calculating nutrient retention value:" + gp.GetMessages(2))
        raise Exception

    # Populate output tables with retention/value
    
    try:
        
        # Create output tables - one for sub-basins, one for basins

        out_table_retention_field = "nut_ret"
        out_table_value_field = "nut_value"

        gp.CreateTable_management(servicews, nut_sws_value_table_name)
        nut_sws_value_table = servicews + nut_sws_value_table_name
        gp.CopyRows_management(sws_export_table, nut_sws_value_table)

        gp.AddField(nut_sws_value_table, out_table_retention_field, "double")
        gp.AddField(nut_sws_value_table, out_table_value_field, "double")

        gp.CreateTable_management(servicews, nut_ws_value_table_name)
        nut_ws_value_table = servicews + nut_ws_value_table_name
        gp.CopyRows_management(ws_export_table, nut_ws_value_table)

        gp.AddField(nut_ws_value_table, out_table_retention_field, "double")
        gp.AddField(nut_ws_value_table, out_table_value_field, "double")

        gp.AddMessage("\nCreating sub-watershed output table...")

        # Sub-watershed table
        
        # Make tables of retention/value to map to sub-watersheds
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, retention, subws_ret_zstat, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, retention_value, subws_val_zstat, "DATA")

        # Zonal stats field name changed in Arc10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = subwshed_id_field
        else:
            zstat_id_field = "Value"

        
        # Add subwatershed values to table

        swsval_rows = gp.UpdateCursor(nut_sws_value_table)
        swsval_row = swsval_rows.Reset
        swsval_row = swsval_rows.Next()

        while(swsval_row):

            retz_rows = gp.SearchCursor(subws_ret_zstat, "", "", "", zstat_id_field + " A")
            retz_row = retz_rows.Reset
            retz_row = retz_rows.Next()

            valz_rows = gp.SearchCursor(subws_val_zstat, "", "", "", zstat_id_field + " A")
            valz_row = valz_rows.Reset
            valz_row = valz_rows.Next()

            while (int(swsval_row.getValue(subwshed_id_field)) <> int(retz_row.getValue(zstat_id_field))):
                retz_row = retz_rows.Next()
                valz_row = valz_rows.Next()

            # Can use mean because there's only one value per sub-watershed

            swsval_row.setValue(out_table_retention_field, float(retz_row.getValue("MEAN")))
            swsval_row.setValue(out_table_value_field, float(valz_row.getValue("MEAN")))
            swsval_rows.UpdateRow(swsval_row)
            swsval_row = swsval_rows.Next()

            del retz_row, retz_rows, valz_row, valz_rows
        
        gp.AddMessage("\n\tCreated sub-watershed output table: \n\t" + str(nut_sws_value_table))
    
        
        gp.AddMessage("\nCreating watershed output table...")

        # Watershed table

        # Find all sub-watersheds within this watershed and sum their retention/value to create watershed output

        # Map sub-watersheds to their corresponding watersheds
        gp.SpatialJoin_analysis(sub_watersheds, watershed, watersheds_sjoin, "JOIN_ONE_TO_ONE", "KEEP_ALL", "#", "IS_WITHIN")

        wsval_rows = gp.UpdateCursor(nut_ws_value_table)
        wsval_row = wsval_rows.Reset
        wsval_row = wsval_rows.Next()

        while (wsval_row):

            # Order watersheds, retention and value by subwatershed id
            sj_rows = gp.SearchCursor(watersheds_sjoin, "", "", "", subwshed_id_field + " A")
            sj_row = sj_rows.Reset
            sj_row = sj_rows.Next()

            retz_rows = gp.SearchCursor(subws_ret_zstat, "", "", "", zstat_id_field + " A")
            retz_row = retz_rows.Reset
            retz_row = retz_rows.Next()

            valz_rows = gp.SearchCursor(subws_val_zstat, "", "", "", zstat_id_field + " A")
            valz_row = valz_rows.Reset
            valz_row = valz_rows.Next()
        
            ret_sum = 0
            val_sum = 0

            while (sj_row):

                if (int(sj_row.getValue(wshed_id_field)) == int(wsval_row.getValue(wshed_id_field))):
                    ret_sum += retz_row.getValue("MEAN")
                    val_sum += valz_row.getValue("MEAN")

                sj_row = sj_rows.Next()
                retz_row = retz_rows.Next()
                valz_row = valz_rows.Next()

            # Fill in table values for this watershed
            wsval_row.setValue(out_table_retention_field, float(ret_sum))
            wsval_row.setValue(out_table_value_field, float(val_sum))
            wsval_rows.UpdateRow(wsval_row)

            wsval_row = wsval_rows.Next()
            
            del sj_rows, sj_row, retz_row, retz_rows, valz_row, valz_rows
        
        gp.AddMessage("\n\tCreated watershed output table: \n\t" + str(nut_ws_value_table))
        
    except:
        gp.AddError("\nError creating output tables: " + gp.GetMessages(2))
        raise Exception
    

    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Output\\Water_Purification_Valuation_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("WATER PURIFICATION 3 - VALUATION MODEL PARAMETERS\n")
        parafile.writelines("_________________________________________________\n\n")

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
        del ws_rows, ws_row, wtable_rows, wtable_row, pv_rows, pv_row, wsval_row, wsval_rows
        gp.Delete_management(interws)
    except:
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception

    
except:
    gp.AddError ("\nError running script.")
    raise Exception
