# ---------------------------------------------------------------------------
# Sediment_2_Valuation.py
# 
# Coded by :
# Stacie Wolny, Bradley Eichelberger
# Written by:
# Driss Ennaanay, Guillermo Mendoza, Marc Conte 
# for the Natural Capital Project
#
# Last edited: 4/22/2011
#
# Calculates the value of the landscape for keeping
#   sediment out of streams
#
# Used for valuing sediment retention as a function of water quality
#   and/or reservoir dredging
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

        # Watersheds shapefile
        watershed = sys.argv[2]
        parameters.append("Watersheds: " + watershed)

        # Sub-watersheds shapefile
        sub_watersheds = sys.argv[3]
        parameters.append("Sub-watersheds: " + sub_watersheds)

        # Value dredging?
        value_dredging = sys.argv[4]

        if value_dredging =='true':
            value_dredging = True
            parameters.append("Value dredging: Yes")
        else :
            value_dredging = False
            parameters.append("Value dredging: No")

        # Sediment retention for dredging (sum) output from biophysical script
        retention_dr = sys.argv[5]
        parameters.append("Sediment retention for dredging (sum): " + retention_dr)
        
        if value_dredging and not gp.Exists(retention_dr):
            gp.AddError("\nError: If dredging is valued, the sediment retention for dredging (sum) raster must be specified in the parameter window.  Exiting.\n")
            raise Exception

        # Value water quality?
        value_wquality = sys.argv[6]

        if value_wquality =='true':
            value_wquality = True
            parameters.append("Value water quality: Yes")
        else:
            value_wquality = False
            parameters.append("Value water quality: No")

        # Sediment retention for water quality (sum) output from biophysical script
        retention_wq = sys.argv[7]
        parameters.append("Sediment retention for water quality (sum): " + retention_wq)

        if value_wquality and not gp.Exists(retention_wq):
            gp.AddError("\nError: If water quality is valued, the sediment retention for water quality (sum) raster must be specified in the parameter window.  Exiting.\n")
            raise Exception

        if value_wquality == '#' and value_dredging == '#':
            gp.AddError("\nError: You must choose to value Dredging and/or Water Quality in the parameter window checkboxes.  Exiting.\n")
            raise Exception

        # Watershed sediment export table (water quality) from biophysical script
        ws_export_table = sys.argv[8]
        parameters.append("Watershed sediment export table: " + ws_export_table)

        # Sub-watershed sediment export table from biophysical script
        sws_export_table = sys.argv[9]
        parameters.append("Sub-watershed sediment export table: " + sws_export_table)

        # Table with model values for the watershed/point of interest
        sed_table = sys.argv[10]
        parameters.append("Sediment watershed parameter table: " + sed_table)

        # Resolution for output rasters
        resolution = sys.argv[11]
        parameters.append("Output resolution: " + str(resolution))

        # Suffix to add to end of output files, as <filename>_<suffix>
        Suffix = sys.argv[12]
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

        # Final output
        dredging_value = servicews + "sed_val_dr"
        wquality_value = servicews + "sed_val_wq"
        sws_value_table_name = "sediment_value_subwatershed" + Suffix + ".dbf"
        ws_value_table_name = "sediment_value_watershed" + Suffix + ".dbf"

        # Watershed/sub-watershed ID fields
        wshed_id_field = "ws_id"
        subwshed_id_field = "subws_id"

        # output table field names
        out_table_wsid_field = "ws_id"
        out_table_wq_value_field = "sed_val_wq"
        out_table_dr_value_field = "sed_val_dr"
        
    except:
        gp.AddError( "\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Verify length of output rasters does not exceed 13 characters
    try:

        Outputnames = [dredging_value, wquality_value]
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
        dredging_value = str(Outputnames[0])
        wquality_value = str(Outputnames[1])
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
        gp.AddError( "\nError setting geoprocessing environment: " + gp.GetMessages(2))
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
                if (table_field.Name.upper() == wshed_id_field.upper()):
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


    # Create output tables - one for sub-basins, one for basins
    try:

        gp.AddMessage("\nCreating output tables...")

        # Table for whole watershed - add value fields to biophysical output table
        gp.CreateTable_management(servicews, ws_value_table_name)
        ws_value_table = servicews + ws_value_table_name
        gp.CopyRows_management(ws_export_table, ws_value_table)

        if value_dredging:
            gp.AddField(ws_value_table, out_table_dr_value_field, "double")

        if value_wquality:
            gp.AddField(ws_value_table, out_table_wq_value_field, "double")            

        # Table for sub-watersheds - add value fields to biophysical output table
        gp.CreateTable_management(servicews, sws_value_table_name)
        sws_value_table = servicews + sws_value_table_name
        gp.CopyRows_management(sws_export_table, sws_value_table)

        if value_dredging:
            gp.AddField(sws_value_table, out_table_dr_value_field, "double")

        if value_wquality:
            gp.AddField(sws_value_table, out_table_wq_value_field, "double")

    except:
        gp.AddError("\nError creating output tables: " + gp.GetMessages(2))
        raise Exception


    # Value sediment retention for water quality

    if (value_wquality):

        try:

            gp.AddMessage("\nProcessing sediment retention value for water quality...")

            # Intermediate output
            watershed_join = interws + "q_ws_join"
            watershed_join2 = interws + "q_ws_join2"
            watershed_join3 = interws + "q_ws_join3"
            watershed_cost = interws + "q_ws_cost"
            watershed_max_csl = interws + "q_ws_max_csl"
            watershed_pv = "q_wshed_pv"
            retention_value1 = interws + "q_ret_val1"
            wshed_ras = interws + "q_wshed_ras"
            subws_ret_zstat = interws + "subws_ret_zstat_wq.dbf"
            subws_val_zstat =  interws + "subws_val_zstat_wq.dbf"
            watersheds_sjoin = interws + "wsheds_sjoin_wq.shp"

            # Sediment input table field names for valuing water quality
            wq_annual_load_field = "wq_annload"
            wq_cost_field = "wq_cost"
            wq_time_field = "wq_time"
            wq_discount_field = "wq_disc"
            

            # Map table values to watersheds
            try:
                
                # Verify input table fields
                checkfields([wshed_id_field, wq_cost_field, wq_discount_field], sed_table)
            
                gp.AddMessage("\n\tMapping table values to watersheds...")
                
                # Get sediment values from table, map to watershed
                # Copy to dbf to work around table OID requirement change in ArcMap 9.3
                dsc = gp.describe(retention_wq)
                cell_size = str(dsc.MeanCellHeight)

                gp.CopyRows_management(sed_table, "sed_tmp_wq.dbf")
                gp.FeatureToRaster_conversion(watershed, wshed_id_field, wshed_ras, cell_size)
                gp.BuildRasterAttributeTable_management(wshed_ras)
                gp.MakeRasterLayer_management(wshed_ras, "wshed_tmp_wq")
                gp.AddJoin_management("wshed_tmp_wq", "Value", "sed_tmp_wq.dbf", wshed_id_field)
                gp.CopyRaster_management("wshed_tmp_wq", watershed_join)
                gp.Lookup_sa(watershed_join, wq_cost_field, watershed_cost)
                
                if (install_info["Version"] == "10.0" or install_info["Version"] == "10.1"):
                    gp.Delete_management("sed_tmp_wq.dbf")
                    gp.Delete_management("wshed_tmp_wq")


            except:
                gp.AddError("\nError mapping table values to watersheds: " + gp.GetMessages(2))
                raise Exception

            # Create temporary output table to hold discount values
            try:
                
                gp.CreateTable_management(interws, "sed_pv_wq.dbf")
                pv_table = interws + "sed_pv_wq.dbf"
                gp.CopyRows_management(sed_table, pv_table)
                
                # new output table field names
                pv_table_pv_field = "pv"
                gp.AddField(pv_table, pv_table_pv_field, "double")

                pv_rows = gp.UpdateCursor(pv_table)
                pv_row = pv_rows.Reset
                pv_row = pv_rows.Next()

            except:
                gp.AddError("\nError creating output table: " + gp.GetMessages(2))
                raise Exception

            # Calculate value for each watershed
            ws_rows = gp.SearchCursor(watershed)
            ws_row = ws_rows.Reset
            ws_row = ws_rows.Next()

            while ws_row:

                stable_rows = gp.SearchCursor(sed_table)
                stable_row = stable_rows.Reset
                stable_row = stable_rows.Next()

                # Find watershed in sediment valuation table that matches watershed shapefile
                while (int(stable_row.getValue(wshed_id_field)) <> int(ws_row.getValue(wshed_id_field))):
                    stable_row = stable_rows.Next()

                gp.AddMessage("\n\tCalculating present value for watershed id " + str(stable_row.getValue(wshed_id_field)) + "...")

                # Present Value
                try:

                    station_time = float(stable_row.getValue(wq_time_field))
                    station_discount = float(stable_row.getValue(wq_discount_field))

                    pv = 0

                    for t in range (0, int(station_time)):
                        pv += 1 / pow(1 + (station_discount / 100), t)

                except:
                    gp.AddError("\n\tError calculating present value: " + gp.GetMessages())
                    raise Exception

                # Add present value to table
            
                try:

                    while (int(stable_row.getValue(wshed_id_field)) <> int(pv_row.getValue(wshed_id_field))):
                        pv_row = pv_rows.Next()

                    pv_row.setValue(pv_table_pv_field, float(pv))
                    pv_rows.UpdateRow(pv_row)

                    pv_rows = gp.UpdateCursor(pv_table)
                    pv_row = pv_rows.Reset
                    pv_row = pv_rows.Next()
                    
                    stable_row = stable_rows.Next()
                
                except:
                    gp.AddError("\n\tError updating output table: " + gp.GetMessages(2))
                    raise Exception

                ws_row = ws_rows.Next()

            # Calculate value of sediment retention

            try:

                gp.AddMessage ("\n\tCalculating sediment retention value...")

                # Map present values to watershed raster
                gp.CopyRows_management(pv_table, "pv_tmp_wq.dbf")
                gp.MakeRasterLayer_management(wshed_ras, "wshed_tmp_wq3")
                gp.AddJoin_management("wshed_tmp_wq3", "Value", "pv_tmp_wq.dbf", wshed_id_field)
                gp.CopyRaster_management("wshed_tmp_wq3", watershed_join3)
                gp.Lookup_sa(watershed_join3, pv_table_pv_field, watershed_pv)

                if (install_info["Version"] == "10.0" or install_info["Version"] == "10.1"):
                    gp.Delete_management("pv_tmp_wq.dbf")
                    gp.Delete_management("wshed_tmp_wq3")

                # Value of landscape for retaining nutrient
                gp.SingleOutputMapAlgebra_sa(watershed_pv +  " * " + retention_wq + " * " + watershed_cost, retention_value1)

                # Can't have negative value, so change negatives to zero
                gp.SingleOutputMapAlgebra_sa("CON( " + retention_value1 + " < 0 , 0, " + retention_value1 + " )", wquality_value )

                gp.AddMessage("\n\tCreated sediment retention value for water quality output file: \n\t" + str(wquality_value))

                del stable_row, stable_rows, ws_row, ws_rows, pv_row, pv_rows

            except:
                gp.AddError ("\nError calculating sediment retention value: " + gp.GetMessages(2))
                raise Exception

            # Add value to output tables

            try:

                gp.AddMessage("\n\tUpdating sub-watershed output table for water quality...")

                # Sub-watershed table
                
                # Aggregate value to sub-watersheds
                gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, wquality_value, subws_val_zstat, "DATA")

                # Zonal stats field name changed in Arc10
                if (install_info["Version"] == "10.0" or install_info["Version"] == "10.1"):
                    zstat_id_field = subwshed_id_field
                else:
                    zstat_id_field = "Value"

                
                # Add subwatershed values to table

                swsval_rows = gp.UpdateCursor(sws_value_table)
                swsval_row = swsval_rows.Reset
                swsval_row = swsval_rows.Next()

                while(swsval_row):

                    valz_rows = gp.SearchCursor(subws_val_zstat, "", "", "", zstat_id_field + " A")
                    valz_row = valz_rows.Reset
                    valz_row = valz_rows.Next()

                    # Match watershed in output table with a row in the value table
                    while (int(swsval_row.getValue(subwshed_id_field)) <> int(valz_row.getValue(zstat_id_field))):
                        valz_row = valz_rows.Next()


                    # Can use mean because there's only one value per sub-watershed
                    swsval_row.setValue(out_table_wq_value_field, float(valz_row.getValue("MEAN")))
                    swsval_rows.UpdateRow(swsval_row)
                    swsval_row = swsval_rows.Next()
                    
                del valz_row, valz_rows

                gp.AddMessage("\n\tUpdated sub-watershed output table for water quality: \n\t" + str(sws_value_table))
            
                
                gp.AddMessage("\n\tUpdating watershed output table for water quality...")

                # Watershed table

                # Find all sub-watersheds within this watershed and sum their value to create watershed output

                # Map sub-watersheds to their corresponding watersheds
                gp.SpatialJoin_analysis(sub_watersheds, watershed, watersheds_sjoin, "JOIN_ONE_TO_ONE", "KEEP_ALL", "#", "IS_WITHIN")

                wsval_rows = gp.UpdateCursor(ws_value_table)
                wsval_row = wsval_rows.Reset
                wsval_row = wsval_rows.Next()

                while (wsval_row):

                    # Order watersheds and value by subwatershed id
                    sj_rows = gp.SearchCursor(watersheds_sjoin, "", "", "", subwshed_id_field + " A")
                    sj_row = sj_rows.Reset
                    sj_row = sj_rows.Next()

                    valz_rows = gp.SearchCursor(subws_val_zstat, "", "", "", zstat_id_field + " A")
                    valz_row = valz_rows.Reset
                    valz_row = valz_rows.Next()
                
                    val_sum = 0

                    while (sj_row):

                        if (int(sj_row.getValue(wshed_id_field)) == int(wsval_row.getValue(wshed_id_field))):
                            val_sum += valz_row.getValue("MEAN")

                        sj_row = sj_rows.Next()
                        valz_row = valz_rows.Next()

                    # Fill in table values for this watershed
                    wsval_row.setValue(out_table_wq_value_field, float(val_sum))
                    wsval_rows.UpdateRow(wsval_row)

                    wsval_row = wsval_rows.Next()

                del sj_row, sj_rows, valz_row, valz_rows
                
                gp.AddMessage("\n\tUpdated watershed output table for water quality: \n\t" + str(ws_value_table))                


            except:
                gp.AddError ("\nError updating output tables for water quality: " + gp.GetMessages(2))
                raise Exception

        except:
            gp.AddError ("\nError processing sediment value for water quality: " + gp.GetMessages(2))
            raise Exception


    # Value sediment retention for reservoir dredging

    if (value_dredging):

        try:
            gp.AddMessage("\nProcessing sediment retention value for dredging...")

            # Intermediate output
            wshed_ras = interws + "d_ws_ras"
            watershed_join = interws + "d_ws_join"
            watershed_join2 = interws + "d_ws_join2"
            watershed_cost = interws + "d_ws_cost"
            watershed_time = interws + "d_ws_time"
            watershed_pv = interws + "d_ws_pv"
            wshed_num_cells = interws + "d_num_cells"
            sed_ret = interws + "d_sed_ret"
            sed_ret0 = interws + "d_sed_ret0"
            max_sed = interws + "d_max_sed"
            subws_ret_zstat = interws + "subws_ret_zstat_dr.dbf"
            subws_val_zstat =  interws + "subws_val_zstat_dr.dbf"
            watersheds_sjoin = interws + "wsheds_sjoin_dr.shp"

            # Sediment input table field names for valuing dredging
            dredge_dead_vol_field = "dr_deadvol"
            dredge_cost_field = "dr_cost"
            dredge_time_field = "dr_time"
            dredge_discount_field = "dr_disc"

     
            # Create temporary output table to hold discount values
            try:

                # Verify input table fields
                checkfields([wshed_id_field, dredge_cost_field, dredge_time_field, dredge_discount_field], sed_table)
                
                gp.CreateTable_management(interws, "sed_pv_dr.dbf")
                pv_table = interws + "sed_pv_dr.dbf"
                gp.CopyRows_management(sed_table, pv_table)
                
                # new output table field names
                pv_table_pv_field = "pv"
                gp.AddField(pv_table, pv_table_pv_field, "double")

                pv_rows = gp.UpdateCursor(pv_table)
                pv_row = pv_rows.Reset
                pv_row = pv_rows.Next()

            except:
                gp.AddError("\nError creating temporary output table: " + gp.GetMessages(2))
                raise Exception

                
            try:

                # Calculate present value

                # Get sediment values from table, map to watershed           
                # Copy to dbf to work around table OID requirement change in ArcMap 9.3
                gp.CopyRows_management(sed_table, "sed_tmp_dr.dbf")

                dsc = gp.describe(retention_dr)
                cell_size = str(dsc.MeanCellHeight)

                gp.FeatureToRaster_conversion(watershed, wshed_id_field, wshed_ras, cell_size)
                gp.BuildRasterAttributeTable_management(wshed_ras)
                gp.MakeRasterLayer_management(wshed_ras, "wshed_tmp_dr")
                gp.AddJoin_management("wshed_tmp_dr", "Value", "sed_tmp_dr.dbf", wshed_id_field)
                gp.CopyRaster_management("wshed_tmp_dr", watershed_join)
                gp.Lookup_sa(watershed_join, dredge_time_field, watershed_time)
                gp.Lookup_sa(watershed_join, dredge_cost_field, watershed_cost)

                if (install_info["Version"] == "10.0" or install_info["Version"] == "10.1"):
                    gp.Delete_management("sed_tmp_dr.dbf")
                    gp.Delete_management("wshed_tmp_dr")

                # Find matching watershed and table values
                ws_rows = gp.SearchCursor(watershed)
                ws_row = ws_rows.Reset
                ws_row = ws_rows.Next()

                while ws_row:
                    
                    stable_rows = gp.SearchCursor(sed_table)
                    stable_row = stable_rows.Reset
                    stable_row = stable_rows.Next()

                    while (int(stable_row.getValue(wshed_id_field)) <> int(ws_row.getValue(wshed_id_field))):
                        stable_row = stable_rows.Next()

                    gp.AddMessage("\n\tCalculating present value for reservoir id " + str(stable_row.getValue(wshed_id_field)) + "...")
                    
                    # Net Present Value calculation
                    try:

                        station_time = float(stable_row.getValue(dredge_time_field))
                        station_discount = float(stable_row.getValue(dredge_discount_field))
                        station_id = int(stable_row.getValue(wshed_id_field))

                        pv = 0

                        for t in range (0, int(station_time)):
                            pv += 1 / pow(1 + (station_discount / 100), t)
                                    
                    except:
                        gp.AddError("\n\tError calculating present value: " + gp.GetMessages(2))
                        raise Exception

                    
                    # Add new field values to discount table
                    
                    try:

                        while (int(stable_row.getValue(wshed_id_field)) <> int(pv_row.getValue(wshed_id_field))):
                            pv_row = pv_rows.Next()
                        
                        pv_row.setValue(pv_table_pv_field, float(pv))
                        pv_rows.UpdateRow(pv_row)

                        pv_row = pv_rows.Next()
                        stable_row = stable_rows.Next()
                    
                    except:
                        gp.AddError("\n\tError updating output table: " + gp.GetMessages(2))
                        raise Exception

                    ws_row = ws_rows.Next()

            except:
                gp.AddError("\n\tError calculating present value: " + gp.GetMessages(2))
                raise Exception

            # Calculate value of sediment retention for dredging
            gp.AddMessage ("\n\tCalculating sediment retention value...")
            try:

                # Map discount values to watershed raster
                gp.CopyRows_management(pv_table, "pv_tmp_dr.dbf")
                gp.MakeRasterLayer_management(wshed_ras, "wshed_tmp_dr2")
                gp.AddJoin_management("wshed_tmp_dr2", "Value", "pv_tmp_dr.dbf", wshed_id_field)
                gp.CopyRaster_management("wshed_tmp_dr2", watershed_join2)
                gp.Lookup_sa(watershed_join2, pv_table_pv_field, watershed_pv)

                if (install_info["Version"] == "10.0" or install_info["Version"] == "10.1"):
                    gp.Delete_management("pv_tmp_dr.dbf")
                    gp.Delete_management("wshed_tmp_dr2")               

                # Value of landscape for retaining sediment
                gp.SingleOutputMapAlgebra_sa(watershed_pv + " * " + retention_dr + " * " + watershed_cost, dredging_value)

                gp.AddMessage("\n\tCreated sediment retention value for dredging output file: \n\t" + str(dredging_value))

                del stable_row, stable_rows, ws_row, ws_rows, pv_row, pv_rows

            except:
                gp.AddError("\nError calculating sediment retention value:" + gp.GetMessages(2))
                raise Exception

            # Add value to output tables

            try:

                gp.AddMessage("\n\tUpdating sub-watershed output table for dredging...")

                # Sub-watershed table
                
                # Aggregate value to sub-watersheds
                gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, dredging_value, subws_val_zstat, "DATA")

                # Zonal stats field name changed in Arc10
                if (install_info["Version"] == "10.0" or install_info["Version"] == "10.1"):
                    zstat_id_field = subwshed_id_field
                else:
                    zstat_id_field = "Value"

                
                # Add subwatershed values to table

                swsval_rows = gp.UpdateCursor(sws_value_table)
                swsval_row = swsval_rows.Reset
                swsval_row = swsval_rows.Next()
                
                while(swsval_row):

                    valz_rows = gp.SearchCursor(subws_val_zstat, "", "", "", zstat_id_field + " A")
                    valz_row = valz_rows.Reset
                    valz_row = valz_rows.Next()

                    while (int(swsval_row.getValue(subwshed_id_field)) <> int(valz_row.getValue(zstat_id_field))):
                        valz_row = valz_rows.Next()

                    # Can use mean because there's only one value per sub-watershed
                    swsval_row.setValue(out_table_dr_value_field, float(valz_row.getValue("MEAN")))
                    swsval_rows.UpdateRow(swsval_row)
                    swsval_row = swsval_rows.Next()

                    del valz_row, valz_rows
                
                gp.AddMessage("\n\tUpdated sub-watershed output table for dredging: \n\t" + str(sws_value_table))
            
                
                gp.AddMessage("\n\tUpdating watershed output table for dredging...")

                # Watershed table

                # Find all sub-watersheds within this watershed and sum their retention/value to create watershed output

                # Map sub-watersheds to their corresponding watersheds
                gp.SpatialJoin_analysis(sub_watersheds, watershed, watersheds_sjoin, "JOIN_ONE_TO_ONE", "KEEP_ALL", "#", "IS_WITHIN")

                wsval_rows = gp.UpdateCursor(ws_value_table)
                wsval_row = wsval_rows.Reset
                wsval_row = wsval_rows.Next()

                while (wsval_row):

                    # Order watersheds and value by subwatershed id
                    sj_rows = gp.SearchCursor(watersheds_sjoin, "", "", "", subwshed_id_field + " A")
                    sj_row = sj_rows.Reset
                    sj_row = sj_rows.Next()

                    valz_rows = gp.SearchCursor(subws_val_zstat, "", "", "", zstat_id_field + " A")
                    valz_row = valz_rows.Reset
                    valz_row = valz_rows.Next()
                
                    val_sum = 0

                    while (sj_row):
                        if (int(sj_row.getValue(wshed_id_field)) == int(wsval_row.getValue(wshed_id_field))):
                            val_sum += valz_row.getValue("MEAN")

                        sj_row = sj_rows.Next()
                        valz_row = valz_rows.Next()

                    # Fill in table values for this watershed
                    wsval_row.setValue(out_table_dr_value_field, float(val_sum))
                    wsval_rows.UpdateRow(wsval_row)

                    wsval_row = wsval_rows.Next()

                    del sj_row, sj_rows, valz_row, valz_rows
                
                gp.AddMessage("\n\tUpdated watershed output table for dredging: \n\t" + str(ws_value_table))                


            except:
                gp.AddError ("\nError creating output tables for dredging: " + gp.GetMessages(2))
                raise Exception
    
        except:
            gp.AddError ("\nError processing sediment retention value for dredging: " + gp.GetMessages(2))
            raise Exception


    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Output\\Sediment_Valuation_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("SEDIMENT 2 - VALUATION MODEL PARAMETERS\n")
        parafile.writelines("_______________________________________\n\n")

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
        del wsval_row, wsval_rows, swsval_row, swsval_rows
        gp.Delete_management(interws)
    except:
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception

    
except:
    gp.AddError ("\nError running script.")
    raise Exception
