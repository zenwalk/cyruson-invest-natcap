#---------------------------------------------------------------------
#
# Hydropower_Scarcity.py
#
# Coded by :
# Stacie Wolny, Bradley Eichelberger
# Written by:
# Driss Ennaanay, Guillermo Mendoza, Marc Conte
# for the Natural Capital Project
#
# Last edit: 2/6/2011
#
# Creates grid of realized supply - water supply minus demand -
# over the landscape, and a table of water supply and demand
# per watershed.
#
# Takes into account only surface water supply
#
#---------------------------------------------------------------------

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
        
        # Water supply volume raster from Yield script
        wyield_vol = gp.GetParameterAsText(1)
        parameters.append("Water yield (volume): " + wyield_vol)

        # Water supply mean raster from Yield script
        wyield_mean = gp.GetParameterAsText(2)
        parameters.append("Water yield (mean): " + wyield_mean)

        # Land use raster
        landuse = gp.GetParameterAsText(3)
        parameters.append("Landcover: " + landuse)

        # Watersheds shapefile
        watersheds = gp.GetParameterAsText(4)
        parameters.append("Watersheds: " + watersheds)

        # Sub-watersheds shapefile
        sub_watersheds = gp.GetParameterAsText(5)
        parameters.append("Sub-watersheds: " + sub_watersheds)

        # Water demand table
        demand_table = gp.GetParameterAsText(6)
        parameters.append("Water demand table: " + demand_table)

        # Hydropower station info table 
        station_table = gp.GetParameterAsText(7)
        parameters.append("Hydropower table: " + station_table)

        # Output resolution
        resolution = gp.GetParameterAsText(8)
        parameters.append("Output resolution: " + str(resolution))
        
        # Results suffx
        Suffix = gp.GetParameterAsText(9)
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
        thefolders=["Output","Intermediate"]
        for folder in thefolders:
            if not gp.exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
                
    except:
        gp.AddError("\nError creating folders: " + gp.GetMessages(2))
        raise Exception


    # Output files/directories

    try:
        # Intermediate and output directories
        outputws = gp.workspace + os.sep + "Output" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep

        # Temporary variables
        consump_rc = interws + "consump_rc"
        consump = interws + "consump"
        consump2 = interws + "consump2"
        consump_zstat = interws + "cons_zstat"
        supply_zstat = interws + "sup_zstat"
        watersheds_join = interws + "wsheds_join"
        watersheds_calib = interws + "wsheds_calib"
        wshed_ras = interws + "wshed_ras"
        wyield_calib_zstat = interws + "wyield_calib.dbf"
        wyield_mean_zstat = interws + "wyield_mean.dbf"
        consump_vol_zstat = interws + "consum_vol.dbf"
        consump_mean_zstat = interws + "consum_mean.dbf"
        rsupply_vol_zstat = interws + "rsupply_vol.dbf"
        rsupply_mean_zstat = interws + "rsupply_mean.dbf"

        # Input table field names
        lucode_field = "lucode"
        demand_field = "demand"
        # Hydropower table field names
        station_id_field = "ws_id"
        station_calib_field = "calib"
        # ID fields for watersheds/sub-watersheds inputs
        wshed_id_field = "ws_id"
        subwshed_id_field = "subws_id"

        # Output files
        wyield_calib = outputws + "cyield_vol"
        consump_vol = outputws + "consum_vol"
        consump_mean = outputws + "consum_mn"
        rsupply_vol = outputws + "rsupply_vol"
        rsupply_mean = outputws + "rsupply_mn"
        out_table_name = "water_scarcity"

    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Verify length of output rasters does not exceed 13 characters
    try:
        Outputnames = [consump_vol, consump_mean, rsupply_vol, rsupply_mean, wyield_calib]
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
        consump_vol =  str(Outputnames[0])
        consump_mean = str(Outputnames[1])
        rsupply_vol = str(Outputnames[2])
        rsupply_mean = str(Outputnames[3])
        wyield_calib = str(Outputnames[4])
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
        gp.AddError( "Error configuring output resolution: " + gp.GetMessages(2))
        raise Exception
    

    # Verify required input table fields' existence and type

    def checkfields (fields, table):
        
        try:

            table_fields = gp.listfields(table, "*", "All")
            table_field = table_fields.next()

            foundfields = []
            while (table_field):
                foundfields.append(table_field.Name.upper())
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
            gp.AddError ("\nError verifying input table fields")
            raise Exception


    # Create realized supply grid (yield minus consumption)

    try:

        # Verify input table fields
        checkfields ([station_id_field, station_calib_field], station_table)
        checkfields ([lucode_field, demand_field], demand_table)
        
        # Set the watersheds input layer to be the mask for subsequent calculations
        gp.mask = watersheds

        gp.AddMessage("\nProcessing water supply...")

        # Map calibration constants to each watershed
        # Copy to dbf to work around table OID requirement change in ArcMap 9.3
        gp.CopyRows_management(station_table, "stable_tmp.dbf")

        desc = gp.describe(wyield_vol)
        cell_size = str(desc.MeanCellHeight)

        gp.FeatureToRaster_conversion(watersheds, wshed_id_field, wshed_ras, cell_size)
        gp.MakeRasterLayer_management(wshed_ras, "wsheds_tmp")
        gp.AddJoin_management("wsheds_tmp", "Value", "stable_tmp.dbf", station_id_field)
        gp.CopyRaster_management("wsheds_tmp", watersheds_join)
        gp.Lookup_sa(watersheds_join, station_calib_field, watersheds_calib)
        # Delete now or Arc won't delete the Intermediate folder later
        gp.Delete_management("wsheds_tmp")
        gp.Delete_management("stable_tmp.dbf")

        # Supply times the calibration constant
        gp.Times_sa(watersheds_calib, wyield_vol, wyield_calib)

        gp.AddMessage("\nProcessing water demand...")
        
        # Reclass by demand per land use class
        gp.ReclassByTable_sa(landuse, demand_table, lucode_field, lucode_field, demand_field, consump_rc, "NODATA")
        gp.SingleOutputMapAlgebra_sa("CON(IsNull(" + consump_rc + "), 0, " + consump_rc + ")", consump)
        gp.SingleOutputMapAlgebra_sa("FLOAT(" + consump + ")", consump2)
        # Aggregate by sub-watershed
        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, consump2, consump_vol, "SUM")
        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, consump2, consump_mean, "MEAN")

        # Yield minus consumption
        gp.Minus_sa(wyield_calib, consump_vol, rsupply_vol)
        gp.Minus_sa(wyield_mean, consump_mean, rsupply_mean)

        gp.AddMessage("\nCreated calibrated yield output grid: \n\t" + str(wyield_calib))
        gp.AddMessage("\nCreated consumption output grids: \n\t" + str(consump_vol) + "\n\t" + str(consump_mean))
        gp.AddMessage("\nCreated realized supply output grids: \n\t" + str(rsupply_vol) + "\n\t" + str(rsupply_mean))

    except:
        gp.AddError("\nError calculating water scarcity: " + gp.GetMessages(2))
        raise Exception

    # Create and populate table with output values for supply and consumption
    try:
        gp.AddMessage("\nCreating output table...")

        # output table field names
        out_table_swsid_field = "subws_id"
        out_table_cyield_vol_field = "cyield_vl"
        out_table_yield_mean_field = "yield_mn"
        out_table_consump_vol_field = "consump_vl"
        out_table_consump_mean_field = "consump_mn"
        out_table_supply_vol_field = "rsupply_vl"
        out_table_supply_mean_field = "rsupply_mn"

        out_table_dbf = out_table_name + Suffix + ".dbf"
        gp.CreateTable_management(outputws, out_table_dbf)
        out_table = outputws + out_table_dbf
        
        gp.AddField(out_table, out_table_swsid_field, "long")
        gp.AddField(out_table, out_table_cyield_vol_field, "double")
        gp.AddField(out_table, out_table_yield_mean_field, "double")
        gp.AddField(out_table, out_table_consump_vol_field, "double")
        gp.AddField(out_table, out_table_consump_mean_field, "double")
        gp.AddField(out_table, out_table_yield_mean_field, "double")
        gp.AddField(out_table, out_table_supply_vol_field, "double")
        gp.AddField(out_table, out_table_supply_mean_field, "double")

        gp.DeleteField_management(out_table, "Field1")

        out_table_rows = gp.InsertCursor(out_table)

        # Supply and consumption totals by sub-watershed
        
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, wyield_calib, wyield_calib_zstat, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, wyield_mean, wyield_mean_zstat, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, consump_vol, consump_vol_zstat, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, consump_mean, consump_mean_zstat, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, rsupply_vol, rsupply_vol_zstat, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, rsupply_mean, rsupply_mean_zstat, "DATA")

        # Zonal stats field name changed in Arc10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = subwshed_id_field
        else:
            zstat_id_field = "Value"

        # Set up cursors for each zonal stat table, sort to ensure matching sub-watersheds
        
        yc_rows = gp.SearchCursor(wyield_calib_zstat, "", "", "", zstat_id_field + " A")
        yc_row = yc_rows.Reset
        yc_row = yc_rows.Next()

        ym_rows = gp.SearchCursor(wyield_mean_zstat, "", "", "", zstat_id_field + " A")
        ym_row = ym_rows.Reset
        ym_row = ym_rows.Next()

        cv_rows = gp.SearchCursor(consump_vol_zstat, "", "", "", zstat_id_field + " A")
        cv_row = cv_rows.Reset
        cv_row = cv_rows.Next()

        cm_rows = gp.SearchCursor(consump_mean_zstat, "", "", "", zstat_id_field + " A")
        cm_row = cm_rows.Reset
        cm_row = cm_rows.Next()

        sv_rows = gp.SearchCursor(rsupply_vol_zstat, "", "", "", zstat_id_field + " A")
        sv_row = sv_rows.Reset
        sv_row = sv_rows.Next()

        sm_rows = gp.SearchCursor(rsupply_mean_zstat, "", "", "", zstat_id_field + " A")
        sm_row = sm_rows.Reset
        sm_row = sm_rows.Next()

        sws_rows = gp.SearchCursor(sub_watersheds)
        sws_row = sws_rows.Reset
        sws_row = sws_rows.Next()

        # Populate output table with supply/consumption values
        while (yc_row):

            subwshed_id = int(yc_row.getValue(zstat_id_field))

            # Since there's only one value across each sub-watershed,
            # the mean will be that value

            out_table_wyield_calib = float(yc_row.getValue("MEAN"))
            out_table_wyield_mean = float(ym_row.getValue("MEAN"))
            out_table_consump_vol = float(cv_row.getValue("MEAN"))
            out_table_consump_mean = float(cm_row.getValue("MEAN"))
            out_table_rsupply_vol = float(sv_row.getValue("MEAN"))
            out_table_rsupply_mean = float(sm_row.getValue("MEAN"))

            new_row = out_table_rows.NewRow()
            new_row.setValue(out_table_swsid_field, subwshed_id)
        
            new_row.setValue(out_table_cyield_vol_field, out_table_wyield_calib)
            new_row.setValue(out_table_yield_mean_field, out_table_wyield_mean)
            new_row.setValue(out_table_consump_vol_field, out_table_consump_vol)
            new_row.setValue(out_table_consump_mean_field, out_table_consump_mean)
            new_row.setValue(out_table_supply_vol_field, out_table_rsupply_vol)
            new_row.setValue(out_table_supply_mean_field, out_table_rsupply_mean)
            
            out_table_rows.InsertRow(new_row)

            yc_row = yc_rows.Next()
            ym_row = ym_rows.Next()
            cv_row = cv_rows.Next()
            cm_row = cm_rows.Next()
            sv_row = sv_rows.Next()
            sm_row = sm_rows.Next()
                
        gp.AddMessage("\nCreated output table: \n\t" + str(out_table))

    except:
        gp.AddError("\nError creating output table: " + gp.GetMessages(2))
        raise Exception
    

    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(outputws + "\\Hydropower_Scarcity_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("HYDROPOWER 2 - SCARCITY MODEL PARAMETERS\n")
        parafile.writelines("________________________________________\n\n")

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
        del yc_row, yc_rows, ym_row, ym_rows, cv_row, cv_rows, cm_row, cm_rows, sv_row, sv_rows, sm_row, sm_rows
        gp.Delete_management(interws)
    except:
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception


except Exception, ErrorDesc:
    gp.AddError("\nError running script")  
    raise Exception

