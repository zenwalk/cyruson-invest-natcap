# ---------------------------------------------------------------------------
# Hydropower_Water_Yield.py
#
# Coded by :
# Stacie Wolny, Bradley Eichelberger
# Written by:
# Driss Ennaanay, Guillermo Mendoza, Marc Conte
# for the Natural Capital Project
#
# Last edit: 7/19/2011
#
# Calculates actual evapotranspiration and surface water yield for each
# sub-basin in a landscape
# ---------------------------------------------------------------------------

# Import system modules
import sys, string, os, arcgisscripting, time, datetime

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out necessary licenses
gp.CheckOutExtension("spatial")

# Allow overwriting of output files
gp.overwriteoutput = 1

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
        
        # Folder where output files will be saved
        gp.workspace = gp.GetParameterAsText(0)
        parameters.append("Workspace: " + gp.workspace)
        
        # Precipitation raster
        precip = gp.GetParameterAsText(1)
        parameters.append("Precipitation: " + precip)

        # Evapotranspiration raster
        eto = gp.GetParameterAsText(2)
        parameters.append("Evapotranspiration: " + eto)

        # Soil depth raster 
        soil_depth = gp.GetParameterAsText(3)
        parameters.append("Soil depth: " + soil_depth)

        # Plant available water fraction raster
        pawf = gp.GetParameterAsText(4)
        parameters.append("Available water content: " + pawf)

        # Landuse raster
        landuse = gp.GetParameterAsText(5)
        parameters.append("Landcover: " + landuse)

        # Watersheds shapefile
        watersheds = gp.GetParameterAsText(6)
        parameters.append("Watersheds: " + watersheds)

        # Sub-watersheds shapefile
        sub_watersheds = gp.GetParameterAsText(7)
        parameters.append("Sub-watersheds: " + sub_watersheds)

        # Biophysical coefficient table
        lu_table = gp.GetParameterAsText(8)
        parameters.append("Biophysical coefficient table: " + lu_table)

        # Zhang constant
        zhang = gp.GetParameterAsText(9)
        parameters.append("Zhang constant: " + zhang)

        # Output resolution
        resolution = gp.GetParameterAsText(10)
        parameters.append("Output resolution: " + str(resolution))

        # Results suffix
        Suffix = gp.GetParameterAsText(11)
        parameters.append("Suffix: " + Suffix)

        if (Suffix == "") or (Suffix == string.whitespace) or (Suffix == "#"):
            Suffix = ""
        else:
            Suffix = "_" + Suffix
        
    except:
        gp.AddMessage("\nError in input arguments: " + gp.GetMessages(2))
        raise Exception


    # Check and create output folders

    try:
        thefolders=["Intermediate", "Output", "Service"]
        for folder in thefolders:
            if not gp.exists(gp.workspace + folder):
                gp.CreateFolder_management(gp.workspace, folder)
        pixelfolder = "Pixel"
        pfparent = gp.workspace + os.sep + "Output"
        if not gp.exists(pfparent + folder):
            gp.CreateFolder_management(pfparent, pixelfolder)

    except:
        gp.AddError("\nError creating output folders: " + gp.GetMessages(2))
        raise Exception


    # Output files
    
    try:
        
        # Intermediate and output directories
        outputws = gp.workspace + os.sep + "Output" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        servicews = gp.workspace + os.sep + "Service" + os.sep
        pixelws = gp.workspace + os.sep + "Output" + os.sep + "Pixel" + os.sep

        # Temporary variables
        tmp_DI = interws + "tmp_DI"
        tmp_eff_depth = interws + "tmp_eff_depth"
        tmp_AWC = interws + "tmp_AWC"
        tmp_Max_AET = interws + "tmp_Max_AET"
        tmp_flt1000et = interws + "tmp_flt1000et"
        tmp_PET = interws + "tmp_PET"
        tmp_Root_Depth = interws + "tmp_RootDep"
        tmp_wodivP = interws + "tmp_wodivP"
        tmp_W = interws + "tmp_W"
        tmp_root_dpth = interws + "tmp_root_dpth"
        tmp_v1000ETk = interws + "tmp_v1000ETk"
        sub_watersheds_areas = interws + "subws_areas.shp"
        subws_areas = interws + "subws_areas"
        ws_precip_table = interws + "ws_precip_zstat.dbf"
        ws_pet_table = interws + "ws_pet_zstat.dbf"
        ws_aet_table = interws + "ws_aet_zstat.dbf"
        ws_wyield_table = interws + "ws_wyield_zstat.dbf"
        sws_precip_table = interws + "sws_precip_zstat.dbf"
        sws_pet_table = interws + "sws_pet_zstat.dbf"
        sws_aet_table = interws + "sws_aet_zstat.dbf"
        sws_wyield_table = interws + "sws_wyield_zstat.dbf"
        watersheds_sjoin = interws + "wsheds_sjoin.shp"
        

        # Input Biophysical Model table field names
        lucode_field = "lucode"
        root_depth_field = "root_depth"
        etk_field = "etk"

        # ID fields for watersheds/sub-watersheds inputs
        wshed_id_field = "ws_id"
        subwshed_id_field = "subws_id"
        
        # Output files
        fractp = pixelws + "fractp"
        aet = pixelws + "aet"
        wyield = pixelws + "wyield"
        fractp_mean = outputws + "fractp_mn"
        aet_mean = outputws + "aet_mn"
        wyield_mean = servicews + "wyield_mn"
        wyield_vol = servicews + "wyield_vol"
        wyield_ha = servicews + "wyield_ha"
        ws_out_table_name = "water_yield_watershed" + Suffix + ".dbf"
        sws_out_table_name = "water_yield_subwatershed" + Suffix + ".dbf"

    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [aet_mean, wyield_vol, wyield_mean, fractp_mean, wyield, aet, fractp, wyield_ha]
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
                aet_mean = str(Outputnames[0])
                wyield_vol = str(Outputnames[1])
                wyield_mean = str(Outputnames[2])
                fractp_mean = str(Outputnames[3])
                wyield = str(Outputnames[4])
                aet = str(Outputnames[5])
                fractp = str(Outputnames[6])
                wyield_ha = str(Outputnames[7])
    except:
        gp.AddError ("\nError validating output filenames: " + gp.GetMessages(2))
        raise Exception
    

    # Check input raster projections - they should all be the same

    try:
        gp.AddMessage("\nChecking input raster projections...")
        precipDesc = gp.describe(landuse)
        precipspatref = precipDesc.SpatialReference
        rasters = (eto, soil_depth, pawf)
        for x in rasters:
            rasterDesc = gp.describe(x)
            spatreflc = rasterDesc.SpatialReference
            if spatreflc.Type <> 'Projected':
                gp.AddMessage(x + " does not appear to be projected.  It is assumed to be in meters")
            elif spatreflc.LinearUnitName <> 'Meter':
                gp.AddMessage("This model assumes that data in " + x + " is projected in meters.  You may get erroneous results")
            if str(precipspatref.name) <> str(spatreflc.name):
                gp.AddError("\nError: " + x + " is not in the same coordinate system as the precipitation raster.  " + x + " is projected in " + spatreflc.name + " while the precipitation layer is in " + precipspatref.name + ".  Please project all rasters in the same projection and rerun this tool.  Exiting tool....")  
                raise Exception
    except:
        gp.AddError("\nError checking input raster projections: " + gp.GetMessages(2)) 
        raise Exception


    # Set the Geoprocessing environment

    try:
        gp.cellSize = resolution
        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
        install_info = gp.GetInstallInfo("desktop")

    except:
        gp.AddError("\nError setting geoprocessing environment: " + gp.GetMessages(2))
        raise Exception


    # Verify required input table fields' existence and type

    def checkfields (fields, table):
        
        try:

            table_fields = gp.listfields(table, "*", "All")
            table_field = table_fields.next()

            foundfields = []
            while (table_field):
                # Allow for upper case field names, especially for .dbf input
                foundfields.append(table_field.Name.upper())
                if (table_field.Name.upper() == lucode_field.upper()) or (table_field.Name.upper() == etk_field.upper()) or (table_field.Name.upper() == root_depth_field.upper()):
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
            

    # Run yield model

    try:

        # Verify input table fields
        checkfields([lucode_field, etk_field, root_depth_field], lu_table)            

        # Set the watersheds input layer to be the mask for subsequent calculations
        gp.mask = watersheds

        gp.AddMessage("\nRunning Yield calculations...")

        # AET/Yield calculations based on the Budyko curve
        
        gp.ReclassByTable_sa(landuse, lu_table, lucode_field, lucode_field, etk_field, tmp_v1000ETk, "DATA")

        gp.Float_sa(tmp_v1000ETk, tmp_flt1000et)

        gp.SingleOutputMapAlgebra_sa(eto + " * " + tmp_flt1000et + " DIV 1000", tmp_PET)

        gp.Divide_sa(tmp_PET, precip, tmp_DI)

        gp.ReclassByTable_sa(landuse, lu_table, lucode_field, lucode_field, root_depth_field, tmp_root_dpth, "DATA")

        gp.Float_sa(tmp_root_dpth, tmp_Root_Depth)
        
        gp.SingleOutputMapAlgebra_sa("MIN( " + soil_depth + "," + tmp_Root_Depth + " )", tmp_eff_depth)

        gp.Times_sa(tmp_eff_depth, pawf, tmp_AWC)

        gp.SingleOutputMapAlgebra_sa(tmp_AWC + " DIV (" + precip + " + 1)", tmp_wodivP)

        gp.SingleOutputMapAlgebra_sa(zhang + " * " + tmp_wodivP, tmp_W)

        gp.SingleOutputMapAlgebra_sa("CON( " + tmp_DI + " > 1, 1, " + tmp_DI + " )", tmp_Max_AET)

        # The fraction of Precipitation that leaves the system as Evapotranspiration 
        gp.SingleOutputMapAlgebra_sa("MIN( " + tmp_Max_AET + ", ( " + tmp_W + " * " + tmp_DI + " + 1) DIV ( ( 1 DIV " + tmp_DI + " ) + " + tmp_W + " * " + tmp_DI + " + 1 ) )", fractp)
        # Aggregate results into sub-watersheds
        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, fractp, fractp_mean, "MEAN", "DATA")
        
        # Surface water yield
        gp.SingleOutputMapAlgebra_sa("( 1 - " + fractp + " ) * " + precip, wyield)
        # Aggregate results into sub-watersheds
        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, wyield, wyield_mean, "MEAN", "DATA")

        # Water yield volume
        gp.CalculateAreas_stats(sub_watersheds, sub_watersheds_areas)
        cell_size = gp.GetRasterProperties_management(wyield, "CELLSIZEX")
        gp.FeatureToRaster_conversion(sub_watersheds_areas, "F_AREA", subws_areas, cell_size)
        gp.SingleOutputMapAlgebra_sa(wyield_mean + " * " + subws_areas + " / 1000", wyield_vol)

        # Water yield per hectare output - 1 meter = .0001 hectare
        gp.SingleOutputMapAlgebra_sa("(" + wyield_mean + " * " + subws_areas + " / 1000) / (.0001 * " + subws_areas + ")", wyield_ha)
        
        # Actual evapotranspiration
        gp.SingleOutputMapAlgebra_sa(fractp + " * " + precip, aet)
        # Aggregate results into sub-watersheds
        gp.ZonalStatistics_sa(sub_watersheds, subwshed_id_field, aet, aet_mean, "MEAN", "DATA")

        gp.AddMessage("\n\tCreated Precipitation Fraction output file:\n\t" + str(fractp_mean))
        gp.AddMessage("\n\tCreated Actual Evapotranspiration output file:\n\t" + str(aet_mean))
        gp.AddMessage("\n\tCreated Water Yield output files:\n\t" + str(wyield_mean) + "\n\t" + str(wyield_vol) + "\n\t" + str(wyield_ha))

    except:
        gp.AddError("\nError running Yield calculations:" + gp.GetMessages(2))
        raise Exception

    # Table to hold generated output values per watershed
    gp.AddMessage("\nCreating watershed yield output table...")
    try:
        # output table field names
        out_table_wshed_id_field = "ws_id"
        out_table_precip_field = "precip_mn"
        out_table_pet_field = "PET_mn"
        out_table_aet_field = "AET_mn"
        out_table_wyield_field = "wyield_mn"
        out_table_wyield_sum_field = "wyield_sum"

        gp.CreateTable_management(outputws, ws_out_table_name)
        ws_out_table = outputws + ws_out_table_name

        gp.AddField(ws_out_table, out_table_wshed_id_field, "long")
        gp.AddField(ws_out_table, out_table_precip_field, "double")
        gp.AddField(ws_out_table, out_table_pet_field, "double")
        gp.AddField(ws_out_table, out_table_aet_field, "double")
        gp.AddField(ws_out_table, out_table_wyield_field, "double")
        gp.AddField(ws_out_table, out_table_wyield_sum_field, "double")
        
        # Remove Field1 - it's added by default and not used
        gp.DeleteField_management(ws_out_table, "Field1")

        out_table_rows = gp.InsertCursor(ws_out_table)

        # Aggregate values by watershed
        gp.ZonalStatisticsAsTable_sa(watersheds, wshed_id_field, precip, ws_precip_table, "DATA")
        gp.ZonalStatisticsAsTable_sa(watersheds, wshed_id_field, eto, ws_pet_table, "DATA")
        gp.ZonalStatisticsAsTable_sa(watersheds, wshed_id_field, aet, ws_aet_table, "DATA")
        gp.ZonalStatisticsAsTable_sa(watersheds, wshed_id_field, wyield, ws_wyield_table, "DATA")

        # Zonal stats field name changed in Arc10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = wshed_id_field
        else:
            zstat_id_field = "Value"

        precip_rows = gp.SearchCursor(ws_precip_table, "", "", "", zstat_id_field + " A")
        precip_row = precip_rows.Reset
        precip_row = precip_rows.Next()

        pet_rows = gp.SearchCursor(ws_pet_table, "", "", "", zstat_id_field + " A")
        pet_row = pet_rows.Reset
        pet_row = pet_rows.Next()

        aet_rows = gp.SearchCursor(ws_aet_table, "", "", "", zstat_id_field + " A")
        aet_row = aet_rows.Reset
        aet_row = aet_rows.Next()

        wyield_rows = gp.SearchCursor(ws_wyield_table, "", "", "", zstat_id_field + " A")
        wyield_row = wyield_rows.Reset
        wyield_row = wyield_rows.Next()

        # Add values to table
        while (precip_row):
            # Zonal stats field name has changed in Arc 10
            if (install_info["Version"] == "10.0"):
                ws_id = precip_row.getValue(wshed_id_field)
            else:
                ws_id = precip_row.getValue("VALUE")
                
            ws_precip = float(precip_row.getValue("MEAN"))
            ws_pet = float(pet_row.getValue("MEAN"))
            ws_aet = float(aet_row.getValue("MEAN"))
            ws_wyield = float(wyield_row.getValue("MEAN"))
            ws_wyield_sum = float(wyield_row.getValue("SUM"))
            
            new_row = out_table_rows.NewRow()
            new_row.setValue(out_table_wshed_id_field, ws_id)
            new_row.setValue(out_table_precip_field, ws_precip)
            new_row.setValue(out_table_pet_field, ws_pet)
            new_row.setValue(out_table_aet_field, ws_aet)
            new_row.setValue(out_table_wyield_field, ws_wyield)
            new_row.setValue(out_table_wyield_sum_field, ws_wyield_sum)

            out_table_rows.InsertRow(new_row)
            precip_row = precip_rows.Next()
            pet_row = pet_rows.Next()
            aet_row = aet_rows.Next()
            wyield_row = wyield_rows.Next()

        del precip_row, precip_rows, pet_row, pet_rows, aet_row, aet_rows, wyield_row, wyield_rows
        del out_table_rows, new_row

        gp.AddMessage("\n\tCreated watershed water yield output table: \n\t" + str(ws_out_table))

    except:
        gp.AddError ("\nError creating watershed water yield output table:" + gp.GetMessages(2))
        raise Exception


    # Table to hold generated output values per sub-watershed
    gp.AddMessage("\nCreating sub-watershed yield output table...")
    try:
        # output table field names
        out_table_wshed_id_field = "ws_id"
        out_table_subwshed_id_field = "subws_id"
        out_table_precip_field = "precip_mn"
        out_table_pet_field = "PET_mn"
        out_table_aet_field = "AET_mn"
        out_table_wyield_field = "wyield_mn"
        out_table_wyield_sum_field = "wyield_sum"

        gp.CreateTable_management(outputws, sws_out_table_name)
        sws_out_table = outputws + sws_out_table_name

        gp.AddField(sws_out_table, out_table_wshed_id_field, "long")
        gp.AddField(sws_out_table, out_table_subwshed_id_field, "long")
        gp.AddField(sws_out_table, out_table_precip_field, "double")
        gp.AddField(sws_out_table, out_table_pet_field, "double")
        gp.AddField(sws_out_table, out_table_aet_field, "double")
        gp.AddField(sws_out_table, out_table_wyield_field, "double")
        gp.AddField(sws_out_table, out_table_wyield_sum_field, "double")
        
        # Remove Field1 - it's added by default and not used
        gp.DeleteField_management(sws_out_table, "Field1")

        out_table_rows = gp.InsertCursor(sws_out_table)

        # Map sub-watersheds to their corresponding watersheds
        gp.SpatialJoin_analysis(sub_watersheds, watersheds, watersheds_sjoin, "JOIN_ONE_TO_ONE", "KEEP_ALL", "#", "IS_WITHIN")

        # Aggregate values by sub-watershed
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, precip, sws_precip_table, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, eto, sws_pet_table, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, aet, sws_aet_table, "DATA")
        gp.ZonalStatisticsAsTable_sa(sub_watersheds, subwshed_id_field, wyield, sws_wyield_table, "DATA")

        # Zonal stats field name changed in Arc10
        if (install_info["Version"] == "10.0"):
            zstat_id_field = subwshed_id_field
        else:
            zstat_id_field = "Value"

        precip_rows = gp.SearchCursor(sws_precip_table, "", "", "", zstat_id_field + " A")
        precip_row = precip_rows.Reset
        precip_row = precip_rows.Next()

        pet_rows = gp.SearchCursor(sws_pet_table, "", "", "", zstat_id_field + " A")
        pet_row = pet_rows.Reset
        pet_row = pet_rows.Next()

        aet_rows = gp.SearchCursor(sws_aet_table, "", "", "", zstat_id_field + " A")
        aet_row = aet_rows.Reset
        aet_row = aet_rows.Next()

        wyield_rows = gp.SearchCursor(sws_wyield_table, "", "", "", zstat_id_field + " A")
        wyield_row = wyield_rows.Reset
        wyield_row = wyield_rows.Next()

        # Add values to table
        while (precip_row):

            sj_rows = gp.SearchCursor(watersheds_sjoin, "", "", "", subwshed_id_field + " A")
            sj_row = sj_rows.Reset
            sj_row = sj_rows.Next()

            # Match subwatershed in spatial join to one in the precip table
            while (int(sj_row.getValue(subwshed_id_field)) <> int(precip_row.getValue(zstat_id_field))):
                sj_row = sj_rows.Next()           

            
            # Zonal stats field name has changed in Arc 10
            if (install_info["Version"] == "10.0"):
                sws_id = precip_row.getValue(subwshed_id_field)
            else:
                sws_id = precip_row.getValue("VALUE")
                
            sws_precip = float(precip_row.getValue("MEAN"))
            sws_pet = float(pet_row.getValue("MEAN"))
            sws_aet = float(aet_row.getValue("MEAN"))
            sws_wyield = float(wyield_row.getValue("MEAN"))
            sws_wyield_sum = float(wyield_row.getValue("SUM"))
            
            new_row = out_table_rows.NewRow()
            new_row.setValue(out_table_wshed_id_field, int(sj_row.getValue(wshed_id_field)))
            new_row.setValue(out_table_subwshed_id_field, sws_id)
            new_row.setValue(out_table_precip_field, sws_precip)
            new_row.setValue(out_table_pet_field, sws_pet)
            new_row.setValue(out_table_aet_field, sws_aet)
            new_row.setValue(out_table_wyield_field, sws_wyield)
            new_row.setValue(out_table_wyield_sum_field, sws_wyield_sum)

            out_table_rows.InsertRow(new_row)
            precip_row = precip_rows.Next()
            pet_row = pet_rows.Next()
            aet_row = aet_rows.Next()
            wyield_row = wyield_rows.Next()

            del sj_row, sj_rows

        del precip_row, precip_rows, pet_row, pet_rows, aet_row, aet_rows, wyield_row, wyield_rows
        del out_table_rows, new_row

        gp.AddMessage("\n\tCreated sub-watershed water yield output table: \n\t" + str(sws_out_table))

    except:
        gp.AddError ("\nError creating sub-watershed water yield output table:" + gp.GetMessages(2))
        raise Exception
    

    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        parafile = open(outputws + "\\Water_Yield" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("WATER YIELD MODEL PARAMETERS\n")
        parafile.writelines("____________________________\n\n")

        for para in parameters:
            parafile.writelines(para + "\n")
            parafile.writelines("\n")
        parafile.close()
    except:
        gp.AddError ("\nError creating parameter file: " + gp.GetMessages(2))
        raise Exception
    
        
    # Clean up temporary files
    gp.AddMessage("\nCleaning up temporary files...\n")
    try:
        gp.Delete_management(interws)
    except:
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception


except:
    gp.AddError("\nError running script")  
    raise Exception
