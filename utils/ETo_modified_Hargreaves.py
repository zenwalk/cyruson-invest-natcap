# ---------------------------------------------------------------------------
# ETo_modified_Hargreaves.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
# 
# Calculates potential ET (ETo) from monthly temperature,
# precipitation and radiation data, using the modified
# Hargreaves method.
#
# Reads temp and precip values from a specified folder,
# calculates monthly ETo values, then combines them into
# a single annual raster.
#
# *** NOTE: Temperature input data is x10 actual values ***
# ***       Temperatures will be divided by 10          ***
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

        # Folder with monthly temperature and precipitation data
        input_dir = gp.GetParameterAsText(1)
        parameters.append("Input folder: " + input_dir)

        # Create dictionary to hold monthly radiation data
        RA = []
        # Fill slot 0 to make subsequent month assignments easier
        RA.append(0)
        
        # Radiation value for January
        RA.append(gp.GetParameterAsText(2))
        parameters.append("Radiation for January: " + str(RA[1]))

        # Radiation value for February
        RA.append(gp.GetParameterAsText(3))
        parameters.append("Radiation for February: " + str(RA[2]))

        # Radiation value for March
        RA.append(gp.GetParameterAsText(4))
        parameters.append("Radiation for March: " + str(RA[3]))

        # Radiation value for April
        RA.append(gp.GetParameterAsText(5))
        parameters.append("Radiation for April: " + str(RA[4]))

        # Radiation value for May
        RA.append(gp.GetParameterAsText(6))
        parameters.append("Radiation for May: " + str(RA[5]))

        # Radiation value for June
        RA.append(gp.GetParameterAsText(7))
        parameters.append("Radiation for June: " + str(RA[6]))

        # Radiation value for July
        RA.append(gp.GetParameterAsText(8))
        parameters.append("Radiation for July: " + str(RA[7]))

        # Radiation value for August
        RA.append(gp.GetParameterAsText(9))
        parameters.append("Radiation for August: " + str(RA[8]))

        # Radiation value for September
        RA.append(gp.GetParameterAsText(10))
        parameters.append("Radiation for September: " + str(RA[9]))

        # Radiation value for October
        RA.append(gp.GetParameterAsText(11))
        parameters.append("Radiation for October: " + str(RA[10]))

        # Radiation value for November
        RA.append(gp.GetParameterAsText(12))
        parameters.append("Radiation for November: " + str(RA[11]))

        # Radiation value for December
        RA.append(gp.GetParameterAsText(13))
        parameters.append("Radiation for December: " + str(RA[12]))

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(14)
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
        thefolders=["Post_Process", "Intermediate"]
        for folder in thefolders:
            if not gp.Exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        gp.AddError("\nError creating folders: " + gp.GetMessages())
        raise Exception


    # Output files
    try:
        # Base output/working directories
        postprocws = gp.workspace + os.sep + "Post_Process" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        
        # Intermediate variables
##        na_all = interws + "na_all"
        
        # Output layers
##        eto_annual = postprocws + "eto" + Suffix + ".tif"
        eto_annual = postprocws + "eto" + Suffix + ".tif"
        
    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    try:

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws

        # Filename prefix for monthly max/min temp and precip
        # Assumption is that the full filenames are
        # tmin_1...tmin_12, tmax_1...tmax_12 prec_1...prec_12
        
        tmin_fname_pre = "tmin_"
        tmax_fname_pre = "tmax_"
        precip_fname_pre = "prec_"

        # Calculate ETo one month at a time
        # range goes to 1 - final value
        for month in range (1, 13):
            try:

                gp.AddMessage("\nCalculating daily ETo for month " + str(month) + "...")
                
                # Create full tmin/max/precip filenames for this month
                tmin_x10 = input_dir + os.sep + tmin_fname_pre + str(month)
                if not gp.exists(tmin_x10):
                    gp.AddError("\nError: Tmin file " + tmin_x10 + " does not exist")
                    raise Exception
                tmin = interws + "tmin_" + str(month)
                gp.Divide_sa(tmin_x10, "10", tmin)
                
                tmax_x10 = input_dir + os.sep + tmax_fname_pre + str(month)
                if not gp.exists(tmax_x10):
                    gp.AddError("\nError: Tmax file " + tmax_x10 + " does not exist")
                    raise Exception
                tmax = interws + "tmax_" + str(month)
                gp.Divide_sa(tmax_x10, "10", tmax)
                
                precip = input_dir + os.sep + precip_fname_pre + str(month)
                if not gp.exists(precip):
                    gp.AddError("\nError: Precip file " + precip + " does not exist")
                    raise Exception

                # Average temp for month
                tavg = interws + "tavg_" + str(month)
                gp.SingleOutputMapAlgebra_sa("(" + tmin + " + " + tmax + ") / 2", tavg)

                # Difference between max and min temps
                tdiff = interws + "tdiff_" + str(month)
                gp.Minus_sa(tmax, tmin, tdiff)

                # ETo calculation for month
                td_p_tmp = interws + "td_p_" + str(month)
                gp.SingleOutputMapAlgebra_sa(tdiff + " - (.0123 * " + precip + ")", td_p_tmp)
                
                eto_day = interws + "eto_d_" + str(month)
                gp.SingleOutputMapAlgebra_sa(".0013 * .408 * " + str(RA[month]) + \
                                             " * (17 + " + tavg + ") * POW(" + td_p_tmp + ", .76)", eto_day)

        
            except:
                gp.AddError ("\nError calculating daily ETo for month " + str(month) + ": " + gp.GetMessages(2))
                raise Exception


        try:

            gp.AddMessage("\nCalculating annual ETo...")

            eto_tmp = interws + "eto_tmp"
            
            for month in range (1, 13):

                eto_month = interws + "eto_m_" + str(month)
                eto_day = interws + "eto_d_" + str(month)

                # 30 days hath September...
                if month == 2:
                    gp.Times_sa(eto_day, "28.2", eto_month)
                elif month in [4, 6, 9, 11]:
                    gp.Times_sa(eto_day, "30", eto_month)
                else:
                    gp.Times_sa(eto_day, "31", eto_month)

                # If first month, running sum just equals that month's value
                if month == 1:
                    gp.CopyRaster_management(eto_month, eto_annual)
                # If subsequent months, add to running sum
                else:
                    gp.Plus_sa(eto_month, eto_annual, eto_tmp)
                    gp.CopyRaster_management(eto_tmp, eto_annual)

            gp.AddMessage("\nCreated annual ETo output: " + eto_annual + "\n")
        
        except:
            gp.AddError ("\nError calculating annual ETo: " + gp.GetMessages(2))
            raise Exception
        
    except:
        gp.AddError ("\nError creating output files" )
        raise Exception




    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Post_Process\\ETo_Modified_Hargreaves_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("ETo MODIFIED HARGREAVES\n")
        parafile.writelines("_______________________\n\n")

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
