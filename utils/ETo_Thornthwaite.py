# ---------------------------------------------------------------------------
# ETo_mThornthwaite.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
# 
# Calculates potential ET (ETo) from mean temperature and
# hours of sun information using the Thornthwaite equation
#
# Reads temp values from a specified folder,
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
        HS = []
        # Fill slot 0 to make subsequent month assignments easier
        HS.append(0)
        
        # Hours of sun value for January
        HS.append(gp.GetParameterAsText(2))
        parameters.append("Hours of sun for January: " + str(HS[1]))

        # Hours of sun value for February
        HS.append(gp.GetParameterAsText(3))
        parameters.append("Hours of sun for February: " + str(HS[2]))

        # Hours of sun value for March
        HS.append(gp.GetParameterAsText(4))
        parameters.append("Hours of sun for March: " + str(HS[3]))

        # Hours of sun value for April
        HS.append(gp.GetParameterAsText(5))
        parameters.append("Hours of sun for April: " + str(HS[4]))

        # Hours of sun value for May
        HS.append(gp.GetParameterAsText(6))
        parameters.append("Hours of sun for May: " + str(HS[5]))

        # Hours of sun value for June
        HS.append(gp.GetParameterAsText(7))
        parameters.append("Hours of sun for June: " + str(HS[6]))

        # Hours of sun value for July
        HS.append(gp.GetParameterAsText(8))
        parameters.append("Hours of sun for July: " + str(HS[7]))

        # Hours of sun value for August
        HS.append(gp.GetParameterAsText(9))
        parameters.append("Hours of sun for August: " + str(HS[8]))

        # Hours of sun value for September
        HS.append(gp.GetParameterAsText(10))
        parameters.append("Hours of sun for September: " + str(HS[9]))

        # Hours of sun value for October
        HS.append(gp.GetParameterAsText(11))
        parameters.append("Hours of sun for October: " + str(HS[10]))

        # Hours of sun value for November
        HS.append(gp.GetParameterAsText(12))
        parameters.append("Hours of sun for November: " + str(HS[11]))

        # Hours of sun value for December
        HS.append(gp.GetParameterAsText(13))
        parameters.append("Hours of sun for December: " + str(HS[12]))

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
        eto_annual = postprocws + "eto" + Suffix
        
    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception


    try:

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws

        # Filename prefix for mean temp inputs
        # Assumption is that the full filenames are from WorldClim, zone 23 (Cali)
        # tmean1_23.tif...tmin12_23.tif
        
        tmean_fname_pre = "tmean"
        tmean_fname_suf = "_23.tif"

        # Calculate ETo one month at a time
        # range needs to go to 13 so that 12 is calculated
        for month in range (1, 13):
            try:

                gp.AddMessage("\nCalculating daily ET for month " + str(month) + "...")
                
                # Create full tmean filename for this month
                tmean_x10 = input_dir + os.sep + tmean_fname_pre + str(month) + tmean_fname_suf
                if not gp.exists(tmean_x10):
                    gp.AddError("\nError: Tmean file " + tmean_x10 + " does not exist")
                    raise Exception
                tmean = interws + "tmean_" + str(month)
                gp.Divide_sa(tmean_x10, "10", tmean)

                # Thornthwaite says that if temp < 0, set it to 0
                tmean_noneg = interws + "tmean0_" + str(month)
                gp.SingleOutputMapAlgebra_sa("CON(" + tmean + " < 0, 0, " + tmean + ")", tmean_noneg)

                # Daily ETo calculation for month (i/5) ^ 1.514
                t_tmp = interws + "t_" + str(month)
                gp.Divide_sa(tmean_noneg, "5", t_tmp)
                
                eto_day = interws + "eto_d_" + str(month)
                gp.SingleOutputMapAlgebra_sa("POW(" + t_tmp + ", 1.514)", eto_day)
        
            except:
                gp.AddError ("\nError calculating daily ET for month " + str(month) + ": " + gp.GetMessages(2))
                raise Exception


        try:

            gp.AddMessage("\nCalculating heat index I...")

            # Heat index I - sum of monthly values calculated above (I = sum((i/5)^1.514)
            heat_index = interws + "hidx"
            heat_index_tmp = interws + "hidx_tmp"
            
            for month in range (1, 13):
                # Value calculated above
                eto_day = interws + "eto_d_" + str(month)
                
                # If first month, running sum just equals that month's value
                if month == 1:
                    gp.CopyRaster_management(eto_day, heat_index)
                # If subsequent months, add to running sum
                else:
                    gp.Plus_sa(eto_day, heat_index, heat_index_tmp)
                    gp.CopyRaster_management(heat_index_tmp, heat_index)

        except:
                gp.AddError ("\nError calculating heat index I: " + gp.GetMessages(2))
                raise Exception


        try:

            # Monthly ET = 16 * (day length/12) * (month days/30)
            #                   * (((10 * mean temp)/I)^alpha)
            # alpha = (6.75x10^-7 * I^3) - (7.71x10^-5 * I^2) + (1.792x10^-2 * I) + 0.49239
            # I = heat index calculated above
            # if mean temp < 0, set it to 0


            alpha = interws + "alpha"

            gp.SingleOutputMapAlgebra_sa("(.000000675 * POW(" + heat_index + ", 3)) - " + \
                                         "(.0000771 * POW(" + heat_index + ", 2)) + " + \
                                         "(.01792 * " + heat_index + ") + 0.49239)", alpha)
                                                     
            for month in range (1, 13):

                gp.AddMessage("\nCalculating ET for month " + str(month) + "...")

                eto_month = interws + "eto_m_" + str(month)
                eto_day = interws + "eto_d_" + str(month)
                tmean_noneg = interws + "tmean0_" + str(month)
                alpha_arg = "alpha_arg"

                # 30 days hath September...
                if month == 2:
                    month_days = "28"
                elif month in [4, 6, 9, 11]:
                    month_days = "30"
                else:
                    month_days = "31"

                gp.SingleOutputMapAlgebra_sa("(10 * " + tmean_noneg + ") / " + heat_index, alpha_arg)
                gp.SingleOutputMapAlgebra_sa("16 * (" + HS[month] + " / 12) * (" + \
                                             month_days + " / 30) * POW(" + alpha_arg + ", " + alpha \
                                             + ")", eto_month)
        except:
            gp.AddError ("\nError calculating ET for month " + str(month) + ": " + gp.GetMessages(2))
            raise Exception
                                                             

        try:

            # Add all months up to get final annual PET value
            
            gp.AddMessage("\nCalculating annual ET...")
            
            for month in range (1, 13):

                eto_month = interws + "eto_m_" + str(month)
                eto_tmp = interws + "eto_tmp"

                # If first month, running sum just equals that month's value
                if month == 1:
                    gp.CopyRaster_management(eto_month, eto_annual)
                # If subsequent months, add to running sum
                else:
                    gp.Plus_sa(eto_month, eto_annual, eto_tmp)
                    gp.CopyRaster_management(eto_tmp, eto_annual)

            gp.AddMessage("\nCreated annual ET output: " + eto_annual + "\n")
        
        except:
            gp.AddError ("\nError calculating annual ET: " + gp.GetMessages(2))
            raise Exception
        
    except:
        gp.AddError ("\nError creating output files" )
        raise Exception




    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Post_Process\\ETo_Thorntwaithe_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("ETo Thornthwaithe\n")
        parafile.writelines("_________________\n\n")

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
