# ---------------------------------------------------------------------------
# Combine_landuses.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
# 
# Last edit: 4/21/2010 
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

        # Fence grid
        fence = gp.GetParameterAsText(1)
        parameters.append("Fence grid: " + fence)

        # Restoration grid
        restore = gp.GetParameterAsText(2)
        parameters.append("Restoration grid: " + restore)

        # Reforestation grid
        reforest = gp.GetParameterAsText(3)
        parameters.append("Reforestation grid: " + reforest)

        # Silvopastoral grid
        silvo = gp.GetParameterAsText(4)
        parameters.append("Silvopastoral grid: " + silvo)

        # Protection grid (optional)
        protect = gp.GetParameterAsText(5)
        parameters.append("Protection grid: " + protect)

        # Original land use grid
        orig_lu = gp.GetParameterAsText(6)
        parameters.append("Original land use grid: " + orig_lu)

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(7)
        parameters.append("Suffix: " + Suffix)

        if (Suffix == "") or (Suffix == string.whitespace) or (Suffix == "#"):
            Suffix = ""
        else:
            Suffix = "_" + Suffix
            
    except:
        gp.AddError("Error in input arguments: " + gp.GetMessages(2))
        raise Exception


    # Check and create output folders
    try:
        thefolders=["Post_Process", "Intermediate"]
        for folder in thefolders:
            if not gp.Exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        gp.AddError("Error creating folders: " + gp.GetMessages())
        raise Exception


    # Output files
    try:
        # Base output/working directories
        postprocws = gp.workspace + os.sep + "Post_Process" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep
        
        # Intermediate variables
        na_all = interws + "na_all"
        
        # Output layers
        lu_port = postprocws + "port_lu"
        
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [lu_port]           
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

        lu_port = str(Outputnames[0])
            
    except:
        gp.AddError ("Error validating output filenames: " + gp.GetMessages(2))
        raise Exception


    try:

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws

        # Remap input rasters so that NoDatas become zeros
        if not gp.Exists(fence):
            gp.AddMessage("no fence")
            fence = interws + "fence_tmp"
            gp.SingleOutputMapAlgebra_sa("CON(" + orig_lu + " > 0, 0, 0 )", fence)
        if not gp.Exists(restore):
            gp.AddMessage("no restore")
            restore = interws + "rest_tmp"
            gp.SingleOutputMapAlgebra_sa("CON(" + orig_lu + " > 0, 0, 0 )", restore)
        if not gp.Exists(reforest):
            gp.AddMessage("no reforest")
            reforest = interws + "refor_tmp"
            gp.SingleOutputMapAlgebra_sa("CON(" + orig_lu + " > 0, 0, 0 )", reforest)
        if not gp.Exists(silvo):
            gp.AddMessage("no silvo")
            silvo = interws + "silvo_tmp"
            gp.SingleOutputMapAlgebra_sa("CON(" + orig_lu + " > 0, 0, 0 )", silvo)
        if not gp.Exists(protect):
            gp.AddMessage("no protect")
            protect = interws + "prot_tmp"
            gp.SingleOutputMapAlgebra_sa("CON(" + orig_lu + " > 0, 0, 0 )", protect)
        
        # Combine new activities into one file
        gp.SingleOutputMapAlgebra_sa(fence + " + " + restore + " + " + reforest + " + " + silvo + " + " + protect, na_all)
        # Combine new activities with original land use to produce final investment portfolio
        gp.SingleOutputMapAlgebra_sa("CON( " + na_all + " > 0, " + na_all + ", " + orig_lu + ")", lu_port)
        gp.AddMessage("\n\tCreated investment portfolio output file: \n\t" + str(lu_port))
        
    except:
        gp.AddError ("Error creating output files: " + gp.GetMessages(2))
        raise Exception




    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Post_Process\\PP_Combine_Landuses_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("POST-PROCESSING - COMBINE LANDUSES\n")
        parafile.writelines("__________________________________\n\n")

        for para in parameters:
            parafile.writelines(para + "\n")
            parafile.writelines("\n")
        parafile.close()
    except:
        gp.AddError ("Error creating parameter file:" + gp.GetMessages(2))
        raise Exception


    # Clean up temporary files
##    gp.AddMessage("\nCleaning up temporary files...\n")
##    try:
##        gp.Delete_management(interws)
##    except:
##        gp.AddError("Error cleaning up temporary files:  " + gp.GetMessages(2))
##        raise Exception


##    # REMOVE THIS
##    gp.AddMessage("\nIGNORE THE FOLLOWING ERROR\n")
##    raise Exception


except:
    gp.AddError ("Error running script")
    raise Exception
