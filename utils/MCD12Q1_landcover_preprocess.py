#-----------------------------------------------------------------
#
# MCD12Q1_landcover_preprocess.py
#
# Coded by :
# Stacie Wolny
# for the Natural Capital Project
#
# Last edit: 2/28/2011
#
# Pre-processes the MCD12Q1 Modis yearly global landcover data,
# which is provided in HDF format.  Extracts only the IGBP
# landcover layer and mosaics the tiles together.
#
#-----------------------------------------------------------------

# Import system modules
import sys, string, os, arcgisscripting, math, time, datetime, glob

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

        # Directory where input HDF tiles are located
        hdf_dir = gp.GetParameterAsText(1)
        parameters.append("HDF source directory: " + hdf_dir)

        # Raster with reference projection info
        ref_proj = gp.GetParameterAsText(2)
        parameters.append("Reference projection file: " + ref_proj)

    except:
        gp.AddError("\nError in specifying arguments: " + gp.GetMessages(2))
        raise Exception


    # Check and create folders

    try:
        output_dir_name = "Extracted_grids_IGBP_landcover"
        
        # Temporary working dir
        if not gp.exists(gp.workspace + "Intermediate"):
            gp.CreateFolder_management(gp.workspace, "Intermediate")
        # Final grid output dir
        if not gp.exists(hdf_dir + output_dir_name):
            gp.CreateFolder_management(hdf_dir, output_dir_name)
                
    except:
        gp.AddError("\nError creating folders: " + gp.GetMessages(2))
        raise Exception


    # Temporary files

    # Base output/working directories
    outputws = hdf_dir + os.sep + output_dir_name + os.sep
    interws = gp.workspace + os.sep + "Intermediate" + os.sep

    # Intermediate variables

    # Output filenames
    mosaic_2001 = "MCD12_IGBP_01"
    mosaic_2002 = "MCD12_IGBP_02"
    mosaic_2003 = "MCD12_IGBP_03"
    mosaic_2004 = "MCD12_IGBP_04"
    mosaic_2005 = "MCD12_IGBP_05"

    # Set the Geoprocessing environment 

    try:
        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws

    except:
        gp.AddError( "\nError configuring output resolution: " + gp.GetMessages(2))
        raise Exception
    


    # Extract IGBP land cover layer from each HDF in the input directory
    try:

        gp.AddMessage("\nExtracting IGBP land cover...")

        # List all HDFs in input dir - ListRasters doesn't work 
        hdf_list = glob.glob('G:\\GIS\\Global\\Input\\MCD12Q1_global_landcover\\2001.01.01\\*.hdf')

        # Extract subdataset 0 - Land_Cover_Type_1 - IGBP classification
        gp.workspace = outputws
        for hdf in hdf_list:
            basename = os.path.basename(hdf)
            filename, ext = os.path.splitext(basename)
            imgname = filename + "_IGBP.img"
            gp.ExtractSubdataset_management(hdf, imgname, "0")
            mosaic_filenames = mosaic_filenames + imgname + ";"
            gp.AddMessage("Extracted = " + str(imgname))
                        
    except:
        gp.AddError("\nError extracting IGBP land cover: " + gp.GetMessages(2))
        raise Exception

    # Final values = energy and value over lifetime of dam

    try:

        gp.AddMessage("\nMosaicking land cover tiles...")

        # Remove trailing semi-colon from mosaic file list
        mosaic_filenames = mosaic_filenames[:-1]

        dsc = gp.describe(ref_proj)
        proj = dsc.SpatialReference
        
        gp.MosaicToNewRaster_management(mosaic_filenames, outputws, mosaic_2001, proj)
        

        
    except:
        gp.AddError("\nError mosaicking land cover tiles: " + gp.GetMessages(2))
        raise Exception




    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(outputws + "\\MCD12Q1_Preprocess_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("MCD12Q1 PRE-PROCESSING MODEL PARAMETERS\n")
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
        gp.Delete_management(interws)
        
    except:
        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception

    # REMOVE THIS
    gp.AddMessage("\nIGNORE THE FOLLOWING ERROR\n")
    raise Exception


except:
    gp.AddError ("\nError running script")
    raise Exception
