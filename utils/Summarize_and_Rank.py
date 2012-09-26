# ---------------------------------------------------------------------------
# Summarize_and_Rank.py
#
# Coded by Stacie Wolny
# for the Natural Capital Project
# 
# Last edit: 3/8/2010 
#
# Summarizes InVEST output by user-defined zone (watershed, parcel, etc)
# Ranks zones by summarized value
#
# ** To do: separate out ranking, add option for ascending/descending
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

        # InVEST data to summarize and rank
        in_raster = gp.GetParameterAsText(1)
        parameters.append("Data to summarize: " + in_raster)

        # Zones to summarize by (watersheds, parcels, etc)
        zones = gp.GetParameterAsText(2)
        parameters.append("Zones: " + zones)

        # Field in zone layer to use as identifiers in the summary
        zone_field = gp.GetParameterAsText(3)
        parameters.append("Zone field: " + zone_field)

##        # Save summary to a table?
##        do_table = gp.GetParameterAsText(4)
##        if do_table =='true':
##            do_table = True
##            parameters.append("Save to table: Yes")
##        else:
##            do_table = False
##            parameters.append("Save to table: No")

        # Suffix to append to output filenames, as <filename>_<suffix>
        Suffix = gp.GetParameterAsText(4)
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
        thefolders=["Intermediate", "Post_Process"]
        for folder in thefolders:
            if not gp.Exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        gp.AddError("Error creating folders: " + gp.GetMessages())
        raise Exception


    # Local variables
    try:
        # Base output/working directories
        postprocws = gp.workspace + os.sep + "Post_Process" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep

        # Intermediate variables
        
        # Output layers
        summary = postprocws + "summary"
        ranking = postprocws + "ranking"
        ranking_table_name = "ranking" + Suffix + ".dbf"
        summary_table_name = "summary_table" + Suffix + ".dbf"
        
    except:
        gp.AddError("Error configuring local variables: " + gp.GetMessages(2))
        raise Exception


    # Add suffix to end of output filenames
    # Constrain length of output raster filenames to 13 characters
    try:
        Outputnames = [summary, ranking]           
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

        summary = str(Outputnames[0])
        ranking = str(Outputnames[1])
            
    except:
        gp.AddError ("Error validating output filenames: " + gp.GetMessages(2))
        raise Exception


    # Summarize input data by zone
    try:

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws
        
        gp.AddMessage ("\nSummarizing by zone...")
        gp.ZonalStatistics_sa(zones, zone_field, in_raster, summary, "SUM", "DATA") 
        gp.AddMessage("\tCreated summary raster: \n\t" + str(summary))
        gp.CreateTable_management(postprocws, summary_table_name)
        summary_table = postprocws + summary_table_name
        gp.AddMessage("\nSummarizing to a table...")
        gp.ZonalStatisticsAsTable_sa(zones, zone_field, in_raster, summary_table, "DATA") 
        gp.AddMessage("\tCreated summary table: \n\t" + str(summary_table))
    except:
        gp.AddError ("Error summarizing by zone: " + gp.GetMessages(2))
        raise Exception


    # Rank summarized zones
    try:
        gp.AddMessage ("\nRanking...")

        # Create empty table, where ranked table will be written
        gp.CreateTable_management(interws, "rank_tmp.dbf", summary_table)
        rank_tmp_table = interws + "rank_tmp.dbf"

        # Fields to preserve for the output table
        summary_table_fields = ["VALUE", "COUNT", "AREA", "SUM"]

        # Sort input table by zonal sum, write to temporary sorted table 
        srows = gp.SearchCursor(summary_table, "", "", "", "SUM A")
        srow = srows.Next()
        irows = gp.InsertCursor(rank_tmp_table)
        while srow:
            irow = irows.newRow()
            for field in summary_table_fields:
                irow.setValue(field, srow.getValue(field))
            irows.insertRow(irow)
            srow = srows.Next()
        del srow, srows, irow, irows

        # Add ranking field to sorted table
        gp.AddField(rank_tmp_table, "rank", "LONG")
        
        urows = gp.UpdateCursor(rank_tmp_table)
        urow = urows.Next()
        rank_val = 1

        # Add ranking value to each row in table
        while urow:
            urow.rank = rank_val
            urows.UpdateRow(urow)
            urow = urows.Next()
            rank_val += 1

        # Remove Zonal Statistics fields that we don't need 
        gp.DeleteField_management(rank_tmp_table, "MIN; MAX; RANGE; MEAN; STD")
        gp.CreateTable_management(postprocws, ranking_table_name)
        ranking_table = postprocws + ranking_table_name
        gp.CopyRows_management(rank_tmp_table, ranking_table)
        gp.Delete_management(rank_tmp_table)
        
        gp.AddMessage("\n\tCreated ranking table output file: \n\t" + str(ranking_table))
    except:
        gp.AddError ("Error calculating rank: " + gp.GetMessages(2))
        raise Exception


##    # Write input parameters to an output file for user reference
##    try:
##        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
##        gp.workspace = gp.GetParameterAsText(0)
##        parafile = open(gp.workspace + "\\Output\\PP_Summarize_Rank_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
##        parafile.writelines("POST-PROCESSING - SUMMARIZE AND RANK\n")
##        parafile.writelines("____________________________________\n\n")
##
##        for para in parameters:
##            parafile.writelines(para + "\n")
##            parafile.writelines("\n")
##        parafile.close()
##    except:
##        gp.AddError ("Error creating parameter file:" + gp.GetMessages(2))
##        raise Exception


    # Clean up temporary files
    gp.AddMessage("\nCleaning up temporary files...\n")
    try:
        gp.Delete_management(interws)
    except:
        gp.AddError("Error cleaning up temporary files:  " + gp.GetMessages(2))
        raise Exception


##    # REMOVE THIS
##    gp.AddMessage("\nIGNORE THE FOLLOWING ERROR\n")
##    raise Exception

except:
    gp.AddError ("Error running script")
    raise Exception
