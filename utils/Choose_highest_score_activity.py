# ---------------------------------------------------------------------------
# Choose_highest_score_activity.py
# Created on: 12/15/2011
# by Stacie Wolny
# for the Natural Capital Project
#   
# Input: a list of water fund activities and associated rank/score rasters,
# where some cells may be included in multiple activity layers.
#
# NOTE: Score/activity layers are entered in increasing order of cost, so
# if there's a tie for score, the lower-cost activity will be chosen.
# Maybe change this to have the order info as inputs instead of hard-coded.
#
# NOTE: Activity layers must have zeroes assigned for non-activity cells
#
# Compare all scores for each cell, choose the activity cell with the highest
# score, de-select the other (coincident) activity layer cells with lower scores.
#
# Output: The modified activity maps, to be the final maps used in budget selection.
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


# Script arguments...
try:
    try:
        # Make parameters array, and later write input parameter values to an output file
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))

        # Dictionaries to hold scores/activities and their order of importance (key)
        activity_rasters = {}
        score_rasters = {}
        
        gp.AddMessage("\nGetting input arguments...")
        
        gp.workspace = sys.argv[1]
        parameters.append("Workspace: " + gp.workspace)

        act_protect = sys.argv[2]
        parameters.append("Protection activity raster: " + act_protect)

        act_enrich = sys.argv[3]
        parameters.append("Enrichment activity raster: " + act_enrich)

        act_fence = sys.argv[4]
        parameters.append("Fencing activity raster: " + act_fence)
        
        act_reforest = sys.argv[5]
        parameters.append("Reforest activity raster: " + act_reforest)

        act_silvo = sys.argv[6]
        parameters.append("Silvopasture activity raster: " + act_silvo)

        score_protect = sys.argv[7]
        parameters.append("Protection score raster: " + score_protect)
               
        score_enrich = sys.argv[8]
        parameters.append("Enrichment score raster: " + score_enrich)

        score_fence = sys.argv[9]
        parameters.append("Fencing score raster: " + score_fence)
        
        score_reforest = sys.argv[10]
        parameters.append("Reforest score raster: " + score_reforest)
        
        score_silvo = sys.argv[11]
        parameters.append("Silvopasture score raster: " + score_silvo)

        watersheds = sys.argv[12]
        parameters.append("Watersheds shapefile: " + watersheds)
        
        Suffix = sys.argv[13]
        parameters.append("Suffix: " + Suffix)

        if (Suffix == "") or (Suffix == string.whitespace) or (Suffix == "#"):
            Suffix = ""
        else:
            Suffix = "_" + Suffix

    except:
        gp.AddError("\nError in input arguments: " + gp.GetMessages(2))
        raise Exception
    

    try:
        gp.AddMessage("\nCreating output folders...")
        thefolders=["Intermediate", "Output"]
        for folder in thefolders:
            if not gp.exists(gp.workspace + folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        gp.AddError("\nError creating output folders: " + gp.GetMessages(2))
        raise Exception

    # Local variables...
    try:

        # Intermediate and output directories
        outputws = gp.workspace + os.sep + "Output" + os.sep
        interws = gp.workspace + os.sep + "Intermediate" + os.sep

        # Temporary variables
        score_fence_act = interws + "sc_fen_act"
        score_protect_act = interws + "sc_pro_act"
        score_enrich_act = interws + "sc_enr_act"
        score_reforest_act = interws + "sc_ref_act"
        score_silvo_act = interws + "sc_sil_act"

        act_fence_ws = interws + "act_fen_ws"
        act_protect_ws = interws + "act_pro_ws"
        act_enrich_ws = interws + "act_enr_ws"
        act_reforest_ws = interws + "act_ref_ws"
        act_silvo_ws = interws + "act_sil_ws"

        sc_fence_ws = interws + "sc_fen_ws"
        sc_protect_ws = interws + "sc_pro_ws"
        sc_enrich_ws = interws + "sc_enr_ws"
        sc_reforest_ws = interws + "sc_ref_ws"
        sc_silvo_ws = interws + "sc_sil_ws"       

        # Output files
        highest_score_pos = outputws + "hisc_pos"
        
    except:
        gp.AddError("\nError configuring local variables: " + gp.GetMessages(2))
        raise Exception

    try:

        gp.AddMessage("\nComparing activities by score...")

        # Make sure all temporary files go in the Intermediate folder
        gp.workspace = interws

        # Include score in comparison only where activities occur
        # Where activities don't occur, set their score to 0 so they
        #   won't skew the HighestPostition results

        # DO THIS MORE ELEGANTLY - TOTAL QUICK HACK RIGHT NOW

        gp.ExtractByMask_sa(act_fence, watersheds, act_fence_ws)
        gp.ExtractByMask_sa(act_protect, watersheds, act_protect_ws)
        gp.ExtractByMask_sa(act_enrich, watersheds, act_enrich_ws)
        gp.ExtractByMask_sa(act_reforest, watersheds, act_reforest_ws)
        gp.ExtractByMask_sa(act_silvo, watersheds, act_silvo_ws)

        gp.ExtractByMask_sa(score_fence, watersheds, sc_fence_ws)
        gp.ExtractByMask_sa(score_protect, watersheds, sc_protect_ws)
        gp.ExtractByMask_sa(score_enrich, watersheds, sc_enrich_ws)
        gp.ExtractByMask_sa(score_reforest, watersheds, sc_reforest_ws)
        gp.ExtractByMask_sa(score_silvo, watersheds, sc_silvo_ws)
        
        gp.SingleOutputMapAlgebra_sa("CON(" + act_fence_ws + " > 0," + sc_fence_ws + ", 0)", score_fence_act)
        gp.SingleOutputMapAlgebra_sa("CON(" + act_protect_ws + " > 0," + sc_protect_ws + ", 0)", score_protect_act)
        gp.SingleOutputMapAlgebra_sa("CON(" + act_enrich_ws + " > 0," + sc_enrich_ws + ", 0)", score_enrich_act)
        gp.SingleOutputMapAlgebra_sa("CON(" + act_reforest_ws + " > 0," + sc_reforest_ws + ", 0)", score_reforest_act)
        gp.SingleOutputMapAlgebra_sa("CON(" + act_silvo_ws + " > 0," + sc_silvo_ws + ", 0)", score_silvo_act)

        score_rasters[1] = score_protect_act
        score_rasters[2] = score_enrich_act
        score_rasters[3] = score_fence_act
        score_rasters[4] = score_reforest_act
        score_rasters[5] = score_silvo_act

        activity_rasters[1] = act_protect_ws
        activity_rasters[2] = act_enrich_ws
        activity_rasters[3] = act_fence_ws
        activity_rasters[4] = act_reforest_ws
        activity_rasters[5] = act_silvo_ws
                    
    

        # Make a list of score rasters to be used in Spatial Analyst tools
        score_ras_vals = score_rasters.values()
        activity_ras_vals = activity_rasters.values()

        score_ras_list = ""
        for score in score_ras_vals:
            score_ras_list = score_ras_list + score + ";"
        # Remove trailing ;
        score_ras_list = score_ras_list[:-1]

        # Choose which score raster has the highest score
        # In case of a tie, HighestPosition takes the first raster in 
        #  the list with the tie score (input list is in order of increasing cost,
        #  so less expensive activities are chosen first, an ROI approach)
        gp.HighestPosition_sa(score_ras_list, highest_score_pos)

        # Loop through all activity rasters
        # Where the activity list position number (dict key) is not equal to highest_score_pos,
        #   set the activity to NoData.  ie, only keep overlapping cells where the activity score
        #   is highest of all.
        
        for akey, araster in activity_rasters.iteritems():
            # Create output filename with suffix
            try:
                # One intermediate filename showing deselected cells - CHANGE THIS TO A BETTER NAME
                out_activity = outputws + os.path.basename(araster)
                # One final filename (deselected cells turned into NoData)
                out_activity0 = interws + os.path.basename(araster)
                Outputnames = [out_activity0]
                num = 0
                iSuffix = "_ds"
                for x in Outputnames:
                    y = x.split("\\")
                    z = y[-1]
                    lenz = len(z + iSuffix)
                    if lenz > 13:
                        x2 = x.rstrip(z)
                        overlen = 13 - lenz 
                        newz = z[0:overlen]
                        x2 = x2 + newz
                        Outputnames[num] = x2 + iSuffix
                    else:
                        Outputnames[num] = x + iSuffix
                    num = num + 1       
                out_activity0 = str(Outputnames[0])
            except:
                gp.AddError ("\nError creating output activity filenames: " + gp.GetMessages(2))
                raise Exception

            # Set changed cells to -999 so it's easy to see which cells were deselected due to score comparison
            gp.SingleOutputMapAlgebra_sa("CON(" + highest_score_pos + " <> " + str(akey) + ", -999, " + araster + ")", out_activity0)
            gp.Reclassify_sa(out_activity0, "VALUE", "0 0 NoData;-999 -999 NoData", out_activity, "DATA")
            gp.AddMessage("\nCreated output activity file: " + out_activity)

    except:
        gp.AddError("\nError comparing activities: " + gp.GetMessages(2))
        raise Exception

    # Write input parameters to an output file for user reference
    try:
        parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
        gp.workspace = gp.GetParameterAsText(0)
        parafile = open(gp.workspace + "\\Output\\Choose_highest_score_activity_" + now.strftime("%Y-%m-%d-%H-%M") + Suffix + ".txt", "w")
        parafile.writelines("CHOOSE HIGHEST SCORE ACTIVITY\n")
        parafile.writelines("_____________________________\n\n")

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
##        del ws_row, ws_rows
##        gp.Delete_management(interws)
##    except:
##        gp.AddError("\nError cleaning up temporary files:  " + gp.GetMessages(2))
##        raise Exception

    

except:
    gp.AddError("\nError running script: " + gp.GetMessages(2))
    raise Exception
