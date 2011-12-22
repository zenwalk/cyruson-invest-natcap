# ---------------------------------------------------------------------------
# pollination.py
# Coded by: Nasser Olwero (nasser.olwero@wwfus.org)
# Contributor: Taylor Ricketts (taylor.ricketts@wwfus.org)
# Contributor: Eric Lonsdorf (EricLonsdorf@lpzoo.org)
# Contributor: Christina Kennedy (christina.marie.kennedy@gmail.com)
# Date: May 26, 2008
# Date: March 4, 2010
# Date: Feb 2, 2011
# Description
# InVEST pollination model uses estimates of the availability of nesting sites and floral resources, 
# as well as flight ranges of bees, to create a map of an index of bee abundance at nest sites across 
# the landscape (i.e. Pollinator supply).  The model then uses the pollinator supply map and flight ranges 
# to predict an index (scaled from o to 1) of bee abundance visiting agricultural parcels on a landscape 
# (i.e. Farm abundance).  In an optional step, the model then calculates a simple index of the value of bees 
# to agricultural production, and attributes this value back to source parcels, creating a pollinator service 
# value map.
# ---------------------------------------------------------------------------

# Import system modules
from __future__ import division
#check if arcpy exists, if so import it and use it else use gp
try:
    import arcpy
    from arcpy.sa import *
    import invest10
    gp = None
except ImportError:
    import arcgisscripting
    gp = arcgisscripting.create()
    import invest
import numpy
from numpy import *    

import time, sys, datetime, os
now = datetime.datetime.now()

msgArguments = "Problem with arguments"
msgRasterDescription = "Error reading raster properties"
msgSeasonsRaster = "Error creating seasons raster"

########################################################################################################################
## Configuration
########################################################################################################################
# Sets the level of reporting on the dialog. True=print reports, False=Do not print
verbose = True
# Sets the option to cleanup intermediate files. True=delete temporary files, False=Do not delete temporary files
cleanup = True
# Toggles debug mode.  In debug mode the parameters are set in this file and model can be run from python.
debug = False    
#initialize cleanup list.  All files to be cleaned after run are added here
cleanuplist=[]
#initialize os cleanup list.  All files to be cleaned after run are added here
cleanuposlist=[]    
#set scaling factor to adjust for fractionals. change this to cater for smaller values.  ArcGIS does not perform lookups correctly with fractions
scaling=100.0
#use square neighborhoods.  If this is set to tru then when creating distance array around a focal cell a square pattern is used
squareneighborhood = True
#the first 4 characters of the species names are similar.  Do not change this flag, used to correct naming
sppdup = 0

########################################################################################################################
## Common Functions - functions common to 9.3 and 10
########################################################################################################################
#function to select unique items from list for pollinator abundance on farms
def unique(seq, idfun=None):
    if idfun is None:
        def idfun(x): return x
    seen={}
    result=[]
    for item in seq:
        marker=idfun(item)
        if marker in seen:continue
        seen[marker]=1
        result.append("frm_"+item[:4])
    return result
#----------------------------------------------------------------------    
# function to check if value is a number stored as string
def sppNameToText(theVar):
    try:
        newnum = int(theVar)
    except ValueError:
        try: 
            newnum = float(theVar)
            return (int(newnum))
        except:
            return False
    else:
        return newnum


#select the correct code to run based on the version on ArcGIS user is running

if gp:
    ######################################################################################################################## 
    ########################################################################################################################
    ## ArcGIS 9.3
    ########################################################################################################################  
    ########################################################################################################################     
    try:
        from osgeo import ogr, gdal
        from osgeo.gdalconst import *
    except:
        raise Exception, "Numpy and GDal python libraries are required.  Please refer to the manual for instructions."
    gdal.AllRegister()
    now = datetime.datetime.now()
    
   
    # Overwrite output
    gp.overwriteoutput=1
    
    # Check out any necessary licenses
    gp.CheckOutExtension("spatial")
    
    msgArguments = "Problem with arguments"
    msgRasterDescription = "Error reading raster properties"
    msgSeasonsRaster = "Error creating seasons raster"
    
    try:
        ########################################################################################################################
        ## Configuration
        ########################################################################################################################
        # Sets the level of reporting on the dialog. True=print reports, False=Do not print
        verbose = True
        # Sets the option to cleanup intermediate files. True=delete temporary files, False=Do not delete temporary files
        cleanup = True
        # Toggles debug mode.  In debug mode the parameters are set in this file and model can be run from python.
        debug = False    
        #initialize cleanup list.  All files to be cleaned after run are added here
        cleanuplist=[]
        #initialize os cleanup list.  All files to be cleaned after run are added here
        cleanuposlist=[]    
        #set scaling factor to adjust for fractionals. change this to cater for smaller values.  ArcGIS does not perform lookups correctly with fractions
        scaling=100.0
        #use square neighborhoods.  If this is set to tru then when creating distance array around a focal cell a square pattern is used
        squareneighborhood = True
        #the first 4 characters of the species names are similar.  Do not change this flag, used to correct naming
        sppdup = 0
        
        ########################################################################################################################
        ## Get parameters
        ########################################################################################################################
        parameters = []
        if debug:
            print "debug mode"
        else:
            
            now = datetime.datetime.now()
            parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
            gp.workspace = gp.GetParameterAsText(0)
            original_workspace = gp.workspace
            parameters.append("Workspace: "+gp.Workspace) 
            lcraster = gp.GetParameterAsText(1)
            parameters.append("Current Landcover: "+ lcraster)
            CDate= gp.GetParameterAsText(2)
            parameters.append("Current Landcover Date: "+ str(CDate))
            if CDate:
                CDate = int(gp.GetParameterAsText(2))
            resolution = gp.GetParameterAsText(3)
            parameters.append("Resolution: "+ str(resolution))
            if resolution:
                resolution = int(gp.GetParameterAsText(3))
            lcrasterf = gp.GetParameterAsText(4)
            parameters.append("Future Landcover: "+lcrasterf)        
            FDate= gp.GetParameterAsText(5)
            parameters.append("Future Landcover Date: "+ str(FDate))
            if FDate:
                FDate = int(gp.GetParameterAsText(5))
            guild = gp.GetParameterAsText(6)
            parameters.append("Guild Table: "+ guild)
            lutable = gp.GetParameterAsText(7)
            parameters.append("Landuse Parameters Table: "+ lutable)
            farmstable = ""
            ag=gp.GetParameterAsText(8)
            parameters.append("Agricultural Classes: "+ ag)
            sat_const = gp.GetParameterAsText(9)
            parameters.append("Saturation Constant: "+ sat_const)
            wild_prop = gp.GetParameterAsText(10)
            parameters.append("Proportion of crop wild pollinated: "+ str(wild_prop))            
            if ag:
                ag = ag.replace(",",";")
                ag = ag.split(';')
            suffix=gp.GetParameterAsText(11)
            if suffix is "#":
                suffix = ""
            parameters.append("Suffix: "+str(suffix))
            if lcrasterf:
                lcovers=["Current","Future"]
            else:
                lcovers=["Current"]
            realized = gp.GetParameterAsText(12)
        if realized == 'true':
            realized = True
        else: 
            realized = False
        
        
    
            
        ########################################################################################################################
        ## Functions
        ########################################################################################################################
        #----------------------------------------------------------------------
        #function to build list of pollinators required by given crop
        def checkpollinators(crop):
            cur=gp.searchcursor(guild)
            row=cur.next()
            pollinators=[]
            while row:
                if row.GetValue(crop):
                    pollinators.append(row.species)
                row=cur.next()
            del cur
            return pollinators
        #----------------------------------------------------------------------
        #check projections
        def checkProjections(thedata):
            dataDesc = gp.describe(thedata)
            spatreflc = dataDesc.SpatialReference
            if spatreflc.Type <> 'Projected':
                gp.Error(thedata +" Your data : "+thedata+" is not projected.  You must use a projected raster and projection must be defined.")
            elif spatreflc.LinearUnitName <> 'Meter':
                gp.AddWarning("This model assumes that "+thedata+" is projected in meters.  You may get erroneous results")
            if thedata<>lcraster:
                if gp.describe(lcraster).SpatialReference.name<>spatreflc.name and spatreflc.Type == 'Projected':
                    gp.AddError("WARNING: "+ thedata +" is in a different projection from the land cover data.  You may get erroneous results.")
        #----------------------------------------------------------------------    
        # function to toggle from main folder to intermediate folder for reading and writing data
        down="down"
        def switchws(direction="up"):
            if direction=="down":
                gp.workspace=original_workspace+"\\intermediate"
            else:
                gp.workspace=original_workspace 
                
        #check and create folders
        invest.createfolders(gp.workspace)

        
        ########################################################################################################################
        ## Process each landcover given in turn
        ## Steps
        ## 1. 
        ########################################################################################################################
        
        for cover in lcovers:
            #reset suffix
            suffix=""
            ##set suffix based on the cover
            if cover=="Future":
                lcraster=lcrasterf
                suffix = "_fut"+suffix
            else:
                suffix = "_cur"+suffix
                
            gp.addmessage("Processing: " + cover + " landcover")
            
            
            #----------------------------------------------------------------------  
            # Get raster properties of the cover
            rasterDesc=gp.describe(lcraster)
            rasterExtent = rasterDesc.extent
            CellSize = rasterDesc.MeanCellHeight
            gp.Extent=rasterExtent
            lcwidth = rasterDesc.width
            lcheight = rasterDesc.height
            gp.mask = lcraster
            gp.CellSize=lcraster
            
            if resolution:
                gp.CellSize=int(resolution)
            checkProjections(lcraster)
            if lcrasterf:
                checkProjections(lcrasterf)
            
            switchws(down)
            #----------------------------------------------------------------------  
            #create agricultural areas map
            gp.CreateConstantRaster_sa ("unity", "1", "INTEGER", gp.Cellsize, gp.extent) 
            cleanuplist.append("unity")
            
            gp.addmessage("Creating agricultural areas map...")
    
            if ag:
                
                agfile = open(gp.workspace+"\\ag.asc","w")
                cur = gp.SearchCursor(lutable)
                row = cur.next()
                while row:
                    agfile.write(str(int(row.LULC))),
                    agfile.write(":"),
                    if str(int(row.LULC)) in ag:
                        agfile.write("1\n")
                    else:
                        agfile.write("0\n")
                    row = cur.next()
                agfile.close()
                #create an agricultural areas map using the text file created above
                filecheck = invest.checkASCIIFile(gp.workspace+"\\ag.asc")
                if filecheck.hasDuplicates():
                    raise Exception, "Your landuse classes file has duplicates"
                if not filecheck.isSorted():
                    raise Exception, "Your landuse file is not sorted"
                gp.ReclassByASCIIFile_sa(lcraster, gp.workspace+"\\ag.asc", "agmap", "NODATA")
            else:
                #create constant raster
                gp.extent = lcraster
                #gp.copy_management(gp.workspace+"\\unity", gp.workspace+"\\agmap")
                gp.SingleOutputMapAlgebra_sa("unity", "agmap")
    
                
               
            ########################################################################################################################
            ## This section calculates the land cover nest suitability and creates a grid for each species
            ########################################################################################################################
            ## Initial validation check and preparation
            ########################################################################################################################
            
            # Its assumed nest type fields will be named N_Nestx where x is an index        
            # ----------------------------------------------------------------------------------------
            #create array of nests from land cover table
            NestFields = gp.listfields(lutable, "N_*", "All")
            field = NestFields.next()
            if not field:
                raise Exception, "Nest data column (N_)required in land cover table"
            nnests=[]
            while field:
                nnests.append(field.name)
                field = NestFields.next()
            # ----------------------------------------------------------------------------------------
            #count number of bees/pollinatios and build list of pollinators
            guildcursor = gp.searchcursor(guild)
            beerow = guildcursor.next()
            numberofguilds = 0
            pollinatorlist = []
            while beerow:
                numberofguilds += 1
                pollinatorlist.append(str(beerow.Species))
                beerow = guildcursor.next()
                
            # ----------------------------------------------------------------------------------------
            #create array of nests from guild table
            NestFields = gp.listfields(guild, "NS_*", "All")
            field = NestFields.next()
            if not field:
                raise Exception, "Nest data column (NS_) required in guild table"
            nsnests=[]
            while field:
                nsnests.append(field.name)
                field = NestFields.next()
            
            # ----------------------------------------------------------------------------------------
            #create array of seasons from land cover table
            SeasonFields = gp.listfields(lutable, "F_*", "All")
            field = SeasonFields.next()
            if not field:
                raise Exception, "Season data column(F_) required in landuse table"
            lcoverseasons=[]
            while field:
                lcoverseasons.append(field.name)
                field = SeasonFields.next()
            
            # ----------------------------------------------------------------------------------------
            #create array of nests from guild table
            SeasonFields = gp.listfields(guild, "FS_*", "All")
            field = SeasonFields.next()
            if not field:
                raise Exception, "Season data column (FS_) required"
            guildseasons=[]
            while field:
                guildseasons.append(field.name)
                field = SeasonFields.next()
            
            # ----------------------------------------------------------------------------------------    
            #ensure that the nest count and nest names match in landuse and guild tables 
            if len(nsnests)!=len(nnests):
                raise Exception, "Nest mismatch"
            else:
                for x in nnests:
                    x = x.upper()[2:]
                    if not x in [x.upper()[3:] for x in nsnests]:
                        raise Exception, "Nest mismatch"
            
            # ----------------------------------------------------------------------------------------
            #ensure that the seasons count and seasons names match in landuse and guild tables 
            if len(lcoverseasons)!=len(guildseasons):
                raise Exception, "Seasons mismatch seasons columns (FS_ and F_) must match in guilds and landuse tables"
            else:
                for x in lcoverseasons:
                    x = x.upper()[2:]
                    if not x in [x.upper()[3:] for x in guildseasons]:
                        raise Exception, "Seasons mismatch seasons columns (FS_ and F_) must match in guilds and landuse tables"
        
            #----------------------------------------------------------------------    
            # check if there are species duplicate names and if first 4 characters are unique
            if invest.ListHasDuplicates(pollinatorlist):
                raise Exception, "Pollinator names must be unique, exiting."
            spp = [x[:4] for x in pollinatorlist]  
            if invest.ListHasDuplicates(spp):
                    sppdup = 1
            ########################################################################################################################
            ## 
            ########################################################################################################################
            # ----------------------------------------------------------------------------------------
            #create nesting availability raster for each nest type
            for nnest in nnests:
                try:
                    nestingfile = open(gp.workspace+"\\"+nnest+"values.asc","w")
                    cleanuposlist.append(nnest+"values.asc")
                    cur = gp.SearchCursor(lutable)
                    row = cur.next()
                    while row:
                        nestingfile.write(str(int(row.LULC))),
                        nestingfile.write(":"),
                        nestingfile.write(str(int(row.GetValue(nnest)*scaling))+"\n")
                        row = cur.next()
                    nestingfile.close()
                    filecheck = invest.checkASCIIFile(gp.workspace+"\\"+nnest+"values.asc")
                    if filecheck.hasDuplicates():
                        raise Exception, "Your landuse classes file has duplicates"
                    if not filecheck.isSorted():
                        raise Exception, "Your landuse file is not sorted"
                    gp.ReclassByASCIIFile_sa(lcraster, gp.workspace+"\\"+nnest+"values.asc", nnest+"temp", "NODATA")
                    del cur
                    
                    # scale back since it was multiplied by scaling factor specified in the configuration
                    expr = nnest+"temp / "+str(scaling)
                    gp.SingleOutputMapAlgebra_sa(expr, nnest)
                    cleanuplist.append(nnest)
                    if cleanup:
                        if gp.exists(nnest+"temp"):
                            gp.delete_management(nnest+"temp")
                except:
                    raise Exception, "Error sup_ting raster"
            
            # ----------------------------------------------------------------------------------------          
            #Create nesting maps
            gp.addmessage("Creating nesting maps...")
            cur=gp.searchcursor(guild)
            row=cur.next()
            while row:
                maxexpr="max("
                if sppNameToText(row.Species):
                    focalspp = 's'+str(int(row.Species))
                    focalsppfullname = 'species ' + str(int(row.Species))
                    focalspp = focalspp[:4]
                else:
                    focalspp = row.Species[:4]
                    focalsppfullname = row.Species
                if sppdup:
                    focalsppfullname = row.Species
                    focalspp = row.Species[:3]+ str(row.OID+1) 
                gp.addmessage("   "+str(focalsppfullname))
            
                for nsnest in nsnests:
                    expr="N_"+nsnest[3:]+" * "+str(row.GetValue(nsnest))
                    gp.SingleOutputMapAlgebra_sa(expr, focalspp+nsnest[3:])
                    cleanuplist.append(focalspp+nsnest[3:])
                    maxexpr+=focalspp+nsnest[3:]+","
                    
                maxexpr=maxexpr[:-1]
                maxexpr+=")"
                gp.SingleOutputMapAlgebra_sa(maxexpr, "hn_"+focalspp+suffix)
                row=cur.next()
            del cur
            
            # ----------------------------------------------------------------------------------------            
            #create floral availability raster for each season
            gp.addmessage("Creating floral resources map...")
            for season in lcoverseasons:
                try:
                    seasonfile = open(gp.workspace+"\\"+season+"values.asc","w")
                    cleanuposlist.append(season+"values.asc")
                    cur = gp.SearchCursor(lutable)
                    row = cur.next()
                    while row:
                        seasonfile.write(str(int(row.LULC))),
                        seasonfile.write(":"),
                        seasonfile.write(str(int(row.GetValue(season)*scaling))+"\n")
                        row = cur.next()
                    seasonfile.close() 
                    filecheck = invest.checkASCIIFile(gp.workspace+"\\"+season+"values.asc")
                    if filecheck.hasDuplicates():
                        raise Exception, "Your landuse classes file has duplicates"
                    if not filecheck.isSorted():
                        raise Exception, "Your landuse file is not sorted"
                    gp.ReclassByASCIIFile_sa(lcraster, gp.workspace+"\\"+season+"values.asc", season+"temp", "DATA")
                    del cur
                    
                    # scale back since it was multiplied by scaling factor specified in the configuration
                    expr = season+"temp / "+str(scaling)
                    gp.SingleOutputMapAlgebra_sa(expr, season)
                    if gp.exists(season+"temp"):
                        if cleanup:
                            gp.delete_management(season+"temp")
                except:
                    raise Exception, msgSeasonsRaster
                
            ########################################################################################################################
            ## process each species or guild table row by row (from the guild table) and for each season
            ########################################################################################################################          
    
            rows=gp.searchcursor(guild)
            row=rows.next()
            PAFarms_expression="" #Used for summung pollinator availability on farms
            Pm_expression="(" #PAFarms expression builder, used later
            countspp=0
            while row: #for each row/species
                if sppNameToText(row.Species):
                    focalspp = 's'+str(int(row.Species))
                    focalsppfullname = 'species ' + str(int(row.Species))
                    focalspp = focalspp[:4]
                else:
                    focalspp = row.Species[:4]
                    focalsppfullname = row.Species
                if sppdup:
                    focalsppfullname = row.Species
                    focalspp = row.Species[:3]+ str(row.OID+1) 
                
                gp.addmessage("   "+focalsppfullname+"...")
                alpha = row.alpha
                #check that species name is a string
                if sppNameToText(row.Species):
                    focalspp = 's'+str(int(row.Species))
                    focalsppfullname = 'species ' + str(int(row.Species))
                    focalspp = focalspp[:4]
                else:
                    focalspp = row.Species[:4]
                    focalsppfullname = row.Species
                if sppdup:
                    focalsppfullname = row.Species
                    focalspp = row.Species[:3]+ str(row.OID+1) 
                hfexpr=""
                for season in guildseasons: #for each season
                    fs = row.GetValue("FS_"+season[3:])
                    if not fs > 0:
                        continue
                    gp.addmessage("     Season: "+season[3:])
                    
                    # ----------------------------------------------------------------------------------------  
                    # Calculate neighborhood based on alpha - used to create the distance matrix
                    nb=float(math.ceil(alpha))*2 
                    neighborhood=int(nb)
                    neighbors =math.ceil(neighborhood/float(gp.CellSize))
                    if not neighbors % 2:
                        neighbors += 1
                    neighbors = int(neighbors)
                    if neighbors > 101:
                        neighbors = 101
                    if neighbors < 5:
                        neighbors = 5
                    nbcells = (neighbors - 1) / 2
        
                    # ----------------------------------------------------------------------------------------  
                    #create a dist array which has weight values for running focal statistics.  The weight is 
                    #used to describe the exponential decline of pollinator visits with distance
                    
                    distplain=[]
                    neighborstring=str(neighbors)+" "
                    #step through the neighborhood array
                    for i in range(0,neighbors):
                        innerlist=[]
                        innerlistplain=[]
                        innerlistdenom = []
                        for j in range(0,neighbors):
                            thedist = math.sqrt((abs(i - nbcells) * float(gp.CellSize)) ** 2 + (abs(j - nbcells) * float(gp.CellSize)) ** 2)
                            if squareneighborhood:
                                weightdistplain = round(math.exp(-thedist / alpha),4)
                            else: 
                                if thedist > neighborhood / 2:
                                    weightdistplain = 0
                                else:
                                    weightdistplain = round(math.exp(-thedist / alpha),4)
                            innerlistplain.append(weightdistplain)
                        distplain.append(innerlistplain)
                    # ----------------------------------------------------------------------------------------                     
                    # write the weight distance files to a text file for focal statistics    
                    weightfile = open(gp.workspace+"\\weightpollen"+focalspp[:4]+".txt","w")
                    cleanuplist.append("weightpollen"+focalspp+".txt")
                    
                    weightfile.write(neighborstring*2+"\n")
                   # ----------------------------------------------------------------------------------------                  
                    #calculate the sum for the weights for normalizing the weight matrix
                    sumweight = 0
                    for x in distplain:
                        for y in x:
                            sumweight += y
                   # ---------------------------------------------------------------------------------------- 
                   #write weight file
                    for x in distplain:
                        for y in x:
                            if sumweight:
                                weightfile.write(str(round(y/float(sumweight), 6))+' ')
                                #if debug: num.write(str(round(y, 6))+' ')
                            else:
                                weightfile.write('0 ')
                        weightfile.write("\n")
                        #if debug: num.write("\n")
                        
                    weightfile.close()
                     
                    # ---------------------------------------------------------------------------------------- 
                    #calculate foraging map for this species
                    #gp.addmessage("  Creating foraging map (hf_*) for "+str(focalsppfullname)+"...")
                    fraster = "f_"+season[3:9]
                    cleanuplist.append("f_"+season[3:9])
                    OutRaster = gp.workspace+"\\hf"+focalspp+season[3:]
                    InNeighborhood = "WEIGHT "+gp.workspace+"\\weightpollen"+focalspp[:4]+".txt CELL"
                    InNoDataOption = "DATA"
                    # FocalStatistics
                    gp.FocalStatistics_sa(fraster, OutRaster, InNeighborhood, "SUM", InNoDataOption)
                    cleanuplist.append("hf"+focalspp+season[3:])
                    
                    #if there are multiple season take the minimum value
                    
                    if len(guildseasons)>1:#if more than one season
                        hfexpr+="hf"+focalspp+season[3:]+" * " + str(fs) + " + "
                    else:
                        hfexpr=""
                if hfexpr:
                    hfexpr = hfexpr[:-3]
                    
                    gp.SingleOutputMapAlgebra_sa(hfexpr, gp.workspace+"\\hf_"+focalspp+suffix) #take the minimum of the seasons
                else:
                    gp.SingleOutputMapAlgebra_sa("hf"+focalspp+guildseasons[0][3:], gp.workspace+"\\hf_"+focalspp+suffix) #just one season so take the only value
                
                # ---------------------------------------------------------------------------------------- 
                #calculate the pollinator abundance for this species (nesting map - determined by nest availability and floral resources availability within range 
                #of this pollinator flight distance - alpha)
                #gp.addmessage("  Creating nesting map (sup_*) for "+str(focalsppfullname)+"...")
                switchws()
                pexpr = gp.workspace+"\\intermediate\\hf_"+focalspp+suffix+" * " +gp.workspace+"\\intermediate\\hn_"+focalspp+suffix
                gp.SingleOutputMapAlgebra_sa(pexpr, "\\intermediate\\sup_"+focalspp[:4]+suffix)
                
                # ---------------------------------------------------------------------------------------- 
                #Calculate pollinator abundance on each cell abundance of each pollinator
                inraster = gp.workspace+"\\intermediate\\sup_"+focalspp[:4]+suffix
                OutRaster = gp.workspace+"\\intermediate\\frmt_"+focalspp[:4]
                InNeighborhood = "WEIGHT "+gp.workspace+"\\intermediate\\weightpollen"+focalspp[:4]+".txt CELL"
                InNoDataOption = "DATA"
                cleanuplist.append("frmt_"+focalspp[:4])
        
                # FocalStatistics
                gp.FocalStatistics_sa(inraster, OutRaster, InNeighborhood, "SUM", InNoDataOption)  
                switchws(down)
                
                #create mask
                gp.SingleOutputMapAlgebra_sa("1 - agmap", "aginvert")
                
                #remove non ag cells
                gp.SetNull_sa ("aginvert", "frmt_"+focalspp[:4], "frm_"+focalspp[:4], "") 
                
                #create an expression for abundance for all species        
                PAFarms_expression+= "frm_"+focalspp[:4]+" + "
    
                #increment species counter
                countspp+=1
        
                #update expression for calculating relative pollinator score
                switchws()
                Pm_expression+= gp.workspace+"\\intermediate\\sup_"+focalspp[:4]+suffix+" + "
                switchws(down)
                row=rows.next() # next species/guild
            del rows
            
            # ---------------------------------------------------------------------------------------- 
            #complete PAFarms expression builder (pollinator abundance on farms)
            PAFarms_expression = PAFarms_expression[:-3]
        
            # ---------------------------------------------------------------------------------------- 
            #complete Pm expression builder and calculate Pm
            Pm_expression = Pm_expression[:-3]+") / "+str(countspp)
            # ---------------------------------------------------------------------------------------- 
            #compute the total nesting map for all species/guilds
            switchws()
            gp.SingleOutputMapAlgebra_sa(Pm_expression, "\\output\\sup_tot"+suffix)
    
            switchws(down)
    
            gp.addmessage("Creating map of pollinator abundance...")
            # ---------------------------------------------------------------------------------------- 
            # Calculate pollinatior abundance on farms for all species/guilds (PAFarms)
            gp.SingleOutputMapAlgebra_sa(PAFarms_expression, "frm_tot")
            switchws()
            # ---------------------------------------------------------------------------------------- 
            #calculate average abundance (first statement muted to test with non ag pixels not removed)
            gp.SingleOutputMapAlgebra_sa(gp.workspace+"\\intermediate\\frm_tot / "+str(numberofguilds), gp.workspace+"\\output\\frm_avg"+suffix)
            if realized:
                # ---------------------------------------------------------------------------------------- 
                #create value map
                gp.addmessage("Creating value map...")
                
                Vm_expression =  "(1 - "+str(wild_prop)+ ") + ("+str(wild_prop) +" * ("+gp.workspace+"\\output\\frm_avg"+suffix+" / ("+sat_const+" + "+gp.workspace+"\\output\\frm_avg"+suffix+")))"
                
                gp.SingleOutputMapAlgebra_sa(Vm_expression, gp.workspace+"\\intermediate\\frm_val"+suffix)
        
                ########################################################################################################################
                ## Calculate the pollinator service value - this part uses gdal and numpy libraries
                ########################################################################################################################
            
                #read the projection info of the value raster
                ds = gdal.Open(gp.workspace+"\\intermediate\\frm_val"+suffix, GA_ReadOnly)   
                geoTransform = ds.GetGeoTransform()
                proj = ds.GetProjection()
                ds = None
                #process each species or guild table row by row
                rows=gp.searchcursor(guild)
                row=rows.next()
                switchws(down)
                gp.addmessage("Calculating pollinator service value...")
                speciesnames = []
                while row:
                    if sppNameToText(row.Species):
                        focalspp = 's'+str(int(row.Species))
                        focalsppfullname = 'species ' + str(int(row.Species))
                        focalspp = focalspp[:4]
                    else:
                        focalspp = row.Species[:4]
                        focalsppfullname = row.Species
                    if sppdup:
                        focalsppfullname = row.Species
                        focalspp = row.Species[:3]+ str(row.OID+1)    
                    gp.addmessage("   "+focalsppfullname)
                    
                    #get the contribution of this species
                    switchws()
                    contexpr = gp.workspace+"\\intermediate\\frm_val"+suffix+" * "+gp.workspace+"\\intermediate\\frm_"+focalspp[:4]+" / "+gp.workspace+"\\intermediate\\frm_tot"
                    gp.SingleOutputMapAlgebra_sa(contexpr, gp.workspace+"\\intermediate\\ctr_"+focalspp[:4]+suffix)
                    cleanuplist.append("ctr_"+str(focalspp[:4])+suffix)
                    
                    switchws(down)
                    val_raster = "ctr_"+focalspp[:4]+suffix
                    sup_raster = "sup_"+focalspp[:4]+suffix
                    frm_raster = "frm_"+focalspp[:4]
                    crop_raster = "agmap"
                    distfile = gp.workspace+"\\weightpollen"+focalspp[:4]+".txt"
                    gp.extent = lcraster
                    speciesnames.append("sp_"+focalspp[:4]+suffix+".img")
        
                    workingdir = gp.workspace
                    matrix = invest.readFileIntoArray(distfile, 1)
                    nst =  invest.convertRastertoArray(workingdir, sup_raster)
                    poll_abund =  invest.convertRastertoArray(workingdir, frm_raster)
                    value =  invest.convertRastertoArray(workingdir, val_raster)
                    agmap =  invest.convertRastertoArray(workingdir, crop_raster)
                    output = zeros(shape=value.shape)
                    outputsum = zeros(shape=value.shape)
                    
                    height = value.shape[0]
                    width = value.shape[1]
                    
                    sizematrix = matrix.shape[0]
                    nb = int(sizematrix / 2)
                    try:
                        for j in arange(height): #process each row from the top working downwards
                            for i in arange(width): #process each column from left to right
                                if agmap[j][i] > 0: #only process if its a cropped pixel, value > 0
                                    if i > nb-1 and i < (width - nb) and j > nb-1 and j < (height - nb): #'ensure the block does not go beyond raster
                                        for jstart in arange(sizematrix): # 0 To nb - 1 'loop through the neighborhood
                                            for istart in arange(sizematrix): #= 0 To nb - 1
                                                if poll_abund[j][i] > 0: # 'this must be greater than 0 to prevent overflow in next line. If its zero then set the result to zero
                                                    result = value[j][i] * matrix[jstart][istart] * nst[j - nb + jstart][i - nb + istart] / poll_abund[j][i]
                                                else: 
                                                    result = 0
                                                temphold = round(output[j - nb + jstart][i - nb + istart] + result, 6)
                                                if temphold < 0:
                                                    temphold = 0
                                                if temphold > 3.402823466e+38:
                                                    temphold = 0                                        
                                                output[j - nb + jstart][i - nb + istart] = temphold 
                                                
                        invest.convertArraytoRaster(output, "sp_"+focalspp[:4]+suffix+".img", workingdir, geoTransform, proj)
                    except:
                        raise Exception, "Error calculating realized supply for species: "+focalsppfullname
                    row=rows.next()
                # ----------------------------------------------------------------------------------------     
                #calculate total realized supply
                switchws()
                exprSupply = ''
                for x in speciesnames:
                    exprSupply += gp.workspace+"\\intermediate\\"+x+ " + "
                exprSupply = exprSupply[:-3]
                gp.addmessage("  Combining pollinator service value maps")
                gp.extent = lcraster
                gp.SingleOutputMapAlgebra_sa(exprSupply, gp.workspace+"\\output\\sup_val"+suffix+".img")
                # ---------------------------------------------------------------------------------------- 
            #write parameters
            switchws()
            #if suffix:
                #suffix = "_"+suffix
            switchws()
            parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
            
            parafile = open(gp.workspace+"\\output\\pollination_"+now.strftime("%Y-%m-%d-%H-%M")+suffix+".txt","w") 
            parafile.writelines("POLLINATION MODEL PARAMETERS\n")
            parafile.writelines("____________________________\n\n")
            
            for para in parameters:
                parafile.writelines(para+"\n")
                parafile.writelines("\n")
            parafile.close()
            switchws(down)
        #check if name speciesnames is defined.  
        try:
            speciesnames
        except NameError:
            speciesnames = []
            
        if cleanup:
            gp.addmessage("Cleaning up...")
            for theraster in cleanuplist+["PAFarmsDenom","PAFarmsNum","unityraster","numerator","denominator","ag.asc","tempo", "agmap", "aginvert"]+speciesnames:
                if gp.exists(theraster):
                    gp.delete_management(theraster)
            #delete non datasets
            os.chdir(gp.workspace)
            for thefile in cleanuposlist+["ag.asc"]:
                if os.path.exists(thefile):
                    os.remove(thefile)
        
    except Exception, ErrorDesc:
        gp.AddError(str(ErrorDesc))


        
        
        
        
        
else:
    ######################################################################################################################## 
    ########################################################################################################################
    ## ArcGIS 10
    ########################################################################################################################  
    ######################################################################################################################## 
    #user is running ArcGIS 10

    # Overwrite output 
    arcpy.env.overwriteOutput = 1
    # Check out any necessary licenses
    arcpy.CheckOutExtension("spatial")
    
    try:
       
        ########################################################################################################################
        ## Get parameters
        ########################################################################################################################
        parameters = []
        if debug:
            now = datetime.datetime.now()
            parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
            arcpy.env.workspace = r"C:\InVEST1005\pollination"
            original_workspace = arcpy.env.workspace
            parameters.append("Workspace: "+arcpy.env.workspace) 
            lcraster = r"C:\InVEST1005\Base_Data\lulc_samp_cur"
            parameters.append("Current Landcover: "+ lcraster)
            CDate= ''
            parameters.append("Current Landcover Date: "+ str(CDate))
            if CDate:
                CDate = int(arcpy.GetParameterAsText(2))
            resolution = '200'
            parameters.append("Resolution: "+ str(resolution))
            if resolution:
                resolution = int('200')
            lcrasterf = ''
            parameters.append("Future Landcover: "+lcrasterf)        
            FDate= ''
            parameters.append("Future Landcover Date: "+ str(FDate))
            if FDate:
                FDate = int('')
            guild = r"C:\InVEST1005\pollination\Input\Guild.dbf"
            parameters.append("Guild Table: "+ guild)
            lutable = r"C:\InVEST1005\pollination\Input\LU.dbf"
            parameters.append("Landuse Parameters Table: "+ lutable)
            farmstable = ""
            ag=r"67;68;71;72;73;74;75;76;78;79;80;81;82;83;84;85;88;90;91;92"
            parameters.append("Agricultural Classes: "+ ag)
            sat_const = "0.125"
            parameters.append("Saturation Constant: "+ sat_const)
            if ag:
                ag = ag.replace(",",";")
                ag = ag.split(';')
            suffix="#"
            if suffix is "#":
                suffix = ""
            parameters.append("Suffix: "+str(suffix))
            if lcrasterf:
                lcovers=["Current","Future"]
            else:
                lcovers=["Current"]
            realized = 'true'
        else:
            
            now = datetime.datetime.now()
            parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
            arcpy.env.workspace = arcpy.GetParameterAsText(0)
            original_workspace = arcpy.env.workspace
            parameters.append("Workspace: "+arcpy.env.workspace) 
            lcraster = arcpy.GetParameterAsText(1)
            parameters.append("Current Landcover: "+ lcraster)
            CDate= arcpy.GetParameterAsText(2)
            parameters.append("Current Landcover Date: "+ str(CDate))
            if CDate:
                CDate = int(arcpy.GetParameterAsText(2))
            resolution = arcpy.GetParameterAsText(3)
            parameters.append("Resolution: "+ str(resolution))
            if resolution:
                resolution = int(arcpy.GetParameterAsText(3))
            lcrasterf = arcpy.GetParameterAsText(4)
            parameters.append("Future Landcover: "+lcrasterf)        
            FDate= arcpy.GetParameterAsText(5)
            parameters.append("Future Landcover Date: "+ str(FDate))
            if FDate:
                FDate = int(arcpy.GetParameterAsText(5))
            guild = arcpy.GetParameterAsText(6)
            parameters.append("Guild Table: "+ guild)
            lutable = arcpy.GetParameterAsText(7)
            parameters.append("Landuse Parameters Table: "+ lutable)
            farmstable = ""
            ag=arcpy.GetParameterAsText(8)
            parameters.append("Agricultural Classes: "+ ag)
            sat_const = arcpy.GetParameterAsText(9)
            parameters.append("Saturation Constant: "+ sat_const)
            wild_prop = arcpy.GetParameterAsText(10)
            parameters.append("Proportion of crop wild pollinated: "+ str(wild_prop))               
            if ag:
                ag = ag.replace(",",";")
                ag = ag.split(';')
            suffix=arcpy.GetParameterAsText(11)
            if suffix is "#":
                suffix = ""
            parameters.append("Suffix: "+str(suffix))
            if lcrasterf:
                lcovers=["Current","Future"]
            else:
                lcovers=["Current"]
            realized = arcpy.GetParameterAsText(12)
        if realized == 'true':
            realized = True
        else: 
            realized = False

            
            
        
        
    
            
        ########################################################################################################################
        ## Functions
        ########################################################################################################################
        #function to select unique items from list for pollinator abundance on farms
        #----------------------------------------------------------------------
        #function to build list of pollinators required by given crop
        def checkpollinators(crop):
            cur=arcpy.SearchCursor(guild)
            
            pollinators=[]
            for row in cur:
                if row.getValue(crop):
                    pollinators.append(row.species)
            del cur
            return pollinators
        #----------------------------------------------------------------------
        #check projections
        def checkProjections(thedata):
            dataDesc = arcpy.Describe(thedata)
            spatreflc = dataDesc.SpatialReference
            if spatreflc.Type <> 'Projected':
                arcpy.AddError(thedata +" Your data : "+thedata+" is not projected.  You must use a projected raster and projection must be defined.")
            elif spatreflc.LinearUnitName <> 'Meter':
                arcpy.AddWarning("This model assumes that "+thedata+" is projected in meters.  You may get erroneous results")
            if thedata<>lcraster:
                if arcpy.Describe(lcraster).SpatialReference.name<>spatreflc.name and spatreflc.Type == 'Projected':
                    arcpy.AddError("WARNING: "+ thedata +" is in a different projection from the land cover data.  You may get erroneous results.")
        #----------------------------------------------------------------------    
        # function to toggle from main folder to intermediate folder for reading and writing data
        down="down"
        def switchws10(direction="up"):
            if direction=="down":
                arcpy.env.workspace=original_workspace+"\\intermediate"
            else:
                arcpy.env.workspace=original_workspace 
                
        #check and create folders
        invest10.createfolders(arcpy.env.workspace)

        
        ########################################################################################################################
        ## Process each landcover given in turn
        ## Steps
        ## 1. 
        ########################################################################################################################
        
        for cover in lcovers:
            #reset suffix
            suffix=""
            ##set suffix based on the cover
            if cover=="Future":
                lcraster=lcrasterf
                suffix = "_fut"+suffix
            else:
                suffix = "_cur"+suffix
                
            arcpy.AddMessage("Processing: " + cover + " landcover")
            
            
            #----------------------------------------------------------------------  
            # Get raster properties of the cover
            rasterDesc=arcpy.Describe(lcraster)
            rasterExtent = rasterDesc.extent
            rasterSpatref = rasterDesc.SpatialReference        
            CellSize = rasterDesc.MeanCellHeight
            arcpy.env.extent=rasterExtent
            lcwidth = rasterDesc.width
            lcheight = rasterDesc.height
            arcpy.env.mask = lcraster
            arcpy.env.cellSize=lcraster
            
            if resolution:
                arcpy.env.cellSize=int(resolution)
            checkProjections(lcraster)
            if lcrasterf:
                checkProjections(lcrasterf)
            
            switchws10(down)
            #----------------------------------------------------------------------  
            #create agricultural areas map
     
            outConstRaster = CreateConstantRaster(1, "INTEGER", arcpy.env.cellSize, arcpy.env.extent)
            arcpy.DefineProjection_management(outConstRaster, rasterSpatref)
            unityraster = outConstRaster
            outConstRaster.save("unity")
    
            cleanuplist.append("unity")
            
            arcpy.AddMessage("Creating agricultural areas map...")
    
            if ag:
                
                agfile = open(arcpy.env.workspace+"\\ag.asc","w")
                cur = arcpy.SearchCursor(lutable)
    
                for row in cur:
                    agfile.write(str(int(row.LULC))),
                    agfile.write(":"),
                    if str(int(row.LULC)) in ag:
                        agfile.write("1\n")
                    else:
                        agfile.write("0\n")
                agfile.close()
                #create an agricultural areas map using the text file created above
                filecheck = invest10.checkASCIIFile(arcpy.env.workspace+"\\ag.asc")
                if filecheck.hasDuplicates():
                    raise Exception, "Your landuse classes file has duplicates"
                if not filecheck.isSorted():
                    raise Exception, "Your landuse file is not sorted"
                #gp.ReclassByASCIIFile_sa(lcraster, arcpy.env.workspace+"\\ag.asc", "agmap", "NODATA")
                outRaster = ReclassByASCIIFile(lcraster, arcpy.env.workspace+"\\ag.asc", "NODATA")
                outRaster.save("agmap")
    
            else:
                #create constant raster
                arcpy.env.extent = lcraster
                #gp.copy_management(arcpy.env.workspace+"\\unity", arcpy.env.workspace+"\\agmap")
                arcpy.gp.SingleOutputMapAlgebra_sa("unity", "agmap")
    
                
               
            ########################################################################################################################
            ## This section calculates the land cover nest suitability and creates a grid for each species
            ########################################################################################################################
            ## Initial validation check and preparation
            ########################################################################################################################
            
            # Its assumed nest type fields will be named N_Nestx where x is an index        
            # ----------------------------------------------------------------------------------------
            #create array of nests from land cover table
            NestFields = arcpy.ListFields(lutable, "N_*", "All")
            if not NestFields:
                raise Exception, "Nest data column (N_)required in land cover table"
            nnests=[]
            for field in NestFields:
                nnests.append(field.name)
    
            # ----------------------------------------------------------------------------------------
            #count number of bees/pollinatios and build list of pollinators
            guildcursor = arcpy.SearchCursor(guild)
    
            numberofguilds = 0
            pollinatorlist = []
            for beerow in guildcursor:
                numberofguilds += 1
                pollinatorlist.append(str(beerow.Species))
                
            # ----------------------------------------------------------------------------------------
            #create array of nests from guild table
            NestFields = arcpy.ListFields(guild, "NS_*", "All")
            if not NestFields:
                raise Exception, "Nest data column (NS_) required in guild table"
            nsnests=[]
            for field in NestFields:
                nsnests.append(field.name)
            
            # ----------------------------------------------------------------------------------------
            #create array of seasons from land cover table
            SeasonFields = arcpy.ListFields(lutable, "F_*", "All")
            if not SeasonFields:
                raise Exception, "Season data column(F_) required in landuse table"
            lcoverseasons=[]
            for field in SeasonFields:
                lcoverseasons.append(field.name)
            
            # ----------------------------------------------------------------------------------------
            #create array of nests from guild table
            SeasonFields = arcpy.ListFields(guild, "FS_*", "All")
            if not SeasonFields:
                raise Exception, "Season data column (FS_) required"
            guildseasons=[]
            for field in SeasonFields:
                guildseasons.append(field.name)
    
            
            # ----------------------------------------------------------------------------------------    
            #ensure that the nest count and nest names match in landuse and guild tables 
            if len(nsnests)!=len(nnests):
                raise Exception, "Nest mismatch"
            else:
                for x in nnests:
                    x = x.upper()[2:]
                    if not x in [x.upper()[3:] for x in nsnests]:
                        raise Exception, "Nest mismatch"
            
            # ----------------------------------------------------------------------------------------
            #ensure that the seasons count and seasons names match in landuse and guild tables 
            if len(lcoverseasons)!=len(guildseasons):
                raise Exception, "Seasons mismatch seasons columns (FS_ and F_) must match in guilds and landuse tables"
            else:
                for x in lcoverseasons:
                    x = x.upper()[2:]
                    if not x in [x.upper()[3:] for x in guildseasons]:
                        raise Exception, "Seasons mismatch seasons columns (FS_ and F_) must match in guilds and landuse tables"
        
            #----------------------------------------------------------------------    
            # check if there are species duplicate names and if first 4 characters are unique
            if invest10.ListHasDuplicates(pollinatorlist):
                raise Exception, "Pollinator names must be unique, exiting."
            spp = [x[:4] for x in pollinatorlist]  
            if invest10.ListHasDuplicates(spp):
                    sppdup = 1
            ########################################################################################################################
            ## 
            ########################################################################################################################
            # ----------------------------------------------------------------------------------------
            #create nesting availability raster for each nest type
            for nnest in nnests:
                try:
                    nestingfile = open(arcpy.env.workspace+"\\"+nnest+"values.asc","w")
                    cleanuposlist.append(nnest+"values.asc")
                    cur = arcpy.SearchCursor(lutable)
                    for row in cur:
                        nestingfile.write(str(int(row.LULC))),
                        nestingfile.write(":"),
                        nestingfile.write(str(int(row.getValue(nnest)*scaling))+"\n")
    
                    nestingfile.close()
                    filecheck = invest10.checkASCIIFile(arcpy.env.workspace+"\\"+nnest+"values.asc")
                    if filecheck.hasDuplicates():
                        raise Exception, "Your landuse classes file has duplicates"
                    if not filecheck.isSorted():
                        raise Exception, "Your landuse file is not sorted"
                    #gp.ReclassByASCIIFile_sa(lcraster, arcpy.env.workspace+"\\"+nnest+"values.asc", nnest+"temp", "NODATA")
                    outRaster = ReclassByASCIIFile(lcraster, arcpy.env.workspace+"\\"+nnest+"values.asc", "NODATA")
                    outRaster.save(nnest+"temp")
    
                    del cur
                    
                    # scale back since it was multiplied by scaling factor specified in the configuration
                    expr = nnest+"temp / "+str(scaling)
                    arcpy.gp.SingleOutputMapAlgebra_sa(expr, nnest)
                    cleanuplist.append(nnest)
                    if cleanup:
                        if arcpy.Exists(nnest+"temp"):
                            arcpy.Delete_management(nnest+"temp")
                except:
                    raise Exception, "Error sup_ting raster"
            
            # ----------------------------------------------------------------------------------------          
            #Create nesting maps
            arcpy.AddMessage("Creating nesting maps...")
            cur=arcpy.SearchCursor(guild)
         
            for row in cur:
                maxexpr="max("
                if sppNameToText(row.Species):
                    focalspp = 's'+str(int(row.Species))
                    focalsppfullname = 'species ' + str(int(row.Species))
                    focalspp = focalspp[:4]
                else:
                    focalspp = row.Species[:4]
                    focalsppfullname = row.Species
                if sppdup:
                    focalsppfullname = row.Species
                    focalspp = row.Species[:3]+ str(row.OID+1) 
                arcpy.AddMessage("   "+str(focalsppfullname))
            
                for nsnest in nsnests:
                    expr="N_"+nsnest[3:]+" * "+str(row.getValue(nsnest))
                    arcpy.gp.SingleOutputMapAlgebra_sa(expr, focalspp+nsnest[3:])
                    cleanuplist.append(focalspp+nsnest[3:])
                    maxexpr+=focalspp+nsnest[3:]+","
                    
                maxexpr=maxexpr[:-1]
                maxexpr+=")"
                arcpy.gp.SingleOutputMapAlgebra_sa(maxexpr, "hn_"+focalspp+suffix)
    
            del cur
            
            # ----------------------------------------------------------------------------------------            
            #create floral availability raster for each season
            arcpy.AddMessage("Creating floral resources map...")
            for season in lcoverseasons:
                try:
                    seasonfile = open(arcpy.env.workspace+"\\"+season+"values.asc","w")
                    cleanuposlist.append(season+"values.asc")
                    cur = arcpy.SearchCursor(lutable)
    
                    for row in cur:
                        seasonfile.write(str(int(row.LULC))),
                        seasonfile.write(":"),
                        seasonfile.write(str(int(row.getValue(season)*scaling))+"\n")
                    seasonfile.close() 
                    filecheck = invest10.checkASCIIFile(arcpy.env.workspace+"\\"+season+"values.asc")
                    if filecheck.hasDuplicates():
                        raise Exception, "Your landuse classes file has duplicates"
                    if not filecheck.isSorted():
                        raise Exception, "Your landuse file is not sorted"
                    #gp.ReclassByASCIIFile_sa(lcraster, arcpy.env.workspace+"\\"+season+"values.asc", season+"temp", "DATA")
                    outRaster = ReclassByASCIIFile(lcraster, arcpy.env.workspace+"\\"+season+"values.asc", "DATA")
                    outRaster.save(season+"temp")
    
                    del cur
                    
                    # scale back since it was multiplied by scaling factor specified in the configuration
                    expr = season+"temp / "+str(scaling)
                    arcpy.gp.SingleOutputMapAlgebra_sa(expr, season)
                    if arcpy.Exists(season+"temp"):
                        if cleanup:
                            arcpy.Delete_management(season+"temp")
                except:
                    raise Exception, msgSeasonsRaster
                
            ########################################################################################################################
            ## process each species or guild table row by row (from the guild table) and for each season
            ########################################################################################################################          
    
            rows=arcpy.SearchCursor(guild)
            PAFarms_expression="" #Used for summung pollinator availability on farms
            Pm_expression="(" #PAFarms expression builder, used later
            countspp=0
            for row in rows: #for each row/species
                if sppNameToText(row.Species):
                    focalspp = 's'+str(int(row.Species))
                    focalsppfullname = 'species ' + str(int(row.Species))
                    focalspp = focalspp[:4]
                else:
                    focalspp = row.Species[:4]
                    focalsppfullname = row.Species
                if sppdup:
                    focalsppfullname = row.Species
                    focalspp = row.Species[:3]+ str(row.OID+1) 
                
                arcpy.AddMessage("   "+focalsppfullname+"...")
                alpha = row.alpha
                #check that species name is a string
                if sppNameToText(row.Species):
                    focalspp = 's'+str(int(row.Species))
                    focalsppfullname = 'species ' + str(int(row.Species))
                    focalspp = focalspp[:4]
                else:
                    focalspp = row.Species[:4]
                    focalsppfullname = row.Species
                if sppdup:
                    focalsppfullname = row.Species
                    focalspp = row.Species[:3]+ str(row.OID+1) 
                hfexpr=""
                for season in guildseasons: #for each season
                    fs = row.getValue("FS_"+season[3:])
                    if not fs > 0:
                        continue
                    arcpy.AddMessage("     Season: "+season[3:])
                    
                    # ----------------------------------------------------------------------------------------  
                    # Calculate neighborhood based on alpha - used to create the distance matrix
                    nb=float(numpy.ceil(alpha))*2 
                    neighborhood=int(nb)
                    neighbors =numpy.ceil(neighborhood/float(arcpy.env.cellSize))
                    if not neighbors % 2:
                        neighbors += 1
                    neighbors = int(neighbors)
                    if neighbors > 101:
                        neighbors = 101
                    if neighbors < 5:
                        neighbors = 5
                    nbcells = (neighbors - 1) / 2
        
                    # ----------------------------------------------------------------------------------------  
                    #create a dist array which has weight values for running focal statistics.  The weight is 
                    #used to describe the exponential decline of pollinator visits with distance
                    
                    distplain=[]
                    neighborstring=str(neighbors)+" "
                    #step through the neighborhood array
                    for i in range(0,neighbors):
                        innerlist=[]
                        innerlistplain=[]
                        innerlistdenom = []
                        for j in range(0,neighbors):
                            thedist = numpy.sqrt((abs(i - nbcells) * float(arcpy.env.cellSize)) ** 2 + (abs(j - nbcells) * float(arcpy.env.cellSize)) ** 2)
                            if squareneighborhood:
                                weightdistplain = round(numpy.exp(-thedist / alpha),4)
                            else: 
                                if thedist > neighborhood / 2:
                                    weightdistplain = 0
                                else:
                                    weightdistplain = round(numpy.exp(-thedist / alpha),4)
                            innerlistplain.append(weightdistplain)
                        distplain.append(innerlistplain)
                    # ----------------------------------------------------------------------------------------                     
                    # write the weight distance files to a text file for focal statistics    
                    weightfile = open(arcpy.env.workspace+"\\weightpollen"+focalspp[:4]+".txt","w")
                    cleanuplist.append("weightpollen"+focalspp+".txt")
                    
                    weightfile.write(neighborstring*2+"\n")
                   # ----------------------------------------------------------------------------------------                  
                    #calculate the sum for the weights for normalizing the weight matrix
                    sumweight = 0
                    for x in distplain:
                        for y in x:
                            sumweight += y
                   # ---------------------------------------------------------------------------------------- 
                   #write weight file
                    for x in distplain:
                        for y in x:
                            if sumweight:
                                weightfile.write(str(round(y/float(sumweight), 6))+' ')
                                #if debug: num.write(str(round(y, 6))+' ')
                            else:
                                weightfile.write('0 ')
                        weightfile.write("\n")
                        #if debug: num.write("\n")
                        
                    weightfile.close()
                     
                    # ---------------------------------------------------------------------------------------- 
                    #calculate foraging map for this species
                    #arcpy.AddMessage("  Creating foraging map (hf_*) for "+str(focalsppfullname)+"...")
                    fraster = "f_"+season[3:9]
                    cleanuplist.append("f_"+season[3:9])
                    OutRaster = arcpy.env.workspace+"\\hf"+focalspp+season[3:]
                    InNeighborhood = "WEIGHT "+arcpy.env.workspace+"\\weightpollen"+focalspp[:4]+".txt CELL"
                    InNoDataOption = "DATA"
                    # FocalStatistics
                    arcpy.gp.FocalStatistics_sa(fraster, OutRaster, InNeighborhood, "SUM", InNoDataOption)
                    cleanuplist.append("hf"+focalspp+season[3:])
                    
                    #if there are multiple season take the minimum value
                    
                    if len(guildseasons)>1:#if more than one season
                        hfexpr+="hf"+focalspp+season[3:]+" * " + str(fs) + " + "
                    else:
                        hfexpr=""
                if hfexpr:
                    hfexpr = hfexpr[:-3]
                    
                    arcpy.gp.SingleOutputMapAlgebra_sa(hfexpr, arcpy.env.workspace+"\\hf_"+focalspp+suffix) #take the minimum of the seasons
                else:
                    arcpy.gp.SingleOutputMapAlgebra_sa("hf"+focalspp+guildseasons[0][3:], arcpy.env.workspace+"\\hf_"+focalspp+suffix) #just one season so take the only value
                
                # ---------------------------------------------------------------------------------------- 
                #calculate the pollinator abundance for this species (nesting map - determined by nest availability and floral resources availability within range 
                #of this pollinator flight distance - alpha)
                #arcpy.AddMessage("  Creating nesting map (sup_*) for "+str(focalsppfullname)+"...")
                switchws10()
                pexpr = arcpy.env.workspace+"\\intermediate\\hf_"+focalspp+suffix+" * " +arcpy.env.workspace+"\\intermediate\\hn_"+focalspp+suffix
                arcpy.gp.SingleOutputMapAlgebra_sa(pexpr, "\\intermediate\\sup_"+focalspp[:4]+suffix)
                
                # ---------------------------------------------------------------------------------------- 
                #Calculate pollinator abundance on each cell abundance of each pollinator
                inraster = arcpy.env.workspace+"\\intermediate\\sup_"+focalspp[:4]+suffix
                OutRaster = arcpy.env.workspace+"\\intermediate\\frmt_"+focalspp[:4]
                InNeighborhood = "WEIGHT "+arcpy.env.workspace+"\\intermediate\\weightpollen"+focalspp[:4]+".txt CELL"
                InNoDataOption = "DATA"
                cleanuplist.append("frmt_"+focalspp[:4])
        
                # FocalStatistics
                #gp.FocalStatistics_sa(inraster, OutRaster, InNeighborhood, "SUM", InNoDataOption)  
                outFocalStatistics = FocalStatistics(inraster, InNeighborhood, "SUM",InNoDataOption)
                outFocalStatistics.save(OutRaster)
    
                switchws10(down)
                
                #create mask
                arcpy.gp.SingleOutputMapAlgebra_sa("1 - agmap", "aginvert")
                
                #remove non ag cells
                outSetNull = SetNull("aginvert", "frmt_"+focalspp[:4], "")
                outSetNull.save("frm_"+focalspp[:4])
    
                #gp.SetNull_sa ("aginvert", "frmt_"+focalspp[:4], "frm_"+focalspp[:4], "") 
                
                #create an expression for abundance for all species        
                PAFarms_expression+= "frm_"+focalspp[:4]+" + "
    
                #increment species counter
                countspp+=1
        
                #update expression for calculating relative pollinator score
                switchws10()
                Pm_expression+= arcpy.env.workspace+"\\intermediate\\sup_"+focalspp[:4]+suffix+" + "
                switchws10(down)
                
            del rows
            
            # ---------------------------------------------------------------------------------------- 
            #complete PAFarms expression builder (pollinator abundance on farms)
            PAFarms_expression = PAFarms_expression[:-3]
        
            # ---------------------------------------------------------------------------------------- 
            #complete Pm expression builder and calculate Pm
            Pm_expression = Pm_expression[:-3]+") / "+str(countspp)
            # ---------------------------------------------------------------------------------------- 
            #compute the total nesting map for all species/guilds
            switchws10()
            arcpy.gp.SingleOutputMapAlgebra_sa(Pm_expression, "\\output\\sup_tot"+suffix)
    
            switchws10(down)
    
            arcpy.AddMessage("Creating map of pollinator abundance...")
            # ---------------------------------------------------------------------------------------- 
            # Calculate pollinatior abundance on farms for all species/guilds (PAFarms)
            arcpy.gp.SingleOutputMapAlgebra_sa(PAFarms_expression, "frm_tot")
            switchws10()
            # ---------------------------------------------------------------------------------------- 
            #calculate average abundance (first statement muted to test with non ag pixels not removed)
            arcpy.gp.SingleOutputMapAlgebra_sa(arcpy.env.workspace+"\\intermediate\\frm_tot / "+str(numberofguilds), arcpy.env.workspace+"\\output\\frm_avg"+suffix)
            if realized:
                # ---------------------------------------------------------------------------------------- 
                #create value map
                arcpy.AddMessage("Creating value map...")
                Vm_expression =  "(1 - "+str(wild_prop)+ ") + ("+str(wild_prop) +" * ("+arcpy.env.workspace+"\\output\\frm_avg"+suffix+" / ("+sat_const+" + "+arcpy.env.workspace+"\\output\\frm_avg"+suffix+")))"
                #Vm_expression =  arcpy.env.workspace+"\\output\\frm_avg"+suffix+" / ("+sat_const+" + "+arcpy.env.workspace+"\\output\\frm_avg"+suffix+")"
                arcpy.gp.SingleOutputMapAlgebra_sa(Vm_expression, arcpy.env.workspace+"\\intermediate\\frm_val"+suffix)
        
                ########################################################################################################################
                ## Calculate the pollinator service value - this part uses gdal and numpy libraries
                ########################################################################################################################
            
                #read the projection info of the value raster
                #ds = gdal.Open(arcpy.env.workspace+"\\intermediate\\frm_val"+suffix, GA_ReadOnly)   
                #geoTransform = ds.GetGeoTransform()
                #proj = ds.GetProjection()
                #ds = None
                #process each species or guild table row by row
                rows=arcpy.SearchCursor(guild)
                
                switchws10(down)
                arcpy.AddMessage("Calculating pollinator service value...")
                speciesnames = []
                for row in rows:
                    if sppNameToText(row.Species):
                        focalspp = 's'+str(int(row.Species))
                        focalsppfullname = 'species ' + str(int(row.Species))
                        focalspp = focalspp[:4]
                    else:
                        focalspp = row.Species[:4]
                        focalsppfullname = row.Species
                    if sppdup:
                        focalsppfullname = row.Species
                        focalspp = row.Species[:3]+ str(row.OID+1)    
                    arcpy.AddMessage("   "+focalsppfullname)
                    
                    #get the contribution of this species
                    switchws10()
                    contexpr = arcpy.env.workspace+"\\intermediate\\frm_val"+suffix+" * "+arcpy.env.workspace+"\\intermediate\\frm_"+focalspp[:4]+" / "+arcpy.env.workspace+"\\intermediate\\frm_tot"
                    arcpy.gp.SingleOutputMapAlgebra_sa(contexpr, arcpy.env.workspace+"\\intermediate\\ctr_"+focalspp[:4]+suffix)
                    cleanuplist.append("ctr_"+str(focalspp[:4])+suffix)
                    
                    switchws10(down)
                    val_raster = "ctr_"+focalspp[:4]+suffix
                    sup_raster = "sup_"+focalspp[:4]+suffix
                    frm_raster = "frm_"+focalspp[:4]
                    crop_raster = "agmap"
                    distfile = arcpy.env.workspace+"\\weightpollen"+focalspp[:4]+".txt"
                    arcpy.env.extent = lcraster
                    speciesnames.append("sp_"+focalspp[:4]+suffix)
        
                    workingdir = arcpy.env.workspace
                    #matrix = arcpy.RasterToNumPyArray(distfile)
                    matrix = invest10.readFileIntoArray(distfile, 1)
                    #nst =  invest.convertRastertoArray(workingdir, sup_raster)
                    nst = arcpy.RasterToNumPyArray(sup_raster)
                    #poll_abund =  invest.convertRastertoArray(workingdir, frm_raster)
                    poll_abund = arcpy.RasterToNumPyArray(frm_raster)
                    #value =  invest.convertRastertoArray(workingdir, val_raster)
                    value = arcpy.RasterToNumPyArray(val_raster)
                    #agmap =  invest.convertRastertoArray(workingdir, crop_raster)
                    agmap = arcpy.RasterToNumPyArray(crop_raster)
                    output = zeros(shape=value.shape)
                    outputsum = zeros(shape=value.shape)
                    
                    height = value.shape[0]
                    width = value.shape[1]
                    
                    sizematrix = matrix.shape[0]
                    nb = int(sizematrix / 2)
                    try:
                        for j in arange(height): #process each row from the top working downwards
                            for i in arange(width): #process each column from left to right
                                if agmap[j][i] > 0: #only process if its a cropped pixel, value > 0
                                    if i > nb-1 and i < (width - nb) and j > nb-1 and j < (height - nb): #'ensure the block does not go beyond raster
                                        for jstart in arange(sizematrix): # 0 To nb - 1 'loop through the neighborhood
                                            for istart in arange(sizematrix): #= 0 To nb - 1
                                                if poll_abund[j][i] > 0: # 'this must be greater than 0 to prevent overflow in Next line. If its zero then set the result to zero
                                                    result = value[j][i] * matrix[jstart][istart] * nst[j - nb + jstart][i - nb + istart] / poll_abund[j][i]
                                                else: 
                                                    result = 0
                                                temphold = round(output[j - nb + jstart][i - nb + istart] + result, 6)
                                                if temphold < 0:
                                                    temphold = 0
                                                if temphold > 3.402823466e+38:
                                                    temphold = 0                                        
                                                output[j - nb + jstart][i - nb + istart] = temphold 
                                                
                        #invest.convertArraytoRaster(output, "sp_"+focalspp[:4]+suffix+".img", workingdir, geoTransform, proj)
                        point = arcpy.Point(rasterDesc.extent.XMin, rasterDesc.extent.YMin)
                        arcpy.env.mask = lcraster
                        newRaster = arcpy.NumPyArrayToRaster(output, point, int(arcpy.env.cellSize), int(arcpy.env.cellSize), 255)
                        newRaster = newRaster * unityraster
                        arcpy.DefineProjection_management(newRaster, rasterSpatref)
                        newRaster.save("sp_"+focalspp[:4]+suffix)
                    except:
                        raise Exception, "Error calculating realized supply for species: "+focalsppfullname
    
                # ----------------------------------------------------------------------------------------     
                #calculate total realized supply
                arcpy.AddMessage("  Combining pollinator service value maps")
                switchws10()
                exprSupply = ''
                for x in speciesnames:
                    exprSupply += arcpy.env.workspace+"\\intermediate\\"+x+ " + "
                exprSupply = exprSupply[:-3]
                arcpy.env.extent = lcraster
                arcpy.env.cellSize = lcraster
                outConstant = CreateConstantRaster(0)
                arcpy.DefineProjection_management(outConstant, rasterSpatref)
                outConstant.save("C:/InVEST1005/pollination/constant")
                for sppoutput in speciesnames:
                    outConstant = outConstant + arcpy.Raster(arcpy.env.workspace+os.sep+'intermediate'+os.sep+sppoutput)
                outConstant.save(arcpy.env.workspace+os.sep+'output'+os.sep+"sup_val"+suffix)
                #arcpy.gp.SingleOutputMapAlgebra_sa(exprSupply, arcpy.env.workspace+"\\output\\sup_val"+suffix)
                # ---------------------------------------------------------------------------------------- 
            #write parameters
            switchws10()
            #if suffix:
                #suffix = "_"+suffix
            switchws10()
            parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
            
            parafile = open(arcpy.env.workspace+"\\output\\pollination_"+now.strftime("%Y-%m-%d-%H-%M")+suffix+".txt","w") 
            parafile.writelines("POLLINATION MODEL PARAMETERS\n")
            parafile.writelines("____________________________\n\n")
            
            for para in parameters:
                parafile.writelines(para+"\n")
                parafile.writelines("\n")
            parafile.close()
            switchws10(down)
        #check if name speciesnames is defined.  
        try:
            speciesnames
        except NameError:
            speciesnames = []
            
        if cleanup:
            arcpy.AddMessage("Cleaning up...")
            for theraster in cleanuplist+["PAFarmsDenom","PAFarmsNum","unityraster","numerator","denominator","ag.asc","tempo", "agmap", "aginvert"]+speciesnames:
                if arcpy.Exists(theraster):
                    arcpy.Delete_management(theraster)
            #delete non datasets
            os.chdir(arcpy.env.workspace)
            for thefile in cleanuposlist+["ag.asc"]:
                if os.path.exists(thefile):
                    os.remove(thefile)
        
    except Exception, ErrorDesc:
        arcpy.AddError(str(ErrorDesc))
