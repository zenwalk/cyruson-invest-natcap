# ---------------------------------------------------------------------------
# biodiversity.py
# Coded by: Nasser Olwero (nasser.olwero@wwfus.org)
# Contributor: Erik Nelson (nels1069@umn.edu)
# Contributor: Derric Pennington (penn0107@umn.edu)
# Created on: Wed June 08 2008 revision Friday December 11 2009
# Update: 01/21/2010
# Description: Computes habitat quality and rarity
# ---------------------------------------------------------------------------

# Import system modules 
import sys, string, os, arcgisscripting, math, datetime, re
now = datetime.datetime.now()

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out any necessary licenses
gp.CheckOutExtension("spatial")

#Set output handling
gp.OverwriteOutput = 1

msgCountThreats = "Error enumerating threats"
msgThreatMismatch = "Threats are mismatched with sensitivity values"
msgRasterDescription = "Error reading raster properties"
msgDistArray = "Error creating distance array"
msgWeightFile = "Error creating weight file"
msgThreatNeighborhood = "Error computing proximity threat values for "
msgThreatProcessing = "Error processing "
msgReclass = "Error creating habitat cover layer"
msgSens = "Sensitivity fields missing"
msgAccess = "Error creating access layer"
msgArguments = "Problem with arguments."
msgDistArray = "Error creating distance array"
msgWeightFile = "Error creating weight file"
msgThreatProcessing = "Error processing "
msgRarityAscii = "Error creating rarity ASCII"
msgRarity = "Error calculating rarity"
msgCleanup = "Error cleaning up(deleting data)"
msgQuality = "Error adjusting quality by habitat"
msgDegradation = "Error calculating degradation for all threats"
msgQualityfromDeg = "Error creating quality layer from degradation"

########################################################################################################################
## Configuration
########################################################################################################################
# Sets the level of reporting on the dialog. True=print reports, False=Do not print
verbose = True
# Sets the option to cleanup intermediate files. True=delete temporary files, False=Do not delete temporary files
cleanup = True
# Toggles debug mode.  In debug mode the parameters are set in this file and model can be run from python.
debug = 0    #initialize cleanup list.  All files to be cleaned after run are added here
cleanuplist=[]
#initialize os cleanup list.  All files to be cleaned after run are added here
cleanuposlist=[]    


try:
    ##############################################################################################################################
    # read in arguments
    # Local variables...
    try:
        parameters = []
        if debug:
            pass
        else:
            parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
            gp.Workspace = gp.GetParameterAsText(0)
            parameters.append("Workspace: "+gp.Workspace) 
            original_workspace = gp.Workspace
            
            lcraster = gp.GetParameterAsText(1)
            parameters.append("Current Landcover: "+ lcraster)
            
            lcrasterf = gp.GetParameterAsText(2)
            parameters.append("Future Landcover: "+ lcrasterf)
                 
            lcrasterbase = gp.GetParameterAsText(3)
            parameters.append("Baseline Landcover: "+ lcrasterbase)
    
            threats_data = gp.GetParameterAsText(4)
            parameters.append("Threats data table: "+ threats_data)
            
            access_lyr = gp.GetParameterAsText(5)
            parameters.append("Accessibility table: "+ access_lyr)
            
            landcover_data = gp.GetParameterAsText(6)
            parameters.append("L: "+ access_lyr)
            
            k = gp.GetParameterAsText(7)
            parameters.append("Half Saturation Constant: "+ k)
            k = int(k)
            
            resolution = gp.GetParameterAsText(8)
            
            parameters.append("Resolution: "+ str(resolution))
            
            if resolution:
                resolution = int(gp.GetParameterAsText(8))
              
            parameters.append("Future Landcover: "+ lcrasterf)
           
            suffix=gp.GetParameterAsText(9)
            
            if suffix is "#":
                suffix = ""
            parameters.append("Suffix: "+ suffix)
            if len(suffix) > 1:
                msgArguments += " Suffix too long, you used '"+suffix+"'. please use just one character"
                raise
     
    except:
        raise Exception, msgArguments
    
    ##############################################################################################################################
    # Functions
    #check fields
    def checkfields(fields, target):
        thefields = gp.listfields(target,"*","All")
        field = thefields.next()
        foundfields = []
        while field:
            foundfields.append(field.Name.upper())
            field = thefields.next()
        for x in fields:
            x = x.upper()
            if not x in foundfields:
                raise Exception, "Expected field:"+x + " in " + target+". Exiting." 
    checkfields(["threat","decay","weight","max_dist"], threats_data)
    
    checkfields(["LULC","habitat"], landcover_data)
    
  
    
    #a function that checks if a file exists
    def fileexists(f):
        try:
            testfile = open(f)
        except IOError:
            exists = 0
        else:
            exists = 1
        return exists
    
    def checkProjections(thedata):
                dataDesc = gp.describe(thedata)
                spatreflc = dataDesc.SpatialReference
                if spatreflc.Type <> 'Projected':
                    gp.AddError(thedata +" does not appear to be projected.  It is assumed to be in meters")
                elif spatreflc.LinearUnitName <> 'Meter':
                    gp.AddError("This model assumes that "+thedata+" is projected in meters.  You may get erroneous results")
                if thedata<>lcraster:
                    if gp.describe(lcraster).SpatialReference.name<>spatreflc.name and spatreflc.Type == 'Projected':
                        gp.AddError("WARNING: "+ thedata +" is in a different projection from the land cover data.  You may get erroneous results.")
    
    def switchws(direction="up"):
        if debug:
            if direction=="down":
                gp.Workspace=r"C:\InVEST20\Biodiversity"+"\\intermediate"
            else:
                gp.Workspace=r"C:\InVEST20\Biodiversity"
        
        else:    
            if direction=="down":
                gp.Workspace=gp.GetParameterAsText(0)+"\\intermediate"
            else:
                gp.Workspace=gp.GetParameterAsText(0)
            
    def reclassfn(table,field,outraster,reclassraster=lcraster,scale=0):
        #this function creates a raster from table values.  However note that the values are inverted.  Made for a specific purpose
        try:
            tempfile = open(gp.Workspace+"\\reclassfile.asc","w")
            cleanuposlist.append(gp.Workspace+"\\reclassfile.asc")
            cur = gp.SearchCursor(table)
            row = cur.next()
            values=[]
            lulc=[]
            while row:
                lulc.append(int(row.LULC))
                values.append(row.GetValue(field))
                row = cur.next()
            if scale>0:
                values=[int(x*scale) for x in values]

            for x in range(len(lulc)):
                tempfile.write(str(lulc[x])),
                tempfile.write(":"),
                tempfile.write(str(values[x])+"\n")
            tempfile.close()
            
            del cur
            sortreclassfile(gp.Workspace+"\\reclassfile.asc")
            gp.ReclassByASCIIFile_sa(reclassraster, gp.Workspace+"\\reclassfile.asc", outraster, "DATA")
            #cleanuplist.append(outraster)
        except:
            raise Exception, msgReclass +" "+ field
        #function for calculating totals
    def summary(theraster):
        summaryfile = original_workspace+os.sep+"intermediate"+os.sep+"summary.dbf"
        gp.cellSize = theraster
        gp.extent = lcraster
        gp.mask = lcraster
        gp.CreateConstantRaster_sa(original_workspace+os.sep+"intermediate"+os.sep+"unitynew", 1, "INTEGER")
        gp.ZonalStatisticsAsTable_sa (original_workspace+os.sep+"intermediate"+os.sep+"unitynew", "Value", theraster, summaryfile, "DATA")
        cur = gp.SearchCursor(summaryfile)
        row = cur.Next()
        valuesum = row.SUM
        del cur
        del row
        if gp.Exists(summaryfile):
            gp.Delete_management(summaryfile)
        if gp.Exists(original_workspace+os.sep+"intermediate"+os.sep+"unitynew"):
            gp.Delete_management(original_workspace+os.sep+"intermediate"+os.sep+"unitynew")
        return valuesum         
    #function for formatting string with commas
    def thous(x):
        return re.sub(r'(\d{3})(?=\d)', r'\1,', str(x)[::-1])[::-1]    
    #function sorts a text file before its used for reclassing raster
    def sortreclassfile(thefile):
        newlines = []
        filetosort = open(thefile)
        lines = filetosort.readlines()
        for line in lines:
            splat = line.split(":")
            splat = map(int, splat)
            newlines.append(splat)
        filetosort.close()
        filetosort = open(thefile, "w")
        newlines = sorted(newlines, key=lambda cover: cover[0])
        for line in newlines:
            filetosort.write(str(line[0])+":"+str(line[1])+"\n")
        filetosort.close()    
    ##############################################################################################################################   
    quality_sum = []
    if access_lyr:
        checkfields(["access"], access_lyr)
    
    lcovers = [lcraster]
    #lcovers  = []
    if lcrasterf:
        lcovers.append(lcrasterf)
    if lcrasterbase:
        lcovers.append(lcrasterbase)


     

    for cover in lcovers:
        
        if cover==lcraster:
            period="cur"
            period1="Current"
        elif cover==lcrasterf:
            period="fut"
            period1="Future"
        elif cover==lcrasterbase:
            period="bas"
            period1="Base"        
        gp.addmessage("Processing "+period1+" landcover")
        
        ##############################################################################################################################
        # get raster properties
        try:
            rasterDesc=gp.describe(cover)
            rasterExtent = rasterDesc.extent
            CellSize = rasterDesc.MeanCellHeight

            gp.mask = cover
            width = rasterDesc.Width
            height = rasterDesc.Height
            cells = width * height
            gp.Extent = cover

            gp.CellSize=cover
            rasterDesc=gp.describe(cover)
            if resolution:
                gp.CellSize=resolution
            rasterDesc=gp.describe(cover)
            CellSize = gp.CellSize



        except:
            raise Exception, msgRasterDescription   
        ##############################################################################################################################
        down="down"
        #if cover==lcraster:

        #check and create folders
        switchws()
        try:
            thefolders=["Output","intermediate"]
            for folder in thefolders:
                if not gp.exists(gp.Workspace+"\\"+folder):
                    gp.CreateFolder_management(gp.Workspace, folder)
        except:
            raise Exception, "Error creating intermediate/output folders. ensure you have write permissions and the folders are not currently being used."
        ##############################################################################################################################
        # create a raster of habitat based on the habitat field

        try:
            switchws("down")
            reclassfn(landcover_data,"habitat",gp.Workspace+"\\habitat_1", cover, 10)
            #scale back habitat map by 10 since it was multiplied by 10 above
            gp.Times_sa(gp.Workspace+"\\habitat_1","0.1",gp.Workspace+"\\habitat")
            cleanuplist.append("habitat_1")
            switchws()
            
        except:
            raise Exception, "Error creating habitat layer"
        
        ##############################################################################################################################
        #check projections
        checkProjections(cover)

        
        ##############################################################################################################################
        switchws(down)
        
        #count threats given
        try:
            cur = gp.SearchCursor(threats_data)
            row = cur.next()
            threats = []
            threatnames=[]

            while row:
                if len(row.threat)>8:
                    msgCountThreats = "Threats names should not be longer than 8 characters. You used: " + row.threat
                    raise Exception, msgCountThreats
                threatnames.append(row.threat)
                innerlist = []
                innerlist.append(row.threat)
                innerlist.append(row.weight)
                innerlist.append(row.decay)
                innerlist.append(row.max_dist)
                
                threats.append(innerlist)
                row = cur.next()
            del cur
            #calculate the sum of w, for weighting
            sumw = 0
            for x in threats:
                sumw += float(x[1])
            #change to upper case for validation
            threatnames=[y.upper() for y in threatnames] 
        except:
            raise Exception, msgCountThreats
        # check that threatcount matches with sensitivity
        try:
            SensFields = gp.listfields(landcover_data, "L_*", "All")
            field = SensFields.next()
            count = 0
            Threats_Sens=[] #use a list to track threat names used in sensitivity table
            while field:
                count +=1
                Threats_Sens.append(field.name[2:])
                field = SensFields.next()

            SensFields.reset()
            field = SensFields.next()
            Threats_Sens=[y.upper() for y in Threats_Sens]
            if not field:
                raise Exception, msgSens
            if count <> len(threats):
                raise Exception, msgThreatMismatch
            if not threatnames.sort() == Threats_Sens.sort():
                msgThreatMismatch = "Threat names used in threats table should match the names used in sensitivity table"
                raise Exception, msgThreatMismatch
        except:
            raise Exception, msgThreatMismatch

        

        ##############################################################################################################################
        #Create access layer
        if access_lyr:
            try:
                if verbose: gp.addmessage("Creating legal access layer...")
                #process access layer
                gp.FeatureToRaster_conversion(access_lyr, "access", "access_1")
                cleanuplist.append("access_1")

                #fill up access
                expr = "Con(IsNull(access_1), 1, access_1)"
                gp.SingleOutputMapAlgebra_sa(expr, "legal_access")
                cleanuplist.append("legal_access")
            except:
                raise Exception, msgAccess

        ##############################################################################################################################
        # Process density layers
        # Process each threat
        try:
            exprdeg=""
            exprqual=""
            gp.addmessage("Processing threats...") 
            for x in threats:
                
                theThreat = x[0]
                w = x[1]
                exprdeg += theThreat+"_fin + "
                
                decay = x[2]
                max_dist = x[3]
                max_dist = max_dist*1000
                densraster = x[0]
                switchws()
                #add a flag to advance to the next landcover
                nextcover = 0
                if cover == lcraster: #processing current landcover
                    if gp.exists(gp.Workspace+"\\Input\\"+densraster):
                        pass
                    elif gp.exists(gp.Workspace+"\\Input\\"+densraster+"_c"):
                        densraster += "_c"
                    elif gp.exists(gp.Workspace+"\\Input\\"+densraster+"_c.img"):
                        densraster += "_c.img"
                    else:
                        gp.addmessage("Density raster for "+period1+ " not found, moving to next period")
                        nextcover = 1
                        break
                    
                elif cover == lcrasterf:
                    if gp.exists(gp.Workspace+"\\Input\\"+densraster+"_f"):
                        densraster += "_f"
                    elif gp.exists(gp.Workspace+"\\Input\\"+densraster+"_f.img"):
                        densraster += "_f.img"
                    else:
                        gp.addmessage("Density raster for "+period1+ " not found, moving to next period")
                        nextcover = 1
                        break
                        
                elif cover == lcrasterbase:
                    
                    if gp.exists(gp.Workspace+"\\Input\\"+densraster+"_b"):
                        densraster += "_b"
                    elif gp.exists(gp.Workspace+"\\Input\\"+densraster+"_b.img"):
                        densraster += "_b.img"
                    else:
                        gp.addmessage("Density raster for "+period1+ " not found, quality will not be calculated for base landcover")
                        
                        nextcover = 1
                        break
                
                gp.addmessage("   "+x[0])        
                switchws()
                checkProjections(gp.Workspace+"\\Input\\"+densraster)
                      
                expr = gp.Workspace+"\\Input\\"+densraster
                
                switchws(down)
                gp.SingleOutputMapAlgebra_sa(expr, theThreat+"_ini")
                cleanuplist.append(theThreat+"_ini")
                
                #original threat density is densraster
                #1. adjust by distance
                #2. adjust by weight
                #3. adjust by protection/access
                #4. adjust by sensitivity
                

                ##############################################################################################################################
                ## 1. adjust threat by distance
                # Calculate neighborhood 
                nbcells = int(max_dist / int(gp.CellSize))
                neighbors = (nbcells*2)+1

                #create a dist array which has weight values for running focal statistics
                try:
                    dist=[]
                    neighborstring=str(neighbors)+" "
                    for i in range(neighbors):
                        innerlist=[]
                        for j in range(neighbors):
                            thedist = math.sqrt((abs(i - nbcells) * int(gp.CellSize)) ** 2 + (abs(j - nbcells) * int(gp.CellSize)) ** 2)
                            if thedist > max_dist:
                                weight = 0
                            else:
                                if decay == 1:
                                    weight = (1 - ((1 / float(max_dist)) * thedist))
                                else:
                                    weight = math.exp((-2.99573/float(max_dist)) * thedist)
                            innerlist.append(weight)
                        dist.append(innerlist)
                except:
                    raise Exception, msgDistArray    

                # Create a weight file for the neighborhood
                switchws(down)
                try:
                    weightfile = open(gp.Workspace+"\\weightthreat"+theThreat+".txt","w")
                    cleanuposlist.append(gp.Workspace+"\\weightthreat"+theThreat+".txt")
                    weightfile.write(neighborstring*2+"\n")
                    for xx in dist:
                        for y in xx:
                            weightfile.write(str(y)+' ')
                        weightfile.write("\n")
                    weightfile.close()
                except:
                    weightfile.close()
                    raise Exception, msgWeightFile               
                #calculate neighborhood using weight file
                try:
                    OutRaster = theThreat+"_prox"
                    cleanuplist.append(theThreat+"_prox")
                    InNeighborhood = "WEIGHT "+gp.Workspace+"\\weightthreat"+theThreat+".txt CELL"
                    # FocalStatistics 
                    gp.FocalStatistics_sa(theThreat+"_ini", OutRaster, InNeighborhood, "SUM", "DATA")
                except:
                    raise Exception, msgThreatNeighborhood + theThreat

                ##############################################################################################################################
                ## 2. adjust threat by weight
                weight_multiplier = w/float(sumw)
                gp.SingleOutputMapAlgebra_sa(theThreat+"_prox * "+ str(weight_multiplier), theThreat+"_wt")
                
                ##############################################################################################################################
                ## 3. adjust threat by protection/access
                
                if access_lyr:
                    gp.Times_sa(theThreat+"_wt","legal_access",theThreat+"_acc")
                else:
                    gp.Times_sa(theThreat+"_wt",1,theThreat+"_acc")
                
                ##############################################################################################################################
                ## 4. adjust threat by sensitivity
                               
                # Create sensitivity raster for this threat
                sensitivityfile = open(gp.Workspace+"\\sensvalues.asc","w")
                cleanuposlist.append(gp.Workspace+"\\sensvalues.asc")
                cur = gp.SearchCursor(landcover_data)
                row = cur.next()
                sens=[]
                lulc=[]
                while row:
                    lulc.append(int(row.LULC))
                    #gp.addmessage(str(row.LULC)+" L_"+theThreat+" "+row.GetValue("L_"+theThreat))
                    sens.append(float(row.GetValue("L_"+theThreat)))
                    row = cur.next()
                if max(sens)>1:
                    sens=[int(z/float(max(sens))*100) for z in sens]
                elif max(sens)<=1:
                    sens=[int(z*100) for z in sens]
                for z in range(len(lulc)):
                    sensitivityfile.write(str(lulc[z])),
                    sensitivityfile.write(":"),
                    sensitivityfile.write(str(sens[z])+"\n")
                sensitivityfile.close()
                sortreclassfile(gp.Workspace+"\\sensvalues.asc")
                gp.ReclassByASCIIFile_sa(cover, gp.Workspace+"\\sensvalues.asc", "senstemp", "DATA")
                cleanuplist.append("senstemp")
                del cur

                # scale back sensitivity by 0.01 since it was multiplied by 100 above.
                gp.Times_sa("senstemp","0.01","sens_"+theThreat)
                
                cleanuplist.append("sens_"+theThreat)
                gp.SingleOutputMapAlgebra_sa(theThreat+"_acc * "+ "sens_"+theThreat, theThreat+"_fin")

                
            exprdeg = exprdeg[:-3]
            
        except:
            raise Exception, msgThreatProcessing + theThreat +" "+ gp.getmessages()
        
        if nextcover: # the current cover had no density raster so advance to next raster
            nextcover = 0 # reset the flag
            continue
        # Compute degradation for all threats
        gp.addmessage("Calculating degradation...")
        try:
            gp.SingleOutputMapAlgebra_sa(exprdeg, "degradation")
            cleanuplist.append("degradation")
            switchws()
            gp.copy_management(gp.Workspace+"\\intermediate\\degradation", gp.Workspace+"\\Output\\degrad_"+period+suffix)
            
        except:
            raise Exception, msgDegradation

        # Compute quality for all threats
        gp.addmessage("Calculating quality...")
        switchws(down)
        try:
            z = 2.5
            
            gp.Power_sa("degradation", z, "degsq")
            
            ksq = k ** z
            gp.SingleOutputMapAlgebra_sa("1 - (degsq / (degsq + "+str(ksq)+"))", "quality_1")
            
            
            #gp.SingleOutputMapAlgebra_sa("100 - degradation", "quality_1")
            switchws(down)
        except:
            raise Exception, msgQualityfromDeg
        
        # Adjust quality by habitat status
        
        try:
            gp.SingleOutputMapAlgebra_sa("quality_1 * habitat", "quality")
            cleanuplist.append("quality")
            cleanuplist.append("quality_1")
            switchws()
            gp.copy_management(gp.Workspace+"\\intermediate\\quality", gp.Workspace+"\\Output\\qual_"+period+suffix)
            quality_sum.append([period1, summary(gp.Workspace+"\\Output\\qual_"+period+suffix)])
            switchws(down)
            
            
        except:
            raise Exception, msgQuality
        
        ## ############################################################################################################
        ## This part computes rarity
        ## ############################################################################################################
        if lcrasterbase and lcrasterbase <> cover: #if user has supplied baseline landcover, calculate rarity
            gp.addmessage("Calculating rarity...")
            #intersect cover and base to get the common area
            try:
                gp.WeightedSum_sa (cover + "; "+ lcrasterbase, "weightsum")
                gp.ExtractByMask_sa (cover, "weightsum", "newcover")
                cover = gp.Workspace+"\\newcover"
                gp.ExtractByMask_sa (lcrasterbase, "weightsum", "newbase")
                lcrasterbase = gp.Workspace+"\\newbase"
            except:
                raise Exception, "Error intersecting cover with base"            
            
            #get the cell sizes
            lcDesc=gp.describe(cover)
            covercellsize = lcDesc.MeanCellHeight
            lcbaseDesc=gp.describe(lcrasterbase)
            basecellsize = lcbaseDesc.MeanCellHeight
            cellratio = covercellsize/float(basecellsize)
    
            #check fields
            def checkfields(fields, target):
                thefields = gp.listfields(target,"*","All")
                field = thefields.next()
                foundfields = []
                while field:
                    foundfields.append(field.Name)
                    field = thefields.next()
                for x in fields:
                    if not x in foundfields:
                        raise Exception, "Expected field:"+x + " in " + target+". Exiting." 
            checkfields(["COUNT","VALUE"], cover)
            checkfields(["COUNT","VALUE"], lcrasterbase)
            def getcover(num,thearray):
                thecount = 0
                for x in thearray:
                    if x[0]==num:
                        thecount = x[1]
                        break
                return thecount

            #read baseline landcover
            try:
                cur = gp.SearchCursor(lcrasterbase)
                row = cur.next()
                basevalues = []
                
                while row:
                    innerlist = []
                    innerlist.append(row.VALUE)
                    innerlist.append(row.COUNT)
                    basevalues.append(innerlist)
                    
                    row = cur.next()
                del cur
            except:
                raise Exception, "Error reading baseline landcover"
            try:
                cur = gp.SearchCursor(cover)
                row = cur.next()
                lulcrarity = []
                while row:
                    basecover = getcover(row.VALUE, basevalues)
                    #if the count of the cells in baseline returned 0, use the count on the current/projected cover
                    if basecover == 0:
                        basecover = row.COUNT
                    innerlist = []
                    innerlist.append(row.VALUE)
                    coverratio=row.COUNT/float(basecover)
                    if coverratio>1: coverratio=1
                    innerlist.append(1 - (coverratio)*cellratio)
                    lulcrarity.append(innerlist)
                    row = cur.next()
                del cur
            except:
                raise Exception, msgCountThreats

            try:
                # create rarity ascii
                rarityfile = open(gp.Workspace+"\\rarity.txt","w")
                cleanuposlist.append(gp.Workspace+"\\rarity.txt")
                for x in lulcrarity:
                    rarityfile.write(str(x[0])),
                    rarityfile.write(":"),
                    rarityfile.write(str(int(round(x[1],4)*10000))+"\n")
                rarityfile.close()
            except:
                raise Exception, msgRarityAscii
            try:
                switchws()
                sortreclassfile(gp.Workspace+"\\intermediate\\rarity.txt")
                gp.ReclassByASCIIFile_sa(cover, gp.Workspace+"\\intermediate\\rarity.txt", gp.Workspace+"\\intermediate\\rarity_1", "DATA")
                gp.Times_sa(gp.Workspace+"\\intermediate\\rarity_1", "0.0001", gp.Workspace+"\\Output\\rarity_"+period+suffix)
                #reset base raster
                if debug:
                    lcrasterbase = r"C:\InVEST\Base_Data\Terrestrial\base_samp"
                else:
                    lcrasterbase = gp.GetParameterAsText(3)
            except:
                raise Exception, msgRarity
            
    
    #write parameters
    switchws()
    if suffix:
        suffix = "_"+suffix
    switchws()
    for qual in quality_sum:
        parameters.append("Quality Sum "+qual[0]+": "+str(thous(round(qual[1],2))))
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    
    parafile = open(gp.Workspace+"\\Output\\Biodiversity_"+now.strftime("%Y-%m-%d-%H-%M")+suffix+".txt","w") 
    parafile.writelines("BIODIVERSITY MODEL PARAMETERS\n")
    parafile.writelines("______________________________\n\n")
    
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()
    
    switchws(down)
    #Final cleanup
    if cleanup:
        try:
            gp.addmessage("Cleaning up...")
            #delete datasets
            for theraster in cleanuplist:
                if gp.exists(theraster):
                    gp.delete_management(theraster)
            #delete non datasets
            os.chdir(gp.Workspace)
            for thefile in cleanuposlist:
                if os.path.exists(thefile):
                    os.remove(thefile)
        except:
            raise Exception, msgCleanup+ gp.getmessages(2)

except Exception, ErrorDesc:
    gp.AddError(str(ErrorDesc))    
