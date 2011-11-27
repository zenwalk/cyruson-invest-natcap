# ---------------------------------------------------------------------------
# Author: Nasser Olwero
# Date: May 21, 2008
# Description: This model calculates standing stock and sequestration or Carbon
# ---------------------------------------------------------------------------
# changelog
# 6/15/08 Modified to take care of null values and resolution
# 11/10/08 When not doing HWP, no need for nullraster, unity raster
# 11/19/08 Current raster used for both current and future giving sequestration as 0 every time
# 12/1/08 valuation of sequestered carbon was reporting a wrong value
# 12/1/08  if discount rate is set as 0, crashes
# 12/4/08 Added calculation of volume and biomass of harvested wood
# 12/4/08 Resampled image is now included if a resolution values is set. These are in the intermediate folder
# 12/11/08 Biomass volume and harvested wood products was not getting calculated correctly due to a problem with counting years when there is harvest in the future 
# 02/03/2011 Summary of carbon storage and sequestration added

# Import system modules
import sys, string, os, arcgisscripting, math, time, datetime, re

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out any necessary licenses
gp.CheckOutExtension("spatial")

#Set output handling
gp.OverwriteOutput = 1

#Configuration
verbose = True
cleanup = True
debug = False

msgArguments = "Problem with arguments."
msgRasterDescription = "Error reading raster properties"
msgPressure = "Error calculating pressure layer"
msgProximity = "Error calculating proximity"
msgReclass = "Error reclassing "
msgNormalizeProximity = "Error normalizing proximity"
msgAccess = "Error creating access layer"
msgFinalCalc = "Error calculating final product"
msgWeight = "Error creating weight file"
msgMaxdist = "Error reading the max distance"
msgReclassCurrent = "Error looking up current C values"
msgReclassFuture = "Error looking up scenario C values"
msgCompute4PoolsC = "Error summing current C values"
msgCompute4PoolsF = "Error summing scenario C values"
msgHarvest = "Error computing harvest rate and time"
msgFrequency = "Error creating frequency raster"
msgHeader = "Error reading in raster header"
msgHWP = "Error calculating C in harvested wood products"
msgHWPF = "Error calculating C in harvested wood products scenario"
msgTotalCurrentC = "Error calculating total current C"
msgTotalCFutureC = "Error calcualting total scenario C"
msgComputehalflifeC = "Error computing current half life"
msgComputehalflifeF = "Error computing scenario half life"
msgEcon = "Error computing economic value"
msgCleanup = "Error cleaning up"
msgTC = "Error computing time since harvest started for current landcover"
msgSequestration = "Error calculating sequestration"
msgProjection = "Your landcover data is not projected.  Exiting"
msgComputeTF = "Error computing time for scenario HWP"
msgFldType = "Field should be of integer type"
msgHLC = "Error creating half life raster"

#initialize summation terms
total_current_carbon = None
total_future_carbon = None
total_sequest_carbon = None
total_value_sequestered = None
total_value_stored_carbon = None


try:
     # Get parameters
     try:

          parameters = []
          now = datetime.datetime.now()
          parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
          
          gp.workspace = gp.GetParameterAsText(0)
          parameters.append("Workspace: "+gp.workspace) 
          originalws = gp.workspace
          scratchws = originalws+os.sep+"scratch"
          
          
          lcraster = gp.GetParameterAsText(1)
          parameters.append("Current Landcover: "+ lcraster)
          
          CDate= int(gp.GetParameterAsText(2))
          parameters.append("Current Landcover Date: "+ str(CDate))
          
          resolution = gp.GetParameterAsText(3)
          parameters.append("Resolution: "+ str(resolution))
          
          if resolution:
               resolution = int(gp.GetParameterAsText(3))
          carbon_data = gp.GetParameterAsText(4)
          parameters.append("Carbon pools data: "+carbon_data)
          
          harv_rate_current = gp.GetParameterAsText(5)
          parameters.append("Current harvest rate map: "+ harv_rate_current)
          
          lcrasterf = gp.GetParameterAsText(6)
          parameters.append("Future Landcover: "+lcrasterf)
          
          FDate = gp.GetParameterAsText(7)
          parameters.append("Future Landcover Date: "+ str(FDate))
          
          if FDate:
               FDate = int(gp.GetParameterAsText(7))
          harv_rate_future = gp.GetParameterAsText(8)
          parameters.append("Future harvest rate map: "+harv_rate_future)
          
          econ = gp.GetParameterAsText(9)
          parameters.append("Compute economic valuation: "+econ)
          
          
          p = gp.GetParameterAsText(10)
          parameters.append("Price of Carbon per metric tonne: "+str(p))
          
          if p:
               p=int(p)
          c = gp.GetParameterAsText(11)
          parameters.append("Carbon discount rate: "+str(c))
          
          if c:
               c = int(c)
               c = c / 100.0
          else:
               c = 0
               
          r = gp.GetParameterAsText(12)
          parameters.append("Market discount rate: "+ str(r))
          
          if r:
               r = int(r)
               r = r / 100.0
          else:
               r = 0

          suffix = gp.GetParameterAsText(13)
          

          if suffix is "#":
               suffix = ""
          if len(suffix)>1:
               msgArguments += " Suffix can only be 1 character."
               raise Exception, msgArguments

          if lcrasterf:
               if not FDate:
                    msgArguments += " Year of future land cover required."
                    raise Exception, msgArguments
          if econ=='true':
               if not p:
                    msgArguments += " Price of carbon required."
                    raise Exception, msgArguments
               if not c:
                    c=0
                    #msgArguments += " Carbon discount rate required."
                    #raise Exception, msgArguments
               if not r:
                    r=0
                    #msgArguments += " Market discount rate required." 
                    #raise Exception, msgArguments

     except:
          raise Exception, msgArguments + gp.GetMessages(2)
     #check and create folders
     try:
          thefolders=["Output","intermediate","scratch"]
          for folder in thefolders:
               if not gp.exists(gp.workspace+folder):
                    gp.CreateFolder_management(gp.workspace, folder)
     except:
          raise Exception, "Error creating folders"
     
     #function for calculating totals
     def summary(theraster):
          summaryfile = scratchws+os.sep+"summary.dbf"
          gp.cellSize = theraster
          gp.extent = lcraster
          gp.mask = lcraster
          gp.CreateConstantRaster_sa(gp.workspace+os.sep+"intermediate"+os.sep+"unitynew", 1, "INTEGER")
          gp.ZonalStatisticsAsTable_sa (gp.workspace+os.sep+"intermediate"+os.sep+"unitynew", "Value", theraster, summaryfile, "DATA")
          cur = gp.SearchCursor(summaryfile)
          row = cur.Next()
          valuesum = row.SUM
          del cur
          del row
          if gp.Exists(summaryfile):
               gp.Delete_management(summaryfile)
          if gp.Exists(gp.workspace+os.sep+"intermediate"+os.sep+"unitynew"):
               gp.Delete_management(gp.workspace+os.sep+"intermediate"+os.sep+"unitynew")
          return valuesum    
     #function for switching between folders. 'down' changes to the intermediate folder. 
     down="down"
     def switchws(direction="up"):
          if direction=="down":
               gp.workspace=gp.GetParameterAsText(0)+"\\intermediate"
          else:
               gp.workspace=gp.GetParameterAsText(0)
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
      
     switchws(down)
     
     #function for formatting string with commas
     def thous(x):
          return re.sub(r'(\d{3})(?=\d)', r'\1,', str(x)[::-1])[::-1]

     #Set raster environment

     try:
          gp.Extent = lcraster
          gp.Mask = lcraster
          gp.CellSize=lcraster
          rasterDesc=gp.describe(lcraster)
          if resolution:
               gp.CellSize=resolution
          rasterDesc=gp.describe(lcraster)

          #check projections
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
          checkProjections(lcraster)
          if lcrasterf:
               checkProjections(lcrasterf)
          unity=[]     
          if harv_rate_current or (harv_rate_future and lcrasterf):
               #create a zero raster for later computation
               gp.CreateConstantRaster_sa("nullraster", 0, "INTEGER")
               
               #create a unity raster for later computation
               gp.CreateConstantRaster_sa("unityraster", 1, "INTEGER")
               
               #create an ascii dataset for reading header info
               gp.RasterToASCII_conversion("nullraster", "lcrasterforheader.asc")
               
               #read the number of columns and rows from the header
               myfile = open(gp.workspace+"\\lcrasterforheader.asc","r")
               cols=myfile.readline().split()[1]
               rows=myfile.readline().split()[1]
     
               #create a unity array for use later
               
               for x in range(int(rows)): 
                    innerlist=[]
                    for y in range(int(cols)): 
                         innerlist.append(1)	
                    unity.append(innerlist)

          #multiply by 0.0001 to change from sq. m. to ha
          AreaofCell = round(((int(gp.CellSize) ** 2) * 0.0001),5)

          #check projections
          #spatreflc = rasterDesc.SpatialReference
          #if harv_rate_current:
               #spatrefharvc = gp.describe(harv_rate_current).SpatialReference
          #if harv_rate_future:
               #spatrefharvf = gp.describe(harv_rate_future).SpatialReference
          #if spatreflc.Type <> 'Projected':
               #gp.addmessage("Your data does not appear to be projected.  It is assumed to be in meters")
               ##raise Exception, msgProjection
          #elif spatreflc.LinearUnitName <> 'Meter':
               #gp.addmessage("This model assumes that data is projected in meters.  You may get erroneous results")
          
          #create resampled landcover data
          if resolution:
               gp.SingleOutputMapAlgebra_sa(lcraster, "lc_res_cur"+suffix)
               if lcrasterf:
                    gp.SingleOutputMapAlgebra_sa(lcrasterf, "lc_res_fut"+suffix)

     except:
          raise Exception, msgRasterDescription

     #check fields
     def checkfields(fields, target):
          thefields = gp.listfields(target,"*","All")
          field = thefields.next()
          foundfields = []
          while field:
               foundfields.append(field.Name.upper())
               if field.Name[:5]=="Decay":
                    if not "Integer" in field.Type:
                         msgFldType = "Field "+ field.Name+" should be of integer type"
                         raise Exception, msgFldType
               field = thefields.next()

          for x in fields:
               x=x.upper()
               if not x in foundfields:
                    raise Exception, "Expected field:"+x + " in " + target+". Exiting." 
               
     def confirmfields(field, target):
          thefields = gp.listfields(target,field,"All")
          field = thefields.next()
          if field:
               return True
          else:
               return False

               
     if harv_rate_current:
          checkfields(["start_date","Decay_cur","Cut_cur"], harv_rate_current)

     if harv_rate_future:
          checkfields(["Cut_fut"], harv_rate_future)
          
     #check Carbon pools in the carbon table
     carbon_pools=[]
     thefields = gp.listfields(carbon_data,"C_*","All")
     field = thefields.next()
     if not field:
          raise Exception, "At least one Carbon pool is expected in the table: "+carbon_data
     while field:
          carbon_pools.append(field.Name)
          field = thefields.next()
    
     if harv_rate_current or (harv_rate_future and lcrasterf):     
          #read ascii header into a variable header
          rasterheader = []
          count = 0
     
          try:
               myfile = open(gp.workspace+"\\lcrasterforheader.asc","r")
               rasterheader = []
               for x in range(6):
                    rasterheader.append(myfile.readline())
               myfile.close()
          except:
               myfile.close()
               raise Exception, msgHeader


     def reclassfn(table,field,outraster,reclassraster=lcraster,scale=0):
          try:
               tempfile = open(scratchws+os.sep+"reclassfile.asc","w")

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
               sortreclassfile(scratchws+os.sep+"reclassfile.asc")
               gp.ReclassByASCIIFile_sa(reclassraster, scratchws+os.sep+"reclassfile.asc", outraster, "DATA")
          except:
               raise Exception, msgReclass +" "+ field 

     # Calculate the first 4 carbon pools (above ground, below ground, soil, dead organic matter)     
     if verbose: gp.addmessage("Calculating current storage...")
     exprtot = "Sum(" #build up expression for summing the 4 pools
     del4=[]
     for pool in carbon_pools:
          
          exprtot += pool+"_cur, "
          reclassfn(carbon_data,pool,pool+"_cur1",lcraster,100)
          
          del4.append(pool+"_cur1")
     exprtot = exprtot[:-2]
     exprtot += ")"
     
     
     #scale back the pools and convert to ha
     for pool in carbon_pools:
          expr = pool+"_cur1 / 100 * " + str(AreaofCell)
          gp.SingleOutputMapAlgebra_sa(expr, pool+"_cur")
     
     try:
          # sum the 4 pools
          gp.SingleOutputMapAlgebra_sa(exprtot, "c4pools")
     except:
          raise Exception, msgCompute4PoolsC

     #Starting HWP computation    
     #try:
          ## Convert the Decay values to a raster
          #if harv_rate_current:
               #reclassfn(carbon_data,"Decay","hlc")
               ##gp.ReclassByTable_sa(lcraster, carbon_data, "LULC", "LULC", "Decay", "hlc", "DATA")
     #except:
          #raise Exception, msgComputehalflifeC

     # ###########################################################################################
     # These functions are used for computing HWP
     # ###########################################################################################   
     def writearraytoascii(asciifile, thelist):
          # this fn writes out array list to ascii
          #Write out array list to ascii
          try:
               outputfile = open(asciifile,"w")
               for item in rasterheader:
                    outputfile.write(item)
               for x in thelist:
                    for y in x:
                         #if y==-9999:
                         outputfile.write(str(y)+' ')
                         #else:
                         #    outputfile.write(y+' ')
                    outputfile.write('\n')
               outputfile.close()
          except:
               msgHWP = "error writing ascii to file"
               outputfile.close()
               raise Exception, msgHWP
     if harv_rate_future:
          startcount = int((FDate - CDate)/2 + 1)    
     if debug: checkfile=open(r"C:\Temp\checkcarboncalc.txt","w")
     def hwpcalc(stage,TimeArray, hwpfile, theraster, DArray, CArray, stopval=0, Frequency=unity,mid=0):
          if debug: checkfile.write("\n"+hwpfile+"\n")
          HWP=[] #initialize HWP list
          try:
               for i in range(len(TimeArray)): #for each row
                    innerlist=[] #initialize a container for result rows
                    for j in range(len(TimeArray[0])): # for each value within the row (sublist)
                         sumhwp = 0 #initialize summation term
                         if CArray[i][j] > 0 and TimeArray[i][j] > 0 and DArray[i][j] > 0:
                              
                              if stage=="current":
                                   for x in range(TimeArray[i][j],0, -Frequency[i][j]): # count from number of years down to 0, or other given value
                                        sumhwp += round((1 - math.exp(-DArray[i][j]))/(DArray[i][j] * math.exp(x*DArray[i][j])),6) #add sum to the summation term
                                   innerlist.append(CArray[i][j] * sumhwp) #append result to the row
                              elif stage=="stage1":
                                   if (FDate-CDate)%2:
                                        startodd = math.ceil((FDate-CDate)/2.0)
                                   else:
                                        startodd = math.floor((FDate-CDate)/2.0)
                                   for x in range(TimeArray[i][j]+(FDate-CDate),startodd, -Frequency[i][j]): # count from number of years down to 0, or other given value
                                        sumhwp += round((1 - math.exp(-DArray[i][j]))/(DArray[i][j] * math.exp(x*DArray[i][j])),6) #add sum to the summation term
                                        if debug: 
                                             if i==19 and j==10: 
                                                  checkfile.write(str(x)+"\t")
                                                  checkfile.write(str(round((1 - math.exp(-DArray[i][j]))/(DArray[i][j] * math.exp(x*DArray[i][j])),6))+"\t")
                                                  checkfile.write(str(sumhwp)+"\n")
                                   innerlist.append(CArray[i][j] * sumhwp) #append result to the row
                              elif stage=="stage2":
                                   if (FDate-CDate)%2:
                                        startodd = math.ceil((FDate-CDate)/2.0)
                                   else:
                                        startodd = math.floor((FDate-CDate)/2.0)
                                   for x in range(startodd,0, -Frequency[i][j]): # count from number of years down to 0, or other given value                                   
                                        sumhwp += round((1 - math.exp(-DArray[i][j]))/(DArray[i][j] * math.exp(x*DArray[i][j])),6) #add sum to the summation term
                                        if debug: 
                                             if i==19 and j==10: 
                                                  checkfile.write(str(x)+"\t")
                                                  checkfile.write(str(round((1 - math.exp(-DArray[i][j]))/(DArray[i][j] * math.exp(x*DArray[i][j])),6))+"\t")
                                                  checkfile.write(str(sumhwp)+"\n")

                                   innerlist.append(CArray[i][j] * sumhwp) #append result to the row
                                   
                         else:
                              innerlist.append(-9999)
                    HWP.append(innerlist) #append completed row to the output
          except:
               msgHWP = "Error computing HWP values"
               raise Exception, msgHWP
          
          #Write HWP to ascii
          try:
               writearraytoascii(hwpfile, HWP)
               del HWP
          except:
               msgHWP = "Error writing hwp to file " + hwpfile 
               raise Exception, msgHWP
          try:
               gp.ASCIIToRaster_conversion(hwpfile, theraster, "FLOAT")
          except:
               msgHWP = "Error converting " +str(hwpfile)+ " to raster"
               raise Exception, msgHWP
     def funcc(theval):
          if float(theval) == -9999.0:
               return theval
          else:
               return float(theval)
     def funcd(theval):
          if theval==-9999 or theval==0:
               return theval
          else:
               return round((math.log(2)/theval),6)


     # ############################################################################
     #This part creates the time raster (TC)which shows the time between when the  
     #wood products harvesting started and date of current landscape                 
     # ############################################################################
     #compute the time since harvest started using the year given and the date of the current raster

     if harv_rate_current:
          
          #create half life (decay) raster
          try:
               InField = "Decay_cur"
               OutRaster = "hlc_1"
               gp.FeatureToRaster_conversion(harv_rate_current, InField, OutRaster)
               #fill half life raster with zeros

               expr="Con(IsNull(hlc_1),0,hlc_1)"
               gp.SingleOutputMapAlgebra_sa(expr, "hlc")
               
          except:
               raise Exception, msgHLC

          try:
               InField = "start_date"
               OutRaster = "start_date_1"
               gp.FeatureToRaster_conversion(harv_rate_current, InField, OutRaster)
               #fill start_date with zeros

               expr="Con(IsNull(start_date_1),0,start_date_1)"
               gp.SingleOutputMapAlgebra_sa(expr, "start_date")

               expr = str(CDate) +" - "+ str("start_date")
               gp.SingleOutputMapAlgebra_sa(expr, "TC_1")
               expr = "Con(TC_1 == "+str(CDate)+", 0, TC_1)"
               gp.SingleOutputMapAlgebra_sa(expr, "TC_2")

               #fill TC with zeros
               expr="Con(IsNull(TC_2),0,TC_2)"
               gp.SingleOutputMapAlgebra_sa(expr, "TC")
          except:
               raise Exception, msgTC
          if confirmfields("C_den_cur", harv_rate_current) and confirmfields("BCEF_cur", harv_rate_current): 
               #create raster of WD
               try:
                    InField = "C_den_cur"
                    OutRaster = "C_den_cur"
                    gp.FeatureToRaster_conversion(harv_rate_current, InField, OutRaster)
     
               except:
                    raise Exception, "Error processing C_den_cur"
               
               #create raster of BCEF
               try:
                    InField = "BCEF_cur"
                    OutRaster = "BCEF_cur"
                    gp.FeatureToRaster_conversion(harv_rate_current, InField, OutRaster)
     
               except:
                    raise Exception, "Error processing BCEF_cur"
          if lcrasterf and confirmfields("C_den_fut", harv_rate_future) and confirmfields("BCEF_fut", harv_rate_future): 
               #create raster of WD
               try:
                    InField = "C_den_fut"
                    OutRaster = "C_den_fut"
                    gp.FeatureToRaster_conversion(harv_rate_future, InField, OutRaster)
     
               except:
                    raise Exception, "Error processing C_den_fut"
               
               #create raster of BCEF
               try:
                    InField = "BCEF_fut"
                    OutRaster = "BCEF_fut"
                    gp.FeatureToRaster_conversion(harv_rate_future, InField, OutRaster)
     
               except:
                    raise Exception, "Error processing BCEF_fut"
       
          try:
               # Set local variables
               InField = "Cut_cur"
               OutRaster = "cut_1"

               # Create raster of the harvest rates
               gp.FeatureToRaster_conversion(harv_rate_current, InField, OutRaster)

               #fill harv rate
               expr="Con(IsNull(cut_1),0,cut_1)"
               gp.SingleOutputMapAlgebra_sa(expr, "cut_cur")

               #Calculate the amount of C harvested per cell            
               expr="cut_cur * " + str(AreaofCell)
               gp.SingleOutputMapAlgebra_sa(expr, "CharvestC")    

               #write the time since harvest started raster to ascii
               gp.RasterToASCII_conversion("TC", "tcascii.asc")
          except:
               raise Exception, msgHarvest
          try:
               # Calculate the frequency
               InField = "Freq_cur"
               OutRaster = "freq_1"

               # Create raster of the frequencies
               gp.FeatureToRaster_conversion(harv_rate_current, InField, OutRaster)

               #fill frequencies
               expr="Con(IsNull(freq_1),0,freq_1)"
               gp.SingleOutputMapAlgebra_sa(expr, "freqc")

               #write the frequencies raster to ascii
               gp.RasterToASCII_conversion("freqc", "freqcascii.asc")
          except:
               raise Exception, msgFrequency

          # ############################################################################
          #This part calculates the harvested wood products carbon pool 
          #using time since harvest, harvest rate and the half life
          # ############################################################################
          
          if harv_rate_current:
               try:
                    #write time since harvest started into an array  
                    inputfile = open(gp.workspace+"\\tcascii.asc","r")
                    TC=[]
                    for line in inputfile.readlines():
                         TC.append([])
                         for i in line.split():
                              TC[-1].append(i)
                    TC=TC[6:]
                    TC=[[(int(y)) for y in x] for x in TC]
                    inputfile.close()

                    inputfile = open(gp.workspace+"\\freqcascii.asc","r")
                    FREQC=[]
                    for line in inputfile.readlines():
                         FREQC.append([])
                         for i in line.split():
                              FREQC[-1].append(i)
                    FREQC=FREQC[6:]
                    FREQC=[[(int(y)) for y in x] for x in FREQC]
                    inputfile.close()               


                    gp.RasterToASCII_conversion("HLC", "hlcascii.asc")
                    #write half life into an array            
                    inputfile = open(gp.workspace+"\\hlcascii.asc","r")
                    HLC=[]
                    for line in inputfile.readlines():
                         HLC.append([])
                         for i in line.split():
                              HLC[-1].append(i)
                    HLC=HLC[6:]
                    HLC=[[(int(y)) for y in x] for x in HLC]
                    inputfile.close()

                    #write d, a term used in calculating HWP into an array. Created from half life               
                    DC=[]
                    DC=[[funcd(y) for y in x] for x in HLC]

                    #write harvested carbon into an array   
                    gp.RasterToASCII_conversion("Charvestc", "charvestascii.asc")
                    inputfile = open(gp.workspace+"\\charvestascii.asc","r")
                    C=[]
                    for line in inputfile.readlines():
                         C.append([])
                         for i in line.split():
                              C[-1].append(i)
                    C=C[6:]
                    C=[[funcc(y) for y in x] for x in C]
                    inputfile.close()
                    #Calculate HWP using function created above               
                    hwpcalc("current",TC, gp.workspace+"\\hwpc.asc", "C_HWP_Cur_1", DC, C,0,FREQC)
                    
               except:
                    raise Exception, msgHWP+ gp.GetMessages(2)

          # ############################################################################
          #This part calculates the total C pool for current landcover 
          #by adding the components from above (4)
          # ############################################################################
     try:
          #fill up HWP raster with zeros
          
          if harv_rate_current:
               expr="Con(IsNull(C_HWP_Cur_1),0,C_HWP_Cur_1)"
               gp.SingleOutputMapAlgebra_sa(expr, "C_HWP_Cur")

          #add first 4 pools to HWP raster
          switchws()
          if harv_rate_current:
               expr="Sum("+gp.workspace+"\\intermediate\\c4pools,"+gp.workspace+"\\intermediate\\C_HWP_Cur)"
          else:
               expr=gp.workspace+"\\intermediate\\c4pools"

          gp.SingleOutputMapAlgebra_sa(expr, gp.workspace+"\\Output\\tot_C_cur"+suffix)
          total_current_carbon = summary(gp.workspace+"\\Output\\tot_C_cur"+suffix)
          

          #Calculate economic value of standing stock and hwp
          if econ=='true':
               if verbose: gp.addmessage("Calculating economic valuation...")
               expr=gp.workspace+"\\Output\\tot_C_cur"+suffix+" * "+ str(p)
               gp.SingleOutputMapAlgebra_sa(expr, gp.workspace+"\\Output\\value_stor"+suffix)
               total_value_stored_carbon = summary(gp.workspace+"\\Output\\value_stor"+suffix)
          
     except:
          raise Exception, msgTotalCurrentC
     # ############################################################################
     #This part calculates Biomass removed from each polygon. 
     #It is only calculated if harvest rate polygon is given
     # ############################################################################
     switchws(down)
     if harv_rate_current:
          if confirmfields("C_den_cur", harv_rate_current) and confirmfields("BCEF_cur", harv_rate_current): 
               if verbose: gp.addmessage("Calculating biomass and volume...")
               expr = "Ceil(TC / Float(freqc))"
               gp.SingleOutputMapAlgebra_sa(expr, "rd_cur")
               expr = "rd_cur * charvestc / C_den_cur"
               gp.SingleOutputMapAlgebra_sa(expr, "bio_hwp_cur"+suffix)
               
               expr = "bio_hwp_cur"+suffix+" / BCEF_cur"
               gp.SingleOutputMapAlgebra_sa(expr, "vol_hwp_cur"+suffix)
                              
     
     # ############################################################################
     #This part calculates the first 4 C pools for the scenario landcover. 
     #This is a repetition of the first step, can me modified into a fn
     # ############################################################################
     
     if lcrasterf: 
          #run this section only if a scenario landcover has been specified.
          #Calculate the first 4 carbon pools (above ground, below ground, soil, dead organic matter)
          
          if verbose: gp.addmessage("Calculating future storage...")
          exprtot = "Sum(" #build up expression for summing the 4 pools
          for pool in carbon_pools:
               
               exprtot += pool+"_fut, "
               reclassfn(carbon_data,pool,pool+"_fut1",lcrasterf,100)
          exprtot = exprtot[:-2]
          exprtot += ")"
          
          
          #scale back the pools and convert to ha
          for pool in carbon_pools:
               expr = pool+"_fut1 / 100 * " + str(AreaofCell)
               gp.SingleOutputMapAlgebra_sa(expr, pool+"_fut")
          
          try:
               # sum the 4 pools
               gp.SingleOutputMapAlgebra_sa(exprtot, "c4poolsf")
          except:
               raise Exception, msgCompute4PoolsF

          # ############################################################################
          #This part, dependent on scenario selection calculates the C stored in harvested 
          #wood products. 
          # ############################################################################
          if harv_rate_future:
               #############################################################################
               ## Prepare Frequency array (FREQF)
               # Calculate the frequency for future
               try:
                    #check if frequency field is present in scenario harvest table.  If yes, create frequency raster otherwise make the 
                    #frequency equal to the current frequency.
                    
                    #check for the three variations of the case
                    thefields1 = gp.listfields(harv_rate_future,"freq_fut","All")
                    thefields2 = gp.listfields(harv_rate_future,"FREQ_FUT","All")
                    thefields3 = gp.listfields(harv_rate_future,"Freq_fut","All")
                    if thefields1.next():
                         nofrequency = 0
                    elif thefields2.next():
                         nofrequency = 0
                    elif thefields3.next():
                         nofrequency = 0
                    else:
                         nofrequency = 1
                    
                    #if frequency present, do calculation
                    if not nofrequency:
                         InField = "Freq_fut"
                         OutRaster = "freq_1"
          
                         # Create raster of the frequencies
                         gp.FeatureToRaster_conversion(harv_rate_future, InField, OutRaster)
          
                         #fill frequencies
                         expr="Con(IsNull(freq_1),0,freq_1)"
                         gp.SingleOutputMapAlgebra_sa(expr, "freqf")
          
                         #write the frequencies raster to ascii
                         gp.RasterToASCII_conversion("freqf", "freqfascii.asc")
                    #otherwise use frequency from before
                    else:
                         #write the frequencies raster to ascii if frequencies exist else use unity raster
                         if gp.exists("freqc"):
                              gp.RasterToASCII_conversion("freqc", "freqfascii.asc")
                         else:
                              gp.RasterToASCII_conversion("unityraster", "freqfascii.asc")
                              
               except:
                    raise Exception, msgFrequency
               #write the frequency to an array
               inputfile = open(gp.workspace+"\\freqfascii.asc","r")
               FREQF=[]
               for line in inputfile.readlines():
                    FREQF.append([])
                    for i in line.split():
                         FREQF[-1].append(i)
               FREQF=FREQF[6:]
               FREQF=[[(int(y)) for y in x] for x in FREQF]
               inputfile.close()  
               
               
               #############################################################################
               ## Prepare Harvested C amounts array CF)
               try:
                    # Create a raster of amount of C harvested from harvest polygons
                    gp.FeatureToRaster_conversion(harv_rate_future, "Cut_fut", "cutf_1")

                    #fill harv rate
                    expr="Con(IsNull(cutf_1),0,cutf_1)"
                    gp.SingleOutputMapAlgebra_sa(expr, "cut_fut")

                    expr="cut_fut * " + str(AreaofCell)
                    gp.SingleOutputMapAlgebra_sa(expr, "CharvestF")

               except:
                    raise Exception, msgComputeTF
               
               gp.RasterToASCII_conversion("CharvestF", "Charvestfascii.asc")
               inputfile = open(gp.workspace+"\\Charvestfascii.asc","r")
               CF=[]
               for line in inputfile.readlines():
                    CF.append([])
                    for i in line.split():
                         CF[-1].append(i)
               CF=CF[6:]
               CF=[[funcc(y) for y in x] for x in CF]
               inputfile.close()
               
               #############################################################################
               ## Prepare Decay/Halflife array HLF)
               if harv_rate_future:
                    try:
                         InField = "Decay_fut"
                         OutRaster = "hlf_1"
                         gp.FeatureToRaster_conversion(harv_rate_future, InField, OutRaster)
                         #fill half life raster with zeros
                         expr="Con(IsNull(hlf_1),0,hlf_1)"
                         gp.SingleOutputMapAlgebra_sa(expr, "hlf")
                    except:
                         raise Exception, msgComputehalflifeF
               
               gp.RasterToASCII_conversion("HLF", "hlfascii.asc")
               inputfile = open(gp.workspace+"\\hlfascii.asc","r")
               HLF=[]
               for line in inputfile.readlines():
                    HLF.append([])
                    for i in line.split():
                         HLF[-1].append(i)
               HLF=HLF[6:]
               HLF=[[(int(y)) for y in x] for x in HLF]
               inputfile.close()
               
               #############################################################################
               ## Prepare Time array (TF) TF is based on TC except mid value is added to stage 1
               mid = (FDate - CDate)/2
               #Create the DF array from HLF
               DF=[[funcd(y) for y in x] for x in HLF]   
         
               
               #use function to calculate hwp
               
               #if there was harvesting in current landscape then calculate stage 1
               if harv_rate_current:
                    #HWP stage 1  
                    hwpcalc("stage1", TC, gp.workspace+"\\hwpfst1.asc", "C_HWP_futst1", DC, C, 0, FREQC,mid)
                    #fill up HWP rasters with zeros
                    expr="Con(IsNull(C_HWP_futst1),0,C_HWP_futst1)"
                    gp.SingleOutputMapAlgebra_sa(expr, "C_HWP_futst11")
               
               else:
                    TC = unity
               
               #calculate biomass and volume for stage 1
               if harv_rate_current:
                    if confirmfields("C_den_cur", harv_rate_current) and confirmfields("BCEF_cur", harv_rate_current): 
                         if verbose: gp.addmessage("Calculating future biomass and volume...")
                         if (FDate-CDate)%2:
                              biomid = math.ceil((FDate-CDate)/2.0)
                         else:
                              biomid = math.floor((FDate-CDate)/2.0)
                         expr = "Ceil((TC + "+str(FDate - CDate - biomid)+") / freqc)"
                         gp.SingleOutputMapAlgebra_sa(expr, "rd_futst1")
                         expr = "rd_futst1 * charvestc / C_den_cur"
                         gp.SingleOutputMapAlgebra_sa(expr, "bio_hwp_fut1")
                         
                         expr = "bio_hwp_fut1"+" / BCEF_fut"
                         gp.SingleOutputMapAlgebra_sa(expr, "vol_hwp_fut1")
                    
               
                    
               
               #HWP stage 2
               hwpcalc("stage2", TC, gp.workspace+"\\hwpfst2.asc", "C_HWP_futst2", DF, CF, 0, FREQF,mid)
               
               #calculate biomas and volume for stage 2
               if confirmfields("C_den_fut", harv_rate_future) and confirmfields("BCEF_fut", harv_rate_future): 
                    if (FDate-CDate)%2:
                         biomid = math.ceil((FDate-CDate)/2.0)
                    else:
                         biomid = math.floor((FDate-CDate)/2.0)
                    expr = "Ceil("+str(biomid)+" / freqf)"
                    gp.SingleOutputMapAlgebra_sa(expr, "rd_futst2")
                    expr = "(rd_futst2 * charvestf) / C_den_fut"
                    gp.SingleOutputMapAlgebra_sa(expr, "bio_hwp_fut21")
                    
                    expr = "bio_hwp_fut21 / BCEF_fut"
                    gp.SingleOutputMapAlgebra_sa(expr, "vol_hwp_fut21")

                    #fill up future biomas rasters with zeros
                    expr="Con(IsNull(bio_hwp_fut21),0,bio_hwp_fut21)"
                    gp.SingleOutputMapAlgebra_sa(expr, "bio_hwp_fut2")  
                    
                    #fill up future volume rasters with zeros
                    expr="Con(IsNull(vol_hwp_fut21),0,vol_hwp_fut21)"
                    gp.SingleOutputMapAlgebra_sa(expr, "vol_hwp_fut2")  
               
               #fill up future HWP rasters with zeros
               expr="Con(IsNull(C_HWP_futst2),0,C_HWP_futst2)"
               gp.SingleOutputMapAlgebra_sa(expr, "C_HWP_futst21")        

               #sum up all HWP. If doesnt apply, will have a value of zero
               if harv_rate_current:
                    expr = "C_HWP_futst11 + C_HWP_futst21"
               else:
                    expr = "C_HWP_futst21"
               gp.SingleOutputMapAlgebra_sa(expr, "C_HWP_fut")          

               #sum up all biomass and volume. If doesnt apply, will have a value of zero
               if harv_rate_current:
                    expr = "bio_hwp_fut1 + bio_hwp_fut2"
               else:
                    expr = "bio_hwp_fut2"
               gp.SingleOutputMapAlgebra_sa(expr, "bio_hwp_fut"+suffix)   
               
               if harv_rate_current:
                    expr = "vol_hwp_fut1 + vol_hwp_fut2"
               else:
                    expr = "vol_hwp_fut2"
               gp.SingleOutputMapAlgebra_sa(expr, "vol_hwp_fut"+suffix)   
               


          #Add the pools (first 4 and the HWP) for scenario
          switchws()
          if harv_rate_future:
               expr="Sum("+gp.workspace+"\\intermediate\\c4poolsf,"+gp.workspace+"\\intermediate\\C_HWP_fut)"
          else:
               expr=gp.workspace+"\\intermediate\\c4poolsf"            
          gp.SingleOutputMapAlgebra_sa(expr, gp.workspace+"\\Output\\tot_C_fut"+suffix)
          total_future_carbon = summary(gp.workspace+"\\Output\\tot_C_fut"+suffix)

          #Compute sequestration
          if verbose: gp.addmessage("Calculating sequestration...")  
          try:
               expr=gp.workspace+"\\Output\\tot_C_fut"+suffix+" - "+gp.workspace+"\\Output\\tot_C_cur"+suffix
               gp.SingleOutputMapAlgebra_sa(expr, gp.workspace+"\\Output\\sequest"+suffix)
               total_sequest_carbon = summary(gp.workspace+"\\Output\\sequest"+suffix)

          except:
               raise Exception, msgSequestration

          try:
               if econ=='true':
                    if verbose: gp.addmessage("Calculating economic valuation...")
                    #Compute economic valuation
                    looptime=FDate-CDate
                    discount=0
                    for x in range(looptime):
                         discount += p / (((1 + r)**x)*((1 + c)**x))
                    expr="("+gp.workspace+"\\Output\\sequest"+suffix+" / " + str(looptime) + ") * " + str(discount)
                    gp.SingleOutputMapAlgebra_sa(expr, gp.workspace+"\\Output\\val_seq"+suffix)
                    total_value_sequestered = summary(gp.workspace+"\\Output\\val_seq"+suffix)
          except:
               raise Exception, msgEcon
          switchws(down)
 
     #write parameters
     if suffix:
          suffix = "_"+suffix
     switchws()
     
     
     if total_current_carbon:
          parameters.append("Total current carbon: "+thous(total_current_carbon)+" Mg")
     if total_future_carbon:
          parameters.append("Total scenario carbon: "+thous(total_future_carbon)+" Mg")  
     if total_sequest_carbon:
          parameters.append("Total sequestered carbon: "+thous(total_sequest_carbon)+" Mg")  
     if total_value_stored_carbon:
          parameters.append("Total value of current stored carbon: "+thous(total_value_stored_carbon)+" Mg")   
     if total_value_sequestered:
          parameters.append("Total value of sequestered carbon: "+thous(total_value_sequestered)+" Mg")
     parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
     parameters.append("Suffix: "+str(suffix))
     
     parafile = open(gp.workspace+"\\Output\\Carbon_"+now.strftime("%Y-%m-%d-%H-%M")+suffix+".txt","w") 
     parafile.writelines("CARBON MODEL PARAMETERS\n")
     parafile.writelines("________________________\n\n")
     
     for para in parameters:
          parafile.writelines(para+"\n")
          parafile.writelines("\n")
     
  
     parafile.close()
     if debug: checkfile.close()
     #Cleanup
     switchws(down)
     if gp.exists(scratchws):
          gp.delete_management(scratchws)     
     if cleanup:
          if verbose: gp.addmessage("Cleaning up...")  

          del1=["unityraster","unitynew","charvestascii.asc","charvestc","charvestf","charvestfascii.asc","cut_cur","cut_fut","hlcascii.asc","hlfascii.asc","hwpc.asc","hwpf.asc","hwpfst1.asc","hwpfst2.asc","reclassfile.asc","C_den_cur","C_den_fut","hlc_1","hlf_1","freqfascii.asc","freqcascii.asc","freqf","freqc","freq_1"]
          del2=["start_date","tc","tc_1","tc_2","tcascii.asc","tconst","tf","tf_stage1","tf_stage2","tfascii.asc","tfst1ascii.asc","tfst2ascii.asc","hlc","hlf","resampled","resampledf","rd_futst1","rd_futst2","rd_cur","vol_hwp_fut1","vol_hwp_fut2","vol_hwp_fut21"]
          del3=["nullraster","lcrasterforheader.asc","c_hwp_cur_1","c_hwp_fut_1","c_hwp_fut_2","c_hwp_futst1","c_hwp_futst11","c_hwp_futst2","c_hwp_futst21","cut_1","cutf_1","start_date_1","c4pools","c4poolsf","bio_hwp_fut1","bio_hwp_fut21","bio_hwp_fut2"]
          deletelist=del1 + del2 + del3 + del4
          for data in deletelist:
               if gp.exists(data):
                    gp.delete_management(data)

except Exception, ErrorDesc:
     gp.addmessage(gp.GetMessages())
     gp.AddError(str(ErrorDesc))
