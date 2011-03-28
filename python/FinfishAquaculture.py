# Marine InVEST: Finfish Aquaculture Model
# Authors: Gregg Verutes, Jodie Toft, Apollo Xi
# Coded for ArcGIS 9.3 and 10
# 02/16/11

# import modules
import sys, string, os, datetime
import arcgisscripting, shlex, numpy
from win32com.client import Dispatch
from math import *

# create the geoprocessor object
gp = arcgisscripting.create()
# set output handling
gp.OverwriteOutput = 1
# check out any necessary licenses
gp.CheckOutExtension("management")

# error messages
msgArguments = "\nProblem with arguments."
msgFarmOp = "\nError reading in farm operation data."
msgGrowthSim = "\nError during growth simulation."
msgPopTables = "\nError while tallying results and creating outputs."

try:
    try:
        # get parameters
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
        gp.workspace = gp.GetParameterAsText(0)
        parameters.append("Workspace: "+ gp.workspace)
        gp.scratchWorkspace = gp.GetParameterAsText(0)
        parameters.append("Scratch Workspace: "+ gp.scratchWorkspace)
        FinFishFC = gp.GetParameterAsText(1)
        parameters.append("Finfish Farm Location: "+ FinFishFC)
        FarmIDField = gp.GetParameterAsText(2)
        parameters.append("Farm Identifier Name: "+ FarmIDField)
        aa = float(gp.GetParameterAsText(3))
        parameters.append("Fish growth paramter (a): "+ str(aa))
        bb = float(gp.GetParameterAsText(4))
        parameters.append("Fish growth paramter (b): "+ str(bb))
        TempData = gp.GetParameterAsText(5)
        parameters.append("Daily Water Temperature at Farm Table: "+ TempData)
        FarmOperations = gp.GetParameterAsText(6)
        parameters.append("Farm Operations Table: "+ FarmOperations)
        ValuationBoolean = gp.GetParameterAsText(7)
        parameters.append("Run Valuation? "+ ValuationBoolean)
        MarketPrice = float(gp.GetParameterAsText(8))
        parameters.append("Market Price Per Kilogram of Processed Fish: "+ str(MarketPrice))
        PctPriceCost = float(gp.GetParameterAsText(9))
        parameters.append("Fraction of Market Price that Accounts for Costs: "+ str(PctPriceCost))
        DiscRateDay = float(gp.GetParameterAsText(10))
        parameters.append("Daily Market Discount Rate: "+ str(DiscRateDay))
    except:
        raise Exception, msgArguments + gp.GetMessages(2)
    
    try:
        thefolders=["Output"]
        for folder in thefolders:
            if not gp.exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        raise Exception, "Error creating folders"
    
    
    # local variables 
    outputws = gp.workspace + os.sep + "Output" + os.sep
    
    FinFishHarv = outputws + "Finfish_Harvest.shp"
    ResultsHTML = outputws + "HarvestResults_"+now.strftime("%Y-%m-%d-%H-%M")+".html"
    hrvwght_kg = outputws + "hrvwght_kg"
    npv_usd_1k = outputws + "npv_usd_1k"

    # growth variables
    F_p = 0.45
    F_f = 0.3
    F_c = 0.07
    C_p = 5650
    C_f = 9450
    C_c = 4100
    P_p = 0.18
    P_f = 0.18
    A_p = 0.89
    A_f = 0.92
    A_c = 0.5
    alpha = 11
    gamma = 0.8
    N_p = 0.166666667
    NtoP = 6
    O_p = 1.89
    O_f = 2.91
    O_c = 1.07
    Ox_Am_Nit = 0.47
    Temp = 11.5
    tau  = 0.08
    
    # make calculations based on input from growth parameters
    delta = F_p*C_p+F_f*C_f+F_c*C_c
    E_p = (F_p*C_p)/delta
    E_f = (F_f*C_f)/delta
    E_c = (F_c*C_c)/delta
    C_fi = (P_p*C_p)+(P_f*C_f)
    FL = (1-A_p)*E_p+(1-A_f)*E_f+(1-A_c)*E_c
    BC = 0.3*A_p*E_p+0.05*(A_f*E_f+A_c*E_c)
    eps = 1-FL-BC
    eps_star = eps-0.15*E_p*A_p
    C_fi_star = 0.85*C_p*P_p+C_f*P_f
    
    
    ##############################################
    ################ CHECK INPUTS ################
    ##############################################

    def ckProjection(data):
        dataDesc = gp.describe(data)
        spatreflc = dataDesc.SpatialReference
        if spatreflc.Type <> 'Projected':
            gp.AddError(data +" does not appear to be projected.  It is assumed to be in meters.")
            raise Exception
        if spatreflc.LinearUnitName <> 'Meter':
            gp.AddError("This model assumes that "+data+" is projected in meters for area calculations.  You may get erroneous results.")
            raise Exception

    ckProjection(FinFishFC)
    
    # if conducting valuation, check that all valuation parameters exist
    if ValuationBoolean == "true":
        inputs = [MarketPrice, PctPriceCost, DiscRateDay]
        for x in inputs:
            if x == "":
                gp.AddError("\nOne or more of the valuation input parameters was not defined.")
                raise Exception
    
    
    #####################################################################
    ################## OPEN FARM OPERATIONS DATA ########################
    #####################################################################
    
    try:
        x1App = Dispatch("Excel.Application")
        x1App.Visible = 0
        x1App.DisplayAlerts = 0
        xlBook = x1App.Workbooks.Open(FarmOperations[:-(1+len(FarmOperations.split("\\")[-1]))])
        FarmOperations_path = FarmOperations.split("\\")
        FarmOperations_sheet = FarmOperations_path[-1]
        xlSheet = xlBook.Worksheets(FarmOperations_sheet[:-1])
        
        # grab variables
        DressWght = xlSheet.Range("b2").Value
        MortRate = xlSheet.Range("b3").Value
        SimLength = xlSheet.Range("b4").Value
        
        # search for number of rows
        row = 11
        col = 1
        
        bottom = row
        while xlSheet.Cells(bottom+1, col).Value not in [None, '']:
            bottom += 1
        
        # find number of farms in operations SS
        FarmNum = xlSheet.Range(xlSheet.Cells(10,1), xlSheet.Cells(bottom,1)).Value
        gp.AddMessage("\nNumber of farms: "+str(len(FarmNum)))
        
        FarmSpecsList = numpy.zeros(len(FarmNum)*6, dtype=numpy.float64)
        FarmSpecsArray = numpy.reshape(FarmSpecsList, (len(FarmNum),6))
        FarmHarvList = []
        FarmHarvArray = []
        FarmTotalsList = numpy.zeros(len(FarmNum)*4, dtype=numpy.float64)
        FarmTotalsArray = numpy.reshape(FarmTotalsList, (len(FarmNum),4))
        
        rowNum = 10
        for row in range(0,len(FarmNum)):
            FarmSpecsArray[row][0] = row+1 # Farm ID
            FarmSpecsArray[row][1] = xlSheet.Range("b"+str(rowNum)).Value # weight of fish at start (kg)
            FarmSpecsArray[row][2] = xlSheet.Range("c"+str(rowNum)).Value # weight of fish at harvest (kg)
            FarmSpecsArray[row][3] = xlSheet.Range("d"+str(rowNum)).Value # number of fish in farm
            FarmSpecsArray[row][4] = xlSheet.Range("e"+str(rowNum)).Value # start day (ID) for growing
            FarmSpecsArray[row][5] = xlSheet.Range("f"+str(rowNum)).Value # length of fallowing period (days)
            rowNum += 1
        
        x1App.ActiveWorkbook.Close(SaveChanges=0)
        x1App.Quit()
    except:
        raise Exception, msgFarmOp
    
    
    #####################################################################
    ####################### OPEN TEMP DATA ##############################
    #####################################################################
    
    x2App = Dispatch("Excel.Application")
    x2App.Visible = 0
    x2App.DisplayAlerts=0
    xlBook = x2App.Workbooks.Open(TempData[:-(1+len(TempData.split("\\")[-1]))])
    TempData_path = TempData.split("\\")
    TempData_sheet = TempData_path[-1]
    x2Sheet = xlBook.Worksheets(TempData_sheet[:-1])
    
    # letters list
    letterList = ['c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','aa','ab','ac','ad','ae','af','ag','ah','ai','aj']
    
    # grab start temp and put in list
    startDayList = [] # day of the year the farm starts (outplanting)
    tempCountList = [] # temperature counter added to tempCountList (for going through the days)
    farmHarvList = [] # [no, yes, hold]
    stopFallowList = [] # fallowing counter --> starts at fallow length and counts down
    HarvCountList = [] # number of full harvests for each farm
    numHarvDaysList = [] # number of days since start of the harvest
    outplantDayList = [] # day of outplanting for each cycle
    yearCountList = [] # calendar year counter
    yearHarvestList = [] # year that outplanting occurred for each farm
    
    for farm in range(0,len(FarmNum)):
        bottom = 6
        while x2Sheet.Cells(bottom, 1).Value <> FarmSpecsArray[farm][4]:  ## removed the farm-1
            bottom += 1
        startDayList.append(bottom)
        tempCountList.append(1)
        farmHarvList.append('no')
        stopFallowList.append(0)
        HarvCountList.append(0)
        numHarvDaysList.append(2)  # starts at "2" because first day happens only once
        outplantDayList.append(bottom-5)
        yearCountList.append(1)
        yearHarvestList.append(1)

    # growth simulation variables    
    DayCount = 1
    dayList = []
    dayArray = []
    
    
    #####################################################################
    ######################### SIMULATION ################################
    #####################################################################
    try:
        gp.AddMessage("\nConducting growth simulation for "+str(int(SimLength))+" years...")
        ## on first day, grabs input temp data (this approach happens only once)
        ## model goes at daily time step, but is at different days of the year for different farms.
        for farm in range(0,len(FarmNum)):
            dayList.append(DayCount) # first day
            dayList.append(x2Sheet.Range(letterList[farm]+str(startDayList[farm])).Value) # start temp 
            dayList.append(exp((x2Sheet.Range(letterList[farm]+str(startDayList[farm])).Value) * tau)) # temp_effect
            dayList.append((FarmSpecsArray[farm][1])*1000) # start weight in grams
        dayArray.append(dayList)
        dayList = []
        DayCount += 1
        
        ## remaining code --> model loops until simulation end day is reach
        while DayCount <= (SimLength*365):
            for farm in range(0,len(FarmNum)):
                # next day number
                dayList.append(DayCount)
                # reset counter for each farm if it reaches the end of the year
                if (startDayList[farm]+tempCountList[farm]) > 370: 
                    startDayList[farm] = 5
                    tempCountList[farm] = 1
                    yearCountList[farm] = yearCountList[farm] + 1
                # next day's temp 
                dayList.append(x2Sheet.Range(letterList[farm]+str(startDayList[farm]+tempCountList[farm])).Value)
                # next day's temp_effect
                dayList.append(exp((x2Sheet.Range(letterList[farm]+str(startDayList[farm]+tempCountList[farm])).Value) * tau))
                ## fallowing period logic: since mass of harvested farm = 0 there will be no growth nor incorrect harvesting
                ## if fallowing counter == 0, outplant and assign start weight
                if stopFallowList[farm] == 0 and farmHarvList[farm] == 'hold':
                    dayList.append((FarmSpecsArray[farm][1])*1000) 
                    farmHarvList[farm] = 'no' # change farmHarvList to 'no'
                    numHarvDaysList[farm] = 1
                    outplantDayList[farm] = (startDayList[farm] + tempCountList[farm])-5
                    yearHarvestList[farm] = yearCountList[farm]
                else:
                    # next day's weight (temp_effect and day before's weight) 
                    dayList.append(aa*(dayArray[(DayCount-2)][((farm*4)+3)])**bb*(exp((x2Sheet.Range(letterList[farm]+str(startDayList[farm]+tempCountList[farm])).Value)*tau))+(dayArray[(DayCount-2)][((farm*4)+3)]))
            dayArray.append(dayList)         

            ## always check if each farm has reached harvest weight --> if so, harvest
            for farm in range(0,len(FarmNum)):
                w = dayList[((farm*4)+3)]
                if w >= (FarmSpecsArray[farm][2]*1000.0): # harvest weight
                    farmHarvList[farm] = 'yes'
                
            ## harvest the farm --> add mass to FarmSpecsArray[row][6] --> also, allow for fallowing time (dormant) before growing for that farm again
            for farm in range(0,len(farmHarvList)):
                if farmHarvList[farm] == 'yes':
                    ## add total weight of processed fish (TPW) to each farm's harvest mass tally (this includes mortality ratio)
                    ## TPW = (weight at harvest * pct of fish wght after processing * (number of fish * EXP(-daily mortality rate * number of days since start of harvest)))
                    TPW = ((dayArray[(DayCount-1)][((farm*4)+3)]/1000)*DressWght*(FarmSpecsArray[farm][3]*numpy.exp(-(MortRate*numHarvDaysList[farm]))))
                    if ValuationBoolean == "true":
                        NetRev = (TPW*(MarketPrice*(1-PctPriceCost)))
                        NPV = (NetRev*(1/((1+DiscRateDay)**DayCount)))
                    else:
                        NetRev = 0
                        NPV = 0

                    # harvest information stored in FarmHarvList & FarmHarvArray
                    FarmHarvList = [farm+1, HarvCountList[farm]+1, int(numHarvDaysList[farm]+FarmSpecsArray[farm][5]), TPW, NetRev, NPV, outplantDayList[farm], yearHarvestList[farm]] # eventually populated into HTML table #2
                    FarmHarvArray.append(FarmHarvList)
                    FarmTotalsArray[farm][3] = FarmTotalsArray[farm][3] + TPW # add this harvests weight to farm-specific tally
                                                  
                    # reset values
                    dayArray[(DayCount-1)][((farm*4)+3)] = 0.0 # set that day's weight in array to 0.0
                    # next two variables: part of fallowing counter
                    farmHarvList[farm] = 'hold'
                    stopFallowList[farm] = (FarmSpecsArray[farm][5]+1)
                    HarvCountList[farm] = HarvCountList[farm] + 1
                    numHarvDaysList[farm] = 0 # restart number of days since last harvest list
        
            dayList = [] # clear day list
            DayCount = DayCount + 1 # move to next day
            if DayCount == int((SimLength*365)*0.25):
                gp.AddMessage("...25% completed")
            elif DayCount == int((SimLength*365)*0.50):
                gp.AddMessage("...50% completed")
            elif DayCount == int((SimLength*365)*0.75):
                gp.AddMessage("...75% completed")
            
            for k in range(0,len(tempCountList)):     
                tempCountList[k] = tempCountList[k] + 1
                stopFallowList[k] = stopFallowList[k] - 1
                numHarvDaysList[k] = numHarvDaysList[k] + 1
    except:
        raise Exception, msgGrowthSim 
            
              
    x2App.ActiveWorkbook.Close(SaveChanges=0)
    x2App.Quit()
    
    
    ############################################################################################################################################################

    try:    
        # populate "FarmTotalsArray" after growth simulation
        for row in range(len(FarmTotalsArray)):
            FarmTotalsArray[row][0] = row+1 # Farm ID
            for rows in range(len(FarmHarvArray)):
                if FarmHarvArray[rows][0] == FarmTotalsArray[row][0]:
                    FarmTotalsArray[row][1] = FarmTotalsArray[row][1] + FarmHarvArray[rows][5]  # sum up all NPV values for each cycle
            FarmTotalsArray[row][2] = HarvCountList[row] # Number of cycles

        # create html file output
        htmlfile = open(ResultsHTML, "w")
        htmlfile.write("<html>\n")
        htmlfile.write("<title>Marine InVEST</title>")
        htmlfile.write("<CENTER><H1>Aquaculture Model (Finfish Harvest)</H1></CENTER><br>")
        htmlfile.write("This page contains results from running the Marine InVEST Finfish Aquaculture model.<p>Cells highlighted in yellow are values that were also populated in the attribute table of the netpens feature class.  Cells highlighted in red should be interpreted as null values since valuation was not selected.<br>")
        htmlfile.write("<br><HR><H2>Farm Operations (input)</H2><table border=\"1\", cellpadding=\"5\"><tr>")
        htmlfile.write("<td align=\"center\"><b>Farm ID Number</b></td><td align=\"center\"><b>Weight of fish at start</b><br>(kg)</td><td align=\"center\"><b>Weight of fish at harvest</b><br>(kg)</td><td align=\"center\"><b>Number of fish in farm</b></td><td align=\"center\"><b>Start day for growing</b><br>(1-365)</td><td align=\"center\"><b>Length of fallowing period</b><br>(days)</td></tr>")
        for i in range(0,len(FarmSpecsArray)):
            htmlfile.write("<tr><td align=\"center\">"+str(int(FarmSpecsArray[i][0]))+"</td><td align=\"center\">"+str(FarmSpecsArray[i][1])+"</td><td align=\"center\">"+str(FarmSpecsArray[i][2])+"</td><td align=\"center\">"+str(int(FarmSpecsArray[i][3]))+"</td><td align=\"center\">"+str(int(FarmSpecsArray[i][4]))+"</td><td align=\"center\">"+str(int(FarmSpecsArray[i][5]))+"</td></tr>")
        htmlfile.write("</table>")
        htmlfile.write("<br><HR><H2>Farm Harvesting (output)</H2><table border=\"1\", cellpadding=\"5\"><tr>")
        htmlfile.write("<td align=\"center\"><b>Farm ID Number</b></td><td align=\"center\"><b>Cycle Number</b></td><td align=\"center\"><b>Days Since Outplanting Date</b><br>(including fallowing period)</td><td align=\"center\"><b>Harvested Weight</b><br>(kg/cycle)</td><td align=\"center\"><b>Net Revenue</b><br>(Thousands of $)</td><td align=\"center\"><b>Net Present Value</b><br>(Thousands of $)</td><td align=\"center\"><b>Outplant Day</b><br>(Julian Day)</td><td align=\"center\"><b>Outplant Year</b></td></tr>")
        if ValuationBoolean == "true":
            for j in range(0,len(FarmHarvArray)):
                htmlfile.write("<tr><td align=\"center\">"+str(FarmHarvArray[j][0])+"</td><td align=\"center\">"+str(FarmHarvArray[j][1])+"</td><td align=\"center\">"+str(FarmHarvArray[j][2])+"</td><td align=\"center\">"+str(int(FarmHarvArray[j][3]))+"</td><td align=\"center\">"+str(int(FarmHarvArray[j][4]/1000))+"</td><td align=\"center\">"+str(int(FarmHarvArray[j][5]/1000))+"</td><td align=\"center\">"+str(FarmHarvArray[j][6])+"</td><td align=\"center\">"+str(FarmHarvArray[j][7])+"</td></tr>")
        else:
            for j in range(0,len(FarmHarvArray)):
                htmlfile.write("<tr><td align=\"center\">"+str(FarmHarvArray[j][0])+"</td><td align=\"center\">"+str(FarmHarvArray[j][1])+"</td><td align=\"center\">"+str(FarmHarvArray[j][2])+"</td><td align=\"center\">"+str(int(FarmHarvArray[j][3]))+"</td><td align=\"center\" bgcolor=\"#ff0000\">"+str(int(FarmHarvArray[j][4]/1000))+"</td><td align=\"center\" bgcolor=\"#ff0000\">"+str(int(FarmHarvArray[j][5]/1000))+"</td><td align=\"center\">"+str(FarmHarvArray[j][6])+"</td><td align=\"center\">"+str(FarmHarvArray[j][7])+"</td></tr>")
        htmlfile.write("</table>")
        htmlfile.write("<br><HR><H2>Farm Result Totals (output)</H2><table border=\"1\", cellpadding=\"5\"><tr>")
        htmlfile.write("<td align=\"center\"><b>Farm ID Number</b></td><td align=\"center\"><b>Net Present Value</b><br>(Thousands of $)<br>(for duration of model run)</td><td align=\"center\"><b>Number of completed<br>harvest cycles</b></td><td align=\"center\"><b>Total Volume Harvested</b><br>(kg) (after processing occurs)</td></tr>")
        if ValuationBoolean == "true":
            for k in range(0,len(FarmTotalsArray)):
                htmlfile.write("<tr><td align=\"center\">"+str(int(FarmTotalsArray[k][0]))+"</td><td align=\"center\" bgcolor=\"#ffff00\">"+str(int(FarmTotalsArray[k][1]/1000))+"</td><td align=\"center\" bgcolor=\"#ffff00\">"+str(int(FarmTotalsArray[k][2]))+"</td><td align=\"center\" bgcolor=\"#ffff00\">"+str(int(FarmTotalsArray[k][3]))+"</td></tr>")
        else:
            for k in range(0,len(FarmTotalsArray)):
                htmlfile.write("<tr><td align=\"center\">"+str(int(FarmTotalsArray[k][0]))+"</td><td align=\"center\" bgcolor=\"#ff0000\">"+str(int(FarmTotalsArray[k][1]/1000))+"</td><td align=\"center\" bgcolor=\"#ffff00\">"+str(int(FarmTotalsArray[k][2]))+"</td><td align=\"center\" bgcolor=\"#ffff00\">"+str(int(FarmTotalsArray[k][3]))+"</td></tr>")    
        htmlfile.write("</table>")
        htmlfile.write("</html>")
        htmlfile.close()

        # copy and populate input FC attribute table
        gp.CopyFeatures_management(FinFishFC, FinFishHarv, "", "0", "0", "0")
        
        # add three new fields
        gp.AddField_management(FinFishHarv, "TOT_CYCLES", "LONG", "5", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        gp.AddField_management(FinFishHarv, "HRVWGHT_KG", "DOUBLE", "19", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        gp.AddField_management(FinFishHarv, "NPV_USD_1K", "DOUBLE", "19", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "") 
        
        for i in range(0,len(FarmNum)):
            # search and activate conditions
            ID = int(FarmTotalsArray[i][0])
            SrchCondition = str(FarmIDField)+" = "+str(ID)
            ActivCondition = str(FarmIDField)+"; TOT_CYCLES; HRVWGHT_KG; NPV_USD_1K"
             
            # translate information from array into FinFishFC
            cur = gp.UpdateCursor(FinFishHarv, SrchCondition, "", ActivCondition)
            row = cur.Next()
            row.SetValue("TOT_CYCLES", int(FarmTotalsArray[i][2])) # Total Number of Cycles
            row.SetValue("HRVWGHT_KG", int(FarmTotalsArray[i][3])) # Total Harvest Weight (tons)
            if ValuationBoolean == "true":
                row.SetValue("NPV_USD_1K", int(FarmTotalsArray[i][1]/1000)) # NPV
            if ValuationBoolean == "false":
                row.SetValue("NPV_USD_1K", 0) # NPV
            cur.UpdateRow(row)
        del cur 
        del row

        # create raster outputs
        gp.FeatureToRaster_conversion(FinFishHarv, "HRVWGHT_KG", hrvwght_kg, "25") # 25 meters
        if ValuationBoolean == "true":
            gp.FeatureToRaster_conversion(FinFishHarv, "NPV_USD_1K", npv_usd_1k, "25") # 25 meters
        
    except:
        raise Exception, msgPopTables

    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile = open(gp.workspace+"\\Output\\parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("FINFISH AQUACULTURE MODEL PARAMETERS\n")
    parafile.writelines("____________________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()
    gp.AddMessage("")    
    del gp
except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())
    gp.AddError(str(ErrorDesc))