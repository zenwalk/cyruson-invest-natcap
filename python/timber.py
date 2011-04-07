# ---------------------------------------------------------------------------
# Author: Nasser Olwero (nasser.olwero@wwfus.org)
# Description: Timber
# Date: June 3, 2008
# ---------------------------------------------------------------------------

# Import system modules
import sys, string, os, arcgisscripting, math, datetime

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out any necessary licenses
gp.CheckOutExtension("spatial")

#Set output handling
gp.OverwriteOutput = 1

msgArguments = "Problem with arguments: " + gp.getmessages()
msgGetPrice = "Error reading the price: " + gp.getmessages()
msgPlantation = "Error reading plantation data: "+ gp.getmessages()

verbose = True
cleanup = True

try:
    try:
        #Set workspace
        parameters = []
        now = datetime.datetime.now()
        parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
        gp.Workspace = gp.GetParameterAsText(0)
        parameters.append("Workspace: "+gp.Workspace) 
        plantationvector = gp.GetParameterAsText(1)
        parameters.append("Managed timber forest parcels: "+plantationvector) 
        plantation_data = gp.GetParameterAsText(2)
        parameters.append("Production table: "+plantation_data) 
        r = int(gp.GetParameterAsText(3))
        parameters.append("Market discount rate: "+str(r)) 
        suffix = gp.GetParameterAsText(4)
        if suffix is "#":
            suffix = ""
        parameters.append("Suffix: "+str(suffix))

    except:
        raise Exception, msgArguments
    down="down"
    def switchws(direction="up"):
        if direction=="down":
            gp.workspace=gp.GetParameterAsText(0)+"\\intermediate"
        else:
            gp.workspace=gp.GetParameterAsText(0)   
    #check and create folders
    try:
        thefolders=["Output","intermediate"]
        for folder in thefolders:
            if not gp.exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace, folder)
    except:
        raise Exception, "Error creating folders"
    switchws()
    #check fields
    def checkfields(fields, target):
        thefields = gp.listfields(target,"*","All")
        field = thefields.next()
        foundfields = []
        while field:
            foundfields.append(field.Name.upper())
            if field.Name=="Decay":
                if not "Integer" in field.Type:
                    msgFldType = "Field "+ field.Name+" should be of integer type"
                    raise Exception, msgFldType
            field = thefields.next()
            
        for x in fields:
            x=x.upper()
            if not x in foundfields:
                raise Exception, "Expected field:"+x + " in " + target+". Exiting." 
    checkfields(["Parcel_ID","T","Parcl_area","Price","Harv_cost","Maint_cost","Perc_harv","Harv_mass","Freq_harv","Immed_harv","BCEF"], plantation_data)
    def getPrice(thespp):
        try:
            rows = gp.SearchCursor(prices)
            row = rows.next()
            while row:
                if row.spp == thespp:
                    return row.price
                row = rows.next()
        except:
            raise Exception, msgGetPrice
        
    def confirmfields(field, target):
        thefields = gp.listfields(target,field,"All")
        field = thefields.next()
        if field:
            return True
        else:
            return False
        
    try:       
        cur = gp.SearchCursor(plantation_data)
        row = cur.next()
    except:
        raise Exception, msgPlantation
    #make a copy of plantation vector
    try:
        out_shape = gp.workspace+"\\Output\\timber"+suffix+".shp"
        if not gp.exists(out_shape):
            gp.copy_management(plantationvector, out_shape)
        if gp.exists(out_shape):
            if not confirmfields("TNPV", out_shape): gp.addfield(out_shape, "TNPV", "DOUBLE")
            if not confirmfields("TBiomass", out_shape): gp.addfield(out_shape, "TBiomass", "DOUBLE")
            if not confirmfields("TVolume", out_shape): gp.addfield(out_shape, "TVolume", "DOUBLE")
        else:
            gp.addmessage(out_shape +" does not exist")
            
    except:
        raise Exception, "Error copying plantation shapefile.  Please ensure the shapefile is not in use and you have write permissions to the workspace"
    try:
        switchws()
    except:
        raise Exception, "Error creating file"
    shapecur = gp.UpdateCursor(out_shape)
    shaperow = shapecur.next()
    try:
        while row:
            n = row.Parcel_ID
            
            Parcl_area = row.Parcl_area
            Perc_harv=row.Perc_harv
            Price = row.Price
            Harv_cost = row.Harv_cost
            Maint_cost = row.Maint_cost
            T = row.T
            Freq_harv = row.Freq_harv
            Immed_harv = row.Immed_harv
            Harv_mass = row.Harv_mass
            BCEF = row.BCEF
            
            if Freq_harv > T:
                gp.addmessage("Frequency for parcel "+ str(n) +" is higher than the total period(T). Frequency will be set equal to T.")
                Freq_harv = T
            
            VH = (Perc_harv/100.0) * ((Price * Harv_mass) - Harv_cost)
            
            if Immed_harv.upper()=="Y" or Immed_harv.upper()=="YES":
                NPV1 = 0
                NPV2 = 0
                for t in range(0,int(math.ceil(T/float(Freq_harv)))):
                    NPV1 += VH / ((1 + float(r/100.0))**(Freq_harv*t))
                    
                for t in range(0,T):
                    NPV2 += Maint_cost / ((1 + float(r/100.0))**t)
                biomass = Parcl_area * (Perc_harv/100.0) * Harv_mass * math.ceil(T/float(Freq_harv))
                    
            elif Immed_harv.upper()=="N" or Immed_harv.upper()=="NO":
                NPV1 = 0
                NPV2 = 0
                for t in range(1,int(math.floor(T/float(Freq_harv))+1)):
                    NPV1 += VH / ((1 + float(r/100.0))**((Freq_harv*t)-1))
                for t in range(0,T):
                    NPV2 += Maint_cost / ((1 + float(r/100.0))**t)
                biomass = Parcl_area * (Perc_harv/100.0) * Harv_mass * math.floor(T/float(Freq_harv))
            volume = biomass / float(BCEF)
            NPV = NPV1 - NPV2    
            TNPV = Parcl_area * NPV

            shaperow.SetValue("TNPV", round(TNPV,3))
            shaperow.SetValue("TBiomass", round(biomass,3))
            shaperow.SetValue("TVolume", round(volume,3))

            shapecur.UpdateRow(shaperow)
            shaperow = shapecur.next()
            row = cur.next()
    except:
        raise Exception, "Error processing plantation table"  + gp.getmessages()
    del cur
    del shapecur
    switchws()
    #write parameters
    if suffix:
        suffix = "_"+suffix
    switchws()
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    
    parafile = open(gp.workspace+"\\Output\\Timber_"+now.strftime("%Y-%m-%d-%H-%M")+suffix+".txt","w") 
    parafile.writelines("TIMBER MODEL PARAMETERS\n")
    parafile.writelines("________________________\n\n")
    
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()
    gp.workspace = gp.workspace+"\\Output"
    
    #Cleanup
    if cleanup:
        if verbose: gp.addmessage("Cleaning up...")  
        for data in []:
            if gp.exists(data):
                gp.delete_management(data)

except Exception, ErrorDesc:
    gp.AddError(str(ErrorDesc))
