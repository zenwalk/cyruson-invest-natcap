# Marine InVEST: Coastal Protection (Profile Builder)
# Authors: Greg Guannel, Gregg Verutes
# 09/07/10

## TO DO ##
## REQUIRE LANDPT TO BE WGS84
## AVOID DOUBLE PRINTING ON PROFILE_TXT

import numpy as num
import CPf_SignalSmooth as SignalSmooth
import string, sys, os, time, datetime, shlex
import fpformat, operator
import arcgisscripting

from win32com.client import Dispatch
from scipy.interpolate import interp1d
from scipy import optimize
from math import *
from matplotlib import *
from pylab import *

# create the geoprocessor object
gp = arcgisscripting.create()

# set output handling
gp.OverwriteOutput = 1
# check out any necessary extensions
gp.CheckOutExtension("management")
gp.CheckOutExtension("analysis")
gp.CheckOutExtension("conversion")
gp.CheckOutExtension("spatial")

# error messages
msgArguments = "Problem with arguments."

try:
    # get parameters
    parameters = []
    now = datetime.datetime.now()
    parameters.append("Date and Time: "+ now.strftime("%Y-%m-%d %H:%M"))
    gp.workspace = gp.GetParameterAsText(0)
    parameters.append("Workspace: "+ gp.workspace)
    LandPoint = gp.GetParameterAsText(1)
    parameters.append("Land Point: "+ LandPoint)
    LandPoly = gp.GetParameterAsText(2)
    parameters.append("Land Polygon: "+ LandPoly)
    InputTable = gp.GetParameterAsText(3)
    parameters.append("Profile Builder Table: "+ InputTable)
    ProfileQuestion = gp.GetParameterAsText(4)
    parameters.append("Do you have a nearshore bathymetry GIS layer?: "+ ProfileQuestion)
    BathyGrid = gp.GetParameterAsText(5)
    parameters.append("IF 1: Bathymetric Grid (DEM): "+ BathyGrid)
    CSProfile = gp.GetParameterAsText(6)
    parameters.append("IF 2: Upload Your Cross-Shore Profile: "+ CSProfile)
    WW3_Pts = gp.GetParameterAsText(7)
    parameters.append("Wave Watch 3 Model Data: "+ WW3_Pts)
except:
    raise Exception, msgArguments + gp.GetMessages(2)

try:
    thefolders=["intermediate","Output"]
    for folder in thefolders:
        if not gp.exists(gp.workspace+folder):
            gp.CreateFolder_management(gp.workspace, folder)
except:
    raise Exception, "Error creating folders"

# variables (hard-coded)
BufferDist = 150
SampInterval = 1
TransectDist = 1.0
BearingsNum = 16
RadLineDist = 100000

# intermediate and output directories
outputws = gp.workspace + os.sep + "Output" + os.sep
interws = gp.workspace + os.sep + "intermediate" + os.sep

PT1 = interws + "PT1.shp"
PT2 = interws + "PT2.shp"
PT1_Z = interws + "PT1_Z.shp"
PT2_Z = interws + "PT2_Z.shp"
LandPoint_Buff = interws + "LandPoint_Buff.shp"
LandPoint_Geo = interws + "LandPoint_Geo.shp"
Shoreline = interws + "Shoreline.shp"
Shoreline_Buff_Clip = interws + "Shoreline_Buff_Clip.shp"
Shoreline_Buff_Clip_Diss = interws + "Shoreline_Buff_Clip_Diss.shp"
PtsCopy = interws + "PtsCopy.shp"
PtsCopy2 = interws + "PtsCopy2.shp"
PtsCopyLR = interws + "PtsCopy2_lineRotate.shp"
Fetch_AOI = interws + "Fetch_AOI.shp"
UnionFC = interws + "UnionFC.shp"
SeaPoly = interws + "SeaPoly.shp"
seapoly_rst = interws + "seapoly_rst"
seapoly_e = interws + "seapoly_e"
PtsCopyEL = interws + "PtsCopy2_eraseLand.shp"
PtsCopyExp = interws + "PtsCopy2_explode.shp"
PtsCopyExp_Lyr = interws + "PtsCopy2_explode.lyr"
WW3_Pts_prj = interws + "WW3_Pts_prj.shp"
costa_ww3 = interws + "costa_ww3"
LandPoint_WW3 = interws + "LandPoint_WW3.shp"

Profile_Txt = outputws + "Profile_Txt.txt"
Profile_Pts = outputws + "Profile_Pts.shp"
Profile_Plot = outputws + "Profile_Plot.png"
Profile_HTML = outputws + "Profile_HTML.html"
Fetch_Vectors = outputws + "Fetch_Vectors.shp"

# various functions and checks
def AddField(FileName, FieldName, Type, Precision, Scale):
    fields = gp.ListFields(FileName, FieldName)
    field_found = fields.Next()
    if field_found:
        gp.DeleteField_management(FileName, FieldName)
    gp.AddField_management(FileName, FieldName, Type, Precision, Scale, "", "", "NON_NULLABLE", "NON_REQUIRED", "")
    return FileName

def getDatum(thedata):
    desc = gp.describe(thedata)
    SR = desc.SpatialReference
    if SR.Type == "Geographic":
        strDatum = SR.DatumName         
    else:
        gp.OutputCoordinateSystem = SR
        strSR = str(gp.OutputCoordinateSystem)
        gp.OutputCoordinateSystem = ""
        n1 = strSR.find("GEOGCS")
        n2 = strSR.find("PROJECTION")
        strDatum = strSR[n1:n2-1]
    return strDatum

def ckProjection(data):
    dataDesc = gp.describe(data)
    spatreflc = dataDesc.SpatialReference
    if spatreflc.Type <> 'Projected':
        gp.AddError(data +" does not appear to be projected.  It is assumed to be in meters.")
        raise Exception
    if spatreflc.LinearUnitName <> 'Meter':
        gp.AddError("This model assumes that "+data+" is projected in meters for area calculations.  You may get erroneous results.")
        raise Exception
    
def grabProjection(data):
    dataDesc = gp.describe(data)
    sr = dataDesc.SpatialReference
    gp.OutputCoordinateSystem = sr
    strSR = str(gp.OutputCoordinateSystem)
    return strSR

def compareProjections(LandPoint, LandPoly):
    if gp.describe(LandPoint).SpatialReference.name <> gp.describe(LandPoly).SpatialReference.name:
        gp.AddError("Projection Error: "+LandPoint+" is in a different projection from the LandPoly data.  The two inputs must be the same projection to calculate depth profile.")
        raise Exception

def Indexed(x,value):
    mylist=abs(x-value);    
    if isinstance(x,num.ndarray):
        mylist=mylist.tolist();
    minval=min(mylist)
    ind=[i for i, v in enumerate(mylist) if v == minval];ind=ind[0];
    return ind

# function to create point transects
def PTCreate(PTType, midx, midy, TransectDist):
    if PTType == 1:
        y1 = midy + TransectDist
        y2 = midy - TransectDist
        x1 = midx
        x2 = midx
    elif PTType == 2:
        y1 = midy
        y2 = midy 
        x1 = midx + TransectDist
        x2 = midx - TransectDist
    elif PTType == 3:
        y1 = NegRecip*(TransectDist) + midy
        y2 = NegRecip*(-TransectDist) + midy
        x1 = midx + TransectDist
        x2 = midx - TransectDist
    elif PTType == 4:
        y1 = midy + TransectDist
        y2 = midy - TransectDist
        x1 = (TransectDist/NegRecip) + midx
        x2 = (-TransectDist/NegRecip) + midx
    elif PTType == 5:
        y1 = midy + TransectDist
        y2 = midy - TransectDist
        x1 = (TransectDist/NegRecip) + midx
        x2 = (-TransectDist/NegRecip) + midx
    elif PTType == 6:
        y1 = NegRecip*(TransectDist) + midy
        y2 = NegRecip*(-TransectDist) + midy
        x1 = midx + TransectDist
        x2 = midx - TransectDist
    return x1, y1, x2, y2


# check that correct inputs were provided based on 'ProfileQuestion'
if ProfileQuestion == "(1) Yes":
    if not BathyGrid:
        gp.AddError("A bathymetry grid input is required.")
        raise Exception
elif ProfileQuestion == "(2) No, but I will upload a cross-shore profile":
    if not CSProfile:
        gp.AddError("A cross-shore profile input is required.")
        raise Exception

# check that three inputs are projected
ckProjection(LandPoint)
ckProjection(LandPoly)
ckProjection(BathyGrid)
# check that two inputs have the same projection
compareProjections(LandPoint, LandPoly)
geo_projection = getDatum(LandPoint)

# import Profile Builder info from Excel file
xlApp = Dispatch("Excel.Application")
xlApp.Visible=0
xlApp.DisplayAlerts=0
xlApp.Workbooks.Open(InputTable)
cell = xlApp.Worksheets("Profile Generator Input")

# check what type of data user has
WaveClimateCheck = cell.Range("e92").Value # 1 model chooses, 2 if enter
DuneCheck = cell.Range("e93").Value # 1 don't know, 2 no, 3 don't know, 4 Yes
Diam = cell.Range("e14").Value # Sediment diam [mm]

# load wave climate info
if WaveClimateCheck == 2:
    He=cell.Range("h25").Value # Effective Wave Height
    Hm=cell.Range("i25").Value # Modal Wave Height
    Tm=cell.Range("j25").Value # Modal Wave Period
    
##else: get data from WW3 [GV]##

# load tide information        
MSL = cell.Range("d30").Value # Mean Sea Level
HT = cell.Range("e30").Value # High tide elevation
HAT = cell.Range("f30").Value # High tide elevation

# foreshore    
Slope = cell.Range("f38").Value # foreshore slope = 1/Slope
m =1.0/Slope; #Bed slope

# beach parameters
A = cell.Range("e94").Value # sediment scale factor
hc=-ceil(1.57*He) # closure depth 

# berm and dune
BermCrest = cell.Range("f48").Value
BermLength = cell.Range("g48").Value

# if user doesn't know dune size, estimate from Short and Hesp
if DuneCheck==1 or DuneCheck==3: 
    Hb=0.39*9.81**(1.0/5)*(Tm*Hm**2)**(2.0/5)
    a=0.00000126
    b=num.sqrt(3.61**2+1.18*(1.56*9.81*(Diam/1000.0)**3/a**2)**(1.0/1.53))-3.61
    ws=(a*b**1.53)/(Diam/1000) #Fall velo
    RTR=HT/Hb
        
    if RTR>3: # in this case, beach is not Wave Dominated, can't know value, so take zero
        DuneCrest=0
        BermLength=50
    else: # else, beach is wave dominated, we read Short and Hesp
        Type=Hb/(ws*Tm)
        if Type<3:
            DuneCrest=5
        elif Type<4:
            DuneCrest=10
        elif Type<5:
            DuneCrest=12
        elif Type<6:
            DuneCrest=20
        else:
            DuneCrest=23
                
elif DuneCheck==2: # no Dunes
    DuneCrest=0
    BermLength=50 # beach has no dune and infinitely long berm
        
elif DuneCheck==4: # user has data
    DuneCrest = cell.Range("j58").Value

# vegetation
do1v=cell.Range("e67").Value # starting depth of vegetation field
Lveg=cell.Range("f67").Value # ending depth of vegetation field

# modifications in case profile is uploaded
SlopeMod1=cell.Range("e77").Value
SlopeMod2=cell.Range("e78").Value
SlopeMod3=cell.Range("e79").Value
OffMod1=cell.Range("f77").Value
OffMod2=cell.Range("f78").Value
OffMod3=cell.Range("f79").Value
ShoreMod1=cell.Range("g77").Value
ShoreMod2=cell.Range("g78").Value
ShoreMod3=cell.Range("g79").Value

xlApp.ActiveWorkbook.Close(SaveChanges=0)
xlApp.Quit()
# put bathy profiles together #
x=num.arange(0,10000.1,.1) # long axis

# nearshore bathy profile
if ProfileQuestion == "(1) Yes": # model extracts value from GIS layers
    gp.AddMessage("\nCreating Point Transects...")
    # create transect and read transect file
    gp.Buffer_analysis(LandPoint, LandPoint_Buff, str(BufferDist)+" Meters", "FULL", "ROUND", "NONE", "")
    gp.Extent = LandPoint_Buff
    gp.PolygonToLine_management(LandPoly, Shoreline)
    gp.Extent = ""    
    gp.Clip_analysis(Shoreline, LandPoint_Buff, Shoreline_Buff_Clip, "")
    # check to make sure that clipped shoreline is not empty FC
    if gp.GetCount_management(Shoreline_Buff_Clip) == 0:
        gp.AddError("Shoreline was not found within "+str(BufferDist)+" meters of 'LandPoint' input.  \
                     Either increase the buffer distance or move the 'LandPoint' input closer to the coastline.")
        raise Exception
    gp.Dissolve_management(Shoreline_Buff_Clip, Shoreline_Buff_Clip_Diss, "", "", "MULTI_PART", "UNSPLIT_LINES")

    # set coordinate system to same projection (in meters) as the shoreline point input
    gp.outputCoordinateSystem = LandPoint
    cur = gp.UpdateCursor(LandPoint)
    row = cur.Next()
    feat = row.Shape
    midpoint = feat.Centroid
    midList = shlex.split(midpoint)
    midList = [float(s) for s in midList]
    midx = midList[0]
    midy = midList[1]
    del cur
    del row

    # grab coordinates of the start and end of the coastline segment
    cur = gp.SearchCursor(Shoreline_Buff_Clip_Diss)
    row = cur.Next()
    counter = 1
    feat = row.Shape
    firstpoint = feat.FirstPoint
    lastpoint = feat.LastPoint
    startList = shlex.split(firstpoint)
    endList = shlex.split(lastpoint)
    startx = float(startList[0])
    starty = float(startList[1])
    endx = float(endList[0])
    endy = float(endList[1])

    # diagnose the type of perpendicular transect to create (PerpTransType)
    PerpTransType = 0
    if starty==endy or startx==endx:
        if starty == endy:
            y1 = midy + TransectDist
            y2 = midy - TransectDist
            x1 = midx
            x2 = midx
            PerpTransType = 1
        if startx == endx:
            y1 = midy
            y2 = midy 
            x1 = midx + TransectDist
            x2 = midx - TransectDist
            PerpTransType = 2
    else:
        # get the slope of the line
        m = ((starty - endy)/(startx - endx))
        # get the negative reciprocal
        NegRecip = -1*((startx - endx)/(starty - endy))

        if m > 0:
            # increase x-values, find y
            if m >= 1:
                y1 = NegRecip*(TransectDist) + midy
                y2 = NegRecip*(-TransectDist) + midy
                x1 = midx + TransectDist
                x2 = midx - TransectDist
                PerpTransType = 3
            # increase y-values, find x
            if m < 1:
                y1 = midy + TransectDist
                y2 = midy - TransectDist
                x1 = (TransectDist/NegRecip) + midx
                x2 = (-TransectDist/NegRecip) + midx
                PerpTransType = 4
        if m < 0:
            # add to x, find y-values
            if m >= -1:
            # add to y, find x-values
                y1 = midy + TransectDist
                y2 = midy - TransectDist
                x1 = (TransectDist/NegRecip) + midx
                x2 = (-TransectDist/NegRecip) + midx
                PerpTransType = 5
            if m < -1:
                y1 = NegRecip*(TransectDist) + midy
                y2 = NegRecip*(-TransectDist) + midy
                x1 = midx + TransectDist
                x2 = midx - TransectDist
                PerpTransType = 6
    del cur
    del row

    # grab projection spatial reference from 'LandPoint'
    dataDesc = gp.describe(LandPoint)
    spatialRef = dataDesc.SpatialReference
    gp.CreateFeatureClass_management(interws, "PT1.shp", "POINT", "#", "#", "#", spatialRef)
    gp.CreateFeatureClass_management(interws, "PT2.shp", "POINT", "#", "#", "#", spatialRef)

    # create two point transects, each point is 1 meter away from the previous    
    cur1 = gp.InsertCursor(PT1)
    cur2 = gp.InsertCursor(PT2)
    while TransectDist <= 10000:
        # call 'PTCreate' function to use the correct perpendicular transect formula based on coastline slope (m)
        x1, y1, x2, y2 = PTCreate(PerpTransType, midx, midy, TransectDist)
        row1 = cur1.NewRow()
        pnt = gp.CreateObject("POINT")
        pnt.x = x1
        pnt.y = y1
        row1.shape = pnt
        cur1.InsertRow(row1)
        row2 = cur2.NewRow()
        pnt = gp.CreateObject("POINT")
        pnt.x = x2
        pnt.y = y2
        row2.shape = pnt
        cur2.InsertRow(row2)
        TransectDist = TransectDist + 1
    del cur1, row1
    del cur2, row2

    # extract depth values from 'BathyGrid' to point transects
    gp.ExtractValuesToPoints_sa(PT1, BathyGrid, PT1_Z, "INTERPOLATE")
    gp.ExtractValuesToPoints_sa(PT2, BathyGrid, PT2_Z, "INTERPOLATE")
    PT1_Z = AddField(PT1_Z, "PT_ID", "LONG", "", "")        
    gp.CalculateField_management(PT1_Z, "PT_ID", "[FID]+1", "VB")
    PT2_Z = AddField(PT2_Z, "PT_ID", "LONG", "", "")        
    gp.CalculateField_management(PT2_Z, "PT_ID", "[FID]+1", "VB")    

    # create depth lists of two point transects
    Dmeas1 = []
    cur = gp.UpdateCursor(PT1_Z)
    row = cur.Next()
    while row:
        Dmeas1.append(row.GetValue("RASTERVALU"))
        cur.UpdateRow(row)
        row = cur.next()
    del cur, row
    Dmeas2 = []
    cur = gp.UpdateCursor(PT2_Z)
    row = cur.Next()
    while row:  
        Dmeas2.append(row.GetValue("RASTERVALU"))
        cur.UpdateRow(row)
        row = cur.next()
    del cur, row

    # find which point transect hits water first
    DepthStart1 = 1
    for DepthValue1 in Dmeas1:
        if DepthValue1 < 0.0 and DepthValue1 <> -9999.0:
            break
        DepthStart1 = DepthStart1 +1
    DepthStart2 = 1
    for DepthValue2 in Dmeas2:
        if DepthValue2 < 0.0 and DepthValue2 <> -9999.0:
            break
        DepthStart2 = DepthStart2 +1

    # create final lists of cross-shore distance (xd) and depth (Dmeas)
    xd = []   
    Dmeas = []
    counter = 0
    if DepthStart1 < DepthStart2:
        for i in range(DepthStart1-1,len(Dmeas1)):
            if Dmeas1[i] < 0.0 and Dmeas1[i] <> -9999.0:
                xd.append(counter)
                Dmeas.append(Dmeas1[i])
                counter = counter + 1
            else:
                break
    else:
        for j in range(DepthStart2-1,len(Dmeas2)):
            if Dmeas2[j] < 0.0 and Dmeas2[j] <> -9999.0:
                xd.append(counter)
                Dmeas.append(Dmeas2[j])
                counter = counter + 1
            else:
                break

    if len(Dmeas1) == 0 and len(Dmeas2) == 0:
        gp.AddError("Neither transect overlaps the seas.  Please check the location of your 'LandPoint' and bathymetry inputs.")
        raise Exception

    # create txt profile for erosion portion
    file = open(Profile_Txt, "a")
    for i in range(0,len(Dmeas)):
        file.writelines(str(xd[i])+"\t"+str(Dmeas[i])+"\n")
    file.close()

    # create final point transect file
    if DepthStart1 < DepthStart2:
        gp.Select_analysis(PT1_Z, Profile_Pts, "\"PT_ID\" > "+str(DepthStart1-1)+" AND \"PT_ID\" < "+str(DepthStart1+counter))
    else:
        gp.Select_analysis(PT2_Z, Profile_Pts, "\"PT_ID\" > "+str(DepthStart2-1)+" AND \"PT_ID\" < "+str(DepthStart2+counter))

# read in user's cross-shore profile
elif ProfileQuestion == "(2) No, but I will upload a cross-shore profile":
    TextData = open(CSProfile,"r") 
    xd = []
    Dmeas = []
    for line in TextData.readlines():
        linelist= [float(s) for s in line.split("\t")] # split the list by tab delimiter
        xd.append(linelist[0])
        Dmeas.append(linelist[1])

if ProfileQuestion <> "(3) No, please assume an equilibrium beach profile":
    # smooth profile and create x axis
    Dx=xd;L=len(Dmeas) # length of original data
    Dmeas=num.array(Dmeas);xd=num.array(xd)
    Dmeas=Dmeas[::-1] # reverse order so deeper values starts at x=0 (for plot only)

    # smooth data
    yd=SignalSmooth.smooth(Dmeas,max(int(len(Dmeas)/20),5),'flat')

    # find A value for erosion model
    Yeq=yd;Xeq=xd
    if abs(Yeq[0])>abs(Yeq[-1]): # make sure that profile starts with deeper point
        Yeq=Yeq[::-1]
        
    la=num.argmin(abs(Yeq-hc)) # locate closure depth
    Yeq=Yeq[0:la];Xeq=Xeq[0:la]

    fitfunc=lambda p, ix: p[0]*(ix)**(2.0/3) # target function for eq. beach profile
    errfunc=lambda p, ix, iy: fitfunc(p, ix) - iy # distance to the target function
    p0=[0.1] # initial guess for the parameters
    p1,success=optimize.leastsq(errfunc, p0[:], args=(Xeq,-Yeq))
    A=p1[0] # sediment scale factor

# equilibrium beach profile; in case we don't have nearshore bathy
if ProfileQuestion == "(3) No, please assume an equilibrium beach profile":
    temp=float(2)/(3)
    yd=x**(temp)# eq. profile
    out=Indexed(yd,-hc)
    yd=num.delete(yd,out,None)# water depths down to hc
    xd=num.delete(x,out,None)
    yd=-yd[::-1] # reverse order so deeper values starts at x=0
    Dmeas=yd;Dx=xd

# prepare whole profile for erosion model
if ProfileQuestion == "(1) Yes" or ProfileQuestion == "(3) No, please assume an equilibrium beach profile": # model extracts value from GIS layer
    # add foreshore, berm and dune
    # foreshore
    yf=float(1)/Slope*x+yd[-1]
    Above=num.nonzero(yf>BermCrest)
    xf=num.delete(x,Above[0],None)
    yf=num.delete(yf,Above[0],None) # remove values that are above MHW

    # berm and dune
    if DuneCheck == 2: # no dunes are present, just berm   # 1 = DK, 2 = No, 3 = Maybe 4 = Have data
        xb=num.arange(0,100,1)
        yb=num.array(len(xb)*[0.0])+BermCrest # horizontal berm 100m long
    elif DuneCheck == 1 and RTR > 3: # user doesn't know, and not wave Dominated: no dunes, just berm
        xb=num.arange(0,100,1)
        yb=num.array(len(xb)*[0.0])+BermCrest # horizontal berm 100m long
    elif DuneCheck == 3 and RTR > 3: # user doesn't know, and not wave Dominated: no dunes, just berm
        xb=num.arange(0,100,1)
        yb=num.array(len(xb)*[0.0])+BermCrest # horizontal berm 100m long          
    else: # dune exists; we'll create it as sinusoid for representation
        xb=num.arange(0,1000.1,1)
        if BermLength <> 0: # there is a berm
            # berm profile
            yb=num.array(len(xb)*[0.0])+BermCrest
            Toe=abs(xb-BermLength).argmin()# locate toe to separate berm and dune
            
            # dune profile
            DuneWidth=3*DuneCrest # width of sinusoid....won't use...for plotting purposes only
            yb[Toe:-1]=float(DuneCrest)*num.sin(2*pi*(xb[Toe:-1]-xb[Toe])/float(DuneWidth))+(BermCrest)
            DunePlotEnd=xb[Toe]+3*DuneWidth
            DunePlotSmall=xb[Toe]+DuneCrest

            out=num.arange(DunePlotEnd,len(yb),1)
            yb[DunePlotSmall:-1]=yb[DunePlotSmall:-1]/10
            yb[DunePlotSmall:-1]=yb[DunePlotSmall:-1]+yb[DunePlotSmall-1]-yb[DunePlotSmall]

            xb=num.delete(xb,out,None)
            yb=num.delete(yb,out,None); yb[-1]=yb[-2];
           
        else: # there's no berm; pretty much same code as above
           # dune profile
            DuneWidth=3*DuneCrest # width of sinusoid....won't use...for plotting purposes only

            yb[Toe:-1]=float(DuneCrest)*num.sin(2*pi*(xb[Toe:-1]-xb[Toe])/float(DuneWidth))+(BermCrest)
            DunePlotEnd=xb[Toe]+3*DuneWidth
            DunePlotSmall=xb[Toe]+DuneCrest
  
            out=num.arange(DunePlotEnd,len(yb),1)
            yb[DunePlotSmall:-1]=yb[DunePlotSmall:-1]/10
            yb[DunePlotSmall:-1]=yb[DunePlotSmall:-1]+yb[DunePlotSmall-1]-yb[DunePlotSmall]

            xb=num.delete(xb,out,None)
            yb=num.delete(yb,out,None); yb[-1]=yb[-2];
  
    # combine all vectors together
    xf=xf+xd[-1];xb=xb+xf[-1] # make one long x-axis
    xd=xd.tolist();yd=yd.tolist() # transform into lists
    xf=xf.tolist();yf=yf.tolist()
    xb=xb.tolist();yb=yb.tolist()
  
    yd.extend(yf);xd.extend(xf) # make one y-axis
    xd.extend(xb);yd.extend(yb)
else:
    TextData=open(r'E:\MarineInVEST\CoastalProtection\Tier1\081711\ProfileBuilder\De.txt',"r") 
    ex=[];Dmeas=[];h=[]

    for line in TextData.readlines():
        linelist = [float(s) for s in line.split("\t")] # split the list by comma delimiter
        ex.append(linelist[0])
        Dmeas.append(linelist[1])
        h.append(linelist[1])

    Dmeas=num.array(Dmeas);ex=num.array(ex);
    
    mylist=abs(ex-ShoreMod1);mylist=mylist.tolist();
    minval=min(mylist)
    sho=[i for i, v in enumerate(mylist) if v == minval];sho=sho[0];

    mylist=abs(ex-OffMod1);mylist=mylist.tolist();
    minval=min(mylist)
    of=[i for i, v in enumerate(mylist) if v == minval];of=of[0];

    m=1./SlopeMod1;
    h[of:sho]=m*ex[of:sho]+Dmeas[of]-m*ex[of]

    if SlopeMod2<>0:
        mylist=abs(ex-ShoreMod2);mylist=mylist.tolist();
        minval=min(mylist)
        sho=[i for i, v in enumerate(mylist) if v == minval];sho=sho[0];

        mylist=abs(ex-OffMod2);mylist=mylist.tolist();
        minval=min(mylist)
        of=[i for i, v in enumerate(mylist) if v == minval];of=of[0];

        m=1./SlopeMod2;
        h[of:sho]=m*ex[of:sho]+Dmeas[of]-m*ex[of]

    if SlopeMod3<>0:
        mylist=abs(ex-ShoreMod3);mylist=mylist.tolist();
        minval=min(mylist)
        sho=[i for i, v in enumerate(mylist) if v == minval];sho=sho[0];

        mylist=abs(ex-OffMod3);mylist=mylist.tolist();
        minval=min(mylist)
        of=[i for i, v in enumerate(mylist) if v == minval];of=of[0];

        m=1./SlopeMod3;
        h[of:sho]=m*ex[of:sho]+Dmeas[of]-m*ex[of]
   
# vegetation variables
gp.AddMessage("\nPlotting Profile...")
if Lveg>0: # fill in seagrass vectors
    temp=num.nonzero(num.array(yd)<=do1v);temp=temp[0];pos1=temp[-1]; # depths shallower than shallowest depth
    Xveg=xd[pos1]-Lveg; #X position where veg. ends
    temp=num.nonzero(num.array(xd)<=Xveg);temp=temp[0];pos2=temp[-1]; # index where veg. ends
    PlntPlot=yd[pos2:pos1+1];Xplant=xd[pos2:pos1+1];do2v=yd[pos2]

# plot and save
if ProfileQuestion <> "(3) No, please assume an equilibrium beach profile": # model extracts value from GIS Layer
    # depth limits for plotting
    VegLow=max([do2v,hc]);# min depth for veg. plot
    AbvHT=num.nonzero(num.array(yd)>HT+2);AbvHT=AbvHT[0]
    if len(AbvHT)==1: VegEnd=xd[AbvHT]
    else: VegEnd=xd[-1]

    AtMSLoc=num.nonzero(num.array(yd)>0);AtMSLoc=AtMSLoc[0];AtMSLoc=AtMSLoc[0]
    AtMSL=xd[AtMSLoc]

    # plot  
    subplot(221)
    plot(num.array(xd),num.array(yd),Xplant,PlntPlot,'-xg');grid();hold
    plot(num.array(xd),num.array(yd)*0,'k',num.array(xd),num.array(yd)*0-MSL,'--k',num.array(xd),num.array(yd)*0+HT,'--k') 
    xlim(Xplant[0]-50,VegEnd);ylim(VegLow-2,HT+2)

    ylabel('Elevation [m]', size='large')

    subplot(222)
    plot(num.array(xd),num.array(yd));grid()
    xlim(AtMSL-10,xd[-1]+10);ylim(-1,max(yd)+2)

    subplot (223)
    plot(num.array(xd),num.array(yd));grid()
    ylim(min(yd)-1,2);

    subplot (224)
    plot(num.array(xd),num.array(yd));grid()
    
else:
    plot(ex,Dmeas,'r');hold;plot(ex,h);grid()
    legend(('Initial Profile','Modified Profile'),'lower right')

# save plot to .PNG
savefig(Profile_Plot, dpi=(640/8))


# create fetch vectors
# copy original point twice and add fields to second copy
gp.CopyFeatures_management(LandPoint, PtsCopy, "", "0", "0", "0")
gp.CopyFeatures_management(LandPoint, PtsCopy2, "", "0", "0", "0")
PtsCopy2 = AddField(PtsCopy2, "DISTANCE", "SHORT", "8", "")
PtsCopy2 = AddField(PtsCopy2, "BEARING", "DOUBLE", "", "")
PtsCopy2 = AddField(PtsCopy2, "BISECTANG", "DOUBLE", "", "")

CopyExpr = PtsCopy
for i in range(0,((BearingsNum*9)+BearingsNum)-2):
    CopyExpr = CopyExpr + ";"+PtsCopy

BiSectAngFullList = []
BiSectAngList = [0.0, 0.15707963267948966, 0.11780972450961724, 0.078539816339744828, 0.039269908169872414, \
                 0.0, 0.039269908169872414, 0.078539816339744828, 0.11780972450961724, 0.15707963267948966, 0.0]

for i in range(0,16):
    for j in range(0,10):
        BiSectAngFullList.append(BiSectAngList[j])
    
gp.Append_management(CopyExpr, PtsCopy2, "NO_TEST", "","")
gp.CalculateField_management(PtsCopy2, "DISTANCE", str(RadLineDist), "PYTHON", "")

# translate information from list into perp transect attribute table
cur = gp.UpdateCursor(PtsCopy2, "", "", "BEARING; FID; BISECTANG")
row = cur.Next()
m = 0
while row:
    FID = float(row.GetValue("FID"))
    Bearing = float((360.000/((BearingsNum*9)+BearingsNum))* FID)
    row.SetValue("Bearing", Bearing)
    row.SetValue("BiSectAng", BiSectAngFullList[m])
    m = m + 1
    cur.UpdateRow(row)
    row = cur.Next()
del cur    
del row

# get the parameters
fc = string.replace(PtsCopy2,"\\","/")
# describe
descfc = gp.describe(fc)
sr = descfc.spatialreference
# process the feature class attributes
lstfc = string.split(fc,"/")
for fl in lstfc:
    fn = fl
strWorkspace = string.replace(fc,fl,"")
gp.workspace = strWorkspace
# shapefile
newfn = string.replace(fl, ".shp", "_lineRotate.shp")
# check for existence
if gp.exists(strWorkspace + newfn):
    gp.delete_management(strWorkspace + newfn )
    gp.refreshcatalog(gp.workspace)
# create the feature class
gp.CreateFeatureClass_management(gp.workspace, newfn, "POLYLINE", fc, "SAME_AS_TEMPLATE", "SAME_AS_TEMPLATE", sr)
addrecs = gp.insertcursor(strWorkspace + newfn)
# refresh the catalog
gp.refreshcatalog(gp.workspace)
  
recs = gp.SearchCursor(fc)
rec = recs.next()
lstFields = gp.listfields(fc)
while rec:
    # get the angle
    rotation = rec.getvalue("BEARING")
    length = rec.getvalue("DISTANCE")
    bearing = math.radians(rotation)
    angle = math.radians((360 - math.degrees(bearing)) + 90)        
    # get the feature and compute the to point
    pt = rec.shape.getpart(0)
    x = operator.add(math.cos(angle) * length, pt.x)
    y = operator.add(math.sin(angle) * length, pt.y)
    # build up the record
    addrec = addrecs.newrow()    
    # create the shape
    newArray = gp.createobject("array")
    newArray.add (pt)
    newPt = gp.createobject("point")
    newPt.x = x
    newPt.y = y
    newArray.add(newPt)
    # maintain the attributes
    lstFields.reset()
    fld = lstFields.next()
    while fld:
        if fld.name <> "FID" and fld.name <> "OBJECTID" and fld.name <> "SHAPE":
            addrec.setvalue(fld.name, rec.getvalue(fld.name))
        fld = lstFields.next()
    # add shape
    addrec.shape = newArray
    addrecs.insertrow(addrec)
    rec = recs.next()

# convert DEM into bathy polygon
gp.Extent = PtsCopyLR

# grab projection spatial reference from 'LandPoly' input
dataDesc = gp.describe(LandPoly)
spatialRef = dataDesc.SpatialReference
gp.CreateFeatureClass_management(interws, "Fetch_AOI.shp", "POLYGON", "#", "#", "#", spatialRef)

# grab four corners from 'PtsCopyLR'
CoordList = shlex.split(gp.Extent)

# when creating a polygon, the coordinates for the starting point must be the same as the coordinates for the ending point
cur = gp.InsertCursor(Fetch_AOI)
row = cur.NewRow()
PolygonArray = gp.CreateObject("Array")
pnt = gp.CreateObject("Point")
pnt.x = float(CoordList[0])
pnt.y = float(CoordList[1])
PolygonArray.add(pnt)
pnt.x = float(CoordList[0])
pnt.y = float(CoordList[3])
PolygonArray.add(pnt)
pnt.x = float(CoordList[2])
pnt.y = float(CoordList[3])
PolygonArray.add(pnt)
pnt.x = float(CoordList[2])
pnt.y = float(CoordList[1])
PolygonArray.add(pnt)
pnt.x = float(CoordList[0])
pnt.y = float(CoordList[1])
PolygonArray.add(pnt)
row.shape = PolygonArray
cur.InsertRow(row)
del row, cur

# erase from 'Fetch_AOI' areas where there is land
LandPoly = AddField(LandPoly, "ERASE", "SHORT", "0", "0")
gp.CalculateField_management(LandPoly, "ERASE", "1", "VB")
UnionExpr = Fetch_AOI+" 1; "+LandPoly+" 2"        
gp.Union_analysis(UnionExpr, UnionFC)

# select features where "ERASE = 0"
gp.Select_analysis(UnionFC, SeaPoly, "\"ERASE\" = 0")

# erase parts of line where it overlaps land (works for ArcView)
gp.Intersect_analysis(PtsCopyLR+" 1;"+SeaPoly+" 2", PtsCopyEL, "ALL", "", "INPUT")
gp.MultipartToSinglepart_management(PtsCopyEL, PtsCopyExp)
# convert to layer to select only lines originating from point source
gp.MakeFeatureLayer_management(PtsCopyExp, PtsCopyExp_Lyr, "", gp.workspace, "")
gp.SelectLayerByLocation_management(PtsCopyExp_Lyr, "WITHIN_A_DISTANCE", LandPoint, "20 Meters", "NEW_SELECTION")
gp.CopyFeatures_management(PtsCopyExp_Lyr, Fetch_Vectors, "", "0", "0", "0")
# add and calculate "LENGTH_M" field
Fetch_Vectors = AddField(Fetch_Vectors, "LENGTH_M", "LONG", "6", "")
gp.CalculateField_management(Fetch_Vectors, "LENGTH_M", "!shape.length@meters!", "PYTHON", "")

### plot fetch on rose
##radians = (num.pi / 180.0)
##pi = num.pi
##theta16 = [0*radians,22.5*radians,45*radians,67.5*radians,90*radians,112.5*radians,135*radians,157.5*radians,180*radians,202.5*radians,225*radians,247.5*radians,270*radians,292.5*radians,315*radians,337.5*radians]
##rc('grid', color='#316931', linewidth=1, linestyle='-')
##rc('xtick', labelsize=0)
##rc('ytick', labelsize=15)
### force square figure and square axes looks better for polar, IMO
##width, height = matplotlib.rcParams['figure.figsize']
##size = min(width, height)
### make a square figure
##plt = figure(figsize=(size, size))
##ax = plt.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='w')
### plot
##bars = ax.bar(theta16, FetchList, width=.35, color='#ee8d18', lw=1)
##for r,bar in zip(FetchList, bars):
##    bar.set_facecolor( cm.YlOrRd(r/10.))
##    bar.set_alpha(.65)
##ax.set_rmax(max(FetchList)+1)
##grid(True)
##ax.set_title("Average Fetch (meters)", fontsize=15)
##plt.savefig(outputws+"Fetch_Plot.png", dpi=(640/8))


# create cost surface based on 'SeaPoly'
gp.Extent = Fetch_AOI
projection = grabProjection(LandPoint)
gp.Project_management(WW3_Pts, WW3_Pts_prj, projection)
SeaPoly = AddField(SeaPoly, "SEA", "SHORT", "", "")
gp.CalculateField_management(SeaPoly, "SEA", "1", "PYTHON", "")
gp.FeatureToRaster_conversion(SeaPoly, "SEA", seapoly_rst, "250")
gp.Expand_sa(seapoly_rst, seapoly_e, "1", "1")
# allocate 'WW3_Pts' throughout cost surface
gp.CostAllocation_sa(WW3_Pts_prj, seapoly_e, costa_ww3, "", "", "FID", "", "")
# determine which point is closest to 'LandPoint'
gp.ExtractValuesToPoints_sa(LandPoint, costa_ww3, LandPoint_WW3, "NONE")
cur = gp.UpdateCursor(LandPoint_WW3)
row = cur.Next()
WW3_FID = row.GetValue("RASTERVALU")
del row
del cur

# populate list with data from closest WW3 point
WW3_ValuesList = []
dirList = [0, 22, 45, 67, 90, 112, 135, 157, 180, 202, 225, 247, 270, 292, 315, 337]
letterList = ['c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t']
SrchCondition = "FID = "+str(WW3_FID)
cur = gp.SearchCursor(WW3_Pts_prj, SrchCondition, "", "")
row = cur.Next()
WW3_ValuesList.append(row.GetValue("LAT")) # 0
WW3_ValuesList.append(row.GetValue("LONG")) # 1
for i in range(0,len(dirList)):
    WW3_ValuesList.append(row.GetValue("V10PCT_"+str(dirList[i]))) # 2 - 17
for i in range(0,len(dirList)):
    WW3_ValuesList.append(row.GetValue("V25PCT_"+str(dirList[i]))) # 18 - 33
for i in range(0,len(dirList)):
    WW3_ValuesList.append(row.GetValue("V_MAX_"+str(dirList[i]))) # 34 - 49
WW3_ValuesList.append(row.GetValue("V_10YR")) # 50
WW3_ValuesList.append(row.GetValue("H_10PCT")) # 51
WW3_ValuesList.append(row.GetValue("T_10PCT")) # 52
WW3_ValuesList.append(row.GetValue("H_25PCT")) # 53
WW3_ValuesList.append(row.GetValue("T_25PCT")) # 54
WW3_ValuesList.append(row.GetValue("H_MAX")) # 55
WW3_ValuesList.append(row.GetValue("T_MAX")) # 56
WW3_ValuesList.append(row.GetValue("H_10YR")) # 57
WW3_ValuesList.append(row.GetValue("He")) # 58
WW3_ValuesList.append(row.GetValue("Hmod")) # 59
WW3_ValuesList.append(row.GetValue("Tmod")) # 60
del row
del cur

# import Profile Builder info from Excel file
xlApp = Dispatch("Excel.Application")
xlApp.Visible=0
xlApp.DisplayAlerts=0
xlApp.Workbooks.Open(InputTable)
# write WW3 results to Excel sheet 'Erosion Model Input'
cell2 = xlApp.Worksheets("Erosion Model Input")
# maximum wave height
cell2.Range("e85").Value = WW3_ValuesList[55]
cell2.Range("f85").Value = WW3_ValuesList[56]
# top 10% wave height
cell2.Range("e86").Value = WW3_ValuesList[51]
cell2.Range("f86").Value = WW3_ValuesList[52]
# top 25% wave height
cell2.Range("e87").Value = WW3_ValuesList[53]
cell2.Range("f87").Value = WW3_ValuesList[54]
# 10-yr wave height
cell2.Range("e88").Value = WW3_ValuesList[57]
# maximum wind speed
for i in range(34,50):
    cell2.Range(letterList[i-34]+"89").Value = WW3_ValuesList[i]
# top 10% wind speed
for i in range(2,18):
    cell2.Range(letterList[i-2]+"90").Value = WW3_ValuesList[i]
# top 25% wind speed
for i in range(18,34):
    cell2.Range(letterList[i-18]+"91").Value = WW3_ValuesList[i]
# 10-yr wind speed
cell2.Range("e92").Value = WW3_ValuesList[50]
xlApp.ActiveWorkbook.Close(SaveChanges=1) # save changes
xlApp.Quit()


# return projected point to geographic (unprojected)
gp.Project_management(LandPoint, LandPoint_Geo, geo_projection)
# grab coordinates for Google Maps plot
cur = gp.UpdateCursor(LandPoint_Geo)
row = cur.Next()
feat = row.Shape
midpoint1 = feat.Centroid
midList1 = shlex.split(midpoint1)
midList1 = [float(s) for s in midList1]
del cur
del row
PtLat = str(midList1[1])
PtLong = str(midList1[0])

# create html file
htmlfile = open(Profile_HTML, "w")
htmlfile.write("<html>\n")
htmlfile.write("<title>Marine InVEST</title>")
htmlfile.write("<CENTER><H1>Visualizing Coastal Protection - Tier 1</H1></CENTER>")
htmlfile.write("<br><HR><H2>Map and Plots</H2>\n")
htmlfile.write("This is a map and plots showing the location and characteristics of where the Profile Builder tool was run. <br>\n")
htmlfile.write("<table border=\"0\"><tr><td>")
htmlfile.write("<iframe width=\"640\" height=\"640\" frameborder=\"0\" scrolling=\"no\" marginheight=\"0\" marginwidth=\"0\"") 
htmlfile.write("src=\"http://maps.google.com/maps/api/staticmap?center=")
htmlfile.write(PtLat+","+PtLong)
htmlfile.write("&zoom=11&size=640x640&maptype=hybrid&markers=color:red%7Ccolor:red%7Clabel:X%7C")
htmlfile.write(PtLat+","+PtLong)
htmlfile.write("&sensor=false\"></iframe><br/><small><a href=\"http://maps.google.com/maps?f=q&amp;source=embed&amp;hl=en&amp;geocode=&amp;q=")
htmlfile.write(PtLat+","+PtLong)
htmlfile.write("&amp;aq=&amp;")
htmlfile.write("sll=37.160317,-95.712891&amp;sspn=48.113934,71.455078&amp;ie=UTF8&amp;z=10&amp;ll=")
htmlfile.write(PtLat+","+PtLong)
htmlfile.write("\" style=\"color:#0000FF;text-align:center\">View Larger Map</a></small><br>\n")
htmlfile.write("</td><td>")
htmlfile.write("<img src=\"Profile_Plot.png\" width=\"640\" height=\"480\">")
htmlfile.write("</td></tr></table>")
htmlfile.write("<img src=\"Profile_Plot.png\" width=\"640\" height=\"480\">")
htmlfile.write("<br>\n")
htmlfile.write("<br><HR><H2>Site Information</H2>\n")
htmlfile.write("<li><u>The site is located at</u> - Latitude: "+PtLat+", Longitude: "+PtLong+"<br>\n")
htmlfile.write("<li><u>The tidal range is</u>: xxx m (High Tide value)<br>\n")
htmlfile.write("<li><u>The foreshore slope is</u>: xxx<br>\n")
htmlfile.write("<li><u>The backshore has a slope that is</u>: xxx m high and yyy m long<br>\n")
htmlfile.write("<li><u>The beach is backed by a sand dune that is</u>: xxx m high<br>\n")
htmlfile.write("<li>There is vegetation in the sub- and inter-tidal area.  The vegetation characteristics are xxx.<br>\n")
htmlfile.write("</html>")
htmlfile.close()

# create parameter file
parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
parafile = open(outputws+"parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
parafile.writelines("PROFILE BUILDER PARAMETERS\n")
parafile.writelines("___________________________________________\n\n")
     
for para in parameters:
    parafile.writelines(para+"\n")
    parafile.writelines("\n")
parafile.close()