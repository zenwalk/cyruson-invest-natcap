# Marine InVEST: Coastal Protection (Profile Generator)
# Authors: Greg Guannel, Gregg Verutes, Jeremy Davies
# 10/04/12

# import libraries
import CPf_SignalSmooth as SignalSmooth
import string, sys, os, time, datetime, shlex, shutil
import fpformat, operator
import arcgisscripting
from math import *

# create the geoprocessor object
gp=arcgisscripting.create()
gp.AddMessage("\nChecking and preparing inputs...")

# set output handling
gp.OverwriteOutput=1
# check out any necessary extensions
gp.CheckOutExtension("management")
gp.CheckOutExtension("analysis")
gp.CheckOutExtension("conversion")
gp.CheckOutExtension("spatial")

# error messages
msgArguments="\nProblem with arguments."
msgCheckInputs = "\nError checking and preparing inputs."
msgPrepWW3Fetch = "\nError preparing Wave Watch III inputs or calculating fetch."
msgReadExcel = "\nError reading Profile Generator Excel file inputs."
msgCreateProfile = "\nError creating profile."
msgModifyProfile = "\nError modifying profile."
msgPlotProfile = "\nError plotting profile."
msgHTMLOutputs = "\nError generating HTML outputs."
msgNumPyNo = "NumPy extension is required to run the Coastal Protection Model.  Please consult the Marine InVEST FAQ for instructions on how to install."
msgSciPyNo = "SciPy extension is required to run the Coastal Protection Model.  Please consult the Marine InVEST FAQ for instructions on how to install."
msgWin32ComNo = "PythonWin extension is required to run the Coastal Protection Model.  Please consult the Marine InVEST FAQ document for instructions on how to install."
msgMatplotlibNo = "Matplotlib extension (version 1.0 or newer) is required to run the Coastal Protection Model.  Please consult the Marine InVEST FAQ document for instructions on how to install."

# import modules
try:
    import numpy as num
except:
    gp.AddError(msgNumPyNo)
    raise Exception

try:
    from scipy.interpolate import interp1d
    from scipy import optimize
except:
    gp.AddError(msgSciPyNo)
    raise Exception

try:
    from win32com.client import Dispatch
except:
    gp.AddError(msgWin32ComNo)
    raise Exception

try:
    from matplotlib import *
    from pylab import *
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
except:
    gp.AddError(msgMatplotlibNo)
    raise Exception

try:
    try:
        # get parameters
        parameters=[]
        now=datetime.datetime.now()
        parameters.append("Date and Time: "+now.strftime("%Y-%m-%d %H:%M"))
        gp.workspace=gp.GetParameterAsText(0)
        parameters.append("Workspace: "+gp.workspace)
        subwsStr=gp.GetParameterAsText(1)
        parameters.append("Label for Profile Generator Run (10 characters max): "+subwsStr)
        LandPoint=gp.GetParameterAsText(2)
        parameters.append("Land Point: "+LandPoint)
        LandPoly=gp.GetParameterAsText(3)
        parameters.append("Land Polygon: "+LandPoly)
        ProfileQuestion=gp.GetParameterAsText(4)
        parameters.append("Do you want us to cut a cross-shore transect in GIS?: "+ProfileQuestion)
        BathyGrid=gp.GetParameterAsText(5)
        parameters.append("IF 1: Bathymetric Grid (DEM): "+BathyGrid)
        HabDirectory=gp.GetParameterAsText(6)
        parameters.append("IF 1: Habitat Data Directory: "+HabDirectory)
        BufferDist=gp.GetParameterAsText(7)
        parameters.append("IF 1: Land Point Buffer Distance: "+BufferDist)
        RadLineDist=gp.GetParameterAsText(8)
        parameters.append("IF 1: Length of your profile [km] "+RadLineDist)
        CSProfile=gp.GetParameterAsText(9)
        parameters.append("IF 2: Upload Your Cross-Shore Profile: "+CSProfile)
        SmoothParameter=float(gp.GetParameterAsText(10))
        parameters.append("Smoothing Percentage (Value of '0' means no smoothing): "+str(SmoothParameter))
        InputTable=gp.GetParameterAsText(11)
        parameters.append("Nearshore Waves Excel Table: "+InputTable)
        WW3_Pts=gp.GetParameterAsText(12)
        parameters.append("Wave Watch III Model Data: "+WW3_Pts)
        WW3_SearchDist=int(gp.GetParameterAsText(13))
        parameters.append("Wave Watch III Search Distance: "+str(WW3_SearchDist))
        FetchQuestion=gp.GetParameterAsText(14)
        parameters.append("Do you wish to calculate fetch for LandPoint?: "+FetchQuestion)

    except:
        raise Exception,msgArguments+gp.GetMessages(2)


    ###############################################           
    ###### CREATE DIRECTORIES AND VARIABLES #######
    ###############################################

    # remove spaces and shorten 'subwsStr' if greater than 10 characters
    subwsStr=subwsStr.replace(" ","")
    subwsStr=subwsStr[0:10]

    # intermediate and output directories
    interws=gp.GetParameterAsText(0)+os.sep+"scratch"+os.sep
    outputws=gp.GetParameterAsText(0)+os.sep+"_ProfileGenerator_Outputs"+os.sep
    gp.scratchworkspace=interws
    subws=outputws+subwsStr+os.sep
    maps=subws+"maps"+os.sep
    html_txt=subws+"html_txt"+os.sep

    try:
        thefolders=["scratch","_ProfileGenerator_Outputs"]
        for folder in thefolders:
            if not gp.exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace,folder)

        thefolders=[subwsStr]
        for folder in thefolders:
            if not gp.exists(outputws+folder):
                gp.CreateFolder_management(outputws,folder)  

        thefolders=["maps","html_txt"]
        for folder in thefolders:
            if not gp.exists(subws+folder):
                gp.CreateFolder_management(subws,folder)  
    except:
        raise Exception,"Error creating folders"

    LandPointLyr=interws+"LandPointLyr.lyr"
    LandPolyLyr=interws+"LandPolyLyr.lyr"
    PT1=interws+"PT1.shp"
    PT2=interws+"PT2.shp"
    PT1_Z=interws+"PT1_Z.shp"
    PT2_Z=interws+"PT2_Z.shp"
    PT_Z_Near=interws+"PT_Z_Near.shp"
    Backshore_Pts=interws+"Backshore_Pts.shp"
    Profile_Pts_Merge=interws+"Profile_Pts_Merge.shp"
    Profile_Pts_Lyr=interws+"Profile_Pts_Lyr.lyr"
    LandPoint_Buff=interws+"LandPoint_Buff.shp"
    LandPoint_Buff50k=interws+"LandPoint_Buff50k.shp"
    LandPoint_Geo=interws+"LandPoint_Geo.shp"
    Shoreline=interws+"Shoreline.shp"
    Shoreline_Buff_Clip=interws+"Shoreline_Buff_Clip.shp"
    Shoreline_Buff_Clip_Diss=interws+"Shoreline_Buff_Clip_Diss.shp"
    PtsCopy=interws+"PtsCopy.shp"
    PtsCopy2=interws+"PtsCopy2.shp"
    PtsCopyLR=interws+"PtsCopy2_lineRotate.shp"
    PtsCopy3=interws+"PtsCopy3.shp"
    PtsCopyLD=interws+"PtsCopy3_lineDist.shp"
    Fetch_AOI=interws+"Fetch_AOI.shp"
    UnionFC=interws+"UnionFC.shp"
    SeaPoly=interws+"SeaPoly.shp"
    seapoly_rst=interws+"seapoly_rst"
    seapoly_e=interws+"seapoly_e"
    PtsCopyEL=interws+"PtsCopy2_eraseLand.shp"
    PtsCopyExp=interws+"PtsCopy2_explode.shp"
    PtsCopyExp_Lyr=interws+"PtsCopy2_explode.lyr"
    WW3_Pts_prj=interws+"WW3_Pts_prj.shp"
    costa_ww3=interws+"costa_ww3"
    LandPoint_WW3=interws+"LandPoint_WW3.shp"

    BathyProfile=html_txt+"BathyProfile_"+subwsStr+".txt"
    CreatedProfile=html_txt+"CreatedProfile_"+subwsStr+".txt"
    ProfileCutGIS=html_txt+"ProfileCutGIS_"+ subwsStr+".txt"
    FetchDistances=html_txt+"FetchDistances_"+ subwsStr+".txt"
    HabitatLocation=html_txt+"HabitatLocation_"+ subwsStr+".txt"
    Profile_HTML=html_txt+"profile.html"
    FetchWindWave_HTML=html_txt+"fetchwindwave.html"
    Wind_Plot=html_txt+"Wind_Plot.png"
    Fetch_Plot=html_txt+"Fetch_Plot.png"

    Profile_Pts=maps+"Profile_Pts.shp"
    Profile_Pts_Hab=maps+"Profile_Pts_Hab.shp"
    Fetch_Vectors=maps+"Fetch_Vectors.shp"
    Fetch_Distances=maps+"Fetch_Distances.shp"
    UploadedProfile=maps+"UploadedProfile.shp"


    ################################         
    ###### VARIOUS FUNCTIONS #######
    ################################

    def AddField(FileName,FieldName,Type,Precision,Scale):
        fields=gp.ListFields(FileName,FieldName)
        field_found=fields.Next()
        if field_found:
            gp.DeleteField_management(FileName,FieldName)
        gp.AddField_management(FileName,FieldName,Type,Precision,Scale,"","","NON_NULLABLE","NON_REQUIRED","")
        return FileName

    def getDatum(thedata):
        desc=gp.describe(thedata)
        SR=desc.SpatialReference
        if SR.Type=="Geographic":
            strDatum=SR.DatumName         
        else:
            gp.OutputCoordinateSystem=SR
            strSR=str(gp.OutputCoordinateSystem)
            gp.OutputCoordinateSystem=""
            n1=strSR.find("GEOGCS")
            n2=strSR.find("PROJECTION")
            strDatum=strSR[n1:n2-1]
        return strDatum

    def ckDatum(thedata):
        desc=gp.describe(thedata)
        SR=desc.SpatialReference
        if SR.Type=="Geographic":
            strDatum=SR.DatumName         
        else:
            gp.OutputCoordinateSystem=SR
            strSR=str(gp.OutputCoordinateSystem)
            gp.OutputCoordinateSystem=""
            n1=strSR.find("DATUM[\'")
            n2=strSR.find("\'",n1+7)
            strDatum=strSR[n1+7:n2]
        if strDatum=="D_WGS_1984":
            pass
        else:
            gp.AddError(thedata+" is not a valid input.\nThe model requires data inputs and a projection with the \"WGS84\" datum.\nSee InVEST FAQ document for how to reproject datasets.")
            raise Exception

    def ckProjection(data):
        dataDesc=gp.describe(data)
        spatreflc=dataDesc.SpatialReference
        if spatreflc.Type <> 'Projected':
            gp.AddError(data+" does not appear to be projected.  It is assumed to be in meters.")
            raise Exception
        if spatreflc.LinearUnitName <> 'Meter':
            gp.AddError("This model assumes that "+data+" is projected in meters for area calculations.  You may get erroneous results.")
            raise Exception

    def checkGeometry(thedata,Type,Message):
        if gp.Describe(thedata).ShapeType <> Type:
            raise Exception,"\nInvalid input: "+thedata+"\n"+Message+" must be of geometry type "+Type+"."

    def checkInteger(thedata):
        if thedata.find("0")==-1 and thedata.find("1")==-1 and thedata.find("2")==-1 and thedata.find("3")==-1 and thedata.find("4")==-1 and thedata.find("5")==-1 and thedata.find("6")==-1 and thedata.find("7")==-1 and thedata.find("8")==-1 and thedata.find("9")==-1:
            gp.AddError(thedata +" must contain an underscore followed by an integer ID at the end of it's name (e.g. filename_1.shp). This is necessary to properly link it with the input table.")
            raise Exception

    def grabProjection(data):
        dataDesc=gp.describe(data)
        sr=dataDesc.SpatialReference
        gp.OutputCoordinateSystem=sr
        strSR=str(gp.OutputCoordinateSystem)
        return strSR

    def PTCreate(PTType,midx,midy,TransectDist): # function to create point transects
        if PTType==1:
            y1=midy+TransectDist
            y2=midy-TransectDist
            x1=midx
            x2=midx
        elif PTType==2:
            y1=midy
            y2=midy 
            x1=midx+TransectDist
            x2=midx-TransectDist
        elif PTType==3:
            y1=NegRecip*(TransectDist)+midy
            y2=NegRecip*(-TransectDist)+midy
            x1=midx+TransectDist
            x2=midx-TransectDist
        elif PTType==4:
            y1=midy+TransectDist
            y2=midy-TransectDist
            x1=(TransectDist/NegRecip)+midx
            x2=(-TransectDist/NegRecip)+midx
        elif PTType==5:
            y1=midy+TransectDist
            y2=midy-TransectDist
            x1=(TransectDist/NegRecip)+midx
            x2=(-TransectDist/NegRecip)+midx
        elif PTType==6:
            y1=NegRecip*(TransectDist)+midy
            y2=NegRecip*(-TransectDist)+midy
            x1=midx+TransectDist
            x2=midx-TransectDist
        return x1,y1,x2,y2

    def compareProjections(LandPoint,LandPoly):
        if gp.describe(LandPoint).SpatialReference.name <> gp.describe(LandPoly).SpatialReference.name:
            gp.AddError("Projection Error: "+LandPoint+" is in a different projection from the LandPoly data.  \
                         The two inputs must be the same projection to calculate depth profile.")
            raise Exception
        
    def LocateHabitat(X,HAB):
        # generate a vector showing location of beg. and end of patches of a certain habitat along the X axis of bathy
        # input: X axis assigned with bathy and HAB vector of same length with 0's and 1's to locate hab
        # output: vector with X location of beg. and end of patches
    
        HabPst=num.nonzero(HAB);HabPst=HabPst[0] # locate non-zero values in HAB; extent of habitation on profile
        HabGap=num.nonzero(diff(HabPst)<>1);HabGap=HabGap[0]  # locate where there's a gap in the habitat     
        if len(HabGap)>0:
            # indices of where patches begin and end
            HabPst1=num.append(HabPst[0],HabPst[HabGap],None)
            HabPst1=num.append(HabPst1,HabPst[HabGap+1],None)
            HabPst1=num.append(HabPst1,HabPst[-1],None)
            HabPst1=sort(HabPst1) # even (odd) numbers indicate loc. of beg (end) of patches
            habx=X[HabPst1]
            
            # empty array for X Location of beg. and end of each patch
            begHABx=num.arange(0,len(HabGap)+1,1)*0
            finHABx=num.arange(0,len(HabGap)+1,1)*0
            
            # fill in those array when there are many patches
            for kk in range(len(habx)/2):
                begHABx[kk]=habx[2*kk]
                finHABx[kk]=habx[2*kk+1]
                
        elif len(HabPst)>0: # no patches, just ID beg and end
            begHABx=[X[HabPst[0]]];
            finHABx=[X[HabPst[-1]]];
        else: # no habitat -> return empty string
            begHABx=[];finHABx=[];
        begHABx=num.array(begHABx);finHABx=num.array(finHABx)
        
        return begHABx,finHABx
    
    def Indexed(x,value): # locates index of point in vector x that has closest value as variable value
        mylist=abs(x-value)    
        if isinstance(x,num.ndarray):
            mylist=mylist.tolist()
        minval=min(mylist)
        ind=[i for i,v in enumerate(mylist) if v==minval]
        ind=ind[0]
        return ind

    def SlopeModif(X,Y,SlopeMod,Offx,Shorex): 
        # this function replaces/adds linear portion to an existing profile 
        # inputs are values of bathy profile along and X axis, as well as information for the linear approx: slope, offshore and shoreward extent
        # outputs are a new bathy profile that includes the new linear slope
    
        # special case when we start with one point: make it two
        if len(X)==1:
            X=[X[0],X[0]+1]
            Y=[Y[0],Y[0]-.5]
        
        X=num.array(X)
        Y=num.array(Y)
        # store initial values of X and Y
        Xo=X
        Yo=Y
        F=interp1d(Xo,Yo)    
    
        # slope 
        if SlopeMod <> 0:
            m=-1.0/SlopeMod # slope
        else:
            m=0.0
        
        # locate footprint of change
        ShorePos=Indexed(X,Shorex) # locate position of shoreward point for linear approx
        if X[ShorePos]>Shorex: # if shoreward bound. outside of existing X vector, add more points
            X=num.arange(Shorex,X[-1]+1,1)
            ShorePos=0 # new start point of linear approx
            
        OffPos=Indexed(X,Offx) # locate position of offshore point for linear approx
        if X[OffPos]<Offx: # if new offshore bound. outside of existing, add more points
            X=num.arange(X[0],Offx+1,1)
            OffPos=len(X)-1 # new end point of linear approx
        
        # create new bathy profile
        Y=X*0+Yo[-1] # new vector Y that will be used to create the linear approx
        temp1=Indexed(X,Xo[0])
        temp2=Indexed(X,Xo[-1])
        Y[temp1:temp2+1]=F(X[temp1:temp2+1]) # map old values of X on new vector Y
        Y[ShorePos:OffPos+1]=m*(X[ShorePos:OffPos+1]-X[OffPos])+Y[OffPos] # create linear approximation in range that user defined
    
        return X,Y

    # wind-wave generation
    def WindWave(U,F,d):
        ds=g*d/U**2.0;Fs=g*F/U**2.0
        A=tanh(0.343*ds**1.14)
        B=tanh(4.14e-4*Fs**0.79/A)
        H=0.24*U**2/g*(A*B)**0.572 # wave height
        A=tanh(0.1*ds**2.01)
        B=tanh(2.77e-7*Fs**1.45/A)
        T=7.69*U/g*(A*B)**0.18 # wave period
        return H,T

    ######################################           
    ###### VARIOUS CHECK TO INPUTS #######
    ######################################

    try:
        # variables (hard-coded)
        SampInterval=1
        TransectDist=0
        BearingsNum=16
        RadLineDist=2.0*int(RadLineDist)*1000

        # check that correct inputs were provided based on 'ProfileQuestion'
        if ProfileQuestion=="(1) Yes":
            if not BathyGrid:
                gp.AddError("A bathymetry grid input is required to create a point transect.")
                raise Exception

        elif ProfileQuestion=="(2) No, but I will upload a cross-shore profile":
            if not CSProfile:
                gp.AddError("A cross-shore profile input is required.")
                raise Exception
            # provide message about not being able to calculate WW3 and fetch if user doesn't provide LandPoint 
            if WW3_Pts or FetchQuestion=='(1) Yes':
                gp.AddWarning("Unable to gather wind and wave statistics using Wave Watch III and/or calculate fetch because you have chosen NOT to cut a cross-shore transect in GIS.")
                WW3_Pts=''
                FetchQuestion='(2) No'

        # check that datum is WGS84
        # need it to be WGS84 datum because Arc geoprocesing won't do transformations on-the-fly
        ckDatum(LandPoint) 

        # check that 'LandPoint' only has one feature in it
        if gp.GetCount_management(LandPoint) <> 1:
            gp.AddError("'Land Point' input should contain only one feature (point).")
            raise Exception
            
        # limit 'BufferDist'
        if int(BufferDist) > 500 or int(BufferDist) < 40:
            gp.AddError("Buffer distance for 'Land Point' must be greater than 40 and less than 500 meters.")
            raise Exception

        # bounds for 'WW3_SearchDist'
        if WW3_SearchDist < 50 or WW3_SearchDist > 500:
            gp.AddError("Wave Watch III Search Distance must be greater than 50 and less than 500 kilometers.")
            raise Exception        

        # check that 'LandPoint' is within 150 meters of coastline
        if ProfileQuestion=="(1) Yes":
            gp.MakeFeatureLayer_management(LandPoint,LandPointLyr,"","","")
            gp.MakeFeatureLayer_management(LandPoly,LandPolyLyr,"","","")
            PointSelect=gp.SelectLayerByLocation_management(LandPointLyr,"WITHIN_A_DISTANCE",LandPolyLyr,str(int(BufferDist)+150)+" Meters","NEW_SELECTION")
            if gp.GetCount_management(PointSelect)==0:
                gp.AddError("Shoreline was not found within "+str(BufferDist)+" meters of 'Land Point' input.\nEither increase the buffer distance or move the 'LandPoint' input closer to the coastline.")
                raise Exception

        # check that three inputs are projected
        ckProjection(LandPoint)
        ckProjection(LandPoly)
        if BathyGrid:
            ckProjection(BathyGrid)
        geo_projection=getDatum(LandPoint) # get datum of 'Land Point'

        # check that smoothing input is within proper range
        if SmoothParameter < 0.0 or SmoothParameter > 100.0:
            gp.AddError("Smoothing parameter value must be between 0 and 100.  Please modify your input and run the tool again.")
            raise Exception
        
        # angle direction list for fetch and WW3
        dirList=[0,22,45,67,90,112,135,157,180,202,225,247,270,292,315,337]

        # modal wave height in case user doesn't enter
        Hm=-1; Tm=-1

        # buffer 'LandPoint' by 50km
        gp.Buffer_analysis(LandPoint,LandPoint_Buff50k,str(WW3_SearchDist*1000)+" Meters","FULL","ROUND","NONE","")
        # convert buffered 'LandPoint' into bathy polygon
        gp.Extent=LandPoint_Buff50k
        # grab projection spatial reference from 'Land Poly' input
        dataDesc=gp.describe(LandPoly)
        spatialRef=dataDesc.SpatialReference
        gp.CreateFeatureClass_management(interws,"Fetch_AOI.shp","POLYGON","#","#","#",spatialRef)
        # grab four corners from 'PtsCopyLR'
        CoordList=shlex.split(gp.Extent)

        # when creating a polygon,the coordinates for the starting point must be the same as the coordinates for the ending point
        cur=gp.InsertCursor(Fetch_AOI)
        row=cur.NewRow()
        PolygonArray=gp.CreateObject("Array")
        pnt=gp.CreateObject("Point")
        pnt.x=float(CoordList[0])
        pnt.y=float(CoordList[1])
        PolygonArray.add(pnt)
        pnt.x=float(CoordList[0])
        pnt.y=float(CoordList[3])
        PolygonArray.add(pnt)
        pnt.x=float(CoordList[2])
        pnt.y=float(CoordList[3])
        PolygonArray.add(pnt)
        pnt.x=float(CoordList[2])
        pnt.y=float(CoordList[1])
        PolygonArray.add(pnt)
        pnt.x=float(CoordList[0])
        pnt.y=float(CoordList[1])
        PolygonArray.add(pnt)
        row.shape=PolygonArray
        cur.InsertRow(row)
        del row,cur

    except:
        gp.AddError(msgCheckInputs)
        raise Exception


    ##############################################           
    ###### WW3, FETCH AND AOI CALCULATIONS #######
    ##############################################

    try:
        if (WW3_Pts or FetchQuestion=='(1) Yes') and ProfileQuestion=="(1) Yes":
            gp.AddMessage("\nPreparing inputs for Wave Watch III and/or fetch calculations...")
            # erase from 'Fetch_AOI' areas where there is land
            LandPoly=AddField(LandPoly,"ERASE","SHORT","0","0")
            gp.CalculateField_management(LandPoly,"ERASE","1","VB")
            UnionExpr=Fetch_AOI+" 1; "+LandPoly+" 2"        
            gp.Union_analysis(UnionExpr,UnionFC)

            # select features where "ERASE=0"
            gp.Select_analysis(UnionFC,SeaPoly,"\"ERASE\"=0")

        # read WW3 Info
        if WW3_Pts:
            gp.AddMessage("...reading Wave Watch III information")
            # create cost surface based on 'SeaPoly'
            gp.Extent=Fetch_AOI
            projection=grabProjection(LandPoint)
            gp.Project_management(WW3_Pts,WW3_Pts_prj,projection)
            WW3_Pts_prj=AddField(WW3_Pts_prj,"PT_ID","LONG","","")
            gp.CalculateField_management(WW3_Pts_prj,"PT_ID","!FID! + 1","PYTHON")

            if gp.GetCount(WW3_Pts_prj) == 0:
                gp.AddError("No Wave Watch III points found within your AOI.  Please increase the 'WW3 Search Distance' input.")
                raise Exception
            elif gp.GetCount(WW3_Pts_prj) == 1:
                cur=gp.UpdateCursor(WW3_Pts_prj)
                row=cur.Next()
                WW3_PT_ID=row.GetValue("PT_ID")
                del row, cur
            else:
                SeaPoly=AddField(SeaPoly,"SEA","SHORT","","")
                gp.CalculateField_management(SeaPoly,"SEA","1","PYTHON","")
                gp.FeatureToRaster_conversion(SeaPoly,"SEA",seapoly_rst,"250")
                gp.Expand_sa(seapoly_rst,seapoly_e,"1","1")
                
                # allocate 'WW3_Pts' throughout cost surface
                gp.CostAllocation_sa(WW3_Pts_prj,seapoly_e,costa_ww3,"","","PT_ID","","")
                # determine which point is closest to 'LandPoint'
                gp.ExtractValuesToPoints_sa(LandPoint,costa_ww3,LandPoint_WW3,"NONE")
                cur=gp.UpdateCursor(LandPoint_WW3)
                row=cur.Next()
                WW3_PT_ID=row.GetValue("RASTERVALU")
                del row, cur

            # populate list with data from closest WW3 point
            WW3_ValuesList=[]
            SrchCondition="PT_ID = "+str(WW3_PT_ID)
            cur=gp.SearchCursor(WW3_Pts_prj,SrchCondition,"","")
            row=cur.Next()            
            WW3_ValuesList.append(row.GetValue("LAT")) # 0
            WW3_ValuesList.append(row.GetValue("LONG_")) # 1
            for i in range(0,len(dirList)):
                WW3_ValuesList.append(row.GetValue("V10PCT_"+str(dirList[i]))) # 2-17
            for i in range(0,len(dirList)):
                WW3_ValuesList.append(row.GetValue("V25PCT_"+str(dirList[i]))) # 18-33
            for i in range(0,len(dirList)):
                WW3_ValuesList.append(row.GetValue("V_MAX_"+str(dirList[i]))) # 34-49
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
            del row, cur
            
            # maximum wave height
            WavMax=[0,0];WavMax[0]=WW3_ValuesList[55];WavMax[1]=WW3_ValuesList[56]
            # top 10% wave height
            Wav10=[0,0];Wav10[0]=WW3_ValuesList[51];Wav10[1]=WW3_ValuesList[52]
            # top 25% wave height
            Wav25=[0,0];Wav25[0]=WW3_ValuesList[53];Wav25[1]=WW3_ValuesList[54]
            # 10-yr wave height
            Wav10yr=WW3_ValuesList[57]
            # maximum wind speed
            WiMax=num.arange(0,16,1)*0      
            for ii in range(0,16):
                WiMax[ii]=WW3_ValuesList[ii+34]
            # top 10% wind speed
            Wi10=num.arange(0,16,1)*0       
            for ii in range(0,16):
                Wi10[ii]=WW3_ValuesList[ii+2]
            # top 25% wind speed
            Wi25=num.arange(0,16,1)*0       
            for ii in range(0,16):
                Wi25[ii]=WW3_ValuesList[ii+18]
            # Hmod and Tmod 
            Hmww3=WW3_ValuesList[59]
            Tmww3=WW3_ValuesList[60]

        # compute fetch
        if FetchQuestion=='(1) Yes':
            # create fetch vectors
            gp.AddMessage("...creating fetch vectors")

            # radials function
            def Radials(PtsName,Suffix):
                fc=string.replace(PtsName,"\\","/")
                descfc=gp.describe(fc)
                sr=descfc.spatialreference
                # process the feature class attributes
                lstfc=string.split(fc,"/")
                for fl in lstfc:
                    fn=fl
                strWorkspace=string.replace(fc,fl,"")
                gp.workspace=strWorkspace
                newfn=string.replace(fl,".shp",Suffix)
                # check for existence
                if gp.exists(strWorkspace+newfn):
                    gp.delete_management(strWorkspace+newfn)
                    gp.refreshcatalog(gp.workspace)
                # create the fc
                gp.CreateFeatureClass_management(gp.workspace,newfn,"POLYLINE",fc,"SAME_AS_TEMPLATE","SAME_AS_TEMPLATE",sr)
                addrecs=gp.insertcursor(strWorkspace+newfn)
                gp.refreshcatalog(gp.workspace)
                  
                recs=gp.SearchCursor(fc)
                rec=recs.next()
                lstFields=gp.listfields(fc)
                while rec:
                    # get the angle
                    rotation=rec.getvalue("BEARING")
                    length=rec.getvalue("DISTANCE")
                    bearing=math.radians(rotation)
                    angle=math.radians((360-math.degrees(bearing))+90)        
                    # get the feature and compute the to point
                    pt=rec.shape.getpart(0)
                    x=operator.add(math.cos(angle)*length,pt.x)
                    y=operator.add(math.sin(angle)*length,pt.y)
                    # build up the record
                    addrec=addrecs.newrow()    
                    # create the shape
                    newArray=gp.createobject("array")
                    newArray.add (pt)
                    newPt=gp.createobject("point")
                    newPt.x=x
                    newPt.y=y
                    newArray.add(newPt)
                    # maintain the attributes
                    lstFields.reset()
                    fld=lstFields.next()
                    while fld:
                        if fld.name <> "FID" and fld.name <> "OBJECTID" and fld.name <> "SHAPE":
                            addrec.SetValue(fld.name,rec.GetValue(fld.name))
                        fld=lstFields.next()
                    # add shape
                    addrec.shape=newArray
                    addrecs.insertrow(addrec)
                    rec=recs.next()
                del rec,recs    
            
            # copy original point twice and add fields to second copy
            gp.CopyFeatures_management(LandPoint,PtsCopy,"","0","0","0")
            gp.CopyFeatures_management(LandPoint,PtsCopy2,"","0","0","0")
            gp.CopyFeatures_management(LandPoint,PtsCopy3,"","0","0","0")
            PtsCopy2=AddField(PtsCopy2,"DISTANCE","SHORT","8","")
            PtsCopy2=AddField(PtsCopy2,"BEARING","DOUBLE","","")
            PtsCopy2=AddField(PtsCopy2,"BISECTANG","DOUBLE","","")

            CopyExpr=PtsCopy
            for i in range(0,(BearingsNum*9)-2):
                CopyExpr=CopyExpr+";"+PtsCopy

            BiSectAngFullList=[]
            BiSectAngList=[1.0,0.999048,0.996195,0.991445,0.984808,0.984808,0.991445,0.996195,0.999048]
            for i in range(0,16):
                for j in range(0,len(BiSectAngList)):
                    BiSectAngFullList.append(BiSectAngList[j])
            gp.Append_management(CopyExpr,PtsCopy2,"NO_TEST","","")
            gp.CalculateField_management(PtsCopy2,"DISTANCE","50000","PYTHON","")

            # translate information from list into perp transect attribute table
            cur=gp.UpdateCursor(PtsCopy2,"","","BEARING; FID; BISECTANG")
            row=cur.Next()
            m=0
            while row:
                FID=float(row.GetValue("FID"))
                Bearing=float((360.000/((BearingsNum*8.0)+BearingsNum))*FID)
                row.SetValue("Bearing",Bearing)
                row.SetValue("BiSectAng",BiSectAngFullList[m])
                m=m+1
                cur.UpdateRow(row)
                row=cur.Next()
            del cur,row

            # create radials
            Radials(PtsCopy2,"_lineRotate.shp")

            # erase parts of line where it overlaps land (works for ArcView)
            gp.Intersect_analysis(PtsCopyLR+" 1;"+SeaPoly+" 2",PtsCopyEL,"ALL","","INPUT")
            gp.MultipartToSinglepart_management(PtsCopyEL,PtsCopyExp)
            # convert to layer to select only lines originating from point source
            gp.MakeFeatureLayer_management(PtsCopyExp,PtsCopyExp_Lyr,"",gp.workspace,"")
            gp.SelectLayerByLocation_management(PtsCopyExp_Lyr,"WITHIN_A_DISTANCE",LandPoint,"20 Meters","NEW_SELECTION")
            gp.CopyFeatures_management(PtsCopyExp_Lyr,Fetch_Vectors,"","0","0","0")
            # add and calculate "LENGTH_M" field
            Fetch_Vectors=AddField(Fetch_Vectors,"LENGTH_M","LONG","6","")
            gp.CalculateField_management(Fetch_Vectors,"LENGTH_M","!shape.length@meters!","PYTHON","")

            # populate fetch distances to a list
            AngleList=[0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5]
            FetchList=[0.0]*16
            # translate information from list into perp transect attribute table
            cur=gp.UpdateCursor(Fetch_Vectors,"","","BEARING; LENGTH_M")
            row=cur.Next()
            while row:
                Angle=float(row.GetValue("BEARING"))
                if Angle in AngleList:
                    indexAngle=AngleList.index(Angle)
                    FetchList[indexAngle]=float(row.GetValue("LENGTH_M"))
                row=cur.Next()
            del cur,row

            binD1=[]; binBiAng1=[]; binD2=[]; binBiAng2=[]; binD3=[]; binBiAng3=[]; binD4=[]; binBiAng4=[]
            binD5=[]; binBiAng5=[]; binD6=[]; binBiAng6=[]; binD7=[]; binBiAng7=[]; binD8=[]; binBiAng8=[]
            binD9=[]; binBiAng9=[]; binD10=[]; binBiAng10=[]; binD11=[]; binBiAng11=[]; binD12=[]; binBiAng12=[]
            binD13=[]; binBiAng13=[]; binD14=[]; binBiAng14=[]; binD15=[]; binBiAng15=[]; binD16=[]; binBiAng16=[]

            cur=gp.UpdateCursor(Fetch_Vectors,"","","BEARING; LENGTH_M; BISECTANG")
            row=cur.Next()    
            while row:
                Bearing=float(row.GetValue("BEARING"))
                if Bearing >= 350.0 or Bearing <= 10.0:
                    binD1.append(row.GetValue("LENGTH_M"))
                    binBiAng1.append(row.GetValue("BiSectAng"))
                elif Bearing >= 12.5 and Bearing <= 32.5:
                    binD2.append(row.GetValue("LENGTH_M"))
                    binBiAng2.append(row.GetValue("BiSectAng"))
                elif Bearing >= 35.0 and Bearing <= 55.0:
                    binD3.append(row.GetValue("LENGTH_M"))
                    binBiAng3.append(row.GetValue("BiSectAng"))
                elif Bearing >= 57.5 and Bearing <= 77.5:
                    binD4.append(row.GetValue("LENGTH_M"))
                    binBiAng4.append(row.GetValue("BiSectAng"))
                elif Bearing >= 80.0 and Bearing <= 100.0:
                    binD5.append(row.GetValue("LENGTH_M"))
                    binBiAng5.append(row.GetValue("BiSectAng"))
                elif Bearing >= 102.5 and Bearing <= 122.5:
                    binD6.append(row.GetValue("LENGTH_M"))
                    binBiAng6.append(row.GetValue("BiSectAng"))
                elif Bearing >= 125.0 and Bearing <= 145.0:
                    binD7.append(row.GetValue("LENGTH_M"))
                    binBiAng7.append(row.GetValue("BiSectAng"))
                elif Bearing >= 147.5 and Bearing <= 167.5:
                    binD8.append(row.GetValue("LENGTH_M"))
                    binBiAng8.append(row.GetValue("BiSectAng"))
                elif Bearing >= 170.0 and Bearing <= 190.0:
                    binD9.append(row.GetValue("LENGTH_M"))
                    binBiAng9.append(row.GetValue("BiSectAng"))
                elif Bearing >= 192.5 and Bearing <= 212.5:
                    binD10.append(row.GetValue("LENGTH_M"))
                    binBiAng10.append(row.GetValue("BiSectAng"))
                elif Bearing >= 215.0 and Bearing <= 235.0:
                    binD11.append(row.GetValue("LENGTH_M"))
                    binBiAng11.append(row.GetValue("BiSectAng"))
                elif Bearing >= 237.5 and Bearing <= 257.5:
                    binD12.append(row.GetValue("LENGTH_M"))
                    binBiAng12.append(row.GetValue("BiSectAng"))
                elif Bearing >= 260.0 and Bearing <= 280.0:
                    binD13.append(row.GetValue("LENGTH_M"))
                    binBiAng13.append(row.GetValue("BiSectAng"))
                elif Bearing >= 282.5 and Bearing <= 302.5:
                    binD14.append(row.GetValue("LENGTH_M"))
                    binBiAng14.append(row.GetValue("BiSectAng"))
                elif Bearing >= 305.0 and Bearing <= 325.0:
                    binD15.append(row.GetValue("LENGTH_M"))
                    binBiAng15.append(row.GetValue("BiSectAng"))
                elif Bearing >= 327.5 and Bearing <= 347.5:
                    binD16.append(row.GetValue("LENGTH_M"))
                    binBiAng16.append(row.GetValue("BiSectAng"))
                cur.UpdateRow(row)
                row=cur.Next()
            del row,cur
            
            # use 'FetchMean' function to summarize bins
            def FetchCalc(binD,binBiAng,index):
                if len(binD) > 0 and binD.count(0) < 3: # and binD has less than 3 zeros in it.
                    numer=0.0
                    denom=0.0
                    for i in range(0,len(binD)):
                        numer=numer+binD[i]*num.cos(binBiAng[i])
                        denom=denom+num.cos(binBiAng[i])
                    FetchList[index]=(numer/denom)
                return FetchList

            FetchList=num.zeros(16,dtype=num.float64)
            FetchCalc(binD4,binBiAng4,0); FetchCalc(binD3,binBiAng3,1); FetchCalc(binD2,binBiAng2,2); FetchCalc(binD1,binBiAng1,3)
            FetchCalc(binD16,binBiAng16,4); FetchCalc(binD15,binBiAng15,5); FetchCalc(binD14,binBiAng14,6); FetchCalc(binD13,binBiAng13,7)
            FetchCalc(binD12,binBiAng12,8); FetchCalc(binD11,binBiAng11,9); FetchCalc(binD10,binBiAng10,10); FetchCalc(binD9,binBiAng9,11)
            FetchCalc(binD8,binBiAng8,12); FetchCalc(binD7,binBiAng7,13); FetchCalc(binD6,binBiAng6,14); FetchCalc(binD5,binBiAng5,15)

            # save fetch distance values to print out in HTML
            FetchFinalList=[]
            for i in range(3,-1,-1):
                FetchFinalList.append(FetchList[i])
            for i in range(len(dirList)-1,3,-1):
                FetchFinalList.append(FetchList[i])           

            # copy original point add fields to third copy
            gp.CopyFeatures_management(LandPoint,PtsCopy3,"","0","0","0")
            PtsCopy3=AddField(PtsCopy3,"DISTANCE","SHORT","8","")
            PtsCopy3=AddField(PtsCopy3,"BEARING","DOUBLE","","")

            CopyExpr=PtsCopy
            for i in range(0,14):
                CopyExpr=CopyExpr+";"+PtsCopy
            gp.Append_management(CopyExpr,PtsCopy3,"NO_TEST","","")

            # translate information from list into perp transect attribute table
            cur=gp.UpdateCursor(PtsCopy3,"","","BEARING; DISTANCE")
            row=cur.Next()
            m=0
            while row:
                row.SetValue("BEARING",AngleList[m])
                row.SetValue("DISTANCE",FetchFinalList[m])
                m=m+1
                cur.UpdateRow(row)
                row=cur.Next()
            del cur,row

            # create radials
            Radials(PtsCopy3,"_lineDist.shp")
            
            # select fetch lines where "DISTANCE > 0"
            gp.Select_analysis(PtsCopyLD,Fetch_Distances,"\"DISTANCE\" > 0")
            
            # save fetch distances as text
            file=open(FetchDistances,"w")
            for i in range(0,16):
                file.writelines(str(AngleList[i])+"\t"+str(FetchFinalList[i])+"\n")
            file.close()
            

        # compute wave height from wind
        if WW3_Pts and FetchQuestion=='(1) Yes':
            gp.AddMessage("...computing locally generated wind-wave characteristics from fetch")
            g=9.81;
            Fd=FetchFinalList;
            WiWavMax=16*[0.0];WiPerMax=16*[0.0]
            for ff in range(16): # compute wave height in each fetch direction
                if WiMax[ff] <> 0:
                    WiWavMax[ff],WiPerMax[ff]=WindWave(WiMax[ff],Fd[ff],100)
            WiWav10=16*[0.0];WiPer10=16*[0.0]
            for ff in range(16): # compute wave height in each fetch direction
                if Wi10[ff] <> 0:
                    WiWav10[ff],WiPer10[ff]=WindWave(Wi10[ff],Fd[ff],100)
            WiWav25=16*[0.0];WiPer25=16*[0.0]
            for ff in range(16): # compute wave height in each fetch direction
                if Wi25[ff] <> 0:
                    WiWav25[ff],WiPer25[ff]=WindWave(Wi25[ff],Fd[ff],100)
                    
    except:
        gp.AddError(msgPrepWW3Fetch)
        raise Exception


    ##################################          
    ###### READ PG EXCEL INPUT #######
    ##################################
                
    try:
        # read Excel file inputs
        gp.AddMessage("\nReading Nearshore Waves Excel file inputs...")
        xlApp=Dispatch("Excel.Application")
        xlApp.Visible=0
        xlApp.DisplayAlerts=0
        xlApp.Workbooks.Open(InputTable)
        cell=xlApp.Worksheets("ModelInput")
        # sediment
        Diam=cell.Range("e15").Value # sediment diameter (mm)
        A=cell.Range("e17").Value # sediment scale factor
        # tide         
        MSL=cell.Range("f5").Value # mean sea level
        HT=cell.Range("g5").Value # high tide elevation

        # backshore type
        BackHelp=cell.Range("i10").Value # 1) beach, 2) mangroves/marshes,
        if BackHelp==1: # read ProfileGenerator information
            # foreshore    
            Slope=cell.Range("j18").Value # foreshore slope=1/Slope
            BermCrest=cell.Range("j17").Value
            BermLength=cell.Range("j16").Value
            DuneCrest=cell.Range("j15").Value # 1) yes,2) no,3) don't know
            
            if Slope==0 :  # user didn't enter enough data
                xlApp.ActiveWorkbook.Close(SaveChanges=0)
                xlApp.Quit()
                gp.AddError('You did not enter beach information')
                raise Exception
            if BermCrest==0  and DuneCrest==0:  # user didn't enter enough data
                xlApp.ActiveWorkbook.Close(SaveChanges=0)
                xlApp.Quit()
                gp.AddError('You did not enter beach information')
                raise Exception
            m=1.0/Slope # bed slope
            
            # estimate dune size from Short and Hesp
            if DuneCrest==-1:  #User doesn't know
                # wave climate data
                if Tm+Hm==0 and WW3_Pts:    Hm=Hmww3;Tm=Tmww3;   #User doesn't know dune height
                     
                if Tm > 0: # estimate dune height for user
                    Hb=0.39*9.81**(1.0/5)*(Tm*Hm**2)**(2.0/5)
                    a=0.00000126
                    b=num.sqrt(3.61**2+1.18*(1.56*9.81*(Diam/1000.0)**3/a**2)**(1.0/1.53))-3.61
                    ws=(a*b**1.53)/(0.0001*Diam) # fall velocity
                    RTR=HT/Hb
                    if RTR > 3: # in this case, beach is not wave dominated, can't know value, so take zero
                        DuneCrest=0;
                    else: # else, beach is wave dominated (we read Short and Hesp)
                        Type=Hb/(ws*Tm)
                        if Type < 3:
                            DuneCrest=5
                        elif Type < 5:
                            DuneCrest=12
                        else:
                            DuneCrest=20 
                else:
                    Tm=-1;RTR=0
                    DuneCrest=0; # beach has no dune and infinitely long berm

        # profile modification information        
        SlopeM=num.array([0,0,0])
        OffX=num.array([0,0,0])
        ShoreX=num.array([0,0,0])
        ModifInfo=cell.Range("e33:g35").Value
        # Modif1
        temp=ModifInfo[0]
        SlopeM[0]=temp[0]
        OffX[0]=temp[1]
        ShoreX[0]=temp[2]
        # Modif2
        temp=ModifInfo[1]
        SlopeM[1]=temp[0]
        OffX[1]=temp[1]
        ShoreX[1]=temp[2]
        # Modif3
        temp=ModifInfo[2]
        SlopeM[2]=temp[0]
        OffX[2]=temp[1]
        ShoreX[2]=temp[2]
        
        # read in habitat IDs
        if HabDirectory: 
            ExcelHabIDList=[]
            HabAbbrevList=[]
            temp=["e","f","g","h","i","j"]
            for i in range(len(temp)):
                if cell.Range(temp[i]+"39").Value not in [None,'']:
                    ExcelHabIDList.append(int(cell.Range(temp[i]+"39").Value)) # input must be an integer
                    ExcelHabName = str(cell.Range(temp[i]+"40").Value).replace(" ", "")
                    HabAbbrevList.append(ExcelHabName[:9])

            Hab1Zip=zip(ExcelHabIDList,HabAbbrevList)
            Hab1Zip.sort()
            ExcelHabIDList,HabAbbrevList=zip(*Hab1Zip)
            HabAbbrevList = list(HabAbbrevList) # convert from tuple to list        
        
        # save changes and close Excel file
        xlApp.ActiveWorkbook.Close(SaveChanges=0) # don't save changes
        xlApp.Quit()
        
    except:
        xlApp.ActiveWorkbook.Close(SaveChanges=0)
        xlApp.Quit()
        gp.AddError(msgReadExcel)
        raise Exception

    try:
        # cut, read, or create nearshore bathy profile
        if ProfileQuestion=="(1) Yes": # model extracts value from GIS layers
            gp.AddMessage("\nCreating profile points from transect...")
            # create transect and read transect file
            gp.Buffer_analysis(LandPoint,LandPoint_Buff,str(BufferDist)+" Meters","FULL","ROUND","NONE","")
            gp.Extent=LandPoint_Buff
            gp.Intersect_analysis(LandPoint_Buff+"; "+LandPoly,Shoreline_Buff_Clip,"NO_FID","","POINT")
            gp.Extent=""
            # check to make sure that clipped shoreline is not empty FC
            if gp.GetCount_management(Shoreline_Buff_Clip)==0:
                gp.AddError("Shoreline was not found within "+str(BufferDist)+" meters of 'LandPoint' input.  \
                             Either increase the buffer distance or move the 'LandPoint' input closer to the coastline.")
                raise Exception
            
            # set coordinate system to same projection (in meters) as the shoreline point input
            gp.outputCoordinateSystem=LandPoint
            cur=gp.UpdateCursor(LandPoint)
            row=cur.Next()
            feat=row.Shape
            midpoint=feat.Centroid
            midList = midpoint.split(' ')
            midList=[float(s) for s in midList]
            midx=midList[0]
            midy=midList[1]
            del cur,row

            # grab coordinates of the start and end of the coastline segment
            cur=gp.SearchCursor(Shoreline_Buff_Clip)
            row=cur.Next()
            counter=1
            feat=row.Shape
            firstpoint=feat.FirstPoint
            lastpoint=feat.LastPoint
            startList = firstpoint.split(' ')
            endList = lastpoint.split(' ')
            startx=float(startList[0])
            starty=float(startList[1])
            endx=float(endList[0])
            endy=float(endList[1])

            # diagnose the type of perpendicular transect to create (PerpTransType)
            PerpTransType=0
            if starty==endy or startx==endx:
                if starty==endy:
                    y1=midy+TransectDist
                    y2=midy-TransectDist
                    x1=midx
                    x2=midx
                    PerpTransType=1
                if startx==endx:
                    y1=midy
                    y2=midy 
                    x1=midx+TransectDist
                    x2=midx-TransectDist
                    PerpTransType=2
            else:
                # get the slope of the line
                m=((starty-endy)/(startx-endx))
                # get the negative reciprocal
                NegRecip=-1*((startx-endx)/(starty-endy))

                if m > 0:
                    # increase x-values,find y
                    if m >= 1:
                        y1=NegRecip*(TransectDist)+midy
                        y2=NegRecip*(-TransectDist)+midy
                        x1=midx+TransectDist
                        x2=midx-TransectDist
                        PerpTransType=3
                    # increase y-values,find x
                    if m < 1:
                        y1=midy+TransectDist
                        y2=midy-TransectDist
                        x1=(TransectDist/NegRecip)+midx
                        x2=(-TransectDist/NegRecip)+midx
                        PerpTransType=4
                if m < 0:
                    # add to x,find y-values
                    if m >= -1:
                    # add to y,find x-values
                        y1=midy+TransectDist
                        y2=midy-TransectDist
                        x1=(TransectDist/NegRecip)+midx
                        x2=(-TransectDist/NegRecip)+midx
                        PerpTransType=5
                    if m < -1:
                        y1=NegRecip*(TransectDist)+midy
                        y2=NegRecip*(-TransectDist)+midy
                        x1=midx+TransectDist
                        x2=midx-TransectDist
                        PerpTransType=6
            del cur, row

            # grab projection spatial reference from 'LandPoint'
            dataDesc=gp.describe(LandPoint)
            spatialRef=dataDesc.SpatialReference
            gp.CreateFeatureClass_management(interws,"PT1.shp","POINT","#","#","#",spatialRef)
            gp.CreateFeatureClass_management(interws,"PT2.shp","POINT","#","#","#",spatialRef)

            x1,y1,x2,y2=PTCreate(PerpTransType,midx,midy,1)
            xDelta = midx-x1
            yDelta = midy-y1

            # create two point transects, each point is 1 meter away from the previous    
            cur1=gp.InsertCursor(PT1)
            cur2=gp.InsertCursor(PT2)

            for bb in range(1,((RadLineDist/2.0)+1)):
                row1=cur1.NewRow()
                pnt=gp.CreateObject("POINT")
                pnt.x=midx+(bb*(xDelta))
                pnt.y=midy+(bb*(yDelta))
                row1.shape=pnt
                cur1.InsertRow(row1)
                row2=cur2.NewRow()
                pnt=gp.CreateObject("POINT")
                pnt.x=midx-(bb*(xDelta))
                pnt.y=midy-(bb*(yDelta))
                row2.shape=pnt
                cur2.InsertRow(row2)
            del cur1,row1
            del cur2,row2            

            # extract depth values from 'BathyGrid' to point transects
            gp.ExtractValuesToPoints_sa(PT1,BathyGrid,PT1_Z,"INTERPOLATE")
            gp.ExtractValuesToPoints_sa(PT2,BathyGrid,PT2_Z,"INTERPOLATE")
            PT1_Z=AddField(PT1_Z,"PT_ID","LONG","","")
            gp.CalculateField_management(PT1_Z,"PT_ID","!FID! + 1","PYTHON")
            PT2_Z=AddField(PT2_Z,"PT_ID","LONG","","")
            gp.CalculateField_management(PT2_Z,"PT_ID","!FID! + 1","PYTHON")      

            # create depth lists of two point transects
            Dmeas1=[]
            cur=gp.UpdateCursor(PT1_Z)
            row=cur.Next()
            while row:
                Dmeas1.append(row.GetValue("RASTERVALU"))
                cur.UpdateRow(row)
                row=cur.next()
            del cur,row
            
            Dmeas2=[]
            cur=gp.UpdateCursor(PT2_Z)
            row=cur.Next()
            while row:  
                Dmeas2.append(row.GetValue("RASTERVALU"))
                cur.UpdateRow(row)
                row=cur.next()
            del cur,row

            if len(Dmeas1)==0 and len(Dmeas2)==0:
                gp.AddError("\nNeither transect overlaps the seas.  Please check the location of your 'LandPoint' and bathymetry inputs.")
                raise Exception

            # find which point transect hits water first
            def locate(l):
                negative=False
                negativeID=-555555
                positiveID=555555
                for i in range(len(l)):
                    if (negative==False) and (l[i] < 0) and (l[i] != -9999):
                        negative=True
                        negativeID=i
                    if (negative==True) and (l[i] > 0) or (l[i] == -9999):
                        positiveID=i
                        break
                if negativeID==-555555:
                    negativeID=len(l)                     
                if positiveID==555555:
                    positiveID=len(l)
                return negativeID,positiveID

            neg1ID,pos1ID=locate(Dmeas1)
            neg2ID,pos2ID=locate(Dmeas2)

            # create final lists of cross-shore distance (Dx) and depth (Dmeas) 
            Dx=[]   
            Dmeas=[]
            counter=0
            TransectChoice=1
       
            if neg1ID < neg2ID:
                # need to add this shift b/c profile used can be offset from cut profile because of the location of the point shapefile on the shoreline
                # this shift value will be used to correct location of the vegetation so it matches profile used x-axis
                shift=neg1ID
                Xorig=[-len(Dmeas2)+ii for ii in range(len(Dmeas2)+len(Dmeas1))]
                Dorig=[Dmeas2[ii] for ii in range(len(Dmeas2))];Dorig.extend(Dmeas1)
                for i in range(neg1ID,pos1ID):
                    Dx.append(counter)
                    Dmeas.append(Dmeas1[i])
                    counter += 1

            elif neg1ID > neg2ID:
                shift=neg2ID
                Xorig=[-len(Dmeas1)+ii for ii in range(len(Dmeas1)+len(Dmeas2))]
                Dorig=[Dmeas1[ii] for ii in range(len(Dmeas1))];Dorig.extend(Dmeas2)
                for i in range(neg2ID,pos2ID):
                    Dx.append(counter)
                    Dmeas.append(Dmeas2[i])
                    counter += 1
                TransectChoice=2

            else:
                if pos1ID > pos2ID:
                    shift=neg1ID
                    Xorig=[-len(Dmeas2)+ii for ii in range(len(Dmeas2)+len(Dmeas1))]
                    Dorig=[Dmeas2[ii] for ii in range(len(Dmeas2))];Dorig.extend(Dmeas1)
                    for i in range(neg1ID,pos1ID):
                        Dx.append(counter)
                        Dmeas.append(Dmeas1[i])
                        counter += 1
                else:
                    shift=neg2ID
                    Xorig=[-len(Dmeas1)+ii for ii in range(len(Dmeas1)+len(Dmeas2))]
                    Dorig=[Dmeas1[ii] for ii in range(len(Dmeas1))];Dorig.extend(Dmeas2)
                    for i in range(neg2ID,pos2ID):
                        Dx.append(counter)
                        Dmeas.append(Dmeas2[i])
                        counter += 1
                    TransectChoice=2

            # save original profile
            Dorig=num.array(Dorig)
            for kk in range(len(Dorig)):
                if abs(Dorig[kk]) > 100:
                    Dorig[kk]=0
            file=open(ProfileCutGIS,"w")
            for i in range(0,len(Dorig)):
                file.writelines(str(Xorig[i])+"\t"+str(Dorig[i])+"\n")
            file.close()
            
            # create txt for bathy portion
            file=open(BathyProfile,"w")
            for i in range(0,len(Dmeas)):
                file.writelines(str(Dx[i])+"\t"+str(Dmeas[i])+"\n")
            file.close()

            # create final point transect file
            if TransectChoice == 1:
                gp.Select_analysis(PT1_Z,Profile_Pts,"\"PT_ID\" > "+str(neg1ID-1)+" AND \"PT_ID\" < "+str(pos1ID))
                if HabDirectory:
                    # add near and backshore to profile for habitat extraction
                    gp.Select_analysis(PT1_Z,PT_Z_Near,"\"PT_ID\" <= "+str(neg1ID-1))
                    PT_Z_Near=AddField(PT_Z_Near,"DEPTH","DOUBLE","","")
                    gp.CalculateField_management(PT_Z_Near,"DEPTH","!RASTERVALU!","PYTHON")
                    gp.DeleteField_management(PT_Z_Near,"RASTERVALU")
                    PT2_Z=AddField(PT2_Z,"DEPTH","DOUBLE","","")
                    gp.CalculateField_management(PT2_Z,"DEPTH","!RASTERVALU!","PYTHON")
                    gp.DeleteField_management(PT2_Z,"RASTERVALU")
                    gp.CalculateField_management(PT2_Z,"PT_ID","!PT_ID! * -1","PYTHON")
                    gp.Select_analysis(PT2_Z,Backshore_Pts,"\"PT_ID\" < 0 AND \"PT_ID\" > -2001")
                    
            else:
                gp.Select_analysis(PT2_Z,Profile_Pts,"\"PT_ID\" > "+str(neg2ID-1)+" AND \"PT_ID\" < "+str(pos2ID))
                if HabDirectory:
                    # add near and backshore to profile for habitat extraction
                    gp.Select_analysis(PT2_Z,PT_Z_Near,"\"PT_ID\" <= "+str(neg2ID-1))
                    PT_Z_Near=AddField(PT_Z_Near,"DEPTH","DOUBLE","","")
                    gp.CalculateField_management(PT_Z_Near,"DEPTH","!RASTERVALU!","PYTHON")
                    gp.DeleteField_management(PT_Z_Near,"RASTERVALU")
                    PT1_Z=AddField(PT1_Z,"DEPTH","DOUBLE","","")
                    gp.CalculateField_management(PT1_Z,"DEPTH","!RASTERVALU!","PYTHON")
                    gp.DeleteField_management(PT1_Z,"RASTERVALU")
                    gp.CalculateField_management(PT1_Z,"PT_ID","!PT_ID! * -1","PYTHON")
                    gp.Select_analysis(PT1_Z,Backshore_Pts,"\"PT_ID\" < 0 AND \"PT_ID\" > -2001")
                    
            # add and calculate field for "DEPTH"
            Profile_Pts=AddField(Profile_Pts,"DEPTH","DOUBLE","","")
            gp.CalculateField_management(Profile_Pts,"DEPTH","!RASTERVALU!","PYTHON")
            gp.DeleteField_management(Profile_Pts,"RASTERVALU")

            # merge 'Profile_Pts' with 'Land Point' and backshore profile up to 2km
            if HabDirectory:
                MergeExpr=Profile_Pts+"; "+PT_Z_Near+"; "+LandPoint+"; "+Backshore_Pts
                gp.Merge_management(MergeExpr,Profile_Pts_Merge,"")

            # smooth profile and create x axis
            lx=len(Dmeas) # length of original data
            Dx=num.array(Dx)
            Dmeas=num.array(Dmeas) # shoreline starts at x=0
            
            yd=[Dmeas[ii] for ii in range(len(Dx))]
            xd=[Dx[ii] for ii in range(len(Dx))]
            SmoothValue=(SmoothParameter/100.0)*len(Dmeas)
            TempY=SignalSmooth.smooth(Dmeas,round(SmoothValue,2),'flat') # smooth function
            TempX=num.array(Dx)

            # remove portions offshore that are shallower than average depth above deepest point
            LocDeep=Indexed(TempY,min(TempY)) # locate deepest point
            davg=mean(TempY);davg=average([davg,min(TempY)])
            out=num.nonzero(TempY<davg);out=out[0] # locate points deeper than average depth

            loc=out[-1]
            if loc>LocDeep: # if point is offshore of deepest point
                out=num.arange(loc,len(yd),1)
                yd=num.delete(yd,out,None)
                xd=num.delete(Dx,out,None);xd=xd-xd[0]
                TempY[out]=0
                TempX[out]=-1         
            lx=len(xd) 
            yd=num.array(yd);xd=num.array(xd)
            
            # insert 'TempY' and 'TempX' in shapefile output
            Profile_Pts=AddField(Profile_Pts,"SM_BATHY_X","DOUBLE","","")
            Profile_Pts=AddField(Profile_Pts,"SM_BATHY_Y","DOUBLE","","")
            cur=gp.UpdateCursor(Profile_Pts)
            row=cur.Next()
            j=0
            while row:
                row.SetValue("SM_BATHY_X",TempX[j])
                row.SetValue("SM_BATHY_Y",TempY[j])
                j+=1
                cur.UpdateRow(row)
                row=cur.Next()
            del cur,row,TempX,TempY

            # read habitat GIS layers to profile
            if HabDirectory:
                gp.AddMessage("...locating each habitat's presence along the profile")
                gp.workspace=HabDirectory
                fcList=gp.ListFeatureClasses("*","all")
                fc=fcList.Next()
                HabLyrList=[]
                HabIDList=[]
                fcCount = 0
                while fc:
                    checkGeometry(fc,"Polygon","Natural Habitat")
                    checkInteger(fc)
                    ckProjection(fc)
                    HabLyrList.append(fc)
                    fc2=fc[::-1]
                    j=fc2.find('_')
                    indexS=len(fc)-j-1
                    indexE=fc.find(".")
                    fc_ID=fc[indexS+1:indexE]
                    HabIDList.append(int(fc_ID))
                    fc=fcList.Next()
                    fcCount =+ 1
                del fc

                # check that habitat data directory has at least one layer
                if fcCount == 0:
                    gp.AddError("\nNo habitat layers were found in the specified 'Habitat Data Directory'.")
                    raise Exception

                Hab2Zip=zip(HabIDList,HabLyrList)
                Hab2Zip.sort()
                HabIDList,HabLyrList=zip(*Hab2Zip)

                # check that 'HabDirectory' and Excel input are consistent
                if ExcelHabIDList <> HabIDList:
                    gp.AddError("There is an inconsistency between the number of habitat layers in the specified directory and the input spreadsheet.")
                    raise Exception

                gp.workspace=gp.GetParameterAsText(0) # reset workspace
                gp.Extent=Profile_Pts_Merge # revise extent
                
                # rasterize the layers and generate zone of influence
                AbbrevList=['Mangroves', 'Marsh', 'SeagrassB', 'SandDunes', 'CoralReef', 'Other']
                ExcludeList=["FID","Shape","Id","PT_ID","DEPTH"]
                HabPresentCount = 0
                IntersectExpr=''
                for i in range(0,len(HabLyrList)):
                    HabVector=HabDirectory+"\\"+HabLyrList[i]
                    HabVector=AddField(HabVector,"VID","SHORT","","")
                    gp.CalculateField_management(HabVector,"VID",1,"PYTHON")
                    gp.FeatureToRaster_conversion(HabVector,"VID",interws+HabAbbrevList[i],"10")
                    gp.BuildRasterAttributeTable_management(interws+HabAbbrevList[i],"OVERWRITE")
                    if gp.GetCount(interws+HabAbbrevList[i]) > 0:
                        pass
                        HabPresentCount += 1
                        gp.Reclassify_sa(interws+HabAbbrevList[i],"VALUE","1 1;NODATA 0",interws+HabAbbrevList[i]+"_rc","DATA")
                        gp.ExtractValuesToPoints_sa(Profile_Pts_Merge,interws+HabAbbrevList[i]+"_rc",interws+HabAbbrevList[i]+".shp","NONE")
                        gp.AddField_management(interws+HabAbbrevList[i]+".shp",HabAbbrevList[i],"SHORT","0","0","","","NON_NULLABLE","NON_REQUIRED","")
                        gp.CalculateField_management(interws+HabAbbrevList[i]+".shp",HabAbbrevList[i],"[RASTERVALU]","VB")
                        gp.DeleteField_management(interws+HabAbbrevList[i]+".shp","RASTERVALU")
                        if i==0:
                            IntersectExpr=IntersectExpr+interws+HabAbbrevList[i]+".shp "+str(i+1)
                        else:
                            IntersectExpr=IntersectExpr+"; "+interws+HabAbbrevList[i]+".shp "+str(i+1)

                if HabPresentCount > 0: # process if at least one habitat layer is within the AOI    
                    # intersect the various habitat profile point plots
                    if HabPresentCount == 1: # if only 1, just make a copy
                        gp.CopyFeatures_management(IntersectExpr[:-2],Profile_Pts_Hab,"","0","0","0")
                    else: # otherwise, intersect them
                        gp.Intersect_analysis(IntersectExpr,Profile_Pts_Hab,"NO_FID","","INPUT")

                    # delete extra fields
                    DeleteFieldList=[]
                    HabFieldList=[]
                    fieldList=gp.ListFields(Profile_Pts_Hab,"*","All")
                    field=fieldList.Next()
                    while field <> None:
                        if field.Name not in HabAbbrevList+ExcludeList:
                            DeleteFieldList.append(field.Name)
                        if field.Name in AbbrevList:
                            HabFieldList.append(field.Name)
                        field=fieldList.Next()
                    del fieldList,field

                    DelExpr=''
                    for i in range(0,len(DeleteFieldList)):
                        DelExpr=DelExpr+DeleteFieldList[i]+";"
                    DelExpr=DelExpr[:-1]
                    if len(DelExpr) > 0:                    
                        gp.DeleteField_management(Profile_Pts_Hab,DelExpr)

                    # add all habitat abbreviation fields, even if not part of inputs
                    for k in range(0,len(AbbrevList)):
                        if AbbrevList[k] not in HabFieldList:
                            Profile_Pts_Hab=AddField(Profile_Pts_Hab,AbbrevList[k],"SHORT","0","0")

                    # read 'Profile_Pts_Hab' to array 'ProfileHabArray'
                    ProfileHabLength = gp.GetCount_management(Profile_Pts_Hab)
                    ProfileHabList = np.zeros(ProfileHabLength*8, dtype=np.float64)
                    ProfileHabArray = np.reshape(ProfileHabList, (ProfileHabLength,8))
                    cur=gp.UpdateCursor(Profile_Pts_Hab)
                    row=cur.Next()
                    j=0
                    while row:
                        # order will be: PT_ID, DEPTH, MG, MR, SG, DN, CR, OT
                        ProfileHabArray[j][0]=row.GetValue("PT_ID")
                        ProfileHabArray[j][1]=row.GetValue("DEPTH")
                        for i in range(0,len(AbbrevList)):
                            ProfileHabArray[j][i+2]=row.GetValue(AbbrevList[i])
                        cur.UpdateRow(row)
                        row=cur.next()
                        j+=1
                    del cur,row
                    
                    # sort the array by 'PT_ID'
                    ProfileHabArray = ProfileHabArray[ProfileHabArray[:,0].argsort()]

                    # export bathy and vegetation information       
                    TextData=open(HabitatLocation,"w")
                    for i in range(0,ProfileHabLength):
                        for j in range(0,8):
                            TextData.write(str(ProfileHabArray[i][j])+" ")
                        TextData.write("\n")
                    TextData.close() 
                    
                    # locate the Xbeg and Xend of each patch
                    temp=transpose(ProfileHabArray)
                    TempX=temp[0]
                    TempY=temp[1]
                    MG1=temp[2];MR1=temp[3];SG1=temp[4];DN1=temp[5];CR1=temp[6];OT1=temp[7]
                    MG1[find(MG1==-999)]=0;SG1[find(SG1==-999)]=0;MR1[find(MR1==-999)]=0
                    DN1[find(DN1==-999)]=0;CR1[find(CR1==-999)]=0;OT1[find(OT1==-999)]=0
                    MG2=MG1*0;MR2=MR1*0;SG2=SG1*0;DN2=DN1*0;CR2=CR1*0;OT2=OT1*0;
                    
                    # mangrove, marsh and dune: switch beginning and end indices b/c they're on land
                    finMGx1,begMGx1=LocateHabitat(TempX,MG1);begMGx2=num.array([]);finMGx2=num.array([]);
                    finMGx1=finMGx1[::-1];begMGx1=begMGx1[::-1]; # switch order to go sea to land
                    temp1=find(finMGx1>0);temp2=find(begMGx1>0); # if there are points in the water, make a note
                    if len(temp1)+len(temp2):
                        gp.AddWarning("......you have mangrove in subtidal areas.  Please double check your habitat input layer.");mangwater=1;
                    del temp1,temp2
                    
                    finMRx1,begMRx1=LocateHabitat(TempX,MR1);begMRx2=num.array([]);finMRx2=num.array([]);
                    finMRx1=finMRx1[::-1];begMRx1=begMRx1[::-1]; # switch order to go sea to land
                    temp1=find(finMRx1>0);temp2=find(begMRx1>0); # if there are points in the water, make a note
                    if len(temp1)+len(temp2):
                        gp.AddWarning("......you have marsh in subtidal areas.  Please double check your habitat input layer.");marwater=1;
                    del temp1,temp2

                    begDNx1,finDNx1=LocateHabitat(TempX,DN1);begDNx2=num.array([]);finDNx2=num.array([]);
                    finDNx1=finDNx1[::-1];begDNx1=begDNx1[::-1]; # switch order to go sea to land
                    temp1=find(finDNx1>0);temp2=find(begDNx1>0); # if there are points in the water, make a note
                    if len(temp1)+len(temp2):
                        gp.AddWarning("......you have a dune in subtidal areas.  Please double check your habitat input layer.");dunewater=1;
                    del temp1,temp2
                    
                    # seagrass and corals are in subtidal areas
                    begSGx1,finSGx1=LocateHabitat(TempX,SG1);begSGx2=num.array([]);finSGx2=num.array([]);
                    temp1=find(finSGx1<0);temp2=find(begSGx1<0); # if there are points on land, make a note
                    if len(temp1)+len(temp2):
                        gp.AddWarning("......you have seagrass on land.  Please double check your habitat input layer.");seagwater=1;
                    del temp1,temp2
                    
                    begCRx1,finCRx1=LocateHabitat(TempX,CR1);begCRx2=num.array([]);finCRx2=num.array([]);
                    temp1=find(finCRx1<0);temp2=find(begCRx1<0); # if there are points on land, make a note
                    if len(temp1)+len(temp2):
                        gp.AddWarning("......you have a coral reef on land.  Please double check your habitat input layer.");coralwater=1;
                    del temp1,temp2
                     
                    # can't say much about other
                    begOTx1,finOTx1=LocateHabitat(TempX,OT1);begOTx2=num.array([]);finOTx2=num.array([]);
                    
                    del ProfileHabArray
                else:
                    gp.AddWarning(".....no habitat was found to be within the area of your profile.")


        # upload user's profile
        elif ProfileQuestion=="(2) No, but I will upload a cross-shore profile":
            gp.AddMessage("\nRetrieving your profile...")
            # read in user's cross-shore profile
            TextData=open(CSProfile,"r") 
            Dx=[];Dmeas=[]
            for line in TextData.read().strip("\n").split("\n"):
                linelist=[float(s) for s in line.split("\t")] # split the list by tab delimiter
                Dx.append(linelist[0])
                Dmeas.append(linelist[1])
            TextData.close()

            Dmeas=num.array(Dmeas)
            Dx=num.array(Dx)
            Lx=len(Dx)
            if Lx > 5:
                F=interp1d(Dx,Dmeas);
                xd=num.arange(Dx[0],Dx[-1],1)
                lx=len(xd)
                temp=F(xd)
                yd=num.array(temp)
            else:
                xd=num.array(Dx)
                yd=num.array(Dmeas)
                lx=len(yd)

            if BackHelp==1:
                keep=num.nonzero(yd<=0);keep=keep[0]
                xd=xd[keep];yd=yd[keep]
                if len(keep)<len(yd):
                    gp.AddMessage("\nA portion of your uploaded profile above 'Mean Lower Low Water' was removed so we can build your backshore profile.")
            yd=num.array(yd);xd=num.array(xd)    
             
        # equilibrium beach profile,in case we don't have nearshore bathy
        elif ProfileQuestion=="(3) No, but please create a theoretical profile for me":
            if Diam < 0.1:
                gp.AddError("\nCannot create an equilibrium profile for cohesive/fine sediments (size smaller than 0.1mm).\
                             \nYou may want to adjust your site's sediment size in the 'Profile Generator Excel Table' input.")
                raise Exception
            else:
                gp.AddMessage("\nCreating theoretical profile for you...")

            x=num.arange(0.0,10001.0,1.0) # long axis
            temp=2.0/3
            Dmeas=-A*x**(temp)# eq. profile
            out=num.nonzero(Dmeas<-20);temp1=out[0]
            out=num.nonzero(Dmeas>-1);temp2=out[0]
            out=num.append(temp1,temp2)
            Dmeas=num.delete(Dmeas,out,None) # water depths down to -20
            Dx=num.delete(x,out,None);lx=len(Dx)
            yd=num.array(Dmeas);xd=num.array(Dx)

    except:
        gp.AddError(msgCreateProfile)
        raise Exception

    try:
        # profile modification
        gp.AddMessage("...customizing depth profile")
        yd=yd-MSL # yd is referenced to MLLW; need to reference to MSL
        
        # smooth the signal
        SmoothValue=round((SmoothParameter/100.0)*len(xd),2)
        yd=SignalSmooth.smooth(num.array(yd),SmoothValue,'flat')

        if BackHelp==1: # create backshore profile for beach systems 
            # add foreshore from MLLW to Berm Elevation
            MLLW_loc=Indexed(yd,-MSL) #Locate MLLW on profile
            x,y=SlopeModif(xd,yd,Slope,xd[MLLW_loc],-2000)
            AbvBerm=num.nonzero(y> BermCrest)
            x=num.delete(x,AbvBerm[0],None)
            y=num.delete(y,AbvBerm[0],None) # remove values that are above BermCrest
            
            # Add Berm and Dune
            if BermLength <>0:
                x,y=SlopeModif(x,y,0,x[0],x[0]-BermLength-1)    
            if DuneCrest<>0:
                y[0]=DuneCrest+BermCrest
                x,y=SlopeModif(x,y,0,x[0],x[0]-50)
            
            xd=num.array(x)
            yd=num.array(y)
            
        # slope modification
        if sum(SlopeM)<>0:
            for oo in range(len(SlopeM)):
                if OffX[oo]-ShoreX[oo]<>0:
                    xd,yd=SlopeModif(xd,yd,SlopeM[oo],OffX[oo],ShoreX[oo])
            # smooth signal again
            SmoothValue=round((2./100.)*len(xd),2)
            yd=SignalSmooth.smooth(num.array(yd),SmoothValue,'flat')
        
        if BackHelp==1 and any(SlopeM):
            gp.AddWarning("......you have asked us to create a sand dune AND modify your profile.  Make sure those demands are not contradictory.")
                
    except:
        gp.AddError(msgModifyProfile)
        raise Exception
    
    # locate X=0 at y=0
    DistZero=xd[Indexed(yd,0)]
    xd=xd-DistZero # x=0 is now located at y=0


    ###################################          
    ########## PLOT OUTPUTS ###########
    ###################################

    try:
        gp.AddMessage("...plotting profile and creating outputs\n")
        Fig3=0;Fig2=0;Fig6=0
        # plot and save
        figure(1)
        ax=subplot(111);plot(xd,yd,xd,yd*0,'k',xd,yd*0-MSL,'--k',xd,yd*0+HT-MSL,'-.k',linewidth=2);grid()
        box=ax.get_position();
        ax.set_position([box.x0, box.y0, box.width*1, box.height])
        ax.legend(('Created Profile','Mean Sea Level','Mean Low Water','Mean High Water'),prop=fontP,loc='upper left',ncol=2)
        ylabel('Depth [m]',size='large',weight='bold')
        xlabel('Cross-Shore Distance [m]',size='large',weight='bold')
        axvline(x=0, linewidth=1, color='k')
        xlim(xd[-1],xd[0])
        title('Created Profile',size='large',weight='bold')
        savefig(html_txt+"CreatedProfile.png",dpi=(640/8));

        Xzoom=Indexed(yd,1) # locate point 1m abv MSL
        figure(2)
        if ProfileQuestion=="(1) Yes": # profile was cut with GIS
            # zoom on nearshore area
            ax=subplot(211);plot(Dx-DistZero,Dmeas-MSL,'b',xd,yd,'r',linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            ax.legend(('Initial Profile','Created Profile'), prop=fontP,loc='upper left',ncol=2)
            title('Bathymetry-Smoothing Factor='+str(SmoothParameter),size='large',weight='bold')
            xlim(Dx[-1],-MSL)
            ylabel('Depth [m]',size='large',weight='bold')
            
        elif ProfileQuestion=="(2) No, but I will upload a cross-shore profile":
            # zoom on nearshore area
            ax=subplot(211);plot(Dx,Dmeas-MSL,'b',xd,yd,'r',linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            ax.legend(('Initial Profile','Created Profile'), prop=fontP,loc='upper left',ncol=2)
            title('Bathymetry-Smoothing Factor='+str(SmoothParameter),size='large',weight='bold')
            xlim(max(Dx[-1],xd[-1]),0);ylim(Dmeas[0]*1.1,MSL)
            ylabel('Depth [m]',size='large',weight='bold')
    
        elif ProfileQuestion=="(3) No, but please create a theoretical profile for me":
            # zoom on nearshore area
            ax=subplot(211);plot(xd,yd,'r',linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            ax.legend(('Created Profile'), prop=fontP,loc='upper left',ncol=2)
            title('Bathymetry-Smoothing Factor='+str(SmoothParameter),size='large',weight='bold')
            xlim(Dx[-1],-MSL)
            ylabel('Depth [m]',weight='bold',size='large')
                 
        if (any(SlopeM+OffX+ShoreX) or BackHelp==1) and Xzoom>0: # plot point abv MLLW
            Xzoom=Indexed(yd,-MSL-5) # locate point 5m below MLLW                
            ax=subplot(212);plot(xd[0:Xzoom],yd[0:Xzoom],xd[0:Xzoom],yd[0:Xzoom]*0,'k',xd[0:Xzoom],yd[0:Xzoom]*0-MSL,'--k',xd[0:Xzoom],yd[0:Xzoom]*0+HT-MSL,'-.k',linewidth=2);grid()
            ax.legend(('Created Profile','Mean Sea Level','Mean Low Water','Mean High Water'),prop=fontP,loc='upper left',ncol=2)
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            xlabel('Cross-Shore Distance [m]',weight='bold')
            ylabel('Elevation [m]',size='large',weight='bold')
            xlim(xd[Xzoom],xd[0]);
            if yd[0]<yd[Xzoom]:
                ylim(yd[0]*1.1,yd[Xzoom]*1.1)
            else:
                ylim(yd[Xzoom]*1.1,yd[0]*1.1)
                
            axvline(x=0, linewidth=1, color='k')
            title('Foreshore and Backshore Areas',size='large',weight='bold')
        savefig(html_txt+"ZoomIns.png",dpi=(640/8))
 
        if ProfileQuestion=="(1) Yes": # profile was cut with GIS
            Dorig=num.array(Dorig);
            Xorig=num.array(Xorig)
            fig=figure(3)
            ax=subplot(211);plot(Xorig-shift+1,Dorig,Dx,Dmeas,'r',Xorig,Dorig*0,'k',linewidth=2);grid()
            ax.legend(('Raw Profile','Extracted Bathy'), prop=fontP,loc='upper left',ncol=2)            
            box=ax.get_position();

            ax.set_position([box.x0, box.y0, box.width*1, box.height])
            title('DEM Profile',size='large',weight='bold')
            ylabel('Elevation [m]',size='large',weight='bold')
            xlabel('Cross-Shore Distance [m]',size='large',weight='bold')
            xlim(Xorig[-1],Xorig[0])
            axvline(x=0, linewidth=1, color='k')
            extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()) # save subplot
            savefig(html_txt+"GISProfile.png", bbox_inches=extent.expanded(1.8, 1.4),dpi=(640/8))

        # plot vegetation 
        if ProfileQuestion=="(1) Yes":
            try:
                HabPrst=int(any(MG1))+int(any(MR1))+int(any(CR1))+int(any(DN1))+int(any(SG1))+int(any(OT1));
            except:
                HabPrst=0
                
            if HabPrst<>0:
                TempY=TempY-MSL # zero is at MSL
                
                # remove points that are too deep (mostly lack of data) and resample to dx=1
                out=num.nonzero(TempY<-500);out=out[0] # points that are in ok water depths
                TempY[out]=0
                leg=1; # only one legend
                
                fig=figure(4);HabCount=1
                if any(MG1):
                    ax=subplot(HabPrst,1,HabCount); temp1=TempY+nan;temp2=TempY+nan
                    la=num.nonzero(MG1);temp1[la]=TempY[la]
                    la=num.nonzero(MG2);temp2[la]=TempY[la]
                    plot(TempX,TempY,xd,yd,'r');plot(TempX-shift,temp1,'or',markersize=15);grid()
                    xlim(TempX[-1],TempX[0])     
                    if leg:
                        ax.legend(('Raw Profile (GIS)','Created Profile','Habitat'), prop=fontP,loc='upper left',ncol=2);leg=0
                    ylabel('Mangrove',size='large',weight='bold');HabCount=HabCount+1
                    axvline(x=0, linewidth=1, color='k')
                if any(DN1):
                    ax=subplot(HabPrst,1,HabCount); temp1=TempY+nan;temp2=TempY+nan
                    la=num.nonzero(DN1);temp1[la]=TempY[la]
                    la=num.nonzero(DN2);temp2[la]=TempY[la]
                    plot(TempX,TempY,xd,yd,'r');plot(TempX-shift,temp1,'or',markersize=15);grid()
                    xlim(TempX[-1],TempX[0])                    
                    if leg:
                        ax.legend(('Raw Profile (GIS)','Created Profile','Habitat'), prop=fontP,loc='upper left',ncol=2);leg=0
                    box=ax.get_position();
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                    ylabel('Dunes',size='large',weight='bold');HabCount=HabCount+1
                    axvline(x=0, linewidth=1, color='k')
                if any(MR1):
                    ax=subplot(HabPrst,1,HabCount); temp1=TempY+nan;temp2=TempY+nan
                    la=num.nonzero(MR1);temp1[la]=TempY[la]
                    la=num.nonzero(MR2);temp2[la]=TempY[la]
                    plot(TempX,TempY,xd,yd,'r');plot(TempX-shift,temp1,'or',markersize=15);grid()
                    xlim(TempX[-1],TempX[0])                    
                    if leg:
                        ax.legend(('Raw Profile (GIS)','Created Profile','Habitat'), prop=fontP,loc='upper left',ncol=2);leg=0
                    ylabel('Marsh',size='large',weight='bold');HabCount=HabCount+1
                    box=ax.get_position();
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                    axvline(x=0, linewidth=1, color='k')
                if any(SG1):
                    ax=subplot(HabPrst,1,HabCount); temp1=TempY+nan;temp2=TempY+nan
                    la=num.nonzero(SG1);temp1[la]=TempY[la]
                    la=num.nonzero(SG2);temp2[la]=TempY[la]
                    plot(TempX,TempY,xd,yd,'r');plot(TempX,temp1,'or',markersize=15);grid()
                    xlim(TempX[-1],TempX[0])                    
                    if leg:
                        ax.legend(('Raw Profile (GIS)','Created Profile','Habitat'), prop=fontP,loc='upper left',ncol=2);leg=0
                    ylabel('Seagrass',size='large',weight='bold');HabCount=HabCount+1
                    box=ax.get_position();
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                    axvline(x=0, linewidth=1, color='k')
                if any(CR1):
                    ax=subplot(HabPrst,1,HabCount); temp1=TempY+nan;temp2=TempY+nan
                    la=num.nonzero(CR1);temp1[la]=TempY[la]
                    la=num.nonzero(CR2);temp2[la]=TempY[la]
                    plot(TempX,TempY,xd,yd,'r');plot(TempX-shift,temp1,'or',markersize=15);grid()
                    xlim(TempX[-1],TempX[0])                    
                    if leg:
                        ax.legend(('Raw Profile (GIS)','Created Profile','Habitat'), prop=fontP,loc='upper left',ncol=2);leg=0
                    ylabel('Coral',size='large',weight='bold');HabCount=HabCount+1
                    box=ax.get_position();
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                    axvline(x=0, linewidth=1, color='k')
                if any(OT1):
                    ax=subplot(HabPrst,1,HabCount); temp1=TempY+nan;temp2=TempY+nan
                    la=num.nonzero(OT1);temp1[la]=TempY[la]
                    la=num.nonzero(OT2);temp2[la]=TempY[la]
                    plot(TempX,TempY,xd,yd,'r');plot(TempX-shift,temp1,'or',markersize=15);grid()
                    xlim(TempX[-1],TempX[0])                    
                    if leg:
                        ax.legend(('Raw Profile (GIS)','Created Profile','Habitat'), prop=fontP,loc='upper left',ncol=2);leg=0
                    ylabel('Other',size='large',weight='bold');HabCount=HabCount+1
                    box=ax.get_position();
                    ax.set_position([box.x0, box.y0, box.width*1, box.height])
                xlabel('Cross-Shore Distance [m]',weight='bold')
                axvline(x=0, linewidth=1, color='k')
                savefig(html_txt+"HabFigure.png",dpi=(640/8))

        # fetch and wind rose plots
        if FetchQuestion=='(1) Yes':
            # plot fetch on rose    
            radians=(num.pi/180.0)
            pi=num.pi
            theta16=[0*radians,22.5*radians,45*radians,67.5*radians,90*radians,112.5*radians,135*radians,157.5*radians,180*radians,202.5*radians,225*radians,247.5*radians,270*radians,292.5*radians,315*radians,337.5*radians]
            rc('grid',color='#316931',linewidth=1,linestyle='-')
            rc('xtick',labelsize=0)
            rc('ytick',labelsize=15)
            # force square figure and square axes looks better for polar,IMO
            width,height=matplotlib.rcParams['figure.figsize']
            size=min(width,height)
            # make a square figure
            fig1=plt.figure()
            ax=fig1.add_axes([0.1,0.1,0.8,0.8],polar=True,axisbg='w')
            # plot
            bars=ax.bar(theta16,FetchList,width=.35,color='#ee8d18',lw=1)
            for r,bar in zip(FetchList,bars):
                bar.set_facecolor(cm.YlOrRd(r/10.))
                bar.set_alpha(.65)
            ax.set_rmax(max(FetchList)+1)
            grid(True)
            ax.set_title("Average Fetch (meters)",fontsize=15,weight='bold')
            plt.savefig(Fetch_Plot,dpi=(640/8))

        if WW3_Pts:
            # modify 'Wi10' list so it conforms to rose order
            Wi10List=[]
            for i in range(3,-1,-1):
                Wi10List.append(Wi10[i])
            for i in range(len(dirList)-1,3,-1):
                Wi10List.append(Wi10[i])
            # plot wind on rose    
            radians=(num.pi/180.0)
            pi=num.pi
            theta16=[0*radians,22.5*radians,45*radians,67.5*radians,90*radians,112.5*radians,135*radians,157.5*radians,180*radians,202.5*radians,225*radians,247.5*radians,270*radians,292.5*radians,315*radians,337.5*radians]
            rc('grid',color='#316931',linewidth=1,linestyle='-')
            rc('xtick',labelsize=0)
            rc('ytick',labelsize=15)
            # force square figure and square axes looks better for polar,IMO
            width,height=matplotlib.rcParams['figure.figsize']
            size=min(width,height)
            # make a square figure
            fig1=plt.figure()
            ax=fig1.add_axes([0.1,0.1,0.8,0.8],polar=True,axisbg='w')
            # plot
            bars=ax.bar(theta16,Wi10List,width=.35,color='#ee8d18',lw=1)
            for r,bar in zip(Wi10List,bars):
                bar.set_facecolor(cm.YlOrRd(r/10.))
                bar.set_alpha(.65)
            ax.set_rmax(max(Wi10List)+1)
            grid(True)
            ax.set_title("Wind Speed (meters/second)",fontsize=15,weight='bold')
            plt.savefig(Wind_Plot,dpi=(640/8))

        # create txt profile for created portion
        file=open(CreatedProfile,"w")
        for i in range(0,len(yd)):
            file.writelines(str(xd[i])+"\t"+str(yd[i])+"\n")
        file.close()

        # return projected point to geographic (unprojected)
        gp.Project_management(LandPoint,LandPoint_Geo,geo_projection)
        # grab coordinates for Google Maps plot
        cur=gp.UpdateCursor(LandPoint_Geo)
        row=cur.Next()
        feat=row.Shape
        midpoint1=feat.Centroid
        midList1=midpoint1.split(' ')
        midList1=[float(s) for s in midList1]
        del cur,row
        PtLat=str(round(midList1[1],4))
        PtLong=str(round(midList1[0],4))
        
        TR=str(HT)
        WaveClimateCheck=0

    except:
        gp.AddError(msgPlotProfile)
        raise Exception


    ##################################          
    ####### CREATE HTML FILES ########
    ##################################

    try:
        htmlfile=open(Profile_HTML,"w")
        htmlfile.write("<html><title>Marine InVEST - Profile Generator</title><CENTER><H1>Coastal Protection - Tier 1</H1><H2>Profile Generator Results ("+subwsStr+")<br></H2><p>")
        if ProfileQuestion=="(1) Yes":
            htmlfile.write("[ PROFILE INFO ]<br><a href=\"fetchwindwave.html\">[ FETCH, WAVE, AND WIND INFO ]</a><br>")
        htmlfile.write("<CENTER><br><HR>")
        htmlfile.write("<table border=\"0\" width=\"800\" cellpadding=\"5\" cellspacing=\"20\"><tr><td>")
        htmlfile.write("<iframe width=\"350\" height=\"325\" frameborder=\"0\" scrolling=\"no\" marginheight=\"0\" marginwidth=\"0\"")
        htmlfile.write("src=\"http://maps.google.com/maps?f=q&amp;source=s_q&amp;hl=en&amp;geocode=&amp;q=")
        htmlfile.write(PtLat+","+PtLong)
        htmlfile.write("&amp;aq=&amp;sspn=0.009467,0.021136&amp;vpsrc=6&amp;ie=UTF8&amp;ll=")
        htmlfile.write(PtLat+","+PtLong)
        htmlfile.write("&amp;spn=0.063714,0.169086&amp;t=h&amp;z=12&amp;output=embed\"></iframe>")
        htmlfile.write("</td><td>")
        htmlfile.write("<H2><u>Site Information</u></H2>")
        htmlfile.write("The site is located at:<br><li> latitude = "+PtLat+"<br><li> longitude = "+PtLong+"<p>")
        htmlfile.write("The tidal range is: "+TR+"m (high tide value)<br>")
        if HT < 2:
            htmlfile.write("Your site is microtidal (Tidal Range < 2m)<br>")
        elif HT <= 4:
            htmlfile.write("Your site is meso-tidal (2 <= Tidal Range <= 4m)<br>")
        else:
            htmlfile.write("Your site is macro-tidal (Tidal Range > 4m)<br>")
        htmlfile.write("</td></tr></table>")
        htmlfile.write("Your average backshore sediment size is: "+str(Diam)+"mm<br><p>")
        if Diam > 1.1:
            htmlfile.write("<i>Your beach is composed of coarse sand/gravel. It won't be eroded during a storm</i><p>")
        elif Diam >= 0.1:
            htmlfile.write("<i>You have a sandy beach.  It can be eroded during a storm</i><p>")
        else:
            htmlfile.write("<i>Your backshore is composed of fines/consolidated sediments. It is not an erodible beach</i><p>")
        
        # backshore information
        if BackHelp==1:
            htmlfile.write("<HR><H2><u>Bathymetry and Backshore Profiles for a Sandy Beach System</u></H2>")
        elif BackHelp==2:
            htmlfile.write("<HR><H2><u>Bathymetry and Backshore Profiles for a Mangrove or Marsh System</u></H2>")
            
        if BackHelp==1 and any(SlopeM):
            htmlfile.write("<b>You have asked us to create a sand dune AND modify your profile.  Make sure those demands are not contradictory.</b>");
            
        htmlfile.write("<table border=\"0\" width=\"1100\" cellpadding=\"5\" cellspacing=\"10\">")
        htmlfile.write("<tr><td>The figure below shows the entire smoothed profile that was created for you.<br>")
        htmlfile.write("<img src=\"CreatedProfile.png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\"></td>")
            
        if Xzoom>0:
            htmlfile.write("<li> The top subplot shows original and smoothed bathymetry.<br><li> The bottom subplot shows a zoom-in on intertidal and backshore profiles.<br>")
            htmlfile.write("<td>The figure below shows a zoom-in on the bathymetry and inter- to supratidal portions.<br>")
            htmlfile.write("<img src=\"ZoomIns.png\" alt=\"Profile Plot #2\" width=\"640\" height=\"480\"></td></tr></table><br>")
            
        htmlfile.write("<u><td>Additional information about your site from your inputs:<br></u>")
        if BackHelp==1:
            Foreshore=str(int(Slope))
            DuneH=str(round(DuneCrest,1))
            BermH=str(round(BermCrest,1))
            BermW=str(round(BermLength,1))
            htmlfile.write("<li> The foreshore slope is: 1/"+Foreshore+"<br>")
            htmlfile.write("<li> The dune at your site is is: "+DuneH+"m high<br>")
        
        if any(SlopeM+OffX+ShoreX):
            htmlfile.write("<td>Slope Modifications:<br>")
            if abs(OffX[0]+ShoreX[0])<>0:
                htmlfile.write("<li>     You have a slope of 1:"+str(SlopeM[0])+" for "+str(abs(ShoreX[0]-OffX[0]))+"m, from "+str(OffX[0])+" to "+str(ShoreX[0])+"m<br>")
            if abs(OffX[1]+ShoreX[1])<>0:
                htmlfile.write("<li>     You have a slope of 1:"+str(SlopeM[1])+" for "+str(abs(ShoreX[1]-OffX[1]))+"m, from "+str(OffX[1])+" to "+str(ShoreX[1])+"m<br>")
            if abs(OffX[2]+ShoreX[2])<>0:
                htmlfile.write("<li>     You have a slope of 1:"+str(SlopeM[2])+" for "+str(abs(ShoreX[2]-OffX[2]))+"m, from "+str(OffX[2])+" to "+str(ShoreX[2])+"m<br>")
        
        if ProfileQuestion=="(1) Yes": # profile was cut with GIS
            htmlfile.write("<p>")
            htmlfile.write("The figure below shows the profile that was cut from GIS,<br>including the bathymetry and topography at your location of interest.<br>")
            htmlfile.write("<img src=\"GISProfile.png\" alt=\"Profile Plot #3\" width=\"640\" height=\"240\">")
        
            if HabDirectory and 'HabCount'  in locals(): 
                htmlfile.write("<hr><H2><u>Location of Natural Habitats</u></H2>")
                htmlfile.write("<table border=\"0\" width=\"1100\" cellpadding=\"5\" cellspacing=\"10\">")
                htmlfile.write("<tr><td><img src=\"HabFigure.png\" alt=\"Location of Natural Habitats\" width=\"640\" height=\"480\"></td><td>")
                htmlfile.write("We indicate below the type, start and end locations of the habitats at your site.  Distances are referenced in meters from the shoreline.\
                                Positive distances are oriented seaward and negative distances landward. In other words, if a distance is positive, your habitat is in the water, and if a distance is negative, your habitat in on land.<p>")
                
                # open Excel sheet to write extent of habitats
                xlApp=Dispatch("Excel.Application")
                xlApp.Workbooks.Open(InputTable)
                cell=xlApp.Worksheets("ModelInput")
                
                if any(MG1):
                    htmlfile.write("You have a mangrove field that starts and ends at the locations indicated below:<br>")
                    if 'mangwater' in locals():
                        htmlfile.write("<b>Your mangrove starts in subtidal areas.  Please double check your habitat input layer.</b>")
                    B=begMGx1;F=finMGx1
                    # need to substract "shift" from coordinate because the profile Dx,Dmeas shifted from Xorig,Dorig by measure shift -> see how it's defined
                    # there's also the same shift with TempX,TempY which we use to plot the location of the vegetation
                    htmlfile.write("<table border=\"1\" width=\"200\" cellpadding=\"0\" cellspacing=\"0\"><tr>")
                    htmlfile.write("<td><b> </b></td>")               
                    for kk in range(len(B)):
                        htmlfile.write("<td>Field "+str(kk)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>Start [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(B[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>End [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(F[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("</table>")
                    htmlfile.write("<p>")
                    
                    # write values in Excel
                    cell.Range("i49").Value=int(B[0]) # beginning of the field
                    cell.Range("j49").Value=int(F[-1]) # end of the field

                if any(DN1):
                    htmlfile.write("You have a dune field that starts and ends at the locations indicated below:<br>")
                    if 'dunewater' in locals():
                        htmlfile.write("<b>Your dune starts in subtidal areas.  Please double check your habitat input layer.</b>")
                    B=begDNx1;F=finDNx1
                    htmlfile.write("<table border=\"1\" width=\"200\" cellpadding=\"0\" cellspacing=\"0\"><tr>")
                    htmlfile.write("<td><b> </b></td>")               
                    for kk in range(len(B)):
                        htmlfile.write("<td>Field "+str(kk)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>Start [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(B[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>End [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(F[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("</table>")
                    htmlfile.write("<p>")
                    
                if any(MR1):
                    htmlfile.write("You have a marsh that starts and ends at the locations indicated below:<br>")
                    if 'marwater' in locals():
                        htmlfile.write("<b>Your marsh starts in subtidal areas.  Please double check your habitat input layer.</b>")
                    B=begMRx1;F=finMRx1
                    htmlfile.write("<table border=\"1\" width=\"200\" cellpadding=\"0\" cellspacing=\"0\"><tr>")
                    htmlfile.write("<td><b> </b></td>")               
                    for kk in range(len(B)):
                        htmlfile.write("<td>Field "+str(kk)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>Start [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(B[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>End [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(F[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("</table>")
                    htmlfile.write("<p>")
                    
                    # write values in Excel
                    cell.Range("i52").Value=int(B[0]) # beginning of the field
                    cell.Range("j52").Value=int(F[-1]) # end of the field                 
    
                if any(SG1): # if seagrass pre-management exists
                    htmlfile.write("You have a seagrass bed that starts and ends at the locations indicated below:<br>")
                    if 'seagwater' in locals():
                        htmlfile.write("<b>Your seagrass is on land.  Please double check your habitat input layer.</b>")
                    B=begSGx1;F=finSGx1
                    htmlfile.write("<table border=\"1\" width=\"200\" cellpadding=\"0\" cellspacing=\"0\"><tr>")
                    htmlfile.write("<td><b> </b></td>")               
                    for kk in range(len(B)):
                        htmlfile.write("<td>Field "+str(kk)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>Start [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(B[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>End [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(F[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("</table>")
                    htmlfile.write("<p>")
                    
                    # write values in Excel
                    cell.Range("i53").Value=int(B[0]) # beginning of the field
                    cell.Range("j53").Value=int(F[-1]) # end of the field      
    
                if any(CR1):
                    htmlfile.write("You have a coral reef that starts and ends at the locations indicated below:<br>")
                    if 'coralwater' in locals():
                        htmlfile.write("<b>Your coral reef is on land.  Please double check your habitat input layer.</b>")
                    B=begCRx1;F=finCRx1
                    htmlfile.write("<table border=\"1\" width=\"200\" cellpadding=\"0\" cellspacing=\"0\"><tr>")
                    htmlfile.write("<td><b> </b></td>")               
                    for kk in range(len(B)):
                        htmlfile.write("<td>Field "+str(kk)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>Start [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(B[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>End [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(F[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("</table>")
                    htmlfile.write("<p>")
                    
                    # write values in Excel
                    cell.Range("e56").Value=int(B[0]) # beginning of the field
                    cell.Range("f56").Value=int(F[-1]) # end of the field                 
                                 
                if any(OT1):
                    htmlfile.write("You have another habitat (OT) that starts and ends at the locations indicated below:<br>")
                    B=begOTx1;F=finOTx1
                    htmlfile.write("<table border=\"1\" width=\"200\" cellpadding=\"0\" cellspacing=\"0\"><tr>")
                    htmlfile.write("<td><b> </b></td>")               
                    for kk in range(len(B)):
                        htmlfile.write("<td>Field "+str(kk)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>Start [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(B[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<td>End [m]</td>")
                    for kk in range(len(B)):
                        htmlfile.write("<td>"+str(F[kk])+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("</table>")
                    htmlfile.write("<p>")
                        
                htmlfile.write("</td></tr></table>")
               
                xlApp.ActiveWorkbook.Close(SaveChanges=1) # save changes in Excel
                xlApp.Quit()
            
            elif HabDirectory and 'HabCount' not in locals():
                htmlfile.write("<hr><H2><u>Location of Natural Habitats</u></H2>")
                htmlfile.write("No habitat was identified near your profile.<p>")
        htmlfile.close() # close HTML
       
        if ProfileQuestion=="(1) Yes":
            htmlfile=open(FetchWindWave_HTML,"w")
            htmlfile.write("<html><title>Marine InVEST - Profile Generator</title><CENTER><H1>Coastal Protection - Tier 1</H1><H2>Profile Generator Results ("+subwsStr+")<br></H2><p>")
            htmlfile.write("<a href=\"profile.html\">[ PROFILE INFO ]</a><br>[ FETCH, WAVE, AND WIND INFO ]<br>")
            htmlfile.write("<CENTER><br><HR><H2><u>Fetch and Wind Roses</u></H2></CENTER>")
            htmlfile.write("<table border=\"0\" width=\"1100\" cellpadding=\"5\" cellspacing=\"0\"><tr><td>")
            htmlfile.write("<iframe width=\"350\" height=\"325\" frameborder=\"0\" scrolling=\"no\" marginheight=\"0\" marginwidth=\"0\"")
            htmlfile.write("src=\"http://maps.google.com/maps?f=q&amp;source=s_q&amp;hl=en&amp;geocode=&amp;q=")
            htmlfile.write(PtLat+","+PtLong)
            htmlfile.write("&amp;aq=&amp;sspn=0.009467,0.021136&amp;vpsrc=6&amp;ie=UTF8&amp;ll=")
            htmlfile.write(PtLat+","+PtLong)
            htmlfile.write("&amp;spn=0.063714,0.169086&amp;t=h&amp;z=12&amp;output=embed\"></iframe>")
            htmlfile.write("</td><td>")
            htmlfile.write("<img src=\"Fetch_Plot.png\" width=\"347\" height=\"260\" alt=\"Fetch Rose: No Fetch Calculation Selected\"></td><td>")
            htmlfile.write("<img src=\"Wind_Plot.png\" width=\"347\" height=\"260\" alt=\"Wind Rose: No Wave Watch III Info Provided\"></td></tr></table>")
            
            temp=0
            if WW3_Pts or FetchQuestion=='(1) Yes':
                htmlfile.write("<HR><H2><u>Wind and Wave Information</u></H2>")
                
                if WW3_Pts:
                    htmlfile.write("<i>From Wave Watch III data, we estimated various wind speed values<br>that we used to generate waves from each of the 16 fetch directions.</i><br>")
                htmlfile.write("<table border=\"1\" width=\"900\" cellpadding=\"6\" cellspacing=\"0\"><tr>")
                htmlfile.write("</td><th colspan=\"1\"></th><th colspan=\"16\"><b>Direction (degrees)</b></th></tr>")
                htmlfile.write("<tr align=\"center\"><td></td>")
                for kk in range(0,16):
                    htmlfile.write("<td><b>"+str(dirList[kk])+"°</b></td>")
                htmlfile.write("</tr><tr align=\"center\"><td><b>Fetch (m)</b></td>")
                
                if FetchQuestion=='(1) Yes':
                    for kk in range(0,16):
                        htmlfile.write("<td>"+str(int(FetchFinalList[kk]))+"</td>")
                else:
                    for kk in range(0,16):
                        htmlfile.write("<td>-</td>")
                if WW3_Pts:
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b><FONT COLOR=\"980000\">Max Wind Speed (m/s)</FONT></b></td>")
                    for kk in range(0,16):
                        htmlfile.write("<td>"+str(WiMax[kk])+"</td>")
                        
                if WW3_Pts and FetchQuestion=='(1) Yes':
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b>Wind-Wave Height (m)</b></td>")
                    for kk in range(0,16):
                        temp1=round(WiWavMax[kk],2)
                        if temp1==0.0:
                            temp1=0 
                        htmlfile.write("<td>"+str(temp1)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b>Wind-Wave Period (s)</b></td>")
                    for kk in range(0,16):
                        temp2=round(WiPerMax[kk],2)
                        if temp2==0.0:
                            temp2=0     
                        htmlfile.write("<td>"+str(temp2)+"</td>")
                        
                if WW3_Pts:            
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b><FONT COLOR=\"C80000\">Top 5% Wind Speed (m/s)</FONT></b></td>")
                    for kk in range(0,16):
                        htmlfile.write("<td>"+str(Wi10[kk])+"</td>")
                        
                if WW3_Pts and FetchQuestion=='(1) Yes':
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b>Wind-Wave Height (m)</b></td>")
                    for kk in range(0,16):
                        temp1=round(WiWav10[kk],2)
                        if temp1==0.0:
                            temp1=0 
                        htmlfile.write("<td>"+str(temp1)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b>Wind-Wave Period (s)</b></td>")
                    for kk in range(0,16):
                        temp2=round(WiPer10[kk],2)
                        if temp2==0.0:
                            temp2=0     
                        htmlfile.write("<td>"+str(temp2)+"</td>")
                        
                if WW3_Pts:            
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b><FONT COLOR=\"FF0000\">Top 10% Wind Speed (m/s)</FONT></b></td>")
                    for kk in range(0,16):
                        htmlfile.write("<td>"+str(Wi25[kk])+"</td>")
                                
                if WW3_Pts and FetchQuestion=='(1) Yes':
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b>Wind-Wave Height (m)</b></td>")
                    for kk in range(0,16):
                        temp1=round(WiWav25[kk],2)
                        if temp1==0.0:
                            temp1=0 
                        htmlfile.write("<td>"+str(temp1)+"</td>")
                    htmlfile.write("</tr>")
                    htmlfile.write("<tr align=\"center\"><td><b>Wind-Wave Period (s)</b></td>")
                    for kk in range(0,16):
                        temp2=round(WiPer25[kk],2)
                        if temp2==0.0:
                            temp2=0     
                        htmlfile.write("<td>"+str(temp2)+"</td>")
                htmlfile.write("</tr></table><p>")
                
                if WW3_Pts and FetchQuestion=='(1) Yes':
                    htmlfile.write("<u>Summary</u>:<br>")
                    temp=[WiPerMax[ii]*WiWavMax[ii]**2 for ii in range(len(WiPerMax))];loc=argmax(num.array(temp))
                    htmlfile.write("<li> The most powerful wave generated by the maximum wind speed is: Ho="+str(round(WiWavMax[loc],2))+"m, with a period of To="+str(round(max(WiPerMax),2))+"s<br>")
                    temp=[WiPer10[ii]*WiWav10[ii]**2 for ii in range(len(WiPer10))];loc=argmax(num.array(temp))
                    htmlfile.write("<li> The most powerful wave generated by the top 5% wind speed is: Ho="+str(round(WiWav10[loc],2))+"m, with a period of To="+str(round(max(WiPer10),2))+"s<br>")
                    temp=[WiPer25[ii]*WiWav10[ii]**2 for ii in range(len(WiPer25))];loc=argmax(num.array(temp))
                    htmlfile.write("<li> The most powerful wave generated by the top 10% wind speed is: Ho="+str(round(WiWav25[loc],2))+"m, with a period of To="+str(round(max(WiPer25),2))+"s<p><br>")
            else:
                htmlfile.write("<p><i>We cannot provide you with wave information at your site since you didn't include Wave Watch III in the analysis.</i><p>")
                temp=1
            
            if WW3_Pts:
                htmlfile.write("<b>Wave Height Data from Wave Watch III </b></br>")
                htmlfile.write("<i>From WaveWatchIII data, we estimated various <br>wave height/period values that you could use as input into the wave model</i>")
                htmlfile.write("<table border=\"1\" width=\"600\" cellpadding=\"6\" cellspacing=\"0\">")
                htmlfile.write("<tr align=\"center\"></td><th colspan=\"1\"></th><td>Wave Height (m)</td><td>Wave Period (s)</td></tr><tr align=\"center\">")
                htmlfile.write("<td><b>Maximum Wave</b></td>")
                htmlfile.write("<td>"+str(WavMax[0])+"</td>")
                htmlfile.write("<td>"+str(WavMax[1])+"</td>")
                htmlfile.write("</tr><tr align=\"center\">")
                htmlfile.write("<td><b>Top 5% Wave</b></td>")
                htmlfile.write("<td>"+str(Wav10[0])+"</td>")
                htmlfile.write("<td>"+str(Wav10[1])+"</td>")
                htmlfile.write("</tr><tr align=\"center\">")
                htmlfile.write("<td><b>Top 10% Wave</b></td>")
                htmlfile.write("<td>"+str(Wav25[0])+"</td>")
                htmlfile.write("<td>"+str(Wav25[1])+"</td>")
                htmlfile.write("</tr>")
                if Tm<>-1:
                    htmlfile.write("<tr align=\"center\"><td><b>Most Frequent</b></td>")
                    htmlfile.write("<td>"+str(Hm)+"</td>")
                    htmlfile.write("<td>"+str(Tm)+"</td>")
                htmlfile.write("</table>")
            elif temp==0:
                htmlfile.write("<p><i>We cannot provide you with wave information at your site since you didn't include Wave Watch III in the analysis.</i><p>")
            htmlfile.close() # close HTML

    except:
        gp.AddError(msgHTMLOutputs)
        raise Exception

    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile=open(subws+"parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("PROFILE GENERATOR PARAMETERS\n")
    parafile.writelines("____________________________\n\n")
         
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

    # delete superfluous intermediate data
    if ProfileQuestion=="(1) Yes":
        del1=[LandPointLyr,LandPolyLyr,LandPoint_Buff,LandPoint_Buff50k,LandPoint_Geo,Shoreline,Shoreline_Buff_Clip,Shoreline_Buff_Clip_Diss]
        del2=[PT1,PT2,PT1_Z,PT2_Z,PT_Z_Near,Backshore_Pts,Profile_Pts_Merge,Profile_Pts_Lyr]
        del3=[PtsCopy,PtsCopy2,PtsCopyLR,PtsCopy3,PtsCopyLD,Fetch_AOI,UnionFC,SeaPoly,seapoly_rst,seapoly_e,PtsCopyEL,PtsCopyExp,PtsCopyExp_Lyr,LandPoint_WW3]
        del4=[]
        if HabDirectory:
            for i in range(0,len(HabAbbrevList)):
                del4.append(interws+HabAbbrevList[i])
                del4.append(interws+HabAbbrevList[i]+"_rc")
                del4.append(interws+HabAbbrevList[i]+".shp")
        deletelist = del1+del2+del3+del4
        for data in deletelist:
            if gp.exists(data):
                gp.delete_management(data)
        del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())