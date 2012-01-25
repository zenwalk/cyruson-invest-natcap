# Marine InVEST: Coastal Protection (Wave and Erosion)
# Authors: Greg Guannel, Gregg Verutes, Apollo Yi
# 12/19/11

import sys,os,string,time,datetime,shlex
from math import *
import fpformat,operator
import arcgisscripting

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
msgReadExcel = "\nError reading Wave Erosion Excel file inputs."
msgReadCSProfile = "\nError reading cross-shore profile."
msgCreateInputs = "\nError creating inputs for model run."
msgManagementAction = "\nError running model for management action."
msgPlotErosion = "\nError plotting erosion."
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
    from scipy.special import erf
    from scipy import optimize
    from scipy.integrate import trapz
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
        parameters.append("Label for Erosion Run (10 characters max): "+subwsStr)
        InputTable=gp.GetParameterAsText(2)
        parameters.append("Nearshore Waves and Erosion Excel Table: "+InputTable)
        CSProfile=gp.GetParameterAsText(3)
        parameters.append("Cross-Shore Profile: "+CSProfile)
        WaveErosionQuestion=gp.GetParameterAsText(4)
        parameters.append("Do you have wave height and wave period values? "+WaveErosionQuestion)
        Ho=gp.GetParameterAsText(5)
        parameters.append("IF 1: Wave Height (meters): "+str(Ho))
        if Ho:
            Ho=float(gp.GetParameterAsText(5))
        To=gp.GetParameterAsText(6)
        parameters.append("IF 1: Wave Period (seconds): "+str(To))
        if To:
            To=float(gp.GetParameterAsText(6))
        Us=gp.GetParameterAsText(7)
        parameters.append("IF 2: Wind Speed (meters per second): "+str(Us))
        if Us:
            Us=float(gp.GetParameterAsText(7))
        Ft=gp.GetParameterAsText(8)
        parameters.append("IF 2: Fetch Distance (meters): "+str(Ft))
        if Ft:
            Ft=float(gp.GetParameterAsText(8))
        depth=gp.GetParameterAsText(9)
        parameters.append("IF 2: Water Depth (meters): "+str(depth))
        if depth:
            depth=float(gp.GetParameterAsText(9))
        StormDur=float(gp.GetParameterAsText(10))
        parameters.append("Storm Duration (hours): "+str(StormDur))
        S=float(gp.GetParameterAsText(11))
        parameters.append("Surge Elevation (meters): "+str(S))
        dx=float(gp.GetParameterAsText(12))
        parameters.append("Model Spatial Resolution (dx): "+str(dx))

    except:
        raise Exception,msgArguments+gp.GetMessages(2)


    ###############################################           
    ###### CREATE DIRECTORIES AND VARIABLES #######
    ###############################################

    # remove spaces and shorten 'subwsStr' if greater than 10 characters
    subwsStr=subwsStr.replace(" ","")
    subwsStr=subwsStr[0:10]

    # intermediate and output directories
    outputws=gp.workspace+os.sep+"_WaveModel_Outputs"+os.sep+subwsStr+os.sep
    interws=gp.workspace+os.sep+"scratch"+os.sep
    gp.scratchworkspace=interws

    try:
        thefolders=["_WaveModel_Outputs"]
        for folder in thefolders:
            if not gp.exists(gp.workspace+folder):
                gp.CreateFolder_management(gp.workspace+os.sep,folder)

        thefolders=[subwsStr]
        for folder in thefolders:
            if not gp.exists(gp.workspace+os.sep+"_WaveModel_Outputs"+os.sep+folder):
                gp.CreateFolder_management(gp.workspace+os.sep+"_WaveModel_Outputs"+os.sep,folder)

    except:
        raise Exception,"Error creating folders"


    #################################################
    ########## VARIOUS FUNCTIONS AND CHECKS #########
    #################################################

    try:
        # check that all associated inputs are specified for 'WaveErosionQuestion'
        Answer1=[Ho,To]
        Answer2=[Us,Ft,depth]
        if WaveErosionQuestion == "(1) Yes, I have these values":
            for x1 in Answer1:
                if x1 == '':
                    gp.AddError("\nOne or more of the required input parameters was not defined.\
                                 \nPlease provide values for either Wave Height and/or Wave Period.")
                    raise Exception
        else:
            for x2 in Answer2:
                if x2 == '':
                    gp.AddError("\nOne or more of the required input parameters was not defined.\
                                 \nPlease provide values for either Wind Speed, Fetch Distance, and/or Water Depth.")
                    raise Exception
    except:
        gp.AddError(msgCheckInputs)
        raise Exception


    def gradient(f,z):
        length=len(f)
        df=length*[0.0]
        df[0]=(f[1]-f[0])/z
        for i in range(1,length-1):
            df[i]=(f[i+1]-f[i-1])/(2.0*z)
        df[length-1]=(f[length-1]-f[length-2])/z
        return df

    def Indexed(x,value): # locates index of point in vector x that has closest value as variable value
        mylist=abs(x-value);    
        if isinstance(x,num.ndarray):
            mylist=mylist.tolist()
        minval=min(mylist)
        ind=[i for i,v in enumerate(mylist) if v==minval]
        ind=ind[0]
        return ind

    def FindRootKD(fun,a,b,tol=1e-16):
        a=float(a);b=float(b)
        assert(sign(fun(a)) != sign(fun(b)))
        c=(a+b)/2
        while math.fabs(fun(c))>tol:
            if a==c or b==c:
                break
            if sign(fun(c))==sign(fun(b)):
                b=c
            else:
                a=c
            c=(a+b)/2
        return c

    # dispersion relationship
    def iterativek(sigma,dh):
        kestimated=(sigma**2)/(g*(sqrt(tanh((sigma**2)*dh/g))))
        kprevious=0.0000001
        count=0
        while (abs(kestimated-kprevious)>0.000005) and (count<1000):
            count += 1
            kh=kestimated*dh
            kcalculated=(sigma**2)/(tanh(kh)*g)
            kprevious=kestimated
            kestimated=kcalculated
        qk=kcalculated
        return qk

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

    # wave transformation model
    def WaveModel(X,h,Ho,To,Seagr,Marsh,Mang,Cf):
        global LagKick
        # constants
        g=9.81;rho=1024.0;B=1.0;Beta=0.05
        lx=len(X);dx=abs(X[1]-X[0])
        # initialize vegetation vectors,zeros everywhere except where veg present
        VegLoc=X*0.0;
        
        # seagrass
        veg=Seagr[0];loc=Seagr[1]
        hs=veg[0];ds=veg[1];Ns=veg[2] #height,diam,dens.
        alphs=hs/h;o=find(alphs>=1);alphs[o]=1 #relative height
        off=Indexed(X,loc[0]);sho=Indexed(X,loc[1]) # locate loc of seagrass
        if sho<>off:
            VegLoc[off:sho+1]=1 # 1 for location of seagrass

        # marsh
        veg=Marsh[0];loc=Marsh[1]
        hr=veg[0];dr=veg[1];Nr=veg[2]#height,diam,dens.
        alphr=hr/h;o=find(alphr>=1);alphr[o]=1#relative height
        off=Indexed(X,loc[0]);sho=Indexed(X,loc[1])# locate loc of marsh
        if sho<>off:
            VegLoc[off:sho+1]=2 #2 for location of marshes

        # mangroves
        veg=Mang[0];loc=Mang[1]
        rt=veg[0];tk=veg[1];cn=veg[2] # roots,trunk and canopy info
        hgr=rt[0];dgr=rt[1];Ngr=rt[2]
        hgt=tk[0];dgt=tk[1];Ngt=tk[2]
        hgc=cn[0];dgc=cn[1];Ngc=cn[2]
        off=Indexed(X,loc[0]);sho=Indexed(X,loc[1])# locate loc of mangr
        if sho<>off:
            VegLoc[off:sho+1]=3 # 3 for location of mangroves
           
        alphgr=hgr/h;alphgt=hgt/h;alphgc=hgc/h
        for kk in range(lx): #Create relative depth values for roots,trunk and canopy
            if alphgr[kk]>1:
                alphgr[kk]=1;alphgt[kk]=0;alphgc[kk]=0 #Roots only
            elif alphgr[kk]+alphgt[kk]>1:
                alphgt[kk]=1-alphgr[kk];alphgc[kk]=0 #Roots and trunk
            elif alphgr[kk]+alphgt[kk]+alphgc[kk]>1:
                alphgc[kk]=1-alphgr[kk]-alphgt[kk] #Roots, trunk and canopy
                
        # drag coefficent for vegetation; mangrove and marsh win over seagrass if they overlap
        Cds=num.arange(0.0,lx,dx)*0.0;o=find(VegLoc==1);Cds[o]=0.1 # seagrass	
        Cdr=num.arange(0.0,lx,dx)*0.0;o=find(VegLoc==2);Cdr[o]=0.1 # marsh	
        Cdg=num.arange(0.0,lx,dx)*0.0;o=find(VegLoc==3);Cdg[o]=1 # mangrove

        # initialize vectors for wave model
        H=lx*[0.0];Eta=lx*[0.0]
        Db=lx*[0.0];Df=lx*[0.0]
        k=lx*[0.0];L=lx*[0.0]
        C=lx*[0.0];n=lx*[0.0];Cg=lx*[0.0]
        Ds=lx*[0.0];Dr=lx*[0.0];Dg=lx*[0.0] # seagrass, marsh, mangrove
        Er=lx*[0.0];Ef=lx*[0.0];Br=lx*[0.0]
        Hs=lx*[0.0];Etas=lx*[0.0]
        Dbs=lx*[0.0];Dfs=lx*[0.0]
        Ers=lx*[0.0];Efs=lx*[0.0];Brs=lx*[0.0]
        
        # wave parameter at 1st grid pt
        ash=[h[ii] for ii in range(lx)] # ash is same as h,but is now an independent variable
        fp=1.0/To; sig=2.0*pi*fp
        k[0]=iterativek(sig,h[0]) # wave number at 1st grid pt
        L[0]=2.0*pi/k[0] # wave length @ 1st grid pt
        n[0]=0.5*(1+(2.0*k[0]*h[0]/sinh(2.0*k[0]*h[0]))) # to compute Cg at 1st grid pt
        C[0]=L[0]/To;Cg[0]=C[0]*n[0] # phase and group velocity at 1st grid pt
        So=Ho/L[0] # deep water wave steepness
        Gam=0.5+0.4*num.tanh(33.0*So) # Gam from Battjes and Stive 85,as per Alsina & Baldock

        # RMS wave height at first grid point; assume no dissipation occurs
        Ho=Ho/sqrt(2.0) # transform significant wave height to rms wave height
        Co=g*To/(2.0*pi); Cgo=Co/2.0 # deep water phase and group speed
        if h[0]>0.5*L[0]: 
            H[0]=Ho # we are in deep water
        else: 
            H[0]=Ho*sqrt(Cgo/Cg[0]) # we are in intermediate water. Assume no brkg occured

        if LagKick==1: #No shoaling if wave goes from coral reef to lagoon
            H[0]=Ho
            
        Ef[0]=0.125*rho*g*H[0]**2*Cg[0];Efs[0]=Ef[0];Hs[0]=H[0] # energy flux @ 1st grid pt

        # begin wave model 
        for xx in range(lx-1): # transform waves,take MWL into account
            # wave in presence of veg.
            Ef[xx]=0.125*rho*g*(H[xx]**2.0)*Cg[xx] # Ef at (xx)      
            Ef[xx+1]=Ef[xx]-dx*(Db[xx]+Df[xx]+Ds[xx]+Dr[xx]+Dg[xx]) # Ef at [xx+1] 

            k[xx+1]=iterativek(sig,h[xx+1])
            n[xx+1]=0.5*(1.0+(2.0*k[xx+1]*h[xx+1]/sinh(2.0*k[xx+1]*h[xx+1])))
            C[xx+1]=sig/k[xx+1];Cg[xx+1]=C[xx+1]*n[xx+1] # phase and group velocity

            H[xx+1]=num.sqrt(8.0*Ef[xx+1]/(rho*g*Cg[xx+1])) # wave height at [xx+1]      
            Br[xx+1]=Br[xx]-dx*(g*Er[xx]*sin(Beta)/C[xx]-0.5*Db[xx]) # roller flux
            Er[xx+1]=Br[xx+1]/(C[xx+1]) # roller energy

            Var=0.25*rho*g*fp*B
            Hb=0.88/k[xx+1]*tanh(Gam*k[xx+1]*h[xx+1]/0.88)		
            temp1=((Hb/H[xx+1])**3.0+1.5*Hb/H[xx+1])*num.exp(-(Hb/H[xx+1])**2.0)
            if isreal(H[xx+1]): temp2=0.75*num.sqrt(pi)*(1-erf(Hb/H[xx+1]))
            else: temp2=0
            Db[xx+1]=Var*H[xx+1]**3/h[xx+1]*(temp1+temp2) # dissipation due to brkg

            Df[xx+1]=rho*Cf[xx+1]/(12.0*pi)*(2.0*pi*fp*H[xx+1]/sinh(k[xx+1]*h[xx+1]))**3.0 # diss due to bot friction 

            CdDN=Cds[xx+1]*ds*Ns
            temp1=rho*CdDN*(k[xx+1]*g/(2.0*sig))**3.0/(2.0*sqrt(pi))
            temp2=sinh(k[xx+1]*alphs[xx+1]*h[xx+1])**3.0+3.0*sinh(k[xx+1]*alphs[xx+1]*h[xx+1])
            temp3=(3.0*k[xx+1]*cosh(k[xx+1]*h[xx+1])**3.0)
            Ds[xx+1]=temp1*temp2/temp3*H[xx+1]**3.0 # diss due to seagrass
            
            CdDN=Cdr[xx+1]*dr*Nr
            temp1=rho*CdDN*(k[xx+1]*g/(2.0*sig))**3.0/(2.0*sqrt(pi))
            temp2=sinh(k[xx+1]*alphr[xx+1]*h[xx+1])**3.0+3.0*sinh(k[xx+1]*alphr[xx+1]*h[xx+1])
            temp3=(3.0*k[xx+1]*cosh(k[xx+1]*h[xx+1])**3.0)
            Dr[xx+1]=temp1*temp2/temp3*H[xx+1]**3.0 # diss due to marshes

            V1=3*sinh(k[xx+1]*alphgr[xx+1]*h[xx+1])+sinh(k[xx+1]*alphgr[xx+1]*h[xx+1])**3 # roots
            V2=(3*sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1])*h[xx+1])-3*sinh(k[xx+1]*alphgr[xx+1]*h[xx+1])+
                sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1])*h[xx+1])**3-
                sinh(k[xx+1]*alphgr[xx+1]*h[xx+1])**3) # trunk
            V3=(3*sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1]+alphgc[xx+1])*h[xx+1])
                -3*sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1])*h[xx+1])+
                sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1]+alphgc[xx+1])*h[xx+1])**3-
                sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1])*h[xx+1])**3) # canopy
            
            CdDN=Cdg[xx+1]*(dgr*Ngr*V1+dgt*Ngt*V2+dgc*Ngc*V3)
            temp1=rho*CdDN*(k[xx+1]*g/(2.0*sig))**3.0/(2.0*sqrt(pi))
            temp3=(3.0*k[xx+1]*cosh(k[xx+1]*h[xx+1])**3)
            Dg[xx+1]=temp1/temp3*H[xx+1]**3
            
            # wave in absence of veg
            Hs[xx+1]=H[xx+1]
            if sum(VegLoc)<>0:
                Efs[xx]=0.125*rho*g*(Hs[xx]**2.0)*Cg[xx] # Ef at (xx)      
                Efs[xx+1]=Efs[xx]-dx*(Dbs[xx]+Dfs[xx]) # Ef at [xx+1] 
            
                Hs[xx+1]=num.sqrt(8.0*Efs[xx+1]/(rho*g*Cg[xx+1])) # wave height at [xx+1]      
                Brs[xx+1]=Brs[xx]-dx*(g*Ers[xx]*sin(Beta)/C[xx]-0.5*Dbs[xx]) # roller flux
                Ers[xx+1]=Brs[xx+1]/(C[xx+1]) # roller energy
        
                temp1=((Hb/Hs[xx+1])**3.0+1.5*Hb/Hs[xx+1])*num.exp(-(Hb/Hs[xx+1])**2.0)
                if isreal(Hs[xx+1]):  temp2=0.75*num.sqrt(pi)*(1-erf(Hb/Hs[xx+1]))
                else: temp2=0
                Dbs[xx+1]=Var*Hs[xx+1]**3/h[xx+1]*(temp1+temp2) # dissipation due to brkg
        
                Dfs[xx+1]=rho*Cf[xx+1]/(12.0*pi)*(2.0*pi*fp*Hs[xx+1]/sinh(k[xx+1]*h[xx+1]))**3.0 # diss due to bot friction 
        
        Ew=lx*[0.0];Ew=[0.125*rho*g*(H[i]**2.0) for i in range(lx)]
        Ews=lx*[0.0];Ews=[0.125*rho*g*(Hs[i]**2.0) for i in range(lx)]

        # estimate MWS
        val=1.0;vals=1.0;dell=1.0;dells=1.0;Os=0.0;O=0.0 # for setup calc
        Sxx=lx*[0.0]; Rxx=lx*[0.0]; Sxxs=lx*[0.0]; Rxxs=lx*[0.0]
        
        # force on plants if they were emergent; take a portion if plants occupy only portion of wc
        Fxs=[rho*g*Cds[i]*ds*Ns*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
        fxs=[-alphs[i]*Fxs[i] for i in range(lx)] # seagrass
        
        Fxr=[rho*g*Cdr[i]*dr*Nr*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
        fxr=[-alphr[i]*Fxr[i] for i in range(lx)] # marsh

        Fxgr=[rho*g*Cdg[i]*dgr*Ngr*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
        Fxgt=[rho*g*Cdg[i]*dgt*Ngt*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
        Fxgc=[rho*g*Cdg[i]*dgc*Ngc*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
        fxg=[-alphgr[i]*Fxgr[i]-alphgt[i]*Fxgt[i]-alphgc[i]*Fxgc[i] for i in range(lx)] # mangrove

        fx=[fxs[ii]+fxr[ii]+fxg[ii] for ii in range(lx)] # all
        
        while dell>1e-10: # iterate until convergence of water level
            val1=val;h=[ash[i]+Eta[i] for i in range(lx)] # water depth
            
            Sxx=[0.5*Ew[i]*(4.0*k[i]*h[i]/sinh(2.0*k[i]*h[i])+1.0) for i in range(lx)] # wave radiation stress
            Rxx=[2.0*Er[i] for i in range(lx)] # roller radiation stress
            # estimate MWL along Xshore transect
            temp1=[Sxx[i]+Rxx[i] for i in range(lx)]
            temp2=gradient(temp1,dx)
            
            Integr=[(-temp2[i]+fx[i])/(rho*g*h[i]) for i in range(lx)]
            Eta[0]=-0.125*H[0]**2.0*k[0]/sinh(2.0*k[0]*h[0])
            Eta[1]=Eta[0]+Integr[0]*dx
            
            #Wave in absence of veg
            Etas[0]=Eta[0];Etas[1]=Eta[1]
            if sum(VegLoc)<>0:
                val1s=vals;hs=[ash[i]+Etas[i] for i in range(lx)] # water depth
                
                Sxxs=[0.5*Ews[i]*(4.0*k[i]*hs[i]/sinh(2.0*k[i]*hs[i])+1.0) for i in range(lx)] # wave radiation stress
                Rxxs=[2.0*Ers[i] for i in range(lx)] # roller radiation stress
                # estimate MWL along Xshore transect
                temp1s=[Sxxs[i]+Rxxs[i] for i in range(lx)]
                temp2s=gradient(temp1s,dx)
                
                Integrs=[(-temp2s[i])/(rho*g*hs[i]) for i in range(lx)]
                Etas[0]=-0.125*Hs[0]**2.0*k[0]/sinh(2.0*k[0]*hs[0])
                Etas[1]=Etas[0]+Integrs[0]*dx

            for i in range(1,lx-2):
                Eta[i+1]=Eta[i-1]+Integr[i]*2*dx
                Etas[i+1]=Eta[i+1]
                if sum(VegLoc)<>0:
                    Etas[i+1]=Etas[i-1]+Integrs[i]*2*dx
            
            Eta[lx-1]=Eta[lx-2]+Integr[lx-1]*dx
            Etas[lx-1]=Eta[lx-1]
            if sum(VegLoc)<>0:
                Etas[lx-1]=Etas[lx-2]+Integrs[lx-1]*dx
        
            temp=gradient(Eta,dx)
            valn=[temp2[i]+rho*g*h[i]*temp[i]-fx[i] for i in range(lx)]
            valn[lx-3]=0;valn[lx-2]=0;valn[lx-1]=0
            dell1=[valn[i]-val1 for i in range(lx)]
            dell1me=sum(dell1)/len(dell1)
            dell1mx=max(dell1)
            temp1=max(dell1me,dell1mx)

            temp2=temp1
            if sum(VegLoc)<>0:
                temp=gradient(Etas,dx)
                valn=[temp2s[i]+rho*g*hs[i]*temp[i] for i in range(lx)]
                valn[lx-3]=0;valn[lx-2]=0;valn[lx-1]=0
                dell1=[valn[i]-val1 for i in range(lx)]
                dell1me=sum(dell1)/len(dell1)
                dell1mx=max(dell1)
                temp2=max(dell1me,dell1mx)
            O=O+1
            if O<3: dell=13
            else: dell=max(temp1,temp2)

        Ubot=[pi*H[ii]/(To*sinh(k[ii]*h[ii])) for ii in range(lx)] # bottom velocity
        H=[H[ii]*sqrt(2) for ii in range(lx)]
        Hs=[Hs[ii]*sqrt(2) for ii in range(lx)]
        
        Ds1=num.array(H);Ds2=num.array(Hs)
        Diss=mean(Ds1[find(VegLoc>0)]**3/(Ds2[find(VegLoc>0)]**3)) # dissipation difference
        if sum(VegLoc)==0:    Diss=1 # if no veg., then ratio=1

        return H,Eta,Hs,Etas,Diss,Ubot,VegLoc

    def WaveModelSimple(X,h,Ho,To):
        global LagKick
        # constants
        g=9.81;rho=1024.0;B=1.0;Beta=0.05
        lx=len(X);dx=X[1]-X[0];
        
        # initialize vectors for wave model
        H=lx*[0.0];Eta=lx*[0.0];L=lx*[0.0]
        Db=lx*[0.0];Df=lx*[0.0]
        Er=lx*[0.0];Ef=lx*[0.0];Br=lx*[0.0]
        C=lx*[0.0];n=lx*[0.0];Cg=lx*[0.0]
        k=lx*[0.0];
        
        # wave parameter at 1st grid pt
        ash=[h[ii] for ii in range(lx)] # ash is same as h,but is now an independent variable.
        fp=1.0/To; sig=2.0*pi*fp
        k[0]=iterativek(sig,h[0]) # wave number at 1st grid pt
        L[0]=2.0*pi/k[0] # wave length @ 1st grid pt
        n[0]=0.5*(1+(2.0*k[0]*h[0]/sinh(2.0*k[0]*h[0]))) # to compute Cg at 1st grid pt
        C[0]=L[0]/To;Cg[0]=C[0]*n[0] # phase and group velocity at 1st grid pt

        # RMS wave height at first grid point; Assume no dissipation occurs
        Ho=Ho/sqrt(2.0) # transform significant wave height to rms wave height
        Co=g*To/(2.0*pi); Cgo=Co/2.0 # deep water phase and group speed
        if h[0]>0.5*L[0]: 
            H[0]=Ho # we are in deep water
        else: 
            H[0]=Ho*sqrt(Cgo/Cg[0]) # we are in intermediate water. Assume no brkg occured

        if LagKick==1: # no shoaling if wave goes from coral reef to lagoon
            H[0]=Ho
            
        Ef[0]=0.125*rho*g*H[0]**2*Cg[0] # energy flux @ 1st grid pt
        So=Ho/L[0] # deep water wave steepness
        Gam=0.5+0.4*num.tanh(33.0*So) # Gam from Battjes and Stive 85, as per Alsina & Baldock
        
        # begin wave model 
        for xx in range(lx-1): # transform waves,take MWL into account
            Ef[xx]=0.125*rho*g*(H[xx]**2.0)*Cg[xx] # Ef at (xx)      
            Ef[xx+1]=Ef[xx]-dx*(Db[xx]+Df[xx]) # Ef at [xx+1] 

            k[xx+1]=iterativek(sig,h[xx+1])
            n[xx+1]=0.5*(1.0+(2.0*k[xx+1]*h[xx+1]/sinh(2.0*k[xx+1]*h[xx+1])))
            C[xx+1]=sig/k[xx+1];Cg[xx+1]=C[xx+1]*n[xx+1] # phase and group velocity

            H[xx+1]=num.sqrt(8.0*Ef[xx+1]/(rho*g*Cg[xx+1])) # wave height at [xx+1]      
            Br[xx+1]=Br[xx]-dx*(g*Er[xx]*sin(Beta)/C[xx]-0.5*Db[xx]) # roller flux
            Er[xx+1]=Br[xx+1]/(C[xx+1]) # roller energy

            Var=0.25*rho*g*fp*B
            Hb=0.88/k[xx+1]*tanh(Gam*k[xx+1]*h[xx+1]/0.88)		
            temp1=((Hb/H[xx+1])**3.0+1.5*Hb/H[xx+1])*num.exp(-(Hb/H[xx+1])**2.0)
            if isreal(H[xx+1]): temp2=0.75*num.sqrt(pi)*(1-erf(Hb/H[xx+1]))
            else: temp2=0
            Db[xx+1]=Var*H[xx+1]**3/h[xx+1]*(temp1+temp2) # dissipation due to brkg

            Df[xx+1]=rho*0.01/(12.0*pi)*(2.0*pi*fp*H[xx+1]/sinh(k[xx+1]*h[xx+1]))**3.0 # diss due to bot friction 

        Ew=lx*[0.0];Ew=[0.125*rho*g*(H[i]**2.0) for i in range(lx)]

        # estimate MWS
        val=1.0;dell=1.0;O=0.0 # for setup calc
        Sxx=lx*[0.0]; Rxx=lx*[0.0]
        
        while dell>1e-10: # iterate until convergence of water level
            val1=val
            h=[ash[i]+Eta[i] for i in range(lx)] # water depth
            
            Sxx=[0.5*Ew[i]*(4.0*k[i]*h[i]/sinh(2.0*k[i]*h[i])+1.0) for i in range(lx)] # wave radiation stress
            Rxx=[2.0*Er[i] for i in range(lx)] # roller radiation stress
            # estimate MWL along Xshore transect
            Tem=[Sxx[i]+Rxx[i] for i in range(lx)]
            Temp=gradient(Tem,dx)
            Terms=[-Temp[i] for i in range(lx)]
            Integr=[Terms[i]/(rho*g*h[i]) for i in range(lx)]
            Eta[0]=-0.125*H[0]**2.0*k[0]/sinh(2.0*k[0]*h[0])
            Eta[1]=Eta[0]+Integr[0]*dx

            for i in range(1,lx-2):
                Eta[i+1]=Eta[i-1]+Integr[i]*2*dx
            
            Eta[lx-1]=Eta[lx-2]+Integr[lx-1]*dx
        
            Temp_eta=gradient(Eta,dx)
            valn=[Temp[i]+rho*g*h[i]*Temp_eta[i] for i in range(lx)]
            valn[lx-3]=0;valn[lx-2]=0;valn[lx-1]=0
            dell1=[valn[i]-val1 for i in range(lx)]
            dell1me=sum(dell1)/len(dell1)
            dell1mx=max(dell1)
            
            O=O+1
            if O<5: dell=13
            else: dell=max(dell1me,dell1mx)

        H=[H[ii]*sqrt(2) for ii in range(lx)]

        return H,Eta

    # wave attenuation by coral reefs
    def WavesCoral(Ho,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf,dx):
        BrkRim=0;BrkFace=0;Kp=0.0
        
        if AlphF+AlphR+he==0: # user doesn't have data
            Kp=mean(kp) # assume waves break on face
        else:
            # reef profile
            Xreef=num.arange(0.0,10000.0,dx)
            Yreef=AlphR*Xreef-he
            loc=find(Yreef>-hr)
            Yreef[loc]=-hr # reef profile
        
            D=To*sqrt(g/hr) # relative subm
            Lo=g*To**2.0/(2.0*pi)
            if D==8:  # first check for breaking location
                if Ho>=0.5*he:
                    BrkFace=1; # wave break on face
                    loc=Indexed(TanAlph,AlphF)
                    Kp=kp[loc] # reef shape factor
        
                    ha=hr # rep. depth
                    
                elif Ho>=0.5*hr:
                    BrkRim=1 # wave break on rim
                    loc=Indexed(TanAlph,AlphR)
                    Kp=kp2[loc] # reef shape factor
        
                    db=0.259*Ho*(tan(AlphR)**2*Ho/Lo)**(-0.17) # breaking wave height
                    if db<hr: db=hr
                    elif db>he: db=he
        
                    xs=(2+1.1*Ho/he)*To*(g*he)**.5 # surf zone width
                    loc1=Indexed(Yreef,-db); # start brkg
                    loc2=Indexed(Xreef,Xreef[loc1]+xs) # end brkg
                    ha=-average(Yreef[loc1:loc2+1]) # rep. depth
            else:  # second check for breaking location
                if Ho>=0.4*he:
                    BrkFace=1 # wave break on face
                    loc=Indexed(TanAlph,AlphF)
                    Kp=kp[loc] # reef shape factor
                    ha=hr
                    
                elif Ho>=0.4*hr:
                    BrkRim=1 # wave break on rim
                    loc=Indexed(TanAlph,AlphR)
                    Kp=kp2[loc] # reef shape factor
        
                    db=0.259*Ho*(tan(AlphR)**2.0*Ho/Lo)**(-0.17) # breaking wave height
                    if db<hr: db=hr
                    elif db>he: db=he
        
                    xs=(2+1.1*Ho/he)*To*(g*he)**.5 # surf zone width
                    loc1=Indexed(Yreef,-db) # start brkg
                    loc2=Indexed(Xreef,Xreef[loc1]+xs) # end brkg
                    ha=-average(Yreef[loc1:loc2+1]) # rep. depth

        # wave transformation on top of coral reef
        Etar=0.0;delta=10;
        Xflat=num.arange(0.0,Wr,dx)
        lx=len(Xflat);
        H_r=lx*[Ho];

        if Kp<> 0: # wave break on reef face or rim
            # wave Setup
            while delta>0.01:
                EtaR=3.0/(64.0*pi)*Kp*Ho**2*To*g**.5/((Etar+ha)**1.5)
                delta=abs(EtaR-Etar);Etar=EtaR
            Etar=Etar[0]    
            H_r[0]=0.42*(hr+Etar) # wave height at the offshore edge of reef
        else:
            H_r[0]=Ho
            
        for xx in range(lx-1):
            H_r[xx+1]=H_r[xx]-Cf*H_r[xx]**2.0/(3.0*pi*(Etar+hr)**2.0)*dx
            
        return H_r,Etar

    # wave attenuation by reef breakwater
    def BreakwaterKt(Hi,To,hi,hc,Cwidth,Bwidth):
        Lo=9.81*To**2.0/(2.0*pi)
        Rc=hc-hi # depth of submergence
        Boff=(Bwidth-Cwidth)/2.0 # base dif on each side
        ksi=(hc/Boff)/sqrt(Hi/Lo)
        
        if Cwidth<>0: # it's not a reef ball
            # van der Meer (2005)
            Kt1=-0.4*Rc/Hi+0.64*(Cwidth/Hi)**(-.31)*(1.0-num.exp(-0.5*ksi)) # transmission coeff: d'Angremond
            Kt1=max(Kt1,0.075);Kt1=min(0.8,Kt1);
        
            Kt2=-0.35*Rc/Hi+0.51*(Cwidth/Hi)**(-.65)*(1.0-num.exp(-0.41*ksi)) # transmission coeff: van der Meer
            Kt2=max(Kt2,0.05);Kt2=min(Kt2,-0.006*Cwidth/Hi+0.93)
        
            if Cwidth/Hi<8.0: # d'Angremond
                Kt=Kt1
            elif Cwidth/Hi>12.0: # van der Meer
                Kt=Kt2
            else: # linear interp
                temp1=(Kt2-Kt1)/4.0;temp2=Kt2-temp1*12.0
                Kt=temp1*Cwidth/Hi+temp2
                
        else: # it's a reef ball
            Kt=1.616-31.322*Hi/(9.81*To**2)-1.099*hc/hi+0.265*hi/Bwidth
            if hc>hi:
                Kt=1;
                gp.AddWarning("The reef balls are emerged. We cannot compute the wave transmission coefficient.")
            
        return Kt

    # K&D Erosion model
    def ErosionKD(A,Ho,TWL,B,D,W,m):
        global BetaKD
        # constants
        BD=D+B;
        if TWL>BD:
            TWL=BD
            gp.AddMessage("Water level is higher than your backshore elevation. You will probably experience flooding.")
            
        # Erosion model 1
        hb=(((Ho**2.0)*g*To/(2*pi))/2.0)**(2.0/5.0)/(g**(1.0/5.0)*0.73**(4.0/5.0))
        xb=(hb/A)**1.5;Hb=0.78*hb
        
        Term1=xb-hb/m
        Rinf=(TWL*Term1)/(BD+hb-TWL/2)-W*(B+hb+0.5*TWL)/(BD+hb-0.5*TWL) # potential erosion distance
        if Rinf<0:
            gp.AddMessage("Berm is so wide that dune is not eroding.")
            Rinf=(TWL*Term1)/(BD+hb-TWL/2)
            
        TS=(320.0*(Hb**(3.0/2.0)/(g**.5*A**3.0))*(1.0/(1.0+hb/BD+(m*xb)/hb)))/3600.0 # response time scale
        BetaKD=2.0*pi*(TS/StormDur)
        expr="num.exp(-2.0*x/BetaKD)-num.cos(2.0*x)+(1.0/BetaKD)*num.sin(2.0*x)" # solve this numerically
        fn=eval("lambda x: "+expr)
        z=FindRootKD(fn,pi,pi/2) # find zero in function,initial guess from K&D
        R01=0.5*Rinf*(1.0-cos(2.0*z)) # final erosion distance
        
        # 2nd method
        x=num.arange(0,10000,1)
        y=A*x**(2.0/3);y=y[::-1]
        y=y[find(y>0.5)];x=x[find(y>0.5)]    
        Hs,Etas=WaveModelSimple(x,y,Ho,To) 
        hb=y[argmin(Etas)]
        xb=(hb/A)**1.5 # surf zone width
        
        # erosion model
        Term1=xb-hb/m
        Rinf=(TWL*Term1)/(BD+hb-TWL/2)-W*(B+hb+0.5*TWL)/(BD+hb-0.5*TWL) # potential erosion distance
        if Rinf<0:
            gp.AddMessage("Berm is so wide that dune is not eroding.")
            Rinf=(TWL*Term1)/(BD+hb-TWL/2)
            
        TS=(320.0*(Hb**(3.0/2.0)/(g**.5*A**3.0))*(1.0/(1.0+hb/BD+(m*xb)/hb)))/3600.0 # response time scale
        BetaKD=2.0*pi*(TS/StormDur)
        expr="num.exp(-2.0*x/BetaKD)-num.cos(2.0*x)+(1.0/BetaKD)*num.sin(2.0*x)" # solve this numerically
        fn=eval("lambda x: "+expr)
        z=FindRootKD(fn,pi,pi/2) # find zero in function,initial guess from K&D
        R02=0.5*Rinf*(1.0-cos(2.0*z)) # final erosion distance
        
        R0=max([R01,R02])
        if R0<0:    R0=0

        return R0,R01,R02

    def MudErosion(Uc,Uw,h,To,me,Cm,values):
        rho=1024.0;nu=1.36e-6;d50=0.03
        ks=2.5*d50;kap=0.4
        
        # current
        Uc=0.0 # current speed before management action
        if Uc<>0:
            us1=0.01;zo1=0.01;dif=10 #Initial value for u* and zo  
            Tc=h*0; # shear stress
            for xx in range(values):
                while dif>1e-4:
                    zo2=ks/30*(1-num.exp(-us1*ks/(27*nu)))+nu/(9*us1)
                    us2=kap*Uc/(num.log(h[xx]/zo2)-1)
                    dif1=abs(us1-us2);dif2=abs(zo1-zo2);dif=dif1+dif2
                    zo1=zo2;us1=us2;
                Tc[xx]=rho*us1**2 # shear stress due to current
        else:    Tc=h*0
        
        # waves
        Rw=Uw**2*To/(2*num.pi)/nu
        fw=0.0521*Rw**(-0.187) #Smooth turbulent flow
        Tw=0.5*rho*fw*Uw**2
        
        # combined Wave and Current
        temp=Tc*(1+1.2*(Tw/(Tc+Tw))**3.2)
        Trms=(temp**2+0.5*Tw**2)**0.5
        
        # erosion 
        Te=h*0+0.0012*Cm**1.2 # erosion threshold
        dmdt=me*(Trms-Te) # erosion rate
        dmdt[find(dmdt<=0)]=0
        Erosion=3600*dmdt/Cm*100 # rate of bed erosion [cm/hr]
        
        return Erosion,Trms,Tc,Tw,Te


    ################################################
    #### READING DEPTH PROFILE AND EXCEL INPUTS ####
    ################################################

    try:
        gp.AddMessage("\nReading Excel inputs and depth profile...")
        # constants
        g=9.81;rho=1024.0
        Fig2=0; # for HTML

        # input from user via Excel SS
        xlApp=Dispatch("Excel.Application")
        xlApp.Visible=0
        xlApp.DisplayAlerts=0
        xlApp.Workbooks.Open(InputTable)
        cell1=xlApp.Worksheets("WaveModelInput")
        cell3=xlApp.Worksheets("ReefShapeFactor")

        # read general information
        MSL=cell1.Range("f16").Value #Mean Sea Level
        MHW=cell1.Range("g16").Value #Mean high water

        temp=cell1.Range("h20").Value
        if temp==1:
            sand=1
            mud=0
        elif temp==2:
            sand=0
            mud=1

        # read muddy shoreline information
        Cm=cell1.Range("h24").Value #dry density
        me=cell1.Range("i24").Value #Erosion constant
        # read beach Information
        d50=cell1.Range("e24").Value # sediment size

        temp=cell1.Range("c41:h41").Value
        temp=temp[0] # all beach information
        Slope=temp[1]
        if Slope<>0:
            m=1.0/Slope # foreshore slope=1/Slope
        else:
            m=0
        A=cell1.Range("f55").Value # sediment scale factor

        B1=temp[2]
        W1=temp[3]
        D1=temp[4]
        Dred=temp[5]
        if sand+mud<>1:
            gp.AddMessage("You didn't specify the sediment type in the backshore area.  We won't be able to estimate amount of erosion.")

        # read Oyster Reef Information
        temp=cell1.Range("d49:g49").value;temp=temp[0] #All oyster information
        Xr=temp[0] # distance from shoreline; need to reverse because user enter with shoreline at X=0
        hc=temp[1] # reef height
        Bw=temp[2] # base width
        Cw=temp[3] # crest width

        # read Vegetation Information
        temp=cell1.Range("f33:k37").Value # all vegetation information

        temp1=temp[2] # mangrove canopy 
        hogc=temp1[0] # height of mangrove canopy
        dogc=temp1[1] # diameter of mangrove canopy
        Nogc=temp1[2] # density of mangrove canopy
        temp1=temp[1] # mangrove trunks
        hogt=temp1[0] # height of mangrove trunks
        dogt=temp1[1] # diameter of mangrove trunks
        Nogt=temp1[2] # density of mangrove trunks
        temp1=temp[0] # mangrove roots
        hogr=temp1[0] # height of mangrove roots
        dogr=temp1[1] # diameter of mangrove roots
        Nogr=temp1[2] # density of mangrove roots
        Xo1g=temp1[4] # offshore edge of mangrove
        Xo2g=temp1[3] # shoreward edge of mangrove
        MangMngt=temp1[5] # management action for mangrove
        if MangMngt is None:
            MangMngt="None"
        if hogc+dogc+Nogc+hogt+dogt+Nogt+hogr+dogr+Nogr==0:
            MangMngt="None"

        temp1=temp[3]
        hos=temp1[0] # height of seagrass
        dos=temp1[1] # diameter of seagrass
        Nos=temp1[2] # density of seagrass
        Xo1s=temp1[4] # offshore edge of seagrass
        Xo2s=temp1[3] # shoreward edge of seagrass
        SeagMngt=temp1[5] # management action for seagrass
        if SeagMngt is None:
            SeagMngt="None"
        if hos+dos+Nos==0:
            SeagMngt="None"

        temp1=temp[4]
        hor=temp1[0] # height of marsh
        dor=temp1[1] # diameter of marsh
        Nor=temp1[2] # density of marsh
        Xo1r=temp1[4] # offshore edge of marsh
        Xo2r=temp1[3] # shoreward edge of marsh
        MarshMngt=temp1[5] # management action for seagrass
        if MarshMngt is None:
            MarshMngt="None"
        if hor+dor+Nor==0:
            MarshMngt="None"

        # read coral info
        temp=cell1.Range("d45:k45").Value;temp=temp[0] # all coral info
        Xco=temp[1] #off. distance
        Xcn=temp[0] # nearshore distance
            
        AlphF=temp[2] # reef face slope
        AlphR=temp[3]    # reef rim slope
        he=temp[4] # reef rim edge
        hr=temp[5] # reef top depth
        Wr=temp[6] # reef top width
        CoralMngt=temp[7]
        TanAlph=cell3.Range("a2:a202").Value;TanAlph=num.array(TanAlph)
        kp=cell3.Range("b2:b202").Value;kp=num.array(kp) # reef shape factor
        kp2=cell3.Range("c2:c202").Value;kp2=num.array(kp2) # reef shape factor
        if CoralMngt is None:
            CoralMngt="None"
        if AlphF+AlphR+he+hr+Wr==0:
            CoralMngt="None"

        xlApp.ActiveWorkbook.Close(SaveChanges=0)
        xlApp.Quit()
        
    except:
        xlApp.ActiveWorkbook.Close(SaveChanges=0)
        xlApp.Quit()
        gp.AddError(msgReadExcel)
        raise Exception

    try:
        # read in user's cross-shore profile
        TextData=open(CSProfile,"r") 
        X_lst=[];h_lst=[]
        for line in TextData.read().strip("\n").split("\n"):
            linelist=[float(s) for s in line.split("\t")] # split the list by comma delimiter
            X_lst.append(linelist[0])
            h_lst.append(linelist[1])
        TextData.close()
        X=num.array(X_lst);h=num.array(h_lst)
        if h[0]>0 or h[0]>h[-1]:
            h=h[::-1] # reverse order if profile starts at shoreline

    except:
        gp.AddError(msgReadCSProfile)
        raise Exception


    #########################################
    ######## CREATING INPUTS FOR RUN ########
    #########################################

    try:
        # compute offshorewave height
        if WaveErosionQuestion=="(2) No, please compute these values from wind speed and fetch distance": 
            Us=Us;Ft=Ft;depth=depth
            Ho,To=WindWave(Us,Ft,depth) # compute wave from wind speed
        gp.AddMessage('\nInput conditions are: Ho = ' +str(round(Ho,2)) +'m and To = ' +str(round(To,2)) +'s')

        if Ho>0.78*(-h[0]):
            Ho=0.78*(-h[0])
            gp.AddWarning("Water depth too shallow for input wave height. That wave broke somewhere in deeper water. We will assume that H = 0.78h")

        # fix veg. start and end point
        if (Xco+Xcn)>0.0:    Xco=X[-1]-Xco;Xcn=X[-1]-Xcn #Coral
        Xo1s=X[-1]-Xo1s;Xo2s=X[-1]-Xo2s # seagrass
        Xo1r=X[-1]-Xo1r;Xo2r=X[-1]-Xo2r # marsh
        Xo1g=X[-1]-Xo1g;Xo2g=X[-1]-Xo2g # mangrove
        Xr=X[-1]-Xr # oyster

        # adjust Bathy
        Zero=Indexed(h,0);ash=-h #Locate MWL 
        h=h-MSL; # adjust water level so that 0 is at MSL
        X=X-X[0]
        # modify the depth profile appropriately
        if mud==1: # if there's a marsh or a mangrove
            h=h-S # add surge level
            keep=find(h<-0.2)
            h=h[keep];X=X[keep] # only keep values that are below water
            dx=X[1]-X[0];m=abs(h[-1]-h[-int(10.0/dx)])/10 # average slope 10m from end of transect
        else: # it's a beach
            keep=find(h<-0.2)
            h=h[keep];X=X[keep] # only keep values that are below water

        # resample bathy to dx
        F=interp1d(X,h);X=num.arange(X[0],X[-1]+dx,dx)
        h=F(X);lx=len(X) #Resample to dx

        # make sure vegetation location is correct
        if Xco>X[-1]:    Xco=X[-1]
        if Xcn>X[-1]:    Xcn=X[-1]
        if Xo1s>X[-1]:    Xo1s=X[-1]
        if Xo2s>X[-1]:    Xo2s=X[-1]
        if Xo1r>X[-1]:    Xo1r=X[-1]
        if Xo2r>X[-1]:    Xo2r=X[-1]
        if Xo1g>X[-1]:    Xo1g=X[-1]
        if Xo2g>X[-1]:    Xo2g=X[-1]
        if Xr>X[-1]:    Xr=X[-1]

        # management actions
        gp.AddMessage("\nEcosystems present at your site...")
        # beach
        D2=(100-Dred)*D1/100 # new dune height
        if A>0 and B1+W1+D1<>0:
            Beach=1 # a sandy beach is present
        else:
            Beach=0
            
        if Beach+sand==2:
            gp.AddMessage("...a sandy beach is present")
        else:
            gp.AddMessage("...you do NOT have a sandy beach\n   The model will NOT compute erosion for the type of sediment that you have.")

        # prepare oyster data
        if hc+Bw<> 0:
            Oyster=1 # oyster info was entered
            if Cw<>0:
                gp.AddMessage("...an oyster reef is present")
            else:
                gp.AddMessage("...reef balls are present")        
        else:
            Oyster=0
            gp.AddMessage("...you do NOT have an oyster reef")

        # prepare mangrove data
        Mang0=[[[0,0,0],[0,0,0],[0,0,0]],[0,0]];Mang=[[[0,0,0],[0,0,0],[0,0,0]],[0,0]]
        Mang1=[[[0,0,0],[0,0,0],[0,0,0]],[0,0]]
        loc=[Xo1g,Xo2g]
        if loc is not None and (hogc+dogc+Nogc+hogt+dogt+Nogt+hogr+dogr+Nogr)<> 0:
            veg=[[hogr,dogr,Nogr],[hogt,dogt,Nogt],[hogc,dogc,Nogc]]
            Mang=[veg,loc] # mangrove characteristics
            
            if MangMngt=="None": # management action
                Mang1=Mang
            elif MangMngt=="Rmv":
                Mang1=Mang0
            elif MangMngt=="Half":
                veg=[[hogr,dogr,Nogr/2.0],[hogt,dogt,Nogt/2.0],[hogc,dogc,Nogc]];#Rdce dens.roots&trees by 1/2
                Mang1=[veg,loc]
                
        MgPst=sum(sum(Mang[0]))+sum(sum(Mang1[0]))
        if MgPst<>0:  
            MgPst=1
            gp.AddMessage("...a mangrove is present")
        else:
            gp.AddMessage("...you do NOT have a mangrove")
            
        # prepare seagrass data
        Seagr0=[[0,0,0],[0,0]];Seagr=[[0,0,0],[0,0]];Seagr1=[[0,0,0],[0,0]]
        loc=[Xo1s,Xo2s]
        if loc is not None and (hos+dos+Nos)<> 0:
            veg=[hos,dos,Nos]
            Seagr=[veg,loc] # seagrass characteristics
            
            if SeagMngt=="None": # management action
                Seagr1=Seagr
            elif SeagMngt=="Rmv":
                Seagr1=Seagr0
            elif SeagMngt=="Half":
                veg=[hos,dos,Nos/2]
                Seagr1=[veg,loc] 
                
        SgPst=sum(Seagr[0])+sum(Seagr[1])+sum(Seagr1[0])+sum(Seagr1[1])
        if SgPst<>0:  
            SgPst=1
            gp.AddMessage("...a seagrass bed is present")
        else:
            gp.AddMessage("...you do NOT have a seagrass bed")

        # prepare marsh data
        Marsh=[[0,0,0],[0,0]];Marsh0=[[0,0,0],[0,0]];Marsh1=[[0,0,0],[0,0]]
        loc=[Xo1r,Xo2r]
        if loc is not None and (hor+dor+Nor)<> 0:
            veg=[hor,dor,Nor]
            Marsh=[veg,loc] # marsh characteristics

            if MarshMngt=="None": # management action
                Marsh1=Marsh
            elif MarshMngt=="Rmv":
                Marsh1=Marsh0
            elif MarshMngt=="Half":
                veg=[hor,dor,Nor/2.0]
                Marsh1=[veg,loc] 
                
        MrPst=sum(Marsh[0])+sum(Marsh[1])+sum(Marsh1[0])+sum(Marsh1[1])
        if MrPst<>0:    
            MrPst=1 
            gp.AddMessage("...a marsh is present")
        else:
            gp.AddMessage("...you do NOT have a marsh")

        # prepare coral data
        if AlphF+AlphR+he+hr+Wr<>0: # check if coral profile is imported
            Coral=1
            gp.AddMessage("...a coral reef is present\n   We will estimate its profile for you.")
        if Xcn+Xco<>0 and AlphF+AlphR+he+hr+Wr==0: # check if coral profile is part of the bathy
            Coral=-1
            gp.AddMessage("...a coral reef is present\n   It is incorporated in the bathy profile you uploaded.")
        elif Xcn+Xco+AlphF+AlphR+he+hr+Wr==0:    
            Coral=0
            gp.AddMessage("...you do NOT have a coral reef")

        Cf_bed=num.arange(0,len(h),.5)*0+0.01 # bottom friction for sand only
        Cf1_bed=num.arange(0,len(h),.5)*0+0.01  # for management action
        Cf =0.2;Cf1=0.2; # friction for live and dead coral

        temp1=find(X>=Xco);temp2=find(X>Xcn)
        Loco=[ii for ii in range(temp1[0],temp2[0],1)] # location of coral reef in profile
        if len(Loco)>2:   Cf_bed[Loco]=0.2 # live coral

        if Coral==1:
            if CoralMngt=="Rmv":        Cf1=.01 # no coral
            elif CoralMngt=="Dead":        Cf1=.1 # dead coral
            elif CoralMngt=="None":        Cf1=.2 # live coral
        elif Coral==-1 and len(Loco)>2:
            if CoralMngt=="Rmv":        Cf1_bed[Loco]=.01 # no coral
            elif CoralMngt=="Dead":        Cf1_bed[Loco]=.1 # dead coral
            elif CoralMngt=="None":        Cf1_bed[Loco]=0.2 # live coral

    except:
        gp.AddError(msgCreateInputs)
        raise Exception

                
    ################################################
    ####### RUN MODEL FOR MANAGEMENT ACTION ########
    ################################################

    try:    
        # BEFORE MANAGEMENT ACTION
        gp.AddMessage("\nComputing wave height profiles...")
        gp.AddMessage("...before management action")

        h=-h # depth is now positive
        Ho1=Ho # offshore wave height
        LagKick=0 # check for lagoon calc
        H_r=[];Xn=[];ho=[];Hsimple1=[];Etasimple1=[]

        # estimate initial wave height if coral reef at beg of transect and no bathy for that reef
        if Xco==0.0 and Xcn==0.0 and Coral==1: # reef top depth
            H_r,Eta_r=WavesCoral(Ho,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf,dx)
            H_r=num.array(H_r);Eta_r=H_r*0+Eta_r # arrays of H_r and Eta_r
            Ho1=H_r[-1];ho=h[0] # wave height at the end of the reef
            Xn=num.arange(-Wr,0,dx)
            LagKick=1 # for future    

        # estimate wave height transformation over bathy profile
        if len(h)>2: 
            H,Eta,Hs,Etas,DissAtn1,Ubot1,VegLoc1=WaveModel(X,h,Ho1,To,Seagr,Marsh,Mang,Cf_bed)

        if Xo1s>Xcn and (MgPst+MrPst+SgPst)<>0: # wave height in lagoon in absence of any vegetation
            Hsimple1=Hs;Etasimple1=Etas        

        # wave height over coral reefs at various locations
        if len(h)>2 and Xco>=0 and Xcn>0 and Coral==1: # coral Reef in middle of transect
            Ho1=H[Loco[0]] # wave height at edge of the coral reef
            H_r,Eta_r=WavesCoral(Ho1,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf,dx)
            H_r=num.array(H_r);Eta_r=H_r*0+Eta_r # arrays of H_r and Eta_r
            LagKick=1;ho=(h[Loco[0]]+h[Loco[-1]])/2
            Xn=num.arange(0,Wr,dx)
            
            Hr=num.array(len(Loco)*[1.0 ])*[H[Loco[0]]] # constant H over reef face/rim
            Hr[-Wr:len(Hr)]=H_r
            H1=num.append(num.array(H[0:Loco[0]]),Hr,None) # all wave profile to edge of reef
            X1=X[0:Loco[-1]+1] # all X profile to edge of reef

            Etar=num.array(len(Loco)*[1.0 ])*[Eta[Loco[0]]] # constant Eta over reef face/rim
            Etar[-Wr:len(Hr)]=Eta_r # MWL profile over reef
            Eta1=num.append(num.array(Eta[0:Loco[0]]),Etar,None) # MWL to edge of reef

            Ho1=H_r[-1] # wave height at the end of the reef
            temp=X[Loco[-1]:-1] # lagoon X
            Hlag,Etalag,Hs,Etas,DissAtn1,temp2,temp=WaveModel(temp,h[Loco[-1]:-1],Ho1,To,Seagr,Marsh,Mang,Cf_bed[Loco[-1]:-1]) # H in lagoon
            Ubot1[Loco[-1]:-1]=temp2 # bottom velocity
            
            if Xo1s>Xcn and (MgPst+MrPst+SgPst)<>0: # wave height in lagoon in case there is vegetation
                Hsimple1=Hs;Etasimple1=Etas        

        elif len(h)>2 and Xco==1 and Xcn==1 and Coral==1: # coral Reef end of transect
            H_r,Eta_r=WavesCoral(H[-1],AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf,dx)
            H_r=num.array(H_r);Eta_r=H_r*0+Eta_r # arrays of H_r and Eta_r
            Ho1=H_r[0];ho=h[-1]; # wave height at offshore edge of reef
            Xn=num.arange(X[-1],X[-1]+Wr,dx)
                
        elif Oyster==1:
            Xloc=find(X>(Xr)) # location of oyster reef
            Ho1=H[Xloc[0]]
            LagKick=1 # wave height at edge of the reef
            hi=h[Xloc[0]] # water depth at reef location
            Kt=BreakwaterKt(Ho1,To,hi,hc,Cw,Bw);Ho1=Kt*Ho1 # transmitted Wave height
            UbotO1=num.array(Ubot1);Ubot1=num.array(Ubot1)
            if Kt<1:
                Hlag,Etalag,Hs,Etas,DissAtn1,temp2,temp=WaveModel(X[Xloc],h[Xloc],Ho1,To,Seagr,Marsh,Mang,Cf_bed[Xloc])
                temp1=num.array(Hlag);temp2=num.array(H)
                Ubot1[Xloc]=UbotO1[Xloc]*temp1/temp2[Xloc] # bottom velocity    

            if Xo1s>Xr and (MgPst+MrPst+SgPst)<>0: # wave height in lagoon in case there is vegetation
                Hsimple1=Hs;Etasimple1=Etas        

        # put all wave characteristics together
        if Xco==0 and Xcn ==0 and Coral==1: # reef at offshore edge of the profile
            X1=num.append(Xn,X,None)
            H1=num.append(H_r,num.array(H),None)
            Eta1=num.append(Eta_r,num.array(Eta),None)
        elif len(h)>2 and Xco>=0 and Xcn>0 and Coral==1: # reef somewhere in profile
            X1=num.append(X1,X[Loco[-1]:-1],None)
            H1=num.append(num.array(H1),num.array(Hlag),None)
            Eta1=num.append(num.array(Eta1),num.array(Etalag),None)
        elif len(h)>2 and Xco==1 and Xcn==1 and Coral==1: # reef at the nearshore edge of the profile
            X1=num.append(X,Xn,None)
            H1=num.append(num.array(H),H_r,None)
            Eta1=num.append(num.array(Eta),Eta_r,None)
        elif Oyster==1:
            XO1=num.array(X);HO1=num.array(H);EtaO1=num.array(Eta)
            HsimpleO1=Hsimple1;EtasimpleO1=Etasimple1
            if Kt<1:
                X1=num.array(X);H1=num.append(H[0:Xloc[0]],Hlag);Eta1=num.append(Eta[0:Xloc[0]],Etalag)
            elif Kt==1:
                X1=XO1;H1=HO1;Eta1=EtaO1
        else:
            X1=num.array(X);H1=num.array(H);Eta1=num.array(Eta)

        # AFTER MANAGEMENT ACTION
        gp.AddMessage("...after management action")

        Ho2=Ho # offshore wave height
        LagKick=0 # check for lagoon calc

        if SgPst+MrPst+MgPst+Coral<>0: # if management action affects vegetation or coral
            Hsimple2=[];Etasimple2=[]
            
            # estimate initial wave height if coral reef at beg of transect and no bathy for that reef
            if Xco==0 and Xcn==0.0 and Coral==1: # reef top depth
                H_r,Eta_r=WavesCoral(Ho2,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf1,dx)    
                H_r=num.array(H_r);Eta_r=H_r*0+Eta_r # arrays of H_r and Eta_r
                Ho2=H_r[-1] # wave height at the end of the reef
                Xn=num.arange(-Wr,0,dx)
                LagKick=1
            
            # estimate wave height transformation over bathy profile
            if len(h)>2 and CoralMngt+MarshMngt+SeagMngt+MangMngt<>"NoneNoneNoneNone": # H over bathy
                H,Eta,Hs,Etas,DissAtn2,Ubot2,VegLoc2=WaveModel(X,h,Ho2,To,Seagr1,Marsh1,Mang1,Cf1_bed)
                
            if Xo1s>Xcn and (MgPst+MrPst+SgPst)<>0: # estimate wave height in lagoon in case there is seagrass
                Hsimple2=Hs;Etasimple2=Etas
                
            if len(h)>2 and Xco>=0 and Xcn>0 and Coral==1: # coral Reef in middle of transect
                Ho2=H[Loco[0]] # wave height at edge of the coral reef
                H_r,Eta_r=WavesCoral(Ho2,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf1,dx)
                H_r=num.array(H_r);Eta_r=H_r*0+Eta_r # arrays of H_r and Eta_r
                Xn=num.arange(0,Wr,dx);LagKick=1
                
                Hr=num.array(len(Loco)*[1.0 ])*[H[Loco[0]]] # constant H over reef face/rim
                Hr[-Wr:len(Hr)]=H_r
                H2=num.append(num.array(H[0:Loco[0]]),Hr,None) # all wave profile to edge of reef
                X2=X[0:Loco[-1]+1] # all X profile to edge of reef
            
                Etar=num.array(len(Loco)*[1.0 ])*[Eta[Loco[0]]] # constant Eta over reef face/rim
                Etar[-Wr:len(Hr)]=Eta_r # MWL profile over reef
                Eta2=num.append(num.array(Eta[0:Loco[0]]),Etar,None) # MWL to edge of reef
                    
                Ho2=H_r[-1] # wave height at the end of the reef
                temp1=X[Loco[-1]:-1];temp2=h[Loco[-1]:-1]
                Hlag,Etalag,Hs,Etas,DissAtn2,temp2,temp=WaveModel(temp1,temp2,Ho2,To,Seagr1,Marsh1,Mang1,Cf1_bed[Loco[-1]:-1]) 
                Ubot2[Loco[-1]:-1]=temp2 # bottom velocity    
                
                if Xo1s>Xcn and (MgPst+MrPst+SgPst)<>0: # wave height in lagoon in case there is vegetation
                    Hsimple2=Hs;Etasimple2=Etas
            
            elif len(h)>2 and Xco==1 and Xcn==1 and Coral==1: # coral Reef at nearshore edge
                H_r,Eta_r=WavesCoral(H[-1],AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf1,dx)
                H_r=num.array(H_r);Eta_r=H_r*0+Eta_r # arrays of H_r and Eta_r
                Ho2=H_r[0] # wave height at offshore edge of reef
                Xn=num.arange(X[-1],X[-1]+Wr,dx)
             
            # oyster Reef
            elif Oyster==1:
                Xloc=find(X>(Xr)) # location of oyster reef
                Ho2=H[Xloc[0]]; LagKick=1 # wave height at edge of the reef
                hi=h[Xloc[0]] # water depth at reef location
                Kt=BreakwaterKt(Ho2,To,hi,hc,Cw,Bw);Ho2=Kt*Ho2 # transmitted wave height
                UbotO2=num.array(Ubot1);Ubot2=num.array(Ubot2)
                if Kt<1:
                    Hlag,Etalag,Hs,Etas,DissAtn2,temp2,temp=WaveModel(X[Xloc],h[Xloc],Ho2,To,Seagr1,Marsh1,Mang1,Cf1_bed[Xloc]) 
                    temp1=num.array(Hlag);temp2=num.array(H)
                    Ubot2[Xloc]=UbotO2[Xloc]*temp1/temp2[Xloc] # bottom velocity    
            
                if Xo1s>Xr and (MgPst+MrPst+SgPst)<>0: # wave height in lagoon in case there is vegetation
                    Hsimple2=Hs;Etasimple2=Etas
            
            # put all wave charact. together
            if Xco==0 and Xcn ==0 and Coral==1: # reef at offshore edge of the profile
                X2=num.append(Xn,X,None)
                H2=num.append(H_r,num.array(H),None)
                Eta2=num.append(Eta_r,num.array(Eta),None)
            elif len(h)>2 and Xco==1 and Xcn==1 and Coral==1: # reef at the nearshore edge of the profile
                X2=num.append(X,Xn,None)
                H2=num.append(num.array(H),H_r,None)
                Eta2=num.append(num.array(Eta),Eta_r,None)
            elif len(h)>2 and Xco>=0 and Xcn>0 and Coral==1: # reef somewhere in profile
                X2=num.append(X2,X[Loco[-1]:-1],None)
                H2=num.append(num.array(H2),num.array(Hlag),None)
                Eta2=num.append(num.array(Eta2),num.array(Etalag),None)
            elif Oyster==1:
                XO2=num.array(X);HO2=num.array(H);EtaO2=num.array(Eta)
                if Kt<1:
                    X2=num.array(X);H2=num.append(H[0:Xloc[0]],Hlag);Eta2=num.append(Eta[0:Xloc[0]],Etalag)
                elif Kt==1:
                    X2=XO2;H2=HO2;Eta2=EtaO2
            else:
                X2=num.array(X);H2=num.array(H);Eta2=num.array(Eta)
        elif SgPst+MrPst+MgPst+Coral+Oyster==0:
            Hsimple2=num.array(Hsimple1);Etasimple2=num.array(Etasimple1)
            H2=num.array(H1);X2=num.array(X1);Eta2=num.array(Eta1)
            Ubot2=num.array(Ubot1)
        elif SgPst+MrPst+MgPst+Coral==0 and Oyster==1:
            Hsimple2=num.array(HsimpleO1);Etasimple2=num.array(EtasimpleO1)
            H2=num.array(HO1);X2=num.array(XO1);Eta2=num.array(EtaO1)
            Ubot2=num.array(UbotO1)

        # save .txt outputs
        Hb4=outputws+"WaveHeightBefore_"+subwsStr+".txt"
        file=open(Hb4,"w")
        for i in range(0,len(X)):
            file.writelines(str(X1[i])+"\t"+str(H1[i])+"\n")
        file.close()
        Haftr=outputws+"WaveHeightAfter_"+subwsStr+".txt"
        file=open(Haftr,"w")
        for i in range(0,len(X)):
            file.writelines(str(X2[i])+"\t"+str(H2[i])+"\n")
        file.close()

    except:
        gp.AddError(msgManagementAction)
        raise Exception

            
    ##############################
    ####### CREATE OUTPUTS #######
    ##############################

    try:
        gp.AddMessage("...plotting wave profiles")
        
        # percent wave attenuation
        lx=len(X1)
        AtnH=H1*0
        for xx in range(lx):
            AtnH[xx]=abs(H2[xx]-H1[xx])/H1[xx]*100 # percent attenuation
        Atn=find(AtnH>0.5) # take Attn > 0.5% to remove any num. diff.
        if len(Atn)>0:
            Atn=trapz(AtnH[Atn],X1[Atn],dx)/(X1[Atn[-1]]-X1[Atn[0]])
        else:
            Atn=0
            
        # save wave outputs
        WaveHeight=outputws+"WaveHeight_"+subwsStr+".txt"
        file=open(WaveHeight,"w")
        for i in range(0,len(X1)):
            file.writelines(str(X1[i])+"\t"+str(H1[i])+"\t"+str(H2[i])+"\n")
        file.close()

        # indices of vegetation location
        keep=find(VegLoc1<>0)
        VegLoc=VegLoc1+nan
        if len(keep)>1:
            VegLoc[keep]=-h[keep]

        # plot
        figure(1)
        ax=subplot(311);plot(X1[::-1],H1,X2[::-1],H2,linewidth=2);grid()
        box=ax.get_position();
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ylabel('Wave Height[m]',size='large')
        ax.legend(('Initial','Mgmt'),loc='center left', bbox_to_anchor=(.95, 0.5))
        ax=subplot(3,1,2);plot(X1[::-1],AtnH,linewidth=2);grid()
        box=ax.get_position();
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
        ylabel('Wave Attn[%]',size='large')
        if sum(VegLoc1)>0 and Coral<>0:
            if Xco+Xcn>3:
                ax=subplot(313);plot(X[-1]-X[Loco],X[Loco]*0.0-ho+S,'ob',X[::-1],VegLoc+S,'xg',X[::-1],-h+S,linewidth=2);grid()
            elif Xco+Xcn==0:
                temp=num.arange(-Wr,0,1);
                ax=subplot(313);plot(X[-1]-temp,temp*0.0-ho+S,'ob',X[::-1],VegLoc+S,'xg',X[::-1],-h+S,linewidth=2);grid()
            elif Xco+Xcn==2:
                temp=num.arange(X[-1],X[-1]+Wr,1);
                ax=subplot(313);plot(X[-1]-temp,temp*0.0-ho+S,'ob',X[::-1],VegLoc,'xg',X[::-1],-h+S,linewidth=2);grid()        
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ax.legend(('Coral','Vegetation','Bed Profile'),loc='center left', bbox_to_anchor=(.95, 0.5))
        elif sum(VegLoc1)>0 and Coral==0:
            ax=subplot(313);plot(X[::-1],VegLoc+S,'xg',X[::-1],-h+S,linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ax.legend(('Vegetation','Bed Profile'),loc='center left', bbox_to_anchor=(.95, 0.5))
        elif sum(VegLoc1)==0 and Coral<>0:
            ax=subplot(313);plot(X[-1]-(X[Loco[0]]+Xn),num.array(Xn)*0.0-ho+S,'ob',X[::-1],-h+S,linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ax.legend(('Coral','Bed Profile'),loc='center left', bbox_to_anchor=(.95, 0.5))
        else:
            ax=subplot(313);plot(X[::-1],-h+S,linewidth=2);grid()  
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            
        ylabel('Depth[m]',size='large')
        xlabel('Cross-Shore Distance from Offshore Boundary',size='large')
        savefig(outputws+"WavePlot_"+subwsStr+".png",dpi=(640/8))

        if Oyster==1:
            ax=subplot(211);plot(X1[::-1],H1,linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ylabel('Wave Height[m]',size='large')
            
            ax=subplot(212);
            Yrf=num.arange((-hi)+S,0.0,0.05);Xrf=X[-1]-(Yrf*0.0+Xloc[0]) # x-axis
            plot(X[::-1],-h+S,'k',Xrf[::-1],Yrf,'r',Xrf[::-1]+.05,Yrf,'r',Xrf[::-1]+.1,Yrf,'r',linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ax.legend(('Bed Profile','Oyster Reef'),loc='center left', bbox_to_anchor=(.95, 0.5))
            ylabel('Depth[m]',size='large')
            xlabel('Cross-Shore Distance from Offshore Boundary [m]',size='large')
            savefig(outputws+"WavePlot_"+subwsStr+".png",dpi=(640/8))
           
        Fig2=0
        gp.AddMessage("\nGenerating outputs...")
        if Beach==1 and sand==1:
            gp.AddMessage("...estimating erosion amount for sandy beach")
            # estimate runup amount
            Lo=g*To**2.0/(2.0*pi)
            
            # before management action
            Rmax1=1.1*(0.35*m*sqrt(Ho1*Lo)+sqrt(Lo*(Ho1*0.563*m**2.0+Ho*0.004))/2.0) # max runup--short and long waves are different (Stockdon adaptation)
            Emax1=0.35*m*num.sqrt(Ho1*Lo)
            
            Etap1=Emax1;Hp1=Ho1;Rnp1=Rmax1
            if len(Etasimple1)<>0 and (MgPst+MrPst+SgPst)<>0: # if there is vegetation
                coef0=max(Etasimple1[-200:-1])/Emax1;Etap1=max(Eta1[-200:-1])/coef0 # corr. MWL at shoreline; take last 200m
                if Etap1<0:
                    Etap1=0 # if Eta with veg. neg,take as zero
                Hp1=(Etap1/(0.35*m))**2/Lo # Hprime to estimate runup with veg
                Rnp1=1.1*(Etap1+num.sqrt(Lo*(Hp1*0.563*m**2+0.004*Ho))/2) # runup with vegetation
                
            # after manamgement action
            Rmax2=1.1*(0.35*m*sqrt(Ho2*Lo)+sqrt(Lo*(Ho2*0.563*m**2.0+Ho*0.004))/2.0) # runup (short and long waves) 
            Emax2=0.35*m*num.sqrt(Ho2*Lo)
            
            Etap2=Emax2;Hp2=Ho2;Rnp2=Rmax2
            if len(Etasimple2)<>0 and (MgPst+MrPst+SgPst)<>0: # if there is vegetation
                coef0=max(Etasimple2[-200:-1])/Emax2;Etap2=max(Eta2[-200:-1])/coef0 # corr. MWL at shoreline; take last 200m
                if Etap2<0:
                    Etap2=0 # if Eta with veg. neg, take as zero
                Hp2=(Etap2/(0.35*m))**2/Lo # Hprime to estimate runup with veg
                Rnp2=1.1*(Etap2+num.sqrt(Lo*(Hp2*0.563*m**2+0.004*Ho))/2) # runup with vegetation
                
            # check if foreshore slope adequate (use worst wave height)
            hb=(((max([Ho2,Ho1])**2.0)*g*To/(2*pi))/2.0)**(2.0/5.0)/(g**(1.0/5.0)*0.73**(4.0/5.0)) # breaking depth
            xb=(hb/A)**1.5 # surf zone width
            temp=xb-hb/m
            if temp<=0:
                mo=floor(xb/hb)
                gp.AddWarning("Your foreshore slope is too flat or sediment size is too high.\nWe'll increase it to 1/" +str(mo) +" to estimate a minimum erosion value.")
            else:
                mo=m

            # before management action
            R_01,temp1,temp2=ErosionKD(A,Ho1,Rmax1+S,B1,D1,W1,mo) # erosion of beach 
            if len(Etasimple1)<>0 and Rnp1<>Rmax1: # if there is vegetation
                temp1=R_01*Rnp1/Rmax1 # scale by runup
                temp2=R_01*DissAtn1 # scale by dissipation
                if Hp1<>0:
                    temp,temp1,temp2=ErosionKD(A,Hp1,Rnp1+S,B1,D1,W1,mo) # erosion of beach 
                else:
                    temp=mean([temp1,temp2])
            else:
                temp1=R_01;temp=R_01;temp2=R_01
            R1=mean([temp,temp1,temp2])
                
            # after manamgement action
            if Ho2+S+D2==Ho1+S+D1:
                R_02=R_01 # if same forcing, no need to rerun
            else:
                R_02,temp1,temp2=ErosionKD(A,Ho2,Rmax2+S,B1,D2,W1,mo) # erosion of beach 
            if len(Etasimple2)<>0 and Rnp2<>Rmax2: # if there is vegetation
                temp1=R_02*Rnp2/Rmax2 # scale by runup
                temp2=R_02*DissAtn2 # scale by dissipation
                if Hp2<>0:
                    temp,temp1,temp2=ErosionKD(A,Hp2,Rnp2+S,B1,D2,W1,mo) # erosion of beach 
                else:
                    temp=mean([temp1,temp2])
            else:
                temp1=R_02;temp=R_02;temp2=R_02
            R2=mean([temp,temp1,temp2])
                
            # plot outputs
            # beach profile plot
            if D1==0:
                W1=100
            Xp=num.arange(0,1000,1)
            Yp=m*Xp
            loc=find(Yp>=B1)
            Yp[loc]=B1
            Lberm=loc[0];Ltoe=Indexed(Xp,Xp[Lberm+W1])
            Yp[Ltoe:-1]=B1+D1
            
            Yp0=m*Xp-1*m*R1 # profile of erosion before management action
            loc=find(Yp0>=B1)
            if R1 <= W1:
                Yp0[loc]=B1
                Lberm0=loc[0]
                W00=W1-(Xp[Lberm0]-Xp[Lberm])
                Ltoe0=Indexed(Xp,Xp[Lberm0+W00])
                Yp0[Ltoe0:-1]=B1+D1
            else:
                Yp0[loc]=B1+D1
                Ltoe0=loc[0]
            
            Ypv=m*Xp-1*m*R2 # profile of erosion after management action
            loc=find(Ypv>=B1)
            if R2 <= W1:
                Ypv[loc]=B1
                Lbermv=loc[0]
                W0v=W1-(Xp[Lbermv]-Xp[Lberm])
                Ltoev=Indexed(Xp,Xp[Lbermv+W0v])
                Ypv[Ltoev:-1]=B1+D2
            else:
                Ypv[loc]=B1+D2
            
            gp.AddMessage("...plotting erosion profiles")
            figure(2)
            ax=subplot(211)    
            if max([R1,R2])>=W1-1:    
                temp2=Indexed(Xp,Xp[Ltoe0]+max([R1,R2])+.5);temp1=0
            else:     temp2=Indexed(Xp,Xp[Lberm]+max([R1,R2])+5);temp1=Indexed(Xp,Xp[Lberm]-10);
            plot(Xp[::-1],Yp,Xp[::-1],Yp0,'--r',Xp[::-1],Ypv,'--g',linewidth=2);grid()
            xlim(Xp[-1]-Xp[temp2],Xp[-1]-Xp[temp1]);ylim(Yp[temp1],Yp[temp2]+.05)
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ax.legend(('Initital','Bef. Mgmt','After Mgmt'),loc='center left', bbox_to_anchor=(.95, 0.5))
            ylabel('Backsh. Elev.[m]',size='large')
            title('Erosion Profiles',size='large',weight='bold')
            
            ax=subplot(212)
            temp1=Indexed(Xp,Xp[Lberm]-10);temp2=Indexed(Xp,Xp[Ltoe0]+10)
            plot(Xp[::-1],Yp,Xp[::-1],Yp0,'--r',Xp[::-1],Ypv,'--g',linewidth=2);grid()
            xlim(Xp[-1]-Xp[temp2]-R1,Xp[-1]-Xp[temp1]);ylim(Yp[temp1],B1+D1+0.5)
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ylabel('Backsh. Elev.[m]',size='large')
            xlabel('Not to Scale',size='large')
            savefig(outputws+"ErosionBed_"+subwsStr+".png",dpi=(640/8))
            Fig2=1; #for HTML
                    
        elif mud==1: # compute erosion amount for consolidated sediments
            gp.AddMessage("...estimating erosion amount for muddy substrate")
            if Oyster:
                R1,Trms1,Tc1,Tw1,Te=MudErosion(0,Ubot1,h,To,me,Cm,Xloc)
                R2,Trms2,Tc2,Tw2,Te=MudErosion(0,Ubot2,h,To,me,Cm,Xloc)
            
            if MgPst+MrPst<>0:
                lx=len(X)
                Ubot1=num.array(Ubot1);Ubot2=num.array(Ubot2)
                values=range(Zero,lx)
                R1,Trms1,Tc1,Tw1,Te=MudErosion(0,Ubot1,h,To,me,Cm,values) # before management action
                R2,Trms2,Tc2,Tw2,Te=MudErosion(0,Ubot2,h,To,me,Cm,values) # after
                
                if MgPst:
                    loc1=Indexed(X,Xo1g);loc2=Indexed(X,Xo2g) # locate the edges of the mangrove
                elif MrPst:
                    loc1=Indexed(X,Xo1r);loc2=Indexed(X,Xo2r) # locate the edges of the marsh              
                
                gp.AddMessage("...plotting erosion profiles")
                figure(2)
                ax=subplot(311)
                plot(X[-1]-X[Zero:-1],Trms1[Zero:-1],X[-1]-X[Zero:-1],Trms2[Zero:-1],X[-1]-X[Zero:-1],Te[Zero:-1],'--k',linewidth=2);grid()
                box=ax.get_position();
                ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
                ax.legend(('Bef. Mgmt','After Mgmt','Mvt Thresh.'),loc='center left', bbox_to_anchor=(.95, 0.5))
                ylabel('Bed Stress[N/m^2]',size='large')
                
                ax=subplot(312)
                plot(X[-1]-X[Zero:-1],R1[Zero:-1],X[-1]-X[Zero:-1],R2[Zero:-1],linewidth=2);grid()
                box=ax.get_position();
                ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
                ylabel('Bed Erosn[cm/hr]',size='large')
                
                ax=subplot(313)
                plot(X[-1]-X[Zero:-1],-h[Zero:-1]+S,X[-1]-X[loc1:loc2],-h[loc1:loc2]+S,'xg',linewidth=2);grid()
                box=ax.get_position();
                ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
                ylabel('Depth[m]',size='large')
                xlabel('Cross-Shore Distance [m]',size='large')
                savefig(outputws+"ErosionBed_"+subwsStr+".png",dpi=(640/8))
                Fig2=1 # for HTML
                
        if Oyster==1:
            temp1=5*(X[-1]-Xr)
            figure(2)
            ax=subplot(211)
            plot(X1[-1]-X1[Xloc[0]-temp1:len(h)],H2[Xloc[0]-temp1:len(h)],X1[-1]-X1[Xloc[0]-temp1:len(h)],H1[Xloc[0]-temp1:len(h)],X1[-1]-X1[Xloc[0]-temp1:len(h)],H2[Xloc[0]-temp1:len(h)],'b',linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ylabel('Wave Height[m]',size='large')
            ax.legend(('No Reef','Reef Present'),loc='center left', bbox_to_anchor=(.95, 0.5))            
            ax=subplot(212)
            plot(X[-1]-X[Xloc[0]-temp1:len(h)],-h[Xloc[0]-temp1:len(h)]+S,'k',Xrf,Yrf,'r',Xrf+.05,Yrf,'r',Xrf+.1,Yrf,'r',linewidth=2);grid()
            box=ax.get_position();
            ax.set_position([box.x0, box.y0, box.width*0.8, box.height])                
            ax.legend(('Bed Profile','Oyster Reef'),loc='center left', bbox_to_anchor=(.95, 0.5))
            ylabel('Depth[m]',size='large')
            xlabel('Cross-Shore Distance from Offshore Boundary [m]',size='large')
            savefig(outputws+"ErosionBed_"+subwsStr+".png",dpi=(640/8))
           
    except:
        gp.AddError(msgPlotErosion)
        raise Exception

                
    #################################
    ####### CREATE HTML FILE ########
    #################################

    try:
        gp.AddMessage("...creating html output")
        # create HTML file
        AddMActionHeader = 'yes'
        htmlfile=open(outputws+"OutputWaveModel_"+subwsStr+".html","w")
        htmlfile.write("<html>\n")
        htmlfile.write("<title>Marine InVEST - Wave, Erosion, and Inundation</title>")
        htmlfile.write("<CENTER><H1>Coastal Protection - Tier 1</H1><H2>Nearshore Waves and Erosion Results ("+subwsStr+")<br></H2></CENTER>")

        htmlfile.write("<HR><H2><u>Site Information</u></H2>")
        htmlfile.write("<table border=\"0\" width=\"1100\" cellpadding=\"5\" cellspacing=\"10\">")
        if sand==1:
            htmlfile.write("You have a sandy beach, and the average sediment size is: "+str(d50)+"mm<br><p>")
        elif mud==1:
            htmlfile.write("Your bed is composed of fines/consolidated sediments. <br><p>")
        else:
            htmlfile.write("<p>")
        htmlfile.write("The tidal range at your site is: "+str(round(MHW,1))+"m (high tide value). ")
        if MHW < 2:
            htmlfile.write("It is <i>microtidal</i> (Tidal Range < 2m)<br>")
        elif MHW <= 4:
            htmlfile.write("It is <i>meso-tidal</i> (2 <= Tidal Range <= 4m)<br>")
        else:
            htmlfile.write("It is <i>macro-tidal</i> (Tidal Range > 4m)<br>")
        htmlfile.write("</td></tr></table>")

        if Slope+B1+D1>0 and sand==1:
            htmlfile.write("<HR><H2><u>Backshore Information for Your Sandy Beach System</u></H2>")
            Foreshore=str(int(Slope))
            DuneH=str(round(D1,1))
            BermH=str(round(B1,1))
            BermW=str(round(W1,1))
            htmlfile.write("<li> The foreshore slope is: 1/"+Foreshore+"<br>")
            if mo<>m:
                htmlfile.write("<li>Your input foreshore slope was too flat or your input sediment size was too high for us to estimate the amount of shoreline erosion at your site.  We increased your foreshore slope to 1/" +str(mo) +" to estimate a minimum erosion value.  Feel free to re-adjust your input paramters.")
            
            htmlfile.write("<li>The berm at your site is is: "+BermH+"m high and "+BermW+"m long<br>")
            htmlfile.write("<li>The dune at your site is is: "+DuneH+"m high<br>")
            if Dred>0:
                htmlfile.write("<HR><H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("You will reduce the dune height at your site by: " +str(Dred) +"% <br>")

        if Coral<>0:
            if AddMActionHeader == 'yes':
                htmlfile.write("<HR><H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
            htmlfile.write("<b><i>A coral reef is present at your site.  </i></b>")
            htmlfile.write("<li> Reef-top depth is " +str(hr) +"m, and reef-top width is " +str(Wr))
            if Xco+Xcn==0:
                htmlfile.write("<li> The reef is located at the offshore edge of the profile")
            elif Xco+Xcn==2:
                htmlfile.write("<li> The reef is located at the shoreward edge of the profile")
                
            if CoralMngt=="Dead":
                if AddMActionHeader == 'yes':
                    htmlfile.write("<H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("<li> You assumed that the reef will die, but will remain structurally intact.  We will assume that the coral is smooth.<br>")
            elif CoralMngt=="Rmv":
                if AddMActionHeader == 'yes':
                    htmlfile.write("<H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("<li>You assumed that the reef will die and break.  It will no longer remain structurally intact.  We will assume that the reef top is covered with sand.<br>")
            elif CoralMngt=="None":
                htmlfile.write("<li>You assumed that the reef will not be affected by a particular management action.</u></H2>")

        if MgPst<>0:    
            htmlfile.write("<b><i>A mangrove forest is present at your site.</i></b>")
            if MangMngt=="Half":
                if AddMActionHeader == 'yes':
                    htmlfile.write("<H2><u> Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("<li>You assumed that half of the mangrove will be removed.  We will reduce its density by half.<br>")
            elif MangMngt=="Rmv":
                if AddMActionHeader == 'yes':
                    htmlfile.write("<H2><u> Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("<li>You assumed that the mangrove forest will be fully removed. <br>")
            elif MangMngt=="None":
                htmlfile.write("<li>You assumed that the mangrove forest will not be affected by a particular management action.<br>")
                
        if SgPst<>0:    
            htmlfile.write("<b><i>A seagrass bed is present at your site.</i></b>")
            if SeagMngt=="Half":
                if AddMActionHeader == 'yes':
                    htmlfile.write("<H2><u> Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("<li>You assumed that half of the seagrass bed will be removed.  We will reduce its density by half.<br>")
            elif SeagMngt=="Rmv":
                if AddMActionHeader == 'yes':
                    htmlfile.write("<H2><u> Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("T<li>You assumed that the seagrass bed will be fully removed.<br>")
            elif SeagMngt=="None":
                htmlfile.write("<li>You assumed that the seagrass bed will not affected by a particular management action.<br>")
                
        if MrPst<>0:    
            htmlfile.write("<b><i>A marsh is present at your site.</i></b>")
            if MarshMngt=="Half":
                if AddMActionHeader == 'yes':
                    htmlfile.write("<H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("You assumed that half of the marsh is removed.  We will reduce its density by half.<br>")
            elif MarshMngt=="Rmv":
                if AddMActionHeader == 'yes':
                    htmlfile.write("<H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
                htmlfile.write("You assumed that the marsh will be fully removed. <br>")
            elif MarshMngt=="None":
                htmlfile.write("You assumed that the marsh will not be affected by a particular management action.<br>")
                
        if Oyster==1:
            htmlfile.write("<b><i>An oyster reef is present at your site.</i></b>") # inputs
            if hi>hc:
                htmlfile.write("<tr>The oyster reef is " +str(X[-1]-Xr) +"m from the shoreline, with a base width of " +str(Bw) +"m, and a crest width of "+str(Cw) +"m. It is "+str(hc) +"m tall, and the water depth is " +str(round(hi,2)) +"m: it is submerged.<br><p>")
            else:
                htmlfile.write("<tr>The oyster reef is " +str(X[-1]-Xr) +"m from the shoreline, with a base width of " +str(Bw) +"m, and a crest width of "+str(Cw) +"m. It is"  +str(hc) +"m tall, and the water depth is " +str(round(hi,2)) +"m: it is emergent...please use results with caution.")            
    
        # figures
        htmlfile.write("<HR><H2><u>Model Outputs</u></H2>")
        htmlfile.write("<table border=\"0\" width=\"1100\" cellpadding=\"5\" cellspacing=\"10\"><tr><td>")
        if WaveErosionQuestion=="(2) No, I need to compute these values from wind speed and fetch distance values": 
            htmlfile.write("We computed offshore wave height and period values based on your input of wind speed (U=" +str(Us)+"m/s), fetch distance (Ft="+ str(Ft)+"km), and average water depth at your site(d=" +str(depth) +"m).<br>")

        if Oyster==1:
            htmlfile.write("<b>The oyster reef attenuated " +str(round((1-Kt)*100,1)) +"% of the incident wave height and " +str(round((1-Kt*Kt)*100,1)) +"% of the incident wave energy.</b><br>")
            htmlfile.write('Offshore wave input conditions are: Ho=' +str(round(Ho,2)) +'m, and To=' +str(round(To,2)) +'s<br>')
            
            
            htmlfile.write("The figure below shows close-ups of the profile of wave height and depth profile in the vicinity of the reef.<br>")
            htmlfile.write("<img src=\"ErosionBed_"+subwsStr+".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")          
            
            htmlfile.write("</td><td>The transmission coefficient is Kt= " +str(round(Kt,2)) +", which means that the wave transmitted shoreward of the reef has a height equal to " +str(round(Kt*100,2)) +"% of  the incident wave height right offshore of the reef.</td></tr>")
            
        else:
            htmlfile.write("The figure below shows profiles of wave height along your profile, before and after your management action.<br>")
            htmlfile.write("<img src=\"WavePlot_"+subwsStr+".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")

        if Fig2:
            htmlfile.write("</td><td>The average percent wave attenuation is " +str(round(Atn,2))+"</td></tr>")
            htmlfile.write("<tr><td>The figure below shows profiles of erosion in the intertidal and backshore regions, before and after your management action.<br>")
            htmlfile.write("<img src=\"ErosionBed_"+subwsStr+".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\"></td>")

            # erosion
            if sand==1:
                htmlfile.write("<td><u>Before</u> your management action, your beach might erode by " +str(round(mean(R1),2)) +"m.<br><u>After</u> your management action, your beach might erode by " +str(round(mean(R2),2)) +"m.")
            elif mud==1:
                htmlfile.write("<td><u>Before</u> your management action, your muddy backshore area might experience a maximum average scour rate of " +str(round(mean(R1),2)) +"m.<br><u>After</u> your management action, your muddy backshore area might experience a maximum average scour rate of " +str(round(mean(R2),2)) +"m.")                    
        elif Oyster==1: # oyster reef case
            htmlfile.write("<tr><td>The figure below shows profiles of wave height along your profile in the presence of an oyster reef. ")
            htmlfile.write("<img src=\"WavePlot_"+subwsStr+".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\"></td>")
        
        htmlfile.write("</td></tr></table></html>")
        htmlfile.close()
        
    except:
        gp.AddError(msgHTMLOutputs)
        raise Exception

    # create parameter file
    parameters.append("Script location: "+os.path.dirname(sys.argv[0])+"\\"+os.path.basename(sys.argv[0]))
    parafile=open(outputws+"parameters_"+now.strftime("%Y-%m-%d-%H-%M")+".txt","w") 
    parafile.writelines("WAVE EROSION MODEL PARAMETERS\n")
    parafile.writelines("_____________________________\n\n")
    for para in parameters:
        parafile.writelines(para+"\n")
        parafile.writelines("\n")
    parafile.close()

    del gp

except Exception, ErrorDesc:
    gp.AddMessage(gp.GetMessages())