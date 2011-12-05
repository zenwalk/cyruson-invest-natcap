# Marine InVEST: Coastal Protection (Wave and Erosion)
# Authors: Greg Guannel, Gregg Verutes, Apollo Xi
# 12/02/11

import numpy as num
import string,sys,os,time,datetime,shlex
import fpformat,operator
import arcgisscripting

from win32com.client import Dispatch
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy import optimize
from scipy.integrate import trapz
from math import *
from matplotlib import *
from pylab import *

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
msgArguments="Problem with arguments."

# get parameters
parameters=[]
now=datetime.datetime.now()
parameters.append("Date and Time: "+now.strftime("%Y-%m-%d %H:%M"))
gp.workspace=gp.GetParameterAsText(0)
parameters.append("Workspace: "+ProfileGeneratorWS)
subwsStr=gp.GetParameterAsText(1)
parameters.append("Label for Erosion Run (10 characters max): "+subwsStr)
InputTable=gp.GetParameterAsText(2)
parameters.append("Nearshore Waves and Erosion Excel Table: "+InputTable)
CSProfile=gp.GetParameterAsText(3)
parameters.append("Cross-Shore Profile: "+CSProfile)
WaveErosionQuestion=gp.GetParameterAsText(4)
parameters.append("Do you have wave height and wave period values? "+WaveErosionQuestion)
Ho=gp.GetParameterAsText(5)
parameters.append("IF 1: Wave Height (meters): "+Ho)
To=gp.GetParameterAsText(6)
parameters.append("IF 1: Wave Period (seconds): "+To)
Us=gp.GetParameterAsText(7)
parameters.append("IF 2: Wind Speed (meters per second): "+Us)
Ft=gp.GetParameterAsText(8)
parameters.append("IF 2: Fetch Distance (meters): "+Ft)
depth=gp.GetParameterAsText(9)
parameters.append("IF 2: Water Depth (meters): "+depth)
StormDur=gp.GetParameterAsText(10)
parameters.append("Storm Duration (hours): "+str(StormDur))
S=gp.GetParameterAsText(11)
parameters.append("Surge Elevation (meters): "+str(S))

##if WaveErosionQuestion=="(1) Yes, I have these values"
##if WaveErosionQuestion=="(2) No, I need to compute these values from wind speed and fetch distance values".

##- Ho: Wave height
##- To: Wave period
##- Us: Wind speed
##- Ft: fetch

# remove spaces and shorten 'subwsStr' if greater than 10 characters
subwsStr=subwsStr.replace(" ","")
subwsStr=subwsStr[0:10]

# intermediate and output directories
outputws=ProfileGeneratorWS+os.sep+"_WaveModel_Outputs"+os.sep
interws=ProfileGeneratorWS+os.sep+"scratch"+os.sep

thefolders=["_WaveModel_Outputs"]
for folder in thefolders:
    if not gp.exists(ProfileGeneratorWS+folder):
        gp.CreateFolder_management(ProfileGeneratorWS+os.sep,folder)

# Various Functions and Checks#######################################################################
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

# wave transformation model
def WaveModel(X,h,Ho,T,Seagr,Marsh,Mang,Cf):
    global LagKick
    # constants
    g=9.81;rho=1024.0;B=1.0;Beta=0.05;
    lx=len(X);dx=X[1]-X[0];
    # initialize vegetation vectors,zeros everywhere except where veg present
    VegLoc=X*0.0;
    
    # seagrass
    veg=Seagr[0];loc=Seagr[1];
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
    veg=Mang[0];loc=Mang[1];
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
            
    # drag coefft for vegetation. mangrove and marsh win over seagrass if they overlap
    Cds=num.arange(0.0,lx,dx)*0.0;o=find(VegLoc==1);Cds[o]=0.1 # seagrass	
    Cdr=num.arange(0.0,lx,dx)*0.0;o=find(VegLoc==2);Cdr[o]=0.1 # marsh	
    Cdg=num.arange(0.0,lx,dx)*0.0;o=find(VegLoc==3);Cdg[o]=1 # mangrove

    # initialize vectors for wave model
    H=lx*[0.0];Eta=lx*[0.0];
    Db=lx*[0.0];Df=lx*[0.0]
    k=lx*[0.0];L=lx*[0.0]
    C=lx*[0.0];n=lx*[0.0];Cg=lx*[0.0]
    Ds=lx*[0.0];Dr=lx*[0.0];Dg=lx*[0.0] # seagrass,marsh,mangrove
    Er=lx*[0.0];Ef=lx*[0.0];Br=lx*[0.0]
    Hs=lx*[0.0];Etas=lx*[0.0];
    Dbs=lx*[0.0];Dfs=lx*[0.0]
    Ers=lx*[0.0];Efs=lx*[0.0];Brs=lx*[0.0]
    
    # wave parameter at 1st grid pt
    ash=[h[ii] for ii in range(lx)] # ash is same as h,but is now an independent variable.
    fp=1.0/T; sig=2.0*pi*fp
    k[0]=iterativek(sig,h[0]) # wave number at 1st grid pt
    L[0]=2.0*pi/k[0] # wave length @ 1st grid pt
    n[0]=0.5*(1+(2.0*k[0]*h[0]/sinh(2.0*k[0]*h[0]))) # to compute Cg at 1st grid pt
    C[0]=L[0]/T;Cg[0]=C[0]*n[0] # phase and group velocity at 1st grid pt
    So=Ho/L[0] # deep water wave steepness
    Gam=0.5+0.4*num.tanh(33.0*So) # Gam from Battjes and Stive 85,as per Alsina & Baldock

    # RMS wave height at first grid point; Assume no dissipation occurs
    Ho=Ho/sqrt(2.0) # transform significant wave height to rms wave height
    Co=g*T/(2.0*pi); Cgo=Co/2.0 # deep water phase and group speed
    if h[0]>0.5*L[0]: 
        H[0]=Ho # we are in deep water
    else: 
        H[0]=Ho*sqrt(Cgo/Cg[0]) # we are in intermediate water. Assume no brkg occured

    if LagKick==1: #No shoaling if wave goes from coral reef to lagoon
        H[0]=Ho
        
    if H[0]>0.78*h[0]:
        H[0]=0.78*h[0]
        gp.addmessage('Water depth too shallow for input wave height. That wave broke somewhere in deeper water. We will assume that H=0.78h')

    Ef[0]=0.125*rho*g*H[0]**2*Cg[0];Efs[0]=Ef[0];Hs[0]=H[0] # energy flux @ 1st grid pt

    # begin wave model 
    for xx in range(lx-1): # transform waves,take MWL into account
        #Wave in presence of veg.
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
        temp2=sinh(k[xx+1]*alphs[xx+1]*h[xx+1])**3.0+3.0*sinh(k[xx+1]*alphs[xx+1]*h[xx+1])
        temp3=(3.0*k[xx+1]*cosh(k[xx+1]*h[xx+1])**3.0)
        Dr[xx+1]=temp1*temp2/temp3*H[xx+1]**3.0 # diss due to marshes

        V1=3*sinh(k[xx+1]*alphgr[xx+1]*h[xx+1])+sinh(k[xx+1]*alphgr[xx+1]*h[xx+1])**3 #Roots
        V2=(3*sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1])*h[xx+1])-3*sinh(k[xx+1]*alphgr[xx+1]*h[xx+1])+
            sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1])*h[xx+1])**3-
            sinh(k[xx+1]*alphgr[xx+1]*h[xx+1])**3) #Trunk
        V3=(3*sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1]+alphgc[xx+1])*h[xx+1])
            -3*sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1])*h[xx+1])+
            sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1]+alphgc[xx+1])*h[xx+1])**3-
            sinh(k[xx+1]*(alphgr[xx+1]+alphgt[xx+1])*h[xx+1])**3) #Canopy
        
        CdDN=Cdg[xx+1]*(dgr*Ngr*V1+dgt*Ngt*V2+dgc*Ngc*V3)
        temp1=rho*CdDN/(2.0*sqrt(pi))*(k[xx+1]*g/(2.0*sig))**3.0
        temp2=(3*k[xx+1]*cosh(k[xx+1]*h[xx+1])**3)
        Dg[xx+1]=temp1*temp2*H[xx+1]**3
        
        #Wave in absence of veg
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

    # Estimate MWS
    val=1.0;vals=1.0;dell=1.0;dells=1.0;Os=0.0;O=0.0 # for setup calc
    Sxx=lx*[0.0]; Rxx=lx*[0.0]; Sxxs=lx*[0.0]; Rxxs=lx*[0.0]
    
    # force on plants if they were emergent; take a portion if plants occupy only portion of wc
    Fxs=[rho*g*Cds[i]*ds*Ns*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
    fxs=[-alphs[i]*Fxs[i] for i in range(lx)] #Seagrass
    
    Fxr=[rho*g*Cdr[i]*dr*Nr*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
    fxr=[-alphr[i]*Fxr[i] for i in range(lx)] #Marsh

    Fxgr=[rho*g*Cdg[i]*dgr*Ngr*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
    Fxgt=[rho*g*Cdg[i]*dgt*Ngt*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
    Fxgc=[rho*g*Cdg[i]*dgc*Ngc*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
    fxg=[-alphgr[i]*Fxgr[i]-alphgt[i]*Fxgt[i]-alphgc[i]*Fxgc[i] for i in range(lx)] #Mangrove

    fx=[fxs[ii]+fxr[ii]+fxg[ii] for ii in range(lx)] #All
    
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

    Ubot=[pi*H[ii]/(T*sinh(k[ii]*h[ii])) for ii in range(lx)] #Bottom velocity
    H=[H[ii]*sqrt(2) for ii in range(lx)]
    Hs=[Hs[ii]*sqrt(2) for ii in range(lx)]

    return H,Eta,Hs,Etas,Ubot,VegLoc

def WaveModelSimple(X,h,Ho,T):
    global LagKick
    # constants
    g=9.81;rho=1024.0;B=1.0;Beta=0.05;
    lx=len(X);dx=X[1]-X[0];
    
    # initialize vectors for wave model
    H=lx*[0.0];Eta=lx*[0.0];L=lx*[0.0]
    Db=lx*[0.0];Df=lx*[0.0]
    Er=lx*[0.0];Ef=lx*[0.0];Br=lx*[0.0]
    C=lx*[0.0];n=lx*[0.0];Cg=lx*[0.0]
    k=lx*[0.0];
    
    # wave parameter at 1st grid pt
    ash=[h[ii] for ii in range(lx)] # ash is same as h,but is now an independent variable.
    fp=1.0/T; sig=2.0*pi*fp
    k[0]=iterativek(sig,h[0]) # wave number at 1st grid pt
    L[0]=2.0*pi/k[0] # wave length @ 1st grid pt
    n[0]=0.5*(1+(2.0*k[0]*h[0]/sinh(2.0*k[0]*h[0]))) # to compute Cg at 1st grid pt
    C[0]=L[0]/T;Cg[0]=C[0]*n[0] # phase and group velocity at 1st grid pt

    # RMS wave height at first grid point; Assume no dissipation occurs
    Ho=Ho/sqrt(2.0) # transform significant wave height to rms wave height
    Co=g*T/(2.0*pi); Cgo=Co/2.0 # deep water phase and group speed
    if h[0]>0.5*L[0]: 
        H[0]=Ho # we are in deep water
    else: 
        H[0]=Ho*sqrt(Cgo/Cg[0]) # we are in intermediate water. Assume no brkg occured

    if LagKick==1: #No shoaling if wave goes from coral reef to lagoon
        H[0]=Ho
        
    if H[0]>0.78*h[0]:
        H[0]=0.78*h[0]
        gp.addmessage('Water depth too shallow for input wave height. That wave broke somewhere in deeper water. We will assume that H=0.78h')

    Ef[0]=0.125*rho*g*H[0]**2*Cg[0] # energy flux @ 1st grid pt
    So=Ho/L[0] # deep water wave steepness
    Gam=0.5+0.4*num.tanh(33.0*So) # Gam from Battjes and Stive 85,as per Alsina & Baldock
    
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
    BrkRim=0;BrkFace=0;Kp=0.0;

    # reef profile
    Xreef=num.arange(0.0,10000.0,dx)
    Yreef=AlphR*Xreef-he
    loc=find(Yreef>-hr);
    Yreef[loc]=-hr # reef profile

    D=T*sqrt(g/hr) # relative subm
    Lo=g*T**2.0/(2.0*pi)
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

            xs=(2+1.1*Ho/he)*T*(g*he)**.5 # surf zone width
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

            xs=(2+1.1*Ho/he)*T*(g*he)**.5 # surf zone width
            loc1=Indexed(Yreef,-db) # start brkg
            loc2=Indexed(Xreef,Xreef[loc1]+xs) # end brkg
            ha=-average(Yreef[loc1:loc2+1]) # rep. depth

    # wave transformation on top of coral reef
    Etar=0.0;delta=10;
    Xflat=num.arange(0.0,Wr,dx)
    lx=len(Xflat);
    H_r=lx*[Ho];

    if Kp<> 0: #Wave break on reef face or rim
        # wave Setup
        while delta>0.01:
            EtaR=3.0/(64.0*pi)*Kp*Ho**2*T*g**.5/((Etar+ha)**1.5)
            delta=abs(EtaR-Etar);Etar=EtaR
        Etar=Etar[0]    
        H_r[0]=0.42*(hr+Etar) #Wave height at the offshore edge of reef
    else:
        H_r[0]=Ho
        
    for xx in range(lx-1):
        H_r[xx+1]=H_r[xx]-Cf*H_r[xx]**2.0/(3.0*pi*(Etar+hr)**2.0)*dx
        
    return H_r,Etar

# wave attenuation by reef breakwater
def BreakwaterKt(Hi,T,hi,hc,Cwidth,Bwidth):
    Lo=9.81*T**2.0/(2.0*pi)
    Rc=hc-hi # depth of submergence
    Boff=(Bwidth-Cwidth)/2.0 # base dif on each side
    ksi=(hc/Boff)/sqrt(Hi/Lo)
    
    if Cwidth<>0: #It's not a reef ball
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
            
    else: #It's a reef ball
        Kt=1.616-31.322*Hi/(9.81*T**2)-1.099*hc/hi+0.265*hi/Bwidth
        
    return Kt

# K&D Erosion model
def ErosionKD(A,Ho,RunupVal,B,D,W,xb,hb):
    # constants
    BD=D+B;TWL=S+RunupVal
    if TWL>B0+D0:
        TWL=B0+D0

    # bkg info 
    if xb+hb==0:
        C0=g*T/(2*pi) # deep water phase speed
        hb=(((Ho**2.0)*C0)/2)**(2.0/5.0)/(g**(1.0/5.0)*0.73**(4.0/5.0));Hb=0.73*hb # breaking depth
        xb=(hb/A)**1.5 # surf zone width

    # erosion model
    Term1=xb-hb/m
    if Term1<0: # profile factors wrong; can't compute erosion
        gp.AddMessage("The waves reaching the beach don't create erosion, foreshore slope too flat or sediment size is too high.")
        Rinf=0 # can't compute
        R0=0
        
    else: # compute erosion
        Rinf=(TWL*Term1)/(BD+hb-TWL/2)-W*(B+hb+0.5*TWL)/(BD+hb-0.5*TWL) # potential erosion distance
        if Rinf<0:
            gp.AddMessage("Berm is so wide that dune is not eroding.")
            Rinf=(TWL*Term1)/(B+hb-TWL/2)
            
        TS=(320.0*(Hb**(3.0/2.0)/(g**.5*A**3.0))*(1.0/(1.0+hb/BD+(m*xb)/hb)))/3600.0 # erosion response time scale
        global BetaKD
        BetaKD=2.0*pi*(TS/StormDur)

        expr="num.exp(-2.0*x/BetaKD)-num.cos(2.0*x)+(1.0/BetaKD)*num.sin(2.0*x)" # solve this numerically
        fn=eval("lambda x: "+expr)
        z=FindRootKD(fn,pi,pi/2) # find zero in function,initial guess from K&D

        R0=0.5*Rinf*(1.0-cos(2.0*z)) # final erosion distance

    return Rinf,R0

# Reading Depth Profile and Excel Inputs############################################################
gp.AddMessage("\nReading depth profile and Excel inputs...")
# constants
g=9.81;rho=1024.0
Fig2=0; # for HTML

# read in user's cross-shore profile
TextData=open(CSProfile,"r") 
X_lst=[];h_lst=[]
for line in TextData.read().strip("\n").split("\n"):
    linelist=[float(s) for s in line.split("\t")] # split the list by comma delimiter
    X_lst.append(linelist[0])
    h_lst.append(linelist[1])
TextData.close()
X=num.array(X_lst);h=num.array(h_lst);
if h[0]>0 or h[0]>h[-1]:    h=h[::-1] # reverse order if profile starts at shoreline

# input from user via Excel SS
xlApp=Dispatch("Excel.Application")
xlApp.Visible=0
xlApp.DisplayAlerts=0
xlApp.Workbooks.Open(InputTable)
cell1=xlApp.Worksheets("WaveModelInput")
cell3=xlApp.Worksheets("ReefShapeFactor")

#Read Beach Information
temp=cell1.Range("c14:h14").value;temp=temp[0] #All beach information
d50=temp[0] # sediment size
Slope=temp[1];
if Slope<>0:
    m=1.0/Slope # foreshore slope=1/Slope# foreshore slope
else:
    m=0
A=cell1.Range("e51").value; #Sediment scale factor

B0=temp[2];
W0=temp[3]
D0=temp[4]
Dred=temp[5]

#Read Muddy shoreline information
Cm=cell1.Range("c18").Value #dry density
me=cell1.Range("d18").Value #Erosion constant

#Read Oyster Reef Information
temp=cell1.Range("c22:f22").value;temp=temp[0] #All oyster information
Xr=X[-1]-temp[0] # distance from shoreline; need to reverse because user enter with shoreline at X=0
hc=temp[1] # reef height
Bw=temp[2] # base width
Cw=temp[3] # crest width

#Read Vegetation Information
temp=cell1.Range("e27:j31").Value #All vegetation information

temp1=temp[2] #Mangrove canopy 
hogc=temp1[0] # height of mangrove canopy
dogc=temp1[1] # diameter of mangrove canopy
Nogc=temp1[2] # density of mangrove canopy
temp1=temp[1] #Mangrove trunks
hogt=temp1[0] # height of mangrove trunks
dogt=temp1[1] # diameter of mangrove trunks
Nogt=temp1[2] # density of mangrove trunks
temp1=temp[0] #Mangrove roots
hogr=temp1[0] # height of mangrove roots
dogr=temp1[1] # diameter of mangrove roots
Nogr=temp1[2] # density of mangrove roots
Xo1g=X[-1]-temp1[4] # offshore edge of mangrove
Xo2g=X[-1]-temp1[3] # shoreward edge of mangrove
MangMngt=temp1[5] #Management action for mangrove

temp1=temp[3]
hos=temp1[0] # height of seagrass
dos=temp1[1] # diameter of seagrass
Nos=temp1[2] # density of seagrass
Xo1s=X[-1]-temp1[4] # offshore edge of seagrass
Xo2s=X[-1]-temp1[3] # shoreward edge of seagrass
SeagMngt=temp1[5] #Management action for seagrass

temp1=temp[4]
hor=temp1[0] # height of marsh
dor=temp1[1] # diameter of marsh
Nor=temp1[2] # density of marsh
Xo1r=X[-1]-temp1[4] # offshore edge of marsh
Xo2r=X[-1]-temp1[3] # shoreward edge of marsh
MarshMngt=temp1[5] #Management action for seagrass

#Read Coral Information
temp=cell1.Range("c35:j35").Value;temp=temp[0] #All coral information
Xco=temp[1] #off. distance
Xcn=temp[0] #nearshore distance
if (Xco+Xcn)>0.0:
    Xco=X[-1]-Xco;Xcn=X[-1]-Xcn
    
AlphF=temp[2] # reef face slope
AlphR=temp[3]    # reef rim slope
he=temp[4] # reef rim edge
hr=temp[5] # reef top depth
Wr=temp[6] # reef top width
CoralMngt=temp[7]
TanAlph=cell3.Range("a2:a202").Value;TanAlph=num.array(TanAlph)
kp=cell3.Range("b2:b202").Value;kp=num.array(kp) # reef shape factor
kp2=cell3.Range("c2:c202").Value;kp2=num.array(kp2) # reef shape factor

xlApp.ActiveWorkbook.Close(SaveChanges=0)
xlApp.Quit()

# Reading Depth Profile and Creating Inputs for Run#################################################
gp.AddMessage("\nPreparing inputs...")

# modify the depth profile appropriately
if hor+hogr+hogt+hogc<> 0: # if there's a marsh or a mangrove
    h=h+S #Add surge level
    out=find(h<-0.2);
    h=h[out];X=X[out];# only keep values that are below water
    dx=X[1]-X[0];m=abs(h[-1]-h[-int(10.0/dx)])/10 # average slope 10m from end of transect
else: #it's a beach
    out=find(h<-0.2);
    h=h[out];X=X[out]; # only keep values that are below water

#Resample bathy to 0.5m
dx=1;F=interp1d(X,h);X=num.arange(X[0],X[-1]+dx,dx);
h=F(X);lx=len(X) #Resample to dx

#Management Actions#
#Beach
D1=(100-Dred)*D0/100 #New dune height
if A>0 and B0+W0+D0<>0:
    Beach=1 #A sandy beach is present
    gp.AddMessage("A sandy beach is present in your system")
else:
    Beach=0
    gp.AddMessage("You do not have a sandy beach in your system. Model will not compute erosion for the type of sediment that you have.")

    #Prepare Oyster data
if hc+Bw<> 0:
    Oyster=1 #Oyster information was entered
    gp.AddMessage("An oyster reef is present in your system")
else:
    Oyster=0
    gp.AddMessage("You do not have an oyster reef in your system")

    #Prepare Mangrove data
Mang0=[[[0,0,0],[0,0,0],[0,0,0]],[0,0]];Mang=[[[0,0,0],[0,0,0],[0,0,0]],[0,0]]
Mang1=[[[0,0,0],[0,0,0],[0,0,0]],[0,0]]
loc=[Xo1g,Xo2g]
if loc is not None and (hogc+dogc+Nogc+hogt+dogt+Nogt+hogr+dogr+Nogr)<> 0:
    veg=[[hogr,dogr,Nogr],[hogt,dogt,Nogt],[hogc,dogc,Nogc]];
    Mang=[veg,loc] # mangrove characteristics
    
    if MangMngt=="None": #Management action
        Mang1=Mang
    elif MangMngt=="Rmv":
        Mang1=Mang0
    elif MangMngt=="Half":
        veg=[[hogr,dogr,Nogr/2.0],[hogt,dogt,Nogt/2.0],[hogc,dogc,Nogc]];#Rdce dens.roots&trees by 1/2
        Mang1=[veg,loc]
        
MgPst=sum(sum(Mang[0]))+sum(sum(Mang1[0]));
if MgPst<>0:  
    MgPst=1
    gp.AddMessage("A mangrove is present in your system")
else:
    gp.AddMessage("You do not have a mangrove in your system")
    
    #Prepare Seagrass data
Seagr0=[[0,0,0],[0,0]];Seagr=[[0,0,0],[0,0]];Seagr1=[[0,0,0],[0,0]]
loc=[Xo1s,Xo2s]
if loc is not None and (hos+dos+Nos)<> 0:
    veg=[hos,dos,Nos];
    Seagr=[veg,loc] # seagrass characteristics
    
    if SeagMngt=="None": #Management action
        Seagr1=Seagr
    elif SeagMngt=="Rmv":
        Seagr1=Seagr0
    elif SeagMngt=="Half":
        veg=[hos,dos,Nos/2];
        Seagr1=[veg,loc] 
        
SgPst=sum(Seagr)+sum(Seagr1)
if SgPst<>0:  
    SgPst=1
    gp.AddMessage("A seagrass bed is present in your system.")
else:
    gp.AddMessage("You do not have a seagrass bed in your system.")

    #Prepare marsh data
Marsh=[[0,0,0],[0,0]];Marsh0=[[0,0,0],[0,0]];Marsh1=[[0,0,0],[0,0]]
loc=[Xo1r,Xo2r]
if loc is not None and (hor+dor+Nor)<> 0:
    veg=[hor,dor,Nor];
    Marsh=[veg,loc] # marsh characteristics

    if MarshMngt=="None": #Management action
        Marsh1=Marsh
    elif MarshMngt=="Rmv":
        Marsh1=Marsh0
    elif MarshMngt=="Half":
        veg=[hor,dor,Nor/2.0];
        Marsh1=[veg,loc] 
        
MrPst=sum(sum(Marsh)+sum(Marsh1))
if MrPst<>0:    
    MrPst=1 
    gp.AddMessage("A marsh is present in your system")
else:
    gp.AddMessage("You do not have a marsh in your system")

    #Prepare Coral data
if AlphF+AlphR+he+hr+Wr<>0: #Check if coral profile is imported o
    Coral=1
    gp.AddMessage("A coral reef is present in your system.  We'll estimate its profile for you")
if Xcn+Xco<>0 and AlphF+AlphR+he+hr+Wr==0: #Check if coral profile is part of the bathy
    Coral=-1
    gp.AddMessage("A coral reef is present in your system.  It is incorporated in the bathy profile you uploaded")
elif Xcn+Xco+AlphF+AlphR+he+hr+Wr==0:    
    Coral=0
    gp.AddMessage("You do not have a coral reef in your system")

Cf_bed=num.arange(0,len(h),.5)*0+0.01 #Bottom friction for sand only
Cf1_bed=num.arange(0,len(h),.5)*0+0.01  #for management action
Cf =0.2;Cf1=0.2; #Friction for live and dead coral

temp1=find(X>=Xco);temp2=find(X>Xcn);
Loco=[ii for ii in range(temp1[0],temp2[0],1)] #location of coral reef in profile
if len(Loco)>2:   Cf_bed[Loco]=0.2; #Live coral

if Coral==1:
    if CoralMngt=="Rmv":        Cf1=.01;#No coral
    elif CoralMngt=="Dead":        Cf1=.1;#Dead coral
    elif CoralMngt=="None":        Cf1=.2 #Live coral
elif Coral==-1 and len(Loco)>2:
    if CoralMngt=="Rmv":        Cf1_bed[Loco]=.01;#No coral
    elif CoralMngt=="Dead":        Cf1_bed[Loco]=.1;#Dead coral
    elif CoralMngt=="None":        Cf1_bed[Loco]=0.2 #Live coral
            
# Run model for management action#################################################################
#BEFORE MANAGEMENT ACTION
gp.AddMessage("\nComputing wave height profiles before management action...")

h=-h-2 #Depth is now positive
Ho1=H0;#Offshore wave height
LagKick=0 #check for lagoon calc
H_r=[];Xn=[];ho=[];Hsimple1=[];Etasimple1=[]

    # Estimate initial wave height if coral reef at beg of transect and no bathy for that reef
if Xco==0.0 and Xcn==0.0 and Coral==1: # reef top depth
    H_r,Eta_r=WavesCoral(H0,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf,dx);
    H_r=num.array(H_r);Eta_r=H_r*0+Eta_r #Arrays of H_r and Eta_r
    Ho1=H_r[-1];ho=h[0]; #Wave height at the end of the reef
    Xn=num.arange(-Wr,0,dx);    LagKick=1; #For future
    

    #Estimate wave height transformation over bathy profile
if len(h)>2: 
    H,Eta,Hs,Etas,Ubot,VegLoc1=WaveModel(X,h,Ho1,T,Seagr,Marsh,Mang,Cf_bed)

if Xo1s>Xcn and (MgPst+MrPst+SgPst)<>0: #Wave height in lagoon in absence of any vegetation
    Hsimple1=Hs;Etasimple1=Etas        

    #Wave height over coral reefs at various locations
if len(h)>2 and Xco>=0 and Xcn>0 and Coral==1: #Coral Reef in middle of transect
    Ho1=H[Loco[0]] #Wave height at edge of the coral reef
    H_r,Eta_r=WavesCoral(Ho1,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf,dx);
    H_r=num.array(H_r);Eta_r=H_r*0+Eta_r #Arrays of H_r and Eta_r
    LagKick=1;ho=(h[Loco[0]]+h[Loco[-1]])/2;
    Xn=num.arange(0,Wr,dx);
    
    temp=num.array((int(X[Loco[-1]]-X[Loco[0]]-Wr)*2)*[H[Loco[0]]]) #Constant H over reef face/rim
    H_r=num.append(temp,H_r,None) #Wave height profile over reef
    H1=num.append(num.array(H[0:Loco[0]]),H_r,None); #All wave profile to edge of reef
    X_r=num.arange(X[Loco[0]],X[Loco[0]]+.5*len(H_r),dx)#X-axis for wave over reef
    X1=num.append(num.array(X[0:Loco[0]]),X_r,None); #All X profile to edge of reef

    temp=num.array((int(X[Loco[-1]]-X[Loco[0]]-Wr)*2)*[Eta[Loco[0]]])#Constant Eta over reef face/rim
    Eta_r=num.append(temp,Eta_r,None) #MWL profile over reef
    Eta1=num.append(num.array(Eta[0:Loco[0]]),Eta_r,None); #MWL to edge of reef

    Ho1=H_r[-1] #Wave height at the end of the reef
    temp=X[Loco[-1]:-1],h[Loco[-1]:-1] #Lagoon X
    Hlag,Etalag,Hs,Etas,Ubot,temp=WaveModel(temp,Ho1,T,Seagr,Marsh,Mang,Cf_bed[Loco[-1]:-1])#H in lagoon
    
    if Xo1s>Xcn and (MgPst+MrPst+SgPst)<>0: #Wave height in lagoon in case there is vegetation
        Hsimple1=Hs;Etasimple1=Etas        

elif len(h)>2 and Xco==1 and Xcn==1 and Coral==1: #Coral Reef end of transect
    H_r,Eta_r=WavesCoral(H[-1],AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf,dx)
    H_r=num.array(H_r);Eta_r=H_r*0+Eta_r #Arrays of H_r and Eta_r
    Ho1=H_r[0];ho=h[-1]; #Wave height at offshore edge of reef
    Xn=num.arange(X[-1],X[-1]+Wr,dx)
        
if Oyster==1:
    Xloc=find(X>=(X[-1]-Xr));#Location of oyster reef
    Ho1=H[Xloc[0]]; LagKick=1#Wave height at edge of the reef
    hi=h[Xloc[0]] #Water depth at reef location
    Kt=BreakwaterKt(Ho1,T,hi,hc,Cw,Bw);Ho1=Kt*Ho1 #Transmitted Wave height
      
    gp.AddMessage("\nComputing wave height profiles in presence of oyster reef...")
    LagKick=1;
    Hlag,Etalag,Hs,Etas,Ubot,temp=WaveModel(X[Xloc],h[Xloc],Ho1,T,Seagr,Marsh,Mang,Cf_bed[Xloc]) 

    if Xo1s>Xr and (MgPst+MrPst+SgPst)<>0: #Wave height in lagoon in case there is vegetation
        Hsimple1=Hs;Etasimple1=Etas        

    #Put all wave charact. together
if Xco==0 and Xcn ==0 and Coral==1: #Reef at offshore edge of the profile
    X1=num.append(Xn,X,None);
    H1=num.append(H_r,num.array(H),None)
    Eta1=num.append(Eta_r,num.array(Eta),None)
elif len(h)>2 and Xco>=0 and Xcn>0 and Coral==1: #Reef somewhere in profile
    X1=num.append(X1,X[Loco[-1]:-1],None);
    H1=num.append(num.array(H1),num.array(Hlag),None)
    Eta1=num.append(num.array(Eta1),num.array(Etalag),None)
elif len(h)>2 and Xco==1 and Xcn==1 and Coral==1: #Reef at the nearshore edge of the profile
    X1=num.append(X,Xn,None);
    H1=num.append(num.array(H),H_r,None)
    Eta1=num.append(num.array(Eta),Eta_r,None)
elif Oyster==1:
    X1=num.array(X);
    H1=num.append(H[0:Xloc[0]],Hlag);
    Eta1=num.append(Eta[0:Xloc[0]],Etalag);
else:
    X1=num.array(X);H1=num.array(H);Eta1=num.array(Eta)

#AFTER MANAGEMENT ACTION
gp.AddMessage("\nComputing wave height profiles after management action...")

Ho2=H0; #Offshore wave height
LagKick=0 #check for lagoon calc

if SgPst+MrPst+MgPst+Coral<>0: #If management action affects vegetation or coral
    Hsimple2=[];Etasimple2=[];
    
        # Estimate initial wave height if coral reef at beg of transect and no bathy for that reef
    if Xco==0 and Xcn==0.0 and Coral==1: # reef top depth
        H_r,Eta_r=WavesCoral(Ho2,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf1,dx)    
        H_r=num.array(H_r);Eta_r=H_r*0+Eta_r #Arrays of H_r and Eta_r
        Ho2=H_r[-1] #Wave height at the end of the reef
        Xn=num.arange(-Wr,0,dx);    LagKick=1;
    
        #Estimate wave height transformation over bathy profile
    if len(h)>2 and CoralMngt+MarshMngt+SeagMngt+MangMngt<>"NoneNoneNoneNone": #H over bathy
        H,Eta,Hs,Etas,Ubot,VegLoc2=WaveModel(X,h,Ho2,T,Seagr1,Marsh1,Mang1,Cf1_bed);
        
    if Xo1s>Xcn and (MgPst+MrPst+SgPst)<>0: #Estimate wave height in lagoon in case there is seagrass
        Hsimple2=Hs;Etasimple2=Etas
        
    if len(h)>2 and Xco>=0 and Xcn>0 and Coral==1: #Coral Reef in middle of transect
        Ho2=H[Loco[0]] #Wave height at edge of the coral reef
        H_r,Eta_r=WavesCoral(Ho2,AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf1,dx);
        H_r=num.array(H_r);Eta_r=H_r*0+Eta_r #Arrays of H_r and Eta_r
    
        Xn=num.arange(0,Wr,dx);LagKick=1;
        temp=num.array((int(X[Loco[-1]]-X[Loco[0]]-Wr)*2)*[H[Loco[0]]]) #H over reef face/rim
        H_r=num.append(temp,H_r,None) #Wave height profile over reef
        temp=num.array((int(X[Loco[-1]]-X[Loco[0]]-Wr)*2)*[Eta[Loco[0]]])
        Eta_r=num.append(temp,Eta_r,None) #MWL profile over reef
        X_r=num.arange(X[Loco[0]],X[Loco[0]]+.5*len(H_r),dx)#X-axis for wave over reef
            
        H2=num.append(num.array(H[0:Loco[0]]),H_r,None); #All wave profile to edge of reef
        Eta2=num.append(num.array(Eta[0:Loco[0]]),Eta_r,None); #All wave profile to edge of reef
        X2=num.append(num.array(X[0:Loco[0]]),X_r,None); #All X profile to edge of reef
    
        Ho2=H_r[-1] #Wave height at the end of the reef
        temp=X[Loco[-1]:-1],h[Loco[-1]:-1]
        Hlag,Etalag,Hs,Etas,Ubot,temp=WaveModel(temp,Ho2,T,Seagr1,Marsh1,Mang1,Cf1_bed[Loco[-1]:-1]) 
        
        if MXo1s>Xcn and (MgPst+MrPst+SgPst)<>0: #Wave height in lagoon in case there is vegetation
            Hsimple2=Hs;Etasimple2=Etas
    
    elif len(h)>2 and Xco==1 and Xcn==1 and Coral==1: #Coral Reef at nearshore edge
        H_r,Eta_r=WavesCoral(H[-1],AlphF ,AlphR,he,hr,Wr,Xco,Xcn,Cf1,dx)
        H_r=num.array(H_r);Eta_r=H_r*0+Eta_r #Arrays of H_r and Eta_r
        Ho2=H_r[0] #Wave height at offshore edge of reef
        Xn=num.arange(X[-1],X[-1]+Wr,dx);
     
        #Oyster Reef
    if Oyster==1:
        Xloc=find(X>=(X[-1]-Xr));#Location of oyster reef
        Ho2=H[Xloc[0]]; LagKick=1#Wave height at edge of the reef
        hi=h[Xloc[0]] #Water depth at reef location
        Kt=BreakwaterKt(Ho2,T,hi,hc,Cw,Bw);Ho2=Kt*Ho2 #Transmitted Wave height
          
        gp.AddMessage("\nComputing wave height profiles in presence of oyster reef...")
        Hlag,Etalag,Hs,Etas,Ubot,temp=WaveModel(X[Xloc],h[Xloc],Ho2,T,Seagr1,Marsh1,Mang1,Cf1_bed[Xloc]) 
        
        if Xo1s>Xr and (MgPst+MrPst+SgPst)<>0: #Wave height in lagoon in case there is vegetation
            Hsimple2=Hs;Etasimple2=Etas
    
        #Put all wave charact. together
    if Xco==0 and Xcn ==0 and Coral==1: #Reef at offshore edge of the profile
        X2=num.append(Xn,X,None);
        H2=num.append(H_r,num.array(H),None)
        Eta2=num.append(Eta_r,num.array(Eta),None)
    elif len(h)>2 and Xco==1 and Xcn==1 and Coral==1: #Reef at the nearshore edge of the profile
        X2=num.append(X,Xn,None);
        H2=num.append(num.array(H),H_r,None)
        Eta2=num.append(num.array(Eta),Eta_r,None)
    elif len(h)>2 and Xco>=0 and Xcn>0 and Coral==1: #Reef somewhere in profile
        X2=num.append(X2,X[Loco[-1]:-1],None);
        H2=num.append(num.array(H2),num.array(Hlag),None)
        Eta2=num.append(num.array(Eta2),num.array(Etalag),None)
    elif Oyster==1:
        X2=num.array(X);
        H2=num.append(H[0:Xloc[0]],Hlag);
        Eta2=num.append(Eta[0:Xloc[0]],Etalag);
    else:
        X2=num.array(X);H2=num.array(H);Eta2=num.array(Eta)
else:
    Hsimple2=num.array(Hsimple1);Etasimple2=num.array(Etasimple1);
    H2=num.array(H1);X2=num.array(X1);Eta2=num.array(Eta1);
        
#CREATE OUTPUTS
# percent wave attenuation
lx=len(X1)
AtnH=H1*0
for xx in range(lx):
    AtnH[xx]=abs(H2[xx]-H1[xx])/H1[xx]*100 # percent attenuation
Atn=find(AtnH>0.5); #Take Attn> 0.5% to remove any num. diff.
if len(Atn)>0:    Atn=trapz(AtnH[Atn],X1[Atn],dx)/(X1[Atn[-1]]-X1[Atn[0]])
else:    Atn=0
    
# save wave outputs
WaveHeight=outputws+"WaveHeight_"+subwsStr+".txt"
file=open(WaveHeight,"w")
for i in range(0,len(X1)):
    file.writelines(str(X1[i])+"\t"+str(H1[i])+"\t"+str(H2[i])+"\n")
file.close()

#indices of vegetation location
keep=find(VegLoc1<>0);
VegLoc=VegLoc1+nan;
if len(keep)>1:  VegLoc[keep]=-h[keep];

# plot
figure(1);
subplot(311);plot(X1,H1,X2,H2);grid()
ylabel('Wave Height [m]',size='large')
legend(('Initial Condition','Management Action'),'upper left')
subplot(3,1,2);plot(X1,AtnH);grid()
ylabel('Wave Attenuation [%]',size='large')
if sum(VegLoc1)>0 and Coral<>0:
    if Xco+Xcn>3:
        subplot(313);plot(X[Loco],X[Loco]*0.0-ho,'ob',X,VegLoc,'xg',X,-h);grid()
    elif Xco+Xcn==0:
        temp=num.arange(-Wr,0,1);
        subplot(313);plot(temp,temp*0.0-ho,'ob',X,VegLoc,'xg',X,-h);grid()
    elif Xco+Xcn==2:
        temp=num.arange(X[-1],X[-1]+Wr,1);
        subplot(313);plot(temp,temp*0.0-ho,'ob',X,VegLoc,'xg',X,-h);grid()
        
    legend(('Coral Reef','Vegetation Field','Depth Profile'),'upper left')
elif sum(VegLoc1)>0 and Coral==0:
    subplot(313);plot(X,VegLoc,'xg',X,-h);grid()
    legend(('Vegetation Field','Depth Profile'),'upper left')
elif sum(VegLoc1)==0 and Coral<>0:
    subplot(313);plot(X[Loco[0]]+Xn,num.array(Xn)*0.0-ho,'ob',X,-h);grid()
    legend(('Coral Reef','Depth Profile'),'upper left')
ylabel('Water Depth [m]',size='large')
xlabel('Cross-Shore Distance from Offshore Boundary')
savefig(outputws+"WavePlot_"+subwsStr+".png",dpi=(640/8))

## ESTIMATE RUNUP AMOUNT
Lo=g*T**2.0/(2.0*pi)

#Before Management Action
Rmax1=1.1*(0.35*m*sqrt(Ho1*Lo)+sqrt(Lo*(Ho1*0.563*m**2.0+H0*0.004))/2.0) # max runup--short and long waves are different (Stockdon adaptation)
Emax1=0.35*m*num.sqrt(Ho1*Lo);

Etap1=Emax1;Hp1=Ho1;Rnp1=Rmax1; 
if len(Etasimple1)<>0 and (MgPst+MrPst+SgPst)<>0: #If there is vegetation
    coef0=max(Etasimple1[-200:-1])/Emax1;Etap1=max(Eta1[-200:-1])/coef0 #Corr. MWL at shoreline; take last 200m
    if Etap1<0:        Etap1=0 # if Eta with veg. neg,take as zero
    Hp1=(Etap1/(0.35*m))**2/Lo # Hprime to estimate runup with veg
    Rnp1=1.1*(Etap1+num.sqrt(Lo*(Hp1*0.563*m**2+0.004*H0))/2) # runup with vegetation
    
#After manamgement action
Rmax2=1.1*(0.35*m*sqrt(Ho2*Lo)+sqrt(Lo*(Ho2*0.563*m**2.0+H0*0.004))/2.0) #Rnp--short and long waves 
Emax2=0.35*m*num.sqrt(Ho2*Lo);

Etap2=Emax2;Hp2=Ho2;Rnp2=Rmax2;
if len(Etasimple2)<>0 and (MgPst+MrPst+SgPst)<>0: #If there is vegetation
    coef0=max(Etasimple2[-200:-1])/Emax2;Etap2=max(Eta2[-200:-1])/coef0 #Corr. MWL at shoreline; take last 200m
    if Etap2<0:        Etap2=0 # if Eta with veg. neg,take as zero
    Hp2=(Etap2/(0.35*m))**2/Lo # Hprime to estimate runup with veg
    Rnp2=1.1*(Etap2+num.sqrt(Lo*(Hp2*0.563*m**2+0.004*H0))/2) # runup with vegetation
    
#ESTIMATE EROSION AMOUNT
if len(Etasimple1)<>0 and Beach==1:
    gp.AddMessage('Estimate erosion amount for sandy beach ...')
    
    #Estimate eq. profile
    temp=argmin(abs(h-20)); #Assume 20 is closure depth
    y=h[temp:-1];y=y[::-1];x=X[temp:-1];x=x-x[0]; #Profile used to estimate sed. scale factor
    fitfunc=lambda p, ix: p[0]*(ix)**(2.0/3) # Target function
    errfunc=lambda p, ix, iy: fitfunc(p, ix) - iy # Distance to the target function
    p0=[0.1] # Initial guess for the parameters
    p1, success=optimize.leastsq(errfunc, p0[:], args=(x,-y))
    A1=abs(p1[0]) #Sediment scale factor
    y=A1*x**(2.0/3);y=y[::-1]
    keep=find(y>0.5);
    y=y[keep];x=x[keep]    
    
        # Before management action
    #Wave height over eq. profile
    Hs,Etas=WaveModelSimple(x,y,Ho1,T) 
    temp=argmin(Etas)
    xb=x[-1]-x[temp];hb=y[temp] #Breaking wave height     
    Rinf,temp1=ErosionKD(A,Ho1,Rmax1,B0,D0,W0,0,0) #Erosion of beach - no segrass (in lagoon)
    Rinf,temp2=ErosionKD(A1,Ho1,Rmax1,B0,D0,W0,0,0) #Erosion of beach - no segrass (in lagoon)
    Rinf,temp=ErosionKD(A1,Ho1,Rmax1,B0,D0,W0,xb,hb) #Erosion of beach - no segrass (in lagoon)
    R01=[temp1,temp2,temp];R_01=mean(R01);
    
    if len(Etasimple1)<>0 and (MgPst+MrPst+SgPst)<>0: #If there is vegetation
        Hs,Etas=WaveModelSimple(x,y,Hp1,T) 
        temp=argmin(Etas)
        xb=x[-1]-x[temp];hb=y[temp] #Breaking wave height     
        Rinf,temp1=ErosionKD(A,Hp1,Rnp1,B0,D0,W0,0,0) #Erosion of beach - no segrass (in lagoon)
        Rinf,temp2=ErosionKD(A1,Hp1,Rnp1,B0,D0,W0,0,0) #Erosion of beach - no segrass (in lagoon)
        Rinf,temp=ErosionKD(A1,Hp1,Rnp1,B0,D0,W0,xb,hb) #Erosion of beach - no segrass (in lagoon)
        R13=R_01*Rmax1/Rnp1
        R_1=[temp1,temp2,temp,R13];R1=mean(R_1);
        
        #After manamgement action
    C0=g*T/(2*pi) # deep water phase speed
    hb=(((Ho2**2.0)*C0)/2)**(2.0/5.0)/(g**(1.0/5.0)*0.73**(4.0/5.0));Hb=0.73*hb # breaking wave depth
    xb=(hb/A)**1.5 # surf zone width
    
    #Wave height over eq. profile
    Hs,Etas=WaveModelSimple(x,y,Ho2,T) 
    temp=argmin(Etas)
    xb=x[-1]-x[temp];hb=y[temp] #Breaking wave height     
    Rinf,temp1=ErosionKD(A,Ho2,Rmax2,B0,D0,W0,0,0) #Erosion of beach - no segrass (in lagoon)
    Rinf,temp2=ErosionKD(A1,Ho2,Rmax2,B0,D0,W0,0,0) #Erosion of beach - no segrass (in lagoon)
    Rinf,temp=ErosionKD(A1,Ho2,Rmax2,B0,D0,W0,xb,hb) #Erosion of beach - no segrass (in lagoon)
    R02=[temp1,temp2,temp];R_02=mean(R02);
        
    if len(Etasimple2)<>0 and (MgPst+MrPst+SgPst)<>0: #If there is vegetation
        Hs,Etas=WaveModelSimple(x,y,Hp2,T) 
        temp=argmin(Etas)
        xb=x[-1]-x[temp];hb=y[temp] #Breaking wave height     
        Rinf,temp1=ErosionKD(A,Hp2,Rnp2,B0,D0,W0,0,0) #Erosion of beach - no segrass (in lagoon)
        Rinf,temp2=ErosionKD(A1,Hp2,Rnp2,B0,D0,W0,0,0) #Erosion of beach - no segrass (in lagoon)
        Rinf,temp=ErosionKD(A1,Hp2,Rnp2,B0,D0,W0,xb,hb) #Erosion of beach - no segrass (in lagoon)
        R13=R_02*Rmax2/Rnp2
        R_2=[temp1,temp2,temp,R13];R2=mean(R_2);
        
        #PLOT OUTPUTS
        # beach profile plot
        if D0==0:
            W0=100
        Xp=num.arange(0,1000,1)
        Yp=m*Xp
        loc=find(Yp>=B0);
        Yp[loc]=B0
        Lberm=loc[0];Ltoe=Indexed(Xp,Xp[Lberm+W0])
        Yp[Ltoe:-1]=B0+D0
        
        Yp0=m*Xp-m*R1 # profile of erosion before management action
        loc=find(Yp0>=B0);
        if R1 <= W0:
            Yp0[loc]=B0
            Lberm0=loc[0]
            W00=W0-(Xp[Lberm0]-Xp[Lberm])
            Ltoe0=Indexed(Xp,Xp[Lberm0+W00])
            Yp0[Ltoe0:-1]=B0+D0
        else:
            Yp0[loc]=B0+D0
            Ltoe0=loc[0]
        
        Ypv=m*Xp-m*R2 # profile of erosion after management action
        loc=find(Ypv>=B0);
        if R2 <= W0:
            Ypv[loc]=B0
            Lbermv=loc[0]
            W0v=W0-(Xp[Lbermv]-Xp[Lberm])
            Ltoev=Indexed(Xp,Xp[Lbermv+W0v])
            Ypv[Ltoev:-1]=B0+D0
        else:
            Ypv[loc]=B0+D0
        
        figure(2)
        subplot(211)
        plot(Xp,Yp,Xp,Yp0,'--r',Xp,Ypv,'--g');grid()
        xlim(0,Xp[Ltoe0]+R1+5);ylim(Yp[0],B0+D0+1)
        ylabel('Backshore Elevation [m]',size='large',weight='demi')
        title('Erosion Profiles',size='large',weight='bold')
        legend(('Initial Profile','No Vegetation','Vegetation Present'),'upper left')
        
        subplot(212)
        temp1=Indexed(Xp,Xp[Lberm]-10);temp2=Indexed(Xp,Xp[Ltoe0]+10)
        plot(Xp,Yp,Xp,Yp0,'--r',Xp,Ypv,'--g');grid()
        xlim(Xp[temp1],Xp[temp2]+R1);ylim(Yp[temp1],B0+D0+0.5)
        ylabel('Backshore Elevation [m]',size='large',weight='demi')
        xlabel('Cross-Shore Distance [m]',size='large',weight='demi')
        savefig(outputws+"ErosionPlot2_"+subwsStr+".png",dpi=(640/8))
        Fig2=1; #for HTML
        
elif (MgPst+MrPst+SgPst)<>0 and Beach==0: #Compute erosion amount for consolidated sediments
    rho=1024.0;nu=1.36e-6;d50=0.03
    ks=2.5*d50;kap=0.4
    
    ##Before management action
    #Uc1=0.5 #Current speed before management action
    #us1=0.01;zo1=0.01;dif=10 #Initial value for u* and zo    
    #while dif>1e-4:
        #zo2=ks/30*(1-num.exp(-us1*ks/(27*nu)))+nu/(9*us1)
        #us2=kap*Uc1/(num.log(h[cc]/zo2)-1)
        #dif1=abs(us1-us2);dif2=abs(zo1-zo2);dif=dif1+dif2;
        #zo1=zo2;us1=us2;
    #Tc=rho*us1**2 #Shear stress due to current
    
    #Rws=Ubot1[cc]**2*T/(2*num.pi)/nu
    #fw=0.0521*Rw**(-0.187) #Smooth turbulent flow
    #Tw=0.5*rho*fw*Ubot1[cc]**2
    
    #temp=Tc*(1+1.2*(Tw/(Tc+Tw))**3.2)
    #Trms=(temp**2+0.5*Tw**2)**0.5
    
    ##Erosion 
    #Te=0.0012*Cm**1.2;dmdt=0 #Erosion threshold
    #if Trms>Te:
        #dmdt=me*(Tmrs-Te) #Erosion rate
    #R1=3600*Dur*dmdt/Cm #Depth of bed eroded [cm]
    
    ##After management action
    #Uc2=0.5 #Current speed before management action
    #us1=0.01;zo1=0.01;dif=10 #Initial value for u* and zo    
    #while dif>1e-4:
        #zo2=ks/30*(1-num.exp(-us1*ks/(27*nu)))+nu/(9*us1)
        #us2=kap*Uc2/(num.log(h[cc]/zo2)-1)
        #dif1=abs(us1-us2);dif2=abs(zo1-zo2);dif=dif1+dif2;
        #zo1=zo2;us1=us2;
    #Tc=rho*us1**2 #Shear stress due to current
    
    #Rws=Ubot2[cc]**2*T/(2*num.pi)/nu
    #fw=0.0521*Rw**(-0.187) #Smooth turbulent flow
    #Tw=0.5*rho*fw*Ubot2[cc]**2
    
    #temp=Tc*(1+1.2*(Tw/(Tc+Tw))**3.2)
    #Trms=(temp**2+0.5*Tw**2)**0.5
    
    ##Erosion 
    #Te=0.0012*Cm**1.2;dmdt=0 #Erosion threshold
    #if Trms>Te:
        #dmdt=me*(Tmrs-Te) #Erosion rate
    #R2=3600*Dur*dmdt/Cm*100 #Depth of bed eroded [cm]
    
    
    

