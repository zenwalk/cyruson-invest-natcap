# Marine InVEST: Coastal Protection (Erosion and Valuation)
# Authors: Greg Guannel, Gregg Verutes, Apollo Yi
# 09/08/10

import numpy as num
import CPf_SignalSmooth as SignalSmooth
import string, sys, os, time, datetime, shlex
import fpformat, operator
import arcgisscripting

from win32com.client import Dispatch
from scipy.interpolate import interp1d
from scipy.special import erf
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
    ProfileGeneratorWS = gp.GetParameterAsText(0)
    parameters.append("Workspace: "+ gp.workspace)
    InputTable = gp.GetParameterAsText(1)
    parameters.append("Profile Generator Excel Table: "+ InputTable)
    StormDur = int(gp.GetParameterAsText(2))
    parameters.append("Storm Duration (in hours): "+ str(StormDur))
    HAT = int(gp.GetParameterAsText(3))
    parameters.append("Surge Level During Storm (in meters): "+ str(HAT))
    VegType = gp.GetParameterAsText(4)
    parameters.append("Vegetation Type: "+ VegType)
except:
    raise Exception, msgArguments + gp.GetMessages(2)


# intermediate and output directories
outputws = ProfileGeneratorWS + os.sep + "Output" + os.sep
interws = ProfileGeneratorWS + os.sep + "intermediate" + os.sep
scratchws = ProfileGeneratorWS + os.sep + "scratch" + os.sep

BathyProfile = outputws + "BathyProfile.txt"
CreatedProfile = outputws + "CreatedProfile.txt"
Erosion_Plot = outputws + "Erosion_Plot.png"
ProfileErosion_HTML = outputws + "ProfileErosion_Results.html"

# various functions and checks
def gradient(f,z):
    length=len(f)
    df=length*[0.0]
    df[0]=(f[1]-f[0])/z
    for i in range(1,length-1):
        df[i]=(f[i+1]-f[i-1])/(2.0*z)
    df[length-1]=(f[length-1]-f[length-2])/z
    return df

def FindRoot( fun, a, b, tol = 1e-16 ):
  a = float(a)
  b = float(b)
  assert(sign(fun(a)) != sign(fun(b)))  
  c = (a+b)/2
  while math.fabs(fun( c )) > tol:
    if a == c or b == c: 
      break
    if sign(fun(c)) == sign(fun(b)):
      b = c
    else:
      a = c
    c = (a+b)/2
  return c

def iterativek(sigma, dh):
    kestimated=(sigma**2)/(g*(sqrt(tanh((sigma**2)*dh/g))))
    kprevious=0.0000001
    count=0
    while (abs(kestimated-kprevious) > 0.000005) and (count < 1000):
        count += 1
        kh=kestimated*dh
        kcalculated=(sigma**2)/(tanh(kh)*g)
        kprevious=kestimated
        kestimated=kcalculated
    qk=kcalculated
    return qk

def WaveModelSimple(h,Ho,T,alph,hv,bv,Nv):
    # constants
    g=9.81;rho=1024.0;B=1.0;
    # initialize vectors
    lx=len(h);L=lx*[0.0]
    H=lx*[0.0];Eta=lx*[0.0]
    Db=lx*[0.0];Dv=lx*[0.0];Df=lx*[0.0]
    Er=lx*[0.0];Ef=lx*[0.0];Br=lx*[0.0];
    C=lx*[0.0];n=lx*[0.0];Cg=lx*[0.0]
    k=lx*[0.0];Gam=lx*[0.0];    
    # wave parameter at 1st grid pt
    ash=h[:] # ash is same as h, but is now an independent variable.
    fp=1.0/T; sig=2.0*pi*fp;
    k[0]=iterativek(sig,h[0]) # wave number at 1st grid pt
    L[0]=2.0*pi/k[0] # wave length @ 1st grid pt
    n[0]=0.5*(1+(2.0*k[0]*h[0]/sinh(2.0*k[0]*h[0]))) # to compute Cg at 1st grid pt
    C[0]=L[0]/T;Cg[0]=C[0]*n[0] # phase and group velocity at 1st grid pt

    So=Ho/L[0] # deep water wave steepness

    # RMS wave height at first grid point; Assume no dissipation occurs
    Ho=Ho/sqrt(2.0) # transform significant wave height to rms wave height
    Co=g*T/(2.0*pi); Cgo=Co/2.0 # deep water phase and group speed
    if h[0]>0.5*L[0]: H[0]=Ho # we are in deep water
    else: H[0]=Ho*sqrt(Cgo/Cg[0]) # we are in intermediate water. Assume no brkg occured

    # wave and roller energy
    Ew=0.125*rho*g*H[0]**2;Ef[0]=Ew*Cg[0] # energy flux @ 1st grid pt
    Gam=0.5+0.4*num.tanh(33.0*So); # Gam from Battjes and Stive 85, as per Alsina & Baldock

    # vegetation parameters
    Cd=.1

    # begin wave model based on Thornton & Guza
    for xx in range(lx-1): # transform waves, take MWL into account
        Ef[xx]=0.125*rho*g*(H[xx]**2.0)*Cg[xx] # Ef at (xx)      
        Ef[xx+1]=Ef[xx]-dx*(Db[xx]+Df[xx]+Dv[xx]) # Ef at [xx+1] 
        Br[xx+1]=Br[xx]-dx*(g*Er[xx]*sin(beta)/C[xx]-0.5*Db[xx]) # roller flux
        
        k[xx+1]=iterativek(sig,h[xx+1])
        n[xx+1]=0.5*(1.0+(2.0*k[xx+1]*h[xx+1]/sinh(2.0*k[xx+1]*h[xx+1])))
        C[xx+1]=sig/k[xx+1];Cg[xx+1]=C[xx+1]*n[xx+1] # phase and group velocity
        
        H[xx+1]=num.sqrt(8.0*Ef[xx+1]/(rho*g*Cg[xx+1])) # wave height at [xx+1]
        Hb=0.88/k[xx+1]*tanh(Gam*k[xx+1]*h[xx+1]/0.88)
        
        Var=0.25*rho*g*fp*B;
        temp1=((Hb/H[xx+1])**3.0+1.5*Hb/H[xx+1])*num.exp(-(Hb/H[xx+1])**2.0);
        if isreal(H[xx+1]): temp2=0.75*num.sqrt(pi)*(1-erf(Hb/H[xx+1]));
        else: temp2=0;
        Db[xx+1]=Var*H[xx+1]**3/h[xx+1]*(temp1+temp2); # dissipation due to brkg

        temp1=rho*Cd*bv[xx+1]*Nv[xx+1]*(k[xx+1]*g/(2.0*sig))**3.0/(2.0*sqrt(pi))
        temp2=sinh(k[xx+1]*alph[xx+1]*h[xx+1])**3.0+3.0*sinh(k[xx+1]*alph[xx+1]*h[xx+1])
        temp3=(3.0*k[xx+1]*cosh(k[xx+1]*h[xx+1])**3.0)
        Dv[xx+1]=temp1*temp2/temp3*H[xx+1]**3.0 # diss due to vegetation
        
        Df[xx+1]=rho*Cf/(16.0*sqrt(pi))*(2.0*pi*fp*H[xx+1]/sinh(k[xx+1]*h[xx+1]))**3.0 # diss due to bot friction 
        Er[xx+1]=Br[xx+1]/(C[xx+1]) # roller energy

    Ew1=lx*[0.0]
    Ew1=[0.125*rho*g*(H[i]**2.0) for i in range(lx)]

    # estimate MWS
    val=1.0;dell=1.0;O=0.0 # for setup calc

    Sxx=lx*[0.0]; Rxx=lx*[0.0]
    # force on plants if they were emergent; take a portion if plants occupy only portion of wc
    Fx=[rho*g*Cd*bv[i]*Nv[i]*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
    fx=[-alph[i]*Fx[i] for i in range(lx)] 
    
    while dell>1e-10: # iterate until convergence of water level
        val1=val
        h=[ash[i]+Eta[i] for i in range(lx)] # water depth
        
        Sxx=[0.5*Ew1[i]*(4.0*k[i]*h[i]/sinh(2.0*k[i]*h[i])+1.0) for i in range(lx)] # wave radiation stress
        Rxx=[2.0*Er[i] for i in range(lx)] # roller radiation stress
        # estimate MWL along Xshore transect
        Tem=[Sxx[i]+Rxx[i] for i in range(lx)]
        Temp=gradient(Tem, dx)
        Terms=[-Temp[i]+fx[i] for i in range(lx)]
        Integr=[Terms[i]/(rho*g*h[i]) for i in range(lx)]
        Eta[0]=-0.125*H[0]**2.0*k[0]/sinh(2.0*k[0]*h[0])
        Eta[1]=Eta[0]+Integr[0]*dx

        for i in range(1,lx-2):
            Eta[i+1]=Eta[i-1]+Integr[i]*2*dx
        
        Eta[lx-1]=Eta[lx-2]+Integr[lx-1]*dx
    
        Temp_eta=gradient(Eta, dx)
        valn=[Temp[i]+rho*g*h[i]*Temp_eta[i]-fx[i] for i in range(lx)]
        valn[lx-3]=0;valn[lx-2]=0;valn[lx-1]=0
        dell1=[valn[i]-val1 for i in range(lx)]
        dell1me=sum(dell1)/len(dell1)
        dell1mx=max(dell1)
        
        O=O+1
        if O<5: dell=13
        else: dell=max(dell1me,dell1mx)
            
    return Eta, H, Db

# erosion model
def ErosionKD(B,D,W):
    # constants
    BD=D+B;

    # bkg info in absence of veg.
    C0=g*T/(2*pi) # deep water phase speed
    hb=(((Ho**2.0)*C0)/2)**(2.0/5.0)/(g**(1.0/5.0)*0.73**(4.0/5.0));Hb=0.73*hb # breaking wave depth
    xb=(hb/A)**1.5 # surf zone width

    # erosion model
    Term1=xb-hb/m;
    if Term1<0: # profile factors wrong; can't compute erosion
        print "'Foreshore slope too flat or Sediment scale factor A is too high'"
        R1=nan # can't compute
        R0=nan
        
    else: # compute erosion
        Rinf=(S*Term1)/(BD+hb-S/2)-W*(B+hb+0.5*S)/(BD+hb-0.5*S);# potential erosion distance
        if Rinf<0:
            print "'Berm is so wide that dune is not eroding'"
            Rinf=(S*Term1)/(B+hb-S/2)
            
        TS=(320.0*(Hb**(3.0/2.0)/(g**.5*A**3.0))*(1.0/(1.0+hb/BD+(m*xb)/hb)))/3600.0;# erosion response time scale
        beta=2.0*pi*(TS/StormDur);

        expr="num.exp(-2*x/beta)-num.cos(2*x)+(1/beta)*num.sin(2*x)"# solve this numerically
        fn = eval( "lambda x: " + expr )
        z=FindRoot(fn,pi,pi/2.0) # find zero in function,initial guess from K&D

        R0=0.5*Rinf*(1.0-cos(2.0*z)); # final erosion distance

    return Rinf, R0

# constants
g=9.81
rho=1024.0

# input from user via Excel SS
xlApp=Dispatch("Excel.Application")
xlApp.Visible=0
xlApp.DisplayAlerts=0
xlApp.Workbooks.Open(InputTable)
cell0=xlApp.Worksheets("Profile Generator Input")
cell1=xlApp.Worksheets("Erosion Model Input")

# data from 'Profile Generator Input'
Slope = cell0.Range("f31").Value # foreshore slope = 1/Slope
m = 1.0/Slope; # bed slope

MSL = cell0.Range("d23").Value # mean sea level
HT = cell0.Range("e23").Value # high tide elevation
##HAT = cell0.Range("f23").Value # high tide elevation
A = cell0.Range("e94").Value # sediment scale factor
He = cell0.Range("h18").Value # effective wave height
Hm = cell0.Range("i18").Value # modal wave height
Tm = cell0.Range("j18").Value # modal wave period
hc = 1.57*He # closure depth 

# data from 'Erosion Model Input'
CheckWaveInput=cell1.Range("b67").Value
if CheckWaveInput==9:
    Ho=cell1.Range("d14").Value # wave height
    T=cell1.Range("e14").Value # wave period
elif CheckWaveInput==10:
    WindD1=cell1.Range("d17:k17").Value # wind direction
    WindU1=cell1.Range("d18:k18").Value # wind speed
    WindD2=cell1.Range("d19:k19").Value # wind direction
    WindU2=cell1.Range("d20:k20").Value # wind speed

    ## insert equation for wave from wind (fetch)
    ## inputs: estimate average depth over profile, use distance of fetch vectors, inputs from that WW3 portion of the Erosion sheet
    ## CHECKS FOR FETCH if CheckWaveInput=5-8 or 10 --> then need accurate GIS fetch file.  Check GIS fetch layer has been calculated properly and is present.

## GV ADD THE CODE TO GRAB DISTANCE MEASURES FROM FETCH VECTORS
## ADD CHECK THAT IF FETCH VECTOR GIS LAYER DOESN'T EXIST, ADD ERROR MESSAGE
    
##    U=u(kk);
##    ds=g*d/U^2;Fs=g*F/U^2;
##    A=tanh(0.343*ds^1.14);
##    B=tanh(4.14e-4*Fs^0.79/A);
##    H(kk)=0.24*U^2/g*(A*B)^0.572;
##    
##    A=tanh(0.1*ds^2.01);
##    B=tanh(2.77e-7*Fs^1.45/A);
##    T(kk)=7.69*U/g*(A*B)^0.187;

BermCrest=cell0.Range("f41").Value
BermLength=cell0.Range("g41").Value
DuneCrest=cell0.Range("j51").Value

# management action input from user
MangmtActn=cell1.Range("c72").Value;
if MangmtActn==1: # vegetation field is modified
    hov=cell1.Range("e48").Value # height of marsh
    dov=cell1.Range("f48").Value # diameter of marsh
    Nov=cell1.Range("g48").Value # density of marsh
    do1v=cell1.Range("h48").Value;do1v=-do1v; # starting depth of marsh
    do2v=cell1.Range("i48").Value;do2v=-do2v; # ending depth of marsh
else: # backshore is modified
    B0=cell1.Range("e53").Value; # initial berm height
    W0=cell1.Range("f53").Value; # initial berm length
    D0=cell1.Range("i53").Value; # initial dune length
    B1=cell1.Range("e54").Value; # modified berm height
    W1=cell1.Range("f54").Value; # modified berm length
    D1=cell1.Range("i54").Value; # modified dune length
   
xlApp.ActiveWorkbook.Close(SaveChanges=0)
xlApp.Quit()


# bathymetry
# load bathy profile
TextData=open(BathyProfile,"r") # assume that it's this profile
xd=[]
Depth=[]
for line in TextData.readlines():
    linelist=[float(s) for s in line.split("\t")] # split the list by tab delimiter
    xd.append(linelist[0])
    Depth.append(linelist[1])
    
depth=SignalSmooth.smooth(array(Depth),int(len(Depth)/20),'flat')   

# eq. profile info
Y=array(depth);X=array(xd)
hc=-20.0;la=num.argmin(abs(Y-hc)) # locate closure depth
Yeq=Y[0:la];Xeq=X[0:la];

fitfunc=lambda p, ix: p[0]*(ix)**(2.0/3) # target function
errfunc=lambda p, ix, iy: fitfunc(p, ix) - iy # distance to the target function
p0=[0.1] # initial guess for the parameters
p1, success=optimize.leastsq(errfunc, p0[:], args=(Xeq,-Yeq))
A=p1[0] # sediment scale factor

# resample x-axis
dx=1;L=len(depth)
x=num.arange(0,X[-1],dx)
F=interp1d(X,depth);h=-F(x)
out=num.nonzero(h<0.5);out=out[0]
h=num.delete(h,out,None)
x=num.delete(x,out,None)
h=h[::-1]
Xshore=x[-1]+0.5/m # location of the zero WL (not -0.5)

# vegetation variables
if MangmtActn==1: # vegetation field is modified
    out=num.nonzero(h<0.5);out=out[0]
    h=num.delete(h,out,None);L=len(h)
    ex=num.delete(x,out,None)
    
    hv=num.arange(0.0,L,1.0)*0.0;dv=num.arange(0.0,L,1.0)*0.0
    Nv=num.arange(0.0,L,1.0)*0.0
    Sv=do1v+do2v
    
    if Sv<>0: # fill in seagrass vectors
        temp1=num.nonzero(h<=do2v) # depths shallower than deepest depth
        temp2=num.nonzero(h<=do1v) # depths shallower than shallowest depth
        temp1=temp1[0];temp2=temp2[0]
        keep=num.setdiff1d(temp1,temp2) # keep depth range that are between shallowest and deepest
        hv[keep]=hov # hs is zero except where there is seagrass
        Nv[keep]=Nov # same with density
        dv[keep]=dov # same with diamter
        PlntPlot=h[keep];Xplant=ex[keep]
    # end of adding vegetation variables--    
    Nv=Nv.tolist();dv=dv.tolist();hv=hv.tolist()
    alph=hv/h

# run wave model
# constants for wave model
gam=0.42 # breaking coef
B=1.0; beta=0.1 # breaking parameter and roller angle
fp=1.0/T; sig=2.0*pi*fp # wave frequency
Cf=1e-3 # sand friction coefficent

# expected runup and max MWL
Lo=g*T**2.0/(2.0*pi)
Rmax=1.1*(0.35*m*sqrt(Ho*Lo)+sqrt(Ho*Lo*(0.563*m**2.0+0.004))/2.0) # max runup
Emax=0.35*m*num.sqrt(Ho*Lo)

# run model
S=1.1; # tide + runup - sea level
Eta0,H0,Db0=WaveModelSimple(h,Ho,T,alph,L*[0.0],L*[0.0],L*[0.0])

# erosion model (Dean & Kriebel)
if MangmtActn==1:
    Eta,H,Db=WaveModelSimple(h,Ho,T,alph,hv,dv,Nv)

    # percent wave attenuation
    xv=num.argmin(abs(h-do2v))
    AtnH=[(H0[i]-H[i])/H0[i]*100 for i in range(xv,len(H0),1)] # percent attenuation

    # correct runup
    coef0=max(Eta0)/Emax;Etap=max(Eta)/coef0 # corrected MWL at shoreline in presence of veg.
    if Etap<0: Etap=0

    Hp=(Etap/(0.35*m))**2/Lo # Hprime to estimate runup with veg
    Rnp=1.1*(Etap+num.sqrt(Lo*(Hp*0.563*m**2+0.004*Ho))/2) # runup with vegetation

    # erosion
    B=BermCrest;W=BermLength;D=DuneCrest;
    Rinf,R0=ErosionKD(B,D,W)
    Rveg=R0*abs(Rnp-Rmax)/Rmax; # scale erosion distance when vegetation is present by ratio of runup

    # test to print
##    print "'Shore retreat in absence of vegetation is %6.2f'" %R0
##    print "'Shore retreat in presence of vegetation is %6.2f'" %Rveg
    
    subplot(221)
    plot(x,H0[::-1], label='No Veg.')
    plot(x,H[::-1], label='Yes Veg.')
    grid();xlim(x[-1],0)
    ylabel('Wave Height [m]', size='large')
    xlabel('Cross-Shore Distance [m]', size='large')
    legend(loc='lower left', bbox_to_anchor=(0.05, 0.25), ncol=1, fancybox=True, shadow=True)
    
    subplot(222)
    plot(x,Eta0[::-1],x,Eta[::-1])
    grid()
    xlim(x[-1],0)
    ylabel('Mean Water Surf. Elev. [m]', size='medium')
    xlabel('Cross-Shore Distance [m]', size='large')
    
    subplot(223)
    plot(x,-h[::-1])
    ylabel('Depth [m]', size='large')
    grid()
    xlim(x[-1],0)
    
    subplot(224)
    plot(AtnH[::-1])
    ylabel('Percent Wave Attenuation', size='large')
    grid()
    xlim(x[-1],0)

else:
    Rinf0,R0=ErosionKD(B0,D0,W0)
    Rinf1,R1=ErosionKD(B1,D1,W1)
##    print "'Shore retreat using initial configuration is %6.2f'" %R0
##    print "'Shore retreat after backshore modification is %6.2f'" %R1

    subplot(311)
    plot(x,H0[::-1], label='No Veg.')
    plot(x,H[::-1], label='Yes Veg')
    grid()
    xlim(x[-1],0)
    ylabel('Wave Height [m]', size='large')
    xlabel('Cross-Shore Distance [m]', size='large')
    legend(loc='lower left', bbox_to_anchor=(0.05, 0.25), ncol=1, fancybox=True, shadow=True)

    subplot(312)
    plot(x,Eta0[::-1],x,Eta[::-1])
    grid();xlim(x[-1],0)
    ylabel('Mean Water Surf. Elev. [m]', size='medium')
    xlabel('Cross-Shore Distance [m]', size='large')
    
    subplot(313)
    plot(x,-h[::-1])
    ylabel('Depth [m]', size='large')
    grid()
    xlim(x[-1],0)

# plot results
##stop=150;la=num.argmin(abs(x-xb0+stop))
##pos=min([la,xv])
savefig(Erosion_Plot, dpi=(640/8))

# re-open and edit html file
htmlfile = open(ProfileErosion_HTML, "a")
if MangmtActn==1:
    htmlfile.write("<li><u>Shore retreat in absence of vegetation is</u>: "+str(R0)+"<br>\n")
    htmlfile.write("<li><u>Shore retreat in presence of vegetation is</u>: "+str(Rveg)+"<br>\n")
else:
    htmlfile.write("<li><u>Shore retreat using initial configuration is</u>: "+str(R0)+"<br>\n")
    htmlfile.write("<li><u>Shore retreat after backshore modification is</u>: "+str(Rveg)+"<br>\n")
htmlfile.write("</html>")
htmlfile.close()