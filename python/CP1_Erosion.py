# Marine InVEST: Coastal Protection (Wave and Erosion)
# Authors: Greg Guannel, Gregg Verutes, Apollo Yi
# 10/09/10

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
    parameters.append("Date and Time: " + now.strftime("%Y-%m-%d %H:%M"))
    ProfileGeneratorWS = gp.GetParameterAsText(0)
    parameters.append("Workspace: " + gp.workspace)
    subwsStr = gp.GetParameterAsText(1)
    parameters.append("Label for Erosion Run (10 characters max): " + subwsStr)
    InputTable = gp.GetParameterAsText(2)
    parameters.append("Profile Generator Excel Table: " + InputTable)
    CSProfile = gp.GetParameterAsText(3)
    parameters.append("Cross Shore Profile: " + CSProfile)
    StormDur = float(gp.GetParameterAsText(4))
    parameters.append("Storm Duration (in hours): " + str(StormDur))
    S = float(gp.GetParameterAsText(5))
    parameters.append("Water Level Elevation During Storm (in meters): " + str(S))
    AvgD = gp.GetParameterAsText(6)
    parameters.append("If Profile is in Shallow Bay or Estuary, Enter Average Water Depth: " + str(AvgD))
    if AvgD:
        AvgD = float(gp.GetParameterAsText(6))
    else: AvgD = 500.00
except:
    raise Exception, msgArguments + gp.GetMessages(2)

# remove spaces and shorten 'subwsStr' if greater than 10 characters
subwsStr = subwsStr.replace(" ", "")
subwsStr = subwsStr[0:10]

# intermediate and output directories
outputws = ProfileGeneratorWS + os.sep + "_WaveModel_Outputs" + os.sep
interws = ProfileGeneratorWS + os.sep + "scratch" + os.sep

thefolders = ["_WaveModel_Outputs"]
for folder in thefolders:
    if not gp.exists(ProfileGeneratorWS + folder):
        gp.CreateFolder_management(ProfileGeneratorWS + os.sep, folder)

# various functions and checks
def gradient(f, z):
    length = len(f)
    df = length * [0.0]
    df[0] = (f[1] - f[0]) / z
    for i in range(1, length - 1):
        df[i] = (f[i + 1] - f[i - 1]) / (2.0 * z)
    df[length - 1] = (f[length - 1] - f[length - 2]) / z
    return df

def Indexed(x, value): # locates index of point in vector x that has closest value as variable value
    mylist = abs(x - value);    
    if isinstance(x, num.ndarray):
        mylist = mylist.tolist()
    minval = min(mylist)
    ind = [i for i, v in enumerate(mylist) if v == minval]
    ind = ind[0]
    return ind

def FindRootKD(fun, a, b, tol=1e-16):
    a = float(a);b = float(b)
    assert(sign(fun(a)) != sign(fun(b)))
    c = (a + b) / 2
    while math.fabs(fun(c)) > tol:
        if a == c or b == c:
            break
        if sign(fun(c)) == sign(fun(b)):
            b = c
        else:
            a = c
        c = (a + b) / 2
    return c

# wind-wave generation
def WindWave(U, F, d):
    ds = g * d / U ** 2.0;Fs = g * F / U ** 2.0
    A = tanh(0.343 * ds ** 1.14)
    B = tanh(4.14e-4 * Fs ** 0.79 / A)
    H = 0.24 * U ** 2 / g * (A * B) ** 0.572 # wave height
    A = tanh(0.1 * ds ** 2.01)
    B = tanh(2.77e-7 * Fs ** 1.45 / A)
    T = 7.69 * U / g * (A * B) ** 0.18 # wave period
    return H, T

# dispersion relationship
def iterativek(sigma, dh):
    kestimated = (sigma ** 2) / (g * (sqrt(tanh((sigma ** 2) * dh / g))))
    kprevious = 0.0000001
    count = 0
    while (abs(kestimated - kprevious) > 0.000005) and (count < 1000):
        count += 1
        kh = kestimated * dh
        kcalculated = (sigma ** 2) / (tanh(kh) * g)
        kprevious = kestimated
        kestimated = kcalculated
    qk = kcalculated
    return qk

# wave transformation model
def WaveModel(X, h, Ho, T, Seagr, Marsh, Mang, Coral):
    # constants
    g = 9.81;rho = 1024.0;B = 1.0;Beta = 0.1
    lx = len(X);dx = X[1] - X[0]
    fp = 1.0 / T; sig = 2.0 * pi * fp
    Cgo = g * T / (4 * pi)
    
    # initialize vegetation vectors, zeros everywhere except where veg present
    VegLoc = num.arange(0.0, lx, 1.0) * 0.0
    # seagrass
    veg = Seagr[0];loc = Seagr[1];sgrs_chk = 0
    hs = veg[0];ds = veg[1];Ns = veg[2]
    alphs = hs / h;o = num.nonzero(alphs >= 1);alphs[o[0]] = 1

    if loc is not None and hs <> 0:
        sgrs_chk = 1
    else: loc = [0.0, 0.0]

    off = Indexed(X, loc[0]);sho = Indexed(X, loc[1]) # locate loc of seagrass
    VegLoc[off:sho + 1] = 1 # 1 for location of seagrass

    # marsh
    veg = Marsh[0];loc = Marsh[1]
    hr = veg[0];dr = veg[1];Nr = veg[2]
    alphr = hr / h;o = num.nonzero(alphr >= 1);alphr[o[0]] = 1
        
    if loc is not None and hr <> 0:
        temp = 1
    else: loc = [0.0, 0.0];

    off = Indexed(X, loc[0]);sho = Indexed(X, loc[1])
    VegLoc[off:sho + 1] = 2 #2 for location of marshes

    # mangroves
    veg = Mang[0];loc = Mang[1];mgrv_chk = 0
    rt = veg[0];hgr = rt[0];tk = veg[1];cn = veg[2] # roots, trunk and canpy info
    hgr = rt[0];dgr = rt[1];Ngr = rt[2]
    hgt = tk[0];dgt = tk[1];Ngt = tk[2]
    hgc = cn[0];dgc = cn[1];Ngc = cn[2]

    if loc is not None and hgr <> 0:
        mgrv_chk = 1
    else: loc = [0.0, 0.0]

    off = Indexed(X, loc[0]);sho = Indexed(X, loc[1])
    VegLoc[off:sho + 1] = 3 # 3 for location of mangroves
       
    alphgr = hgr / h;alphgt = hgt / h;alphgc = hgc / h
    for kk in range(lx):
        if alphgr[kk] > 1:
            alphgr[kk] = 1;alphgt[kk] = 0;alphgc[kk] = 0
        elif alphgr[kk] + alphgt[kk] > 1:
            alphgt[kk] = 1 - alphgr[kk];alphgc[kk] = 0
        elif alphgr[kk] + alphgt[kk] + alphgc[kk] > 1:
            alphgc[kk] = 1 - alphgr[kk] - alphgt[kk]
                
    # drag coefft for vegetation. mangrove and marsh win over seagrass if X overlap
    Cds = num.arange(0.0, lx, 1.0) * 0.0;o = num.nonzero(VegLoc == 1);Cds[o[0]] = 0.1 # seagrass	
    Cdr = num.arange(0.0, lx, 1.0) * 0.0;o = num.nonzero(VegLoc == 2);Cdr[o[0]] = 0.1 # marsh	
    Cdg = num.arange(0.0, lx, 1.0) * 0.0;o = num.nonzero(VegLoc == 3);Cdg[o[0]] = 1 # mangrove

    # initialize vectors for wave model
    H = lx * [0.0];Eta = lx * [0.0];L = lx * [0.0]
    Db = lx * [0.0];Df = lx * [0.0]
    Ds = lx * [0.0];Dr = lx * [0.0];Dg = lx * [0.0] # seagrass, marsh, mangrove
    Er = lx * [0.0];Ef = lx * [0.0];Br = lx * [0.0]
    C = lx * [0.0];n = lx * [0.0];Cg = lx * [0.0]
    k = lx * [0.0];Gam = lx * [0.0]
    
    # wave parameter at 1st grid pt
    ash = [h[ii] for ii in range(lx)] # ash is same as h, but is now an independent variable.
    fp = 1.0 / T; sig = 2.0 * pi * fp
    k[0] = iterativek(sig, h[0]) # wave number at 1st grid pt
    L[0] = 2.0 * pi / k[0] # wave length @ 1st grid pt
    n[0] = 0.5 * (1 + (2.0 * k[0] * h[0] / sinh(2.0 * k[0] * h[0]))) # to compute Cg at 1st grid pt
    C[0] = L[0] / T;Cg[0] = C[0] * n[0] # phase and group velocity at 1st grid pt

    # RMS wave height at first grid point; Assume no dissipation occurs
    Ho = Ho / sqrt(2.0) # transform significant wave height to rms wave height
    Co = g * T / (2.0 * pi); Cgo = Co / 2.0 # deep water phase and group speed
    if h[0] > 0.5 * L[0]: H[0] = Ho # we are in deep water
    else: H[0] = Ho * sqrt(Cgo / Cg[0]) # we are in intermediate water. Assume no brkg occured

    Ef[0] = 0.125 * rho * g * H[0] ** 2 * Cg[0] # energy flux @ 1st grid pt
    So = Ho / L[0] # deep water wave steepness
    Gam = 0.5 + 0.4 * num.tanh(33.0 * So) # Gam from Battjes and Stive 85, as per Alsina & Baldock

    # begin wave model based on Thornton & Guza
    for xx in range(lx - 1): # transform waves, take MWL into account
        Ef[xx] = 0.125 * rho * g * (H[xx] ** 2.0) * Cg[xx] # Ef at (xx)      
        Ef[xx + 1] = Ef[xx] - dx * (Db[xx] + Df[xx] + Ds[xx] + Dr[xx] + Dg[xx]) # Ef at [xx+1] 

        k[xx + 1] = iterativek(sig, h[xx + 1])
        n[xx + 1] = 0.5 * (1.0 + (2.0 * k[xx + 1] * h[xx + 1] / sinh(2.0 * k[xx + 1] * h[xx + 1])))
        C[xx + 1] = sig / k[xx + 1];Cg[xx + 1] = C[xx + 1] * n[xx + 1] # phase and group velocity

        H[xx + 1] = num.sqrt(8.0 * Ef[xx + 1] / (rho * g * Cg[xx + 1])) # wave height at [xx+1]      
        Br[xx + 1] = Br[xx] - dx * (g * Er[xx] * sin(Beta) / C[xx] - 0.5 * Db[xx]) # roller flux
        Er[xx + 1] = Br[xx + 1] / (C[xx + 1]) # roller energy

        Var = 0.25 * rho * g * fp * B
        Hb = 0.88 / k[xx + 1] * tanh(Gam * k[xx + 1] * h[xx + 1] / 0.88)		
        temp1 = ((Hb / H[xx + 1]) ** 3.0 + 1.5 * Hb / H[xx + 1]) * num.exp(-(Hb / H[xx + 1]) ** 2.0)
        if isreal(H[xx + 1]): temp2 = 0.75 * num.sqrt(pi) * (1 - erf(Hb / H[xx + 1]))
        else: temp2 = 0
        Db[xx + 1] = Var * H[xx + 1] ** 3 / h[xx + 1] * (temp1 + temp2) # dissipation due to brkg

        A = 0.5 * H[xx + 1] * cosh(k[xx + 1] * h[xx + 1]) / sinh(k[xx + 1])
        if Coral == 1: Cf = min(0.3, num.exp(-6.0 + 5.2 * (A / 0.16) ** (-0.19))) * 0 + .15
        else: Cf = min(0.3, num.exp(-6.0 + 5.2 * (A / .001) ** (-0.19))) #Sandy systems

        Df[xx + 1] = rho * Cf / (12.0 * pi) * (2.0 * pi * fp * H[xx + 1] / sinh(k[xx + 1] * h[xx + 1])) ** 3.0 # diss due to bot friction 

        CdDN = Cds[xx + 1] * ds * Ns
        temp1 = rho * CdDN * (k[xx + 1] * g / (2.0 * sig)) ** 3.0 / (2.0 * sqrt(pi))
        temp2 = sinh(k[xx + 1] * alphs[xx + 1] * h[xx + 1]) ** 3.0 + 3.0 * sinh(k[xx + 1] * alphs[xx + 1] * h[xx + 1])
        temp3 = (3.0 * k[xx + 1] * cosh(k[xx + 1] * h[xx + 1]) ** 3.0)
        Ds[xx + 1] = temp1 * temp2 / temp3 * H[xx + 1] ** 3.0 # diss due to seagrass
        
        CdDN = Cdr[xx + 1] * dr * Nr
        temp1 = rho * CdDN * (k[xx + 1] * g / (2.0 * sig)) ** 3.0 / (2.0 * sqrt(pi))
        temp2 = sinh(k[xx + 1] * alphs[xx + 1] * h[xx + 1]) ** 3.0 + 3.0 * sinh(k[xx + 1] * alphs[xx + 1] * h[xx + 1])
        temp3 = (3.0 * k[xx + 1] * cosh(k[xx + 1] * h[xx + 1]) ** 3.0)
        Dr[xx + 1] = temp1 * temp2 / temp3 * H[xx + 1] ** 3.0 # diss due to marshes

        V1 = 3 * sinh(k[xx + 1] * alphgr[xx + 1] * h[xx + 1]) + sinh(k[xx + 1] * alphgr[xx + 1] * h[xx + 1]) ** 3
        V2 = (3 * sinh(k[xx + 1] * (alphgr[xx + 1] + alphgt[xx + 1]) * h[xx + 1]) - 
            3 * sinh(k[xx + 1] * alphgr[xx + 1] * h[xx + 1]) + 
            sinh(k[xx + 1] * (alphgr[xx + 1] + alphgt[xx + 1]) * h[xx + 1]) ** 3 - 
            sinh(k[xx + 1] * alphgr[xx + 1] * h[xx + 1]) ** 3)
        V3 = (3 * sinh(k[xx + 1] * (alphgr[xx + 1] + alphgt[xx + 1] + alphgc[xx + 1]) * h[xx + 1]) - 
            3 * sinh(k[xx + 1] * (alphgr[xx + 1] + alphgt[xx + 1]) * h[xx + 1]) + 
            sinh(k[xx + 1] * (alphgr[xx + 1] + alphgt[xx + 1] + alphgc[xx + 1]) * h[xx + 1]) ** 3 - 
            sinh(k[xx + 1] * (alphgr[xx + 1] + alphgt[xx + 1]) * h[xx + 1]) ** 3)
        
        CdDN = Cdg[xx + 1] * (dgr * Ngr * V1 + dgt * Ngt * V2 + dgc * Ngc * V3)
        temp1 = rho * CdDN * (k[xx + 1] * g / (2.0 * sig)) ** 3.0 / (2.0 * sqrt(pi))
        temp2 = (3 * k[xx + 1] * cosh(k[xx + 1] * h[xx + 1]) ** 3) * H[xx + 1] ** 3
        Dg[xx + 1] = temp1 * temp2
        Dg[xx + 1] = (rho * Cdg[xx + 1] / (2 * sqrt(pi)) * (k[xx + 1] * g / (2 * sig)) ** 3 * 
                  (dgr * Ngr * V1 + dgt * Ngt * V2 + dgc * Ngc * V3) / 
                  (3 * k[xx + 1] * cosh(k[xx + 1] * h[xx + 1]) ** 3) * H[xx + 1] ** 3) # diss. due to mangrove forest

    Ew = lx * [0.0];Ew = [0.125 * rho * g * (H[i] ** 2.0) for i in range(lx)]

    if mgrv_chk == 0:  # compute setup for no veg and seagrass systems only
        # estimate MWS
        val = 1.0;dell = 1.0;O = 0.0 # for setup calc
        Sxx = lx * [0.0]; Rxx = lx * [0.0]
        
        # force on plants if they were emergent; take a portion if plants occupy only portion of wc
        Fxs = [rho * g * Cds[i] * ds * Ns * H[i] ** 3.0 * k[i] / (12.0 * pi * tanh(k[i] * ash[i])) for i in range(lx)]
        fxs = [-alphs[i] * Fxs[i] for i in range(lx)] 
        
        Fxr = [rho * g * Cdr[i] * dr * Nr * H[i] ** 3.0 * k[i] / (12.0 * pi * tanh(k[i] * ash[i])) for i in range(lx)]
        fxr = [-alphr[i] * Fxr[i] for i in range(lx)]

        fx = [fxs[ii] + fxr[ii] for ii in range(lx)]         
        
        while dell > 1e-10: # iterate until convergence of water level
            val1 = val
            h = [ash[i] + Eta[i] for i in range(lx)] # water depth
            
            Sxx = [0.5 * Ew[i] * (4.0 * k[i] * h[i] / sinh(2.0 * k[i] * h[i]) + 1.0) for i in range(lx)] # wave radiation stress
            Rxx = [2.0 * Er[i] for i in range(lx)] # roller radiation stress
            # estimate MWL along Xshore transect
            Tem = [Sxx[i] + Rxx[i] for i in range(lx)]
            Temp = gradient(Tem, dx)
            Terms = [-Temp[i] + fx[i] for i in range(lx)]
            Integr = [Terms[i] / (rho * g * h[i]) for i in range(lx)]
            Eta[0] = -0.125 * H[0] ** 2.0 * k[0] / sinh(2.0 * k[0] * h[0])
            Eta[1] = Eta[0] + Integr[0] * dx

            for i in range(1, lx - 2):
                Eta[i + 1] = Eta[i - 1] + Integr[i] * 2 * dx
            
            Eta[lx - 1] = Eta[lx - 2] + Integr[lx - 1] * dx
        
            Temp_eta = gradient(Eta, dx)
            valn = [Temp[i] + rho * g * h[i] * Temp_eta[i] - fx[i] for i in range(lx)]
            valn[lx - 3] = 0;valn[lx - 2] = 0;valn[lx - 1] = 0
            dell1 = [valn[i] - val1 for i in range(lx)]
            dell1me = sum(dell1) / len(dell1)
            dell1mx = max(dell1)
            
            O = O + 1
            if O < 5: dell = 13
            else: dell = max(dell1me, dell1mx)

    Ubot = [pi * H[ii] / (T * sinh(k[ii] * h[ii])) for ii in range(lx)] #Bottom velocity

    return H, Eta, Ubot, VegLoc

# wave attenuation by coral reefs
def WavesCoral():
    global Slope
    BrkRim = 0;BrkFace = 0
    # reef profile
    Xreef = num.arange(0.0, 1000.0, 1.0)
    Yreef = AlphR * Xreef - he
    loc = num.nonzero(Yreef > -hr);loc = loc[0]
    Yreef[loc] = -hr # reef profile
    
    kp = cell3.Range("c2:c202").Value;kp = num.array(kp) # reef shape factor
    Kp = 0.0

    D = T * sqrt(g / hr) # relative subm
    Lo = g * T ** 2.0 / (2.0 * pi)
    Hrms = Ho / sqrt(2) # Hrms wave height
    if D == 8:  # first check for breaking location
        if Hrms >= 0.5 * he:
            BrkFace = 1; # wave break on face
            loc = Indexed(TanAlph, AlphF)
            Kp = kp[loc] # reef shape factor

            ha = hr # rep. depth
            
        elif Hrms >= 0.5 * hr:
            BrkRim = 1 # wave break on rim
            loc = Indexed(TanAlph, AlphR)
            Kp = kp[loc] # reef shape factor

            db = 0.259 * Hrms * (tan(AlphR) ** 2 * Hrms / Lo) ** (-0.17) # breaking wave height
            if db < hr: db = hr
            elif db > he: db = he

            xs = (2 + 1.1 * Hrms / he) * T * (g * he) ** .5 # surf zone width
            loc1 = Indexed(Yreef, -db); # start brkg
            loc2 = Indexed(Xreef, Xreef[loc1] + xs) # end brkg
            ha = -average(Yreef[loc1:loc2 + 1]) # rep. depth
    else:  # second check for breaking location
        if Hrms >= 0.4 * he:
            BrkFace = 1 # wave break on face
            loc = Indexed(TanAlph, AlphF)
            Kp = kp[loc] # reef shape factor
            ha = hr
            
        elif Hrms >= 0.4 * hr:
            BrkRim = 1 # wave break on rim
            loc = Indexed(TanAlph, AlphR)
            Kp = kp[loc] # reef shape factor

            db = 0.259 * Hrms * (tan(AlphR) ** 2.0 * Hrms / Lo) ** (-0.17) # breaking wave height
            if db < hr: db = hr
            elif db > he: db = he

            xs = (2 + 1.1 * Hrms / he) * T * (g * he) ** .5 # surf zone width
            loc1 = Indexed(Yreef, -db) # start brkg
            loc2 = Indexed(Xreef, Xreef[loc1] + xs) # end brkg
            ha = -average(Yreef[loc1:loc2 + 1]) # rep. depth

    # wave transformation on top of coral reef
    Etar = 0.0;delta = 10
    dx = 0.1;Xflat = num.arange(0.0, Wr, dx)
    lx = len(Xflat); H_r = lx * [0.0]
    if Kp <> 0:
        # wave Setup
        while delta > 0.01:
            EtaR = 3.0 / (64.0 * pi) * Kp * Hrms ** 2 * T * g ** .5 / ((Etar + ha) ** 1.5)
            delta = abs(EtaR - Etar);Etar = EtaR
        Etar = Etar[0]
        
        # wave height profile over reef
        H_r[0] = 0.42 * (hr + Etar)
    else:
        H_r[0] = Hrms
        
    for xx in range(lx - 1):
        H_r[xx + 1] = H_r[xx] - 0.15 * H_r[xx] ** 2.0 / (3.0 * pi * (Etar + hr) ** 2.0) * dx
        
    # process-based method
    if AlphF < 0.8:
        Seagr0 = [[0, 0, 0], [0, 0]];Marsh0 = [[0, 0, 0], [0, 0]]
        Mang0 = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [0, 0]]
        H, Eta, temp, temp = WaveModel(X, h, Ho, T, Seagr0, Marsh0, Mang0, 1)
    else:
        H = [0 for ii in range(len(X))]
        Eta = H[:]
    
    return H, Eta, H_r, Etar

# wave attenuation by reef breakwater
def BreakwaterKt(Hi, T, hi, hc, Cwidth, Bwidth):
    Lo = 9.81 * T ** 2.0 / (2.0 * pi)
    Rc = hc - hi # depth of submergence
    Boff = (Bwidth - Cwidth) / 2.0 # base dif on each side
    ksi = (hc / Boff) / sqrt(Hi / Lo)

    # van der Meer (2005)
    Kt1 = -0.4 * Rc / Hi + 0.64 * (Cwidth / Hi) ** (-.31) * (1.0 - num.exp(-0.5 * ksi)) # transmission coeff: d'Angremond
    Kt1 = max(Kt1, 0.075);Kt1 = min(0.8, Kt1);

    Kt2 = -0.35 * Rc / Hi + 0.51 * (Cwidth / Hi) ** (-.65) * (1.0 - num.exp(-0.41 * ksi)) # transmission coeff: van der Meer
    Kt2 = max(Kt2, 0.05);Kt2 = min(Kt2, -0.006 * Cwidth / Hi + 0.93)

    if Cwidth / Hi < 8.0: # d'Angremond
        Kt = Kt1
    elif Cwidth / Hi > 12.0: # van der Meer
        Kt = Kt2
    else: # linear interp
        temp1 = (Kt2 - Kt1) / 4.0;temp2 = Kt2 - temp1 * 12.0
        Kt = temp1 * Cwidth / Hi + temp2
    return Kt

# K&D Erosion model
def ErosionKD(Ho, RunupVal, B, D, W):
    # constants
    BD = D + B;TWL = S + RunupVal
    if TWL > B0 + D0:
        TWL = B0 + D0

    # bkg info in absence of veg.
    C0 = g * T / (2 * pi) # deep water phase speed
    hb = (((Ho ** 2.0) * C0) / 2) ** (2.0 / 5.0) / (g ** (1.0 / 5.0) * 0.73 ** (4.0 / 5.0));Hb = 0.73 * hb # breaking wave depth
    xb = (hb / A) ** 1.5 # surf zone width

    # erosion model
    Term1 = xb - hb / m
    if Term1 < 0: # profile factors wrong; can't compute erosion
        gp.AddMessage("Foreshore slope too flat or sediment size is too high.")
        Rinf = 0 # can't compute
        R0 = 0
        
    else: # compute erosion
        Rinf = (TWL * Term1) / (BD + hb - TWL / 2) - W * (B + hb + 0.5 * TWL) / (BD + hb - 0.5 * TWL) # potential erosion distance
        if Rinf < 0:
            gp.AddMessage("Berm is so wide that dune is not eroding.")
            Rinf = (TWL * Term1) / (B + hb - TWL / 2)
            
        TS = (320.0 * (Hb ** (3.0 / 2.0) / (g ** .5 * A ** 3.0)) * (1.0 / (1.0 + hb / BD + (m * xb) / hb))) / 3600.0 # erosion response time scale
        global BetaKD
        BetaKD = 2.0 * pi * (TS / StormDur)

        expr = "num.exp(-2.0*x/BetaKD)-num.cos(2.0*x)+(1.0/BetaKD)*num.sin(2.0*x)" # solve this numerically
        fn = eval("lambda x: " + expr)
        z = FindRootKD(fn, pi, pi / 2) # find zero in function,initial guess from K&D

        R0 = 0.5 * Rinf * (1.0 - cos(2.0 * z)) # final erosion distance

    return Rinf, R0

# code starts
gp.AddMessage("\nReading inputs...")
# constants
g = 9.81;rho = 1024.0
# input from user via Excel SS
xlApp = Dispatch("Excel.Application")
xlApp.Visible = 0
xlApp.DisplayAlerts = 0
xlApp.Workbooks.Open(InputTable)
cell = xlApp.Worksheets("ProfileGeneratorInput")
cell1 = xlApp.Worksheets("WaveModelInput")
cell2 = xlApp.Worksheets("CreateModifyProfile")
cell3 = xlApp.Worksheets("ReefShapeFactor")
d50 = cell.Range("e8").Value # sediment size

# wave height and period
CheckWaveInput = cell1.Range("b86").Value
if CheckWaveInput < 5:
    Ho = cell1.Range("h13").Value # wave height
    T = cell1.Range("i13").Value # wave period

elif CheckWaveInput >= 5 and CheckWaveInput < 9:
    WindU1 = cell1.Range("d23:k23").Value;WindU1 = WindU1[0] # wind speed
    WindU2 = cell1.Range("d25:k25").Value;WindU2 = WindU2[0]
    U = [WindU1[ii] for ii in range(len(WindU1))];V = [WindU2[ii] for ii in range(len(WindU2))]
    U.extend(V) # wind vector 16 directions

    Fetch = cell1.Range("e100:t100").Value;Fetch = Fetch[0];Fetch = num.array(Fetch) # fetch distance
    Hf = 16 * [0.0];Tf = 16 * [0.0]
    for ff in range(15): # compute wave height in each fetch direction
        if U[ff] <> 0:
            Hf[ff], Tf[ff] = WindWave(U[ff], Fetch[ff], AvgD)

    Pow = [Tf[kk] * Hf[kk] ** 2.0 for kk in range(len(Hf))]
    loc = Indexed(Pow, max(Pow))
    Ho = Hf[loc];T = Tf[loc] # take wave height/period that has highest power
    cell1.Range("e104:t104").Value = Hf
    cell1.Range("e105:t105").Value = Tf
        
elif CheckWaveInput == 9:
    Ho = cell1.Range("d13").Value # wave height
    T = cell1.Range("e13").Value # wave period
    
elif CheckWaveInput == 10:
    WindU1 = cell1.Range("d17:k17").Value;WindU1 = WindU1[0] # wind speed
    WindU2 = cell1.Range("d19:k19").Value;WindU2 = WindU2[0]
    U = [WindU1[ii] for ii in range(len(WindU1))];V = [WindU2[ii] for ii in range(len(WindU2))]
    U.extend(V) # wind vector 16 directions
    
    Fetch = cell1.Range("e100:t100").Value;Fetch = Fetch[0];Fetch = num.array(Fetch) # fetch distance
    Hf = 16 * [0.0];Tf = 16 * [0.0];
    for ff in range(15): # compue wave height in each fetch direction
        if U[ff] <> 0:    Hf[ff], Tf[ff] = WindWave(U[ff], Fetch[ff], AvgD)
        
    Pow = [Tf[kk] * Hf[kk] ** 2.0 for kk in range(len(Hf))]
    loc = Indexed(Pow, max(Pow))
    Ho = Hf[loc];T = Tf[loc] # take wave height/period that has highest power
    cell1.Range("e104:t104").Value = Hf
    cell1.Range("e105:t105").Value = Tf

# foreshore slope input from user
m = cell2.Range("f8").Value
Fig2 = 0; # for HTML

# management action input from user
MangmtActn = cell1.Range("b88").Value

if MangmtActn == 1: # vegetation field is modified
    gp.AddMessage("\nReading depth profile...")
    # read in user's cross-shore profile
    TextData = open(CSProfile, "r") 
    X_lst = [];h_lst = []
    for line in TextData.read().strip("\n").split("\n"):
        linelist = [float(s) for s in line.split("\t")] # split the list by comma delimiter
        X_lst.append(linelist[0])
        h_lst.append(linelist[1])
    TextData.close()
    
    X = num.array(X_lst);h = num.array(h_lst);
    keep = num.nonzero(h < -0.5)
    h = -h[keep];X = X[keep];lx = len(X)
    if h[0] < h[-1]:
       h = h[::-1] # reverse order if profile starts at shoreline

    hos = cell1.Range("e42").Value # height of seagrass
    dos = cell1.Range("f42").Value # diameter of seagrass
    Nos = cell1.Range("g42").Value # density of seagrass
    Xo1s = cell1.Range("h42").Value # offshore edge of seagrass
    Xo2s = cell1.Range("i42").Value # shoreward edge of seagrass
    veg = [hos, dos, Nos];loc = [Xo1s, Xo2s]
    Seagr = [veg, loc] # seagrass characteristics
    Seagr0 = [[0, 0, 0], [0, 0]]

    hor = cell1.Range("e43").Value # height of marsh
    dor = cell1.Range("f43").Value # diameter of marsh
    Nor = cell1.Range("g43").Value # density of marsh
    Xo1r = cell1.Range("h43").Value # offshore edge of marsh
    Xo2r = cell1.Range("i43").Value # shoreward edge of marsh
    veg = [hor, dor, Nor];loc = [Xo1r, Xo2r]
    Marsh = [veg, loc] # marsh characteristics
    Marsh0 = [[0, 0, 0], [0, 0]]

    hogr = cell1.Range("e39").Value # height of mangrove roots
    dogr = cell1.Range("f39").Value # diameter of mangrove roots
    Nogr = cell1.Range("g39").Value # density of mangrove roots
    hogt = cell1.Range("e40").Value # height of mangrove trunks
    dogt = cell1.Range("f40").Value # diameter of mangrove trunks
    Nogt = cell1.Range("g40").Value # density of mangrove trunks
    hogc = cell1.Range("e41").Value # height of mangrove canopy
    dogc = cell1.Range("f41").Value # diameter of mangrove canopy
    Nogc = cell1.Range("g41").Value # density of mangrove canopy
    Xo1g = cell1.Range("h39").Value # offshore edge of mangrove
    Xo2g = cell1.Range("i39").Value # shoreward edge of mangrove

    veg = [[hogr, dogr, Nogr], [hogt, dogt, Nogt], [hogc, dogc, Nogc]];loc = [Xo1g, Xo2g]
    Mang = [veg, loc] # mangrove characteristics
    Mang0 = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [0, 0]]

    # modify the depth profile appropriately
    if hor <> 0 or hogr + hogt + hogc <> 0: # if there's a marsh or a mangrove
        h = h + S #Add surge level
        out = num.nonzero(h > 0.5);out = out[0]
        h = h[out];X = X[out] # only keep values that are below water
        m = abs(h[-1] - h[-10]) / 10 # average slope 10m from end of transect
    else:
        out = num.nonzero(h > 0.5);out = out[0]
        h = h[out];X = X[out] # only keep values that are below water

    # run wave model in presence and absence of vegetation
    gp.AddMessage("\nComputing wave height profiles...")
    H, Eta, Ubot, temp = WaveModel(X, h, Ho, T, Seagr0, Marsh0, Mang0, 0) # vegetation absent
    Hv, Etav, Ubotv, VegLoc = WaveModel(X, h, Ho, T, Seagr, Marsh, Mang, 0) # vegetation present

    # percent wave attenuation
    o = num.nonzero(VegLoc <> 0);o = o[0]
    AtnH = num.arange(0, len(o), 1) * 0
    for xx in range(len(o)):
        pp = o[xx]
        AtnH[xx] = (H[pp] - Hv[pp]) / H[pp] * 100 # percent attenuation
    Atn = average(AtnH)
    VegLoc = VegLoc + nan;VegLoc[o] = -h[o] # location of vegetation
    Eta0 = max(Eta);Eta1 = max(Etav) # mean water level attenuation
        
    # save outputs
    WaveHeight = outputws + "VegetationWaveHeight_" + subwsStr + ".txt"
    file = open(WaveHeight, "w")
    for i in range(0, len(X)):
        file.writelines(str(X[i]) + "\t" + str(H[i]) + "\t" + str(Hv[i]) + "\n")
    file.close()

    # plot
    figure(1);
    subplot(311);plot(X, H, X, Hv);grid()
    ylabel('Wave Height [m]', size='large')
    legend(('No Vegetation', 'Vegetation Present'), 'upper left')
    subplot(312);plot(X, Ubot, X, Ubotv);grid()
    ylabel('Bot. Veloc. [m/s]', size='large')
    subplot(313);plot(X, VegLoc, 'xg', X, -h);grid()
    if len(o) > 0:
        legend(('Vegetation Field', ''), 'upper left')
    ylabel('Water Depth [m]', size='large')
    xlabel('Cross-Shore Distance from Offshore Boundary')
    savefig(outputws + "ErosionPlot1_" + subwsStr + ".png", dpi=(640 / 8))
    
    # estimate erosion in presence or absence of vegetation
    if hogr + hogt + hogc + hor == 0.0:
        gp.AddMessage('Estimate erosion amount...')
        A = cell1.Range("b91").Value;
        Slope = cell1.Range("f49").Value;m = 1.0 / Slope # foreshore slope=1/Slope
        B0 = cell1.Range("g49").Value
        W0 = cell1.Range("h49").Value
        D0 = cell1.Range("i49").Value

        # expected runup and max MWL
        Lo = g * T ** 2.0 / (2.0 * pi)
        Rmax = 1.1 * (0.35 * m * sqrt(Ho * Lo) + sqrt(Ho * Lo * (0.563 * m ** 2.0 + 0.004)) / 2.0) # max runup
        Emax = 0.35 * m * num.sqrt(Ho * Lo)

        # correct runup
        coef0 = max(Eta) / Emax;Etap = max(Etav) / coef0 # corrected MWL at shoreline in presence of veg.
        if Etap < 0:
            Etap = 0 # if Eta with veg. neg, take as zero

        Hp = (Etap / (0.35 * m)) ** 2 / Lo # Hprime to estimate runup with veg
        Rnp = 1.1 * (Etap + num.sqrt(Lo * (Hp * 0.563 * m ** 2 + 0.004 * Ho)) / 2) # runup with vegetation
                
        # erosion
        if d50 < 2 and d50 > 0.06:
            Rinf, R0 = ErosionKD(Ho, Rmax, B0, D0, W0)
            Rinf, Rv1 = ErosionKD(Hp, Rnp, B0, D0, W0)
            Rv2 = R0 * Rnp / Rmax # scale erosion distance when vegetation is present by ratio of runup
            Rveg = min(Rv1, Rv2)
        else:
            R0 = 0;Rveg = 0
            gp.AddMessage("Model can't compute erosion for the type of sediment that you have.")

        # beach profile plot
        if D0 == 0:
            W0 = 100
        Xp = num.arange(0, 1000, 1)
        Yp = m * Xp
        loc = num.nonzero(Yp >= B0);loc = loc[0]
        Yp[loc] = B0
        Lberm = loc[0];Ltoe = Indexed(Xp, Xp[Lberm + W0])
        Yp[Ltoe:-1] = B0 + D0

        Yp0 = m * Xp - m * R0 # profile of erosion w/o vegetation
        loc = num.nonzero(Yp0 >= B0);loc = loc[0]
        if R0 <= W0:
            Yp0[loc] = B0
            Lberm0 = loc[0]
            W00 = W0 - (Xp[Lberm0] - Xp[Lberm])
            Ltoe0 = Indexed(Xp, Xp[Lberm0 + W00])
            Yp0[Ltoe0:-1] = B0 + D0
        else:
            Yp0[loc] = B0 + D0
            Ltoe0 = loc[0]

        Ypv = m * Xp - m * Rveg # profile of erosion with vegetation
        loc = num.nonzero(Ypv >= B0);loc = loc[0]
        if Rveg <= W0:
            Ypv[loc] = B0
            Lbermv = loc[0]
            W0v = W0 - (Xp[Lbermv] - Xp[Lberm])
            Ltoev = Indexed(Xp, Xp[Lbermv + W0v])
            Ypv[Ltoev:-1] = B0 + D0
        else:
            Ypv[loc] = B0 + D0

        figure(2)
        subplot(211)
        plot(Xp, Yp, Xp, Yp0, '--r', Xp, Ypv, '--g');grid()
        xlim(0, Xp[Ltoe0] + R0 + 5);ylim(Yp[0], B0 + D0 + 1)
        ylabel('Backshore Elevation [m]', size='large', weight='demi')
        title('Erosion Profiles', size='large', weight='bold')
        legend(('Initial Profile', 'No Vegetation', 'Vegetation Present'), 'upper left')

        subplot(212)
        temp1 = Indexed(Xp, Xp[Lberm] - 10);temp2 = Indexed(Xp, Xp[Ltoe0] + 10)
        plot(Xp, Yp, Xp, Yp0, '--r', Xp, Ypv, '--g');grid()
        xlim(Xp[temp1], Xp[temp2] + R0);ylim(Yp[temp1], B0 + D0 + 0.5)
        ylabel('Backshore Elevation [m]', size='large', weight='demi')
        xlabel('Cross-Shore Distance [m]', size='large', weight='demi')
        savefig(outputws + "ErosionPlot2_" + subwsStr + ".png", dpi=(640 / 8))
        Fig2 = 1; #for HTML

elif  MangmtActn == 2: # backshore is modified
    gp.AddMessage("\nReading depth profile...")
    # read in user's cross-shore profile
    if CSProfile:
        TextData = open(CSProfile, "r") 
        X_lst = [];h_lst = []
        for line in TextData.read().strip("\n").split("\n"):
            linelist = [float(s) for s in line.split("\t")] # split the list by comma delimiter
            X_lst.append(linelist[0])
            h_lst.append(linelist[1])
        TextData.close()

        X = num.array(X_lst);h = num.array(h_lst);
        keep = num.nonzero(h < -0.5)
        h = -h[keep];X = X[keep];lx = len(X)
        if h[0] < h[-1]:
           h = h[::-1] # reverse order if profile starts at shoreline

    A = cell1.Range("b92").Value
    Slope = cell1.Range("f54").Value;m = 1.0 / Slope # foreshore slope=1/Slope
    B0 = cell1.Range("g54").Value # initial berm height
    W0 = cell1.Range("h54").Value # initial berm length
    D0 = cell1.Range("i54").Value # initial dune height
    B1 = cell1.Range("g55").Value # modified berm height
    W1 = cell1.Range("h55").Value # modified berm length
    D1 = cell1.Range("i55").Value # modified dune height

    gp.AddMessage("\nComputing wave profile...")
    Seagr0 = [[0, 0, 0], [0, 0]];Marsh0 = [[0, 0, 0], [0, 0]]
    Mang0 = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [0, 0]]
    H, Eta, Ubot, temp = WaveModel(X, h, Ho, T, Seagr0, Marsh0, Mang0, 0) # vegetation absent
    # save outputs
    WaveHeightExplicit = outputws + "BeachChangeWaveHeight_" + subwsStr + ".txt"
    file = open(WaveHeightExplicit, "w")
    for i in range(0, len(X)):
        file.writelines(str(X[i]) + "\t" + str(H[i]) + "\n")
    file.close()

    # plot
    figure(1)
    subplot(311);plot(X, H);grid()
    ylabel('Wave Height [m]', size='large')
    subplot(312);plot(X, Ubot);grid()
    ylabel('Bot. Veloc. [m/s]', size='large')
    subplot(313);plot(X, -h);grid()
    ylabel('Water Depth [m]', size='large')
    xlabel('Cross-Shore Distance from Offshore Boundary')
    savefig(outputws + "ErosionPlot1_" + subwsStr + ".png", dpi=(640 / 8))
    
    gp.AddMessage("\nCalculating amount of erosion...")
    Lo = g * T ** 2.0 / (2.0 * pi)
    Rmax = 1.1 * (0.35 * m * sqrt(Ho * Lo) + sqrt(Ho * Lo * (0.563 * m ** 2.0 + 0.004)) / 2.0) # max runup
    Rinf, R0 = ErosionKD(Ho, Rmax, B0, D0, W0)
    Rinf1, R1 = ErosionKD(Ho, Rmax, B1, D1, W1)

    # beach profile plot (original profile)
    if D0 == 0:
        W0 = 100
    Xp = num.arange(0, 1000, 1)
    Yp = m * Xp
    loc = num.nonzero(Yp >= B0);loc = loc[0]
    Yp[loc] = B0
    Lberm = loc[0];Ltoe = Indexed(Xp, Xp[Lberm + W0])
    Yp[Ltoe:-1] = B0 + D0

    Yp0 = m * Xp - m * R0 # profile of erosion 
    loc = num.nonzero(Yp0 >= B0);loc = loc[0]
    Lberm0 = loc[0]
    if R0 <= W0:
        Yp0[loc] = B0        
        W00 = W0 - (Xp[Lberm0] - Xp[Lberm])
        Ltoe0 = Indexed(Xp, Xp[Lberm0 + W00])
        Yp0[Ltoe0:-1] = B0 + D0
    else:
        Yp0[loc] = B0 + D0
        Ltoe0 = loc[0]

    # plot
    figure(2)    
    subplot(221)
    plot(Xp, Yp, Xp, Yp0, '--r');grid()
    legend(('Initial Profile', 'Eroded Profile'), 'upper left')
    xlim(0, Xp[Ltoe0] + 2);ylim(Yp[0], B0 + D0 + 0.5)
    ylabel('Backshore Elevation [m]', size='large', weight='demi')
    title('Erosion of Initial Profile', size='large', weight='bold')
    subplot(222)
    temp1 = Indexed(Xp, Xp[Lberm] - 2);temp2 = Indexed(Xp, Xp[Ltoe0] + 2)
    plot(Xp, Yp, Xp, Yp0, '--r',);grid()
    xlim(Xp[temp1], Xp[temp2]);ylim(Yp[temp1], B0 + D0 + 0.5)

    # beach profile plot (modified profile)
    if D1 == 0:
        W1 = 100
    Xp = num.arange(0, 1000, 1)
    Yp = m * Xp
    loc = num.nonzero(Yp >= B1);loc = loc[0]
    Yp[loc] = B1
    Lberm = loc[0];Ltoe = Indexed(Xp, Xp[Lberm + W1])
    Yp[Ltoe:-1] = B1 + D1

    Yp1 = m * Xp - m * R1 # profile of erosion
    loc = num.nonzero(Yp1 >= B1);loc = loc[0]
    Lberm1 = loc[0]
    if R1 <= W1:
        Yp1[loc] = B1
        W10 = W1 - (Xp[Lberm1] - Xp[Lberm])
        Ltoe1 = Indexed(Xp, Xp[Lberm1 + W10])
        Yp1[Ltoe1:-1] = B1 + D1
    else:
        Yp1[loc] = B1 + D1
        Ltoe1 = Indexed(Xp, Xp[Lberm1 + R0])

    figure(2)    
    subplot(223)
    plot(Xp, Yp, Xp, Yp1, '--r');grid()
    legend(('Initial Profile', 'Eroded Profile'), 'upper left')
    xlim(0, Xp[Ltoe1] + R0 + 2);ylim(Yp[0], B1 + D1 + .5)
    ylabel('Backshore Elevation [m]', size='large', weight='demi')
    xlabel('Cross-Shore Distance [m]', size='large', weight='demi')
    title('Erosion of Modified Profile', size='large', weight='bold')

    subplot(224)
    temp1 = Indexed(Xp, Xp[Lberm] - 2)
    temp2 = Indexed(Xp, Xp[Ltoe1] + 2)
    plot(Xp, Yp, Xp, Yp1, '--r', Xp);grid()
    xlim(Xp[temp1], Xp[temp2] + R0 + 2);ylim(Yp[temp1], B1 + D1 + 0.5)
    xlabel('Cross-Shore Distance [m]', size='large', weight='demi')    
    savefig(outputws + "ErosionPlot2_" + subwsStr + ".png", dpi=(640 / 8))
    Fig2 = 1;
    
elif MangmtActn == 3: # coral reef
    gp.AddMessage("\nCoral reef calculations...")
    AlphF = cell1.Range("d59").Value # reef face slope
    AlphR = cell1.Range("e59").Value	# reef rim slope
    he = cell1.Range("f59").Value # reef rim edge
    hr = cell1.Range("g59").Value # reef top depth
    Wr = cell1.Range("h59").Value # reef top width
    hl = cell1.Range("i59").Value # lagoon depth
    TanAlph = cell3.Range("a2:a202").Value;TanAlph = num.array(TanAlph)
    Lo = g * T ** 2.0 / (2.0 * pi);sig = 2.0 * pi / T

    # read in user's cross-shore profile
    X_built = num.arange(0, 10000.0, 1.0)
    h_built = AlphF * X_built - 100.0;xx = X_built[:]; hh = h_built[:] # reef face
    ind = Indexed(h_built, -he)
    h_built[ind:-1] = AlphR * X_built[ind:-1] + (h_built[ind] - AlphR * X_built[ind]) # reef rim
    o = num.nonzero(h_built > -hr);o = o[0];Xloc = X_built[o[0]]
    h_built[o] = -hr # reef top
    o = num.nonzero(X_built <= Xloc + Wr);o = o[0]
    X_built = X_built[o];h_built = h_built[o];lx = len(X_built)
                  
    if CSProfile:
        TextData = open(CSProfile, "r") 
        X_lst = [];h_lst = []
        for line in TextData.read().strip("\n").split("\n"):
            linelist = [float(s) for s in line.split("\t")] # split the list by comma delimiter
            X_lst.append(linelist[0])
            h_lst.append(linelist[1])
        TextData.close()

        X = num.array(X_lst);h = num.array(h_lst);
        keep = num.nonzero(h < -0.5)
        h = -h[keep];X = X[keep];lx = len(X)
        if h[0] < h[-1]:
           h = h[::-1] # reverse order if profile starts at shoreline
    else:
        h = h_built;X = X_built
        
    H, Eta, H_r, Etar = WavesCoral()

    # reconfigure outputs so we can plot them
    Eta_r = [Eta[ii] * 0 + Etar for ii in range(len(Eta))]
    X = [ii for ii in range(len(H))];Xr = [ii * Wr / len(H_r) for ii in range(len(H_r))]

    figure(2)
    subplot(311);plot(X, H, Xr, H_r);grid()
    title('Coral Reef', size='large', weight='bold')
    ylabel('Wave Height [m]', size='large', weight='demi')
    ## legend(('Explicit Solution','Gourlay''s Method'),'lower right')
    savefig(outputws + "ErosionPlot1_" + subwsStr + ".png", dpi=(640 / 8))
    subplot(312);plot(X, Eta, X, Eta_r);grid()
    ylabel('Water Level on Reef[m]', size='large', weight='demi')
    subplot(313);plot(X, -h, X_built, h_built);grid()
    ylabel('Reef Profile')
    xlabel('Cross-Shore Distance from Offshore Boundary [m]', size='large', weight='demi')    
    savefig(outputws + "ErosionPlot1_" + subwsStr + ".png", dpi=(640 / 8))

    # save outputs
    WaveHeightExplicit = outputws + "CoralWaveHeight_" + subwsStr + ".txt"
    file = open(WaveHeightExplicit, "w")
    for i in range(0, len(X)):
        file.writelines(str(X[i]) + "\t" + str(H[i]) + "\n")
    file.close()
    
    CoralDepth = outputws + "CoralDepth_" + subwsStr + ".txt"
    file = open(CoralDepth, "w")
    for i in range(0, len(X_built)):
        file.writelines(str(X_built[i]) + "\t" + str(h_built[i]) + "\n")
    file.close()

    # erosion at beach
    if hl <> 0: # there is a lagoon
        Hreef = H_r[-1] # wave height in lagoon
    else:
        Hreef = H_r[0]

    Er_chk = cell1.Range("b89").Value
    if Er_chk == 1:
        A = cell1.Range("b93").Value
        Slope = cell1.Range("f65").Value;m = 1.0 / Slope # foreshore slope=1/Slope
        B0 = cell1.Range("g65").Value;D0 = 0 # berm height

        if hl <> 0: # there is a lagoon
            Lo = g * T ** 2 / (2 * pi);Cgo = 0.5 * Lo / T
            Rnp = 1.1 * (0.35 * m * sqrt(Hreef * Lo) + sqrt(0.563 * Hreef * Lo * m ** 2.0 + 0.004 * Ho * Lo) / 2.0) # max runup
        else: # fringing reef
            ksir = m / num.sqrt(H_r[0] / (g * T ** 2 / (2 * pi)))
            Rnp = 1.1 * (0.35 * m * sqrt(Hreef * Lo) + sqrt(0.563 * Hreef * Lo * m ** 2.0 + 0.004 * Ho * Lo) / 2.0) # max runup

        # expected runup and max MWL
        Rmax = 1.1 * (0.35 * m * sqrt(Ho * Lo) + sqrt(Ho * Lo * (0.563 * m ** 2.0 + 0.004)) / 2.0) # max runup
        # erosion
        Rinf, R0 = ErosionKD(Ho, Rmax, B0, 0, 0)
        Rinf, R1 = ErosionKD(Hreef, Rnp, B0, 0, 0)

        # beach profile plot
        D0 = 0;W0 = 100
        Xp = num.arange(0, 1000, 1)
        Yp = m * Xp
        loc = num.nonzero(Yp >= B0);loc = loc[0]
        Yp[loc] = B0;Lberm = Xp[loc[0]]

        Yp0 = m * Xp - m * R0 # profile of erosion w/t coral
        loc = num.nonzero(Yp0 >= B0);loc = loc[0]
        Yp0[loc] = B0

        Yp1 = m * Xp - m * R1 # profile of erosion with coral
        loc = num.nonzero(Yp1 >= B0);loc = loc[0]
        Yp1[loc] = B0

        figure(3);subplot(221)
        plot(Xp, Yp, Xp, Yp0, '--r', Xp, Yp1, '--g');grid()
        xlim(0, Xp[Lberm] + R0 + 10);ylim(Yp[0] - .5, B0 + 0.5)
        ylabel('Backshore Elevation [m]', size='large', weight='demi')
        xlabel('Cross-Shore Distance [m]', size='large', weight='demi')
        ## legend(('Initial Profile','Erosion Profile in Absence of Coral','Erosion Profile in Presence of Coral'),'upper left')
        savefig(outputws + "ErosionPlot2_" + subwsStr + ".png", dpi=(640 / 8))
        Fig2 = 1

elif MangmtActn == 4: # it is an oyster reef
    gp.AddMessage("\nPerforming oyster reef calculations...")
    TextData = open(CSProfile, "r") 
    X_lst = [];h_lst = [];
    for line in TextData.read().strip("\n").split("\n"):
        linelist = [float(s) for s in line.split("\t")] # split the list by comma delimiter
        X_lst.append(linelist[0])
        h_lst.append(linelist[1])
    TextData.close()     
    h_lst = [-h_lst[ii] for ii in range(len(h_lst))]
    if h_lst[0] < h_lst[-1]:
        h_lst = h_lst[::-1] # reverse order if profile starts at shoreline
    X = num.array(X_lst);h = num.array(h_lst);lx = len(X)

    # get data from Excel
    Xr = cell1.Range("d69").Value # distance from shoreline
    hc = cell1.Range("e69").Value # reef height
    Bw = cell1.Range("f69").Value # base width
    Cw = cell1.Range("g69").Value # crest width
   
    # wave height w/t reef
    Seagr0 = [[0, 0, 0], [0, 0]]
    Marsh0 = [[0, 0, 0], [0, 0]]
    Mang0 = [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [0, 0]]
    h = h + S # add surge level
    out = num.nonzero(h > 0.5);out = out[0]
    h = h[out];X = X[out] # only keep values that are below water

    # locate reef 
    h_lst = [-h[i] for i in range(len(h))]
    h_ar = num.array(h_lst) # water depth
    Xval = X[-1] - Xr;Xloc = Indexed(X, Xval) # location of reef

    dx = X[1] - X[0]
    H, Eta, Ubot, temp = WaveModel(X, h, Ho, T, Seagr0, Marsh0, Mang0, 0) # wave Profile
    hi = np.abs(h_ar[Xloc])
    
    Hi = H[Xloc] # wave height and water depth at edge of reef
    # compute transmission coefficient
    Kt = BreakwaterKt(Hi, T, hi, hc, Cw, Bw)
    if Bw - Cw > 1:
        KtB1 = BreakwaterKt(Hi, T, hi, hc, Cw, Bw + 1);KtB0 = BreakwaterKt(Hi, T, hi, hc, Cw, Bw - 1) # base width+-1
        KtC1 = BreakwaterKt(Hi, T, hi, hc, Cw + 1, Bw);KtC0 = BreakwaterKt(Hi, T, hi, hc, Cw - 1, Bw) # crest width+-1
    else:
        KtB1 = 'Bw-Cw=1';KtB0 = 'Bw-Cw=1';KtC1 = 'Bw-Cw=1';KtC0 = 'Bw-Cw=1'
    Kth1 = BreakwaterKt(Hi, T, hi, hc + 1, Cw, Bw);Kth0 = BreakwaterKt(Hi, T, hi, hc - 1, Cw, Bw) # reef height+-1
    # compute wave height behind breakwater
    Htemp, temp, Utemp, temp = WaveModel(X[Xloc:len(h)], h[Xloc:len(h)], H[Xloc] * Kt, T, Seagr0, Marsh0, Mang0, 0)
    H1_lst = [H[i] for i in range(len(h))];U1 = [Ubot[ii] for ii in range(len(h))]
    H1_lst[Xloc:len(h)] = Htemp;U1[Xloc:len(h)] = Utemp
    # plotting
    Yrf = num.arange((-hi), 0.0, 0.05);Xrf = Yrf * 0.0 + Xval # x-axis
    figure(1)
    subplot(211);
    plot(X, H, X, H1_lst);grid()
    ylabel('Wave Height [m]', size='large', weight='demi')
    legend(('Oyster Reef Absent', 'Oyster Reef Present'), 'upper left')
    subplot(212);
    plot(X, -h, Xrf, Yrf, 'r', Xrf + .05, Yrf, 'r', Xrf + .1, Yrf, 'r');grid()
    ylabel('Water Depth [m]', size='large', weight='demi')
    xlabel('Cross-Shore Distance from Offshore Boundary [m]', size='large', weight='demi')
    legend(('Depth Profile', 'Oyster Reef'), 'upper left')
    savefig(outputws + "ErosionPlot1_" + subwsStr + ".png", dpi=(640 / 8))

    figure(2)
    subplot(311);
    plot(X[Xloc - 50:len(h)], H[Xloc - 50:len(h)], X[Xloc - 50:len(h)], H1_lst[Xloc - 50:len(h)]);grid()
    ylabel('Wave Height', size='large', weight='demi')
    legend(('Oyster Reef Absent', 'Oyster Reef Present'), 'upper right')
    subplot(312);
    plot(X[Xloc - 50:len(h)], Ubot[Xloc - 50:len(h)], X[Xloc - 50:len(h)], U1[Xloc - 50:len(h)]);grid()
    ylabel('Bed Velo. [m/s]', size='large', weight='demi')
    subplot(313);
    plot(X[Xloc - 50:len(h)], -h[Xloc - 50:len(h)], Xrf, Yrf, 'r', Xrf + .05, Yrf, 'r', Xrf + .1, Yrf, 'r');grid()
    ylabel('Depth', size='large', weight='demi')
    xlabel('Cross-Shore Distance from Offshore Boundary [m]', size='large', weight='demi')
    legend(('Depth Profile', 'Oyster Reef'), 'upper left')
    savefig(outputws + "ErosionPlot2_" + subwsStr + ".png", dpi=(640 / 8))

    # save profile
    WaveHeightOyster = outputws + "OysterWaveHeight_" + subwsStr + ".txt"
    file = open(WaveHeightOyster, 'w')
    for ii in range(0, len(X)):
        file.writelines(str(X[ii]) + "\t" + str(H[ii]) + "\t" + str(H1_lst[ii]) + "\n")
    file.close()

    # output results in Excel
    cell4 = xlApp.Worksheets("OysterReefsBreakwater") # will go in HTML
    cell4.Range("f2").Value = hi # water depth
    cell4.Range("f8").Value = Hi # incident wave height
    cell4.Range("f10").Value = Kt # transmission coef
    cell4.Range("e15").Value = KtB1
    cell4.Range("e17").Value = KtB0
    cell4.Range("j15").Value = KtC1
    cell4.Range("j17").Value = KtC0
    cell4.Range("o15").Value = Kth1
    cell4.Range("o17").Value = Kth0

xlApp.ActiveWorkbook.Close(SaveChanges=1)
xlApp.Quit()

# create html file
gp.AddMessage("\nCreating outputs...")
htmlfile = open(outputws + "OutputWaveModel_" + subwsStr + ".html", "w")
htmlfile.write("<html>\n")
htmlfile.write("<title>Marine InVEST - Wave, Erosion, and Inundation</title>")
htmlfile.write("<CENTER><H1>Coastal Protection - Tier 1</H1><H2>Erosion Results<br></H2></CENTER>")
htmlfile.write("<br><HR><H2>Plot</H2>\n")

if MangmtActn == 1:
    if Fig2:
        htmlfile.write("This plot shows how vegetation affects waves height and beach erosion. <br>\n")
        htmlfile.write("<img src=\"ErosionPlot1_" + subwsStr + ".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")
        htmlfile.write("</td><td>")
        htmlfile.write("<img src=\"ErosionPlot2_" + subwsStr + ".png\" alt=\"Profile Plot #2\" width=\"640\" height=\"480\">")
        htmlfile.write("</td><td></td></tr></table><br>\n")
        htmlfile.write("<br><HR><H2>More Information</H2>\n")
        if S + Rmax > B0 + D0:
            htmlfile.write("<li>Without vegetation, the water is: " + str(round(S + Rmax, 1)) + "m above the highest point in the backshore<br>\n")
            htmlfile.write("<li>Without vegetation, the beach inundated, and it is going to erode by more than: " + str(round(R0, 1)) + "m <br>\n")
        else:
            htmlfile.write("<li>Without vegetation, the beach is going to erode by: " + str(round(R0, 1)) + "m <br>\n")
        if S + Rnp > B0 + D0:
            htmlfile.write("<li>With vegetation present, the water is: " + str(round(S + Rnp, 1)) + "m above the highest point in the backshore<br>\n")
            htmlfile.write("<li>With vegetation present, the beach inundated, and it is going to erode by more than: " + str(round(Rveg, 1)) + "m <br>\n")
        else:
            htmlfile.write("<li>With vegetation present, the beach is going to erode by: " + str(round(Rveg, 1)) + "m <br>\n")
    else:
        htmlfile.write("This plot shows how vegetation affects wave height. <br>\n")
        htmlfile.write("<img src=\"ErosionPlot1.png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")
        htmlfile.write("</td><td></td></tr></table><br>\n")
        htmlfile.write("<br><HR><H2>More Information</H2>\n")
        
    htmlfile.write("<li>Average wave attenuation in presence of vegetation is: " + str(round(Atn)) + "% <br>\n")
    
if MangmtActn == 2:
    htmlfile.write("This plot shows how changes in backshore configuration affect beach erosion. <br>\n")
    htmlfile.write("<img src=\"ErosionPlot1_" + subwsStr + ".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")
    htmlfile.write("</td><td>")
    htmlfile.write("<img src=\"ErosionPlot2_" + subwsStr + ".png\" alt=\"Profile Plot #2\" width=\"640\" height=\"480\">")
    htmlfile.write("</td><td></td></tr></table><br>\n")
    htmlfile.write("<br><HR><H2>More Information</H2>\n")
    if S + Rmax > B0 + D0:
        htmlfile.write("<li>Based on your initial profile characteristics, the water is: " + str(round(S + Rmax, 1)) + "m above the highest point in the backshore<br>\n")
        htmlfile.write("<li>The beach inundated, and it is going to erode by more than: " + str(round(R0, 1)) + "m <br>\n")
    else:
        htmlfile.write("<li>Based on your initial profile characteristics, the beach is going to erode by: " + str(round(R0, 1)) + "m <br>\n")
    if S + Rmax > B1 + D1:
        htmlfile.write("<li>Based on your modified profile characteristics, the water is: " + str(round(S + Rmax, 1)) + "m above the highest point in the backshore<br>\n")
        htmlfile.write("<li>The beach inundated, and it is going to erode by more than: " + str(round(R1, 1)) + "m <br>\n")
    else:
        htmlfile.write("<li>Based on your modified profile characteristics, the beach is going to erode by: " + str(round(R1, 1)) + "m <br>\n")

if MangmtActn == 3:
    if Fig2:
        htmlfile.write("This plot shows how coral reefs affect waves height and beach erosion. <br>\n")
        htmlfile.write("<img src=\"ErosionPlot1_" + subwsStr + ".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")
        htmlfile.write("</td><td>")
        htmlfile.write("<img src=\"ErosionPlot2_" + subwsStr + ".png\" alt=\"Profile Plot #2\" width=\"640\" height=\"480\">")
        htmlfile.write("</td><td></td></tr></table><br>\n")
        if S + Rmax > B0:
            htmlfile.write("<li>Beach inundated, without the reef, the beach is going to erode more than: " + str(round(R0, 1)) + "m <br>\n")
        else:
            htmlfile.write("<li>Without the reef, the shoreline will erode by: " + str(round(R0, 1)) + "m <br>\n")
        if S + Rnp > B0:
            htmlfile.write("<li>Beach inundated, with the reef, the beach is going to erode more than: " + str(round(R1, 1)) + "m <br>\n")
        else:
            htmlfile.write("<li>With the reef present, the shoreline will erode by: " + str(round(R1, 1)) + "m <br>\n")
        htmlfile.write("<br><HR><H2>More Information</H2>\n")    
    else:
        htmlfile.write("This plot shows how coral reefs affect waves height. <br>\n")
        htmlfile.write("<img src=\"ErosionPlot1_" + subwsStr + ".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")
        htmlfile.write("</td><td></td></tr></table><br>\n")
        htmlfile.write("<br><HR><H2>More Information</H2>\n")
    htmlfile.write("<li>Wave height on the reef is: " + str(round(Hreef, 2)) + "<br>\n")
    htmlfile.write("<li>Percent wave attenuation is: " + str(round(Hreef / Ho * 100.0, 1)) + "<br>\n")
    htmlfile.write("<li>Water level at shoreline in absence of reef is: " + str(round(S + Rmax, 1)) + "<br>\n")
    htmlfile.write("<li>Percent water level attenuation is: " + str(round((S + Rnp) / (S + Rmax) * 100.0, 1)) + "% <br>\n")
    htmlfile.write("<li>Runup in absence of reef is: " + str(round(Rmax, 1)) + "<br>\n")
    htmlfile.write("<li>Percent runup attenuation is: " + str(round(Rnp / Rmax * 100.0, 1)) + "% <br>\n")

if MangmtActn == 4:
    htmlfile.write("This plot shows how oyster reefs affect wave height. <br>\n")
    htmlfile.write("<img src=\"ErosionPlot1_" + subwsStr + ".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")
    htmlfile.write("<img src=\"ErosionPlot2_" + subwsStr + ".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")
    htmlfile.write("</td><td></td></tr></table><br>\n")

    htmlfile.write("<br><HR><H2>More Information</H2>\n")
    htmlfile.write("<li>Transmission coefficient is: " + str(round(Kt * 100.0, 2)) + "% <br>\n")

htmlfile.write("<li>Inputs conditions are: Ho=" + str(round(Ho, 2)) + ", T=" + str(round(T, 1)) + "<br>\n")
htmlfile.write("</html>")
htmlfile.close()

# create parameter file
parameters.append("Script location: " + os.path.dirname(sys.argv[0]) + "\\" + os.path.basename(sys.argv[0]))
parafile = open(outputws + "parameters_" + now.strftime("%Y-%m-%d-%H-%M") + ".txt", "w") 
parafile.writelines("WAVE / EROSION MODEL PARAMETERS\n")
parafile.writelines("_______________________________\n\n")
for para in parameters:
    parafile.writelines(para + "\n")
    parafile.writelines("\n")
parafile.close()
