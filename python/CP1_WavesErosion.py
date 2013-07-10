# Marine InVEST: Coastal Protection (Wave and Erosion)
# Authors: Greg Guannel, Gregg Verutes, Joe Faries, Apollo Xi
# Coded for ArcGIS 9.3, 10, 10.1
# 10/19/12

import sys, os, string, time, datetime
import CPf_SignalSmooth as SignalSmooth
from math import *
import fpformat, operator
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
msgReadExcel = "\nError reading Erosion Protection Excel file inputs."
msgTwoHabitats = "\nYou cannot have marsh and mangrove at your site. Please remove one of the two."
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
		parameters.append("Erosion Protection Excel Table: "+InputTable)
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
		EconCP=gp.GetParameterAsText(13)
		parameters.append("Compute Economic Valuation? "+EconCP)
		Longshore=gp.GetParameterAsText(14)
		parameters.append("Longshore Extent: "+str(Longshore))
		if Longshore:
			Longshore=float(gp.GetParameterAsText(14))
		PropValue=gp.GetParameterAsText(15)
		parameters.append("Property Value: "+str(PropValue))
		if PropValue:
			PropValue=float(gp.GetParameterAsText(15))
		Tr=gp.GetParameterAsText(16)
		parameters.append("Return Period: "+str(Tr))
		if Tr:
			Tr=float(gp.GetParameterAsText(16))
		disc=gp.GetParameterAsText(17)
		parameters.append("Discount Rate: "+str(disc))
		if disc:
			disc=float(gp.GetParameterAsText(17))
		TimeHoriz=gp.GetParameterAsText(18)
		parameters.append("Time Horizon: "+str(TimeHoriz))
		if TimeHoriz:
			TimeHoriz=float(gp.GetParameterAsText(18))
		parameters=[]
		now=datetime.datetime.now()
	except:
		raise Exception,msgArguments+gp.GetMessages(2)


	###############################################
	###### CREATE DIRECTORIES AND VARIABLES #######
	###############################################

	# remove spaces and shorten 'subwsStr' if greater than 10 characters
	subwsStr=subwsStr.replace(" ","")
	subwsStr=subwsStr[0:10]

	#intermediate and output directories
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
		Answer1=[Ho,To] # if user provides offshore height and period
		Answer2=[Us,Ft,depth] # if user would like starting wave conditions to be calculated from wind speed, fetch and depth
		if WaveErosionQuestion == "(1) Yes, I have these values":
			for Xaxis in Answer1:
				if Xaxis == '':
					gp.AddError("\nOne or more of the required input parameters was not defined.\nPlease provide values for either Wave Height and/or Wave Period.")
					raise Exception
		else:
			for x2 in Answer2:
				if x2 == '':
					gp.AddError("\nOne or more of the required input parameters was not defined.\nPlease provide values for either Wind Speed, Fetch Distance, and/or Water Depth.")
					raise Exception

		# check that all economic valuation inputs have been provided
		if EconCP=='true':
			econInputs = [Longshore, PropValue, TimeHoriz, Tr, disc]
			for ii in range(0,len(econInputs)):
				if econInputs[ii] == '':
					gp.AddError("\nOne or more of the required economic valuation input parameters was not defined.  We cannot compute a value for your habitats.\n")
					raise Exception
			valua=1
	except:
		gp.AddError(msgCheckInputs)
		raise Exception

	def gradient(f,z): # calculates the differences between consecutive entries in a list (f) scaled by a factor (z) 
		length=len(f)
		df=length*[0.0]
		df[0]=(f[1]-f[0])/z
		for i in range(1,length-1):
			df[i]=(f[i+1]-f[i-1])/(2.0*z)
		df[length-1]=(f[length-1]-f[length-2])/z
		return df

	def Indexed(x,value): # locates index of point in vector x (type = array) that has closest value as variable value
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

	# calculates the wave number given an angular frequency (sigma) and local depth (dh)
	def iterativek(sigma,dh):
		g=9.81;rho=1024;
		kestimated=(sigma**2)/(g*(sqrt(tanh((sigma**2)*dh/g)))) # intial guess for the wave number
		tol = 1 # initialize the tolerence that will determine if the iteration should continue
		count=0
		while (tol>0.00005) and (count<1000): # iterate over the dispersion relationship
			count += 1
			kh=kestimated*dh
			tol = kestimated-(sigma**2)/(tanh(kh)*g) # check the difference between the previous and the current iterations
			kcalculated=(sigma**2)/(tanh(kh)*g) # k values for current iteration
			kestimated=kcalculated # set the current value for k as the previous value for the subsequent iteration
		qk=kcalculated # return the value of k when accuracy tolerance is achieved
		return qk

	# wind-wave generation
	def WindWave(U,F,d):
		# calculates wind waves which are a function of Wind Speed (U), Fetch length (F), and local depth (d) using empirical equations
		g=9.81;rho=1024;
		ds=g*d/U**2.0;Fs=g*F/U**2.0
		A=tanh(0.343*ds**1.14)
		B=tanh(4.14e-4*Fs**0.79/A)
		H=0.24*U**2/g*(A*B)**0.572 # wave height
		A=tanh(0.1*ds**2.01)
		B=tanh(2.77e-7*Fs**1.45/A)
		T=7.69*U/g*(A*B)**0.18 # wave period
		return H,T

	# wave transformation model
	def WaveModel(X,h,Ho,To,Etao,Roots,Trunk,Canop,VegXloc,TrigRange):
		# X: Vector of consecutive cross-shore positions going shoreward
		# n: Vector of water depths at the corresponding cross-shore position
		# Ho: Initial wave height to be applied at the first cross-shore position
		# To: Wave Period
		# Etao: Mean Water Level increase due wave set-up at the first cross-shore position
		# Roots: An array of physical properties (density, diameter, and height) of mangrove roots at the corresponding cross-shore position (all zeros if mangroves don't exist at that location)
		# Trunk: An array of physical properties (density, diameter, and height) of mangrove trunks or the marsh or seagrass plants at the corresponding cross-shore position (all zeros if vegetation don't exist at that location)
		# Canop: An array of physical properties (density, diameter, and height) of the mangrove canopy at the corresponding cross-shore position (all zeros if mangroves don't exist at that location)
		# VegXloc: A vector with a numeric code indicate what (if any) natural habitat is present at the cross-shore location.  0 = No Habitat, 1 = Mangrove, 2 = Marsh, 3 = Seagrass, 4 = Coral, 5 = Oyster Reef
		# TrigRange: Defines the segment of the cross-shore domain where a given model is valid.  Where VegXloc is 0, 1, 2, or 3 the same wave model is applicable.  If a coral or oyster reef is present the model is interupted, the reef model is run and that output is carried onto the next segment of the cross-shore domain.

		# constants
		g=9.81;rho=1024.0;B=1.0;Beta=0.05;Cf=0.01

		# extract vegetation characteristics
		NRoots=Roots[0];dRoots=Roots[1];hRoots=Roots[2]
		NTrunk=Trunk[0];dTrunk=Trunk[1];hTrunk=Trunk[2]
		NCanop=Canop[0];dCanop=Canop[1];hCanop=Canop[2]

		# input characteristics in segment of interest
		NRoots=NRoots[TrigRange];dRoots=dRoots[TrigRange];hRoots=hRoots[TrigRange]
		NTrunk=NTrunk[TrigRange];dTrunk=dTrunk[TrigRange];hTrunk=hTrunk[TrigRange]
		NCanop=NCanop[TrigRange];dCanop=dCanop[TrigRange];hCanop=hCanop[TrigRange]
		Vegxloc=VegXloc[TrigRange];
		lx=len(X);dx=abs(X[1]-X[0])

		# create relative depth values for roots, trunk and canopy
		alphr=hRoots/h;alpht=hTrunk/h;alphc=hCanop/h
		for kk in range(lx): 
			if alphr[kk]>1:
				alphr[kk]=1;alpht[kk]=0;alphc[kk]=0 # roots only
			elif alphr[kk]+alpht[kk]>1:
				alpht[kk]=1-alphr[kk];alphc[kk]=0 # roots and trunk
			elif alphr[kk]+alpht[kk]+alphc[kk]>1:
				alphc[kk]=1-alphr[kk]-alpht[kk] # roots, trunk and canopy

		# drag coefficent for vegetation; mangrove and marsh win over seagrass if they overlap
		CdVeg=num.arange(0.0,lx,dx)*0.0;
		CdVeg[Vegxloc==1]=1 # mangrove        
		CdVeg[Vegxloc==2]=0.1 # marsh        
		CdVeg[Vegxloc==3]=0.1 # seagrass
		CdVeg=SignalSmooth.smooth(num.array(CdVeg),len(CdVeg)*0.01,'hanning') 

		# initialize vectors for wave model
		H=lx*[0.0]; # RMS Wave Height
		Db=lx*[0.0];Df=lx*[0.0];Dveg=lx*[0.0] 
		k=lx*[0.0];L=lx*[0.0] # wave number; wave length
		C=lx*[0.0];n=lx*[0.0];Cg=lx*[0.0] # wave celerity; shoaling factor group velocity (C*n)
		Er=lx*[0.0]; Ef=lx*[0.0]; Br=lx*[0.0] # roller energy; energy flux; roller flux 
		Hs=lx*[0.0]; Etas=lx*[0.0] # wave height; setup in the absence of vegetation
		Dbs=lx*[0.0]; Dfs=lx*[0.0]; Ers=lx*[0.0] # dissipation due to breaking; dissipation due to bottom friction; roller energy
		Efs=lx*[0.0]; Brs=lx*[0.0] # energy flux in the absence of veg; roller flux in the absence of vegetation

		# wave parameter at 1st grid pt
		ash=[h[ii] for ii in range(lx)] # ash is same as h, but is now an independent variable
		fp=1.0/To; sig=2.0*pi*fp # wave frequency and angular frequency
		k[0]=iterativek(sig,h[0]) # wave number at 1st grid pt
		L[0]=2.0*pi/k[0] # wave length @ 1st grid pt
		n[0]=0.5*(1+(2.0*k[0]*h[0]/sinh(2.0*k[0]*h[0]))) # to compute Cg at 1st grid pt
		C[0]=L[0]/To;Cg[0]=C[0]*n[0] # phase and group velocity at 1st grid pt
		So=Ho/L[0] # deep water wave steepness
		Gam=0.5+0.4*num.tanh(33.0*So) # Gam from Battjes and Stive 85 (as per Alsina & Baldock)

		# RMS wave height at first grid point; assume no dissipation occurs
		H[0]=Ho/sqrt(2.0) # transform significant wave height to rms wave height         
		Ef[0]=0.125*rho*g*H[0]**2*Cg[0];Efs[0]=Ef[0];Hs[0]=H[0] # energy flux @ 1st grid pt

		# begin wave model 
		for xx in range(lx-1): # transform waves,take MWL into account
			# wave in presence of veg.
			Ef[xx]=0.125*rho*g*(H[xx]**2.0)*Cg[xx] # Ef at (xx)      
			Ef[xx+1]=Ef[xx]-dx*(Db[xx]+Df[xx]+Dveg[xx]) # Ef at [xx+1] (subtract dissipation due to: breaking, bottom friction, and vegetation)

			# phase and group velocity
			k[xx+1]=iterativek(sig,h[xx+1]) # compute wave number
			n[xx+1]=0.5*(1.0+(2.0*k[xx+1]*h[xx+1]/sinh(2.0*k[xx+1]*h[xx+1])))
			C[xx+1]=sig/k[xx+1];Cg[xx+1]=C[xx+1]*n[xx+1] # phase and group velocity

			# roller info
			H[xx+1]=num.sqrt(8.0*Ef[xx+1]/(rho*g*Cg[xx+1])) # wave height at [xx+1]      
			Br[xx+1]=Br[xx]-dx*(g*Er[xx]*sin(Beta)/C[xx]-0.5*Db[xx]) # roller flux
			Er[xx+1]=Br[xx+1]/(C[xx+1]) # roller energy

			# dissipation due to brkg            
			Var=0.25*rho*g*fp*B
			Hb=0.88/k[xx+1]*tanh(Gam*k[xx+1]*h[xx+1]/0.88)                
			temp1=((Hb/H[xx+1])**3.0+1.5*Hb/H[xx+1])*num.exp(-(Hb/H[xx+1])**2.0)
			if isreal(H[xx+1]): temp2=0.75*num.sqrt(pi)*(1-erf(Hb/H[xx+1]))
			else: temp2=0
			Db[xx+1]=Var*H[xx+1]**3/h[xx+1]*(temp1+temp2) # dissipation due to brkg

			# dissipation due to bot friction 
			Df[xx+1]=rho*Cf/(12.0*pi)*(2.0*pi*fp*H[xx+1]/sinh(k[xx+1]*h[xx+1]))**3.0 

			# dissipation due to vegetation
			V1=3*sinh(k[xx+1]*alphr[xx+1]*h[xx+1])+sinh(k[xx+1]*alphr[xx+1]*h[xx+1])**3 # roots
			V2=(3*sinh(k[xx+1]*(alphr[xx+1]+alpht[xx+1])*h[xx+1])-3*sinh(k[xx+1]*alphr[xx+1]*h[xx+1])+
				sinh(k[xx+1]*(alphr[xx+1]+alpht[xx+1])*h[xx+1])**3-
				sinh(k[xx+1]*alphr[xx+1]*h[xx+1])**3) # trunk
			V3=(3*sinh(k[xx+1]*(alphr[xx+1]+alpht[xx+1]+alphc[xx+1])*h[xx+1])
				-3*sinh(k[xx+1]*(alphr[xx+1]+alpht[xx+1])*h[xx+1])+
				sinh(k[xx+1]*(alphr[xx+1]+alpht[xx+1]+alphc[xx+1])*h[xx+1])**3-
				sinh(k[xx+1]*(alphr[xx+1]+alpht[xx+1])*h[xx+1])**3) # canopy

			CdDN=CdVeg[xx+1]*(dRoots[xx+1]*NRoots[xx+1]*V1+dTrunk[xx+1]*NTrunk[xx+1]*V2+dCanop[xx+1]*NCanop[xx+1]*V3)
			temp1=rho*CdDN*(k[xx+1]*g/(2.0*sig))**3.0/(2.0*sqrt(pi))
			temp3=(3.0*k[xx+1]*cosh(k[xx+1]*h[xx+1])**3)
			Dveg[xx+1]=temp1/temp3*H[xx+1]**3 # dissipation due to vegetation

			# wave in absence of vegetation
			Hs[xx+1]=H[xx+1]
			if sum(Vegxloc)<>0: # compute if there's vegetation (time saver)
				Efs[xx]=0.125*rho*g*(Hs[xx]**2.0)*Cg[xx] # Ef at (xx)      
				Efs[xx+1]=Efs[xx]-dx*(Dbs[xx]+Dfs[xx]) # Ef at [xx+1] 

				Hs[xx+1]=num.sqrt(8.0*Efs[xx+1]/(rho*g*Cg[xx+1])) # wave height at [xx+1]      
				Brs[xx+1]=Brs[xx]-dx*(g*Ers[xx]*sin(Beta)/C[xx]-0.5*Dbs[xx]) # roller flux
				Ers[xx+1]=Brs[xx+1]/(C[xx+1]) # roller energy

				temp1=((Hb/Hs[xx+1])**3.0+1.5*Hb/Hs[xx+1])*num.exp(-(Hb/Hs[xx+1])**2.0)
				if isreal(Hs[xx+1]):  temp2=0.75*num.sqrt(pi)*(1-erf(Hb/Hs[xx+1]))
				else: temp2=0
				Dbs[xx+1]=Var*Hs[xx+1]**3/h[xx+1]*(temp1+temp2) # dissipation due to brkg

				Dfs[xx+1]=rho*Cf/(12.0*pi)*(2.0*pi*fp*Hs[xx+1]/sinh(k[xx+1]*h[xx+1]))**3.0 # dissipation due to bottom friction 

		Ew=lx*[0.0];Ew=[0.125*rho*g*(H[i]**2.0) for i in range(lx)] # energy density
		Ews=lx*[0.0];Ews=[0.125*rho*g*(Hs[i]**2.0) for i in range(lx)] # energy density in the absence of vegetation

		# force on plants if they were emergent; take a portion if plants occupy only portion of wc
		Fxgr=[rho*g*CdVeg[i]*dRoots[i]*NRoots[i]*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
		Fxgt=[rho*g*CdVeg[i]*dTrunk[i]*NTrunk[i]*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
		Fxgc=[rho*g*CdVeg[i]*dCanop[i]*NCanop[i]*H[i]**3.0*k[i]/(12.0*pi*tanh(k[i]*ash[i])) for i in range(lx)]
		fx=[-alphr[i]*Fxgr[i]-alpht[i]*Fxgt[i]-alphc[i]*Fxgc[i] for i in range(lx)] # scale by height of indiv. elements

		# estimate MWS
		dx=1;Xi=num.arange(X[0],X[-1]+dx,dx);lxi=len(Xi) # use smaller dx to get smoother result
		F=interp1d(X,ash);ashi=F(Xi);
		F=interp1d(X,k);ki=F(Xi);
		F=interp1d(X,Ew);Ewi=F(Xi);
		F=interp1d(X,Er);Eri=F(Xi);
		F=interp1d(X,fx);fxi=F(Xi);
		Sxx=lxi*[0.0];Rxx=lxi*[0.0];Eta=lxi*[0.0];O=0;

		while O<8: # iterate until convergence of water level
			hi=[ashi[i]+Eta[i] for i in range(lxi)] # water depth        
			Sxx=[0.5*Ewi[i]*(4.0*ki[i]*hi[i]/sinh(2.0*ki[i]*hi[i])+1.0) for i in range(lxi)] # wave radiation stress
			Rxx=[2.0*Eri[i] for i in range(lxi)] # roller radiation stress
			# estimate MWL along Xshore transect
			temp1=[Sxx[i]+Rxx[i] for i in range(lxi)]
			temp2=gradient(temp1,dx)

			Integr=[(-temp2[i]+fxi[i])/(rho*g*hi[i]) for i in range(lxi)]
			Eta[0]=Etao
			Eta[1]=Eta[0]+Integr[0]*dx
			for i in range(1,lxi-2):
				Eta[i+1]=Eta[i-1]+Integr[i]*2*dx
			Eta[lxi-1]=Eta[lxi-2]+Integr[lxi-1]*dx
			O=O+1
		F=interp1d(Xi,Eta);Eta=F(X);Etas=F(X);

		if sum(Vegxloc)<>0: # compute if there's vegetation (time saver)
			Sxxs=lx*[0.0]; Rxxs=lx*[0.0];Etas=lx*[0.0];dx=1;O=0;
			while O<5: # iterate until convergence of water level
				h=[ash[i]+Etas[i] for i in range(lx)] # water depth        
				Sxxs=[0.5*Ews[i]*(4.0*k[i]*h[i]/sinh(2.0*k[i]*h[i])+1.0) for i in range(lx)] # wave radiation stress
				Rxxs=[2.0*Ers[i] for i in range(lx)] # roller radiation stress
				# estimate MWL along Xshore transect
				temp1=[Sxxs[i]+Rxxs[i] for i in range(lx)]
				temp2=gradient(temp1,dx)

				Integr=[(-temp2[i])/(rho*g*h[i]) for i in range(lx)]
				Etas[0]=Etao
				Etas[1]=Eta[0]+Integr[0]*dx
				for i in range(1,lx-2):
					Etas[i+1]=Etas[i-1]+Integr[i]*2*dx
				Etas[lx-1]=Etas[lx-2]+Integr[lx-1]*dx
				O=O+1    

		Ubot=[pi*H[ii]/(To*sinh(k[ii]*h[ii])) for ii in range(lx)] # bottom velocity
		H=[H[ii]*sqrt(2) for ii in range(lx)]
		Hs=[Hs[ii]*sqrt(2) for ii in range(lx)]

		# Ds1=num.array(H); Ds2=num.array(Hs)
		Ds1=gradient(Ef,dx);Ds1=[-Ds1[ii] for ii in range(len(Ef))]
		if sum(Vegxloc)<>0: # compute if there's vegetation (time saver)
			Ds2=gradient(Efs,dx);Ds2=[-Ds2[ii] for ii in range(len(Efs))]
		else:
			Ds2=Ds1;

		# Diss=mean(Ds1[find(num.array(Dveg)<>0)]**3/(Ds2[find(num.array(Dveg)<>0)]**3)) # dissipation difference
		Diss=[Ds1,Ds2] # dissipation difference
		if num.isnan(mean(Diss)) == True:
			Diss=[Ds1*0.0,Ds1*0.0]

		return H,Eta,Hs,Etas,Diss,Ubot # returns: wave height, wave setup, wave height w/o veg., wave setup w/o veg, wave dissipation, bottom wave orbital velocity over the cross-shore domain

	# a wave model independent of habitat with dissipation only due to breaking and bottom friction 
	def WaveModelSimple(X,h,Ho,To,Cf):
		# constants
		g=9.81;rho=1024.0;B=1.0;Beta=0.05
		lx=len(X);dx=X[1]-X[0];

		if diff(h).all()==0:
			flat=1;
		else:
			flat=0

		# initialize vectors for wave model
		H=lx*[0.0];Eta=lx*[0.0];L=lx*[0.0]
		Db=lx*[0.0];Df=lx*[0.0]
		Er=lx*[0.0];Ef=lx*[0.0];Br=lx*[0.0]
		C=lx*[0.0];n=lx*[0.0];Cg=lx*[0.0]
		k=lx*[0.0];

		# wave parameter at 1st grid pt
		ash=[h[ii] for ii in range(lx)] # ash is same as h, but is now an independent variable
		fp=1.0/To; sig=2.0*pi*fp
		H[0]=Ho/sqrt(2.0) # transform significant wave height to rms wave height
		Ef[0]=0.125*rho*g*H[0]**2*Cg[0] # energy flux @ 1st grid pt
		k[0]=iterativek(sig,h[0]) # wave number at 1st grid pt
		L[0]=2.0*pi/k[0] # wave length @ 1st grid pt
		n[0]=0.5*(1+(2.0*k[0]*h[0]/sinh(2.0*k[0]*h[0]))) # to compute Cg at 1st grid pt
		C[0]=L[0]/To;Cg[0]=C[0]*n[0] # phase and group velocity at 1st grid pt

		So=Ho/L[0] # deep water wave steepness
		Gam=0.5+0.4*num.tanh(33.0*So) # Gam from Battjes and Stive 85 (as per Alsina & Baldock)

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

			if flat:
				Er[xx+1]=0.0
				Db[xx+1]=0.0            

			Df[xx+1]=1.0*rho*Cf/(12.0*pi)*(2.0*pi*fp*H[xx+1]/sinh(k[xx+1]*h[xx+1]))**3.0 # dissipation due to bottom friction 

		Ew=lx*[0.0];Ew=[0.125*rho*g*(H[i]**2.0) for i in range(lx)]

		# estimate MWS
		val=1.0;dell=1.0;O=0.0 # for setup calculations
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
			if flat:
				Eta[0]=0
			else:
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
		return H,Eta # returns wave height and setup over the cross-shore domain

	# wave attenuation by coral reefs
	def WavesCoral(Ho,AlphF ,AlphR,he,hr,Wr,Cf,dx):
		BrkRim=0;BrkFace=0;Kp=0.0

		if AlphF+AlphR+he==0: # user doesn't have data
			Kp=mean(kp) # assume waves break on face
			ha=hr
		else:
			# reef profile
			Xreef=num.arange(0.0,10000.0,dx)
			Yreef=AlphR*Xreef-he
			loc=find(Yreef>-hr)
			Yreef[loc]=-hr # reef profile

			D=To*sqrt(g/hr) # relative subm
			Lo=g*To**2.0/(2.0*pi)
			if D==8: # first check for breaking location
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
		H_r=num.arange(0,lx)*0.0

		if Kp<> 0: # wave break on reef face or rim
			# wave Setup
			while delta>0.01:
				EtaR=3.0/(64.0*pi)*Kp*Ho**2*To*g**.5/((Etar+ha)**1.5)
				delta=abs(EtaR-Etar);Etar=EtaR
			# Etar=Etar[0]    
			H_r[0]=0.42*(hr+Etar) # wave height at the offshore edge of reef
		else:
			Etar=0
			H_r[0]=Ho

		depth_flat=num.arange(0,lx)*0.0+Etar+hr;X=num.arange(0,lx)
		H_r,eta=WaveModelSimple(X,depth_flat,H_r[0],To,Cf)
		H_r=num.array(H_r)

		Etar=H_r*0.0+Etar # arrays of wave height and setup on the reef                   
		return H_r,Etar

	# wave attenuation by reef breakwater
	def BreakwaterKt(Hi,To,hi,hc,Cwidth,Bwidth,OysterReefType):
		Lo=9.81*To**2.0/(2.0*pi)
		Rc=hc-hi # depth of submergence
		Boff=(Bwidth-Cwidth)/2.0 # base dif on each side
		ksi=(hc/Boff)/sqrt(Hi/Lo)

		if hc/hi<0.5: # the reef is too submerged
			Kt=1
			gp.AddError("\nYour reef is too small and we cannot compute the wave transmission coefficient. Please increase your reef height so it's closer to the water surface, or assume that it does not protect your shoreline.")
			raise Exception
		if hc/hi>1.25: # the reef is too submerged
			Kt=1
			gp.AddError("\nYour reef is too high above the water level and we cannot compute the wave transmission coefficient. Please assume that no waves pass through and the transmission coefficient is zero.")
			raise Exception
		if abs(Rc/Hi)>5:
			Kt=1
			gp.AddError("\nYour reef is too small and outside the range of validity of our model. Please increase your reef height so it's closer to the water surface.")
			raise Exception
		if OysterReefType=="Trapezoidal": # it's not a reef ball
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

		elif OysterReefType=="Dome": # it's a reef ball
			Kt=1.616-31.322*Hi/(9.81*To**2)-1.099*hc/hi+0.265*hc/Bwidth # D'Armono and Hall
			if hc>hi:
				Kt=1
				gp.AddError("\nThe reef balls are emerged. We cannot compute the wave transmission coefficient.")
				raise Exception       
			elif hc<0.1*hi:
				Kt=1
				gp.AddError("\nThe reef is smaller than 10% of your water depth.  We cannot compute a transmission coeffcient and you can assume that they do not protect your shoreline.")
				raise Exception       
			elif Kt<0:
				Kt=1  # no transmission; breakwater fully emerged
				gp.AddError("\nWe were not able to compute a transmission coefficient. Please create a reef that is closer to the water level, but not too emerged.")
				raise Exception       
			elif Kt>1:
				Kt=1
		return Kt

	# K&D Erosion model
	def ErosionKD(A,Ho,TWL,B,D,W,m,count):
		# A: Shape factor for the Equilibrium Beach Profile associated with the sediment size
		# Ho: Deep Water Wave Height
		# TWL: Total Water Level (Surge plus Setup)
		# B: Berm Elevation
		# D: Dune Height
		# W: Dune Width
		# m: Foreshore Slope
		global BetaKD
		# constants
		BD=D+B; 
		msg=0 # tracks whether or not certain messages should be displayed
		B_adj1=0 # The berm may need to be elevated if TWL value is too great (first method)
		B_adj2=0 # The berm may need to be elevated if TWL value is too great (second method)
		inundation=0 # tracks whether or not the backshore is inundated 
		if TWL>BD:
			# TWL=BD
			if count == 0:
				count += 1
				gp.AddMessage("...Water level is higher than your backshore elevation. You will probably experience flooding.") # only display this message if it is the first time the backshore is inundated
			inundation=1

		# Erosion model 1
		hb=(((Ho**2.0)*g*To/(2*pi))/2.0)**(2.0/5.0)/(g**(1.0/5.0)*0.73**(4.0/5.0)) # breaking depth
		xb=(hb/A)**1.5;Hb=0.78*hb # cross-shore breaking location
		if B + hb <= TWL/2: # if this condition is true, the K&D model with blow up or given erroneous values 
			gp.AddWarning("...The surge and wave setup are too great to run erosion on your beach!  The berm elevation has been increased to generate erosion profiles!")
			testval = -9999
			while testval < .1:
				B_adj1 += 0.5
				testval = B_adj1 + B + hb - TWL/2
			BD1 = B+B_adj1+D
			B1 = B+B_adj1
			gp.AddWarning("...The berm has been increased by"+str(B_adj1)+"meters in elevation.") 
		else:
			BD1=BD
			B1=B

		Term1=xb-hb/m
		Rinf1=(TWL*Term1)/(B1+hb-TWL/2.) # erosion without taking width into account   
		if D>0:
			Rinf=(TWL*Term1)/(BD1+hb-TWL/2.)-W*(B1+hb-0.5*TWL)/(BD1+hb-0.5*TWL) # potential erosion distance
			if Rinf < 0:
				Rinf=(TWL*Term1)/(B1+hb-TWL/2.); msg+=1
		else:
			Rinf=Rinf1

		TS=(320.0*(Hb**(3.0/2.0)/(g**.5*A**3.0))*(1.0/(1.0+hb/BD1+(m*xb)/hb)))/3600.0 # response time scale
		BetaKD=2.0*pi*(TS/StormDur)
		expr="num.exp(-2.0*x/BetaKD)-num.cos(2.0*x)+(1.0/BetaKD)*num.sin(2.0*x)" # solve this numerically
		fn=eval("lambda x: "+expr)
		z=FindRootKD(fn,pi,pi/2) # find zero in function,initial guess from K&D
		R01=0.5*Rinf*(1.0-num.cos(2.0*z)) # final erosion distance


		# 2nd method
		x=num.arange(0,10000,1)
		y=A*x**(2.0/3);y=y[::-1]
		y=y[find(y>0.5)];x=x[find(y>0.5)]    
		Hs,Etas=WaveModelSimple(x,y,Ho,To,0.01) 
		hb=y[argmin(Etas)]
		xb=(hb/A)**1.5 # surf zone width
		if B + hb <= TWL/2: # if this condition is true, the K&D model with blow up or given erroneous values 
			gp.AddWarning("...The surge and wave setup are too great to run erosion on your beach!  The berm elevation as been increased to generate erosion profiles!")
			testval = -9999
			while testval < 1:
				B_adj2 += 0.5
				testval = B_adj2 + B + hb - TWL/2
			BD2 = B+B_adj2+D
			B2 = B+B_adj2
		else:
			BD2 = BD
			B2 = B
		Term1=xb-hb/m
		Rinf1=(TWL*Term1)/(B2+hb-TWL/2) # erosion without taking width into account
		if D>0:
			Rinf=(TWL*Term1)/(BD2+hb-TWL/2)-W*(B2+hb-0.5*TWL)/(BD2+hb-0.5*TWL) # potential erosion distance
			if Rinf < 0:
				Rinf=(TWL*Term1)/(B2+hb-TWL/2); msg+=1
		else:
			Rinf = Rinf1
		TS=(320.0*(Hb**(3.0/2.0)/(g**.5*A**3.0))*(1.0/(1.0+hb/BD2+(m*xb)/hb)))/3600.0 # response time scale
		BetaKD=2.0*pi*(TS/StormDur)
		expr="num.exp(-2.0*x/BetaKD)-num.cos(2.0*x)+(1.0/BetaKD)*num.sin(2.0*x)" # solve this numerically
		fn=eval("lambda x: "+expr)
		z=FindRootKD(fn,pi,pi/2) # find zero in function,initial guess from K&D
		R02=0.5*Rinf*(1.0-num.cos(2.0*z)) # final erosion distance        
		R0=max([R01,R02])
		if R0<0:    R0=0
		if msg==2:
			gp.AddMessage("...Berm is so wide that the dune is not eroding.")
			msg=2;
		B_adj = min([B_adj1,B_adj2]) # return the minimum of the required berm elevation adjustments ('0' if one or both methods did not require adjusting 
		# returns: the average retreat value, the retreat values by method 1, by method 2, berm increase value, whether or not inundation occured, and a variable tracking what message to display later. 
		return R0,R01,R02,B_adj,inundation,count

	def MudErosion(Uc,Uw,h,To,me,Cm):
		rho=1024.0;nu=1.36e-6;d50=0.00003
		ks=2.5*d50;kap=0.4

		# current
		if max(Uc)<>0:
			us1=0.01;zo1=0.01;dif=10 # initial value for u* and zo  
			Tc=h*0; # shear stress
			while dif>1e-4:
				zo2=ks/30*(1-num.exp(-us1*ks/(27*nu)))+nu/(9*us1)
				us2=kap*Uc/(num.log(h/zo2)-1)
				dif1=abs(us1-us2);dif2=abs(zo1-zo2);dif=mean(dif1+dif2);
				zo1=zo2;us1=us2;
			Tc=rho*us1**2 # shear stress due to current
		else:    Tc=h*0

		# waves
		Rw=Uw**2*To/(2*num.pi)/nu
		fw=0.0521*Rw**(-0.187) # smooth turbulent flow
		Tw=0.5*rho*fw*Uw**2

		# combined wave and current
		temp=Tc*(1+1.2*(Tw/(Tc+Tw))**3.2)
		Trms=(temp**2+0.5*Tw**2)**0.5

		# erosion 
		Te=h*0+0.0012*Cm**1.2 # erosion threshold
		dmdt=me*(Trms-Te) # erosion rate
		dmdt[find(dmdt<=0)]=0
		Erosion=3600*dmdt/Cm*100 # rate of bed erosion [cm/hr]

		# erosion: bed erosion rate, Trms: wave and current shear stress, Tc: current shear stress, Tw: wave shear stress, Te: erosion threshold
		return Erosion,Trms,Tc,Tw,Te

	def ParseSegment(X):
		# this function returns location of beg. and end of segments that have same value in a vector
		# input: X vector 
		# output: 
		# - Beg and End_ar: indices of beginning and end of segments that have same value
		# - L_seg, Values: length of segment and value of terms in that segment

		X=num.array(X)    
		End_ar=nonzero(diff(X));End_ar=End_ar[0] # ID locations BEFORE number changes

		Beg_ar=num.append([0],End_ar+1) # index of bef. of segments
		End_ar=num.append(End_ar,len(X)-1) # index of end of segments
		L_seg=End_ar-Beg_ar+1 # length of the segments
		Values=X[Beg_ar]

		return Values,Beg_ar,End_ar,L_seg


	################################################
	#### READING DEPTH PROFILE AND EXCEL INPUTS ####
	################################################

	try:
		gp.AddMessage("\nReading Erosion Protection Excel file and depth profile...")
		# constants
		g=9.81;rho=1024.0
		Fig2=0; # for HTML

		# input from user via Excel table
		xlApp=Dispatch("Excel.Application")
		xlApp.Visible=0
		xlApp.DisplayAlerts=0
		xlApp.Workbooks.Open(InputTable)
		cell1=xlApp.Worksheets("ModelInput")
		cell3=xlApp.Worksheets("ReefShapeFactor")

		# read general information
		MSL=cell1.Range("f5").Value # mean sea level
		MHW=cell1.Range("g5").Value # mean high water

		XelVal=cell1.Range("i10").Value
		if XelVal==1:
			sand=1
			mud=0
		elif XelVal==2:
			sand=0
			mud=1
		try:
			# read in user's cross-shore profile
			TextData=open(CSProfile,"r") 
			X_lst=[];h_lst=[]
			for line in TextData.read().strip("\n").split("\n"):
				linelist=[float(s) for s in line.split("\t")] # split the list by comma delimiter
				X_lst.append(linelist[0])
				h_lst.append(linelist[1])
			TextData.close()
			Xinit = num.array(X_lst)
			X=num.array(X_lst);h=num.array(h_lst)
	
			# Adjust Bathy
			if any(h[0:100]>0) or h[0]>h[-1]:
				flip=1;
			else:
				h=h[::-1] # reverse order if profile starts offshore
				flip=0 # flip later to profile starts offshore
	
			h=h-MSL; # adjust water level so that 0 is at MSL
	
			if mud: # if there's a marsh or a mangrove
				gp.AddMessage("...your backshore is a *marsh/mangrove*")            
				h=h-S # add surge level by decreasing depth (increasing the absolute value)
				out=find(h>-0.1)
				h=h[out[-1]+1:-1];X=X[out[-1]+1:-1] # only keep values that are below water
				dx=abs(X[1]-X[0]);m=abs(h[-1]-h[-int(10.0/dx)])/10 # average slope 10m from end of transect
	
			elif sand: # it's a beach
				gp.AddMessage("...your backshore is a *sandy beach*")            
				keep=find(h<-0.1)
				h=h[keep];X=X[keep] # only keep values that are below water
			h=-h # depth is now positive
			if flip: # original profile starts onshore
				h=h[::-1] # reverse order if profile starts at onshore
	
		except:
			gp.AddError(msgReadCSProfile)
			raise Exception		

		# read muddy shoreline information
		Cm=cell1.Range("e25").Value # dry density
		me=cell1.Range("f25").Value # erosion constant
		# read beach Information
		d50=cell1.Range("e15").Value # sediment size
		B1=cell1.Range("j17").Value # berm elevation (relative to MSL)
		W1=cell1.Range("j16").Value # berm width
		D1=cell1.Range("j15").Value # dune height
		Dred=cell1.Range("e45").Value # percent reduction of dune height due to management action
		D2=(100-Dred)*D1/100 # new dune height        
		Slope=cell1.Range("j18").Value
		if Slope<>0:
			m=1.0/Slope # foreshore slope=1/Slope
		else:
			m=0
		A=cell1.Range("e17").Value # sediment scale factor

		if sand+mud<>1:
			gp.AddWarning("...You didn't specify a backshore type.  We won't be able to estimate amount of erosion.")

		# read vegetation information
		XelVal=cell1.Range("f49:m53").Value # all vegetation information
		temp1=XelVal[0] # mangrove roots
		hogr=temp1[0] # height of mangrove roots
		dogr=temp1[1] # diameter of mangrove roots
		Nogr=temp1[2] # density of mangrove roots 
		Xsg=temp1[3] # shoreward edge of mangrove
		Xog=temp1[4] # offshore edge of mangrove
		XsMAg=temp1[5] # shoreward edge after management action
		XoMAg=temp1[6] # offshore edge after management action          
		densdeltagr = temp1[7]
		if XsMAg == 0 and XoMAg == 0: # if there is no habitat footprint, then density must be '0' as well
			del densdeltagr
			densdeltagr = 100.0
		dNogr=(1-densdeltagr *1.0/100.0) # density change in mangrove roots 


		temp1=XelVal[2] # mangrove canopy 
		hogc=temp1[0] # height of mangrove canopy
		dogc=temp1[1] # diameter of mangrove canopy
		Nogc=temp1[2] # density of mangrove canopy       
		densdeltagc = temp1[7]
		if XsMAg == 0.0 and XoMAg == 0.0: # if there is no habitat footprint, then density must be '0' as well
			del densdeltagc
			densdeltagc = 100.0      
		dNogc=(1-densdeltagc *1.0/100.0)# density change in mangrove canopy

		temp1=XelVal[1] # mangrove trunks
		densdeltagt = temp1[7]
		if XsMAg == 0.0 and XoMAg == 0.0: # if there is no habitat footprint, then density must be '0' as well
			del densdeltagt
			densdeltagt = 100.0     
		hogt=temp1[0] # height of mangrove trunks
		dogt=temp1[1] # diameter of mangrove trunks
		Nogt=temp1[2] # density of mangrove trunks
		dNogt=(1-densdeltagt *1.0/100.0) # density change in mangrove trunks

		if Xog>0:
			gp.AddWarning("...Mangrove landward edge should be above mean sea level.  We'll move it for you.")            
			Xog=-Xog
		if Xsg<Xog:
			gp.AddWarning("...You switched mangrove edge distances.  We'll change them for you.")
			temp=Xog;Xog=Xsg;Xsg=temp
		if XsMAg<XoMAg:
			gp.AddWarning("...You switched mangrove edge distances for the management action.  We'll change them for you.")
			temp=XoMAg;XoMAg=XsMAg;XsMAg=temp
		if Xog<min(Xinit) and Xsg<min(Xinit):
			gp.AddWarning("...The mangrove footprint you applied lies completely outside the extent of your topo/bathy profile. The mangrove has been excluded in the analysis. Check your inputs.")
			Xog = 0; Xsg = 0; XoMAg =0; XsMAg = 0;
		elif Xog<min(Xinit):
			gp.AddWarning("...The mangrove footprint you applied extends beyond the limits of your topo/bathy profile. The inland limit of the mangrove has been set to the inland limit of the topo/bathy profile. Check your inputs.")
			Xog = min(Xinit)
		if XoMAg<min(Xinit):
			gp.AddWarning("...The mangrove footprint you applied for the managment action extends beyond the limits of your topo/bathy profile. The inland limit of the mangrove post management action has been set to the inland limit of the topo/bathy profile. Check your inputs.")
			XoMAg = min(Xinit)			
		
		if XsMAg+XoMAg<>0: # if there is vegetation in management action
			if XsMAg>Xsg:
				gp.AddWarning("...Your impacted mangrove footprint has to be within your initial conditon footprint.  We'll change it for you.")
				XsMAg=Xsg
			if XoMAg<Xog:
				gp.AddWarning("...Your impacted mangrove footprint has to be within your initial conditon footprint.  We'll change it for you.")
				XoMAg=Xog

		temp1=XelVal[3]
		hor=temp1[0] # height of marsh
		dor=temp1[1] # diameter of marsh
		Nor=temp1[2] # density of marsh
		Xsr=temp1[3] # shoreward edge of marsh
		Xor=temp1[4] # offshore edge of marsh
		XsMAr=temp1[5] # shoreward edge of marsh
		XoMAr=temp1[6] # offshore edge of marsh
		densdeltar = temp1[7]
		if XsMAr == 0 and XoMAr == 0: # if there is no habitat footprint, then density must be '0' as well
			del densdeltar
			densdeltar = 100.0        
		dNor=(1-densdeltar *1.0/100.0) # density change in marsh
		if Xor>0:
			gp.AddWarning("...Marshes landward edge should be above mean sea level.  We'll move it for you.")  
			Xor=-Xor;
		if Xsr<Xor:
			gp.AddWarning("...You switched marsh offshore and shoreward edge.  We'll change them for you.")
			temp=Xor;Xor=Xsr;Xsr=temp
		if XsMAr<XoMAr:
			gp.AddWarning("...You switched marsh offshore and shoreward edge for the management action.  We'll change them for you.")
			temp=XoMAr;XoMAr=XsMAr;XsMAr=temp
		if Xor<min(Xinit) and Xsr<min(Xinit):
			gp.AddWarning("...The marsh footprint you applied lies completely outside the extent of your topo/bathy profile. The marsh has been excluded in the analysis. Check your inputs.")
			Xor = 0; Xsr = 0; XoMAr =0; XsMAr = 0;
		elif Xor<min(Xinit):
			gp.AddWarning("...The marsh footprint you applied extends beyond the limits of your topo/bathy profile. The inland limit of the marsh has been set to the inland limit of the topo/bathy profile. Check your inputs.")
			Xor = min(Xinit)
		if XoMAr<min(Xinit):
			gp.AddWarning("...The marsh footprint you applied post management extends beyond the limits of your topo/bathy profile. The inland limit of the marsh post management has been set to the inland limit of the topo/bathy profile. Check your inputs.")
			XoMAr = min(Xinit)			
		if XsMAr+XoMAr<>0: # if there is vegetation in management action
			if XsMAr>Xsr:
				gp.AddWarning("...Your impacted marsh footprint has to be within your initial conditon footprint.  We'll change it for you.")
				XsMAr=Xsr
			if XoMAr<Xor:
				gp.AddWarning("...Your impacted marsh footprint has to be within your initial conditon footprint.  We'll change it for you.")
				XoMAr=Xor

		temp1=XelVal[4]
		hos=temp1[0] # height of seagrass
		dos=temp1[1] # diameter of seagrass
		Nos=temp1[2] # density of seagrass
		Xss=temp1[3] # shoreward edge of seagrass
		Xos=temp1[4] # offshore edge of seagrass
		XsMAs=temp1[5] # shoreward edge after management action
		XoMAs=temp1[6] # offshore edge after management action
		densdeltas = temp1[7]
		if XsMAs == 0 and XoMAs == 0: # if there is no habitat footprint, then density must be '0' as well
			del densdeltas
			densdeltas = 100.0        
		dNos=(1-densdeltas *1.0/100.0) # density change in seagrass        
		if Xss>Xos:
			gp.AddWarning("...You switched seagrass offshore and shoreward edge.  We'll change them for you.")
			temp=Xos;Xos=Xss;Xss=temp

		if XsMAs>XoMAs:
			gp.AddWarning("...You switched seagrass offshore and shoreward edge for the management action.  We'll change them for you.")
			temp=XoMAs;XoMAs=XsMAs;XsMAs=temp
		if Xos>max(Xinit) and Xss>max(Xinit):
			gp.AddWarning("...The seagrass footprint you applied lies completely outside the extent of your topo/bathy profile. The seagrass has been excluded in the analysis. Check your inputs.")
			Xos = 0; Xss = 0; XoMAs =0; XsMAs = 0;
		elif Xos>max(Xinit):
			gp.AddWarning("...The seagrass footprint you applied extends beyond the limits of your topo/bathy profile. The offshore limit of the seagrass has been set to the offshore limit of the topo/bathy profile. Check your inputs.")
			Xos = max(Xinit)
		if XoMAs>max(Xinit):
			gp.AddWarning("...The seagrass footprint you applied post management action extends beyond the limits of your topo/bathy profile. The offshore limit of the seagrass post management action has been set to the offshore limit of the topo/bathy profile. Check your inputs.")
			XoMAs = max(Xinit)			
		if XsMAs+XoMAs<>0: # if there is vegetation in management action
			if XsMAs<Xss:
				gp.AddWarning("...Your impacted seagrass footprint has to be within your initial conditon footprint.  We'll change it for you.")
				XsMAs=Xss
			if XoMAs>Xos:
				gp.AddWarning("...Your impacted seagrass footprint has to be within your initial conditon footprint.  We'll change it for you.")
				XoMAs=Xos		

		# read coral info
		CoralType=cell1.Range("d56").Value
		if CoralType<>"Barrier" and CoralType <> "Fringe" and CoralType <> "Fringe Lagoon":
			CoralType = "None"        
		XelVal=cell1.Range("e56:l56").Value;XelVal=XelVal[0] # all coral info
		Xco=XelVal[1] # offshore distance
		Xcn=XelVal[0] # nearshore distance

		AlphF=XelVal[2] # reef face slope
		AlphR=XelVal[3] # reef rim slope
		he=XelVal[4] # reef rim edge
		hr=XelVal[5] # reef top depth
		Wr=XelVal[6] # reef top width
		CoralMA=XelVal[7]
		TanAlph=cell3.Range("a2:a202").Value;TanAlph=num.array(TanAlph)
		kp=cell3.Range("b2:b202").Value;kp=num.array(kp) # reef shape factor
		kp2=cell3.Range("c2:c202").Value;kp2=num.array(kp2) # reef shape factor
		if CoralMA is None:
			CoralMA="None"
		if AlphF+AlphR+he+hr+Wr==0:
			CoralMA="None"        
		if Xcn>Xco:
			gp.AddWarning("...You switched coral offshore and shoreward edge.  We'll change them for you.")
			temp=Xco;Xco=Xcn;Xcn=temp
		if Xco>max(Xinit):
			gp.AddWarning("...The coral reef footprint you applied is offshore of the limits of your profile. The reef will be placed at the offshore limit of the profile (Barrier Reef).")
			Xcn=-1.0; Xco = -1.0; CoralType = "Barrier"

		# read oyster reef information
		OysterReefType=cell1.Range("d59").value;
		if OysterReefType<>"Trapezoidal" and OysterReefType<>"Dome":
			OysterReefType = "None"
		XelVal=cell1.Range("e59:i59").value;XelVal=XelVal[0] # all oyster information
		Xr=XelVal[0] # distance from shoreline; need to reverse because user enter with shoreline at X=0
		hc=XelVal[1] # reef height
		Bw=XelVal[2] # base width
		Cw=XelVal[3] # crest width
		OysterMA=XelVal[4]
		if OysterMA<>"Rmv":
			OysterMA="None"
		if Xr+hc+Bw+Cw==0:
			OysterMA="None"; OysterReefType="None"

		xlApp.ActiveWorkbook.Close(SaveChanges=0)
		xlApp.Quit()

		if Xsr+Xor+XsMAr+XoMAr<>0 and Xsg+Xog+XsMAg+XoMAg<>0:
			gp.AddError(msgTwoHabitats)
			raise Exception

	except:
		xlApp.ActiveWorkbook.Close(SaveChanges=0)
		xlApp.Quit()
		gp.AddError(msgReadExcel)
		raise Exception


	#########################################
	######## CREATING INPUTS FOR RUN ########
	#########################################

	# open HTML file
	AddMActionHeader = 'yes'
	htmlfile=open(outputws+"OutputWaveModel_"+subwsStr+".html","w")
	htmlfile.write("<html>\n")
	htmlfile.write("<title>Marine InVEST - Wave, Erosion, and Inundation</title>")
	htmlfile.write("<CENTER><H1>Coastal Protection - Tier 1</H1><H2>Nearshore Waves and Erosion Results ("+subwsStr+")<br></H2></CENTER>")

	htmlfile.write("<HR><H2><u>Site Information</u></H2>")
	htmlfile.write("<table border=\"0\" width=\"1100\" cellpadding=\"5\" cellspacing=\"10\">")
	if sand==1:
		htmlfile.write("<li>You have a sandy beach, and the average sediment size is: "+str(d50)+"mm<br><p>")
	elif mud==1:
		htmlfile.write("<li>Your bed is composed of fines/consolidated sediments. <br>")

	htmlfile.write("<li>The tidal range at your site is: "+str(round(MHW,1))+"m (high tide value).<br>")
	if MHW < 2:
		htmlfile.write("<li> It is <i>microtidal</i> (Tidal Range < 2m)<br>")
	elif MHW <= 4:
		htmlfile.write("<li> It is <i>meso-tidal</i> (2 <= Tidal Range <= 4m)<br>")
	else:
		htmlfile.write("<li> It is <i>macro-tidal</i> (Tidal Range > 4m)<br>")
	htmlfile.write("</td></tr></table>")

	try:
		# compute offshore wave height
		if WaveErosionQuestion=="(2) No, please compute these values from wind speed and fetch distance": 
			Us=Us;Ft=Ft;depth=depth
			Ho,To=WindWave(Us,Ft,depth) # compute wave from wind speed

		# fill up vegetation information
		VegXloc=num.array(X)*0.0;VegXlocMA=num.array(X)*0.0; # x-axis locating habitats
		HabTrigg=num.array(X)*0.0+100.0; # habitats function triggers; assume it's all wave model for now
		NRoots=num.array(X)*0.0;NRootsMA=num.array(X)*0.0; # density of roots
		dRoots=num.array(X)*0.0;dRootsMA=num.array(X)*0.0; # diameter of roots
		hRoots=num.array(X)*0.0;hRootsMA=num.array(X)*0.0; # height of roots
		NTrunk=num.array(X)*0.0;NTrunkMA=num.array(X)*0.0; # density of trunks
		dTrunk=num.array(X)*0.0;dTrunkMA=num.array(X)*0.0; # diameter of trunks
		hTrunk=num.array(X)*0.0;hTrunkMA=num.array(X)*0.0; # height of trunks
		NCanop=num.array(X)*0.0;NCanopMA=num.array(X)*0.0; # density of canopy
		dCanop=num.array(X)*0.0;dCanopMA=num.array(X)*0.0; # diamter of canopy
		hCanop=num.array(X)*0.0;hCanopMA=num.array(X)*0.0; # height of canopy

		# prepare seagrass data
		if Xos+Xss<>0: # there's a seagrass bed
			temp1=X[-1]-X[0]-Indexed(X,Xos) # offshore boundary of habitat
			temp2=X[-1]-X[0]-Indexed(X,Xss) # shoreward boundary of habitat
			VegXloc[temp1:temp2+1]=3 # '3' for seagrass
			NTrunk[temp1:temp2+1]=Nos 
			dTrunk[temp1:temp2+1]=dos
			hTrunk[temp1:temp2+1]=hos

		if XoMAs+XsMAs<>0: # there's a seagrass bed after management action
			temp1=X[-1]-X[0]-Indexed(X,XoMAs) # offshore boundary of habitat
			temp2=X[-1]-X[0]-Indexed(X,XsMAs) # shoreward boundary of habitat
			VegXlocMA[temp1:temp2+1]=3 # '3' for seagrass
			NTrunkMA[temp1:temp2+1]=Nos *dNos
			dTrunkMA[temp1:temp2+1]=dos
			hTrunkMA[temp1:temp2+1]=hos

		if Xos+Xss+XoMAs+XsMAs<>0: # there's a seagrass bed 
			gp.AddMessage("...a *seagrass bed* is present")
			if AddMActionHeader == 'yes':
				htmlfile.write("<HR><H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
			htmlfile.write("<b>A *seagrass bed* is present at your site</b>")
			if XoMAs+XsMAs==0:
				htmlfile.write("...you assumed that the whole seagrass bed will be removed.<br>")
			elif XoMAs+XsMAs<>Xos+Xss:
				if densdeltas == 0.0:
					htmlfile.write("...you assumed that the seagrass bed footprint will be reduced by "+str(abs(abs(Xos-Xss)-abs(XoMAs-XsMAs)))+"m.<br>")
				else:
					htmlfile.write("...you assumed that the seagrass bed footprint will be reduced by "+str(abs(abs(Xos-Xss)-abs(XoMAs-XsMAs)))+"m and that the density of the seagrass will reduce by " +str(densdeltas)+"%.<br>")
			elif XoMAs+XsMAs==Xos+Xss:
				if densdeltas == 0.0:
					htmlfile.write("...you assumed that the seagrass bed will not be affected by a particular management action.<br>")
				else:
					htmlfile.write("...you assumed that the seagrass bed footprint will be unaffected but the density will reduce by " +str(densdeltas)+"%"+".<br>")
		else:
			gp.AddMessage("...you do NOT have a seagrass bed")

		# prepare marsh data
		if Xor+Xsr<>0: # there's a marsh 
			temp2=X[-1]-X[0]-Indexed(X,Xor) # landward boundary of habitat
			temp1=X[-1]-X[0]-Indexed(X,Xsr) # shoreward boundary of habitat
			VegXloc[temp1:temp2+1]=2 # '2' for marsh
			NTrunk[temp1:temp2+1]=Nor 
			dTrunk[temp1:temp2+1]=dor
			hTrunk[temp1:temp2+1]=hor

		if XoMAr+XsMAr<>0: # there's a marsh after management action
			temp2=X[-1]-X[0]-Indexed(X,XoMAr) # landward boundary of habitat
			temp1=X[-1]-X[0]-Indexed(X,XsMAr) # shoreward boundary of habitat
			VegXlocMA[temp1:temp2+1]=2 # '2' for marsh
			NTrunkMA[temp1:temp2+1]=Nor *dNor
			dTrunkMA[temp1:temp2+1]=dor
			hTrunkMA[temp1:temp2+1]=hor

		if Xor+Xsr+XoMAr+XsMAr<>0: # there's a marsh 
			gp.AddMessage("...a *marsh* is present")
			if AddMActionHeader == 'yes':
				htmlfile.write("<HR><H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
			htmlfile.write("<b>A *marsh* is present at your site</b>")
			if XoMAr+XsMAr==0:
				htmlfile.write("...you assumed that the whole marsh will be removed.<br>")
			elif XoMAr+XsMAr<>Xor+Xsr:
				if densdeltar == 0.0:
					htmlfile.write("...you assumed that the marsh footprint will be reduced by "+str(abs(abs(Xor-Xsr)-abs(XoMAr-XsMAr)))+"m. <br>")
				else:
					htmlfile.write("...you assumed that the marsh footprint will be reduced by "+str(abs(abs(Xor-Xsr)-abs(XoMAr-XsMAr)))+"m and that the density of the marsh will reduce by " +str(densdeltar)+"%. <br>")
			elif XoMAr+XsMAr==Xor+Xsr:
				if densdeltar == 0.0:
					htmlfile.write("...you assumed that the marsh will not be affected by a particular management action.<br>")
				else:
					htmlfile.write("...you assumed that the marsh footprint will be unaffected but the density will reduce by " +str(densdeltar)+"%"+".<br>")
		else:
			gp.AddMessage("...you do NOT have a marsh")

		# prepare mangrove data
		if Xog+Xsg<>0: # there's a mangrove
			temp2=X[-1]-X[0]-Indexed(X,Xog) # landward boundary of habitat
			temp1=X[-1]-X[0]-Indexed(X,Xsg) # shoreward boundary of habitat
			VegXloc[temp1:temp2+1]=1 # '1' for mangroves
			NRoots[temp1:temp2+1]=Nogr
			dRoots[temp1:temp2+1]=dogr
			hRoots[temp1:temp2+1]=hogr 
			NTrunk[temp1:temp2+1]=Nogt
			dTrunk[temp1:temp2+1]=dogt
			hTrunk[temp1:temp2+1]=hogt
			NCanop[temp1:temp2+1]=Nogc
			dCanop[temp1:temp2+1]=dogc 
			hCanop[temp1:temp2+1]=hogc

		if XoMAg+XsMAg<>0: # there's a mangrove
			temp2=X[-1]-X[0]-Indexed(X,XoMAg) # landward boundary of habitat
			temp1=X[-1]-X[0]-Indexed(X,XsMAg) # shoreward boundary of habitat
			VegXlocMA[temp1:temp2+1]=1 # '1' for mangroves
			NRootsMA[temp1:temp2+1]=Nogr*dNogr
			dRootsMA[temp1:temp2+1]=dogr
			hRootsMA[temp1:temp2+1]=hogr 
			NTrunkMA[temp1:temp2+1]=Nogt*dNogt
			dTrunkMA[temp1:temp2+1]=dogt
			hTrunkMA[temp1:temp2+1]=hogt
			NCanopMA[temp1:temp2+1]=Nogc*dNogc
			dCanopMA[temp1:temp2+1]=dogc 
			hCanopMA[temp1:temp2+1]=hogc

		if XoMAg+XsMAg+Xog+Xsg<>0: # there's a mangrove
			gp.AddMessage("...a *mangrove* is present")
			if AddMActionHeader == 'yes':
				htmlfile.write("<HR><H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
			htmlfile.write("<b>A *mangrove forest* is present at your site</b>")
			if XoMAg+XsMAg==0:
				htmlfile.write("...you assumed that the whole mangrove forest will be removed.<br>")
			elif XoMAg+XsMAg<>Xog+Xsg:
				if densdeltagc+densdeltagr+densdeltagt == 0.0:
					htmlfile.write("...you assumed that the mangrove forest footprint will be reduced by "+str(abs(abs(Xog-Xsg)-abs(XoMAg-XsMAg)))+"m. <br>")
				else:
					htmlfile.write("...you assumed that the mangrove forest footprint will be reduced by "+str(abs(abs(Xog-Xsg)-abs(XoMAg-XsMAg)))+"m and that the density of the mangrove will reduce by " +str((densdeltagc+densdeltagr+densdeltagt)/3)+"%. <br>")
			elif XoMAg+XsMAg==Xog+Xsg:
				if densdeltagc+densdeltagr+densdeltagt == 0.0 :
					htmlfile.write("...you assumed that the mangrove forest will not be affected by a particular management action.<br>")
				else:
					htmlfile.write("...you assumed that the mangrove forest footprint will be unaffected but the density will reduce by " +str((densdeltagc+densdeltagr+densdeltagt)/3)+"%"+".<br>")           
		else:
			gp.AddMessage("...you do NOT have a mangrove")

		# prepare coral data
		Coral=Xco+Xcn
		Cf_coral = 0.2
		Cf_coralMA = 0.2
		if CoralType=="Barrier":
			temp=1
		elif CoralType=="Fringe":
			temp=1
		elif CoralType=="Fringe Lagoon": 
			temp=1
		else:
			Xco=0;Xcn=0.0        

		if Xco+Xcn<>0: # there's a coral reef
			gp.AddMessage("...a *coral reef* is present")
			if AddMActionHeader == 'yes':
				htmlfile.write("<HR><H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
			htmlfile.write("<b>A *coral reef* is present at your site</b>")
			htmlfile.write("...reef-top depth is " +str(hr) +"m, and reef-top width is " +str(Wr))
			if Xco+Xcn==-2: # reef is at offshore edge
				VegXloc[0]=4 # '4' for coral reefs
				HabTrigg[0]=400. # coral function will be called
				htmlfile.write("...the reef is located at the offshore edge of the profile.")
				VegXlocMA[0]=4
			else: # reef is somewhere on profile                
				if Xco-Xcn<Wr: # in case footprint is smaller than reef width; take footprint=width
					Xcn=Xco+Wr

				temp1=X[-1]-X[0]-Indexed(X,Xco) # offshore boundary of habitat
				temp2=X[-1]-X[0]-Indexed(X,Xcn) # shoreward boundary of habitat
				VegXloc[temp1:temp2+1]=4 # '4' for coral reefs
				HabTrigg[temp1:temp2+1]=400. # coral function will be calledfunction will be called
				VegXlocMA[temp1:temp2+1]=4 # '4' for coral reefs          
			Cf_coral=0.2 # friction parameter for live coral reef       
			if CoralMA=="None":
				Cf_coralMA=Cf_coral
				htmlfile.write("<li> You assumed that the reef will not be affected by a particular management action./br>")
			elif CoralMA=="Dead": # dead but structurally intact
				Cf_coralMA=0.1
				htmlfile.write("<li> You assumed that the reef will die, but will remain structurally intact.  We will assume that the coral is smooth.<br>")
			elif CoralMA=="Gone": # gone
				Cf_coralMA=0.01       
				htmlfile.write("<li> You assumed that the reef will die and will no longer remain structurally intact.  We will assume that the reef top is covered with sand.<br>")
		else:    
			gp.AddMessage("...you do NOT have a coral reef")

		HabTriggMA=num.array(HabTrigg) # Trigger for habitats is the same so far, except that we'll have to change it for oyster reef

		# prepare oyster data
		Oyster=Xr
		if OysterReefType<>"None": # there's an oyster reef    
			if OysterReefType=="Trapezoidal":
				gp.AddMessage("...a trapezoidal *oyster reef* is present")
				htmlfile.write("<b>A trapezoidal *oyster reef* is present at your site</b>")
				if AddMActionHeader == 'yes':
					htmlfile.write("<HR><H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'
			elif OysterReefType=="Dome":
				gp.AddMessage("...reef dome is present")        
				htmlfile.write("<b>A *reef dome* is present at your site</b>")
				if AddMActionHeader == 'yes':
					htmlfile.write("<HR><H2><u>Management Action</u></H2>"); AddMActionHeader = 'no'

			OysterX=X[-1]-X[0]-Indexed(X,Xr) # shoreward boundary of habitat
			VegXloc[OysterX]=5 # '5' for oyster
			HabTrigg[OysterX]=500. # oyster function will be called

			if OysterMA<>"Rmv":
				VegXlocMA[OysterX]=5 # '5' for oyster
				HabTriggMA[OysterX]=500. # oyster function will be called
				htmlfile.write("<li> Your reef will not be affected by the management action.<br>")
			else: # oyster reef is removed
				VegXlocMA[OysterX]=0 # '5' for oyster
				HabTriggMA[OysterX]=100. # oyster function will be called
				htmlfile.write("<li> You assumed that the reef will be removed.<br>")
		else:
			gp.AddMessage("...you do NOT have an oyster reef")

		# concatenate all habitat characteristics
		Roots=[NRoots,dRoots,hRoots]
		Trunk=[NTrunk,dTrunk,hTrunk]
		Canop=[NCanop,dCanop,hCanop]

		RootsMA=[NRootsMA,dRootsMA,hRootsMA]
		TrunkMA=[NTrunkMA,dTrunkMA,hTrunkMA]
		CanopMA=[NCanopMA,dCanopMA,hCanopMA]

		# locate beginning and end of each habitats for wave transformation
		Values,BegTrig,EndTrig,L_seg=ParseSegment(HabTrigg)
		ValuesMA,BegTrigMA,EndTrigMA,L_segMA=ParseSegment(HabTriggMA)
		if Values[-1]<>400 or ValuesMA[-1]<>400:
			if CoralType=="Fringe":
				CoralType="Fringe Lagoon"        
				gp.AddMessage("...You have a sandy bottom after your reef.  It's a reef with lagoon.")

	# determine what scenario we have:
	# A -> No habitats OR only a dune habitat w/o management on the dune. Only requires the simple wave model and one wave and erosion profile.
	# B -> Only a dune with a reduction in dune height. Only requires the simple wave model and one wave and two erosion profiles.
	# C -> Non-dune habitats (w/ or w/o a dune as well) w/o any management. Requires 1 run of the complex wave model and one wave and erosion profile.
	# D -> Non-dune habitats w/ a dune AND management only on the dune. Requires 1 run of the complex wave model and one wave and two erosion profiles.
	# E -> Non-dune habitats (w/ or w/o a dune as well) w/ management action on some non-dune habitat (w/ or w/o reduction in dune height). Requires two runs of the complex wave model and two wave and erosion profiles.

		if not any(VegXloc) and D1 == D2 and OysterReefType == "None" and CoralType == "None":
			scenario = "A" 
		elif not any(VegXloc) and D1 <> D2 and OysterReefType == "None" and CoralType == "None":
			scenario = "B"
		elif any(VegXloc) and list(VegXloc)==list(VegXlocMA) and D1 == D2 and Cf_coralMA==Cf_coral and dNogc+dNogr+dNogt+dNor+dNos==5 and OysterMA == "None":
			scenario = "C"       
		elif any(VegXloc) and list(VegXloc)==list(VegXlocMA) and D1 <> D2 and Cf_coralMA==Cf_coral and dNogc+dNogr+dNogt+dNor+dNos==5 and OysterMA == "None":
			scenario = "D"
		elif any(VegXloc) and list(VegXloc)<>list(VegXlocMA) or Cf_coralMA<>Cf_coral or dNogc+dNogr+dNogt+dNor+dNos<>5 or OysterMA <> "None":
			scenario = "E"

	except:
		gp.AddError(msgCreateInputs)
		raise Exception


	################################################
	####### RUN MODEL FOR MANAGEMENT ACTION ########
	################################################

	try:    
		# effect of management action
		gp.AddMessage("\nComputing wave height profiles...")

		# wave and setup at first grid point
		if Values[0]==100: # wave height at first grid point if it's not a barrier reef
			fp=1.0/To; sig=2.0*pi*fp # wave frequency and angular frequency
			ko=iterativek(sig,h[0]) # wave number at 1st grid pt
			Lo=2.0*pi/ko # wave length @ 1st grid pt
			no=0.5*(1+(2.0*ko*h[0]/sinh(2.0*ko*h[0]))) # to compute Cg at 1st grid pt
			Cgo=Lo/To*no # phase and group velocity at 1st grid pt
			# wave height at first grid point
			if Ho>0.78*(h[0]): # check that the wave height is supported by the starting water depth (according to depth limited breaking)
				Ho1=0.78*(h[0])
				gp.AddWarning("...Water depth too shallow for input wave height; that wave broke somewhere in deeper water.  We will assume that H = 0.78h");
				shal=1;
			else: # wave not broken; heck if it's in intermediate water
				shal=0
				if h[0]>0.5*Lo: 
					Ho1=Ho; # we are in deep water
				else: 
					Ho1=Ho*sqrt(0.5*g*To/(2.0*pi)/Cgo) # we are in intermediate water; assume no brkg occurred                
			# wave setup at first grid point
			Etao1=-0.125*(Ho1/sqrt(2))**2.0*ko/sinh(2.0*ko*h[0])

		else:
			shal = 0
			# barrier reef
			fp=1.0/To; sig=2.0*pi*fp # wave frequency and angular frequency
			ko=iterativek(sig,h[0]) # wave number at 1st grid pt
			Lo=2.0*pi/ko # wave length @ 1st grid pt
			no=0.5*(1+(2.0*ko*h[0]/sinh(2.0*ko*h[0]))) # to compute Cg at 1st grid pt
			Cgo=Lo/To*no # phase and group velocity at 1st grid pt
			# wave height at first grid point                
			Ho1=Ho
			Etao1=-0.125*(Ho1/sqrt(2))**2.0*ko/sinh(2.0*ko*h[0])

		Ho1MA=Ho1;EtaoMA=Etao1 # values for management option 

		# initialize variables
		H_r=[];Xn=[];ho=[];
		Retreat2 = -9999
		MErodeLen2= -9999
		EqualRetreat = "Null"
		Xaxis=[];Depth=[]; # vector of X distances and depth will be filled with coral and oyster if applicable as we move along profile
		Wave=[];WaveMA=[];Setup=[];SetupMA=[]; Dis1=[];DisSimple1=[];DisMA=[];DisSimpleMA=[] # vector of wave height and setup will be filled as we move along the profile        
		VegLoc=[];VegLocMA=[]; # vector of vegetation presence/absence will be filled
		BottVelo=[];BottVeloMA=[]; # vector of bottom velocity to be used if muddy substrate
		Hshortwave=Ho;HshortwaveMA=Ho; # short wave value to use in runup equation

		if scenario == "A":
			if D1==0:
				gp.AddMessage("...There are NO habitats on your profile.")
			else:
				gp.AddMessage("...The only habitat on your profile is a *dune*.  No management action has been taken on the dune.")
		if scenario == "B":
			gp.AddMessage("...The only habitat on your profile is a *dune*.  You have taken a management action which reduces your dune elevation.")

		if scenario == "A" or scenario == "B":
			# just need to run the simple wave model
			# Depth=num.array(h)
			# Xaxis=num.array(X)          
			# Wave,Eta=WaveModelSimple(X,h,Ho1,To,0.01)
			for rr in range(len(Values)):
				TrigVal=Values[rr];
				TrigRange=range(BegTrig[rr],EndTrig[rr]+1);
				Xtrig=X[TrigRange];htrig=h[TrigRange];

				if TrigVal==100: # bare bed or vegetated bed
					# wave model over bathy
					if Ho1>0.78*htrig[0]:
						Ho1=0.78*htrig[0]
					H,Eta,Hs,Etas,DissAtn1,Ubot=WaveModel(Xtrig,htrig,Ho1,To,Etao1,Roots,Trunk,Canop,VegXloc,TrigRange)
					Etasimple1=Etas # to modify runup  
					Wave=num.append(Wave,H);Setup=num.append(Setup,Eta);BottVelo=num.append(BottVelo,Ubot); #New wave and setup axes
					Xaxis=num.append(Xaxis,Xtrig);Depth=num.append(Depth,htrig);
					# input for next habitat
					Ho1=H[-1];
					Etao1=Eta[-1]
					EtasEnd=Etas[-1];HsEnd=Hs[-1]

				elif TrigVal==400: # coral reef
					H_r,Eta_r1=WavesCoral(Ho1,AlphF,AlphR,he,hr,Wr,Cf_coral,dx)
					if Xco-Xcn>Wr: # if the reef is trapezoidal (the foot print is bigger than the width of the reef top)
						H_reef=num.arange(0,Xco-Xcn+1)*0.0+Ho1 # initialize the wave height across the reef footprint as the incoming wave height value
						H_reef[-Wr:len(H_reef)]=H_r # place the wave height values calculated over the reef top at the end of the reef footprint
					else:
						H_reef=H_r
					Etasimple1=Eta_r1

					# input for next habitat
					Ho1=H_r[-1];
					EtasEnd=Eta_r1[-1];HsEnd=H_r[-1]
					if CoralType=="Barrier":
						Hshortwave=H_r[-1] # use the wave height value at the lagoon for run-up calculation
						Etao1=0
					elif CoralType=="Fringe":
						Hshortwave=H_r[0] # use the wave height value at the top of the reef for run-up calculation
						Etao1=Eta[-1]
					elif CoralType=="Fringe Lagoon":
						Hshortwave=H_r[-1] # use the wave height value at the at the lagoon for run-up calculation
						Etao1=Eta[-1]

					# X,h,wave and vegetation location vectors
					Wave=num.append(Wave,H_reef);Setup=num.append(Setup,Eta_r1);BottVelo=num.append(BottVelo,htrig*0); # new wave and setup axes  
					if Xco==-1 and  Xcn==-1: # reef at edge of profile
						Depth=num.append(num.arange(0,Wr)*0.0+hr,Depth)
						BottVelo=num.append(num.arange(0,Wr-1)*0.0,BottVelo)
					else:
						Xaxis=num.append(Xaxis,Xtrig)
						Depth=num.append(Depth,htrig*0.0-hr)

				elif TrigVal==500: # oyster reef
					# insert oyster reef
					TrigRange=range(BegTrig[rr+1],EndTrig[rr+1]+1);
					Kt=BreakwaterKt(Ho1,To,htrig,hc,Cw,Bw,OysterReefType);
					Ho1=float(Kt*Ho1) # transmitted wave height
					Etao1=0 
					Hshortwave=Ho1 # use the wave height value just shoreward of the reef for runup calculations
					Wave=num.append(Wave,Ho1);Setup=num.append(Setup,Etao1);BottVelo=num.append(BottVelo,0);
					Xaxis=num.append(Xaxis,Xtrig);Depth=num.append(Depth,0);

			if Xco==-1 and  Xcn==-1: # reef at edge of profile
				Xaxis=num.append(Xaxis,num.arange(max(Xaxis),max(Xaxis)+Wr)); 
			Zero=max(Xaxis); # zero is stored as the location of the shoreline at mean low water for plotting purposes.

			# estimate erosion
			if sand == 1:
				gp.AddMessage("...estimating erosion amount for sandy beach")
				# check if foreshore slope adequate (use worst wave height)
				hb=(((Ho**2.)*g*To/(2*pi))/2.)**(2./5.)/(g**(1./5.)*0.73**(4./5.)) # breaking depth
				xb=(hb/A)**1.5 # surf zone width
				ErosionTerm=xb-hb/m # erosion term in Kriebel and Dean; has to be > 0!
				if ErosionTerm<=0:
					mo=floor(xb/hb)
					gp.AddWarning("Your foreshore slope is too flat or sediment size is too high.\nWe'll increase it to 1/" +str(mo) +" to estimate a minimum erosion value.")
				else:
					mo=m

				# html file    
				htmlfile.write("<HR><H2><u>Backshore Information for Your Sandy Beach</u></H2>")
				Foreshore=str(int(Slope))
				DuneH=str(round(D1,1))
				BermH=str(round(B1,1))
				BermW=str(round(W1,1))
				htmlfile.write("<li> Your input foreshore slope is: 1/"+Foreshore+"<br>")
				if mo<>m:
					htmlfile.write("<li>Your input foreshore slope was too flat or your input sediment size was too high for us to estimate the amount of shoreline erosion at your site.  We increased your foreshore slope to 1/" +str(mo) +" to estimate a minimum erosion value.  Feel free to re-adjust your input paramters.")
				htmlfile.write("<li>The berm at your site is is: "+BermH+"m high and "+BermW+"m long<br>")
				if D1 <> 0:
					htmlfile.write("<li>The dune at your site is is: "+DuneH+"m high<br>")
				if "A":
					if D1==0:
						htmlfile.write("<li>There are no habitats on your profile.<br>")
					else:
						htmlfile.write("<li>The only habitat on your profile is a dune.  No management action has been taken on the dune.<br>")
				elif "B":
					htmlfile.write("<li>The only habitat on your profile is a dune.  You have taken a management action which reduces your dune elevation.<br>")                
				if Dred>0:
					htmlfile.write("<u>Management Action:</u>You will reduce the dune height at your site by: " +str(Dred) +"% <br>")     

				# estimate runup amount        
				Lo=g*To**2.0/(2.0*pi)

				Rnp1=1.1*(0.35*m*sqrt(Hshortwave*Lo)+sqrt(Lo*(Hshortwave*0.563*m**2.0+0.004*Ho))/2.0) # runup; short and long waves
				inundationcount = 0
				Retreat1,temp1,temp2,B_adj_pre,inundation_pre,inundationcount=ErosionKD(A,Ho,S+Rnp1,B1,D1,W1,mo,inundationcount) # erosion of beach
				if scenario == "B": # if erosion reaches the dune and the dune height is reduced, erosion must be run again to investigate the effect of the dune height reduction on the erosion distance
					Retreat2,temp1,temp2,B_adj_post,inundation_post,inundationcount=ErosionKD(A,Ho,S+Rnp1,B1,D2,W1,mo,inundationcount) # erosion of beach
				if B_adj_pre>0:
					htmlfile.write("<u>Warning</u>: Your berm elevation had to be increased by: " +str(B_adj_pre)+"m.   This means that your berm elevation was too low for your Surge and Wave Conditions to compute beach retreat.  Please check your berm elevation and forcing conditions.<br>")

			if mud == 1:
				gp.AddMessage("...estimating erosion amount for muddy substrate")
				lx=len(X);
				Retreat1,Trms1,Tc1,Tw1,Te=MudErosion(BottVelo*0,BottVelo,Depth,To,me,Cm)
				ErodeLoc=find(Trms1>Te[0]); ErodeLoc=ErodeLoc[ErodeLoc>=Zero] # indices where erosion rate greater than Threshold
				MErodeLen1=len(ErodeLoc)*dx # erosion rate greater than threshold at each location shoreward of the shoreline (pre-management)
				if any(ErodeLoc)>0:
					MErodeVol1=trapz(Retreat1[ErodeLoc]/100.0,Xaxis[ErodeLoc],dx)* StormDur # volume of mud eroded shoreward of the shoreline (m^3/m)
				else:
					MErodeVol1=0

			# flip all vectors
			Depth=Depth[::-1];Depth=Depth-MSL
			Wave=Wave[::-1];
			if mud==1:
				Depth=Depth-S
				Retreat1=Retreat1[::-1]
				Trms1=Trms1[::-1]     
			# save wave outputs
			WaveHeight=outputws+"WaveHeight_"+subwsStr+".txt"
			file=open(WaveHeight,"w")
			for i in range(0,len(Xaxis)):
				file.writelines(str(Xaxis[i])+"\t"+str(Wave[i])+"\n")
			file.close()   

		elif scenario=="C" or scenario=="D":
			if scenario == "C":
				gp.AddMessage("...There has been no managment action taken on your natural habitat.")
			if scenario == "D":
				gp.AddMessage("...The only management action was reduction of the dune elevation.  Other habitats are assumed to be unchanged.")

			# for "C" and "D" run the complex wave model (only "pre-management" required)
			for rr in range(len(Values)):
				TrigVal=Values[rr];
				TrigRange=range(BegTrig[rr],EndTrig[rr]+1);
				Xtrig=X[TrigRange];htrig=h[TrigRange];

				if TrigVal==100: # bare bed or vegetated bed
					# wave model over bathy
					if Ho1>0.78*htrig[0]:
						Ho1=0.78*htrig[0]
					H,Eta,Hs,Etas,DissAtn1,Ubot=WaveModel(Xtrig,htrig,Ho1,To,Etao1,Roots,Trunk,Canop,VegXloc,TrigRange)
					Etasimple1=Etas # to modify runup  
					Wave=num.append(Wave,H);Setup=num.append(Setup,Eta);BottVelo=num.append(BottVelo,Ubot); # new wave and setup axes
					#Dis1=num.append(Dis1,DissAtn1[0]);DisSimple1=num.append(DisSimple1,DissAtn1[1])
					temp1=[H[ii]**3 for ii in range(len(H))];Dis1=num.append(Dis1,temp1);
					temp2=[Hs[ii]**3 for ii in range(len(Hs))];DisSimple1=num.append(DisSimple1,temp2)
					Xaxis=num.append(Xaxis,Xtrig);Depth=num.append(Depth,htrig);
					# input for next habitat
					Ho1=H[-1];
					Etao1=Eta[-1]
					EtasEnd=Etas[-1];HsEnd=Hs[-1]

				elif TrigVal==400: # coral reef
					H_r,Eta_r1=WavesCoral(Ho1,AlphF,AlphR,he,hr,Wr,Cf_coral,dx)
					ECg=[-0.125*1024*9.81*H_r[ii]**2*sqrt(9.81*(hr+Eta_r1[-1])) for ii in range(len(H_r))]
					Diss=gradient(ECg,dx);Diss=H_r[1:-1]**3;Dis1=num.append(Dis1,Diss);

					H_temp,temp=WavesCoral(Ho1,AlphF,AlphR,he,hr,Wr,0.01,dx)
					ECg=[-0.125*1024*9.81*H_temp[ii]**2*sqrt(9.81*(hr+temp[-1])) for ii in range(len(H_temp))]
					Diss=gradient(ECg,dx);Diss=H_temp[1:-1]**3;DisSimple1=num.append(DisSimple1,Diss);

					if Xco-Xcn>Wr: # if the reef is trapezoidal (the foot print is bigger than the width of the reef top)
						H_reef=num.arange(0,Xco-Xcn+1)*0.0+Ho1 # initialize the wave height across the reef footprint as the incoming wave height value
						H_reef[-Wr:len(H_reef)]=H_r # place the wave height values calculated over the reef top at the end of the reef footprint
					else:
						H_reef=H_r
					Etasimple1=Eta_r1

					# input for next habitat
					Ho1=H_r[-1];
					EtasEnd=Eta_r1[-1];HsEnd=H_r[-1]
					if CoralType=="Barrier":
						Hshortwave=H_r[-1] # use the wave height value at the lagoon for run-up calculation
						Etao1=0
					elif CoralType=="Fringe":
						Hshortwave=H_r[0] # use the wave height value at the top of the reef for run-up calculation
						Etao1=Eta[-1]
					elif CoralType=="Fringe Lagoon":
						Hshortwave=H_r[-1] # use the wave height value at the at the lagoon for run-up calculation
						Etao1=Eta[-1]

					# X,h,wave and vegetation location vectors
					Wave=num.append(Wave,H_reef);Setup=num.append(Setup,Eta_r1);BottVelo=num.append(BottVelo,htrig*0);#New wave and setup axes  
					if Xco==-1 and  Xcn==-1: #Reef at edge of profile
						Depth=num.append(num.arange(0,Wr)*0.0+hr,Depth)
						BottVelo=num.append(num.arange(0,Wr-1)*0.0,BottVelo)
					else:
						Xaxis=num.append(Xaxis,Xtrig)
						Depth=num.append(Depth,htrig*0.0-hr)

				elif TrigVal==500: # oyster reef
					# insert oyster reef
					TrigRange=range(BegTrig[rr+1],EndTrig[rr+1]+1);
					Kt=BreakwaterKt(Ho1,To,htrig,hc,Cw,Bw,OysterReefType);
					Ho1=float(Kt*Ho1) # transmitted Wave height
					Etao1=0 
					Hshortwave=Ho1 # use the wave height value just shoreward of the reef for run-up calculations
					Wave=num.append(Wave,Ho1);Setup=num.append(Setup,Etao1);BottVelo=num.append(BottVelo,0);
					temp=Ho1/Kt;DisSimple1=num.append(DisSimple1,temp**3);Dis1=num.append(Dis1,Ho1**3)
					Xaxis=num.append(Xaxis,Xtrig);Depth=num.append(Depth,0);
			if Xco==-1 and  Xcn==-1: # reef at edge of profile
				Xaxis=num.append(Xaxis,num.arange(max(Xaxis),max(Xaxis)+Wr)); 

			Zero=max(Xaxis); # zero is stored as the location of the shoreline at mean low water for plotting purposes          
			temp=[Dis1[ii]/DisSimple1[ii] for ii in range(len(Dis1))]
			#for ii in range(len(Dis1)):
				#if Dis1[ii]==DisSimple1[ii]:
					#temp[ii]=NaN
			Diss=[]
			for ii in range(len(temp)):
				if isnan(temp[ii])<>1:
					Diss=num.append(Diss,temp[ii])
			if len(Diss)<>0:
				DissBF=mean(Diss)
			else:
				DissBF=-1                

			if sand == 1:
				gp.AddMessage("...estimating erosion amount for sandy beach")
				# check if foreshore slope adequate (use worst wave height)
				hb=(((Ho**2.)*g*To/(2*pi))/2.)**(2./5.)/(g**(1./5.)*0.73**(4./5.)) # breaking depth
				xb=(hb/A)**1.5 # surf zone width
				ErosionTerm=xb-hb/m # erosion term in Kriebel and Dean; has to be > 0!
				if ErosionTerm<=0:
					mo=floor(xb/hb)
					gp.AddWarning("Your foreshore slope is too flat or sediment size is too high.\nWe'll increase it to 1/" +str(mo) +" to estimate a minimum erosion value.")
				else:
					mo=m

				# html file    
				htmlfile.write("<HR><H2><u>Backshore Information for Your Sandy Beach</u></H2>")
				Foreshore=str(int(Slope))
				DuneH=str(round(D1,1))
				BermH=str(round(B1,1))
				BermW=str(round(W1,1))
				htmlfile.write("<li> Your input foreshore slope is: 1/"+Foreshore+"<br>")
				if mo<>m:
					htmlfile.write("<li>Your input foreshore slope was too flat or your input sediment size was too high for us to estimate the amount of shoreline erosion at your site.  We increased your foreshore slope to 1/" +str(mo) +" to estimate a minimum erosion value.  Feel free to re-adjust your input paramters.")
				htmlfile.write("<li>The berm at your site is is: "+BermH+"m high and "+BermW+"m long<br>")
				if D1 <> 0:
					htmlfile.write("<li>The dune at your site is is: "+DuneH+"m high<br>")
				if "C":
					htmlfile.write("<li>There has been no managment action taken on your natural habitat.<br>")
				elif "D":
					htmlfile.write("<li>The only management action taken is a reduction in dune elevation.  All offshore habitats are assumed to be unaffected.<br>")                                  
				if Dred>0:
					htmlfile.write("<u>Management Action:</u>You will reduce the dune height at your site by: " +str(Dred) +"% <br>")

				# estimate runup amount        
				Lo=g*To**2.0/(2.0*pi)

				Rnp1=1.1*(0.35*m*sqrt(Hshortwave*Lo)+sqrt(Lo*(Hshortwave*0.563*m**2.0+0.004*Ho))/2.0) # Runup -short and long waves w/o vegetation effects 
				inundationcount = 0
				R1,temp1,temp2,B_adj_pre,inundation_pre,inundationcount=ErosionKD(A,Ho,S+Rnp1,B1,D1,W1,mo,inundationcount) # erosion without vegetation

				TrigRange=range(BegTrig[-1],EndTrig[-1]+1); # last portion of the profile
				#if Values[-1]==100 and any(VegXloc[TrigRange]>0): # last portion of the profile is not reef and is vegetated
				Emax=0.35*m*num.sqrt(Hshortwave*Lo) # setup at the beach
				coef0=max(Etasimple1[-len(TrigRange):-1])/Emax; # correction factor for MWL
				Etap=max(Setup[-len(TrigRange):-1])/coef0 # corrected MWL at shoreline in presence of habitat
				if Etap<0:
					Etap=0 # if Eta with veg. neg,take as zero
				Hp=(Etap/(0.35*m))**2/Lo # Hprime to estimate runup with veg
				Rnpveg=1.1*(Etap+num.sqrt(Lo*(Hp*0.563*m**2+Ho*0.004))/2) # runup with vegetation
				R_rnp=R1*Rnpveg/Rnp1 # scale beach retreat by runup
				if DissBF<>-1: # wave model did not crash
					R_dissip=R1*DissBF # scale beach retreat by dissipation due to vegetation
				else:
					R_dissip=R_rnp # scale beach retreat by dissipation due to vegetation
				Retreat1=mean([R_rnp,R_dissip]) # take the mean of both approaches  

				if scenario == "D": # if the dune height is reduced, erosion must be run again to investigate the effect of the dune height reduction on the erosion distance
					R2,temp1,temp2,B_adj_post,inundation_post,inundationcount=ErosionKD(A,Ho,Rnp1+S,B1,D2,W1,mo,inundationcount) # erosion of after dune reduction and without vegetation beach
					R_rnp=R2*Rnpveg/Rnp1 # scale beach retreat by runup
					if DissBF<>-1: # wave model did not crash
						R_dissip=R2*DissBF # scale beach retreat by dissipation due to vegetation
					else:
						R_dissip=R_rnp # scale beach retreat by dissipation due to vegetation
					Retreat2=mean([R_rnp,R_dissip]) # take the mean of both approaches    

				if B_adj_pre>0:
					htmlfile.write("<u>Warning</u>: Your berm elevation had to be increased by: " +str(B_adj_pre)+"m.  This means that your berm elevation was too low for your surge and wave conditions to compute beach retreat.  Please check your berm elevation and forcing conditions.<br>")

			if mud == 1:
				gp.AddMessage("...estimating erosion amount for muddy substrate")
				lx=len(X);
				Retreat1,Trms1,Tc1,Tw1,Te=MudErosion(BottVelo*0,BottVelo,Depth,To,me,Cm) # before management action
				ErodeLoc=find(Trms1>Te[0]); ErodeLoc=ErodeLoc[ErodeLoc>=Zero]# Indices where erosion rate greater than Threshold
				MErodeLen1=len(ErodeLoc)*dx # erosion rate greater than threshold at each location shoreward of the shoreline (pre-managament)
				if any(ErodeLoc)>0:
					MErodeVol1=trapz(Retreat1[ErodeLoc]/100.0,Xaxis[ErodeLoc],dx)* StormDur # volume of mud eroded shoreward of the shoreline (m^3/m)
				else:
					MErodeVol1=0

			# flip all vectors
			Depth=Depth[::-1];Depth=Depth-MSL
			Wave=Wave[::-1];
			VegXloc=VegXloc[::-1]
			if mud==1:
				Depth=Depth-S
				Retreat1=Retreat1[::-1]
				Trms1=Trms1[::-1]     
			# save wave outputs
			WaveHeight=outputws+"WaveHeight_"+subwsStr+".txt"
			file=open(WaveHeight,"w")
			for i in range(0,len(Wave)):
				file.writelines(str(Xaxis[i])+"\t"+str(Wave[i])+"\n")
			file.close()  

		elif scenario=="E":
			# before management action
			for rr in range(len(Values)):
				TrigVal=Values[rr];
				TrigRange=range(BegTrig[rr],EndTrig[rr]+1);
				Xtrig=X[TrigRange];htrig=h[TrigRange];

				if TrigVal==100: # bare bed or vegetated bed
					# wave model over bathy
					H,Eta,Hs,Etas,DissAtn1,Ubot=WaveModel(Xtrig,htrig,Ho1,To,Etao1,Roots,Trunk,Canop,VegXloc,TrigRange)
					Etasimple1=Etas # to modify runup  
					Wave=num.append(Wave,H);Setup=num.append(Setup,Eta);BottVelo=num.append(BottVelo,Ubot);
					#Dis1=num.append(Dis1,DissAtn1[0]);DisSimple1=num.append(DisSimple1,DissAtn1[1])
					temp1=[H[ii]**3 for ii in range(len(H))];Dis1=num.append(Dis1,temp1);
					temp2=[Hs[ii]**3 for ii in range(len(Hs))];DisSimple1=num.append(DisSimple1,temp2)
					# new wave and setup axes
					Xaxis=num.append(Xaxis,Xtrig);Depth=num.append(Depth,htrig);
					# input for next habitat
					Ho1=H[-1];Etao1=Eta[-1]
					EtasEnd=Etas[-1];HsEnd=Hs[-1]

				elif TrigVal==400: # coral reef
					H_r,Eta_r1=WavesCoral(Ho1,AlphF,AlphR,he,hr,Wr,Cf_coral,dx)
					ECg=[-0.125*1024*9.81*H_r[ii]**2*sqrt(9.81*(hr+Eta_r1[-1])) for ii in range(len(H_r))]
					Diss=gradient(ECg,dx);Diss=H_r[1:-1]**3;Dis1=num.append(Dis1,Diss);

					H_temp,temp=WavesCoral(Ho1,AlphF,AlphR,he,hr,Wr,0.01,dx)
					ECg=[-0.125*1024*9.81*H_temp[ii]**2*sqrt(9.81*(hr+temp[-1])) for ii in range(len(H_temp))]
					Diss=gradient(ECg,dx);Diss=H_temp[1:-1]**3;DisSimple1=num.append(DisSimple1,Diss);

					if Xco-Xcn>Wr: # if the reef is trapezoidal (the foot print is bigger than the width of the reef top)
						H_reef=num.arange(0,Xco-Xcn+1)*0.0+Ho1 # initialize the wave height across the reef footprint as the incoming wave height value
						H_reef[-Wr:len(H_reef)]=H_r # place the wave height values calculated over the reef top at the end of the reef footprint
					else:
						H_reef=H_r
					Etasimple1=Eta_r1

					# input for next habitat
					Ho1=H_r[-1];Etao1=Eta_r1[-1]
					EtasEnd=Eta_r1[-1];HsEnd=H_r[-1]
					if CoralType=="Barrier":
						Hshortwave=H_r[-1] # use the wave height value at the lagoon for run-up calculation
						Etao1=0
					elif CoralType=="Fringe":
						Hshortwave=H_r[0] # use the wave height value at the top of the reef for run-up calculation
					elif CoralType=="Fringe Lagoon":
						Hshortwave=H_r[-1] # use the wave height value at the at the lagoon for run-up calculation

					Wave=num.append(Wave,H_reef);Setup=num.append(Setup,Eta_r1);BottVelo=num.append(BottVelo,htrig*0); # new wave and setup axes  
					if Xco==-1 and  Xcn==-1: # reef at edge of profile                  
						Depth=num.append(num.arange(0,Wr)*0.0+hr,Depth)    
						BottVelo=num.append(num.arange(0,Wr-1)*0.0,BottVelo)
					else:
						Xaxis=num.append(Xaxis,Xtrig)
						Depth=num.append(Depth,htrig*0.0-hr)

				elif TrigVal==500: # oyster reef
					# insert oyster reef
					TrigRange=range(BegTrig[rr+1],EndTrig[rr+1]+1);
					Kt=BreakwaterKt(Ho1,To,htrig,hc,Cw,Bw,OysterReefType); 
					Ho1=float(Kt*Ho1) # transmitted wave height
					Etao1=0 
					Hshortwave=Ho1 # use the wave height value just shoreward of the reef for runup calculations
					Wave=num.append(Wave,Ho1);Setup=num.append(Setup,Etao1);BottVelo=num.append(BottVelo,0);
					temp=Ho1/Kt;DisSimple1=num.append(DisSimple1,temp**3);Dis1=num.append(Dis1,Ho1**3)
					Xaxis=num.append(Xaxis,Xtrig);Depth=num.append(Depth,0);

			if Xco==-1 and  Xcn==-1: # reef at edge of profile
				Xaxis=num.append(Xaxis,num.arange(max(Xaxis),max(Xaxis)+Wr));        

			temp=[Dis1[ii]/DisSimple1[ii] for ii in range(len(Dis1))]
			#for ii in range(len(Dis1)):
				#if Dis1[ii]==DisSimple1[ii]:
					#temp[ii]=NaN
			Diss=[]
			for ii in range(len(temp)):
				if isnan(temp[ii])<>1:
					Diss=num.append(Diss,temp[ii])
			if len(Diss)<>0:
				DissBF=mean(Diss)
			else:
				DissBF=-1                

			gp.AddMessage("\nRunning Wave Model after habitat modification...")            
			# after management action:
			for rr in range(len(ValuesMA)):
				TrigVal=ValuesMA[rr];
				TrigRange=range(BegTrigMA[rr],EndTrigMA[rr]+1);
				Xtrig=X[TrigRange];htrig=h[TrigRange];

				if TrigVal==100: # bare bed or vegetated bed
					# wave model over bathy
					H,Eta,Hs,Etas,DissAtnMA,Ubot=WaveModel(Xtrig,htrig,Ho1MA,To,EtaoMA,RootsMA,TrunkMA,CanopMA,VegXlocMA,TrigRange)
					EtasimpleMA=Etas # to modify runup  
					Ho1MA=H[-1];EtaoMA=Eta[-1] # input for next habitat
					WaveMA=num.append(WaveMA,H);SetupMA=num.append(SetupMA,Eta); # new wave and setup axes
					BottVeloMA=num.append(BottVeloMA,Ubot);  
					#DisMA=num.append(DisMA,DissAtnMA[0]);DisSimpleMA=num.append(DisSimpleMA,DissAtnMA[1])
					temp1=[H[ii]**3 for ii in range(len(H))];DisMA=num.append(DisMA,temp1);
					temp2=[Hs[ii]**3 for ii in range(len(Hs))];DisSimpleMA=num.append(DisSimpleMA,temp2)

				elif TrigVal==400: # coral reef
					H_r,Eta_rMA=WavesCoral(Ho1MA,AlphF,AlphR,he,hr,Wr,Cf_coralMA,dx)
					ECg=[-0.125*1024*9.81*H_r[ii]**2*sqrt(9.81*(hr+Eta_rMA[-1])) for ii in range(len(H_r))]
					Diss=gradient(ECg,dx);Diss=H_r[1:-1]**3;DisMA=num.append(DisMA,Diss);

					H_temp,temp=WavesCoral(Ho1MA,AlphF,AlphR,he,hr,Wr,0.01,dx)
					ECg=[-0.125*1024*9.81*H_temp[ii]**2*sqrt(9.81*(hr+temp[-1])) for ii in range(len(H_temp))]
					Diss=gradient(ECg,dx);Diss=H_temp[1:-1]**3;DisSimpleMA=num.append(DisSimpleMA,Diss);

					if Xco-Xcn>Wr:
						H_reef=num.arange(0,Xco-Xcn+1)*0.0+Ho1MA
						H_reef[-Wr:len(H_reef)]=H_r
					else:
						H_reef=H_r

					EtasimpleMA=Eta_rMA          
					Ho1MA=H_r[-1];EtaoMA=Eta_rMA[-1] # input for next habitat
					if CoralType=="Barrier":
						HshortwaveMA=H_r[-1] # use the wave height value at the lagoon for run-up calculation
						EtaoMA=0
					elif CoralType=="Fringe":
						HshortwaveMA=H_r[0] # use the wave height value at the top of the reef for run-up calculation
					elif CoralType=="Fringe Lagoon":
						HshortwaveMA=H_r[-1] # use the wave height value at the at the lagoon for run-up calculation

					WaveMA=num.append(WaveMA,H_reef);SetupMA=num.append(SetupMA,Eta_rMA); # new wave and setup axes  
					BottVeloMA=num.append(BottVeloMA,htrig*0);   
					if Xco==-1 and  Xcn==-1: # reef at edge of profile
						BottVeloMA=num.append(num.arange(0,Wr-1)*0.0,BottVeloMA)                                    

				elif TrigVal==500: # oyster reef
					# insert oyster reef
					Kt=BreakwaterKt(Ho1MA,To,htrig,hc,Cw,Bw,OysterReefType);
					Ho1MA=float(Kt*Ho1MA) # transmitted wave height
					EtaoMA=0 
					HshortwaveMA=Ho1MA # wave in lagoon shoreward of the reef
					WaveMA=num.append(WaveMA,Ho1MA);SetupMA=num.append(SetupMA,EtaoMA);
					BottVeloMA=num.append(BottVeloMA,0);
					temp=Ho1/Kt;DisSimpleMA=num.append(DisSimpleMA,temp**3);DisMA=num.append(DisMA,Ho1**3)

			temp=[DisMA[ii]/DisSimpleMA[ii] for ii in range(len(DisMA))]
			#for ii in range(len(DisMA)):
				#if DisMA[ii]==DisSimpleMA[ii]:
					#temp[ii]=NaN
			Diss=[]
			for ii in range(len(temp)):
				if isnan(temp[ii])<>1:
					Diss=num.append(Diss,temp[ii])
			if len(Diss)<>0:
				DissMA=mean(Diss)
			else:
				DissMA=-1  
			if DissMA<DissBF:
				DissBF=DissMA
				#gp.AddWarning("The wave dissapation pre and post management action was set equal due to a model error. There likely is more erosion post management than actually computed. Please use results with caution.")
			

			Zero=max(Xaxis); # zero is stored as the location of the shoreline at mean low water for plotting purposes
			# estimate runup amount        
			Lo=g*To**2.0/(2.0*pi)
			# before management action
			Rnp1=1.1*(0.35*m*sqrt(Hshortwave*Lo)+sqrt(Lo*(Hshortwave*0.563*m**2.0+Ho*0.004))/2.0) # run-up - short and long waves (run-up w/o vegetation effect)

			TrigRange=range(BegTrig[-1],EndTrig[-1]+1); # last portion of the profile
			if Values[-1]==100 : # last portion of the profile is not reef 
				if CoralType<>"Fringe":
					Emax=0.35*m*num.sqrt(Hshortwave*Lo) # setup at the beach
					coef0=max(Etasimple1[-len(TrigRange):-1])/Emax; # correction factor for MWL
					Etap=max(Setup[-len(TrigRange):-1])/coef0 # corrected MWL at shoreline in presence of habitat
				else:
					Etap=Eta_r1[-1]

				if Etap<0:
					Etap=0 # if Eta with vegetation is negative, take as zero
				Hp=(Etap/(0.35*m))**2/Lo # Hprime to estimate run-up with vegetation
				Rnpveg1=1.1*(Etap+num.sqrt(Lo*(Hp*0.563*m**2+Ho*0.004))/2.0) # run-up with vegetation

			# after management action
			RmaxMA=1.1*(0.35*m*sqrt(HshortwaveMA*Lo)+sqrt(Lo*(HshortwaveMA*0.563*m**2.0+Ho*0.004))/2.) # run-up - short and long waves if reefs are located on your profile, otherwise RmaxMA = Rnp (RnpMA = Rnp1)
			RnpMA=RmaxMA # run-up w/o vegetation effect

			TrigRange=range(BegTrigMA[-1],EndTrigMA[-1]+1); # last portion of the profile
			if ValuesMA[-1]==100 : # last portion of the profile is not reef 
				if CoralType<>"Fringe":
					Emax=0.35*m*num.sqrt(HshortwaveMA*Lo) # setup at the beach
					coef0=max(EtasimpleMA[-len(TrigRange):-1])/Emax; # correction factor for MWL
					Etap=max(SetupMA[-len(TrigRange):-1])/coef0 # corrected MWL at shoreline in presence of habitat
				else:
					Etap=Eta_rMA[-1]

				if Etap<0:
					Etap=0 # if Eta with veg. neg,take as zero
				Hp=(Etap/(0.35*m))**2/Lo # Hprime to estimate run-up with vegetation
				RnpvegMA=1.1*(Etap+num.sqrt(Lo*(Hp*0.563*m**2+Ho*0.004))/2) # run-up with vegetation    

			# estimate erosion
			if sand==1:
				gp.AddMessage("...estimating erosion amount for sandy beach")
				# check if foreshore slope adequate (use worst wave height)
				hb=(((Ho**2.)*g*To/(2*pi))/2.)**(2./5.)/(g**(1./5.)*0.73**(4./5.)) # breaking depth
				xb=(hb/A)**1.5 # surf zone width
				ErosionTerm=xb-hb/m # erosion term in Kriebel and Dean; has to be > 0!
				if ErosionTerm<=0:
					mo=floor(xb/hb)
					gp.AddWarning("Your foreshore slope is too flat or sediment size is too high.\nWe'll increase it to 1/" +str(mo) +" to estimate a minimum erosion value.")
				else:
					mo=m

				# html file    
				htmlfile.write("<HR><H2><u>Backshore Information for Your Sandy Beach</u></H2>")
				Foreshore=str(int(Slope))
				DuneH=str(round(D1,1))
				BermH=str(round(B1,1))
				BermW=str(round(W1,1))
				htmlfile.write("<li> Your input foreshore slope is: 1/"+Foreshore+"<br>")
				if mo<>m:
					htmlfile.write("<li>Your input foreshore slope was too flat or your input sediment size was too high for us to estimate the amount of shoreline erosion at your site.  We increased your foreshore slope to 1/" +str(mo) +" to estimate a minimum erosion value.  Feel free to re-adjust your input paramters.")
				htmlfile.write("<li>The berm at your site is is: "+BermH+"m high and "+BermW+"m long<br>")
				if D1 <> 0:
					htmlfile.write("<li>The dune at your site is is: "+DuneH+"m high<br>")
				if Dred>0:
					htmlfile.write("<u>Management Action:</u>You will reduce the dune height at your site by: " +str(Dred) +"% <br>")

				# before management action
				inundationcount = 0 # tracks whether or not the profile was inundation pre MA to avoid redundant messages.

				R1,temp1,temp2,B_adj_pre,inundation_pre,inundationcount=ErosionKD(A,Ho,Rnp1+S,B1,D1,W1,mo,inundationcount) # erosion of beach 
				if Rnpveg1>Rnp1:
					Rnpveg1=Rnp1
					R_rnp=R1
				else:
					R_rnp=R1*Rnpveg1/Rnp1 # scale beach retreat by runup change due to vegetation

				if DissBF<>-1: # wave model did not crash
					R_dissip=R1*DissBF # scale beach retreat by dissipation due to vegetation
				else:
					R_dissip=R_rnp # scale beach retreat by dissipation due to vegetation

				temp1=mean([R_rnp,R_dissip]) # take the mean of both approaches
				temp2=mean([R_rnp]) # take the mean of both approaches

				# after manamgement action
				if HshortwaveMA+S+D2==Hshortwave+S+D1:
					R2 = R1 # if same forcing, no need to re-run
					B_adj_post = B_adj_pre
					inundation_post = inundation_pre
				else:
					R2,temp1,temp2,B_adj_post,inundation_post,inundationcount=ErosionKD(A,Ho,RnpMA+S,B1,D2,W1,mo,inundationcount) # erosion of beach

				if RnpvegMA>RnpMA:
					RnpvegMA=RnpMA                        
				R_rnpMA=R2*RnpvegMA/RmaxMA # scale by runup
				if DissMA<>-1:
					R_dissipMA=R2*DissMA # scale by dissipation
				else:
					R_dissipMA=R_rnpMA # scale by dissipation
				temp1MA=mean([R_rnpMA,R_dissipMA]) # take the mean of both approaches       
				temp2MA=mean([R_rnpMA]) # take the mean of both approaches  
				
				d1=abs(temp1-temp1MA);d2=abs(temp2-temp2MA);
				if d1>d2:
					Retreat1=temp1;Retreat2=temp1MA
				else:
					Retreat1=temp2;Retreat2=temp2MA

				if B_adj_pre>0:
					htmlfile.write("<u>Warning</u>: Your berm elevation had to be increased by: " +str(B_adj_pre)+"m.   This means that your berm elevation was too low for your Surge and Wave Conditions to compute beach retreat.  Please check your berm elevation and forcing conditions.<br>")
				elif B_adj_post>0:
					htmlfile.write("<u>Warning</u>: Your berm elevation had to be increased by: " +str(B_adj_post)+"m <b> after you removed habitat </b>.   This means that your berm elevation was too low for your Surge and Wave Conditions to compute beach retreat after your managment action.  Please check your berm elevation and forcing conditions.<br>")
				else:
					gp.AddMessage("...Berm elevation is appropriate for forcing condition.")           

			elif mud==1: # compute erosion amount for consolidated sediments
				gp.AddMessage("...estimating erosion amount for muddy substrate")
				lx=len(Xaxis);msg=0
				#values=range(Zero,lx)
				Retreat1,Trms1,Tc1,Tw1,Te=MudErosion(BottVelo*0,BottVelo,Depth,To,me,Cm) # before management action
				Retreat2,Trms2,Tc2,Tw2,Te=MudErosion(BottVeloMA*0,BottVeloMA,Depth,To,me,Cm) # after                        
				ErodeLoc=find(Trms1>Te[0]); ErodeLoc=ErodeLoc[ErodeLoc>=Zero]# Indices where erosion rate greater than Threshold
				MErodeLen1=len(ErodeLoc)*dx # Erosion rate greater than Threshold at each location shoreward of the shoreline (pre managament)
				if any(ErodeLoc)>0:
					MErodeVol1=trapz(Retreat1[ErodeLoc]/100.0,Xaxis[ErodeLoc],dx)* StormDur #Volume of mud eroded shoreward of the shoreline (m^3/m)
				else:
					MErodeVol1=0

				ErodeLoc=find(Trms2>Te[0]); ErodeLoc=ErodeLoc[ErodeLoc>=Zero]# Indices where erosion rate greater than Threshold
				MErodeLen2=len(ErodeLoc)*dx # Erosion rate greater than Threshold at each location shoreward of the shoreline (post managament)
				if any(ErodeLoc)>0:
					MErodeVol2=trapz(Retreat2[ErodeLoc]/100.0,Xaxis[ErodeLoc],dx)* StormDur #Volume of mud eroded shoreward of the shoreline (m^3/m)
				else:
					MErodeVol2=0

				gp.AddMessage("......Erosion Length (Pre-MA): " +str(round(MErodeLen1,1)) + "\n......Erosion Length (Post-MA): " + str(round(MErodeLen2,1)) + "\n......Erosion Volume (Pre-MA): " +str(round(MErodeVol1)) + "\n......Erosion Volume (Post-MA): " + str(round(MErodeVol2)))

			# flip all vectors
			Depth=Depth[::-1];Depth=Depth-MSL
			Wave=Wave[::-1];WaveMA=WaveMA[::-1];
			VegXloc=VegXloc[::-1];VegXlocMA=VegXlocMA[::-1];
			if mud==1:
				Depth=Depth-S
				Retreat1=Retreat1[::-1];Retreat2=Retreat2[::-1];
				Trms2=Trms2[::-1];Trms1=Trms1[::-1]

			# save wave outputs
			WaveHeight=outputws+"WaveHeight_"+subwsStr+".txt"
			file=open(WaveHeight,"w")
			for i in range(0,min(len(Wave),len(WaveMA),len(Xaxis))):
				file.writelines(str(Xaxis[i])+"\t"+str(Wave[i])+"\t"+str(WaveMA[i])+"\n")
			file.close()      

	except:
		gp.AddError(msgManagementAction)
		raise Exception


	##############################
	####### CREATE OUTPUTS #######
	##############################

	try:
		gp.AddMessage("\nPlotting wave profiles...")
		from matplotlib.font_manager import FontProperties
		fontP = FontProperties()
		fontP.set_size('small')
		# percent wave attenuation
		if scenario == "E": # only if management actions were performed on non-dune habitats will the wave height profile change ("scenario" E)
			AtnH=Wave*0;AtnE=Wave*0;
			AtnH=abs(WaveMA-Wave)/WaveMA*100 # percent attenuation
			WAtn=mean(AtnH[AtnH>0.01]) # average wave attenuation
			AtnE=abs(WaveMA**2-Wave**2)/WaveMA**2*100 # percent attenuation
			EAtn=mean(AtnE[AtnE>0.01]) # average energy attenuation
		elif list(VegXloc)<>list(VegXlocMA) and OysterMA<>"None" and CoralMA<>"None":
			gp.AddWarning("There was a management action taken along your profile, but it was not captured by the wave model.")

		# plots of wave height
		figure(1)
		if scenario == "E": # only if management actions were performed on non-dune habitats will the wave height profile change ("scenario" E)
			ax=subplot(211);plot(Xaxis,WaveMA, 'r-', Xaxis,Wave, 'g.', linewidth=2);grid()
			box=ax.get_position();
			ax.set_position([box.x0, box.y0, box.width*1, box.height])
			ylabel('Wave Height[m]',size='large')
			xlim(X[-1],X[0])
			ax.legend(('After Mgmt','Before Mgmt'),prop=fontP,loc='upper left',ncol=2)
		else:
			ax=subplot(211);plot(Xaxis,Wave, 'g', linewidth=2);grid()
			box=ax.get_position();
			ax.set_position([box.x0, box.y0, box.width*1, box.height])
			ylabel('Wave Height[m]',size='large')
			xlim(X[-1],X[0])
			ax.legend(('Wave Ht. Prof.'),prop=fontP,loc='upper left',ncol=2)  

		ax=subplot(212);plot(Xaxis,-Depth,linewidth=2);grid()  
		xlim(X[-1],X[0])
		vegleg=['Bed']  
		if scenario == "E" or scenario == "D" or scenario == "C": # only if there are non-dune habitats
			if any(VegXloc==1): # if mangrove is present
				veg=Depth*0.+NaN;veg[VegXloc==1]=-Depth[VegXloc==1]
				plot(Xaxis,veg,'.r',linewidth=1)
				vegleg.append('Mangrove')
			if any(VegXloc==2): # if marsh is present
				veg=Depth*0.+NaN;veg[VegXloc==2]=-Depth[VegXloc==2]
				plot(Xaxis,veg,'xr',linewidth=1)
				vegleg.append('Marsh')
			if any(VegXloc==3): # if Seagrass is present
				veg=Depth*0.+NaN;veg[VegXloc==3]=-Depth[VegXloc==3]
				plot(Xaxis,veg,'+g',linewidth=1)
				vegleg.append('Seagrass')
			if any(VegXloc==4): # if coral is present
				veg=Depth*0.+NaN;veg[VegXloc==4]=-Depth[VegXloc==4]
				plot(Xaxis,veg,'og',linewidth=1)
				vegleg.append('Coral Reef')
			if any(VegXloc==5): # if an oyster reef is present
				hi=-Depth[VegXloc==5]
				Yrf=num.arange((hi),0.0,0.05);Xrf=(Yrf*0.0+Xr) # x-axis
				plot(Xrf,Yrf,'g',Xrf+.05,Yrf,'g',Xrf+.1,Yrf,'g',linewidth=2);grid()
				vegleg.append('Oyster Reef')
		box=ax.get_position();
		ax.set_position([box.x0, box.y0, box.width*1, box.height])                
		xlabel('Cross-Shore Distance from Shoreline',size='large')
		ax.legend((vegleg),prop=fontP,loc='upper left',ncol=2)
		axvline(x=0, linewidth=1, color='k')
		annotate('Shoreline',xy=(1000,0))
		ylabel('Depth[m]',size='large')
		savefig(outputws+"WavePlot_"+subwsStr+".png",dpi=(640/8))

		figure(2) # zoom in of figure 1
		ax=subplot(211)
		vegleg=['Bed'] 
		if scenario == "D" or scenario == "E" or scenario == "C":    
			VegEnd=find(VegXloc>0); VegEnd=min(VegEnd[-1]+300,len(Wave)) # location of most offshore vegetation stem
		else:
			VegEnd=min(1000,len(Wave))
		if scenario == "E": # only if some management action has been taken to effect the wave height profile. ("Scenario" E)
			plot(Xaxis,WaveMA, 'r-', Xaxis,Wave, 'g.', linewidth=2);grid()
			box=ax.get_position();
			ax.set_position([box.x0, box.y0, box.width*1, box.height])
			ylabel('Wave Height[m]',size='large')
			ax.legend(('After Mgmt','Before Mgmt'),prop=fontP,loc='upper left',ncol=2)
			xlim(VegEnd+100,Xaxis[0]);
			if VegEnd == len (Wave):
				xlim(Xaxis[-1],Xaxis[0]);  
				ylim(0,Wave[-1]+.3)                
			else:
				xlim(Xaxis[VegEnd],Xaxis[0]);  
				ylim(0,Wave[VegEnd]+.3)          
		else: # if no management action was taken on non-dune habitat
			plot(Xaxis,Wave, 'g', linewidth=2);grid()
			box=ax.get_position();
			ax.set_position([box.x0, box.y0, box.width*1, box.height])
			ylabel('Wave Height[m]',size='large')
			ax.legend(('Wave Ht. Prof.'),prop=fontP,loc='upper left',ncol=2)
			if VegEnd == len (Wave):
				xlim(Xaxis[-1],Xaxis[0]);  
				ylim(0,Wave[-1]+.3)                
			else:
				xlim(Xaxis[VegEnd],Xaxis[0]);  
				ylim(0,Wave[VegEnd]+.3)
		ax=subplot(212)
		plot(Xaxis,-Depth,linewidth=2);grid()  
		if any(VegXloc==1): # if mangrove is present
			veg=Depth*0.+NaN;veg[VegXloc==1]=-Depth[VegXloc==1]
			plot(Xaxis,veg,'.r',linewidth=1)
			vegleg.append('Mangrove')
		if any(VegXloc==2): # if marsh is present
			veg=Depth*0.+NaN;veg[VegXloc==2]=-Depth[VegXloc==2]
			plot(Xaxis,veg,'xr',linewidth=1)
			vegleg.append('Marsh')
		if any(VegXloc==3): # if Seagrass is present
			veg=Depth*0.+NaN;veg[VegXloc==3]=-Depth[VegXloc==3]
			plot(Xaxis,veg,'+g',linewidth=1)
			vegleg.append('Seagrass')
		if any(VegXloc==4): # if coral is present
			veg=Depth*0.+NaN;veg[VegXloc==4]=-Depth[VegXloc==4]
			plot(Xaxis,veg,'og',linewidth=1)
			vegleg.append('Coral Reef')
		if any(VegXloc==5): # if an oyster reef is present
			hi=-Depth[VegXloc==5]
			Yrf=num.arange((hi),0.0,0.05);Xrf=(Yrf*0.0+Xr) # x-axis
			plot(Xrf,Yrf,'g',Xrf+.05,Yrf,'g',Xrf+.1,Yrf,'g',linewidth=2);grid()
			vegleg.append('Oyster Reef')        
		box=ax.get_position();
		ax.set_position([box.x0, box.y0, box.width*1, box.height])                
		ylabel('Depth[m]',size='large')
		xlabel('Cross-Shore Distance from Shoreline',size='large')
		ax.legend((vegleg),prop=fontP,loc='upper left',ncol=2)
		if VegEnd == len (Wave):
			xlim(Xaxis[-1],Xaxis[0]);                
		else:
			xlim(Xaxis[VegEnd],Xaxis[0]);        
		axvline(x=0, ymin=0.75, linewidth=1, color='k') 
		savefig(outputws+"WavePlotZoomin_"+subwsStr+".png",dpi=(640/8))    
		
		# erosion plot
		if sand==1:
			# beach profile plot
			# create a profile that is a linear slope up to the berm elevation, flat across the berm, and a dune is placed at the end of the berm (if there is one).
			Xp=num.arange(0,500+W1,1)
			Yp=m*Xp
			loc=find(Yp>=B1)
			Yp[loc]=B1

			Lberm=loc[0];Ltoe=Indexed(Xp,Xp[Lberm+W1])
			Yp[Ltoe:-1]=B1+D1           
			# initialize message declaring whether dune or inland erosion occurs for output later
			preduneretreatmsg = "Null" 
			preinlanderrormsg = "Null"           
			postduneretreatmsg = "Null" 
			postinlanderrormsg = "Null"
			noaffect = "Null"
			mgmtbermgone = "Null"

			Yp0=m*Xp-1*m*Retreat1 # profile of erosion before management action
			loc=find(Yp0>=B1)
			if Retreat1 <= W1:
				Yp0[loc]=B1
				Lberm0=loc[0]
				W00=W1-(Xp[Lberm0]-Xp[Lberm])
				Ltoe0=Indexed(Xp,Xp[Lberm0+W00])
				Yp0[Ltoe0:-1]=B1+D1
			elif Retreat1 > W1 and D1<>0:
				preduneretreatmsg = "Your entire berm has been eroded and your dune will be retreated."
				Yp0[loc]=B1+D1
				Ltoe0=loc[0]
			elif Retreat1 > W1 and D1==0:
				preinlanderrormsg = "The erosion has removed your entire berm. Areas immediately inland are severely threatened."
				Yp0[loc]=B1
				Lberm0=loc[0] ## CRASHES HERE
				W00_inlanderosion=-1*(W1-(Xp[Lberm0]-Xp[Lberm]));

			if scenario == "B" or scenario == "D" or scenario == "E":
				Ypv=m*Xp-1*m*Retreat2 # profile of erosion after management action
				loc=find(Ypv>=B1)
				if Retreat1==Retreat2:
					noaffect = "Your dune height reduction does not change beach retreats under these forcing conditions."
				elif Retreat2 <= W1:
					Ypv[loc]=B1
					Lbermv=loc[0]
					W0v=W1-(Xp[Lbermv]-Xp[Lberm])
					Ltoev=Ltoe0
					Ypv[Ltoev:-1]=B1+D2
				elif Retreat2 > W1 and D2<>0:
					postduneretreatmsg = "Due to your management action, your entire berm has been eroded and your dune will be retreated. The entire dune may erode as well, threatening inland areas."
					Ypv[loc]=B1+D2
					Ltoev=loc[0]  
					W0v_inlanderosion=-1*(W1-(Xp[Ltoev]-Xp[Lberm]))
				elif Retreat2 > W1 and D2==0:
					postinlanderrormsg = "Due to your management action, erosion has removed your entire berm. Areas immediately inland are severely threatened"
					Ypv[loc]=B1
					Lbermv=loc[0]
					W0v_inlanderosion=-1*(W1-(Xp[Lbermv]-Xp[Lberm]))

			figure(3)
			ax=subplot(211)
			if scenario == "B" or scenario == "D"  or scenario == "E" and Retreat1 <> Retreat2:
				if Retreat2<=W1: # this means that the berm is wider than the erosion distance both pre- and post-management action (MA)
					axes_beg = max(Lberm-Retreat2, Lberm/2) # the start of the plot will begin offshore of the berm location by the same distance as post MA erosion
					axes_end = Lberm+W1 + 10 # the plot will only show the true berm length plus 10 meters
					plot(Xp,Yp,Xp,Yp0,'-.g',Xp,Ypv,'--r',linewidth=2);grid()
					xlim(Xp[axes_beg],Xp[axes_end]);ylim(0,B1+D1+.2)
					box=ax.get_position();
					ax.set_position([box.x0, box.y0, box.width*1,box.height])                
					ax.legend(('Initital','Before Mgmt','After Mgmt'),prop=fontP,loc='upper left',ncol=2)
					axvline(x=Lberm+W1, linewidth=2, color='k')
					annotate('Berm Limit',xy=(Lberm+W1,B1+.1))                
					ylabel('Backsh. Elev.[m]',size='large')
					xlabel('Cross-Shore Distance [m] from Shoreline',size='large')
					title('Erosion Profiles',size='large',weight='bold') 
				elif Retreat1 > W1: # this means that erosion for both pre- and post-MA has removed the entire berm
					if D1 <> 0 and D2 <> 0: # the dune will retreat both pre and post
						axes_beg = Lberm/2. # the start of the plot will be half way between the berm location and the shoreline
						axes_end = Lberm+Retreat2 + 10 # the profile plot will extend 10 meters past the limit of dune retreat
						plot(Xp,Yp,Xp,Yp0,'-.g',Xp,Ypv,'--r',linewidth=2);grid()
						xlim(Xp[axes_beg],Xp[axes_end]);ylim(0,B1+D1+.2)
						box=ax.get_position();
						ax.set_position([box.x0, box.y0, box.width*1,box.height])
						axvline(x=Lberm+W1, linewidth=2, color='k')
						annotate('Berm Limit',xy=(Lberm+W1,B1+.1))             
						ax.legend(('Initital','Before Mgmt','After Mgmt'),prop=fontP,loc='upper left',ncol=2)
						ylabel('Backsh. Elev.[m]',size='large')
						xlabel('Cross-Shore Distance [m] from Shoreline',size='large')
						title('Erosion Profiles',size='large',weight='bold')

					elif D1 == 0 or D2 == 0: # there is no dune pre- and/or post-management so erosion will extend past erodable portion and the inland area will be severly threatened.
						axes_beg = Lberm/2. # the start of the plot will be half way between the berm location and the shoreline
						axes_end = Lberm+Retreat2 + 10 # the profile plot will extend 10 meters past the erosion limit
						if D1 == 0:
							plot(Xp,Yp,Xp,Yp0,'-.g',Xp,Ypv,'--r',linewidth=2);grid()
							xlim(Xp[axes_beg],Xp[axes_end]);ylim(0,B1+D1+.2)
							fill_between(Xp[Lberm+W1:Lberm+W1+W00_inlanderosion],0,10,facecolor = 'g', alpha = 0.5)
							fill_between(Xp[Lberm+W1:Lberm+W1+ W0v_inlanderosion],0,10,facecolor = 'r', alpha = 0.5)
							axvline(x=Lberm+W1, linewidth=2, color='k')
							annotate('Berm Limit',xy=(Lberm+W1,B1+.1))                              
							ylabel('Backsh. Elev.[m]',size='large')
							xlabel('Cross-Shore Distance [m] from Shoreline',size='large')                            
							box=ax.get_position();
							ax.set_position([box.x0, box.y0, box.width*1,box.height])                
							ax.legend(('Initital','Before Mgmt','After Mgmt'),prop=fontP,loc='upper left',ncol=2)
						else:
							plot(Xp,Yp,Xp,Yp0,'-.g',Xp,Ypv,'--r',linewidth=2);grid()
							xlim(Xp[axes_beg],Xp[axes_end]);ylim(0,B1+D1+.2)
							fill_between(Xp[Lberm+W1:Lberm+W1+ W0v_inlanderosion],0,10,facecolor = 'r', alpha = 0.5)
							ylabel('Backsh. Elev.[m]',size='large')
							xlabel('Cross-Shore Distance [m] from Shoreline',size='large')
							box=ax.get_position();
							axvline(x=Lberm+W1, linewidth=2, color='k')
							annotate('Berm Limit',xy=(Lberm+W1,B1+.1))               
							ax.set_position([box.x0, box.y0, box.width*1,box.height])     
							ax.legend(('Initital','Before Mgmt','After Mgmt'),prop=fontP,loc='upper left',ncol=2) 
				else: 
					axes_beg = Lberm/2. # the start of the plot will be half way between the berm location and the shoreline
					axes_end = Lberm+Retreat2 + 10 # the profile plot will extend 10 meters past the limit of dune retreat                    
					plot(Xp,Yp,Xp,Yp0,'-.g',Xp,Ypv,'--r',linewidth=2);grid()
					xlim(Xp[axes_beg],Xp[axes_end]);ylim(0,B1+D1+.2)
					fill_between(Xp[Lberm+W1:Lberm+W1+ W0v_inlanderosion],0,10,facecolor = 'r', alpha = 0.5)
					ylabel('Backsh. Elev.[m]',size='large')
					xlabel('Cross-Shore Distance [m] from Shoreline',size='large')
					box=ax.get_position();
					axvline(x=Lberm+W1, linewidth=2, color='k')
					annotate('Berm Limit',xy=(Lberm+W1,B1+.1))                      
					ax.set_position([box.x0, box.y0, box.width*1,box.height])                
					ax.legend(('Initital','Before Mgmt','After Mgmt'),prop=fontP,loc='upper left',ncol=2)   
					mgmtbermgone = "Your management action has caused the your entire berm to be lost and inland areas are severely threatened under these conditions."

			else: # for all other scenarios or for scenarios 'B' and 'D' when then dune reduction has no impact on erosion, only one erosion plot is required
				if D1 <> 0 or Retreat1<=W1: # either the erosion does not completely remove the berm or, if it does, a berm is present to retreat
					axes_beg = Lberm/2. # start of the plot will be half way between the berm location and the shoreline
					axes_end = Lberm+max(Retreat1,W1)+ 10 # profile plot will extend 10 meters past the limit of berm erosion or dune retreat
					plot(Xp,Yp,Xp,Yp0,'-.g',linewidth=2);grid()
					xlim(Xp[axes_beg],Xp[axes_end]);ylim(0,B1+D1+.2)
					box=ax.get_position();
					ax.set_position([box.x0, box.y0, box.width*1,box.height])                
					ax.legend(('Initital','After Storm'),prop=fontP,loc='upper left',ncol=2)
					ylabel('Backsh. Elev.[m]',size='large')
					xlabel('Cross-Shore Distance [m] from Shoreline',size='large')
					title('Erosion Profile',size='large',weight='bold')
				else: # erosion has removed the berm and the lack of dune means that inland areas are severely threatened.
					axes_beg = Lberm/2. # start of the plot will be half way between the berm location and the shoreline
					axes_end = Lberm+max(Retreat1,W1)+ 10 # profile plot will extend 10 meters past the limit of berm erosion or dune retreat
					plot(Xp,Yp,Xp,Yp0,'-.g',linewidth=2);grid()
					xlim(Xp[axes_beg],Xp[axes_end]);ylim(0,B1+D1+.2)
					fill_between(Xp[Lberm+W1:Lberm+W1+W00_inlanderosion],0,10,facecolor = 'g', alpha = 0.5)
					box=ax.get_position();
					ax.set_position([box.x0, box.y0, box.width*1,box.height])                
					ax.legend(('Initital','After Storm'),prop=fontP,loc='upper left',ncol=2)
					ylabel('Backsh. Elev.[m]',size='large')
					xlabel('Cross-Shore Distance [m] from Shoreline',size='large')
					title('Erosion Profile',size='large',weight='bold')                    

			savefig(outputws+"ErosionBed_"+subwsStr+".png",dpi=(640/8))

		elif mud==1: # compute erosion amount for consolidated sediments
			figure(2)
			Xplot = Xaxis   
			Shorel=find(abs(Depth)==min(abs(Depth)) )
			if scenario == "E":
				ax=subplot(311)               
				plot(Xplot,Retreat1,'--g',Xplot,Retreat2,'r',linewidth=2);grid()
				box=ax.get_position();
				ax.set_position([box.x0, box.y0, box.width*1, box.height])     
				axvline(x=0, linewidth=1, color='k')
				xlim(0,Xplot[0])
				ylim(-0.,max(max(Retreat1[0:Shorel+10]),max(Retreat2[0:Shorel+10]))*1.1)
				ylabel('Erosion[cm/hr]',size='large')
				ax=subplot(312)
				plot(Xplot,Trms1,'--g',Xplot,Trms2,'r',Xplot,Te,'.k',linewidth=2);grid()
				box=ax.get_position();
				ax.set_position([box.x0, box.y0, box.width*1, box.height])                
				ax.legend(('Before Mgmt', 'After Mgmt', 'Mvt Thresh.'),prop=fontP,loc='upper left',ncol=2)
				xlim(0,Xplot[0]) 
				ylim(-0.,max([max(Trms1[0:Shorel+10]),max(Trms2[0:Shorel+10]),max(Te[0:Shorel+10])])*1.1)
				axvline(x=0, linewidth=1, color='k')
				ylabel('Stress[Nm-2]',size='large')
			else:
				ax=subplot(311) 
				plot(Xplot,Retreat1,'--g',linewidth=2);grid()
				box=ax.get_position();
				ax.set_position([box.x0, box.y0, box.width*1, box.height])     
				axvline(x=0, linewidth=1, color='k')
				xlim(10,Xplot[0]-10)
				ylabel('Erosion[cm/hr]',size='large')

				ax=subplot(312)
				plot(Xplot,Trms1,'--g',Xplot,Te,'.k',linewidth=2);grid()
				box=ax.get_position();
				ax.set_position([box.x0, box.y0, box.width*1, box.height])                
				ax.legend(('Bed Shear Stress', 'Mvt Thresh.'),prop=fontP,loc='upper left',ncol=2)
				xlim(10,Xplot[0]-10)  
				ylim(-0.,max(max(Trms1[0:Shorel+10]),max(Te[0:Shorel+10]))*1.1)
				axvline(x=0, linewidth=1, color='k')
				ylabel('Stress[Nm-2]',size='large')  

			ax=subplot(313)
			plot(Xplot,-Depth,'b',linewidth=2);grid()
			vegleg=['Bed']   

			if any(VegXloc==1): # if mangrove is present
				veg=Depth*0.+NaN;veg[VegXloc==1]=-Depth[VegXloc==1]
				plot(Xaxis,veg,'.r',linewidth=1)
				vegleg.append('Mangrove')
			if any(VegXloc==2): # if marsh is present
				veg=Depth*0.+NaN;veg[VegXloc==2]=-Depth[VegXloc==2]
				plot(Xaxis,veg,'xr',linewidth=1)
				vegleg.append('Marsh')
			if any(VegXloc==3): # if Seagrass is present
				veg=Depth*0.+NaN;veg[VegXloc==3]=-Depth[VegXloc==3]
				plot(Xaxis,veg,'+g',linewidth=1)
				vegleg.append('Seagrass')
			if any(VegXloc==4): # if coral is present
				veg=Depth*0.+NaN;veg[VegXloc==4]=-Depth[VegXloc==4]
				plot(Xaxis,veg,'og',linewidth=1)
				vegleg.append('Coral Reef')
			if any(VegXloc==5): # if an oyster reef is present
				hi=-Depth[VegXloc==5]
				Yrf=num.arange((hi),0.0,0.05);Xrf=(Yrf*0.0+Xr) # x-axis
				plot(Xrf,Yrf,'g',Xrf+.05,Yrf,'g',Xrf+.1,Yrf,'g',linewidth=2);grid()
				vegleg.append('Oyster Reef')               
			box=ax.get_position();
			ax.set_position([box.x0, box.y0, box.width*1, box.height])           
			ax.legend((vegleg),prop=fontP,loc='upper left',ncol=2)
			ylabel('Depth[m]',size='large')
			xlabel('Cross-Shore Distance [m]',size='large')
			xlim(10,Xaxis[0]-10)
			axvline(x=0, linewidth=1, color='k')
			ylim(-.25,-Depth[0]+1)
			savefig(outputws+"ErosionBed_"+subwsStr+".png",dpi=(640/8))
			Fig2=1 # for HTML

	except:
		gp.AddError(msgPlotErosion)
		raise Exception


	#################################
	########## VALUATION ############
	#################################

	if EconCP=='true':
		gp.AddMessage("\nComputing economic valuation...")

		# function to split thousands with commas
		def splitthousands(s, sep=','):
			if len(s) <= 3:
				return s
			return splitthousands(s[:-3], sep) + sep + s[-3:]

		if valua==1:
			if MErodeLen2==-9999 and Retreat2==-9999: # there's no management action
				gp.AddWarning("\nYou haven't defined a management action.  We cannot compute an avoided damage value.")
				valua=0;
			else:
				if sand==1:
					E1=max([Retreat1,Retreat2])
					E2=min([Retreat1,Retreat2])
				elif mud==1:
					E1=max([MErodeLen1,MErodeLen2])
					E2=min([MErodeLen1,MErodeLen2])

				if E2<0:
					gp.AddWarning("\nWe cannot compute an avoided damage cost because the retreat amount for one of your scenario is negative.  The biophysical model did not run appropriately.")
					valua=0
				else:
					E1=E1*Longshore;
					E2=E2*Longshore;
					Eav=E1-E2; # avoided erosion
					D1=E1*PropValue;
					D2=E2*PropValue;
					Dav=D1-D2; # avoided erosion
					p=1.0/Tr # return frequency
					temp1=1.0+disc;
					temp2=[(1.0/temp1)**t for t in range(1,int(TimeHoriz)+1)]
					EPV=p*Dav*sum(temp2)
					gp.AddMessage("...Avoided Erosion between scenarios is "+str(round(Eav))+" meters squared.\n...Avoided Damage Value is $"+splitthousands(str(int(Dav)))+" (in your local currency)\n...Expected Projected Value of habitat is $"+splitthousands(str(int(EPV)))+" (in your local currency)")


	#################################
	####### CREATE HTML FILE ########
	#################################

	try:
		# figures
		gp.AddMessage("\nWriting results to HTML output...\nLocation: "+outputws+"OutputWaveModel_"+subwsStr+".html\n")
		htmlfile.write("<HR><H2><u>Model Outputs</u></H2>")
		if WaveErosionQuestion=="(2) No, please compute these values from wind speed and fetch distance": 
			htmlfile.write("We computed offshore wave height and period values based on your input of wind speed (U=" +str(Us)+"m/s), fetch distance (Ft="+ str(round(Ft))+"km), and average water depth at your site(d=" +str(round(depth,1)) +"m).<br>")
		htmlfile.write('Offshore wave input conditions are:<br>Ho=' +str(round(Ho,2)) +'m, and To=' +str(round(To,2)) +'s<br>')

		if EconCP=='true':
			if valua==1:
				htmlfile.write("<p><i>Because of the presence of habitats</i>, the erosion difference between your scenarios is "+str(round(Eav))+" square meters.<p>We represent your local currency symbol with '$':<br><li> <b>Avoided damages from erosion</b> is $"+splitthousands(str(int(Dav)))+"<br><li> <b>Expected present value of the habitat lost or gained</b> in your scenario is $"+splitthousands(str(int(Dav)))+"<br>")
				#htmlfile.write("<p><center>Details of the valuation inputs and outputs:</center><p>")
				#htmlfile.write("<table border=\"1\" width=\"400\" cellpadding=\"2\" cellspacing=\"0\" align=\"center\"><tr align=\"center\">")
				#htmlfile.write("<td align=\"left\"><b>Inputs:</b></td>")
				#htmlfile.write("<th> </td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("<td align=\"center\" bgcolor=\"#E3170D\">Area of Erosion (E_1, sq. meters)<br>[no habitat]</td>")
				#htmlfile.write("<td align=\"center\" bgcolor=\"#E3170D\">"+str(round(E1))+"</td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("<td align=\"center\" bgcolor=\"#00FF00\">Area of Erosion (E_2, sq. meters)<br>[yes habitat]</td>")
				#htmlfile.write("<td align=\"center\" bgcolor=\"#00FF00\">"+str(round(E2))+"</td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("<td align=\"center\">Avoided Erosion (E_A, sq. meters)</td>")
				#htmlfile.write("<td align=\"center\">"+str(round(Eav))+"</td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("<td align=\"left\"><b>Outputs:</b></td>")
				#htmlfile.write("<th> </td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("<td align=\"center\" bgcolor=\"#E3170D\">Damage (D_1, $)<br>[no habitat]</td>")
				#htmlfile.write("<td align=\"center\" bgcolor=\"#E3170D\">$"+splitthousands(str(int(D1)))+"</td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("<td align=\"center\" bgcolor=\"#00FF00\">Damage (D_2, $)<br>[yes habitat]</td>")
				#htmlfile.write("<td align=\"center\" bgcolor=\"#00FF00\">$"+splitthousands(str(int(D2)))+"</td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("<td align=\"center\">Avoided Damage (D_A, $)</td>")
				#htmlfile.write("<td align=\"center\">$"+splitthousands(str(int(Dav)))+"</td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("<td align=\"center\">Expected Damage Reduction (EPV, $)</td>")
				#htmlfile.write("<td align=\"center\">$"+splitthousands(str(int(EPV)))+"</td>")
				#htmlfile.write("</tr>")
				#htmlfile.write("</table>")

			elif valua==0 and scenario == "C":
				htmlfile.write("We cannot compute an avoided damage cost because no management action has been taken on your habitat.<br>")
			elif valua == 0 and scenario <> "C":
				htmlfile.write("We cannot compute an avoided damage cost because the retreat amount for one of your scenario is negative.  The biophysical model did not run appropriately.<br>")
			elif valua==-1:
				htmlfile.write("One or more of the required economic valuation input parameters was not defined.  We cannot compute a value for your habitats.<br>")

		htmlfile.write("<p>")

		if shal:
			htmlfile.write("Offshore water depth was too shallow for input wave height. We will assume that it broke somewhere in deeper water.")
		if scenario == "E":
			htmlfile.write("The loss of habitat caused the total water level to increase by "+str(round((RnpvegMA-Rnpveg1)/RnpvegMA*100))+"%.<br>")
			if sand == 1:
				if inundation_pre <> 0:
					htmlfile.write("The storm surge plus wave runup is higher than your berm and dune system.  Inland areas will experience flooding.")
				elif HshortwaveMA+S+D2<>Hshortwave+S+D1 and inundation_post <> 0:
					htmlfile.write("Your management action has led to inundation of your berm and dune system.  Inland areas will experience flooding due to your management action.")
				if noaffect <> "Null":
					htmlfile.write(noaffect+"<br>")
				if mgmtbermgone <> "Null":
					htmlfile.write(mgmtbermgone+"<br>")
				if preduneretreatmsg == "Null" and postduneretreatmsg == "Null" and preinlanderrormsg == "Null" and postinlanderrormsg == "Null":
					if D1 <> 0:
						htmlfile.write("Your berm is so wide that your dune is not affected. <br>")
					else:
						htmlfile.write("Your berm is wide enough such that some of the berm remains after erosion under these forcing conditions.<br>")
				elif preduneretreatmsg <> "Null":
					htmlfile.write(preduneretreatmsg+"<br>")
				elif postduneretreatmsg <> "Null":
					htmlfile.write(postduneretreatmsg+"<br>")
				elif preinlanderrormsg <> "Null":
					htmlfile.write(preinlanderrormsg+"<br>")
				elif postinlanderrormsg <> "Null":
					htmlfile.write(postinlanderrormsg+"<br>")  

		htmlfile.write("<table border=\"0\" width=\"1280\" cellpadding=\"5\" cellspacing=\"10\"><tr><td>")
		htmlfile.write("The figure below shows profiles of wave height along your profile,<br>before and after your management action.<br>")
		htmlfile.write("<img src=\"WavePlotZoomin_"+subwsStr+".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\">")
		if scenario == "E":
			htmlfile.write("</td><td>Because of your management action, wave height increased on average by " +str(round(WAtn,2))+"%  and wave energy increased by " +str(round(EAtn,2))+"% near areas where you had natural habitats.")
		if Oyster:
			htmlfile.write("<br><b>The oyster reef attenuated " +str(round((1-Kt)*100,1)) +"% of the wave height and " +str(round((1-Kt*Kt)*100,1)) +"% of the wave energy right offshore of the reef..</b></td></tr>")
		else:
			htmlfile.write("</td></tr>")
		htmlfile.write("<tr><td>The figure below shows profiles of erosion in the backshore region,<br>before and after your management action.<br>")
		htmlfile.write("<img src=\"ErosionBed_"+subwsStr+".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\"></td>")

		# erosion
		if sand==1:
			htmlfile.write("<td>")
			if Retreat2<>-9999:
				htmlfile.write("<u>Before</u> your management action, your beach might erode by " +str(round(Retreat1, 1)) +"m.<br><u>After</u> your management action, your beach might erode by " +str(round(Retreat2, 1)) +"m.")
				if S+Rnp1 > B1 and D1 <> 0:
					htmlfile.write("<br>Under these forcing conditions, your berm is submerged by surge and runup causing your dune to be exposed to wave attack.  This will lead to scarping or more severe dune erosion than shown.")
			else:
				htmlfile.write("Under these forcing conditions, your beach might erode by " +str(round(Retreat1, 1)) +"m.")
				if S+Rnp1 > B1 and D1 <> 0:
					htmlfile.write("<br>Under these forcing conditions, your berm is submerged by surge and runup causing your dune to be exposed to wave attack.  This will lead to scarping or more severe dune erosion than shown.")  
			if preinlanderrormsg <> "Null":
				htmlfile.write("<br><u>Under this forcing condtion</u> your berm has been completely eroded.  The green shading in the erosion plot indicates the area inland of the berm which will be severly threatened under these conditions.")
			elif postinlanderrormsg <> "Null":
				htmlfile.write("<br><u>Due to your management action</u> your berm has been completely eroded and there is no dune to protect inland areas.  Inland areas are severely threatened.")
			htmlfile.write("</td></tr>")
		elif mud==1:
			htmlfile.write("<td>")
			if scenario=="E":
				htmlfile.write("<u>Before</u> your management action, your muddy backshore area might experience a volumetric loss of " +str(round(MErodeVol1,2)) +" cubic meters.  The length of the mud bed where sediment mobilization occurs is " +str(round(MErodeLen1,2)) +" meters.<br><u>After</u> your management action, your muddy backshore area might experience a volumetric loss of " +str(round(MErodeVol2,2)) +" cubic meters.  The length of the mud bed where sediment mobilization occurs is " +str(round(MErodeLen2,2)) +" meters.</td></tr>")
			else:
				htmlfile.write("Under these forcing conditions, your muddy backshore area might experience a volumetric loss of " +str(round(MErodeVol1,2)) +" cubic meters.  The length of the mud bed where sediment mobilization occurs is " +str(round(MErodeLen1,2)) +" meters.</td><tr>")

		htmlfile.write("<tr><td>The figure below shows profiles of wave height along your whole profile.  It's a zoom out of the first figure.<br>")
		htmlfile.write("<img src=\"WavePlot_"+subwsStr+".png\" alt=\"Profile Plot #1\" width=\"640\" height=\"480\"></td></tr>")
		htmlfile.write("</table></html>")
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