from ROOT import *
from math import *
from array import *
import numpy as np
from scipy.optimize import fmin,minimize
import os
import sys
global Filedat,Ebinbefore,Ebin
# my condition
number_simulation=5
initialguesspar=[16000,2.849,2.716,336,0.000215] #initial guess Norm,gamma1,gamm2,Ebreak
model='SPLwHe.f' #choose model SPL,BPL,SPLwHe,BPLwHe
def Fluxcompute(A,gamma1,gamm2,Ebreak,normAll):
	RunFlux='./test1.out %f %f %f %f %f'%(A,gamma1,gamm2,Ebreak,normAll)
	os.system(RunFlux)
def SumlogPois(dummy):
	#print dummy
	A=dummy[0]
	gamma1=dummy[1]
	gamma2=dummy[2]
	Ebreak=dummy[3]
	normAll=dummy[4]
	Fluxcompute(A,gamma1,gamma2,Ebreak,normAll)
	file=open('0.dat')
	data=np.genfromtxt('0.dat')
	file.close()
	x,y=data[:,0],data[:,1]
	sumlogpois=0
	for i in range(len(x)):
		measurement=Sim_Flux275.Eval(x[i],0,'S') # Cubic spline interpolate
		model=y[i]*(x[i]**2.75)
		if TMath.Poisson(measurement,model)==0:
			sumlogpois+=308.
		if TMath.Poisson(measurement,model)!=0:
			sumlogpois+=-log(TMath.Poisson(measurement,model))
	return sumlogpois
def SimulateFlux(flux275): # include
	flux275=[]
	dNsb,Eavgbin,flxlimb275=Filedat[:,0],Filedat[:,1],Filedat[:,2]
	# simulate random coutn (Statistical error)
	simstat10GeV_flux275=(flxlimb275[0]/dNsb[0])*gRandom.PoissonD(dNsb[0])
	simstat100GeV_flux275=(flxlimb275[24]/dNsb[24])*gRandom.PoissonD(dNsb[24])
	simstat1000GeV_flux275=(flxlimb275[49]/dNsb[49])*gRandom.PoissonD(dNsb[49])
	# simulate Systematic error (Aeff err.)
	simtot10GeV_flux275=gRandom.Gaus(simstat10GeV_flux275,simstat10GeV_flux275*0.05)
	simtot100GeV_flux275=gRandom.Gaus(simstat100GeV_flux275,simstat100GeV_flux275*0.05)
	simtot1000GeV_flux275=gRandom.Gaus(simstat1000GeV_flux275,simstat1000GeV_flux275*0.15)
	flux275.append(simtot10GeV_flux275)
	flux275.append(simtot100GeV_flux275)
	flux275.append(simtot1000GeV_flux275)
	return flux275
if __name__ == "__main__":
	# Initialize model
	os.system('gfortran %s frag.f -o test1.out' %(model))
	# open dat file
	Filedat=np.genfromtxt('alldat.olo')
	Eavgbin=Filedat[:,1] # GOT Emidbin
	# choose Emidbin only 3 point
	Eavgbin_simulate=[Eavgbin[0],Eavgbin[24],Eavgbin[49]]
    # open to write output parameters
	foutput=open('outputTotal.dat','w')
	for i in range(number_simulation):
		Flux275=[] # create global variable
		Flux275=SimulateFlux(Flux275) # simulate new flux (Random Error stat.)
		Sim_Flux275=TGraph(3,array('d',Eavgbin_simulate),array('d',Flux275))
		bestfit=fmin(SumlogPois,initialguesspar)
		print bestfit
		foutput.write('%f %f %f \n' %(bestfit[1],bestfit[2],bestfit[3]))
# close dat file
foutput.close()
