from ROOT import *
from math import *
from array import *
import numpy as np
from scipy.optimize import fmin,minimize
import os
import sys
global Filedat,Ebinbefore,Ebin
# my condition
number_simulation=1
initialguesspar=[34554,2.9,2.7,330.,0.000215] #initial guess Norm,gamma1,gamm2,Ebreak
model='BPLwHe.f' #choose model SPL,BPL,SPLwHe,BPLwHe
livetime=70761348.6153 # waiting for Warit's exposure map
Zmin=110.0
Zmax=111.6
solidangle=(cos(Zmin*(pi/180.))-cos(Zmax*(pi/180.)))*2.*pi

def Fluxcompute(A,gamma1,gamm2,Ebreak,normAll):
	os.system('gfortran %s frag.f -o test1.out' %(model))
	RunFlux='./test1.out %f %f %f %f %f'%(A,gamma1,gamm2,Ebreak,normAll)
	os.system(RunFlux)
def SumlogPois(dummy):
        print dummy
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
		measurement=FinalFlux_sim275.Eval(x[i],0,'S') # Cubic spline interpolate
		model=y[i]*(x[i]**2.75)
		if TMath.Poisson(measurement,model)==0:
			sumlogpois+=308
		if TMath.Poisson(measurement,model)!=0:
			sumlogpois+=-log(TMath.Poisson(measurement,model))
	return sumlogpois
def SimulateFlux(flux): # Simlate Random count (Stat. err.)
	flux=[]
    dNsb,Eavgbin,expmap=Filedat[:,0],Filedat[:,1],Filedat[:,2]
	for i in range(len(flux)):
		dNsb[i]=gRandom.PoissonD(dNsb[i]) # Random new dNsb
        binwidth=Ebin[i+1]-Ebin[i]
		flux.append(dNsb[i]/(binwidth*solidangle*expmap[i]))
	return flux
if __name__ == "__main__":
	# Declare energy bin
	Ebinbefore=[(10**((float(i)/25)+1)) for i in range(51)]
    Ebin=array('d',Ebinbefore)
	# Open dat file
	Filedat=np.genfromtxt('alldat.dat')
	Eavgbin=Filedat[:,1] # GOT Emidbin
    # open to write output parameters
	foutput=open('outputStat.dat','w')
	for i in range(number_simulation):
		Flux=[] # create variable
		Flux=SimulateFlux(Flux) # simulate new flux (Random Error stat.)
		# let Flux to E^{2.75}Flux
		Flux_simulate275=[]
		for i in range(len(Flux)):
			Flux_simulate275.append(Flux[i]*(Eavgbin[i]**2.75))
		FinalFlux_sim275=TGraph(50,array('d',Eavgbin),array('d',Flux_simulate275))
		bestfit=fmin(SumlogPois,[32000,2.8,2.6,300.0,0.000215])
		foutput.write('%f %f %f \n'%(bestfit[1],bestfit[2],bestfit[3]))

# close dat file
foutput.close()
