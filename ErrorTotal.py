import pyLikelihood
import pyIrfLoader
from ROOT import *
from math import *
from array import *
import numpy as np
from scipy.optimize import fmin,minimize
import os
import sys
import matplotlib.pyplot as plt ##try

# my condition
number_simulation=1
initialguesspar=[34554,2.9,2.7,330.,0.000215] #initial guess Norm,gamma1,gamm2,Ebreak
model='BPLwHe.f' #choose model SPL,BPL,SPLwHe,BPLwHe
livetime=70761348.6153 # waiting for Warit's exposure map
Zmin=110.0
Zmax=111.6
solidangle=(cos(Zmin*(pi/180.))-cos(Zmax*(pi/180.)))*2.*pi

# setting LAT performance
pyIrfLoader.Loader_go()
myFactory=pyIrfLoader.IrfsFactory_instance()
irfs_f=myFactory.create("P8R2_ULTRACLEANVETO_V6::FRONT")
irfs_b=myFactory.create("P8R2_ULTRACLEANVETO_V6::BACK")
aeff_f=irfs_f.aeff()
aeff_b=irfs_b.aeff()

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
def SimulateFlux(flux):
	flux=[]
	for i in range(50):
		flux.append(0)
        Fileflux=np.genfromtxt('dNsb.dat')
        dNsb,Eavgbin=Fileflux[:,0],Fileflux[:,1]
	for i in range(len(flux)):
		dNsb[i]=gRandom.PoissonD(dNsb[i]) # Random new dNsb
		exposure=(aeff_f.value(Eavgbin[i]*1000.,61.5,0.)+aeff_b.value(Eavgbin[i]*1000.,61.5,0.))/10000.
                binwidth=Ebin[i+1]-Ebin[i]
		flux[i]=(((dNsb[i]/binwidth)/exposure)/livetime)/solidangle
	return flux
if __name__ == "__main__":
	# declare energy bin
    Ebinbefore=[(10**((float(i)/25)+1)) for i in range(51)]
    Ebin=array('d',Ebinbefore)
	# d
	Fileflux=np.genfromtxt('dNsb.dat')
	Eavgbin=Fileflux[:,1] # GOT Emidbin
	# choose Emidbin only 3 point
	Eavgbin_simulate=[Eavgbin[0],Eavgbin[24],Eavgbin[49]]
        # open to write output parameters
        foutput=open('outputTotal.dat','w')
	for i in range(number_simulation):
		Flux=[] # create global variable
		Flux=SimulateFlux(Flux) # simulate new flux (Random Error stat.)
		# random new flux only 3 point (Error sys.)
		Flux_simulate1=gRandom.Gaus(Flux[0],Flux[0]*0.05)
		Flux_simulate2=gRandom.Gaus(Flux[24],Flux[24]*0.05)
		Flux_simulate3=gRandom.Gaus(Flux[49],Flux[49]*0.15)
		Flux_simulate=[Flux_simulate1,Flux_simulate2,Flux_simulate3]
                Flux_simulate275=[]
                for i in range(3):
                        Flux_simulate275.append(Flux_simulate[i]*(Eavgbin_simulate[i]**2.75))
		FinalFlux_sim275=TGraph(3,array('d',Eavgbin_simulate),array('d',Flux_simulate275))#cubic spline
                bestfit=fmin(SumlogPois,[32000,2.8,2.6,300.0,0.000215])
		foutput.write('%f %f %f \n' %(bestfit[1],bestfit[2],bestfit[3]))

# close dat file
foutput.close()
