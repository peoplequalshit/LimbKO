from ROOT import *
from math import *
from array import *
import numpy as np
oV=[(10**((float(i)/25)+1)) for i in range(51)]
V=array('d',oV)
E1=TH1F('E1','E1',len(V)-1,V)
Ebg=TH1F('Ebg','Ebg',len(V)-1,V)
Zmin=110.0
Zmax=111.6
Zbgmin=100.
Zbgmax=108.
solidangle=(cos(Zmin*(pi/180.))-cos(Zmax*(pi/180.)))*2*pi
solidanglebg=(cos(Zbgmin*(pi/180.))-cos(Zbgmax*(pi/180.)))*2*pi
dat=np.genfromtxt('alldat.olo')
dN,Eavgbin,expmap=dat[:,0],dat[:,1],dat[:,2]/10000.
dNbg,Eavgbinbg,expmapbg=dat[:,3],dat[:,4],dat[:,5]
Flux=[]
Fluxbg=[]
for i in range(50):
    Flux.append(dN[i]/(E1.GetBinWidth(i+1)*solidangle*expmap[i]))
    Fluxbg.append(dNbg[i]/(Ebg.GetBinWidth(i+1)*solidanglebg*expmapbg[i]))
    print Flux[i]
    E1.SetBinContent(i+1,Flux[i])
    Ebg.SetBinContent(i+1,Fluxbg[i])
C=TCanvas('C','C',800,600)
E1.SetMarkerStyle(20)
E1.Draw('P')
Ebg.SetMarkerStyle(32)
Ebg.Draw('Psame')
C.SetLogx()
C.SetLogy()
raw_input()
