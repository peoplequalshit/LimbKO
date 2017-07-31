from math import *
from ROOT import *
from array import *
import numpy as np
# our condition (Limb Peak at nadir 68.02)
Zmin=110.0 # nadir 70.0
Zmax=111.6 # nadir 68.4
Zbgmin=100. # nadir 80
Zbgmax=108. # nadir 72
# Energy bin
oV=[(10**((float(i)/25)+1)) for i in range(51)]
V=array('d',oV)
# Build count map
Fcntmap=TFile('CntMap.root','RECREATE')
cntmap=[]
# open expmap
Fexpmap=TFile('ExpMap_P8R2_ULTRACLEANVETO_V6_w010-w399.root')
name_expmap=[]# declare map name in expma
for i in range(50):
    namecntmap='cntmap%.3d'%(i)
    cntmap.append(TH2F(namecntmap,namecntmap,72,0.,360.,800,0.,80.))
    if i<10:
            name_expmap.append('expmap00%i'%i)
    if i>=10:
            name_expmap.append('expmap0%i'%i)
# import data from tree
ev=TChain('Data of photon and spacecraft')
ev.Add('finaltree.root')
# Declare variable
dN=[] # raw count
dNbg=[] # background count
EdN=[]
EavgdN=[] # Energy average in each limb bin
EavgdNbg=[] # Energy average in each bg bin
dNsb=[] # dN subtract background
for i in range(51): #have 0-50 but interest just 1-50
    dN.append(0)
    dNbg.append(0)
    EdN.append(0)
    EavgdN.append(0)
    EavgdNbg.append(0)
    dNsb.append(0)
# process data
E=TH1F('E','E',len(V)-1,V)
for event in ev:
#    print event.EVENTS
#    energy=ev.ENERGY
#    if np.searchsorted(V,energy)>0 and np.searchsorted(V,energy)<51:
#       if event.ZENITH>Zmin and event.ZENITH<Zmax and event.THETA<70.:
#           dN[np.searchsorted(V,energy)]+=1.
#           E.Fill(energy)
#print dN
#exit()
    print event.EVENTS
    energy=0.963*ev.ENERGY # bias energy
    if np.searchsorted(V,energy)>0 and np.searchsorted(V,energy)<51:
        if event.ZENITHSHIFT>Zmin and event.ZENITHSHIFT<Zmax and event.THETA < 70.:# limb count
            dN[np.searchsorted(V,energy)]+=1.
            EavgdN[np.searchsorted(V,energy)]+=energy # sum before, average in next for-loop
        if event.ZENITHSHIFT>Zbgmin and event.ZENITHSHIFT<Zbgmax and event.THETA < 70.:# bg count
            dNbg[np.searchsorted(V,energy)]+=1.
            EavgdNbg[np.searchsorted(V,energy)]+=energy
# write dNsb,EavgdN to file
f1=file('alldat.olo','w')

for i in range(len(V)-1):
    # get expmap in each bin
    expmap=Fexpmap.Get(name_expmap[i])
    # int expmap limb
    expmap.GetYaxis().SetRangeUser(180.-Zmax,180.-Zmin)
    expmap.GetXaxis().SetRangeUser(0.,360.)
    intexplimb=expmap.Integral()/10000. #cm^2->m^2
    # int expmap bg
    expmap=Fexpmap.Get(name_expmap[i])
    expmap.GetYaxis().SetRangeUser(180.-Zbgmax,180.-Zbgmin)
    expmap.GetXaxis().SetRangeUser(0.,360.)
    intexpbg=expmap.Integral()/10000. #cm^2->m^2
    i=i+1
    dNsb[i]=dN[i]-dNbg[i]*((Zmin-Zmax)/(Zbgmin-Zbgmax)) # weight str bg ti str limb
    EavgdN[i]=EavgdN[i]/dN[i]
    EavgdNbg[i]=EavgdNbg[i]/dNbg[i]
    f1.write('%f %f %f %f %f %f\n'%(dNsb[i],EavgdN[i],intexplimb,dNbg[i],EavgdNbg[i],intexpbg))
    # write count map in root file
f1.close()
