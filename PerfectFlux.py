from math import *
from ROOT import *
from array import *
import numpy as np
# our condition
Zmin=110.0
Zmax=111.6
Zbgmin=100. # nadir 80
Zbgmax=108. # nadir 72
# Energy bin
oV=[(10**((float(i)/25)+1)) for i in range(51)]
V=array('d',oV)
# open expmap
Fexpmap=TFile('ExpMap_P8R2_ULTRACLEANVETO_V6_w010-w399.root')
name_expmap=[]# declare map name in expmap
for i in range(50):
    if i<10:
            name_expmap.append('expmap00%i'%i)
    if i>=10:
            name_expmap.append('expmap0%i'%i)
# import data from tree
ev=TChain('Data of photon and spacecraft')
ev.Add('tree.root')
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
for event in ev:
    print event.EVENTS
    energy=0.936*ev.ENERGY # bias energy shift
    if np.searchsorted(V,energy)>0 and np.searchsorted(V,energy)<51:
        if event.ZENITHSHIFT>Zmin and event.ZENITHSHIFT<Zmax and event.THETA < 70:# raw count
            dN[np.searchsorted(V,energy)]+=1.
            EavgdN[np.searchsorted(V,energy)]+=energy # sum before, average in next for-loop
        if event.ZENITHSHIFT>Zbgmin and event.ZENITHSHIFT<Zbgmax and event.THETA < 70:# bg count
            dNbg[np.searchsorted(V,energy)]+=1.
            EavgdNbg[np.searchsorted(V,energy)]+=energy
# write dNsb,EavgdN to file
f1=file('alldat.olo','w')
for i in range(len(V)-1):
    # get expmap in each bin
    expmap=Fexpmap.Get(name_expmap[i])
    # int expmap limb
    x1=expmap.GetXaxis().FindBin(0.)
    x2=expmap.GetXaxis().FindBin(360.)
    y1=expmap.GetYaxis().FindBin(180.-Zmax)
    y2=expmap.GetYaxis().FindBin(180.-Zmin)
    intexplimb=expmap.Integral(x1,x2,y1,y2)
    # int expmap bg
    x1=expmap.GetXaxis().FindBin(0.)
    x2=expmap.GetXaxis().FindBin(360.)
    y1=expmap.GetYaxis().FindBin(180.-Zbgmax)
    y2=expmap.GetYaxis().FindBin(180.-Zbgmin)
    intexpbg=expmap.Integral(x1,x2,y1,y2)
    i=i+1
    dNsb[i]=dN[i]-dNbg[i]*((Zmin-Zmax)/(Zbgmin-Zbgmax)) # weight str bg ti str limb
    EavgdN[i]=EavgdN[i]/dN[i]
    EavgdNbg[i]=EavgdNbg[i]/dNbg[i]
    print EavgdN[i],EavgdNbg[i] # fortest
    f1.write('%f %f %f %f %f %f\n'%(dNsb[i],EavgdN[i],intexplimb,dNbg[i],EavgdNbg[i],intexpbg))
f1.close()
