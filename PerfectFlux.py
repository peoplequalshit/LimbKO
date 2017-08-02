from math import *
from ROOT import *
from array import *
import numpy as np
# our condition (Limb Peak at nadir 68.02)
Zmin=110.0 # nadir 70.0
Zmax=111.6 # nadir 68.4
Zbgmin=100. # nadir 80
Zbgmax=108. # nadir 72
solidangle=(cos(Zmin*(pi/180.))-cos(Zmax*(pi/180.)))*2.*pi
solidanglebg=(cos(Zbgmin*(pi/180.))-cos(Zbgmax*(pi/180.)))*2.*pi
# Energy bin
oV=[(10**((float(i)/25)+1)) for i in range(51)]
V=array('d',oV)
# Build count map
cntmap=[]
# open expmap
Fexpmap=TFile('ExpMap_P8R2_ULTRACLEANVETO_V6_w010-w399.root')
name_expmap=[]# declare map name in expmap
for i in range(50):
    namecntmap='cntmap%03d'%(i)
    cntmap.append(TH2F(namecntmap,namecntmap,180,0.,360.,800,0.,80.))
    name_expmap.append('expmap%03d'%(i))
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
for i in range(50): #have 0-50 but interest just 1-50
    dN.append(0)
    dNbg.append(0)
    EdN.append(0)
    EavgdN.append(0)
    EavgdNbg.append(0)
    dNsb.append(0)
# process data
for event in ev:
#    print event.EVENTS
#    energy=ev.ENERGY
#    if np.searchsorted(V,energy)>0 and np.searchsorted(V,energy)<51:
#       if event.ZENITH>Zmin and event.ZENITH<Zmax and event.THETA<70.:
#           dN[np.searchsorted(V,energy)]+=1.
#           E.Fill(energy)
#print dN
    #print event.EVENTS
    energy=0.963*ev.ENERGY # bias energy
    if np.searchsorted(V,energy)>0 and np.searchsorted(V,energy)<51:
        cntmap[np.searchsorted(V,energy)-1].Fill(event.PHI,180.-event.ZENITHSHIFT)
        if event.ZENITHSHIFT>Zmin and event.ZENITHSHIFT<Zmax and event.THETA < 70.:# limb count
            dN[np.searchsorted(V,energy)-1]+=1.
            EavgdN[np.searchsorted(V,energy)-1]+=energy # sum before, average in next for-loop
        if event.ZENITHSHIFT>Zbgmin and event.ZENITHSHIFT<Zbgmax and event.THETA < 70.:# bg count
            dNbg[np.searchsorted(V,energy)-1]+=1.
            EavgdNbg[np.searchsorted(V,energy)-1]+=energy
#print dN
# create strMap
strmap=TH2F('strmap','strmap',180,0.,360.,800,0.,80.)
dphi=(360./strmap.GetNbinsX())*pi/180.
dtheta=(80./strmap.GetNbinsY())*pi/180.
for j in range(strmap.GetNbinsY()):
    strbin=dphi*dtheta*sin((j+0.5)*dtheta) # d(phi)d(theta)sin(theta)
    for i in range(strmap.GetNbinsX()):
        strmap.SetBinContent(i+1,j+1,strbin)
# define flxmap
flxmap=[]
for i in range(50):
    flxmapadd=cntmap[i].Clone()
    flxmap.append(flxmapadd)
flxvallimb=[]
flxvalbg=[]
# write dNsb,EavgdN to file
f1=file('alldat.olo','w')
for i in range(len(V)-1):
    # dE
    dE=V[i+1]-V[i]
    # Flxmap : flxmap=(cntmap/expmap)/dE/dOmega
    expmap=Fexpmap.Get(name_expmap[i])
    #expmap.Scale(1./10000.) # cm^2->m^2
    #
    flxmap[i].Divide(cntmap[i],expmap)
    flxmap[i].Divide(flxmap[i],strmap)
    # get flux value limb
    flxmap[i].GetXaxis().SetRangeUser(0.,360.)
    flxmap[i].GetYaxis().SetRangeUser(180.-Zmax,180.-Zmin)
    flxvallimb.append(flxmap[i].Integral()/(dE))
    # get flux value bg
    flxmap[i].GetXaxis().SetRangeUser(0.,360.)
    flxmap[i].GetYaxis().SetRangeUser(180.-Zbgmax,180.-Zbgmin)
    flxvalbg.append(flxmap[i].Integral()/dE)
    #
    dNsb[i]=dN[i]-dNbg[i]*((Zmin-Zmax)/(Zbgmin-Zbgmax)) # weight str bg ti str limb
    EavgdN[i]=EavgdN[i]/dN[i]
    EavgdNbg[i]=EavgdNbg[i]/dNbg[i]
    f1.write('%f %f %e %f %f %e\n'%(dNsb[i],EavgdN[i],flxvallimb[i]*(EavgdN[i]**2.75),dNbg[i],EavgdNbg[i],flxvalbg[i]*(EavgdNbg[i]**2.75)))
    # write count map in root file
# Close all file
Fexpmap.Close()
f1.close()
