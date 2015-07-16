# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:07:35 2015

@author: jmwilson
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import h5py

#Schorlemmer 07 RELM testing area
polylats = [430, 430, 394, 357, 343, 329, 322, 317, 315, 319, 328, 337, 342, 377, 402, 405]
polylons = [-1252, -1190, -1190, -1140, -1131, -1135, -1136, -1145, -1171, -1179, -1184, -1210, -1216, -1238, -1254, -1254]
calipoly = zip(polylons, polylats)

def isinPoly(pt, polyvecs=None):
    # be careful; this needs to be a polyvec, not just a polygon...
    typ = type
    minim = min
    maxim = max
    flt = float
    if polyvecs==None: polyvecs=getPolyVecs(calipoly)
    if typ(polyvecs[0][0])==typ(0) or typ(polyvecs[0][0]) == typ(1.0):
        polyvecs=getPolyVecs(polyvecs)
    #
    # now, how many "left of" (or "right of") vectors (where y1 < y_pt < y2)
    x0=pt[0]
    y0=pt[1]
    nlefts=0
    isin=False
    #
    for vec in polyvecs:
        x1, x2 = vec[0][0], vec[1][0]
        y1, y2 = vec[0][1], vec[1][1]
        #if (y0>y1 and y0<y2) or (y0<y1 and y0>y2):
        if (y0<minim(y1, y2)) or (y0>=maxim(y1, y2)): continue
        #
        #we're in the y-domain of this vector.
        if x0<x1 and x0<x2:
            nlefts+=1
        elif x0>=x1 and x0>=x2:
            continue
        else:
            #might still be left...
            slope=(y2-y1)/flt(x2-x1)
            y_line = y1 + slope*(x0-x1)
            if slope>0 and y0>=y_line:
                nlefts+=1
            elif slope<0 and y0<y_line:
                nlefts+=1
            #
        #
    #
    if nlefts%2==0: isin=False
    if nlefts%2==1: isin=True
    #
    return isin


def getPolyVecs(poly=None):
    if poly==None:
        poly = calipoly
    polyvecs=[]
    for i in xrange(len(poly)):
        polyvecs+=[[poly[i], poly[(i+1)%(len(poly))]]]
        if polyvecs[-1][0]==polyvecs[-1][1]: polyvecs.pop()
    return polyvecs


def makePolyMask(binpositions, polyvecs=None, maskshape=()):
    if polyvecs==None:
        polyvecs = getPolyVecs()
    #
    binpositions = np.array(binpositions)
    masklist = np.zeros((np.shape(binpositions)[:-1])  + maskshape)
    #
    for j in xrange(len(binpositions)):
        for i in xrange(len(binpositions[0])):
            inarea = isinPoly(binpositions[j][i], polyvecs=polyvecs)
            if inarea == True:
                pass
            else:
                masklist[j][i] = np.ones(maskshape).astype(int)
    #
    return np.array(masklist)


def getFaultCoords(faultfile="../ALLCAL2_Faults/nocreep/ALLCAL2_VQ_faults_3km.h5"):
    faults = h5py.File(faultfile,'r')
    vert_coords = []
    for vert in faults['vertices']:
        coord = (int(round(vert[2]*10)),int(round(vert[1]*10)))
        if coord not in vert_coords:
            vert_coords.append(coord)
    vert_coords = np.array(vert_coords)
    return vert_coords

def makeFaultMask(binpositions, vertcoords=None, maskshape=()):
    #
    if vertcoords == None:
        vertcoords = getFaultCoords()
    minlon = np.min(binpositions[...,0])
    minlat = np.min(binpositions[...,1])
    masklist = np.ones((np.shape(binpositions)[:-1])  + maskshape)
    #
    for vert in vertcoords:
        masklist[vert[1]-minlat][vert[0]-minlon] = np.zeros(maskshape).astype(int)

    return masklist
    
def sphericalDistNP(singleloc, otherlocs, distType='rad'):
    Rearth = 6371.0    # km
    deg2rad = math.pi/180.0
    #
    # phi = longitude, lambda = latitude
    phis = singleloc[0]*deg2rad    
    lams  = singleloc[1]*deg2rad
    phif = otherlocs[..., 0]*deg2rad
    lamf  = otherlocs[..., 1]*deg2rad
    dlam = (lamf - lams)
    dphi = (phif - phis)
    #
    #
    if distType=='rad': dsig = vincentyNP(phis, lams, phif, lamf, dphi)
    if distType=='rect':
        lamav = (lamf+lams)/2.0
        xsigs = vincentyNP(phis, lamav, phif, lamav*np.ones(phif.shape), dphi)
        xsigs[dphi<0]*=-1
        ysigs = vincentyNP(phis, lams, phis*np.ones(lamf.shape), lamf, 0.0)
        ysigs[dlam<0]*=-1
        dsig = np.concatenate((xsigs[...,np.newaxis], ysigs[...,np.newaxis]), axis=xsigs.ndim)
    #
    Dist = Rearth * dsig
    #
    return Dist

def vincentyNP(phis, lams, phif, lamf, dphi):
    #this one is supposed to be bulletproof (Vincenty formula):
    dsigma = np.arctan2( np.sqrt((np.cos(lamf)*np.sin(dphi))**2.0 + (np.cos(lams)*np.sin(lamf) - np.sin(lams)*np.cos(lamf)*np.cos(dphi))**2.0 ) , (np.sin(lams)*np.sin(lamf) + np.cos(lams)*np.cos(lamf)*np.cos(dphi))  )
    #
    return np.absolute(dsigma)




#==============================================================================
# Script
#==============================================================================

#==============================================================================
# polyvectors = getPolyVecs()
# 
# 
# edge = 2
# lonmin = min(polylons)-edge
# lonmax = max(polylons)+edge
# latmin = min(polylats)-edge
# latmax = max(polylats)+edge
# 
# lonlist = np.arange(lonmin, lonmax+1, 1)
# latlist = np.arange(lonmax, latmax+1, 1)
# 
# mglon, mglat = np.meshgrid([round(x/10.0,2) for x in lonlist], [round(x/10.0,2) for x in latlist])
# 
# binpos = np.array([[(x, y) for x in lonlist] for y in latlist])
# 
# polymask = makePolyMask(binpos, polyvectors)
# numtest = 100
# testlons = np.random.rand(numtest)*(lonmax-lonmin)+lonmin
# testlats = np.random.rand(numtest)*(latmax-latmin)+latmin
# 
# testeqs = zip(testlons,testlats)
# goodeqs = []
# 
# for eq in testeqs:
#     if isinPoly(eq, polyvecs=polyvectors) == True:
#         goodeqs.append(eq)
#     
# goodlons = np.array([eq[0] for eq in goodeqs])
# goodlats = np.array([eq[1] for eq in goodeqs])
# 
# #==============================================================================
# # Display polymask bins and vectors on top of Basemap
# #==============================================================================
# plt.figure(figsize=(14, 10.5))
# m = Basemap(projection='cyl', llcrnrlat=latmin/10.0, urcrnrlat=latmax/10.0, llcrnrlon=lonmin/10.0, urcrnrlon=lonmax/10.0, resolution='i')
# m.drawcoastlines()
# m.drawmapboundary(fill_color='PaleTurquoise')
# m.fillcontinents(color='LemonChiffon',lake_color='PaleTurquoise', zorder=0)
# m.drawstates()
# m.drawcountries()
# m.drawparallels(np.arange(latmin/10.0, (latmax+1)/10.0, 0.1))#,labels=[1,0,0,0])
# m.drawmeridians(np.arange(lonmin/10.0, (lonmax+1)/10.0, 0.1))#,labels=[0,0,0,1])
# m.pcolor(mglon, mglat, polymask, alpha=0.5)
# #m.plot(testlons/10.0, testlats/10.0, 'm.')
# for vect in polyvectors:
#     m.plot([vect[0][0]/10.0, vect[1][0]/10.0], [vect[0][1]/10.0, vect[1][1]/10.0], 'g-')
# plt.show()
#==============================================================================
