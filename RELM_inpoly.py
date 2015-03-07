# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:07:35 2015

@author: jmwilson
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

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


def makePolyMask(binpositions, polyvecs, maskshape):
    #
    masklist = np.zeros((np.shape(binpositions)[:-1])  + maskshape)
    #
    for j in xrange(len(binpositions)):
        for i in xrange(len(binpositions[0])):
            inarea = isinPoly(binpositions[j][i], polyvecs=polyvecs)
            if inarea == True:
                pass
            else:
                masklist[j][i] = np.ones(maskshape)
    #
    return np.array(masklist)
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
# #==============================================================================
# # binpos = np.array([[(x, y) for x in lonlist] for y in latlist])
# # 
# # polymask = makePolyMask(binpos)
# #==============================================================================
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
# #m.pcolor(mglon, mglat, polymask, alpha=0.5)
# m.plot(testlons/10.0, testlats/10.0, 'm.')
# for vect in polyvectors:
#     m.plot([vect[0][0]/10.0, vect[1][0]/10.0], [vect[0][1]/10.0, vect[1][1]/10.0], 'g-')
#     
# plt.figure(figsize=(14, 10.5))
# m = Basemap(projection='cyl', llcrnrlat=latmin/10.0, urcrnrlat=latmax/10.0, llcrnrlon=lonmin/10.0, urcrnrlon=lonmax/10.0, resolution='i')
# m.drawcoastlines()
# m.drawmapboundary(fill_color='PaleTurquoise')
# m.fillcontinents(color='LemonChiffon',lake_color='PaleTurquoise', zorder=0)
# m.drawstates()
# m.drawcountries()
# m.drawparallels(np.arange(latmin/10.0, (latmax+1)/10.0, 0.1))#,labels=[1,0,0,0])
# m.drawmeridians(np.arange(lonmin/10.0, (lonmax+1)/10.0, 0.1))#,labels=[0,0,0,1])
# m.plot(goodlons/10.0, goodlats/10.0, 'm.')
# for vect in polyvectors:
#==============================================================================
#==============================================================================
#     m.plot([vect[0][0]/10.0, vect[1][0]/10.0], [vect[0][1]/10.0, vect[1][1]/10.0], 'g-')
# 
# plt.show()
#==============================================================================
