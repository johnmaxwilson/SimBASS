# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:27:37 2015

@author: jmwilson

    Test of the mechanisms to quickly (numpy) assign each bin in a grid an 
earthquake rate based on Omori behavior for a list of quake epicenters.
    Does not use real lat/lons. Uses single magnitude and strike for all events
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.basemap import Basemap
import RELM_inpoly as poly

def omori_el(quakes_f, binpos, abratio=2, strike=30, q=1.5, mag = 7.0, delm=1.0, mc=5.5, b=1.0, lam=1.76):
    flt = float
    convert = deg2km
    normalize = np.linalg.norm
    cosine = np.cos
    sine = np.sin
    
    histog = np.zeros((len(binpos), len(binpos[0])))
    
    angledeg = strike
    anglerad = math.pi/180.0*angledeg
    
    for quakeloc in quakes_f:
        
        diff = binpos - quakeloc
        diff = [[convert(pairs) for pairs in _] for _ in diff]
                
        transform = [[((coord[0]*cosine(anglerad)-coord[1]*sine(anglerad))*abratio,(coord[0]*sine(anglerad)+coord[1]*cosine(anglerad))) for coord in _] for _ in diff]        
        
        rad = normalize(transform, axis=2)
        
        Nom = 10**(b*(mag-delm-mc))    
        r0 = 0.5*10**(0.5*mag-lam)
        chi = r0**(1-q)/Nom/(q-1)
        dist = 1.0/chi/(r0 + rad)**q/(2*pi*rad)
        
        #dist = Nom/np.sum(dist)*dist
        dist = np.ma.fix_invalid(dist, fill_value=0).data
        
        histog += dist
    
    return histog


def deg2km(degs):
    rE = 6371.0
    latkm = math.pi/180.0*rE*degs[1]/10.0
    lonkm = math.pi/180.0*rE*np.cos(math.pi/180.0*degs[1]/10.0)*degs[0]/10.0
    kms = (lonkm, latkm)
    return kms


edge = 2
lonmin = -1254
lonmax = -1131
latmin = 315
latmax = 430

lonlist = np.arange(lonmin-edge, lonmax+edge+1, 1)
latlist = np.arange(latmin-edge, latmax+edge+1, 1)

mglon, mglat = np.meshgrid(lonlist/10.0, latlist/10.0)

binlocs = np.array([[(x, y) for x in lonlist] for y in latlist])

numquakes = 100

eqlons = np.random.rand(numquakes)*(lonmax-lonmin)+lonmin#[0.5*(lonmax-lonmin)+lonmin+0.5]#
eqlats = np.random.rand(numquakes)*(latmax-latmin)+latmin#[0.5*(latmax-latmin)+latmin+0.5]#
eqlocs = zip(eqlons, eqlats)


t0=time.time()
hist_grid = omori_el(eqlocs, binlocs)
print time.time()-t0

plt.figure(figsize=(14, 10.5))
plt.pcolor(mglon, mglat, hist_grid)
plt.show()

#==============================================================================
# polyvectors = poly.getPolyVecs()
# polymask = poly.makePolyMask(binlocs, polyvectors, ())
# hist_grid_ma = np.ma.masked_array(hist_grid, mask=polymask)
# 
# plt.figure(figsize=(14, 10.5))
# m = Basemap(projection='cyl', llcrnrlat=latmin/10.0, urcrnrlat=latmax/10.0, llcrnrlon=lonmin/10.0, urcrnrlon=lonmax/10.0, resolution='i')
# #m.drawmapboundary(fill_color='PaleTurquoise')
# #m.fillcontinents(color='lemonchiffon',lake_color='PaleTurquoise', zorder=0)
# m.drawcoastlines()
# m.drawstates()
# m.drawcountries()
# m.drawparallels(np.arange(latmin/10.0, (latmax+1)/10.0, 2),labels=[1,0,0,0])
# m.drawmeridians(np.arange(lonmin/10.0, (lonmax+1)/10.0, 2),labels=[0,0,0,1])
# m.pcolor(mglon, mglat, hist_grid_ma)
#==============================================================================
#==============================================================================
# plt.show()
#==============================================================================
