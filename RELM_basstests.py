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

def omori_el(quakes_f, binposit, size, abratio=2, strike=30, q=1.5, mag = 7.0, delm=1.0, mc=5.5, b=1.0, lam=1.76):
    flt = float
    convert = deg2km
    normalize = np.linalg.norm
    cosine = np.cos
    sine = np.sin
    
    histog = np.zeros((len(binposit), len(binposit[0])))
    
    angledeg = -strike
    anglerad = math.pi/180.0*angledeg
    
    for quakeloc in quakes_f:
        
        diff = binposit - quakeloc
        diff = [[convert(pairs) for pairs in _] for _ in diff]
                
        transform = [[((coord[0]*cosine(anglerad)+coord[1]*sine(anglerad))*abratio/1,( -coord[0]*sine(anglerad)+coord[1]*cosine(anglerad))/1) for coord in _] for _ in diff]        
                
        rad = normalize(transform, axis=2)
        
        Nom = 10**(b*(mag-delm-mc))    
        r0 = 0.5*10**(0.5*mag-lam)
        chi = r0**(1-q)/Nom/(q-1)
        dist = 1.0/chi/(r0 + rad)**q#/(2*pi*rad)
        
        #dist = Nom/np.sum(dist)*dist
        #dist = np.ma.
        
        histog += dist
    
    return histog


def deg2km(degs):
    rE = 6371.0
    latkm = math.pi/180.0*rE*degs[0]
    lonkm = math.pi/180.0*rE*np.cos(math.pi/180.0*degs[0])*degs[1]
    kms = (latkm, lonkm)
    return kms


gridsize = 200
numquakes = 1

quakes = [(gridsize/2+0.3, gridsize/2+0.3)]#np.random.random_sample((numquakes, 2))*gridsize#np.random.randint(0, gridsize,(numquakes,2))#
binlocs = np.array([[(x,y) for x in xrange(gridsize)] for y in xrange(gridsize)])

hist_grid = omori_el(quakes, binlocs, gridsize)

plt.figure()
#plt.imshow(hist_grid, interpolation='nearest')
plt.pcolor(hist_grid)
plt.colorbar()
plt.show()