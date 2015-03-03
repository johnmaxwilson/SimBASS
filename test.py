
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 14:55:53 2014

@author: jmwilson
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.basemap import Basemap

sidesize = 200
hist_grid = np.zeros((sidesize, sidesize))
num = 1

abratio = 2
strike = 30
angledeg = -strike
anglerad = math.pi/180.0*angledeg

mag=7.0
q=1.5
delm=1.0
mc=5.5
b=1.0
lam=1.76

centers = [(sidesize/2+0.3, sidesize/2+0.3)]#np.random.random_sample((num, 2))*sidesize#np.random.randint(0, sidesize,(num,2))#
binlocs = np.array([[(x,y) for x in xrange(sidesize)] for y in xrange(sidesize)])


for center in centers:
    diff = binlocs - center
    transform = [[((coord[0]*math.cos(anglerad)+coord[1]*math.sin(anglerad))*abratio/1,( -coord[0]*math.sin(anglerad)+coord[1]*math.cos(anglerad))/1) for coord in _] for _ in diff]
    rad = np.linalg.norm(transform, axis=2)
    #dist = 1.0/(2*np.pi*10**2)*np.exp(-(rad**2)/(2.0*10**2))
    
    
    Nom = 10**(b*(mag-delm-mc))    
    r0 = 0.5*10**(0.5*mag-lam)
    chi = r0**(1-q)/Nom/(q-1)
    dist = 1/chi*(r0 + rad)**-q/(2*pi*rad)
    
    #dist = Nom/np.sum(dist)*dist
    dist = np.ma.    
    
    hist_grid += dist
    print Nom

plt.figure()
#plt.imshow(hist_grid, interpolation='nearest')
plt.pcolor(hist_grid)
plt.colorbar()
plt.show()
print np.sum(hist_grid)

#==============================================================================
# mag=9.0
# q=1.5
# delm=1.0
# mc=4.5
# b=1.0
# lam=1.76
# r0 = 0.5*10**(0.5*mag-lam)
# Nom = 10**(b*(mag-delm-mc))
# 
# rad = np.arange(1,100,0.1)
# dist = (r0 + rad)**-q/rad
# dist = Nom/np.sum(dist)*dist
# 
# plt.figure()
# plt.loglog(rad, dist)
#==============================================================================

plt.show()