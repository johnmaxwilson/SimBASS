# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 14:41:45 2015

@author: jmwilson
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import time

def omori(quakelist, binposit, histog, mag, q=1.5, delm=1.0, mc=5.5, b=1.0, lam=1.76):
    
    convert = deg2km
    normalize = np.linalg.norm
    
    for quakeloc in quakelist:
        diff = binposit - quakeloc
        
        diff = [[convert(pairs) for pairs in _] for _ in diff]
                
        rad = normalize(diff, axis=2)        
        
        Nom = 10**(b*(mag-delm-mc))    
        r0 = 0.5*10**(0.5*mag-lam)
        chi = r0**(1-q)/Nom/(q-1)
        dist = 1.0/chi/(r0 + rad)**q/(2*pi*rad)
        
        #==============================================================================
        # dist = 2**(1-q)*(q-1)*10**((b+q*0.5-0.5)*mag+(1-q)*lam-b*delm-b*mc)*(0.5*10**(0.5*mag-lam)+rad)**-q
        #==============================================================================
                
        histog += dist
    return None
    
#==============================================================================
# def omorilin(rad, mag, q=1.5, delm=1.0, mc=1.0, b=1.0, lam=1.76):
#     Nom = 10**(b*(mag-delm-mc))
#     r0 = 0.5*10**(0.5*mag-lam)
#     chi = r0**(1-q)/Nom/(q-1)
#     dist = 1.0/chi/(r0 + rad)**q
#     
#     return dist
# 
# x = np.linspace(0,5000,1000)
# y = omori(x,9.0, mc=4.5)
# 
# plt.figure()
# plt.loglog(x, y)
# plt.show()
#==============================================================================
    
def deg2km(degs):
    rE = 6371.0
    latkm = pi/180.0*rE*degs[0]
    lonkm = pi/180.0*rE*np.cos(pi/180.0*degs[0])*degs[1]
    kms = (latkm, lonkm)
    return kms

numquakes = 10
gridsize = 500

hist = np.zeros((gridsize,gridsize))


binpos = np.array([[(x,y) for x in xrange(gridsize)] for y in xrange(gridsize)])
quakes = np.random.random_sample((numquakes,2))*gridsize
t0=time.time()

omori(quakes, binpos, hist, 7.1, mc=3.0)

print time.time() - t0

plt.figure()
plt.imshow(hist, interpolation="nearest")
plt.colorbar()
plt.show()