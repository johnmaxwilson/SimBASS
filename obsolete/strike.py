# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:57:39 2015

@author: jmwilson
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
import time
import SimBASS_tools as tools


#==============================================================================
# TRY TO MAKE THIS THING VECTORIZED
#==============================================================================

deg2rad = math.pi/180.0


eqlons = np.random.rand(400)
eqlats = np.random.rand(400)
quakes = np.array(zip(eqlons, eqlats))
for quakeloc in quakes[:1]:    
    #==============================================================================
    # calculate the strike of a single quake "quakeloc"
    #==============================================================================
    ruptlen = 0.5#"CURRENT EQ RUPTURE LENGTH" #r0ssim*2.0
    
    
    tmpX = []
    tmpY = []
    
    distresults = tools.sphericalDistNP(quakeloc, eqlons, eqlats)

    for row in distresults:
        if row[0] <= 5.0*ruptlen:
            tmpX += [row[1]]
            tmpY += [row[2]]
    #
    if len(tmpX)<=1:
        abratio=1.0
        strike=0.0
    #    
    if len(tmpX)>=2:
        abratio = 2.0
        fitparams = tools.linfit(p=scipy.array([0.,0.]), X=tmpX, Y=tmpY, full_output=False)
        strike = 90.0 - math.atan(fitparams[0][1])/deg2rad

    print quakeloc, strike



#plt.close("all")
plt.figure()
plt.plot(eqlons, eqlats, 'k.')
plt.plot(quakeloc[0],quakeloc[1], 'mo')
