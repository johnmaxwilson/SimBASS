# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 14:55:53 2014

@author: jmwilson
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import time
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
#==============================================================================
#==============================================================================
# plt.show()
#==============================================================================

test = [[1,2],[3,4]]
tmask = [[0,1],[0,1]]

masked = np.ma.masked_array(test, mask = tmask)