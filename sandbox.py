# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 14:41:45 2015

@author: jmwilson
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import time

#==============================================================================
# Masked arrays
#==============================================================================
import RELM_inpoly as poly

edge = 2
lonmin = -1254
lonmax = -1131
latmin = 315
latmax = 430



lonlist = np.arange(lonmin-edge, lonmax+edge+1, 1)
latlist = np.arange(latmin-edge, latmax+edge+1, 1)

mglon, mglat = np.meshgrid(lonlist/10.0, latlist/10.0)

binlocs = np.array([[(x,y) for x in lonlist] for y in latlist])
unmask = np.array([[1.5 for x in lonlist] for y in latlist])

polyvectors = poly.getPolyVecs()
polymask = poly.makePolyMask(binlocs, polyvectors, np.shape(unmask[0][0]))

q = np.ma.masked_array(unmask, mask=polymask)
 

t0 = time.time()

#==============================================================================
# for j in xrange(len(q)):
#     for i in xrange(len(q[0])):
#         if polymask[j][i] == 0:
#             31**4
#==============================================================================
        
for row in q:
    for el in row:
        31**4

print time.time()-t0

#==============================================================================
# plotarray = []
# for row in q[~q.mask]:
#     for coord in row:
#         plotarray.append(coord)
#         
# 
# plt.figure(figsize=(16,12))
# plt.plot([coord[0] for coord in plotarray], [coord[1] for coord in plotarray], 'k.')
# #plt.pcolor(mglon, mglat, polymask)
#==============================================================================
#==============================================================================
# plt.show()
#==============================================================================
