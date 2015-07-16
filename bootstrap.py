# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:06:58 2015

@author: jmwilson
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import SimBASS_tools as tools

simind = 3
etasBool = 1

sim_fold = "/home/jmwilson/Desktop/raw_output/allcal2"
folderext=['/allcal/version_1a', '/virtcal/z_01', '/pollitz/version_1a', '/rsqsim/version_1a']
sim_loc = sim_fold+folderext[simind]

cat_loc = "/home/jmwilson/Desktop/SimBASS/cats"
boot_total = []

if etasBool == True:
    simdat = open(sim_loc + '/ETAS_files/hist_grid_omori6.p', 'r')
    catdat = open(cat_loc+'/hist_dict6.p', 'r')#
else:
    simdat = open(sim_loc + '/ETAS_files/hist_grid_noETAS6.p', 'r')
    catdat = open(cat_loc+'/hist_dict6_shifted_creep.p', 'r')#    

sim_info = pickle.load(simdat)
simdat.close()
#
cat_hist_dict = pickle.load(catdat)
catdat.close()

lonlist = np.array(sim_info[0])
latlist = np.array(sim_info[1])

sim_grid = sim_info[2]

binlocs = np.array([[(x,y) for x in lonlist] for y in latlist])

if etasBool == True:
    binmask = tools.makePolyMask(binlocs)
else:
    binmask = tools.makeFaultMask(binlocs)

#==============================================================================
# Sample ETAS rate field to produce new quake catalog
#==============================================================================
#==============================================================================
# areasize = 0
# sim_flat = []
# for j in xrange(len(binmask)):
#     for i in xrange(len(binmask[0])):
#         if binmask[j][i] == 0:
#             areasize += 1
#             sim_flat.append([(binlocs[j][i][0], binlocs[j][i][1]), sim_grid[j][i]])
# 
# sim_sorted = sorted(sim_flat, key=lambda x: x[1], reverse = True)
# 
# sim_sorted = np.array(sim_sorted)
# 
# flat = sim_sorted[...,1]
# pdf = flat/np.sum(flat)
# cdf = np.cumsum(pdf)
# 
# for q in range(1000):
#     loopboot = {}
#     for p in range(sum(cat_hist_dict.values())):#range(len(sim_flat)):#
#         randnum = np.random.rand()
#         #bestdiff = 2.0
#         for i in range(len(cdf)):
#             #diff = abs(randnum-cdf[i])
#             if cdf[i] >= randnum:#diff < bestdiff:#
#                 #bestdiff = diff
#                 besti = i
#                 break
#             #
#         coordkey = sim_sorted[besti][0]
#         if coordkey not in loopboot:
#             loopboot[coordkey] = 1
#         else:
#             loopboot[coordkey]+=1
#         #
#     boot_total.append(loopboot)
# 
# f = open(sim_loc + "/ETAS_files/hist_dict_bootstrap.p", "w")
# pickle.dump(boot_total, f)
# f.close()
#==============================================================================


#==============================================================================
# Bootstrap on ANSS catalog
#==============================================================================
allquakes = []
for key in cat_hist_dict:
    for i in range(cat_hist_dict[key]):
        allquakes.append(key)
    #

for q in range(1000):
    loopboot = {}
    for p in range(sum(cat_hist_dict.values())):
        randind = np.random.randint(len(allquakes))
        coordkey = allquakes[randind]
        if coordkey not in loopboot:
            loopboot[coordkey] = 1
        else:
            loopboot[coordkey] += 1
        #
    boot_total.append(loopboot)
#
f = open(sim_loc + "/ETAS_files/hist_dict_bootstrap2.p", "w")
pickle.dump(boot_total, f)
f.close()


