# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 14:41:45 2015

@author: jmwilson
"""
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import time
import cPickle
import datetime as dtm
import h5py
import SimBASS_tools as tools
from scipy.spatial import KDTree
import matplotlib.animation as animation

#==============================================================================
# #==============================================================================
# # Using Rtree for spatial binning
# #==============================================================================
# from rtree import index
# from scipy.spatial import KDTree
# 
# binsize = 1
# 
# lonbounds=[0,10] #[-125.4, -113.1]
# latbounds=[0,10] #[31.5, 43.0]
# 
# lonborders = np.arange(lonbounds[0],lonbounds[1]+binsize, binsize)
# latborders = np.arange(latbounds[0],latbounds[1]+binsize, binsize)
# 
# borders = []
# for i in range(len(lonborders)-1):
#     for j in range(len(latborders)-1):
#         borders.append((lonborders[i], lonborders[i+1], latborders[j], latborders[j+1]))
# 
# bins = index.Index(interleaved=False)
# count = 0
# for box in borders:
#     bins.insert(count, box, obj=count)
#     count += 1
#     
# hits = bins.intersection((0,1.5))
# 
# for ind in hits:
#     print ind.obj()
# 
# #==============================================================================
# #  P(A)*A = const?
# #==============================================================================
# #==============================================================================
# # simind = 1
# # etasBool = 1
# # bootstrapping = 1
# # 
# # filesuff = ['noETAS6.p','omori6.p']
# # 
# # sim_fold = "/home/jmwilson/Desktop/raw_output/allcal2"
# # folderext=['/allcal/version_1a', '/virtcal/z_01', '/pollitz/version_1a', '/rsqsim/version_1a']
# # sim_loc = sim_fold+folderext[simind]
# # 
# # simdat = open(sim_loc + '/ETAS_files/hist_grid_'+filesuff[etasBool], 'r')
# # sim_info = pickle.load(simdat)
# # simdat.close()
# # 
# # lonlist = np.array(sim_info[0])
# # latlist = np.array(sim_info[1])
# # 
# # sim_grid = np.array(sim_info[2])/np.sum(sim_info[2])
# # 
# # n, bins = np.histogram(sim_grid, bins=np.linspace(0,sim_grid.max(), 1001))
# # plt.close()
# # plt.bar(bins[:-1], n*bins[:-1], width = (bins[1]-bins[0]))#(bins[:-1]+(bins[1]-bins[0]/2.0)), n*(bins[:-1]+(bins[1]-bins[0]/2.0)))
#==============================================================================
# plt.show
#==============================================================================

fig = plt.figure()
heights = np.ones((25,25))

areas = np.ones(heights.shape())
volume = heights/areas
slopes = np.zeros(heights.shape()+(2))
momentum = momconst*slopes
velocity = momentum/volume

for step in range(100):
    for i in range(1,len(heights)-1):
        for j in range(1,len(heights[0])-1):
            slopex = heights[i,j+1] - heights[i,j-1]
            slopey = heights[i+1,j] - heights[i-1,j]
            slopes[i][j] = (slopex, slopey)
    
    momentum = momconst*slopes/volume
    velocity = momentum/volume
    
    im = plt.imshow(heights)
    ims.append([im])
    
    