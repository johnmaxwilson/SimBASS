# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 14:41:45 2015

@author: jmwilson
"""
#==============================================================================
# TODO: 
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pickle
import time

#==============================================================================
# Using nearest 100th degree to bin events
#   - lat and lon are kept as int(coord*100)
#==============================================================================

def dict_filler(quakes_f, binsize):
    hist_dict_f = {}
    offsets = np.arange(-binsize/2, binsize/2 + 1, 1)
    offsets = [(x, y) for x in offsets for y in offsets]
    flt = float
    intgr = int
    for quake in quakes_f:
        lon = intgr(flt(quake['hypo_lon'])*100)
        lat = intgr(flt(quake['hypo_lat'])*100)
        for offset in offsets:
            ind = (lon+offset[0], lat+offset[1])
            if ind not in hist_dict_f:
                hist_dict_f[ind] = [quake]
            else:
                hist_dict_f[ind].append(quake)
    return hist_dict_f

file_loc = "/home/jmwilson/Desktop/SimBASS/cats"
sim_rec = np.load(file_loc+'/sim_output/sim_rec.p')
#cat_rec = np.load('../cats/cat_rec_Feb-05-2015.p')

quakes_rec = sim_rec


for size in xrange(10, 11, 2):
    print size
    hist_dict = dict_filler(quakes_rec, size)
    

    
    f = open(file_loc+ "/opt/hist/hist_dict_%i.p"%size, "w")
    #f = open("../cats/cat_hist_dict.p", "w")
    pickle.dump(hist_dict, f)
    f.close()
    
    #==============================================================================
    # Full grid filling for contour plot
    #   -everything is in integers until put into mesh
    #==============================================================================
    
    #lonmin = int(min(sim_rec['hypo_lon'].astype(np.float))*100)
    lonmin = -12493
    lonmax = -11457
    latmin = 3156
    latmax = 4220
    
    lonlist = np.arange(lonmin-size, lonmax+size+1, 1)
    latlist = np.arange(latmin-size, latmax+size+1, 1)
    mglon, mglat = np.meshgrid([round(x/100.0,2) for x in lonlist], [round(x/100.0,2) for x in latlist])
    
    hist_grid = np.zeros((len(latlist),len(lonlist)), dtype=np.int)
    for key in hist_dict:
        hist_grid[key[1]-(latmin-size)][key[0]-(lonmin-size)] = len(hist_dict[key])
    
    rates_array = np.array([mglon, mglat, hist_grid])
    rates_array.dump(file_loc+'/opt/hist/hist_grid_%i.p'%size)
