# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 12:17:06 2015

@author: jmwilson
"""
#==============================================================================
# TODO: Explicitly convert from km to lat/lon degrees
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pickle
import time

#==============================================================================
# Using nearest 100th degree to bin events
#   - lat and lon are kept as int(coord*100)
#==============================================================================

def dict_filler(quakes_f, extend):
    hist_dict_f = {}
    flt = float
    intgr = int
    rnd = round
    xrng = xrange
    if extend == True:
        for quake in quakes_f:
            lon0 = intgr(flt(quake['hypo_lon'])*100)
            lat0 = intgr(flt(quake['hypo_lat'])*100)
            mag = flt(quake['magnitude'])
            #Since for CA 1 km ~ 0.01 degrees, calc km rupture distance and leave as is.  could be improved!
            radius = intgr(rnd(10**(-3.22+0.69*mag)*0.5))
            for lonoffset in xrng(-radius, radius+1):
                for latoffset in xrng(-radius, radius+1):
                    if lonoffset*lonoffset + latoffset*latoffset <= radius*radius:                
                        ind = (lon0+lonoffset, lat0+latoffset)
                        if ind in hist_dict_f:
                            hist_dict_f[ind].append(quake)
                        else:
                            hist_dict_f[ind] = [quake]
    else:
        for quake in quakes_f:
            lon0 = intgr(flt(quake['hypo_lon'])*100)
            lat0 = intgr(flt(quake['hypo_lat'])*100)
            mag = flt(quake['magnitude'])
            for lonoffset in xrng(-5, 6):
                for latoffset in xrng(-5, 6):             
                    ind = (lon0+lonoffset, lat0+latoffset)
                    if ind in hist_dict_f:
                        hist_dict_f[ind].append(quake)
                    else:
                        hist_dict_f[ind] = [quake]
    return hist_dict_f

spacialextent = True

file_loc = "/home/jmwilson/Desktop/RELM/sim_output"
quakes_rec = np.load(file_loc+'/sim_rec.p')
#quakes_rec = np.load('../cats/cat_rec_Feb-05-2015.p')


hist_dict = dict_filler(quakes_rec, spacialextent)



f = open(file_loc+ "/opt/hist/hist_dict_variable.p", "w")
pickle.dump(hist_dict, f)
f.close()

#==============================================================================
# Full grid filling for contour plot
#   -everything is in integers until put into mesh
#==============================================================================
size = 50
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
rates_array.dump(file_loc+'/opt/hist/hist_grid_variable.p')
