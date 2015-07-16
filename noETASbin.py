# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:09:01 2015

@author: jmwilson
"""

import numpy as np
import math
import time

#==============================================================================
# bin quakes without ETAS
#==============================================================================

def binning(quakes_f, binpos, mc=3.0):
    flt = float
    rnd = round
    intgr = int
    
    histog = np.zeros((len(binpos), len(binpos[0])))
    
    eqlons = np.array(quakes_f['hypo_lon']).astype(float)
    eqlats = np.array(quakes_f['hypo_lat']).astype(float)
    eqcoords = np.array(zip(eqlons, eqlats))
    
    lonoffset = np.min(binpos[...,0])
    latoffset = np.min(binpos[...,1])
    
#==============================================================================
#     offenders={}
#==============================================================================
    for quake in quakes_f:
        mag = flt(quake['magnitude'])
        quakeloc = (flt(quake['hypo_lon']), flt(quake['hypo_lat']))
        
        if mag >= mc and quakeloc != (0.0, 0.0):
        
            lon_ind = intgr(rnd(flt(quakeloc[0])*10)) - lonoffset
            lat_ind = intgr(rnd(flt(quakeloc[1])*10)) - latoffset
            
            histog[lat_ind][lon_ind] += 1
#==============================================================================
#             if quakeloc == (0.0, 0.0):
#                 if quake['sid'] not in offenders:
#                     offenders[quake['sid']] = [(quake['event_id'], quake['hypo_das'])]
#                 else:
#                     offenders[quake['sid']].append((quake['event_id'], quake['hypo_das']))
#    print offenders
#==============================================================================
        
        
    return histog
    


simind = 3
folderext=['/allcal/version_1a', '/virtcal/z_01', '/pollitz/version_1a', '/rsqsim/version_1a']

sim_loc = "/home/jmwilson/Desktop/raw_output/allcal2" + folderext[simind]
sim_rec = np.load(sim_loc+'/ETAS_files/sim_rec3.p')

edge = 2
lonmin = -1254
lonmax = -1131
latmin = 315
latmax = 430

lonlist = np.arange(lonmin-edge, lonmax+edge+1, 1)
latlist = np.arange(latmin-edge, latmax+edge+1, 1)

binpos = np.array([[(x, y) for x in lonlist] for y in latlist])

t0=time.time()
hist_grid = binning(sim_rec, binpos, mc=3.0)
print 'time in minutes: ', (time.time()-t0)/60.0


rates_array = np.array([lonlist, latlist, hist_grid])
rates_array.dump(sim_loc+'/ETAS_files/hist_grid_noETAS6.p')