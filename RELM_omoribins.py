# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 12:17:06 2015

@author: jmwilson

    Script to assign Omori-like earthquake halos in bins around earthquakes.
Also has function for populating a dictionary of quake bin coords
"""
#==============================================================================
# TODO: -Include ellipsoidal work from Mark's BASScast and RELM_basstests
#       -Try to use masked arrays to restrict calculations to RELM test region
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pickle
import time

#==============================================================================
# Using nearest 10th degree to bin events
#   - lat and lon are kept as int(coord*10)
#==============================================================================


def omori(quakes_f, binposit, histog, q=1.5, delm=1.0, mc=5.5, b=1.0, lam=1.76):
    flt = float
    convert = deg2km
    normalize = np.linalg.norm
            
    for quake in quakes_f:
        quakeloc = (flt(quake['hypo_lat'])*10, flt(quake['hypo_lon'])*10)
        mag = flt(quake['magnitude'])
        
        diff = binposit - quakeloc #CONVERT (LAT, LON) TO (KM, KM) BEFORE FINDING NORM
        diff = [[convert(pairs) for pairs in _] for _ in diff]
        rad = normalize(binposit - quakeloc, axis=2)
        
        Nom = 10**(b*(mag-delm-mc))    
        r0 = 0.5*10**(0.5*mag-lam)
        chi = r0**(1-q)/Nom/(q-1)
        dist = 1.0/chi/(r0 + rad)**q/(2*pi*rad)
        
        #==============================================================================
        # dist = 2**(1-q)*(q-1)*10**((b+q*0.5-0.5)*mag+(1-q)*lam-b*delm-b*mc)*(0.5*10**(0.5*mag-lam)+rad)**-q
        #==============================================================================
        
        histog += dist
    return None


def deg2km(degs):
    rE = 6371.0
    latkm = math.pi/180.0*rE*degs[0]
    lonkm = math.pi/180.0*rE*np.cos(math.pi/180.0*degs[0])*degs[1]
    kms = (latkm, lonkm)
    return kms
    
    
def dict_filler(save_loc = "/home/jmwilson/Desktop/RELM/cats", quakes_f):
    cat_rec = np.load('../cats/cat_rec_Feb-05-2015.p')
    hist_dict_f = {}
    flt = float
    intgr = int
    for quake in quakes_f:
        lon = intgr(flt(quake['hypo_lon'])*10)
        lat = intgr(flt(quake['hypo_lat'])*10)
        ind = (lon, lat)
        if ind not in hist_dict_f:
            hist_dict_f[ind] = [quake]
        else:
            hist_dict_f[ind].append(quake)
            
            
    f = open(save_loc + "/opt/hist/hist_dict_omori.p", "w")
    pickle.dump(cat_dict, f)
    f.close()
    
    return hist_dict_f

#==============================================================================
# Omori distribution for forecast generation
#   -everything is in integers until put into mesh
#==============================================================================
sim_loc = "/home/jmwilson/Desktop/RELM/sim_output"
sim_rec = np.load(file_loc+'/sim_rec.p')

edge = 50
lonmin = -12493
lonmax = -11457
latmin = 3156
latmax = 4220

lonlist = np.arange(lonmin-edge, lonmax+edge+1, 1)
latlist = np.arange(latmin-edge, latmax+edge+1, 1)

mglon, mglat = np.meshgrid([round(x/100.0,2) for x in lonlist], [round(x/100.0,2) for x in latlist])

binpos = np.array([[(y, x) for x in lonlist] for y in latlist])

hist_grid = np.zeros((len(latlist),len(lonlist)))

omori(sim_rec, binpos, hist_grid, mc=1.0)

rates_array = np.array([mglon, mglat, hist_grid])
rates_array.dump(file_loc+'/opt/hist/hist_grid_omori.p')

#==============================================================================
# Just fill dictionary with ANSS events
#==============================================================================
#==============================================================================
# cat_dict = dict_filler(cat_rec)
#==============================================================================
