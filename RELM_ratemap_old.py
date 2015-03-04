# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 12:11:32 2015

@author: jmwilson
"""

#==============================================================================
# TODO: *Speed up
#==============================================================================
import numpy as np

events_rec = np.load('/home/jmwilson/Desktop/raw_output/allcal2/allcal/version_1a/events_rec.p')
cat_rec = np.load('../cats/cat_rec_30yr.p')

quakes_rec = events_rec
quake_coords = zip(quakes_rec['hypo_lon'].astype(np.float), quakes_rec['hypo_lat'].astype(np.float))

lonmin = min(events_rec['hypo_lon'].astype(np.float))
lonmax = max(events_rec['hypo_lon'].astype(np.float))
latmin = min(events_rec['hypo_lat'].astype(np.float))
latmax = max(events_rec['hypo_lat'].astype(np.float))

for binsize in [0.1]:#np.linspace(0.05, 0.25, 21):
    #binsize = 0.1
    staggering = 10
    binsperdeg = staggering / binsize
    
    lonsize = int((lonmax-lonmin)*binsperdeg)
    latsize = int((latmax-latmin)*binsperdeg)
    
    lonlist = np.linspace(lonmin, lonmax, lonsize)
    latlist = np.linspace(latmin, latmax, latsize)
    mglon, mglat = np.meshgrid(lonlist, latlist)
    
    
    #Binning
    quake_rates = np.zeros((latsize,lonsize))
    for i in range(lonsize):
        print "lon: %i/%i"%(i+1, lonsize)
        for j in range(latsize):
            #print "%ith lat"%j
            for coord in quake_coords:
                londiff = abs(coord[0] - lonlist[i])
                latdiff = abs(coord[1] - latlist[j])
                if londiff <= binsize/2.0 and latdiff <= binsize/2.0:
                    #print "add a quake"
                    quake_rates[j][i] += 1
    
    rates_array = np.array([mglon, mglat, quake_rates])
    #rates_rec = np.core.records.fromarrays(rates_array.transpose(), names="mglon, mglat, quake_rates")
    rates_array.dump('/home/jmwilson/Desktop/raw_output/allcal2/allcal/version_1a/sim_rates_%i.p'%(int(binsize*100)))