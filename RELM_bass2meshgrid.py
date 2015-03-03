# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 10:53:10 2015

@author: jmwilson
"""

import numpy as np

fileloc = '/home/jmwilson/Desktop/RELM/BASScast/kml'
simulator = 'virtcal'

z    = []
lons = []
lats = []

f = open(fileloc+'/'+simulator+'-ETASgrid.xyz', 'r')
for line in f:
    linelist = line.split()
    if linelist[0][0] != '#':
        z.append(float(linelist[2]))
        lon = float(linelist[0])
        lat = float(linelist[1])
        if lon not in lons:
            lons.append(lon)
        if lat not in lats:
            lats.append(lat)
f.close()

#==============================================================================
# lons = np.array(lons)
# lats = np.array(lats)
#==============================================================================

count = 0
hist_grid = np.zeros((len(lats), len(lons)))
for j in xrange(len(lats)):
    for i in xrange(len(lons)):
        hist_grid[j][i] = z[count]
        count += 1
        
        
mglon, mglat = np.meshgrid(lons, lats)

rates_array = np.array([mglon, mglat, hist_grid])
rates_array.dump(fileloc+'/'+simulator+'-hist_grid_bassomori.p')