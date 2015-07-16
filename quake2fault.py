# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:16:33 2015

@author: jmwilson
"""

import numpy as np
import math
import time
import h5py
import SimBASS_tools as tools


cat_rec = np.load('../cats/cat_rec_mc60_1980.p')
faults = h5py.File("../ALLCAL2_Faults/nocreep/ALLCAL2_VQ_faults_3km.h5",'r')

#==============================================================================
# Collect all Fault Vertices
#==============================================================================
vert_coords = []
for vert in faults['vertices']:
    vert_coords.append((vert[2],vert[1]))
vert_coords = np.array(vert_coords)

t0=time.time()
#==============================================================================
# Loop through quakes, get distance to all fault vertices, sort by distance, pick corresponding coordinate
#==============================================================================
for quake in cat_rec:
    quake_coord = np.array((quake['hypo_lon'], quake['hypo_lat']))
    
    dist = tools.sphericalDistNP(quake_coord, vert_coords)
    coords_n_dists = zip(vert_coords, dist)
    coords_n_dists = sorted(coords_n_dists, key=lambda x: x[1], reverse = False)
    best_coords = coords_n_dists[0][0]
        
    quake['hypo_lon'] = best_coords[0]
    quake['hypo_lat'] = best_coords[1]
        
    
print time.time()-t0

cat_rec.dump('../cats/cat_rec_mc60_1980_shifted_nocreep.p')