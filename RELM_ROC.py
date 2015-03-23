# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 10:02:55 2015

@author: jmwilson
"""
#==============================================================================
# TODO: Consider restrictions on simulation bins (have UCERF2 faults?)
#       -restrict to calipoly
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pickle
import RELM_inpoly as poly

def ROC_fnc(sim_bins, cat_dict, areasize):
    A=0.0
    B=0.0
    AplusC = float(len(cat_dict))
    BplusD = float(areasize-AplusC)
    Hits_f = []
    Falsies_f = []
    h_apnd = Hits_f.append
    f_apnd = Falsies_f.append
    
    for sbin in sim_bins:
        #print sim_bins[N_bins][0]
        if sbin[2] == 1:
            A += 1.0
        else:
            B += 1.0
        
        hit = A/AplusC
        false = B/BplusD
        h_apnd(hit)
        f_apnd(false)
    
    Hits_f = np.array(Hits_f)
    Falsies_f = np.array(Falsies_f)
    return Hits_f, Falsies_f

sim_loc = "/home/jmwilson/Desktop/RELM/sim_output"
cat_loc = "/home/jmwilson/Desktop/RELM/cats"
finaldat = []


simdat = open(sim_loc + '/opt/hist/hist_grid_omori.p', 'r')
sim_info = pickle.load(simdat)
simdat.close()

catdat = open(cat_loc+'/opt/hist/hist_dict_omori.p', 'r')
cat_hist_dict = pickle.load(catdat)
catdat.close()

#lonmin = int(round(min(sim_info[0][0])*10))
#latmin = int(round(min(sim_info[1][0])*10))
lonlist = np.array(sim_info[0])
#lonlist = lonlist.round()
#lonlist = lonlist.astype(int)
latlist = np.array(sim_info[1])
#latlist = latlist.round()
#latlist = latlist.astype(int)

sim_grid = sim_info[2]

binlocs = np.array([[(x,y) for x in lonlist] for y in latlist])
polymask = poly.makePolyMask(binlocs)

#sim_masked = np.ma.masked_array(sim_grid, mask=polymask)

#==============================================================================
# make flattened list of bins in calipoly, sort by forecast value and whether it's a hit
#==============================================================================
areasize = 0
sim_flat = []
for j in xrange(len(polymask)):
    for i in xrange(len(polymask[0])):
        if polymask[j][i] == 0:
            areasize += 1
            sim_flat.append([(binlocs[j][i][0], binlocs[j][i][1]), sim_grid[j][i]])

#sim_flat = [[(i+min(lonlist), j+min(latlist)), sim_grid[j][i]] for j in xrange(len(sim_grid)) for i in xrange(len(sim_grid[0]))]# if sim_grid[j][i]>0 or (i+lonmin, j+latmin) in cat_hist_dict]

sim_flat = [line+[1] if line[0] in cat_hist_dict else line+[0] for line in sim_flat]

sim_sorted = sorted(sim_flat, key=lambda x: x[2], reverse = True)
sim_sorted = sorted(sim_sorted, key=lambda x: x[1], reverse = True)


#ROC curve creation
Hits, Falsies = ROC_fnc(sim_sorted, cat_hist_dict, areasize)

#ROC area skill test
df = [Falsies[i+1]-Falsies[i] for i in xrange(len(Falsies)-1)]
df = np.array(df)
score = sum(Hits[:-1]*df)-0.5
print score

finaldat.append({'Hits':Hits, 'Falsies':Falsies, 'score':score})

#Dumping data
f = open(sim_loc+"/virtcal-ROC_data_omori.p", "w")
pickle.dump(finaldat, f)
f.close