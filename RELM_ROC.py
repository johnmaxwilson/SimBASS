# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 10:02:55 2015

@author: jmwilson
"""
#==============================================================================
# TODO: Consider restrictions on simulation bins (have UCERF2 faults?)
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pickle

def ROC_fnc(sim_bins, cat_dict, maxbins):
    A=0.0
    B=0.0
    AplusC = float(len(cat_dict))
    BplusD = float(len(sim_bins)-AplusC) #Here simbins should be considered
    Hits_f = []
    Falsies_f = []
    h_apnd = Hits_f.append
    f_apnd = Falsies_f.append
    
    for N_bins in xrange(maxbins):
        #print sim_bins[N_bins][0]
        if sim_bins[N_bins][2] == 1:
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

sim_loc = "/home/jmwilson/Desktop/RELM/BASScast/kml"
cat_loc = "/home/jmwilson/Desktop/RELM/cats"
finaldat = []


simdat = open(sim_loc + '/virtcal-hist_grid_bassomori.p', 'r')
sim_info = pickle.load(simdat)
lonmin = int(min(sim_info[0][0])*10)
latmin = int(min(sim_info[1][0])*10)
sim_grid = sim_info[2]
simdat.close()

catdat = open(cat_loc+'/opt/hist/hist_dict_omori.p', 'r')
cat_hist_dict = pickle.load(catdat)
catdat.close()

# Sort simulation by number of events in bin, only include bins with either sim or cat events
sim_flat = [[(i+lonmin, j+latmin), sim_grid[j][i]] for j in xrange(len(sim_grid)) for i in xrange(len(sim_grid[0]))]# if sim_grid[j][i]>0 or (i+lonmin, j+latmin) in cat_hist_dict]
sim_flat = [line+[1] if line[0] in cat_hist_dict else line+[0] for line in sim_flat]

sim_sorted = sorted(sim_flat, key=lambda x: x[2], reverse = True)
sim_sorted = sorted(sim_sorted, key=lambda x: x[1], reverse = True)


N_binsmax = len(sim_sorted)

#ROC curve creation
Hits, Falsies = ROC_fnc(sim_sorted, cat_hist_dict, N_binsmax)

#ROC area skill test
df = [Falsies[i+1]-Falsies[i] for i in xrange(len(Falsies)-1)]
df = np.array(df)
score = sum(Hits[:-1]*df)-0.5

finaldat.append({'Hits':Hits, 'Falsies':Falsies, 'score':score})
    
#Dumping data
f = open(sim_loc+"/virtcal-ROC_data_omori.p", "w")
pickle.dump(finaldat, f)
f.close
