# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 10:02:55 2015

@author: jmwilson
"""
#==============================================================================
# TODO: 
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pickle
import SimBASS_tools as tools

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

simind = 3
etasBool = 1
bootstrapping = 0

filesuff = ['noETAS6.p','omori6_abrat4.p']

sim_fold = "/home/jmwilson/Desktop/raw_output/allcal2"
folderext=['/allcal/version_1a', '/virtcal/z_01', '/pollitz/version_1a', '/rsqsim/version_1a']
sim_loc = sim_fold+folderext[simind]

simdat = open(sim_loc + '/ETAS_files/hist_grid_'+filesuff[etasBool], 'r')
sim_info = pickle.load(simdat)
simdat.close()


if bootstrapping == False:
    cat_loc = "/home/jmwilson/Desktop/SimBASS/cats"
    if etasBool == True:
        catdat = open(cat_loc+'/hist_dict6.p', 'r')#
    else:
        catdat = open(cat_loc+'/hist_dict6_shifted_creep.p', 'r')#
    catalogs = []
    catalogs.append(pickle.load(catdat))
    catdat.close()
else:
    cat_loc = sim_loc+'/ETAS_files/hist_dict_bootstrap.p'
    catdat = open(cat_loc, 'r')
    catalogs = pickle.load(catdat)
    catdat.close()
    



finaldat = []

lonlist = np.array(sim_info[0])
latlist = np.array(sim_info[1])

sim_grid = sim_info[2]

binlocs = np.array([[(x,y) for x in lonlist] for y in latlist])

#==============================================================================
# Ifrastructure for different bin sizes (using distance instead of rounding)
#==============================================================================
#==============================================================================
# binsize = 1 #in 10th of degree; we convert to degree after making the list
# lonrange = (min(lonlist), max(lonlist))
# latrange = (min(latlist), max(latlist))
# lonlist = np.arange(lonrange[0], lonrange[1]+1, binsize)/10.0
# latlist = np.arange(latrange[0], latrange[1]+1, binsize)/10.0
# 
# 
# binlocs = np.array([[(x,y) for x in lonlist] for y in latlist])
#==============================================================================

#------------------------------------------------------------------------------



if etasBool == True:
    binmask = tools.makePolyMask(binlocs)
else:
    binmask = tools.makeFaultMask(binlocs)

#sim_masked = np.ma.masked_array(sim_grid, mask=binmask)

#==============================================================================
# make flattened list of bins in calipoly, sort by forecast value and whether it's a hit
#==============================================================================
areasize = 0
sim_flat = []
for j in xrange(len(binmask)):
    for i in xrange(len(binmask[0])):
        if binmask[j][i] == 0:
            areasize += 1
            sim_flat.append([(binlocs[j][i][0], binlocs[j][i][1]), sim_grid[j][i]])

for cat_hist_dict in catalogs:
    
    sim_flat_scored = [line+[1] if line[0] in cat_hist_dict else line+[0] for line in sim_flat]
    
    sim_sorted = sorted(sim_flat_scored, key=lambda x: x[2], reverse = False)
    sim_sorted = sorted(sim_sorted, key=lambda x: x[1], reverse = True)
    
    
    #ROC curve creation
    Hits, Falsies = ROC_fnc(sim_sorted, cat_hist_dict, areasize)
    
    #ROC area skill test
    df = [Falsies[i+1]-Falsies[i] for i in xrange(len(Falsies)-1)]
    df = np.array(df)
    score = sum(Hits[:-1]*df)-0.5
    #print score
    
    finaldat.append({'Hits':Hits, 'Falsies':Falsies, 'score':score})

#Dumping data
if bootstrapping == False:
    f = open(sim_loc+"/ETAS_files/ROC_data_"+filesuff[etasBool], "w")
else:
    f = open(sim_loc+"/ETAS_files/ROC_data_bootstrap.p", "w")
pickle.dump(finaldat, f)
f.close()
