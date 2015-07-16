# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 12:22:53 2015

@author: jmwilson
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import SimBASS_tools as tools
import pickle
from scipy.stats import norm


#==============================================================================
# event_names = 'event_id, magnitude, time, duration, sid, depth_lo, depth_hi, das_lo, das_hi, hypo_depth, hypo_das, area, mean_slip, moment, shear_before, shear_after, normal_before, normal_after, hypo_lat, hypo_lon'
# 
# cat_names = 'date, hypo_lat, hypo_lon, magnitude, depth'
#
# rate array values = 0-mglon, 1-mglat, 2-quake_rates
#
# ROC dict keys = 'Hits', 'Falsies', 'size', 'score'
#
# UCERF2 bounding coords: lon=[-124.9322, -114.5771], lat=[31.5676, 42.2008]
# RELM:   lon=[-125.4, -113.1],       lat=[31.5, 43.0]
# polylats = [430, 430, 394, 357, 343, 329, 322, 317, 315, 319, 328, 337, 342, 377, 402, 405]
# polylons = [-1252, -1190, -1190, -1140, -1131, -1135, -1136, -1145, -1171, -1179, -1184, -1210, -1216, -1238, -1254, -1254]
#==============================================================================

#==============================================================================
# #Basemap plotting
#==============================================================================
def basemap_plot(file_loc, etasBool):
    sim_rec = np.load(file_loc+'/ETAS_files/sim_rec3.p')    
    if etasBool == True:
        cat_rec = np.load('../cats/cat_rec_mc60_1980.p') #
        rates_array = np.load(file_loc+'/ETAS_files/hist_grid_omori6.p') 
    else:
        cat_rec = np.load('../cats/cat_rec_mc60_1980_shifted_creep.p') #
        rates_array = np.load(file_loc+'/ETAS_files/hist_grid_noETAS6.p') 
    
    
    
    simx = sim_rec['hypo_lon']
    simy = sim_rec['hypo_lat']
    realx = cat_rec['hypo_lon']
    realy = cat_rec['hypo_lat']
    
    #rates_array = np.load(file_loc+'/home/jmwilson/Desktop/SimBASS/BASScast/kml')
    binlocs = np.array([[(x,y) for x in rates_array[0]] for y in rates_array[1]])
    mglon, mglat = np.meshgrid(rates_array[0]/10.0, rates_array[1]/10.0)
    #rates_grid = np.ma.masked_array(np.log10(rates_array[2]), mask = tools.makePolyMask(binlocs))
    if etasBool == True:
        binmask = tools.makePolyMask(binlocs)
        rates_grid = np.ma.masked_array(log10(rates_array[2]), mask = binmask)
    else:
        binmask = tools.makeFaultMask(binlocs)
        rates_grid = np.ma.masked_array(rates_array[2], mask = binmask)
    
    #fig = plt.figure(figsize=(14,10.5))
    m = Basemap(projection='cyl', llcrnrlat=31.3, urcrnrlat=43.2, llcrnrlon=-125.6, urcrnrlon=-112.9, resolution='i')
    m.drawcoastlines()
    m.drawmapboundary(fill_color='PaleTurquoise')
    m.fillcontinents(color='lemonchiffon',lake_color='PaleTurquoise', zorder=0)
    m.drawstates()
    m.drawcountries()
    m.drawparallels(np.arange(30,44,2), labels=[1,0,0,0])
    m.drawmeridians(np.arange(-124,-115,2), labels=[0,0,0,1])
    
    maxtick = rates_grid[~rates_grid.mask].max()
    mintick = rates_grid[~rates_grid.mask].min()
    contlevels = np.linspace(mintick, 1.3, 50)
    
    #m.contourf(mglon, mglat, rates_grid, levels=contlevels)
    m.pcolor(mglon, mglat, rates_grid)
    #m.plot(simx[::], simy[::], 'k.')
    #m.plot(realx, realy, 'm.')
    #plt.legend()
    plt.colorbar()
    #plt.xlabel("Longitude")
    #plt.ylabel("Latitude")
    
    #polyvectors = tools.getPolyVecs()
    #for vect in polyvectors:
    #     m.plot([vect[0][0]/10.0, vect[1][0]/10.0], [vect[0][1]/10.0, vect[1][1]/10.0], 'g-')
    
    #plt.show()
    
    return None

#==============================================================================
# #ROC plotting
#==============================================================================
def ROC_plot(file_loc, etasBool):
    #-----
    g = open(file_loc + '/ETAS_files/ROC_data_bootstrap2.p', 'r')#
    roc_data = pickle.load(g)
    g.close()
    scores = []
    for i in range(len(roc_data)):
        scores.append(roc_data[i]['score'])
    sortscores = sorted(scores)
    
    #plt.subplot(2,1,1)
    plt.title('ViscoSim ROC curves')
    for roc_dict in roc_data:
        Hits = roc_dict['Hits']
        Hits = Hits[::]
        Falsies = roc_dict['Falsies']
        Falsies = Falsies[::]
        plt.plot(Falsies, Hits, 'y', alpha=0.2)
    #-----
    
    if etasBool == True:
        f = open(file_loc + '/ETAS_files/ROC_data_omori6.p', 'r')
    else:
        f = open(file_loc + '/ETAS_files/ROC_data_noETAS6.p', 'r')
    
    roc_data = pickle.load(f)
    f.close()
    
    #fig = plt.figure(figsize=(12,9))
    #ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
    for roc_dict in roc_data:
        Hits = roc_dict['Hits']
        Hits = Hits[::]
        Falsies = roc_dict['Falsies']
        Falsies = Falsies[::]
        score = roc_dict['score']
        plt.plot(Falsies, Hits, 'r', label = "score: %f"%(score), linewidth=2)
    
    plt.plot(Falsies, np.linspace(0, 1, len(Falsies)), 'k--')
    plt.legend(loc = 'center right')#, bbox_to_anchor=(1, 1))
    plt.xlabel('False Alarms')
    plt.ylabel('Hits')
    
    #-----
    plt.figure()
    #plt.subplot(2,1,2)
    plt.title('ViscoSim Boostrap ROC scores histogram')
    mu, std = norm.fit(scores)
    n, bins, patches = plt.hist(scores, bins=50, normed = True, alpha=0.6)#
    ymin, ymax = plt.ylim()
    xmin, xmax = plt.xlim()
    scoresig = abs(score-mu)/std
    print scoresig
    x = np.linspace(xmin, xmax, 500)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=1.6)
    plt.plot([mu-2*std, mu-2*std], [0, 16], 'k--', linewidth=2.0, label='2 $\sigma$')
    plt.plot([mu+2*std, mu+2*std], [0, 16], 'k--', linewidth=2.0)
    plt.plot([score, score], [0, 16], 'r-', linewidth=2.0, label='%.3f $\sigma$'%scoresig)
    plt.legend()
    plt.xlabel('ROC Score')
    plt.ylabel('Normalized Count')
    #-----
    
    
    
    return None


simind = 1
etasBool = 0
file_locs = ["/home/jmwilson/Desktop/raw_output/allcal2/allcal/version_1a", "/home/jmwilson/Desktop/raw_output/allcal2/virtcal/z_01", "/home/jmwilson/Desktop/raw_output/allcal2/pollitz/version_1a", "/home/jmwilson/Desktop/raw_output/allcal2/rsqsim/version_1a"]
titles = ["AllCal", "Virtual Quake", "ViscoSim", "RSQSim"]
#for simind in range(4):
plt.close('all')
fig = plt.figure()#figsize=(9,14))
#plt.subplot(2,1,1)
#plt.title(titles[simind])
basemap_plot(file_locs[simind], etasBool)
#plt.subplot(2,1,2)    
#ROC_plot(file_locs[simind], etasBool)

plt.show()