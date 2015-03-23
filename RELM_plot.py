# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 12:22:53 2015

@author: jmwilson
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import RELM_inpoly as poly
import pickle


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

file_loc = "/home/jmwilson/Desktop/RELM"

sim_rec = np.load(file_loc+'/sim_output/sim_rec.p')
cat_rec = np.load('../cats/cat_rec_Mar-05-2015.p')

simx = sim_rec['hypo_lon']
simy = sim_rec['hypo_lat']
realx = cat_rec['hypo_lon']
realy = cat_rec['hypo_lat']

#==============================================================================
# #Basemap plotting
#==============================================================================
def basemap_plot(file_loc, realx, realy):
    rates_array = np.load(file_loc+'/sim_output/opt/hist/hist_grid_omori.p')
    #rates_array = np.load(file_loc+'/home/jmwilson/Desktop/RELM/BASScast/kml')
    binlocs = [[(x,y) for x in rates_array[0]] for y in rates_array[1]]
    mglon, mglat = np.meshgrid(rates_array[0]/10.0, rates_array[1]/10.0)
    rates_grid = np.ma.masked_array(rates_array[2], mask = poly.makePolyMask(binlocs))
    
    fig = plt.figure(figsize=(14,10.5))
    m = Basemap(projection='cyl', llcrnrlat=31.3, urcrnrlat=43.2, llcrnrlon=-125.6, urcrnrlon=-112.9, resolution='i')
    m.drawcoastlines()
    m.drawmapboundary(fill_color='PaleTurquoise')
    m.fillcontinents(color='lemonchiffon',lake_color='PaleTurquoise', zorder=0)
    m.drawstates()
    m.drawcountries()
    m.drawparallels(np.arange(30,44,2),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-124,-115,2),labels=[0,0,0,1])
    maxtick = np.log10(rates_grid[~rates_grid.mask].max())
    contlevels = np.linspace(0, 1.7, 50)
    #m.pcolor(mglon, mglat, np.log10(rates_grid))
    m.contourf(mglon, mglat, np.log10(rates_grid), levels=contlevels)
    #m.plot(simx, simy, 'm.')
    m.plot(realx, realy, 'm.')
    #plt.legend()
    plt.colorbar()
    #plt.xlabel("Longitude")
    #plt.ylabel("Latitude")
    
    #polyvectors = poly.getPolyVecs()
    #for vect in polyvectors:
    #     m.plot([vect[0][0]/10.0, vect[1][0]/10.0], [vect[0][1]/10.0, vect[1][1]/10.0], 'g-')
    
    plt.show()


#==============================================================================
# #ROC plotting
#==============================================================================
def ROC_plot():
    f = open(file_loc + '/sim_output/virtcal-ROC_data_omori.p', 'r')
    roc_data = pickle.load(f)
    f.close()
    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
    for roc_dict in roc_data:
        Hits = roc_dict['Hits']
        Hits = Hits[::5]
        Falsies = roc_dict['Falsies']
        Falsies = Falsies[::5]
        score = roc_dict['score']
        ax.plot(Falsies, Hits, label = "score:%f"%(score))
        
    plt.plot(Falsies, np.linspace(0, 1, len(Falsies)), 'k--')
    ax.legend(loc = 'upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel('False Alarms')
    ax.set_ylabel('Hits')
    plt.show()
    
basemap_plot(file_loc, realx, realy)