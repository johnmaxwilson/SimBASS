# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 12:22:53 2015

@author: jmwilson
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
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
#==============================================================================

file_loc = "/home/jmwilson/Desktop/RELM"


sim_rec = np.load(file_loc+'/sim_output/sim_rec.p')
cat_rec = np.load('../cats/cat_rec_Feb-05-2015.p')
rates_array = np.load(file_loc+'/BASScast/kml/virtcal-hist_grid_bassomori.p')

simx = sim_rec['hypo_lon']
simy = sim_rec['hypo_lat']
realx = cat_rec['hypo_lon']
realy = cat_rec['hypo_lat']

#==============================================================================
# #Basemap plotting
#==============================================================================
fig = plt.figure(figsize=(14,10.5))

m = Basemap(projection='cyl', llcrnrlat=31.5676, urcrnrlat=42.2008, llcrnrlon=-124.9322, urcrnrlon=-114.5771, resolution='i')
m.drawcoastlines()
m.drawmapboundary(fill_color='PaleTurquoise')
#m.fillcontinents(color='lemonchiffon',lake_color='PaleTurquoise', zorder=0)
m.drawstates()
m.drawcountries()
m.drawparallels(np.arange(31,43,2),labels=[1,0,0,0])
m.drawmeridians(np.arange(-124,-115,2),labels=[0,0,0,1])

maxtick = rates_array[2].max()
#contlevels = np.arange(0, maxtick+1, maxtick/50)
m.contourf(rates_array[0], rates_array[1], rates_array[2], 50)
#m.plot(simx, simy, 'r.')
m.plot(realx, realy, 'm.')
#plt.legend()
plt.colorbar()
#plt.xlabel("Longitude")
#plt.ylabel("Latitude")
plt.show()



#==============================================================================
# #ROC plotting
#==============================================================================
#==============================================================================
# f = open(file_loc + '/BASScast/kml/virtcal-ROC_data_omori.p', 'r')
# roc_data = pickle.load(f)
# f.close()
# 
# fig = plt.figure(figsize=(12,9))
# ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
# for roc_dict in roc_data:
#     Hits = roc_dict['Hits']
#     Hits = Hits[::5]
#     Falsies = roc_dict['Falsies']
#     Falsies = Falsies[::5]
#     score = roc_dict['score']
#     ax.plot(Falsies, Hits, label = "score:%f"%(score))
#     
# plt.plot(Falsies, np.linspace(0, 1, len(Falsies)), 'k--')
# ax.legend(loc = 'upper left', bbox_to_anchor=(1, 1))
# ax.set_xlabel('False Alarms')
#==============================================================================
#==============================================================================
# ax.set_ylabel('Hits')
# plt.show()
#==============================================================================

#==============================================================================
# #sizelist=[]
# #scorelist=[]
# for roc_dict in roc_data:
#     sizelist.append(roc_dict['size']/100.0 + 0.01)
#     scorelist.append(roc_dict['score'])
# sizelist = np.array(sizelist)
# scorelist = np.array(scorelist)
# arealist = sizelist**2
# reducedscore = scorelist/arealist
# redscorelist = scorelist/sizelist
# 
# 
# plt.figure()
# plt.subplot(2,1,1)
# plt.plot(sizelist, scorelist,'b.-')
# plt.xlabel('Bin length (degrees)')
# plt.ylabel('ROC area skill score')
# plt.subplot(2,1,2)
# plt.plot(arealist, scorelist,'b.-')
# plt.xlabel('Bin area (degrees^2)')
# plt.ylabel('ROC area skill score')
#==============================================================================
#plt.show()