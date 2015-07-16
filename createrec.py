# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 12:19:28 2015

@author: jmwilson
"""
#==============================================================================
# TODO: *(?)Output appropriate data as floats and ints instead of strings
#==============================================================================

import numpy as np
import time
import cPickle
import datetime as dtm

placeholderdate = '1985/8/4'#dtm.date.today()

def collect_events(simdat_f, minmag):
    events_f = []
    anss_style = []
    eapnd = events_f.append
    anssapnd = anss_style.append
    intgr = int
    strng = str
    flt = float
    abslt = abs
    rnge = range
    length = len
    
    faultdat = open("../ALLCAL2_Faults/allcal2_1-7-11/ALLCAL2_1-7-11_Geometry.dat","r")
    faults = [line.split() for line in faultdat]
    faultdat.close()
    
    #lines = simdat_f.readlines()
    
    for line in simdat_f:#lines[8875:8890]:
        simline = line.split()
        
        if simline[0] == "200" and flt(simline[2])>=minmag and not np.isinf(flt(simline[2])) and simline[11] != "NA":# and simline[5]!= "11" and simline[5]!= "12":
            sim_mag = simline[2]
            sim_sid = simline[5]
            sim_hypodepth = simline[10]
            sim_hypodas = simline[11]
            
            # Find lat and lon of hypocenter
            record = False
            vertices = []
            vapnd = vertices.append
            
            for faultline in faults:
                #VirtCal needs a -1 after sim_sid (indexing offset)
                if faultline[:2] == ['201', strng(intgr(sim_sid))]:#-1)]:#
                    if sim_hypodas<faultline[12] or sim_hypodas>faultline[13]:
                        break
                    else:
                        record = True                    
                elif record == True and faultline[0] == '202':
                    vapnd(faultline)
                elif record == True and faultline[0] != '202':
                    break
            
            if length(vertices) > 0:
                bestdiff = 1e5
                bestlat = 0.0
                bestlon = 0.0
                hdas_flt = flt(sim_hypodas)
                for vert in vertices:
                    diff = abslt(hdas_flt - flt(vert[5]))
                    if diff < bestdiff:
                        bestdiff = diff
                        bestlat = vert[2]
                        bestlon = vert[3]
                eapnd(simline[1:] + [bestlat, bestlon])
                anssapnd([placeholderdate, flt(bestlat), flt(bestlon), flt(sim_mag), flt(sim_hypodepth)])
    return events_f, anss_style

simind = 0

event_names = 'event_id, magnitude, time, duration, sid, depth_lo, depth_hi, das_lo, das_hi, hypo_depth, hypo_das, area, mean_slip, moment, shear_before, shear_after, normal_before, normal_after, hypo_lat, hypo_lon'
# anss_elements = 'date, hypo_lat, hypo_lon, magnitude, depth'

folderext=['/allcal/version_1a', '/virtcal/z_01', '/pollitz/version_1a', '/rsqsim/version_1a']
filenames=['/ALLCAL2-30k-output[3-24-11].converted', '/ALLCAL2_eqsim_out_v2.d', '/Fred-allcal2-nofric-dip.ge.30-17300yr.txt', '/eqs.ALLCAL2_RSQSim.barall']

file_loc = "/home/jmwilson/Desktop/raw_output/allcal2"+ folderext[simind]
#!!VirtCal needs a -1 in the fault ID index!!

simdat = open(file_loc  + filenames[simind], "r")
t0= time.time()

events, ansslike = collect_events(simdat, 3.0)

t1=time.time()
print t1-t0
simdat.close()

events_array = np.array(events)
events_rec = np.core.records.fromarrays(events_array.transpose(), names=event_names)
events_rec.dump(file_loc + '/ETAS_files/sim_rec3.p')

f = open(file_loc+'/ETAS_files/ansslike3.p', 'w')
cPickle.dump(ansslike, f)
f.close()