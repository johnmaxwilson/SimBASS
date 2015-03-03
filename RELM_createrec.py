# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 12:19:28 2015

@author: jmwilson
"""
#==============================================================================
# TODO: *(?)Output appropriate data as floats and ints instead of strings
#       *Speed up
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
    
    faultdat = open("../UCERF2_Faults/ALLCAL2_1-7-11_no-creep_Geometry.dat","r")
    faults = [line.split() for line in faultdat]
    faultdat.close()
    
    for line in simdat_f:
        simline = line.split()
        if simline[0] == '200' and simline[2]>=minmag:
            # Find lat and lon of hypocenter
            record = False
            vertices = []
            vapnd = vertices.append
            hypodas = simline[11]
            
            if hypodas != "NA":
                for faultline in faults:
                    #VirtCal needs a -1 after simline[5] (indexing offset)
                    if faultline[:2] == ['201', strng(intgr(simline[5])-1)]:
                        record = True                    
                    elif record == True and faultline[0] == '202':
                        vapnd(faultline)
                    elif record == True and faultline[0] != '202':
                        break
                
                bestdiff = 1e5
                bestlat = 0.0
                bestlon = 0.0
                hdas_flt = flt(hypodas)
                for vert in vertices:
                    diff = abslt(hdas_flt - flt(vert[5]))
                    if diff < bestdiff:
                        bestdiff = diff
                        bestlat = vert[2]
                        bestlon = vert[3]
                eapnd(simline[1:] + [bestlat, bestlon])
                anssapnd([placeholderdate, flt(bestlat), flt(bestlon), flt(simline[2]), flt(simline[10])])

    return events_f, anss_style
    
event_names = 'event_id, magnitude, time, duration, sid, depth_lo, depth_hi, das_lo, das_hi, hypo_depth, hypo_das, area, mean_slip, moment, shear_before, shear_after, normal_before, normal_after, hypo_lat, hypo_lon'
# anss_elements = 'date, hypo_lat, hypo_lon, magnitude, depth'

file_loc = "/home/jmwilson/Desktop/RELM/sim_output"
#!!VirtCal needs a -1 in the fault ID index!!
simdat = open(file_loc + "/ALLCAL2_eqsim_out_v2.d","r")
t0= time.time()

events, ansslike = collect_events(simdat, 5.5)

t1=time.time()
print t1-t0
simdat.close()

events_array = np.array(events)
events_rec = np.core.records.fromarrays(events_array.transpose(), names=event_names)
events_rec.dump(file_loc + '/sim_rec.p')

f = open(file_loc+'/ansslike.p', 'w')
cPickle.dump(ansslike, f)
f.close()