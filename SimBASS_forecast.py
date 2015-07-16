# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 11:45:29 2015

@author: jmwilson
"""

import numpy as np
import math
import time
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pickle
import cPickle
import SimBASS_tools as tools
import scipy.optimize
from scipy.stats import norm
import datetime as dtm

class ratemap(object):
    #createrec, 
    
    def __init__(self, params):
        
        return initiate(params)
    
    def initiate(self, params):
        self.params = params
        return
    
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
        
    def createrec():
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
    
    def dict_filler(quakes_f, save_loc = "/home/jmwilson/Desktop/SimBASS/cats"):
        hist_dict_f = {}
        flt = float
        intgr = int
        rnd = round
        for quake in quakes_f:
            lon = intgr(rnd(flt(quake['hypo_lon'])*10))
            lat = intgr(rnd(flt(quake['hypo_lat'])*10))
            ind = (lon, lat)
            if ind not in hist_dict_f:
                #hist_dict_f[ind] = [quake]
                hist_dict_f[ind] = 1
            else:
                #hist_dict_f[ind].append(quake)
                hist_dict_f[ind]+=1
                
                
        f = open(save_loc + "/hist_dict_shifted_creep.p", "w")
        pickle.dump(hist_dict_f, f)
        f.close()
        
        return hist_dict_f
    
    def omori(quakes_f, binpos, pca = False, q=1.5, delm=1.0, mc=3.0, b=1.0, lam=1.76, D=1.5):
        flt = float
        normalize = np.linalg.norm
        cosine = np.cos
        sine = np.sin
        
        lr0const = -((D+2*b)/(2.0+D))*(delm + mc) + (2.0/(2.0+D))*math.log10((2.0+D)/2.0) - lam + math.log10(q-1.0) - math.log10(2.0)
        
        histog = np.zeros((len(binpos), len(binpos[0])))
        
        eqlons = np.array(quakes_f['hypo_lon']).astype(float)
        eqlats = np.array(quakes_f['hypo_lat']).astype(float)
        eqcoords = np.array(zip(eqlons, eqlats))
        
        for quake in quakes_f:
            quakeloc = (flt(quake['hypo_lon']), flt(quake['hypo_lat']))
            mag = flt(quake['magnitude'])
            
            if mag >= mc and quakeloc != (0.0,0.0):
                            
                lr0ssim = mag*(((3+2*b)*D+2+4*b)/(4.0+2*D)) + lr0const
                r0ssim = 10**lr0ssim
                ruptlen = 10**(mag/2.0-lam)
                
                diffgridkm = tools.sphericalDistNP(quakeloc, binpos/10.0, distType='rect') #grid of (x,y) distances to eq in km
                #Calculate strike angle based on locations of nearby events
                fitX = []
                fitY = []
                
                distresults = zip(tools.sphericalDistNP(quakeloc, eqcoords, distType='rad'), tools.sphericalDistNP(quakeloc, eqcoords, distType='rect'))
            
                for row in distresults:
                    if row[0] <= 5.0*(ruptlen):
                        fitX += [row[1][0]]
                        fitY += [row[1][1]]
                
                if len(fitX)<=1:
                    abratio=1.0
                    strikeRad=0.0
                    
                else:
                    if pca == False:
                        abratio = 2.0
                        fitparams = tools.linfit(p=scipy.array([0.,0.]), X=fitX, Y=fitY, full_output=False)
                        strikeRad = math.pi/2.0-math.atan(fitparams[0][1]) #Strike measured clockwise from North
                    else:
                        covcalc = np.cov([fitX, fitY])
                        eig_vals, eig_vecs = np.linalg.eig(covcalc)
                        eigs = zip(eig_vals, eig_vecs.transpose())
                        eigsort = sorted(eigs, key=lambda x: x[0], reverse = True)
                        
                        strikeRad = math.pi/2.0-math.atan(eigsort[0][1][1]/eigsort[0][1][0])
                        
                        abratio = 2.0
                
                transform = [[((coord[0]*cosine(strikeRad)-coord[1]*sine(strikeRad))*abratio, (coord[0]*sine(strikeRad)+coord[1]*cosine(strikeRad))) for coord in row] for row in diffgridkm] #squash and rotate
                
                radprime = normalize(transform, axis=2)
                
                Nom = 10**(b*(mag-delm-mc))
                
                dist = Nom*(q-1.0)*(r0ssim**(q-1.0))/(r0ssim + radprime)**(q)/(2.0*math.pi*radprime) #Omori distribution
                
                dist = np.ma.fix_invalid(dist, fill_value=0).data #Remove any infinities from radii of 0
                dist = Nom/np.sum(dist)*dist #Normalize to N_omori
                total = np.sum(dist)
                
                if np.isnan(total):
                    print quake
                else:
                    histog += dist
        
        return histog
    
    def omoribinning():
        simind = 3
        folderext=['/allcal/version_1a', '/virtcal/z_01', '/pollitz/version_1a', '/rsqsim/version_1a']
        
        sim_loc = "/home/jmwilson/Desktop/raw_output/allcal2" + folderext[simind]
        sim_rec = np.load(sim_loc+'/ETAS_files/sim_rec3.p')
        
        edge = 2
        lonmin = -1254
        lonmax = -1131
        latmin = 315
        latmax = 430
        
        lonlist = np.arange(lonmin-edge, lonmax+edge+1, 1)
        latlist = np.arange(latmin-edge, latmax+edge+1, 1)
        
        binpos = np.array([[(x, y) for x in lonlist] for y in latlist])
        
        t0=time.time()
        hist_grid = omori(sim_rec, binpos, pca=True, mc=6.0)
        print 'time in minutes: ', (time.time()-t0)/60.0
                
        rates_array = np.array([lonlist, latlist, hist_grid])
        rates_array.dump(sim_loc+'/ETAS_files/hist_grid_omori6.p')        
        
    def catdictionary():
        cat_rec = np.load('../cats/cat_rec_mc60_1980_shifted_creep.p')
        cat_dict = dict_filler(cat_rec)
    
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

    def ROC():
        simind = 3
        etasBool = 1
        bootstrapping = 1
        
        filesuff = ['noETAS6.p','omori6.p']
        
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
