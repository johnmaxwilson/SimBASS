# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 12:17:06 2015

@author: jmwilson

    Script to assign Omori-like earthquake halos in bins around earthquakes.
Also has function for populating a dictionary of quake bin coords
"""
#==============================================================================
# TODO: -Include ellipsoidal work from Mark's BASScast and RELM_basstests
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pickle
import time

#==============================================================================
# Using nearest 10th degree to bin events
#   - lat and lon are kept as int(coord*10)
#==============================================================================


def omori(quakes_f, binpos, abratio=2.0, q=1.5, delm=1.0, mc=5.5, b=1.0, lam=1.76, D=1.5):
    flt = float
    convert = decideg2km
    normalize = np.linalg.norm
    cosine = np.cos
    sine = np.sin
    
    lr0const = -((D+2*b)/(2.0+D))*(delm + mc) + (2.0/(2.0+D))*math.log10((2.0+D)/2.0) - lam + math.log10(q-1.0) - math.log10(2.0)
    
    histog = np.zeros((len(binpos), len(binpos[0])))
            
    for quake in quakes_f:
        quakeloc = (flt(quake['hypo_lon'])*10, flt(quake['hypo_lat'])*10)
        mag = flt(quake['magnitude'])
        
        strike=-30 
        strikeRad = math.pi/180.0*strike
        
        lr0ssim = mag*(((3+2*b)*D+2+4*b)/(4.0+2*D)) + lr0const
        r0ssim = 10**lr0ssim
        
        diff = binpos - quakeloc
        diff = [[convert(pairs) for pairs in _] for _ in diff]
                
        transform = [[((coord[0]*cosine(strikeRad)-coord[1]*sine(strikeRad))*abratio, (coord[0]*sine(strikeRad)+coord[1]*cosine(strikeRad))) for coord in _] for _ in diff]        
        
        rad = normalize(transform, axis=2)
        
        Nom = 10**(b*(mag-delm-mc))
        
#==============================================================================
#         lrupture=10**((mag/2.0) - lam)
#         Nas=10.0**((2.0/(2.0+D))*math.log10(1.+D/2.) + (D/(2.+D))*(mag - delm - mc))
#         Nprimemax = 10.0**(lNas - lrupture + math.log10(2.0))
#         r0 = Nom*(q-1.0)/Nprimemax
#         chi = (r0**(1-q))/(Nom*(q-1.0))
#         dist = 1.0/chi/(r0 + rad)**q/(2*math.pi*rad)
#         #dist = 2**(1-q)*(q-1)*10**((b+q*0.5-0.5)*mag+(1-q)*lam-b*delm-b*mc)*(0.5*10**(0.5*mag-lam)+rad)**-q
#==============================================================================
        
        dist = Nom*(q-1.0)*(r0ssim**(q-1.0))/(r0ssim + rad)**(q)/(2.0*math.pi*rad)
        
        dist = np.ma.fix_invalid(dist, fill_value=0).data
        
        histog += dist
    
    return histog


def decideg2km(degs):
    rE = 6371.0
    latkm = math.pi/180.0*rE*degs[0]/10.0
    lonkm = math.pi/180.0*rE*np.cos(math.pi/180.0*degs[0]/10.0)*degs[1]/10.0
    kms = (lonkm, latkm)
    return kms
    
    
    
def dict_filler(quakes_f, save_loc = "/home/jmwilson/Desktop/RELM/cats"):
    hist_dict_f = {}
    flt = float
    intgr = int
    rnd = round
    for quake in quakes_f:
        lon = intgr(rnd(flt(quake['hypo_lon'])*10))
        lat = intgr(rnd(flt(quake['hypo_lat'])*10))
        ind = (lon, lat)
        if ind not in hist_dict_f:
            hist_dict_f[ind] = [quake]
        else:
            hist_dict_f[ind].append(quake)
            
            
    f = open(save_loc + "/opt/hist/hist_dict_omori.p", "w")
    pickle.dump(hist_dict_f, f)
    f.close()
    
    return hist_dict_f

#==============================================================================
# Omori distribution for forecast generation
#   -everything is in integers (coord*10) until put into mesh
#==============================================================================
sim_loc = "/home/jmwilson/Desktop/RELM/sim_output"
sim_rec = np.load(sim_loc+'/sim_rec.p')

edge = 2
lonmin = -1254
lonmax = -1131
latmin = 315
latmax = 430

lonlist = np.arange(lonmin-edge, lonmax+edge+1, 1)
latlist = np.arange(latmin-edge, latmax+edge+1, 1)

mglon, mglat = np.meshgrid(lonlist/10.0, latlist/10.0)

binpos = np.array([[(x, y) for x in lonlist] for y in latlist])

hist_grid = omori(sim_rec, binpos)

rates_array = np.array([mglon, mglat, hist_grid])
rates_array.dump(sim_loc+'/opt/hist/hist_grid_omori.p')

#==============================================================================
# Just fill dictionary with ANSS events
#==============================================================================
#==============================================================================
# cat_rec = np.load('../cats/cat_rec_Mar-05-2015.p')
# cat_dict = dict_filler(cat_rec)
#==============================================================================
