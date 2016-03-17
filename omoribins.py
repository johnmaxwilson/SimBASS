# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 12:17:06 2015

@author: jmwilson

    Script to assign Omori-like earthquake halos in bins around earthquakes.
Also has function for populating a dictionary of quake bin coords
"""
#==============================================================================
# TODO: 
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import pickle
import time
import SimBASS_tools as tools

#==============================================================================
# Using nearest 10th degree to bin events
#   - binpos lat and lon are int(round(coord*10))
#==============================================================================
    
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
    norm = np.linalg.norm
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
                    #=============================================================
                    #try:
                    #    abratio = eigsort[0][0]/eigsort[1][0]
                    #    print abratio
                    #except:
                    #    abratio = 2.0
                    #=============================================================
                    abratio = 4.0
            
            transform = [[((coord[0]*cosine(strikeRad)-coord[1]*sine(strikeRad))*abratio, (coord[0]*sine(strikeRad)+coord[1]*cosine(strikeRad))) for coord in row] for row in diffgridkm] #squash and rotate
                    
            #diffgridkm = [[tools.deg2km((binpos[j][i]/10.0-quakeloc), (binpos[j][i][1]/10.0+quakeloc[1])/2.0) for i in xrange(len(binpos[0]))] for j in xrange(len(binpos))]
            
            radprime = norm(transform, axis=2)
            
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
            
            dist = Nom*(q-1.0)*(r0ssim**(q-1.0))/(r0ssim + radprime)**(q)/(2.0*math.pi*radprime) #Omori distribution
            
            dist = np.ma.fix_invalid(dist, fill_value=0).data #Remove any infinities from radii of 0
            dist = Nom/np.sum(dist)*dist #Normalize to N_omori
            total = np.sum(dist)
            
            if np.isnan(total):
                print quake
            else:
                histog += dist
    
    return histog



#==============================================================================
# Omori distribution for forecast generation
#   -coordinates are in integers (coord*10) until used by omori function
#==============================================================================
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
hist_grid = omori(sim_rec[::int(len(sim_rec)/8000)], binpos, pca=True, mc=6.0)
print 'time in minutes: ', (time.time()-t0)/60.0


rates_array = np.array([lonlist, latlist, hist_grid])
rates_array.dump(sim_loc+'/ETAS_files/hist_grid_omori6_abrat4.p')

#==============================================================================
# Just fill dictionary with ANSS events
#==============================================================================
#==============================================================================
# cat_rec = np.load('../cats/cat_rec_mc60_1980_shifted_creep.p')
# cat_dict = dict_filler(cat_rec)
#==============================================================================
