# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 12:17:06 2015

@author: jmwilson

    Script to assign Omori-like earthquake halos in bins around earthquakes.
Also has function for populating a dictionary of quake bin coords
"""
#==============================================================================
# TODO: -Include ellipsoidal work from Mark's BASScast
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import pickle
import time

#==============================================================================
# Using nearest 10th degree to bin events
#   - binpos lat and lon are int(round(coord*10))
#==============================================================================


def omori(quakes_f, binpos, q=1.5, delm=1.0, mc=5.5, b=1.0, lam=1.76, D=1.5):
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
        
        lr0ssim = mag*(((3+2*b)*D+2+4*b)/(4.0+2*D)) + lr0const
        r0ssim = 10**lr0ssim
        ruptlen = 10**(mag/2.0-lam)
        
        #Calculate strike angle based on locations of nearby events
        fitX = []
        fitY = []
        
        distresults = zip(sphericalDistNP(quakeloc, eqcoords, distType='rad'), sphericalDistNP(quakeloc, eqcoords, distType='rect'))
    
        for row in distresults:
            if row[0] <= 5.0*(ruptlen):
                fitX += [row[1][0]]
                fitY += [row[1][1]]
        
        if len(fitX)<=1:
            abratio=1.0
            strikeRad=0.0
            
        if len(fitX)>=2:
            abratio = 2.0
            fitparams = linfit(p=scipy.array([0.,0.]), X=fitX, Y=fitY, full_output=False)
            strikeRad = math.pi/2.0-math.atan(fitparams[0][1]) #Strike measured clockwise from North
        
        #diffgridkm = [[deg2km((binpos[j][i]/10.0-quakeloc), (binpos[j][i][1]/10.0+quakeloc[1])/2.0) for i in xrange(len(binpos[0]))] for j in xrange(len(binpos))]
        
        diffgridkm = sphericalDistNP(quakeloc, binpos/10.0, distType='rect') #grid of (x,y) distances to eq in km
        
        transform = [[((coord[0]*cosine(strikeRad)-coord[1]*sine(strikeRad))*abratio, (coord[0]*sine(strikeRad)+coord[1]*cosine(strikeRad))) for coord in row] for row in diffgridkm] #squash and rotate
        
        radprime = normalize(transform, axis=2)
        
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
        
        dist = np.ma.fix_invalid(dist, fill_value=0).data #Remove any infinities from events happening on top of bins
        dist = Nom/np.sum(dist)*dist #Normalize to N_omori
        
        histog += dist
    
    return histog

def sphericalDistNP(singleloc, otherlocs, distType='rad'):
    Rearth = 6371.0    # km
    deg2rad = math.pi/180.0
    #
    # phi = longitude, lambda = latitude
    phis = singleloc[0]*deg2rad    
    lams  = singleloc[1]*deg2rad
    phif = otherlocs[..., 0]*deg2rad
    lamf  = otherlocs[..., 1]*deg2rad
    dlam = (lamf - lams)
    dphi = (phif - phis)
    #
    #
    if distType=='rad': dsig = vincentyNP(phis, lams, phif, lamf, dphi)
    if distType=='rect':
        lamav = (lamf+lams)/2.0
        xsigs = vincentyNP(phis, lamav, phif, lamav*np.ones(phif.shape), dphi)
        xsigs[dphi<0]*=-1
        ysigs = vincentyNP(phis, lams, phis*np.ones(lamf.shape), lamf, 0.0)
        ysigs[dlam<0]*=-1
        dsig = np.concatenate((xsigs[...,np.newaxis], ysigs[...,np.newaxis]), axis=xsigs.ndim)
    #
    Dist = Rearth * dsig
    #
    return Dist

def vincentyNP(phis, lams, phif, lamf, dphi):
    #this one is supposed to be bulletproof (Vincenty formula):
    dsigma = np.arctan2( np.sqrt((np.cos(lamf)*np.sin(dphi))**2.0 + (np.cos(lams)*np.sin(lamf) - np.sin(lams)*np.cos(lamf)*np.cos(dphi))**2.0 ) , (np.sin(lams)*np.sin(lamf) + np.cos(lams)*np.cos(lamf)*np.cos(dphi))  )
    #
    return np.absolute(dsigma)

def deg2km(deltadegs, lat):
    Rearth = 6371.0
    deg2rad = math.pi/180.0
    latkm = Rearth*(deg2rad*deltadegs[1]) # Rearth * dlat
    lonkm = Rearth*np.cos(deg2rad*lat)*(deg2rad*deltadegs[0])# Rearth * cos(lat) * dlon
    kms = (lonkm, latkm)
    return kms    

def linres(p, y, x, w=None):
    if w==None:
        w=scipy.ones(len(y))
    err = w*(y - (p[0] + p[1]*x))
    return err
		
def linfit(p, X, Y, W=None, full_output=False):
    # consider replacing this with numpy.linalge.linfit() (i think... something like that) or just scipy.optimize.curve_fit()
    if W==None:
        W=scipy.ones(len(Y))
    #p=scipy.array([0.0, 0.0])
    plsq=scipy.optimize.leastsq(linres, scipy.array(p), args=(scipy.array(Y), scipy.array(X), scipy.array(W)), full_output=full_output)
    return plsq

    
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
            #hist_dict_f[ind] = [quake]
            hist_dict_f[ind] = 1
        else:
            #hist_dict_f[ind].append(quake)
            hist_dict_f[ind]+=1
            
            
    f = open(save_loc + "/opt/hist/hist_dict_omori.p", "w")
    pickle.dump(hist_dict_f, f)
    f.close()
    
    return hist_dict_f

#==============================================================================
# Omori distribution for forecast generation
#   -everything is in integers (coord*10) until put into mesh
#==============================================================================
sim_loc = "/home/jmwilson/Desktop/raw_output/allcal2/rsqsim/version_1a"
sim_rec = np.load(sim_loc+'/sim_rec.p')

edge = 2
lonmin = -1254
lonmax = -1131
latmin = 315
latmax = 430

lonlist = np.arange(lonmin-edge, lonmax+edge+1, 1)
latlist = np.arange(latmin-edge, latmax+edge+1, 1)

binpos = np.array([[(x, y) for x in lonlist] for y in latlist])

t0=time.time()
hist_grid = omori(sim_rec[:], binpos)
print 'time in minutes: ', (time.time()-t0)/60.0


rates_array = np.array([lonlist, latlist, hist_grid])
rates_array.dump(sim_loc+'/hist_grid_omori.p')

#==============================================================================
# Just fill dictionary with ANSS events
#==============================================================================
#==============================================================================
# cat_rec = np.load('../cats/cat_rec_Mar-05-2015.p')
# cat_dict = dict_filler(cat_rec)
#==============================================================================
