# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 14:41:45 2015

@author: jmwilson
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import time
from geographiclib.geodesic import Geodesic as ggp

#==============================================================================
# Geo distances
#==============================================================================

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
        ysigs = vincentyNP(0.0, lams, 0.0, lamf, 0.0)
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


quakeloc = np.array([-116.1332, 33.7054])
binloc = np.array([[-120.6, 36.3]])

#==============================================================================
# print np.linalg.norm(deg2km((binloc[0]-quakeloc), (binloc[0][1]+quakeloc[1])/2.0))
# print np.linalg.norm(sphericalDistNP(quakeloc, binloc, distType='rect')[0])
# print sphericalDistNP(quakeloc, binloc, distType='rad')[0]
# print np.linalg.norm((ggp.WGS84.Inverse((quakeloc[1]+binloc[0][1])/2.0, quakeloc[0], (quakeloc[1]+binloc[0][1])/2.0, binloc[0][0])['s12']/1000.0, ggp.WGS84.Inverse(quakeloc[1], 0.0, binloc[0][1], 0.0)['s12']/1000.0))
# print ggp.WGS84.Inverse(quakeloc[1], quakeloc[0], binloc[0][1], binloc[0][0])['s12']/1000.0
#==============================================================================

print deg2km((binloc[0]-quakeloc), (binloc[0][1]+quakeloc[1])/2.0)
print sphericalDistNP(quakeloc, binloc, distType='rect')[0]
print (ggp.WGS84.Inverse((quakeloc[1]+binloc[0][1])/2.0, quakeloc[0], (quakeloc[1]+binloc[0][1])/2.0, binloc[0][0])['s12']/1000.0, ggp.WGS84.Inverse(quakeloc[1], 0.0, binloc[0][1], 0.0)['s12']/1000.0)
