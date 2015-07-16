# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 13:58:52 2015

@author: jmwilson
"""
#==============================================================================
# TODO: 
#==============================================================================

import ANSStools
import datetime as dtm
import pytz
import numpy as np
import SimBASS_tools as tools

tzutc=pytz.timezone('UTC')

cat_names = 'date, hypo_lat, hypo_lon, magnitude, depth'
#==============================================================================
# UCERF2: lon=[-124.9322, -114.5771], lat=[31.5676, 42.2008]
# RELM:   lon=[-125.4, -113.1],       lat=[31.5, 43.0]
#==============================================================================
eq_list = ANSStools.catfromANSS(lon=[-125.4, -113.1], lat=[31.5, 43.0], minMag=6.0, dates0=[dtm.datetime(1980,01,01, tzinfo=tzutc), None], fout='../cats/mycat_mc60_1980.cat')

polyvectors = tools.getPolyVecs()

goodeqs = []
for eq in eq_list:
    if tools.isinPoly((eq[2]*10, eq[1]*10), polyvecs = polyvectors) == True:
        goodeqs.append(eq)

eq_ar = np.array(goodeqs)
eq_rec = np.core.records.fromarrays(eq_ar.transpose(), names=cat_names)
eq_rec.dump('../cats/cat_rec_mc60_1980.p')