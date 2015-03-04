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

tzutc=pytz.timezone('UTC')

cat_names = 'date, hypo_lat, hypo_lon, magnitude, depth'
#==============================================================================
# UCERF2: lon=[-124.9322, -114.5771], lat=[31.5676, 42.2008]
#==============================================================================
eq_list = ANSStools.catfromANSS(lon=[-124.9322, -114.5771], lat=[31.5676, 42.2008], minMag=5.5, dates0=[dtm.datetime(1980,01,01, tzinfo=tzutc), None], fout='../cats/mycat_Feb-05-2015.cat')

eq_ar = np.array(eq_list)
eq_rec = np.core.records.fromarrays(eq_ar.transpose(), names=cat_names)
eq_rec.dump('../cats/cat_rec_Feb-05-2015.p')