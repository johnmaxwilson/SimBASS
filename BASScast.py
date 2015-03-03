import datetime as dtm
import pytz
import operator
import math
import random
import numpy
import scipy
import scipy.optimize as spo
import os
from PIL import Image as ipp
import multiprocessing as mpp
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mpd
import matplotlib.mpl as mpl
#
#import shapely.geometry as sgp
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from geographiclib.geodesic import Geodesic as ggp
#
import ANSStools as atp
#
days2secs = 60.*60.*24.
year2secs = 60.*60.*24.*365.
deg2km=111.12
deg2rad = 2.0*math.pi/360.
alpha_0 = 2.28
#
tzutc = pytz.timezone('UTC')
calipoly = [[-129.0, 42.75], [-121.5, 29.5], [-113.5, 29.5], [-114.0, 35.6], [-119.3, 39.4], [-119.5, 42.75], [-129.0, 42.75]]

class BASScast(object):
	# BASS/ETAS type forecast container.
	sites=[]
	quakes=[]
	conts=[]			# contour (local) objects (from sites x,y,z data).
	latrange=[]
	lonrange=[]		# lat/lon range (rectangular) of catalog
	gridsize=.1		# degrees
	midlat=None		# average latitude of the catalog, though i think we'll employ a more accurate calc. that uses both lat coords (dx = cos(l2)*x2 - cos(l1)*x1)
	mc=2.0			# minimum magnitude (at least for BASS/ETAS calculations)
	contres=2
	contfact=5
	fcdate=None
	fcdatef=None	# forecast date, float(forecast date)
	deg2km=111.12	# degrees (latitute) to km conversion
	zresolution=9	# decimal resolution of Z_ij array from which contours are calculated.
	#
	reffMag=5.0	# equivalent number of this magnitude/year earthquakes.
							# so, the contours can be read as " n m>5/km^2. nominally
							# we can integrate the contours as well if we like...
							# but from this we should be able to talk about critical rates et al.
							# n -> n0*10^(mc-rateMag.)
	mapres='i'		# basemap map resolution.
	#
	def __init__(self, incat=[], fcdate=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=2, mc=2.0, eqtheta=0, eqeps=1.0, fitfactor=5.0, contour_intervals=None, lats=None, lons=None, doBASScast=True, rtype='ssim', p_quakes=None, p_map=None):
		return self.initialize(incat=incat, fcdate=fcdate, gridsize=gridsize, contres=contres, mc=mc, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor, contour_intervals=contour_intervals, lats=lats, lons=lons, doBASScast=doBASScast, rtype=rtype, p_quakes=p_quakes, p_map=p_map)
		#
	#
	def initialize(self, incat=[], fcdate=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=2, mc=2.0, eqtheta=0.0, eqeps=1.0, fitfactor=5.0, contour_intervals=None, lats=None, lons=None, doBASScast=True, rtype='ssim', p_quakes=None, p_map=None, map_projection='cyl'):
		#
		# input parameters notes:
		# lats, lons: define the map lat, lon ranges. if None, we estimate them from the catalog. if provided, we explicitly set the map prams.
		#
		epsDefault = math.sqrt(2.0)	# some explanation: this should be, as per Kagan2002, eps=a/b. when
		#epsDefault = 2.0					# we calculate r_{effective}, we do something like
												# r_eff ^2 = (eps*x)^2 + (y/eps)^2, so the aspect ratio is:
												# abratio = eps^2, aka, eps/(1/eps). so, for abratio=2.0, eps=sqrt(2)
												# and of course 1/eps equivalently.
		if eqtheta!=None and abs(float(eqtheta))>0.0 and eqeps==None: eqeps = epsDefault
		# incat: [ [dtm, lat, lon, mag, (depth?)], [], ... ] (consistent with ypp.eqcatalog() format).
		# note: catalog dates/times must be in seconds (or convert all the scaling bits to days-- yuck).
		self.mc = mc
		self.contres = contres
		self.gridsize = gridsize
		#
		# so simulator catalogs will probably produce dates that are out of datetime range. hijack this and pivot towards a more
		# float-centric date management.
		#
		self.fcdate = fcdate
		if isinstance(fcdate, dtm.datetime):
			# and note that if we pass fcdate as a float, we need to convert from days to secs??? sort of. times
			# get converted from days to seconds in the location() object initilizations. let's do that here as well;
			# pass in a "days" and we'll convert here.
			fcdatef = mpd.date2num(fcdate)*days2secs
		else:
			fcdatef = fcdate*days2secs
												# and since days2secs ~ ~ 10**7.5, we can probably guess if seconds have been passed
												# (if we want to). a more likely error would be to pass years, not days, and that's
												# harder to trap.
		self.fcdatef = fcdatef
			#
		#
		self.fitfactor = fitfactor
		#if rtype not in ('ssim', 'omorisat', 'satpl'): rtype='ssim'
		self.rtype = rtype
		self.p_quakes = p_quakes
		self.p_map = p_map
		self.map_projection = map_projection
		#
		self.quakes = []
		while len(self.quakes)>0: self.quakes.pop()
		#
		self.contour_intervals=contour_intervals
		#
		# read a catalog
		if type(incat)==type('astring'):
			# it's a filename. assume it's formatted for a catalog object...
			# ... and we're not sure what to do with this just yet. for now, require a list
			# [ [dtm, lat, lon, mag, (depth?)], [], ... ] (consistent with ypp.eqcatalog() format.
			#
			# load catalog into [quakes] and determine X, Y range from catalog.
			print "we need a list-catalog. files are not supported at this time."
			#
			return None
		#
		# make sure, as best we can, that incat is a list type object (we might start seeing recarrays).
		# recarrays would be better, but theyr'e not what whe have, and there are various operations that 
		# will need to be recoded -- namely indexing and maybe some slicing.
		#
		if hasattr(incat, 'tolist'): incat=incat.tolist()
		#
		# trap some weird possible errors:
		if len(incat[0])==1:
			# this can happen in recarray conversions...
			incat = [rw[0] for rw in incat]	# row got nested in a list or tuple.
		
		#
		#print "incat type: ", type(incat)
		self.catalog=incat	# dt, lat, lon, mag
		#
		# permit an empty catalog:
		if len(self.catalog)==0:
			# nothing to do, but we can still use the member functions.
			return None
			
		# initialize the lat/lon range arrays:
		# have lat/lon range been provided? if so then set.
		if lats==None or lats==[] or len(lats)<2: lats=None
		if lons==None or lons==[] or len(lons)<2: lons=None	# the "None" bit looks stupid, but it prevents an error on len(None)
		self.latrange = lats
		self.lonrange = lons
		#... otherwise, we'll determine lat/lon range from the catalog (also see below).
		if lats==None: self.latrange = [incat[0][1], incat[0][1]]
		if lons==None: self.lonrange = [incat[0][2], incat[0][2]]
		#
		#print "initializing. latrange: ", self.latrange
		#		
		# aftershock convolution...
		# so, we want an indexed catalog (sorted by lats, then lons or vice versa) so we can quickly fit the events around each quake
		# also, of course, these need to be linear distances, noting that the exact distances around the sub-catalog are not terribly
		# important, as the idea is just to get a reasonable linear fit for direction and epsilon (width).
		# then, we need a catalog of loc, sig, theta (though we can store these vals. in quakes[] and write a script to retrieve
		# them in various forms.
		print "eqtheta, eqeps: ", eqtheta, eqeps
		if eqtheta==None or eqeps==None:
			# in this case, do a dynamic transform.
			#xy=map(operator.itemgetter(1,2), incat)	# so this will be lat, lon
			#xy.sort()	# and it's now sorted by lat, then lon.
			#incat.sort(key=operator.itemgetter(1))	# incat is now sorted by lat
			#print "incat[0]: ", incat[0]
			incat.sort(key = lambda x: (x[1], x[2]))
			#
			# should this be sorted by lat, lon?:
			# i think:
			# incat.sort(key = lambda x: x[1], x[2])
		#
		# add earthquakes and other catalog related bits.
		catindex=0
		for rw in incat:
			# ... and to accomplish mpp, the whole contents of this loop needs to be wrapped into a single function.
			# it may also be necessary (or at least desirable) to mitigate the dependency on class scope (aka,
			# modularize the function -- pass all values directly, rather than pull them from "self"). 1) classes
			# seem to not pickel into mpp properly, 2) if possible, we want to avoid having to pickle the whole
			# class anyway.
			#
			# we'll want to introduce the elliptical convolution with fault map or local seismicity. for now, fudge it.
			# each rw: [dtm, lat, lon, mag, ...]
			# note: loc[] for earthquake is [lon, lat]; we receive incat as lat,lon
			self.quakes+=[earthquake(mag=rw[3], loc=[rw[2], rw[1]], evtime=rw[0], mc=mc, rtype=self.rtype, p=p_quakes, p_map=p_map)]
			#
			thiseq=self.quakes[-1]
			#
			# start with "default global" type rupture angle (strike) and epsilon[elliptical eccentricity]
			# provided as parameters):
			thiseq.eqtheta=eqtheta
			thiseq.eqeps=eqeps
			#
			# if available...
			# set earthquake orientation (if available). keep an eye on this feature; we might be, in some cases,
			# passing long catalog entries where cols [4,5] are NOT theta,epsilon.
			# this is a problem in many existing scripts since rw[4] is often depth. this functionality was put in place for one of the
			# scenario exercises, right? maybe we can (in the short term) correct it there and expand this to at least rw[5], rw[6].
			# eventually, replace list-rows with dictionary-rows... or a structured array.
			# note also that a structured array (recarray, or "record array" -- all numpy)
			# can (i think) be shared directly in mpp processing.
			#
			# for now, use len()>6, etc.
			if len(rw)>=6:
				if rw[5] > 0.: thiseq.eqeps = rw[5]
			if len(rw)>=6:
				if rw[6]>=(-360.) and rw[6]<=360.:
					thiseq.eqtheta = rw[6] # some simple validation.
			#
			# and after all that, did we get an epsilon?
			if thiseq.eqeps==None: thiseq.eqeps=epsDefault	# so, at this point eqeps is either what was given or default.
			#
			# rupture length, to modify catalog dimensions (we set the grid to the lat/lon from the catalog + rupture-len)
			rl=thiseq.getruptlen()
			#
			# now, get "local" quakes and fit to a lat/lon line
			# the catalog is indexed so this should go quickly and we can use a "while"
			#if eqtheta==None or eqeps==None:
			
			if thiseq.eqtheta==None:
				#thisR = self.dispFrom(loc=self.loc, coords='rad')[0]
				tmpX=[]
				tmpY=[]	# temporary lists for fitting
				#tmpW=[]
				tmpW=None
				thisRLat = 0.0
				#theseLats=[thiseq.loc[1]-(rl/2.)*, thiseq.loc[1]+rl/2.]
				#theseLons=[thiseq.loc[0]-rl/2., thiseq.loc[0]+rl/2.]
				i=0
				iup=i
				idn=i
				rfactor=rl*float(fitfactor)
				while (thisRLat<rfactor and (iup>=0 or idn<len(incat))):
					# be sure these are returning the same units. we have two schools, one running km, another m -- and arguably another 
					# in lat/lon degrees or even radians.
					# now, what sort of transform do we use for the fitting?
					#
					# we'll go one before and one after the event in question
					iup = catindex-i
					idn = catindex+i
					#if catindex-i>0:	# we haven't retreated past the begining of the catalog...
					for newi in [idn, iup]:
						if newi>=0 and newi<len(incat):
							testrw = incat[newi]
							thisR = thiseq.sphericalDist([testrw[2], testrw[1]])
							thisRlat = thiseq.sphericalDist([thiseq.loc[0], testrw[1]])	# "latitude distance"
						# not necessary:
						#else:
						#	thisR=None
						#	break
								
						if thisR!=None and thisR<=rfactor:
							# this earthquake is inside our test-circle, so we'll try again later...
							#getfitters=True
							# and of course, these require some sort of linear transformation...
							# let's use a simple trig. transform from the latitute of the source event.
							#print "***", thiseq.loc[0], testrw[2]
							tmpX+=[(thiseq.loc[0] - testrw[2])*self.deg2km*math.cos(thiseq.loc[1]*deg2rad)]	# R * Delta-Lon * latfact
							tmpY+=[(thiseq.loc[1] - testrw[1])*self.deg2km]	# R*\Delta-Lat
							# tmpW+=[thiseq.mag]
					i+=1
					
					# all of this if-logic can be cleaned up a bit, but at least it's correct (i think) for now.
					#
					#thiseq.eqtheta = 0.0
					#thiseq.eqeps=eqeps
					#if thiseq.eqeps == None:
					#	thiseq.eqeps=epsDefault	# if there are no data to fit for direction, make it a circle
					#							# (unless otherwise specified from input prams).
					#
				
				# note the first entry in tmpX,Y is its distance from itself. this was (originally) by design, and is
				# easy enough to work around...	
				if len(tmpX)<=2:
					if eqeps==None:
						thiseq.eqeps=1.0
						thiseq.eqtheta=0.0
				#
				#print len(tmpX), thiseq.eqeps, eqeps, thiseq.eqtheta, thiseq.mag, thiseq.loc[0], thiseq.loc[1]
				#print tmpX
				#print tmpY
				#
				if len(tmpX)>=3:
					#thiseq.epeqs=eqeps
					if eqtheta==None:
						# if eqtheta is defined, don't redefine it.
						#if tmpW==None: tmpW=scipy.ones(len(tmpY))
						fitprams=self.linfit(p=scipy.array([0.,0.]), X=scipy.array(tmpX), Y=scipy.array(tmpY), W=tmpW, full_output=False)
						#print len(self.errs), sum(self.errs), sum(self.errs)/float(len(self.errs))
						#print "fitprams: " , fitprams, len(tmpX), len(tmpY)
						#
						thiseq.eqtheta = math.atan(fitprams[0][1])/deg2rad
					#if thiseq.eqeps==None:
					#	thiseq.eqeps=epsDefault
			#	
			dlon=rl/(self.deg2km*math.cos(rw[1]*2.0*math.pi/360.0))		# noting the deg -> rad and the km -> deg. conversions.
			dlat=rl/self.deg2km
			# adjust the forecast extents according to catalog events and their rupture lens.
			if lats==None:
				if self.latrange[0]>(rw[1]-dlat): self.latrange[0] = (rw[1]-dlat)		# expand min-lat?
				if self.latrange[1]<(rw[1]+dlat): self.latrange[1] = (rw[1]+dlat)		# expand max-lat?
			if lons==None:
				if self.lonrange[0]>(rw[2]-dlon): self.lonrange[0] = (rw[2]-dlon)
				if self.lonrange[1]<(rw[2]+dlon): self.lonrange[1] = (rw[2]+dlon)
			#
			catindex+=1
			#
		#####
		#
		self.midlat = .5*(sum(self.latrange))
		#
		# for now, let's use gridsize x gridsize (in degrees) size elements (rather than project into a proper rectangular coord-system,
		# though this is potentially not terribly difficult for this application since the grid is not strictly enforced (aka, we
		# have a set of center-points for which we calc. z-values based on distance between [(x1,y1), (x2, y2)]. the potential difficulty arises
		# from the fact that the array from which contours are calculated must be square... so maybe it is difficult.
		#
		# set up forecast sites:
		x0 = self.lonrange[0]-self.lonrange[0]%gridsize
		xmax = gridsize + self.lonrange[1]-self.lonrange[1]%gridsize
		self.nLon = int(math.ceil((xmax-x0)/gridsize))
		#
		y0 = self.latrange[0]-self.latrange[0]%gridsize
		ymax = gridsize + self.latrange[1]-self.latrange[1]%gridsize
		self.nLat = int(math.ceil((ymax-y0)/gridsize))
		#
		#print "during initialization, lonrange: ", self.lonrange, " and latrange: ", self.latrange
		y=y0
		#
		self.sites=[]
		for i in xrange(self.nLat*self.nLon):
			self.sites+=[forecastsite(loc=[x0+gridsize*(i%self.nLon), y0+gridsize*(i/int(self.nLon))], dxdy=[gridsize, gridsize], evtime=self.fcdatef, mc=mc)]
		#		
		#print "lengths: %d, %d, %d" % (self.nLon, self.nLat, len(self.sites))
		#
		if doBASScast==True:
			# and now, this object should br ready to calculate forecasts...
			#bc=self.calcBASScast()
			#bc=self.plotContours()
			print "calculating BASScast with rtype=%s" % self.rtype
			bcast=self.calcBASScast()
			#self.conts = self.getContourSet(X_i=self.X_i, Y_i=self.Y_i, Z_ij=self.Z2d.round(self.zresolution), contres=self.contfact*self.contres)
			#self.conts = self.getContourSet(X_i=bcast[0], Y_i=bcast[1], Z_ij=bcast[2].round(4), contres=self.contfact*self.contres)
			self.conts = self.getContourSet(X_i=bcast[0], Y_i=bcast[1], Z_ij=bcast[2], contres=self.contfact*self.contres)
		#
		# initialize a basemap:
		cm_cntr = [self.lonrange[0] + .5*(self.lonrange[1]-self.lonrange[0]), self.latrange[0] + .5*(self.latrange[1]-self.latrange[0])]
		self.cm = Basemap(llcrnrlon=self.lonrange[0], llcrnrlat=self.latrange[0], urcrnrlon=self.lonrange[1], urcrnrlat=self.latrange[1], resolution=self.mapres, projection=self.map_projection, lon_0=cm_cntr[0], lat_0=cm_cntr[1])
		#
		return None
	
	def linres(self, p, y, x, w=None):
		if w==None:
			w=scipy.ones(len(y))
		err = w*(y - (p[0] + p[1]*x))
		self.errs = err
		return err
		
	def linfit(self, p, X, Y, W=None, full_output=False):
		# consider replacing this with numpy.linalge.linfit() (i think... something like that) or just scipy.optimize.curve_fit()
		# do we use self.errs for anything? (doesn't look like it).
		if W==None:
			W=scipy.ones(len(Y))
		#p=scipy.array([0.0, 0.0])
		plsq=spo.leastsq(self.linres, scipy.array(p), args=(scipy.array(Y), scipy.array(X), scipy.array(W)), full_output=full_output)
		return plsq
	#
	def resetzvals(self):
		# set all z to 0 or None.
		for i in xrange(len(self.sites)):
			self.sites[i].z=None
	#
	def calcBASScast(self, grid=None, quakes=None, nx=None, ny=None, gridsize=None):
		if grid==None: grid = self.sites
		if quakes==None: quakes = self.quakes
		if nx==None: nx = self.nLon
		if ny==None: ny = self.nLat
		if gridsize==None: gridsize=self.gridsize
		#
		#print quakes
		fnow=mpd.date2num(dtm.datetime.now(pytz.timezone('UTC')))
		#for q in quakes:
		#	#print q.eventDtm, q.eventftime, q.mag, q.mc, q.timesince(734974.4474623842)
		#	print q.eventDtm, q.eventftime, q.mag, q.mc, q.timesince(fnow)
		#
		X=[]
		Y=[]
		Z=[]	
		for site in grid:
			site.setz(equakes=quakes)
			# note, z values (eq counts) have been summed at this point
			#
			X+=[site.loc[0]]
			Y+=[site.loc[1]]
			if site.z!=0: Z+=[math.log10(site.z)]
			if site.z==0: Z+=[None]
			#
		#	
		#Z=map(operator.itemgetter(2), cdata)
		Z2d=scipy.array(Z)
		#
		#Z2d.shape=(len(Y), len(X))
		Z2d.shape=(ny, nx)		# aka, (ny rows, of nx elements) such that ny*nx=N
		#X1, Y1=numpy.meshgrid(X, Y)
		#
		X_i=numpy.array(map(float, range(nx)))
		X_i*=gridsize
		X_i+=grid[0].loc[0]
		#
		Y_i=numpy.array(map(float, range(ny)))
		Y_i*=gridsize
		Y_i+=grid[0].loc[1]
		#
		self.X_i = X_i
		self.Y_i = Y_i
		self.Z2d=Z2d
		#
		return [X_i, Y_i, Z2d]
	
	def calcBASScastinPoly(self, grid=None, quakes=None, nx=None, ny=None, gridsize=None):
		if grid==None: grid = self.sites
		if quakes==None: quakes = self.quakes
		if nx==None: nx = self.nLon
		if ny==None: ny = self.nLat
		if gridsize==None: gridsize=self.gridsize
		#
		#print quakes
		fnow=mpd.date2num(dtm.datetime.now(pytz.timezone('UTC')))
		#
		#X, Y, Y = [], [], []
		xyz=[]
		#
		polyvecs = self.getPolyVecs(calipoly)
		for site in grid:
			x=site.loc[0]
			y=site.loc[1]
			if self.isinPoly(pt=[x,y], polyvecs=polyvecs)==False:
				# note: we might want to just skip the z-cacl and return a full grid. this would facilitate easy mpl contouring.
				# for now, just move on and return a compressed data set. it's not too hard to fill in the missing bits on the return side.
				continue
			z=site.getz(equakes=quakes)
			# note, z values (eq counts) have been summed at this point
			#
			#X+=[site.loc[0]]
			#Y+=[site.loc[1]]
			if z!=0.0: lz=math.log10(z)
			if z==0.0: lz=None
			#
			xyz+=[[site.loc[0], site.loc[1], lz]]
			#
		#
		return xyz
	#
	def isinPoly(self, pt, polyvecs=None):
		# be careful; this needs to be a polyvec, not just a polygon...
		if polyvecs==None: polyvecs=self.getPolyVecs(calipoly)
		if type(polyvecs[0][0])==type(0) or type(polyvecs[0][0]) == type(1.0):
			polyvecs=self.getPolyVecs(polyvecs)
		#
		# now, how many "left of" (or "right of") vectors (where y1 < y_pt < y2)
		x0=pt[0]
		y0=pt[1]
		nlefts=0
		isin=False
		#
		for vec in polyvecs:
			x1, x2 = vec[0][0], vec[1][0]
			y1, y2 = vec[0][1], vec[1][1]
			#if (y0>y1 and y0<y2) or (y0<y1 and y0>y2):
			if (y0<min(y1, y2)) or (y0>=max(y1, y2)): continue
			#
			# we're in the y-domain of this vector.
			if x0<x1 and x0<x2:
				nlefts+=1
			else:
				# might still be left...
				b=(y2-y1)/(x2-x1)
				y_line = y1 + b*(x0-x1)
				if y0<=y_line: nlefts+=1
				#
			#
		#
		if nlefts%2==0: isin=False
		if nlefts%2==1: isin=True
		#
		return isin
	
	def getPolyVecs(self, poly):
		polyvecs=[]
		for i in xrange(len(poly)):
			polyvecs+=[[poly[i], poly[(i+1)%(len(poly))]]]
			if polyvecs[-1][0]==polyvecs[-1][1]: polyvecs.pop()
		return polyvecs

	#
	def gridplot(self, ncontours=25):
		plt.figure(0)
		plt.clf()
		minval=self.Z2d.min()
		maxval=self.Z2d.max()
		zrange=maxval-minval
		#
		hexColorMin='0x000000'
		hexColorMax='0xFFFFFF'
		colorRange=int('0xFFFFFF', 0)
		#
		#print "gridlens: %d, %d" % (len(self.Z2d), len(self.Z2d[0]))
		#		
		#dz=(maxval-minval)/float(ncontours)
		
		for y in xrange(len(self.Z2d)):
			for x in xrange(len(self.Z2d[0])):
				colorint=int(colorRange*(self.Z2d[y,x]-minval)/zrange)
				#colorhex=self.fillHexString(hex(colorint), 6)
				colorhex=fillHexString(hex(colorint), 6)
				colorhexstr='#' + colorhex.split('x')[1]
				plt.plot([x], [y], '.', color=colorhexstr, ms=10)
				
	
	'''
	def plotContours(self, grid=None, quakes=None, nx=None, ny=None, gridsize=None):
		XYZ=self.calcBASScast(grid=grid, quakes=quakes, nx=nx, ny=ny, gridsize=gridsize)
		return self.BASScastContours(XYZ[0], XYZ[1], XYZ[2])
	
	def mapContours(self, grid=None, quakes=None, nx=None, ny=None, gridsize=None):
		XYZ=self.calcBASScast(grid=grid, quakes=quakes, nx=nx, ny=ny, gridsize=gridsize)
		return self.BASScastContourMap(XYZ[0], XYZ[1], XYZ[2])
		
	def BASScastContours(self, X_i, Y_i, Z2d, fignum=0):
		#
		if fignum!=None: plt.figure(fignum)
		#plt.clf()
		cnts = self.getContourSet(X_i, Y_i, Z2d, self.contfact*self.contres)
		#
		self.conts = cnts
		#
		return None
	'''
		
	def BASScastContourMap(self, X_i=None, Y_i=None, Z2d=None, fignum=1, maxNquakes=250.0, alpha=.75):
		#
		if X_i==None: X_i = self.X_i
		if Y_i==None: Y_i = self.Y_i
		if Z2d==None: Z2d = self.Z2d
		#
		plt.figure(fignum)
		plt.clf()
		plt.ion()
		#
		cntr = [self.lonrange[0] + .5*(self.lonrange[1]-self.lonrange[0]), self.latrange[0] + .5*(self.latrange[1]-self.latrange[0])]
		try:
			cm = self.cm
		except:
			cm = Basemap(llcrnrlon=self.lonrange[0], llcrnrlat=self.latrange[0], urcrnrlon=self.lonrange[1], urcrnrlat=self.latrange[1], resolution=self.mapres, projection=self.map_projection, lon_0=cntr[0], lat_0=cntr[1])
			self.cm=cm
		cm.drawcoastlines(color='gray', zorder=1)
		cm.drawcountries(color='gray', zorder=1)
		cm.drawstates(color='gray', zorder=1)
		cm.drawrivers(color='gray', zorder=1)
		cm.fillcontinents(color='beige', zorder=0)
		# drawlsmask(land_color='0.8', ocean_color='w', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution='l', grid=5, **kwargs)
		#cm.drawlsmask(land_color='0.8', ocean_color='c', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution=self.mapres, grid=5)


		print "lat, lon ranges: ", self.latrange, self.lonrange
		cm.drawmeridians(range(int(self.lonrange[0]), int(self.lonrange[1])), color='k', labels=[0,0,1,1])
		cm.drawparallels(range(int(self.latrange[0]), int(self.latrange[1])), color='k', labels=[1, 1, 0, 0])
		#
		# get X,Y for contours:
		#
		X,Y=cm(*numpy.meshgrid(X_i, Y_i))
		#X,Y=cm(X_i, Y_i)
		#		
		cnts = self.getContourSet(X, Y, Z2d, contres=self.contfact*self.contres)
		plt.colorbar()
		plt.spectral()
		#
		#maxNquakes=250.0	# max num. equakes to plot (approximately)
		Nquakes = float(len(self.quakes))
		#
		if maxNquakes==0: 
			mthresh=27.
		else:
			mthresh=math.log10(Nquakes/maxNquakes) + self.mc	# we can take  out the "if" bc for N<Nmax, mthresh<mc
		#
		# now, plot the earthquakes:
		for q in self.quakes:
			if q.mag<mthresh: continue
			qx, qy, qm = q.loc[0], q.loc[1], q.mag
			qxp, qyp = cm(qx, qy)
			cm.plot([qxp], [qyp], 'ro', ms=qm, alpha=alpha)
		#
		#self.conts = cnts
		#
		return None
	
	def getContourSet(self, X_i=None, Y_i=None, Z_ij=None, contres=None, zorder=4, alpha=.3):
		#
		if X_i==None: X_i=self.X_i
		if Y_i==None: Y_i=self.Y_i
		if Z_ij==None: Z_ij=self.Z2d
		if contres==None: contres=self.contfact*self.contres
		#
		# X_i, Y_i are x,y coordinates (aka, like range(100), range(200) ), Z_ij is a 2D array that is actually contoured.
		# contres is the resolution of the contoursl; numConts -> 5*contres
		#
		resolution = contres
		#LevelsNumber = self.contfact * resolution		# note: this is probably an accidents, since we multiply by contfact twice.
		LevelsNumber = resolution
		warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#
		if self.contour_intervals!=None:
			try:
				LevelsNumber = self.contour_intervals	# note this is intended to be a list of contour values.
			except:
				# do nothing...
				LevelsNumber = resolution
		#
		# cnts=plt.contourf(gridsize*numpy.array(range(nelements)), gridsize*numpy.array(range(nelements)), Z2d,10)
		# retrieve the collections() object from the contourf() function which returns a matplotlib.contour.ContourSet
		#
		#cs = plt.contourf(X_i, Y_i, Z_ij, LevelsNumber, cm=plt.spectral()).collections
		
		cs = plt.contourf(X_i, Y_i, Z_ij, LevelsNumber, cm=plt.spectral(), alpha=alpha, zorder=zorder)
		#cs = plt.contourf(X_i, Y_i, Z_ij, LevelsNumber, alpha=.3)
	
		return cs

	def arrayFromSet(self, cs=None, alpha='9d'):
		# get an array (list) from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
		# format:
		# [ [z0, [[x,y], [x,y],...], [z1, [paths1]], 
		#
		if cs==None: cs=self.conts
		if alpha==None: alpha='9d'
		if type(alpha)==type(1) or type(alpha)==type(1.0): alpha='%02x' % alpha	# if alpha is a number, convert to string-hex
		#																									# (aka, '%02x' % 220 --> 'dc')
		#
		levels=cs.levels	# lower bound of contours
		layers=cs.layers	# upper bound of contours (at least for contours <0)
		dz= layers[0] - levels[0]
		collects=cs.collections
		carray = []
		for i in xrange(0, len(collects)):
			bgra_array = 255*collects[i].get_facecolor()[0]
			#strclr = '7d%02x%02x%02x' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
			strclr = '%s%02x%02x%02x' % (alpha, bgra_array[2] , bgra_array[1] , bgra_array[0] )
			#
			for trace in collects[i].get_paths():
				#carray+=[[levels[i], layers[i], '%s' % strclr, []]]
				# objContour(low=None, high=None, alpha=None, R=None, G=None, B=None, verts=None, RGB=None)
				carray += [objContour(low=cs.levels[i], high=cs.layers[i], RGB=strclr, verts=[])]
			
				for lng, lat in trace.vertices:
					#carray[-1][-1]+=[[lng, lat]]
					carray[-1].verts+=[[lng, lat]]
					#
				#
			#
		#
		#
		return carray
	#
	def fixedContsToFile(self, cset=None, file_out='contour_shapes.txt', top=1.0, bottom=0.0):
		contours_string = self.fixedContsToStr(cset=cset, top=top, bottom=bottom)
		#
		fout = open(file_out, 'w')
		fout.write(contours_string)
		fout.close()
		#
		return contours_string
	#	
	def fixedContsToStr(self, cset=None, top=1.0, bottom=0.0):
		if cset==None: cset=self.conts
		levels = cset.levels
		fixed_conts = self.contsToPlotLists(cset=cset)
		#				# contours are then like [ level0, level1, ...]
		#				# level: [poly0, poly1, poly2...]
		#				# poly: [X, Y]
		#
		n_levels = len(fixed_conts)-1	# and, of course, not exactly number of levels, but the max index.
		#top      = max(n_levels - int(math.ceil(top*float(n_levels))), 0)
		#bottom   = min(int(math.ceil(bottom*float(n_levels))), n_levels)	# (and just in case we accidentally calc <0 or >len)
		#top      = int(math.ceil(top*float(n_levels)))
		#bottom   = max(n_levels-int(math.ceil(bottom*float(n_levels))), 0)	# (and just in case we accidentally calc <0 or >len)
		top      = int(math.ceil(top*float(n_levels)))
		bottom   = int(math.floor(bottom*float(n_levels)))	# (and just in case we accidentally calc <0 or >len)
		#
		#
		outstr='# contour coordinates\n'
		outstr+='# this functionality is adapted somewhat hastily from some diagnostic functions, so it is a bit rough.\n'
		outstr+='# #!contour indicates a contour level (change). #!trace indicates a trace, or contour shape (borrowed \n'
		outstr+='# from the matplotlib vernacular). matplotlib produces contour shapes in a funny way, so there is a process \n'
		outstr+='# to fix this. in the process, some of the shapes are grouped into sets where there were inner/outer sets and \n'
		outstr+='# spatially independent shapes. this accounts for the notation "trace n,m", in which m indicates shapes within \n'
		outstr+='# a group of shapes n. all shapes within #!contour blocks belong to the same contour; aka, the n,m indices could \n'
		outstr+='# be flattened (and probably will in later iterations of this code.\n\n'
		#
		# levels:
		outstr+='#!contour levels:\t'
		for x in levels: outstr+='%f\t' % x
		outstr=outstr[:-1]+'\n'
		#
		#for i in xrange(len(cs)):
		for i in xrange(bottom, top):
			#
			outstr+='#!contour\t%d\n' % i
			#itrace=0
			for itraces, traces in enumerate(fixed_conts[i]):
			#for trace in cs[i]:
				for itrace, trace in enumerate(traces):
					#print trace
					#print len(trace), len(trace[0])
					#continue
					#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
					outstr+='#!trace\t%d,%d\n' % (itraces, itrace)
					#
					for i_ll in xrange(len(trace[0])):
					#for lng,lat in trace.vertices:
						lng = trace[0][i_ll]
						lat = trace[1][i_ll]
					#for lng, lat in zip(*trace):
						#kmlstr+='%s,%s,0\n' % (round(lng,7),round(lat,7))
						outstr+='%f\t%f,0\n' % (lng, lat)    # and z=0
						#tmpLL+=[[lng, lat]]
					#if tmpLL[0]!=tmpLL[-1]:
					#	print "completing a polygon, %d." % fixedpolys
					#	fixedpolys+=1
					#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])
		#
		return outstr


	#
	def contsToStr(self, cset=None, top=1.0, bottom=0.0):
		'''
		# top/bottom: max/min contours to be used. "top" is the highest level.
		'''
		# contsToPlotLists(self, cset=None)
		if cset==None: cset=self.conts
		cs=cset.collections
		levels = cset.levels
		#
		# what's our range?
		if top==None: top=1.0
		if bottom==None: bottom=0.0
		#
		# if top, bottom are >1, assume they're the number/index of the contour range being selected.
		# otherwise, it's a percentage.
		# first contour is the highest value.
		#
		n_levels = len(levels)-1
		#top      = max(n_levels - int(math.ceil(top*float(n_levels))), 0)
		#bottom   = min(int(math.ceil(bottom*float(n_levels))), n_levels)	# (and just in case we accidentally calc <0 or >len)
		#top      = int(math.ceil(top*float(n_levels)))
		#bottom   = max(n_levels-int(math.ceil(bottom*float(n_levels))), 0)	# (and just in case we accidentally calc <0 or >len)
		top      = int(math.ceil(top*float(n_levels)))
		bottom   = int(math.floor(bottom*float(n_levels)))	# (and just in case we accidentally calc <0 or >len)

		#
		outstr='# contour coordinates\n'
		#
		# levels:
		outstr+='#!contour levels:\t'
		for x in levels: outstr+='%f\t' % x
		outstr=outstr[:-1]+'\n'
		#
		#for i in xrange(len(cs)):
		
		for i in xrange(bottom, top):
			outstr+='#contour\t%d\n' % i
			itrace=0
			for itrace, trace in enumerate(cs[i].get_paths()):
			#for trace in cs[i]:
				#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				outstr+='#trace\t%d\n' % itrace
				for lng,lat in trace.vertices:
				#for lng, lat in zip(*trace):
					#kmlstr+='%s,%s,0\n' % (round(lng,7),round(lat,7))
					outstr+='%f\t%f,0\n' % (lng, lat)    # and z=0
					#tmpLL+=[[lng, lat]]
				#if tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
				#	fixedpolys+=1
				#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])
		#
		return outstr

	#self.plotcolor=plotcolor
	
	def Z2str(self):
		# this is a fast but brutally unindexed and ambiguous format -- just a raw copy of the grid.
		outstr=''
		for rw in self.Z2d:
			for elem in rw:
				outstr+='%f\t' % elem
			outstr=outstr[:-1] + '\n'
		#
		return outstr
	#
	def xyztofile(self, outfile='BASSxyz.dat', minlogz=None):
		# this writes a format like lon \t lat \t z \n, ... to a file. do we need a version to write to a list, or can we avoid risking the memory?
		fout=open(outfile,'w')
		fout.write('#BASScast xyz:\n#lon, lat, log_10(z)\n#!dx=%f\tdy=%f\tlog(z)_min=%s, m_c=%f\n' % (self.sites[0].dxdy[0], self.sites[0].dxdy[1], str(minlogz), self.mc))
		fout.write('# this is a grid representation of Yoder et al. (2014)/NASA E-DECIDER ETAS. z ~ number of earthquakes/km**2 (noting that this is for m>m_c\n')
		for s in self.sites:
			lz=math.log10(s.z)
			if minlogz!=None and lz<minlogz: continue
			fout.write('%f\t%f\t%f\n' % (s.loc[0], s.loc[1], lz))
		fout.close()
		#
		return None
	#
	def innerouterpolys(self, polylist):
		#print "do nothing yet."
		# but, eventually: separate polygons that are "inside" two two polys (aka, and "outside" poly), and
		# associate their inner polys with them.
		# use shapely.geometry (sgp) module?? maybe not, as per compatiblity.
		#
		# note, pollys are coming like [ [[Xs], [Ys]], [] ], 
		# so the jth vertex of the ith poly is (x,y) = (polylist[i][0][j], polylist[i][1][j])
		polylistplus=[]		# indexed entries: [polyindex, [list of inners], [verts] ]
		#outerpolys=[]			# these will be lists. each entry is like: [[outer],[true-inner],[true-inner],..]
		#
		#for ply in polylist:
		for i in xrange(len(polylist)):
			polylistplus += [[i, [], polylist[i]]]
			#
			# shorthand:
			# which poly each polygon it is inside, and how many polys it is inside (aka, len(that list)).
			#x0,y0=poly[i][0][0], poly[i][1][0]	# we only need to test one point since we're starting with contours (don't cross).
			#x0,y0=polylistplus[-1][2][0][0], polylistplus[-1][2][1][0]	# we only need to test one point since we're starting with contours (don't cross).
			x0, y0 = polylist[i][0][0], polylist[i][1][0]
			#print "x0, y0: ", x0, y0
			# in general, we'd need to test all points to be inside.
			#
			# for each polygon in this level:
			# is the ith ("top") polygon inside the j'th poly?
			for j in xrange(len(polylist)):
				if j==i: continue 
				X,Y = polylist[j][0][:], polylist[j][1][:]
				#if x0>=max(X) or x0<=min(X) or y0>max(Y) or y0<min(Y): 
				#	print "outside max/min..."
				#	continue
				#
				if X[0]!=X[-1]:
					X+=[X[0]]
					Y+=[Y[0]]	# complete the poly...
				#
				N=len(X)
				ncrossings = 0
				# how many poly boundaries do we cross if we draw a line out of the poly in one direction.
				# equivalently (and in computer language), how many segments at y1 < y <y2 (or upside down)
			# are to the right of the point (or to the left, or up/down -- pick one).
				for k in xrange(1,N):
					k1 = k-1
					#k2 = (k+1)%N	# note the k%N: the poly does not have to be closed
					k2 = k	# but it should be, or you can count a crossing twice and get a bogus answer.
					x1, y1 = X[k1], Y[k1]
					x2, y2 = X[k2], Y[k2]
					'''
					x1,y1 = polylist[j][0][k1], polylist[j][1][k1]
					print "j,len(polylist):", j, len(polylist)
					print "k2, len(polylist[j][0])", k2, len(polylist[j][0])
					print "k2, len(polylist[j][1])", k2, len(polylist[j][1])
					x2,y2 = polylist[j][0][k2], polylist[j][1][k2]
					'''
					#if y0>=min(y1, y2) and y0<=max(y1, y2) and x0<max(x1, x2):
					#
					if x0>=min(x1, x2) and x0<max(x1, x2):	# note, one must be <= and the other < or, if we're on a grid -- and
						fx = (y2-y1)/(x2-x1)*(x0-x1)			# we're always on a grid, we'll count each crossing (left seg then right).
						if fx>(y0-y1):							# that's why the test script worked (x1,y1 != x2, y2) and the production
							ncrossings += 1						# failed...
				# clean up a bit:
				X=None
				Y=None
				#
				#print i,j,j,ncrossings, ncrossings%2
				if ncrossings%2==1:
					# i'th poly is inside j'th poly...
					polylistplus[-1][1] += [j]	
				#
			#
			# so now, we have a list of polygons and the polys they are inside.
			# for outerPolys, len(innersList)%2==0. if polyA is inside polyB and nA-nB=1, polyA is an inner-poly to polyB
			#
			#for rw in polylistplus:
			#	print rw
			outerpolys=[]
			for ply1 in polylistplus:
				#print "inner-len: ", len(ply1[1]), len(ply1[1])%2
				if len(ply1[1])%2==0:
					#print "***", len(ply1[1])
					# it's an outer poly...
					outerpolys+=[[ply1[2]]]
					for ply2 in polylistplus:
						# find its inners:
						if ply2==ply1: continue
						if len(ply2[1])%2==0: continue	# skip outers...
						#
						#print len(ply2[1]), (len(ply1[1])+1)
						if ply1[0] in ply2[1] and len(ply2[1])==(len(ply1[1])+1):
							# the outer poly's index is in ply2's "inside-list", then ply2 is inside ply...
							# AND, ply2 is one "deeper" (inside exactly one more poly) thatn ply
							outerpolys[-1]+=[ply2[2]]
		#
		#return polylist
		return outerpolys
	#
	def contsToPlotLists(self, cset=None):
		# this corrects a mpl/kml mismatch. mpl likes to plot inner polygons by drawing a line from the outside
		# ring to the inner poly, draw the inenr poly, then back to the outer ring over the first line. mpl
		# interprets the line as infinitely thin and basically ignores it. KML will not ignore that line, and so
		# you get a big fat mess. This script separates the outer/inner groups into outer and inner poly rings.
		#
		# BUT, (2012 09 05) this needs to be revisited. inner-inner polys are not plotting correctly. aka,
		# the z-profile like: /-\_/-\_/-\ (like an s0 within an s1) need to be considered. in such a case,
		# we get: <outer><inner><inner></inner></inner></outer>, but there should be 1 outer/inner poly then a
		# separate poly entirely. the second "inner" needs to be pulled out as a separate poly. (see similar note below)
		#
		# so, polys that are inside an odd number of polys are "inner" polys; polys inside an even number are "outside" 
		# any poly (as per knot theory), or otherwise constitute an "outer" poly. poly1 is an "inner" poly of poly2 if:
		#   1: it is inside poly2
		#   2: its "inner index n2 = n1 - 1
		if cset==None: cset=self.conts
		cs=cset.collections
		outlist=[]	# the length of this will be the number of levels: (or contours?)
						# outlist -> [ [level0: [[cont0x], [cont0y]], [[cont1x], [cont1y]] ], [level1: ] ???
		#levels=[]
		#contlevel=0
		#
		for i in xrange(len(cs)):
			# level-level:
			outlist+=[[]]
			#contcount=0
			# each level will have multiple polygons:
			for trace in cs[i].get_paths():
				# each "trace" will be a polygon (there might be multiple polygons per level), which might have "internal" polygons.
				tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				#outlist[-1]+=[[[],[]]]
				newpoly=[ [[], []] ]	# newpoly[poly-index][x=0 or y=1][x_i or y_i]. note, outer poly is first element.
				# draw the polygons (first the outside one, then the inside ones -- presumably)
				llindex=0
				for lng,lat in trace.vertices:
					#
					# fix mangled contours:
					#
					newpoly[-1][0]+=[lng]
					newpoly[-1][1]+=[lat]
					#
					# have we closed the polygon (in the middle?):
					# what is happening is that the inner polys are being drawn screwy ways.
					# the errors are not a total mess; basically, the link between then inner and outer poly is drawn
					# from the last/first point rather than clostest two points.
					# the first poly is probably the outer-most poly. subsequent closed sequences are inside polys.
					# however, see comments above. the inner-inner->outer polys must be separated. do the following procedure
					# first, then pass that list of polys for further parsing.
					#
					# SO, make a list of all closed objects in a single poly object, then insert them into the outer sequence
					# (noting the inner poly has to close and then return to the starting point on the outer ring.
					# is newpoly[-1] == newpoly[0] ?
					if newpoly[-1][0][-1] == newpoly[-1][0][0] and newpoly[-1][1][-1] == newpoly[-1][1][0] and len(newpoly[-1][0])>1:
						# the polygon has closed.
						# but has it really, or is this just a "bay" or a "nub" for which the closing boundaries intersect.
						# if this poly is ended and a new poly folows, the next point is a line-segment connecting the next
						# (inner) poly. if it is a bay or a nub, the shape will close on itself again. I THINK, that we can say
						# this poly is truly closed IF this end coordinate does not appear again in the sequence.
						# there is probably a smart, compiled way to do this, but for now let's just loop it:
						ispinched=True
						for i in xrange(llindex+1, len(trace.vertices)):
							if trace.vertices[i][0]==lng and trace.vertices[i][1]==lat:
								# we've found our "end" point at least one mor time, so keep building the poly.
								# obviously, there is a faster way to do this, but let's just see if it works first.
								ispinched=False
								break
						#
						if ispinched==True:
							# pinch it off and start a new poly.
							newpoly+=[ [[],[]] ]
							tmpLL=[]
							#contcount+=1
					#
					llindex+=1
				#
				# now, clean up a bit:
				# if it's a bogus poly -- probably a line back to a prev. poly (each poly is like: poly0, line (segment) to poly1,
				#		poly1, either: {line back to poly0, line to poly2, nothing?}
				# i think the problem we're having is multi-inner polys aka, profile like 1 - 0 - 1. (ring=1, dip-ring=0, bump in middle=1)
				# ___|-\_/-\_/-|___ . the inner-inner bit is causing a problem. to fix this, i think, we have to draw all the polys
				# and explicitly determine which polys they are inside. (<outer><inner><inner></innter></inner></outer> ??)
				# note: a poly-ring inside an even  number of polys is an outer ring; a poly inside an odd number of rings is an inner ring.
				if len(newpoly[-1][0])<2: newpoly.pop()
				if newpoly[-1][0][-1] != newpoly[-1][0][0] and newpoly[-1][1][-1] != newpoly[-1][1][0]:
					newpoly.pop()
				#
				# at this point, newpoly is like: [[likely-outermost],[maybe-inner], [maybe-inner], ... ]
				# where we infer, or have in the past, that the first row represents the outer-most poly; subsequent polys are inners.
				# this is not the cast. find the 2xinner (aka, outer) polys within as per comments above.
				#
				# break up newpoly into inner-outer groups
				newpolylist=self.innerouterpolys(newpoly)
				#if len(newpoly)>2:
				#	print "len(newpoly), len(newpolylist): %d, %d" % (len(newpoly), len(newpolylist))
					#print newpoly
					#print newpolylist
				#newpolylist = [newpoly]
				# return the full list and interpret i=0 as the outer ring; i>0 as the inner boundaries.
				#
				#outlist[-1]+=[newpoly]	# add this trace(set) to the last level of outlist.
				for ply in newpolylist:
					#outlist[-1]+=[newpoly]
					outlist[-1]+=[ply]
				#
				#contcount+=1
			#contlevel+=1
		return outlist
	
	#
	def contsToPlotListsRaw(self, cset=None):
		# this corrects a mpl/kml mismatch. mpl likes to plot inner polygons by drawing a line from the outside
		# ring to the inner poly, draw the inenr poly, then back to the outer ring over the first line. mpl
		# interprets the line as infinitely thin and basically ignores it. KML will not ignore that line, and so
		# you get a big fat mess. This script separates the outer/inner groups into outer and inner poly rings.
		#
		# BUT, (2012 09 05) this needs to be revisited. inner-inner polys are not plotting correctly. aka,
		# the z-profile like: /-\_/-\_/-\ (like an s0 within an s1) need to be considered. in such a case,
		# we get: <outer><inner><inner></inner></inner></outer>, but there should be 1 outer/inner poly then a
		# separate poly entirely. the second "inner" needs to be pulled out as a separate poly. (see similar note below)
		#
		# so, polys that are inside an odd number of polys are "inner" polys; polys inside an even number are "outside" 
		# any poly (as per knot theory), or otherwise constitute an "outer" poly. poly1 is an "inner" poly of poly2 if:
		#   1: it is inside poly2
		#   2: its "inner index n2 = n1 - 1
		if cset==None: cset=self.conts
		cs=cset.collections
		outlist=[]	# the length of this will be the number of levels: (or contours?)
						# outlist -> [ [level0: [[cont0x], [cont0y]], [[cont1x], [cont1y]] ], [level1: ] ???
		#levels=[]
		#contlevel=0
		#
		for i in xrange(len(cs)):
			# level-level:
			outlist+=[[]]
			#contcount=0
			# each level will have multiple polygons:
			for trace in cs[i].get_paths():
				# each "trace" will be a polygon (there might be multiple polygons per level), which might have "internal" polygons.
				tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				#outlist[-1]+=[[[],[]]]
				newpoly=[ [[], []] ]	# newpoly[poly-index][x=0 or y=1][x_i or y_i]. note, outer poly is first element.
				# draw the polygons (first the outside one, then the inside ones -- presumably)
				for lng,lat in trace.vertices:
					#
					# fix mangled contours:
					#
					newpoly[-1][0]+=[lng]
					newpoly[-1][1]+=[lat]
					#
				#
				outlist[-1]+=[newpoly]	# add this trace(set) to the last level of outlist.
				#
				#contcount+=1
			#contlevel+=1
		return outlist	
	#
	def contsToPlotListsOuter(self, cset=None):
		# simple polygons, not distinguishing inner/outer.
		if cset==None: cset=self.conts
		cs=cset.collections
		outlist=[]	# the length of this will be the number of levels: (or contours?)
						# outlist -> [ [level0: [[cont0x], [cont0y]], [[cont1x], [cont1y]] ], [level1: ] ???
		contlevel=0
		#
		for i in xrange(len(cs)):
			# level-level:
			outlist+=[[]]
			contcount=0
			for trace in cs[i].get_paths():
				tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				outlist[-1]+=[[[],[]]]
				for lng,lat in trace.vertices:
					#
					# fix mangled contours:
					#
					outlist[-1][-1][0]+=[lng]
					outlist[-1][-1][1]+=[lat]
					tmpLL+=[[lng, lat]]
					#
					# have we closed the polygon (in the middle?):
					# what is happening is that the inner polys are being drawn screwy ways.
					# the errors are not a total mess; basically, the link between then inner and outer poly is drawn
					# from the last/first point rather than clostest two points.
					# the first poly (i think) is the outer poly. subsequent closed sequences are inner polys (excluded regions).
					# SO, make a list of all closed objects in a single poly object, then insert them into the outer sequence
					# (noting the inner poly has to close and then return to the starting point on the outer ring.
					if outlist[-1][-1][0][-1] == outlist[-1][-1][0][0] and outlist[-1][-1][1][-1] == outlist[-1][-1][1][0] and len(outlist[-1][-1][0])>1:
						# the polygon has closed. pinch it off and start a new poly.
						outlist[-1]+=[[[],[]]]
						tmpLL=[]
						contcount+=1

				if len(outlist[-1][-1][0])<2: outlist[-1].pop()
				if outlist[-1][-1][0][-1] != outlist[-1][-1][0][0] and outlist[-1][-1][1][-1] != outlist[-1][-1][1][0]:
					outlist[-1].pop()
				# (outlist[contour level][polynum][x or y][x,y index] )
				
				
				'''
				#if tmpLL[0]!=tmpLL[-1]:
				while tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
					#fixedpolys+=1
					#outlist[-1][-1][0] += [tmpLL[0][0]]
					#outlist[-1][-1][1] += [tmpLL[0][1]]
					outlist[-1][-1][0].pop()
					outlist[-1][-1][1].pop()
					tmpLL.pop()
				'''
					#
				contcount+=1
			contlevel+=1
		return outlist

	def plotPolyList(self, polys=None, markerstr='.-', fignum=3, pllevels=None):
		if polys==None: polys=self.contsToPlotLists()
		if pllevels==None: pllevels = range(len(polys))
		if type(pllevels).__name__ in ('float', 'int'): pllevels=[pllevels]
		plt.figure(fignum)
		plt.clf()
		plt.ion()
		nlevels=len(polys)
		
		#
		ilevel=0
		for level in polys:
			#lvlclr=self.plotcolor(ilevel, nlevels)
			if ilevel not in pllevels:
				ilevel+=1
				continue
				
				 
			lvlclr=plotcolor(ilevel, nlevels*2)
			print ilevel, lvlclr
			for xy in level:
				for ring in xy:
				#fixedpolys = checkPoly(inpoly=xy, fignum=None)
				#for fxy in fixedpolys:
					plt.plot(ring[0], ring[1], markerstr, color=lvlclr)
				#plt.plot(xy[0], xy[1], markerstr, color=lvlclr)
			ilevel+=1
		return polys			
	
	def plotPolyListOuter(self, polys=None, markerstr='.-', fignum=3, pllevels=None):
		# plot the simpler "just outers" version of the polygon list.
		if polys==None: polys=self.contsToPlotListsOuter()
		if pllevels==None: pllevels = range(len(polys))
		if type(pllevels).__name__ in ('float', 'int'): pllevels=[pllevels]
		plt.figure(fignum)
		plt.clf()
		plt.ion()
		nlevels=len(polys)
		
		#
		ilevel=0
		for level in polys:
			#lvlclr=self.plotcolor(ilevel, nlevels)
			if ilevel not in pllevels:
				ilevel+=1
				continue
				
				 
			lvlclr=plotcolor(ilevel, nlevels*2)
			print ilevel, lvlclr
			for xy in level:
				#fixedpolys = checkPoly(inpoly=xy, fignum=None)
				#for fxy in fixedpolys:
				#	plt.plot(fxy[0], fxy[1], markerstr, color=lvlclr)
				plt.plot(xy[0], xy[1], markerstr, color=lvlclr)
			ilevel+=1
		return polys
	
	def getmarkerKMLstr(self, markers):
		#kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
		kmlstr=''
		for ev in markers:
			thisdt=str(mpd.num2date(ev[0]))
			datestr=thisdt.split('.')[0]
			kmlstr+='<Placemark><name>%s, %.2f</name><Point><coordinates>%f,%f,0</coordinates></Point></Placemark>' % (datestr, ev[3], ev[2], ev[1])
		
		return kmlstr
		
	def writeKMLfile(self, cset=None, fout='conts.kml', colorbarname='scale.png', equakes=None, top=1.0, bottom=0.0):
		#kmlstr=self.KMLstrFromConts(cset)
		kmlstr=self.KMLstrFromConts(cset, colorbarname=colorbarname, openfile=True, closefile=False, top=top, bottom=bottom)
		if equakes !=None:
			# we've been given some earthquakes to plot as well.
			kmlstr+='\n'
			kmlstr+=self.getmarkerKMLstr(equakes)
		kmlstr+='</Document>\n'
		kmlstr+='</kml>'
		#
		f=open(fout, 'w')
		f.write(kmlstr)
		f.close()
		return None
	
	def KMLstrFromConts(self, cset=None, colorbarname='scale.png', openfile=True, closefile=True, warnings=None, top=1.0, bottom=0.0):
		# get a KML string from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
		# (will it be necessary to write directly to file? is there a string length limit problem?)
		#
		if cset==None:
			cset=self.conts
			#cset=plt.contourf(self.X_i, self.Y_i, self.Z2d, self.contres*self.contfact, cm=plt.spectral())
		cs=cset.collections
		#
		polys = self.contsToPlotLists(cset=cset)	# this function fixes multi-closed polys. use for actual kml contours.
																# still has same polys[cont. level][contour index][x=0,y=1][x,y index] format
																# the polys array will differ from cset.collections in the number of polys
																# per level; the number of levels should be the same.
		#
		# contour range to render:
		if top==None: top=1.0
		if bottom==None: bottom=0.0
		#
		# if top, bottom are >1, assume they're the number/index of the contour range being selected.
		# otherwise, it's a percentage.
		# first contour is the highest value.
		#
		levels = cset.levels
		n_levels = len(levels)-1
		kml_top      = int(math.ceil(top*float(n_levels)))
		kml_bottom   = int(math.floor(bottom*float(n_levels)))	# (and just in case we accidentally calc <0 or >len)
		#
		print "kml bottom, top: ", kml_bottom, kml_top
		#polys  = self.contsToPlotLists(cset=cset[top:bottom])
		#
		#resolution = 5
		#LevelsNumber = 5 * resolution
		if warnings==None: warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#contoursStart=int(.2*len(cs))
		#resolution = int(len(cs)/len(warnings))
		resolution = int(len(polys)/len(warnings))
		#startindex=resolution
		startindex=0
		kmlstr=''
		#
		if openfile==True:
			kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
		#
		# styles come from the contours collection (cs):
		#for i in xrange(0,len(cs)):
		for i in xrange(kml_bottom, kml_top):
			bgra_array = 255*cs[i].get_facecolor()[0]
			kmlstr+='<Style id="l%d">\n' % i
			kmlstr+='<LineStyle><color>00000000</color></LineStyle>\n'
			kmlstr+='<PolyStyle><color>7d%02x%02x%02x</color></PolyStyle>\n' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
			kmlstr+='</Style>\n'
		#
		#
		kmlstr+='<ScreenOverlay id="scale">\n'
		kmlstr+='<name>Color Scale</name>\n'
		kmlstr+='<Icon><href>scale.png</href></Icon>\n'
		kmlstr+='<Icon><href>%s</href></Icon>\n' % colorbarname
		kmlstr+='<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n'
		kmlstr+='</ScreenOverlay>\n'
		#
		#fixedpolys=0
		# now, stake out the contours from the polys array.
		# note: len(cs) = len(polys) (same number of levels)
		#	   len(cs[i]) != len(polys[i]) (in fact, len(polys[i]>=len(cs[i]) ), which is to say
		#	   polys may have more distinct polygons per contour level. of course, the individual polys
		#	   will have differend length as well.
		#for i in xrange(startindex, len(cs)):
		#for i in xrange(startindex, len(polys)):
		for i in xrange(kml_bottom, kml_top):
			# each i is a contour level.
			kmlstr+='<Placemark>\n'
			#print i, resolution, len(warnings), len(polys)
			warningindex=int(float(i)/float(resolution))
			if warningindex>(len(warnings)-1): warningindex=len(warnings)-1	# in the event that we have off-integer numbers.
			kmlstr+='<name>%s Risk</name>\n' % warnings[warningindex]
			kmlstr+='<styleUrl>#l%d</styleUrl>\n' % i
			kmlstr+='<MultiGeometry>\n'
			 
			#for trace in cs[i].get_paths():
			for ii in xrange(len(polys[i])):
				# each ii is a polygon (set).
				kmlstr+='<Polygon>\n'
				kmlstr+='<extrude>0</extrude>\n'
				kmlstr+='<altitudeMode>clampToGround</altitudeMode>\n'
				kmlstr+='<outerBoundaryIs>\n'
				kmlstr+='<LinearRing>\n'
				kmlstr+='<coordinates>\n'
				
				#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				#for lng,lat in trace.vertices:
				for ill in xrange(len(polys[i][ii][0][0])):
					# first set is the outerBoundary
					# noting that each polygon is stored as [[x], [y]], so len(polys[i][ii])=2 always.
					# len(polys[i][ii][0])=len(polys[i][ii][1]) is the length or size (num. verts.) of the polygon.
					lng=polys[i][ii][0][0][ill]
					lat=polys[i][ii][0][1][ill]
					#
					kmlstr+='%f,%f,0\n' % (lng, lat)
					#tmpLL+=[[lng, lat]]
				#if tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
				#	fixedpolys+=1
				#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])

				kmlstr+='</coordinates>\n</LinearRing>\n</outerBoundaryIs>\n'
				#
				# inner polys?
				for iInner in xrange(1,len(polys[i][ii])):
					thispoly=polys[i][ii][iInner]
					# if any, these will be inner polys, each like [[x], [y]]
					kmlstr+='<innerBoundaryIs>\n<LinearRing>\n<coordinates>\n'
					# coords like 'lon,lat,alt\n'
					# -77.05668055019126,38.87154239798456,100
					for ill in xrange(len(thispoly[0])):
						kmlstr+='%f,%f,0\n' % (thispoly[0][ill], thispoly[1][ill])
					kmlstr+='</coordinates>\n</LinearRing>\n</innerBoundaryIs>\n'
				#
				kmlstr+='</Polygon>\n'
			#
			kmlstr+='</MultiGeometry>\n'
			kmlstr+='</Placemark>\n'
		#
		if closefile==True:
			kmlstr+='</Document>\n'
			kmlstr+='</kml>'
		#
		return kmlstr

	#
	def KMLstrFromConts_raw(self, cset=None, colorbarname='scale.png'):
		# get a KML string from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
		# (will it be necessary to write directly to file? is there a string length limit problem?)
		#
		if cset==None:
			cset=self.conts
			#cset=plt.contourf(self.X_i, self.Y_i, self.Z2d, self.contres*self.contfact, cm=plt.spectral())
		cs=cset.collections
		#
		#resolution = 5
		#LevelsNumber = 5 * resolution
		warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#contoursStart=int(.2*len(cs))
		resolution = int(len(cs)/self.contfact)
		#startindex=resolution
		startindex=0
		#
		#
		kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
		#
		for i in range(0,len(cs)):
			bgra_array = 255*cs[i].get_facecolor()[0]
			kmlstr+='<Style id="l%d">\n' % i
			kmlstr+='<LineStyle><color>00000000</color></LineStyle>\n'
			kmlstr+='<PolyStyle><color>7d%02x%02x%02x</color></PolyStyle>\n' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
			kmlstr+='</Style>\n'
		#
		#
		kmlstr+='<ScreenOverlay id="scale">\n'
		kmlstr+='<name>Color Scale</name>\n'
		#kmlstr+='<Icon><href>scale.png</href></Icon>\n'
		kmlstr+='<Icon><href>%s</href></Icon>\n' % colorbarname
		kmlstr+='<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n'
		kmlstr+='</ScreenOverlay>\n'
		#
		fixedpolys=0
		for i in xrange(startindex, len(cs)):
			kmlstr+='<Placemark>\n'
			kmlstr+='<name>%s Risk</name>\n' % warnings[i/resolution]
			kmlstr+='<styleUrl>#l%d</styleUrl>\n' % i
			kmlstr+='<MultiGeometry>\n'
			 
			for trace in cs[i].get_paths():
				kmlstr+='<Polygon>\n'
				kmlstr+='<extrude>0</extrude>\n'
				kmlstr+='<altitudeMode>clampToGround</altitudeMode>\n'
				kmlstr+='<outerBoundaryIs>\n'
				kmlstr+='<LinearRing>\n'
				kmlstr+='<coordinates>\n'
				
				tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				for lng,lat in trace.vertices:
					#kmlstr+='%s,%s,0\n' % (round(lng,7),round(lat,7))
					
					# fixing contours mangled by pyplot:
					# an older attempted fixing routine that didn't quite work out.
					# keep this version as the "raw" form; see other KML function that 
					# uses a repaired set of polygons.
					#
					thislng=round(lng, self.zresolution)
					thislat=round(lat, self.zresolution)
					kmlstr+='%f,%f,0\n' % (thislng, thislat)
					tmpLL+=[[thislng, thislat]]
				#if tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
				#	fixedpolys+=1
				#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])

				kmlstr+='</coordinates>\n'
				kmlstr+='</LinearRing>\n'
				kmlstr+='</outerBoundaryIs>\n'
				kmlstr+='</Polygon>\n'
			#
			kmlstr+='</MultiGeometry>\n'
			kmlstr+='</Placemark>\n'
		#
		kmlstr+='</Document>\n'
		kmlstr+='</kml>'
		#
		return kmlstr
	
	def makeColorbar(self, cset=None, colorbarname=None, reffMag=5.0, fnum=None, fontcolor='k'):
		#
		# reffMag:
		if reffMag==None: reffMag=self.reffMag
		#
		if colorbarname==None: colorbarname='scale.png'
		if cset==None: cset=self.conts
		cs=cset.collections
		#resolution = 5
		#LevelsNumber = 5 * resolution
		warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#contoursStart=int(.2*len(cs))
		resolution = int(len(cs)/self.contfact)
		#startindex=resolution
		startindex=0
		#
		plt.figure(fnum)
		plt.close()
		#
		fig = plt.figure(num=fnum, figsize=(1,3))
		plt.clf()
		axe = fig.add_axes([0.05, 0.05, 0.2, 0.9])
		#
		fig.figurePatch.set_alpha(0)
		#
		matplotlib.rc('font', family='monospace', weight='black')
		#
		cmap = mpl.cm.spectral
		#norm = mpl.colors.Normalize(vmin=Z_i.min(), vmax=Z_i.max())
		#ratefactor=10**(self.mc-self.reffMag)
		ratefactorExp=self.mc-reffMag	# so the raw rate probably should be with mc_catalog = mc_etas.
													# but, we'll want to report some rate of m>reffMag
		norm = mpl.colors.Normalize(vmin=self.Z2d.min(), vmax=self.Z2d.max())
		timefactExp = math.log10(year2secs)
		#
		#print self.Z2d.min(), self.Z2d.max(), norm
		#tics = [0,10,20,40,80]
		#tics = [round(self.Z2d.min(),3), round(self.Z2d.max(),3)]
		#t1='%.2f' % (self.Z2d.min() + timefactExp + ratefactorExp)
		#t2='%.2f' % (self.Z2d.max() + timefactExp + ratefactorExp)
		#
		# let's add intermediate ticks: half and half-half:
		ticdiff0=self.Z2d.max()-self.Z2d.min()
		ticdiff=ticdiff0/4.0
		#t3='%.2f' % (self.Z2d.min() + timefactExp + ratefactorExp + 1.0*ticdiff
		#t4='%.2f' % (self.Z2d.min() + timefactExp + ratefactorExp + 2.0*ticdiff
		#t5='%.2f' % (self.Z2d.min() + timefactExp + ratefactorExp + 3.0*ticdiff
		#
		print "max set: ", self.Z2d.max(), timefactExp, ratefactorExp
		#tics = [self.Z2d.min(), self.Z2d.max()]
		tics = [self.Z2d.min()]
		tcklbls=['%.2f'% (self.Z2d.min() + timefactExp + ratefactorExp)]
		while tics[-1]<=self.Z2d.max():
			tics+=[tics[-1]+ticdiff]
			tcklbls+=['%.2f'% (tics[-1]+ timefactExp + ratefactorExp)]
		#print t1, t2
		#
		#cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%g%%", orientation='vertical')
		cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%.2f", orientation='vertical')
		#cb1.set_ticklabels([t1, t2], update_ticks=True)
		cb1.set_ticklabels(tcklbls, update_ticks=True)
		#
		cb1.set_label('ETAS rate: $\\log({N(m \geq %.1f)}/{yr \cdot km^2})$' % reffMag, color=fontcolor)
		#cb1.set_label('ETAS rate')
		plt.savefig(colorbarname)	# vertical
		#
		#suffix=''
		i=len(colorbarname)-1
		#while colorbarname[i]!='.':
		#	#suffix.insert(0, colorbarcopy[-i])
		#	# this will get us the last '.'
		#	i-=1
		i=colorbarname.rfind('.')
		hname=colorbarname[:i] + '-h' + colorbarname[i:]
		#
		im = ipp.open(colorbarname, 'r')
		im.load()							# just to be sure. open() only creates a pointer; data are not loaded until analyzed.
		im=im.rotate(-90)
		#im.show()
		im.save(hname)
		#
		# depricated method:
		#os.system('convert -rotate 90 %s %s' % (colorbarname, hname))
		#os.system('cp %s %s' % (colorbarname, hname))
		#
		return [colorbarname, hname]
	#
	def writeCat(self, fout='cat.cat'):
		# dt, lat, lon, mag
		f=open(fout, 'w')
		f.write('#BASScast catalog\n#dtm, lat, lon, mag\n')
		f.write('#!writedate\t%s\n' % str(dtm.datetime.now(pytz.timezone('UTC'))))
		f.write('#!lats=%s\tlons=%s\tfcdate=%s\t%f\n' % (str(self.latrange), str(self.lonrange), self.fcdate, self.fcdatef))
		for rw in self.catalog:
			#row_string = '%f\t%f\t%f\t%f' % (rw[0], rw[1], rw[2], rw[3])
			row_string = ''
			for x in rw:
				row_string += '%s\t' % str(x)
			row_string = row_string[:-1] + '\n'
			#
			f.write(row_string)
		f.close()

#class locbase(object):
class locbase(mpp.Process):
	loc = [0,0]		# [lon, lat]
	#eventDtm = None
	eventftime = None # float version of time
	mc=2.0
	#
	def __init__():
		mpp.Process.__init__()
	#
	# location related functions:
	def timesince(self, intime):
		# return time since the in-event (in seconds):
		#print type(self)
		# this will be a guess until "type()" can be properly implemented, but...
		if type(intime)==type(self):
			# we've been asked to compare with an earthquake object.
			# but the way we SHOULD do this is to check for a member function/variable...
			intime=intime.eventftime
		#
		# otherwise, assume float ops will work.
		#
		tsince = intime - self.eventftime
		#
		return tsince
		
	def dispFrom(self, inloc, coordtype='cart'):
		# displacement from inloc...
		# inloc is a vector [lon, lat]
		# return a vector [dLon, dLat] or [r, theta]
		# return distances in km.
		#
		#meanlat=.5*(inloc[1]+self.loc[1])*2.0*math.pi/360.0
		#latfactSelf=math.cos(self.loc[1]*2.0*math.pi/360.0)
		#latfactin = math.cos(inloc[1]*2.0*math.pi/360.0)
		latfactSelf=math.cos(self.loc[1]*deg2rad)
		#latfactin = math.cos(inloc[1]*deg2rad)
		#meanlat=.5*(self.loc[1] + inloc[1])
		#meanlatfact=math.cos(meanlat*2.0*math.pi/360.0)
		#
		dx = (inloc[0] - self.loc[0])*111.12*latfactSelf		# the dist. along self's latitude.
		dy = (inloc[1] - self.loc[1])*111.12
		rvect=[dx, dy]	# return in rect. coords... but this is, i think, not terribly accurate...
		Rsph = self.sphericalDist(inloc)		# this will give the proper spherically calculated distance.
		# using geographic-lib (which is very computationally expensive):
		#g1=ggp.WGS84.Inverse(self.loc[1], self.loc[0], inloc[1], inloc[0])	# lat, lon, lat, lon (??)
		#Rgeolib=g1['s12']/1000.0	# and converting from m to km.
		#
		Rrect=math.sqrt(dx*dx + dy*dy)		# not totally accurate, but quick and easy
		#
		if coordtype=='radial' or coordtype=='rad':
			#Rsph = self.sphericalDist(inloc)		# this will give the proper spherically calculated distance. (and now it's easy).
			# use one of these two...
			R=Rsph
			#R=Rgeolib
			if R==0:
				# should not happen often except in demos. just return 0.
				theta=0.0
			else:
				theta = math.acos(dx/Rrect)	# and we'll need a better measurement of theta as well (later -- i'm not sure we use it
												# anywhere yet).
				if dy<0: theta=-theta
			#rvect=[R,theta]
			rvect = [R, theta]	# return in radial coords; theta in radians
		#
		return rvect
	
	def sphericalDist(self, inloc):
		# displacement from inloc...
		# inloc is a vector [lon, lat]
		# return a vector [dLon, dLat] or [r, theta]
		# return distances in km.
		#
		# also, we need to get the proper spherical angular displacement (from the parallel)
		#
		Rearth = 6378.1	# km
		deg2rad=2.0*math.pi/360.
		#
		# note: i usually use phi-> longitude, lambda -> latitude, but i accidentally copied a source where this is
		# reversed. oops. so just switch them here.
		#
		#phif  = inloc[0]*deg2rad
		#lambf = inloc[1]*deg2rad
		#phis  = self.loc[0]*deg2rad
		#lambs = self.loc[1]*deg2rad
		phif  = inloc[1]*deg2rad
		lambf = inloc[0]*deg2rad
		phis  = self.loc[1]*deg2rad
		lambs = self.loc[0]*deg2rad
		#
		#print 'phif: ', phif
		#print 'lambf: ', lambf
		#
		dphi = (phif - phis)
		dlambda = (lambf - lambs)
		#this one is supposed to be bulletproof:
		sighat3 = math.atan( math.sqrt((math.cos(phif)*math.sin(dlambda))**2.0 + (math.cos(phis)*math.sin(phif) - math.sin(phis)*math.cos(phif)*math.cos(dlambda))**2.0 ) / (math.sin(phis)*math.sin(phif) + math.cos(phis)*math.cos(phif)*math.cos(dlambda))  )
		R3 = Rearth * sighat3
		#
		return R3
	
	def distFrom(self, inloc, coordtype='cart'):
		return self.dispFrom(inloc, coordtype)

	def setEventTime(self, evtime):
			#print type(evtime).__name__
		#
		#print "setting EventTime: ", evtime, type(evtime), isinstance(evtime, float)
		#if type(evtime).__name__ == 'datetime': # == type(dtm.datetime):
		if isinstance(evtime, dtm.datetime):
			#self.eventDtm = evtime
			self.eventftime = mpd.date2num(evtime) * days2secs	# catalog must be in seconds for scaling to work.
			#
		if isinstance(evtime, int) or isinstance(evtime, float):
		#if type(evtime) == type(1.0) or type(evtime) == type(1):
			self.eventftime = float(evtime)
			#
			#if self.eventftime < (mpd.date2num(dtm.datetime.now(pytz.timezone('UTC')))*100.0):
			if self.eventftime < 10**6.5:
				# date is in days, not seconds
				self.eventftime*=days2secs
				#
			# this fails for large dates (in simulated catalogs) and we don't seem to ever use it...
			#self.eventDtm = mpd.num2date(self.eventftime/days2secs)	# assuming it is properly in seconds, as it should be.
			# 
		#if type(evtime) == type('abc'):
		if isinstance(evtime, str):
			# convert the string, etc.
			#self.eventDtm = self.datetimeFromString(strDtin=evtime)
			self.eventftime = mpd.datestr2num(evtime)*days2secs
			#self.eventDtm = mpd.num2date(self.eventftime/days2secs)			
			#
		return None

	# utils we probably don't need any longer (use matplotlib.date.{various date conversions}

class objContour(object):
	lwr=None
	upr=None
	#clr='7d000000'	# white? (alpha - R - G - B)
	alpha = '7d'
	R='00'
	G='00'
	B='00'
	
	verts=[]
	
	def __init__(self, low=None, high=None, alpha=None, R=None, G=None, B=None, verts=None, RGB=None):
		return self.initialize(low=low, high=high, alpha=alpha, R=R, G=G, B=B, verts=verts, RGB=RGB)
		#
	#
	def initialize(self, low=None, high=None, alpha=None, R=None, G=None, B=None, verts=None, RGB=None):
		#print "RGBinit: %s" % RGB
		if RGB!=None:
			self.setRGB2(RGB=RGB)
		else:
			self.setRGB(R=R, G=G, B=B, alpha=alpha)
		self.lwr=low
		self.upr=high
		#
		self.verts = verts
		#
		return None
		#
	#
	def getRGB(self):
		return '%s%s%s%s' % (self.alpha, self.R, self.G, self.B)
	
	def setRGB(self, R='00', G='00', B='00', alpha='7d'):
		R='00' + str(R)
		G='00' + str(G)
		B='00' + str(B)
		alpha='00' + str(alpha)
		#
		R=R[-2:]
		G=G[-2:]
		B=B[-2:]
		alpha=alpha[-2:]
		#
		self.R=R
		self.G=G
		self.B=B
		self.alpha=alpha
		#
		return None
	
	def setRGB2(self, RGB='7d000000'):
		#if type(RGB).__name__ != 'string': return None
		#print 'RGB: %s' % RGB
		#
		self.alpha=RGB[0:2]
		self.R=RGB[2:4]
		self.G=RGB[4:6]
		self.B=RGB[6:8]
		#
		return None
	
	def aryRow(self):
		return [self.lwr, self.upr, self.getRGB(), self.verts]

class forecastsite(locbase):
	# site, or grid-cell, for a forecast.
	# basically, it needs to hava location (position and time), size, and z-value (forecast metric).
	# nominally, it should also have some catalog prams, like mc to determine an earthquake count.
	# it does beg the question as to whether scaling prams like rho, alpha, etc. are properties of the
	# earthquake or the local earth (aka, earthquake object or forecastsite object).
	#
	#loc = [0,0]		# [lon, lat]
	#eventDtm = None
	#eventftime = None # float version of time
	#mc=2.0
	#z=1.0
	z=None
	#
	#
	def __init__(self, loc=[0.0,0.0], dxdy=[.1, .1], evtime=1.0, mc=2.0):
		return self.initialize(loc=loc, dxdy=dxdy, evtime=evtime, mc=mc)
	
	def initialize(self, loc=[0.0,0.0], dxdy=[.1, .1], evtime=1.0, mc=2.0):
		# note scaling exponents are pre-declared.
		self.loc=loc
		self.latfactor = math.cos(2.0*math.pi*loc[1]/360.0)
		#
		for i in xrange(len(self.loc)):
			self.loc[i]=float(self.loc[i])
		# element size:
		if type(dxdy).__name__==type(.1).__name__:
			dxdy=[dxdy, dxdy]
		for i in xrange(len(dxdy)):
			dxdy[i]=float(dxdy[i])
		self.dxdy=dxdy	# note: these values are in lat/lon at this point.
		#
		self.mc=float(mc)
		#
		self.setEventTime(evtime)
		#
		#
		return None
	
	def getIntensity(self, equake):
		# dN/dt for this site. localIntensity() returns a value like n/[sec][km^2]. here, we itegrate that over area.
		dz=equake.localIntensity(inloc=self.loc, intime=self.eventftime)
		#
		# note this is an approximation in which dz/dA is calculated from circular coords. and now we interpolate using rect. coords.
		dZ = dz*self.dxdy[0]*self.dxdy[1]*self.latfactor # though a more accurate calc. can be made by doing the integral (where the
																			# the approximation dA = drdr,
																			# Ndot = (orate/2*pi*r0^(-(1+rho))) * (1/(rho*(rho-1)) * [r2^{1-rho} - r1^{1-rho}]
		return dZ
	
	def getIntensities(self, equakes):
		# for this individual site (presumably, for all sites but righ here for this site), calculate
		# contribution from all earthquakes.
		# note that this is not integrated -- z ~ n/[sec][km^2]
		#dzs=[]
		# add mpp...
		#pool = mpp.Pool()
		#dzs0 = []
		#
		# to accomplish MPPP:
		# 1) create a Process() from each equake (either subclass it, create a Process() that contains the object,
		#    or create a Process() that references that object's.localIntensity() function as its run-target. some,
		#    if not all, of these approaches may have problems pickling the various bits to the sub-processes,
		#    and it will be necessary to include Pipe()s to return the data.
		# 2) reorganize so that we explicitly pass both locations and times (for equake and site) to a
		# localIntensity() function. this can be done with a Pool() and apply_async().
		#
		# the Process() approach is probably fastest, so long as object pickling can be minimized.
		#
		return [equake.localIntensity(inloc=self.loc, intime=self.eventftime) for equake in equakes]
	
	def getz(self, equakes):
		#print 'setting z.'
		z = sum(self.getIntensities(equakes))
		#self.z=z
		return z
	
	def setz(self, equakes):
		#print 'setting z.'
		z = sum(self.getIntensities(equakes))
		self.z=z
		return z
	
class earthquake(locbase):
	# (probably should have made forecastsite() first and then derived earthquake from it).
	#
	# earthquake parameters (purely observable):
	loc = [0.0,0.0]		# [lon, lat]
	mag=None
	#eventDtm = None
	eventftime = None # float version of time
	#mc=2.0	# obviously, this is not an earthquake propertie, but we will need it to calculate rates, etc.
				# another approach might be to return self-rates, aka, the rate at which an earthquake produces
				# itself. this could be converted to catalog rates, of course, by subtracting mc.
	#
	# scaling exponents:
	#rho=1.37				# spatial scaling exponent
	rho = 1.5			# let's see what this looks like...
	#p = 1.05  			# Omori scaling exponent
	sig = 1.68			# aftershock surface density exponent.
	#alpha = 2.28		# rupture duration scaling exponent, nominally 
	b=1.0
	#dmstar=1.2
	dl0=1.0
	#drho=5.25
	
	# calculated earthquake properties:
	ltau = None  # Omori pram (log[tau]) in seconds
	lrupttime=None	# duration of rupture (log[rupture duration]) in seconds
	lt0 = None	# Omori pram derrived from rupt. duration (typically log(t0) = lrupttime + 1.0)
	lruptlen = None	# log rupt. len (in km).
	#
	def __init__(self, mag=5.0, loc=[0.0,0.0], evtime=1.0, mc=2.0, drho=5.25, p=1.05, dmstar=1.2, dm=1.0, alpha = 2.28, eqtheta=0.0, eqeps=1.0, rtype='ssim', p_map=None):
		return self.initialize(mag=mag, loc=loc, evtime=evtime, mc=mc, drho=drho, p=p, dmstar=dmstar, alpha=alpha, eqtheta=eqtheta, eqeps=eqeps, rtype=rtype, p_map=p_map)
	
	def initialize(self, mag=5.0, loc=[0.0,0.0], evtime=100.0, mc=2.0, drho=5.25, p=1.05, dmstar=1.2, dm=1.0, alpha = 2.28, eqtheta=0.0, eqeps=1.0, rtype='ssim', p_map=None):
		# note scaling exponents are pre-declared.
		self.drho=float(drho)
		if p==None: p=1.05
		self.p=float(p)
		self.p_map = p_map	# use this parameter to produce time-indepentent plots. aka, initialize ETAS with regular p~1.05, then use
							# p_map=0.0 to produce a time-independent (aka, accumulate all events) map.
		self.dmstar=dmstar	# mag difference for a primary sequence (a mainshock and its aftershocks)
		self.dm=dm				# mag difference for a full (compound) aftershock sequence
		self.mag=float(mag)
		self.loc=loc
		self.alpha=alpha
		for i in xrange(len(self.loc)):
			self.loc[i]=float(self.loc[i])
		#
		self.mc=float(mc)
		self.eqtheta=eqtheta
		self.eqeps=eqeps		# earthquake dipole moment (aka, direction and magnitude for elliptical transform).
		self.rtype=rtype
		#
		# note: evtime can be a datetime, days, seconds. eventually, it will convert to seconds. small times (in seconds)
		# will be interpreted as days and converted to seconds.
		#
		#print "initializing: ", evtime
		nothin=self.setEventTime(evtime)
		#print "initializing2: ", self.eventftime, self.eventDtm
		#
		# 13 oct 2014: comment
		#self.dmp=1.0-math.log10(self.p-1.0)/self.b	# scaling exponent for initial rate 1/tau?
		#
		# and of course, here is the debate. how do we calc. tau/t0? this method might work actually...
		# the main problem from before was that the catalog was not being converted from days -> seconds.
		# we do, however, introduce a new method of calculating t0/tau. rather than use the theoretical maximum
		# (GR) rupture rate, we use the event-resolution time (the time when the omori-interval is greater than
		# the rupture duration (or a factor thereof).
		#
		#self.ltau=self.ltauYoder(m=mag, alpha=self.alpha, p=self.p, b=self.b, dmp=self.dmp, mc=mc)
		#self.lt0 = self.getlt0(m=mag, alpha=self.alpha)
		self.lt0 = self.getlt02 (mstar=mag, alpha=self.alpha, dmstar=self.dmstar, drho=self.drho)
		self.ltau = self.getltau2(mstar=mag, mc=self.mc, alpha=self.alpha, p=self.p, drho=self.drho, thislt0=self.lt0)
		self.lrupttime=self.lDeltat(m=mag, alpha=self.alpha)
		self.lruptlen=self.getlruptLen(m=mag, dl0=None)
		#
		return None
	#
	# earthquake utility functions:
	def omoriRate(self, t, t0=None, tau=None, p=None):
		if t0==None: t0=10.0**self.lt0
		if tau==None: tau=10.0**self.ltau
		if p==None: p=self.p
		#
		rinv=tau*(t0+t)**p	# "interval" form of Omori
		tfact=1.0				# don't remember what this is for... some factor of the rate.
		#tfact=5.0
		return tfact/rinv
	#
	def omoriN(self, t1=0., t2=None, t0=None, tau=None, p=None):
		# integral of omoriRate (N occurred since t0) for p!=1.0 (p>1.0).
		if t0==None: t0=10.0**self.lt0
		if tau==None: tau=10.0**self.ltau
		if p==None: p=self.p
		if t2==None: t2=t1+t0	# or should we let this error off?
		#
		if p==1.0:
			return omoriNln(t=t, t0=t0, tau=tau)
		#
		N=( (t0 + t2)**(1.0-p) - (t0 + t1)**(1.0-p) )/(tau*(1.0-p))
		#
		return N
	#			
	def omoriNln(self, t1=0., t2=None, t0=None, tau=None):
		# integral of omoriRate (N occurred since t0) with p=1.0
		if t0==None: t0=10.0**self.lt0
		if tau==None: tau=10.0**self.ltau
		if t2==None: t2=t1+t0	# or should we let this error off?
		p=1.0
		#
		#N=(math.log((1.0 + t/t0))/tau	# note: log(x) = ln(x)
		N=(math.log((t2 + t0)/(t1 + t0)))/tau
		#
		return N
	#	
	def radialDensity(self, m=None, r=0.0, lrupt=None, rho=None, dm=None):
		# lrupt -> l_rupture (rupture length). l_rupture = 10**lruptlen
		if lrupt==None: 
			lrupt=10.0**self.lruptlen
			if dm==None: dm=self.dm
		#
		if m==None: m=self.mag
		if rho==None: rho=self.rho
		if dm==None: dm=self.dm
		Nom=10.0**(m-dm-self.mc)
		#
		#note: is this r=0 at the edge of the rupt. or r=0 at the center?
		# the (r+r0) formulation implies the "scaling" distance. we want just plain distance, and assume
		# the rate is constant for 0 < r < lrupt
		#
		#
		# omori-like formulation. r is the "scaling" distance; r=0 indicates the edge of the rupture.
		# the "scaling" distance r' = r - lrupt
		#rate = (rh0 - 1.0)*(lrupt**(rho-1.0))*Nom*(lrupt + r)**(-rho)
		#
		# equivalently, if r is the dist. from the epicenter,
		if r<lrupt: r=lrupt		# note: this will produce an upper-truncated distribution (which is good).
		rate = (rho-1.0)*(lrupt**(rho-1.0))*Nom*r**(-rho)
		# noting that scaling will start at r=lrupt, so effectively r = lrupt + r'
		#
		return rate	
	#
	'''	
	# i think this is left over from an older (attempted) model... and can be chucked.
	def dmprime(self, p=None, b=None):
		# scaling decay exponent for initial omori rate (?), 1/tau:
		if p==None: p=self.p
		if b==None: b=self.b
		dmp = 1.0-math.log10(p-1.0)/b	# usually dmp~2.0
		#
		return dmp
	'''
	#
	def getruptlen(self, m=None, dl0=None):
		if m==None: m=self.mag
		if dl0==None: dl0=self.dl0
		#
		return 10.0**self.getlruptLen(m,dl0)
	
	def getlruptLen(self, m=None, dl0=None):
		if dl0==None: dl0=self.dl0
		if m==None:
			m=self.mag
		return (m/2.0) - 1.755 - dl0

	#	method1: assuming some maximum rate at the end of rupture rate = N_GR/Delta-t_R * (some factor)
	def gett0(self, lt0=None):
		if lt0==None:
			lt0=self.lt0
		return 10**lt0
	
	def gettau(self, ltau=None):
		if ltau==None:
			ltau=self.ltau
		return 10**ltau
	#
	'''
	def ltauYoder(self, m=None, alpha=None, p=None, b=1.0, dmp=None, mc=None):
		if mc==None: mc=self.mc
		if p==None: p=self.p
		if dmp==None: dmp=self.dmp
		if m==None: m=self.mag
		if alpha==None: alpha=self.alpha
		if b==None: b=self.b
		if mc==None: mc=self.mc
		#
		# actualy log(tau)
		ltau=alpha*(p-1.0) + b*(dmp+mc) + m*((1.0-p)/2.0 - b)
		#
		return ltau
	'''
	#
	def lDeltat(self, m=None, alpha=None):
		if m==None: m=self.mag
		if alpha==None: alpha=self.alpha
		
		return (m/2.0) - alpha

	def getlt0(self, m=None, alpha=None):
		# experimentally, 10*deltaT seem to make nice t0 values. let's explicitly make the distinction between
		# these two parameters and see where it gets us...
		if m==None: m=self.mag
		if alpha==None: alpha=self.alpha
		#
		return self.lDeltat(m=m, alpha=alpha) + 1.0
	###
	##
	# method 2: upper rate limit is where events become resolvable (dt_omori > dt_rupt.)
	# ultimately, i think these two approaches both work (except it is better to estimate one
	# parameter from intrinsic/intensive physics (max rate, etc.) and then constrain the second
	# by integrating Omori.
	#
	# ... and of course, now all of this needs to be cleaned up a bit.
	def getlt02 (self, mstar=None, alpha=alpha_0, dmstar=None, drho=None):
		# the value drho=5.25 appears to be good within about +/- 0.5 for medium earthquakes.
		# ... but with respect to drho -- is this our new dmprime? i think so, more or less. note that we're discovering
		# that dmprime -> dmprime-mc (aka, the initial rate is not magnitude dependent. basically, if we decrease mc 
		# which nominally would increase the number of expected (observable) events, we increast the number of events lost in
		# in the early coda (or that might not be the right language -- we increase the number of events obscured by overlapping
		# larger events).
		if mstar==None: mstar = self.mag
		if alpha==None: alpha=self.alpha
		if dmstar==None: dmstar=self.dmstar
		if drho==None: drho=self.drho
		#
		x=(.5*mstar) - alpha - dmstar - 1.0 + drho
		return x

	def getltau2(self, mstar=None, mc=None, alpha=alpha_0, p=None, drho=None, thislt0=None):
		if mstar==None: mstar = self.mag
		if alpha==None: alpha=self.alpha
		if drho==None: drho=self.drho
		if p==None: p=self.p
		if thislt0==None:
			thislt0=getlt02(mstar=mstar, alpha=alpha, dmstar=1.2, drho=drho)
		x = mc-.5*mstar-alpha - p*thislt0 + drho
		return x
	
	def lruptLen(self, m, dl0=1.0):
		return (m/2.0) - 1.755 - dl0
	
	def localIntensity(self, inloc, intime, rho=None, coordtype='rad', thetaconv=0, abratio=1.0, eqtheta=None, eqeps=None, rtype=None):
		#
		# self-similar formulation.
		#
		# some hard-code valuse for now:
		#
		# eqtheta, eqeps are the local dipole moment of the aftershock distribution. nominally, this is an elliptical
		# convolution along the fault structure or is determined by fitting local seismicity to a line; the uncertainty
		# can be used as the eccentricity. epsilon=a*sigma, or something like that.
		if eqtheta==None: eqtheta = self.eqtheta	# let's take theta in degrees
		if eqeps==None: eqeps = self.eqeps
		#
		if rho==None: rho=self.rho	# radial density scaling exponent, which before too long we'll call "q".
		if rho==None: rho=1.5
		dmtau=0.0	# from HM, Toh., Pf., ElM, Coa., SanSim <dmtau>~7, but for rates 0 is better.
		dtau = 2.28
		dmstar=1.0
		D=1.5	# fractal dimension	# from paper, <D>~1.35, but that seems too low. the best fits are D~1.5.
		
		#self.p=1.05
		p=self.p
		#rho=1.5
		q=rho
		#
		# rtype can be 'ssim', 'omorisat', 'satpl' (see below)
		if rtype==None:
			rtype=self.rtype
			#print "settting rtype=",self.rtype
		if rtype==None:
			rtype='ssim'
		#
		#rtype='omorisat'	# saturated omori-type (aka, saturated then Omori-like with m'=m-dm) "aftershock formulation"
		#rtype='ssim'		# fully self-sim omori-type														"self-similar formulation"
		#rtype='satpl'		# saturated, then PL
		#
		# return estimated current rate of seismicity based on Omori and spatial decay.
		#
		# specifically, report dN/(dt * dA), or d(omori rate)/dA.
		#
		# this is a little bit tricky with the modified Omoir type formats.
		# the easiest way is probably to get the first (Omori) intensity directly.
		# then, calculate the relative spatial intensity by comparing N'(R) to N'(R->l_r).
		# otherwise, we have to solve for the spatial (or Omori) parameters based on the 
		# decayed rate (aka, new t0, tau parameters).
		# also, assume the decay rate of the full (compound) sequence is p=1.0 (not 1.05),
		# and something similar for spatial decay (though the 1.37 numbers from Frolich probably
		# represent the full sequence.
		# do we integrate over the area here, or is that done by the calling function? (which
		# probably assumes linearity across the bin).
		t=self.timesince(intime)
		if t<0.0: return 0.0	# earthquake has not yet happened.
		#
		R=self.dispFrom(inloc=inloc, coordtype='rad')	# displacement...
		#r=R[0]	# for coordtype='rad', we get [R, theta]
		
		# g1=ggp.WGS84.Inverse(lambs0, phis0, rw[1], rw[2])
		g1=ggp.WGS84.Inverse(self.loc[1], self.loc[0], inloc[1], inloc[0])
		r=g1['s12']/1000.0
		
		
		lNom = self.b*(self.mag - self.dm - self.mc)
		Nom = 10.0**lNom
		#
		# for a strict PL type model, use something like this: 
		'''
		tcrit = 10.0**(1.0+self.lrupttime)	
		#
		if t<tcrit:
			#t=10.0**self.lt0	# flat bit of the Omori distribution.
			t = 0.0
		else:
			t=t-tcrit
		'''
		lt0ssim = 7.0*self.mag/6.0 - (2.0*self.mc/3.0) - dtau - dmstar - (2./3.)*math.log10(1.5) + (dmtau/3.0) + math.log10(p-1.0)
		lr0ssim = self.mag*((6.0+D)/(4.0+2*D)) - (2.0/(2.0+D))*(self.dm + self.mc - math.log10((2.0+D)/2.0)) + math.log10(q-1.0) - 1.76 - math.log10(2.0)
		r0ssim = 10**lr0ssim
		#
		t0ssim = 10.**lt0ssim
		taussim = (t0ssim**(1.0-p))/(Nom*(p-1.0))
		#
		lt0 = lt0ssim
		t0 = t0ssim
		tau = taussim
		#
		# spatial parts (this is sort of the long way, but computers are fast):
		lrupture=10**self.lruptlen		# now this is not actually the rupture-length, but we used to think it (more or less) was,
												# and it makes a nice parameter for "ssimbump"
		#
		# self-similar formulation:
		#if self.mag>7.6: print "bigmag: %f" % self.mag
		lNas = (2.0*(2.0+D))*math.log10(1.+D/2.) + (D/(2.+D))*(self.mag - self.dm - self.mc)
		lNprimemax = lNas - lrupture + math.log10(2.0)	# log(2) from: 0 <r < L/2 -- assume mainshock is centered.
		Nas=10.0**lNas
		Nprimemax = 10.0**lNprimemax
		r0 = Nom*(q-1.0)/Nprimemax
		chi = (r0**(1.0-q))/(Nom*(q-1.0))
		#
		#
		# try (eventually) convolving r with the local fault moment. for now, just use an angle and the ratio a/b of the major/minor axes.
		if eqtheta!=None and eqeps!=None:
			#print eqtheta*deg2rad, R[1]
			thetaconv=eqtheta*deg2rad - R[1]
			#abratio = 3.0
			abratio = eqeps	# or maybe 1/eqeps
			#yfact = abratio*r*math.sin(R[1]-thetaconv)
			#xfact = r*math.sin(R[1]-thetaconv)/abratio
			xprime = r * math.cos(thetaconv)
			yprime = r * math.sin(thetaconv)
			rprime = ((abratio*yprime)**2 + (xprime/abratio)**2)**.5
			#print eqtheta, R[1], thetaconv, r, rprime, r/rprime
			r=rprime
			#
			coordtype='rect'
			if coordtype=='rect2':
				R=self.dispFrom(inloc=inloc, coordtype='rect')	# displacement...
				#
				#thistheta=math.atan(R[1]/R[0])
				#print thistheta, 40.0*deg2rad, thistheta-40.0*deg2rad, abs(math.cos(thistheta-40.0*deg2rad))
				# for rect. coords:
				#r=math.sqrt(R[0]*R[0] + R[1]*R[1])
				#
				# trying a convolution...
				#thistheta = math.acos(R[0]/r) - thetaconv
				xprime = R[0]*math.cos(thetaconv) - R[1]*math.sin(thetaconv)
				yprime = R[0]*math.sin(thetaconv) + R[1]*math.cos(thetaconv)
				#
				r = math.sqrt(abratio*yprime*yprime + xprime*xprime/abratio) # this gives an equal-area transformation.
			
		#
		# this would be for explicit saturation...
		#if r<lrupture:
		#	r=lrupture		# this implies that the spatial-rate distribution saturates at (some factor of) lrupture
		#
		# the basic idea is to
		#1) calculate the absolute Omori or spatial rate/density
		#2) then, calculate the expected decay in the other.
		# aka, total omori rate at t=t', then the relative decay in the spatial dimension
		#
		# Omori rate: full rate of aftershocks.
		#orate=self.omoriRate(t=t, t0=self.gett0(), tau=self.gettau(), p=1.0)
		#orate=self.omoriRate(t=t, t0=10**self.lt0, tau=10**self.ltau, p=self.p)
		#
		# 13 oct 2014:
		#self.p_map=0.0
		this_p = self.p
		try:
			if self.p_map != None: this_p = self.p_map
		except:
			# there is no special "map p" defined. use regular Omori p.
			pass
		#
		#print "this_p: ", this_p
		orate = 1.0/(tau * (t0 + t)**this_p)
		#orate = 1.0/(tau * (t0 + t)**self.p)
		#
		# and normalizing the radial distribution, we get:
		# (note, this uses a straight omori-like distribution).
		# N'/N'_0 = (r0/(r0 + r))**q
		#spatialdensity = ((orate * (rho-1.0)*lrupture**(rho-1.0)))/( 2.0*math.pi*r**(1.0+rho))
		#radialDens = self.radialDensity(r=r, lrupt=None, rho=None, dm=None)	# all events at distance r
		#
		#Rrad=self.dispFrom(inloc=inloc, coordtype='rad')
		#r=Rrad[0]
		rprime=r
		#if rdenom<lrupture:
			#rdenom=lrupture
		#
		# remember, these have to be normalized correctly so that integral(f(x)) -> N_omori or 1.0.
		# ... and right now, the normalization is not right, but it is probably good enough for demonstrative purposes.
		# self-similar omori-like:
		rcrit = .5*lrupture	# this is the older "core" formulation (i think) log(thisL) = m/2 -2.76, and maybe a 1/2 
									# as well.
		rcrit = .5*10.0**(.5*self.mag - 1.76)	# this one is worth investigating, particularly for large earthquakes
														# 
		if rtype=='ssim':
			# 
			#radialDens = (r0/(r0+r))**q	#where did this come from? to be fair, it does make a nice picture.
			#i think the correct formulation is:
			#
			#radialDens = (q-1.0)*(r0**(q-1.0))/((r0+rprime)**q)
			
			#radialDens = (q-1.0)*(r0ssim**(q-1.0))*(r0ssim + rprime)**(-q)
			radialDens = (q-1.0)*(r0ssim**(q-1.0))*(r0ssim + rprime)**(-q)
		#
		if rtype=='ssimbump':
			# a hypothetical distribution accounting for Tahir2012, ssim with reduced probability
			# inside the rupture.
			# *math.exp((bumpfact*(x-t0))/(x + t0))
			#print "ssimbumping..."
			ssim_BumpFactor=3.0
			ssim_rimfact = 1.25
			r0prime = ssim_rimfact*lrupture	# lrupture ~ Lr/10
			#
			#bumpPart = math.exp(ssim_BumpFactor*(rprime-r0prime)/(rprime+r0prime))
			bumpPart = math.exp(ssim_BumpFactor*(rprime-r0prime)/(rprime+.15*r0prime))
			radialDens1 = ((q-1.0)*(r0ssim**(q-1.0))*(r0ssim + rprime)**(-q))*bumpPart
			
			#radialDens2 = .5*(q-1.0)*(r0ssim**(q-1.0))*(r0ssim + rprime)**(-q)
			radialDens0 = .5*(q-1.0)*(r0ssim**(q-1.0))*(r0ssim + 0.0)**(-q)
			
			radialDens = radialDens1 + .15*radialDens0
			
			#radialDens = radialDens1 + radialDens2
			#radialDens = ((q-1.0)*(r0ssim**(q-1.0))*(r0ssim + rprime)**(-q))*(math.exp(ssim_BumpFactor*(rprime-rcrit)/(rprime)))
		#
		if rtype=='omorisat':
			#"saturated omori-type"
			#
			# "Aftershock" decay formulation
			# 1) inside rupture, assume full saturation.
			# 2) outside, assume omoriX-like decay but with r0 determined from largest aftershock; chi effectively
			#    determined by the boundary condition that f_inside = f_outside at rupture boundary. note r begins at the 
			# 	  rupture boundary.
			# note that the calculation of chi2 (for the r>L_r part) might be a bit messed up, since it's based on
			# normalization... of course, here we're not using any of those formulations. we can normalize by setting
			# N'(r_crit) = N'0, but there is some question as to whether this is equivalent to the initial analysis
			# (where chi is determined from ssim-framework). i guess, so long as decay rate r0 is the same, this should
			# be ok; chi is an initial rate.
			#
			# normalizing: 1 = N'0*rcrit + int( [1/chi]
			# so chi = 1/(N'0*r**q)
			# and: N0' = (q-1)/(r_crit*(q-1) + r_2)
			r02=.5*10.0**((self.mag-1.0)/2. - 1.76)
			normfact = (q-1.0)/(rcrit*(q-1.0) + r02)
			chi2inv = (r02**(q))*normfact
			#
			if r<rcrit:
				#radialDens=1.0
				#radialDens = 1.0*normfact
				
				#rprime=0.0
				#radialDens = chi2inv/((r02 + rprime)**q)
				r=rcrit*1.001	# this factor of 1.001 helps to mitigate some serious anomalies related to contouring.
				
				#radialDens = chi2inv/((r02 + 0.0)**q)
				
				#print "rd=%f" % radialDens
				#
			if r>=rcrit:
				rprime = r-rcrit
				#chi2 = (r02**(-q))/Nprimemax
				#
				radialDens = chi2inv/((r02 + rprime)**q)
				
				#radialDens=(r02/(r02+rprime))**q		# this normalization also makes for a nice map.
				
				#radialDens = (q-1.0)*(r02**(q-1.0))/((r02+rprime)**q)
				#
		if rtype=='omorisatbump':
			#"saturated omori-type"
			# (see also tahir type below).
			#
			# "Aftershock" decay formulation (with tahir2012, zalohar2013 largest-aftershock formulation)...
			#
			# 1) inside rupture, assume full saturation.
			# 2) outside, assume omoriX-like decay but with r0 determined from largest aftershock; chi effectively
			#    determined by the boundary condition that f_inside = f_outside at rupture boundary. note r begins at the 
			# 	  rupture boundary.
			# note that the calculation of chi2 (for the r>L_r part) might be a bit messed up, since it's based on
			# normalization... of course, here we're not using any of those formulations. we can normalize by setting
			# N'(r_crit) = N'0, but there is some question as to whether this is equivalent to the initial analysis
			# (where chi is determined from ssim-framework). i guess, so long as decay rate r0 is the same, this should
			# be ok; chi is an initial rate.
			#
			# normalizing: 1 = N'0*rcrit + int( [1/chi]
			# so chi = 1/(N'0*r**q)
			# and: N0' = (q-1)/(r_crit*(q-1) + r_2)
			r02=.5*10.0**((self.mag-1.0)/2. - 1.76)
			normfact = (q-1.0)/(rcrit*(q-1.0) + r02)
			
			ssim_BumpFactor=3.0
			ssim_rimfact = 2.0
			r0prime = ssim_rimfact*r02	# note that the normalization for this is screwed up, so it's only for demonstration.
			
			#
			if r<rcrit:
				#radialDens=1.0
				#radialDens = 1.0*normfact
				#print "here?"
				radialDens = normfact*(r/rcrit)
				r=r*1.001
				#if r<(.1*rcrit): radialDens = normfact*(.1)**.5
				#
			if r>=rcrit:
				rprime = r-rcrit
				#chi2 = (r02**(-q))/Nprimemax
				chi2inv = (r02**(q))*normfact
				#
				radialDens = chi2inv/((r02 + rprime)**q)
				
				#radialDens=(r02/(r02+rprime))**q		# this normalization also makes for a nice map.
				
				#radialDens = (q-1.0)*(r02**(q-1.0))/((r02+rprime)**q)
				#
		#
		#
		if rtype in ('tz', 'tahir'):
			# tahir-zalohar type distribution, where hazard is focused on largest aftershock, presumably
			# at a rupture length distance (for now. later versions might be more sophisticated).
			#
			# start with an omori-like distribution, combine with poisson dist. for now, assume k=1, l_r.
			omori_fact   = .25
			poisson_fact = .75
			op_sum = omori_fact + poisson_fact
			#
			# normalize these factors (which we will one day parameterize)
			omori_fact/=op_sum
			poisson_fact/=op_sum
			#
			radialDens_omori   = (q-1.0)*(r0ssim**(q-1.0))*(r0ssim + rprime)**(-q)
			#poisson_mean_r = r0ssim		# or maybe rupture len?
			poisson_mean_r = .5*lrupture	#
			#radialDens_poisson = (1.0/poisson_mean_r)*(rprime/poisson_mean_r)*math.exp(-rprime/poisson_mean_r)
			radialDens_poisson = (rprime/poisson_mean_r)*math.exp(-rprime/poisson_mean_r)
			#
			radialDens = omori_fact*radialDens_omori + poisson_fact*radialDens_poisson
			#
		#		
		#
		if rtype=='satpl':
			# fully saturated PL.
			# the normalization can be done fully explicitly fron N = N'0*r0 + N'0*r0^q*int(r^-q dr)
			# for the simplified, relative intensity: N->1 and we solve for N'0
			# N'0 = (q-1)/(r0*q)
			normfact = (q-1.0)/(rcrit*q)
			if r<rcrit:
				radialDens=1.0*normfact
				#radialDens=Nprimemax/Nom
			if r>=rcrit:
				#rprime = r-rcrit
				#radialDens=(Nprimemax/Nom)*(rcrit**q)*rprime**(-q)
				#radialDens = (rcrit**q)*1.0*rprime**(-q)
				radialDens = (normfact*rcrit**q)*rprime**(-q)		# note same as above but normfact != 1
		#
			
			#
		# avoid divide by zero...
		#if r==0.0: return 0.0
		if rprime==0.0: return 0.0
		#
		spatialdensity = radialDens/(2.0*math.pi*rprime)
		#
		#localIntensity = orate*spatialDecayFactor
		localIntensity = spatialdensity*orate
		
		#
		# so this is the total rate of aftershock production.
		# now, calculate the spatial distribuiton. the total number of remaining (expected) aftershocks is:
		#
		#lruptTime = self.lDeltat(m, alpha)
		
		#return 1.0+lruptTime	
		
		#if self.mag>7.0: print "m=%f, dI=%f (%f)" % (self.mag, math.log10(localIntensity), self.eventftime/(3600.24))
		
		return localIntensity

	def localIntensityLRUPT(self, inloc, intime, rho=None, coordtype='rad', thetaconv=0, abratio=1.0, eqtheta=None, eqeps=None):
		#
		# eqtheta, eqeps are the local dipole moment of the aftershock distribution. nominally, this is an elliptical
		# convolution along the fault structure or is determined by fitting local seismicity to a line; the uncertainty
		# can be used as the eccentricity. epsilon=a*sigma, or something like that.
		if eqtheta==None: eqtheta = self.eqtheta	# let's take theta in degrees
		if eqeps==None: eqeps = self.eqeps
		#
		if rho==None: rho=self.rho	# radial density scaling exponent, which before too long we'll call "q".
		#
		# this is the real thing. return estimated current rate of seismicity based on
		# Omori and spatial decay.
		#
		# specifically, report dN/(dt * dA), or d(omori rate)/dA.
		#
		# this is a little bit tricky with the modified Omoir type formats.
		# the easiest way is probably to get the first (Omori) intensity directly.
		# then, calculate the relative spatial intensity by comparing N'(R) to N'(R->l_r).
		# otherwise, we have to solve for the spatial (or Omori) parameters based on the 
		# decayed rate (aka, new t0, tau parameters).
		# also, assume the decay rate of the full (compound) sequence is p=1.0 (not 1.05),
		# and something similar for spatial decay (though the 1.37 numbers from Frolich probably
		# represent the full sequence.
		# do we integrate over the area here, or is that done by the calling function? (which
		# probably assumes linearity across the bin).
		t=self.timesince(intime)
		if t<0.0: return 0.0	# earthquake has not yet happened.
		#
		#if t<10.0**self.lt0:
		#
		# t is the actual time since the (begining of?) the earthquake. convert to t', the "scaling time".
		# 
		tcrit = 10.0**(1.0+self.lrupttime)	# but i think i just made up this number. let's try something more
														# scientific. how'bout dt_omori(t=0) = tau * t0**p
														# that said, this does not seem to be terribly important, except as 
														# an mechanism to handle the singularity.
		#tcrit = 10**(self.drho+self.mc-self.alpha-self.mag/2.0)
		#print "tcrit: %f, %f" % (tcrit/days2secs, (10.0**(1.0+self.lrupttime))/days2secs)
		if t<tcrit:
			#t=10.0**self.lt0	# flat bit of the Omori distribution.
			t = 0.0
		else:
			t=t-tcrit
		#
		# print "t: %f" % t
		#
		lrupture=10**self.lruptlen
		#
		# rect-lin approximation. note, a minor error gives the best results here, specifically when we (accidentally)
		# omit the latitude displacement (R = dx, aka R=rvect[0] -- this stemmed from a mistake where it was mistakenly
		# interpretedt that dispFrom -> [R,theta] always, not just for the polar-coords option).
		R=self.dispFrom(inloc=inloc, coordtype='rad')	# displacement...
		#r=self.sphericalDist(inloc=inloc)
		r=R[0]	# for coordtype='rad', we get [R, theta]
		#theta=R[1]
		#
		# try convolving r with the local fault moment. for now, just use an angle and the ratio a/b of the major/minor axes.
		if eqtheta!=None and eqeps!=None:
			#print eqtheta*deg2rad, R[1]
			thetaconv=eqtheta*deg2rad - R[1]
			#abratio = 3.0
			abratio = eqeps	# or maybe 1/eqeps
			#yfact = abratio*r*math.sin(R[1]-thetaconv)
			#xfact = r*math.sin(R[1]-thetaconv)/abratio
			xprime = r * math.cos(thetaconv)
			yprime = r * math.sin(thetaconv)
			rprime = ((abratio*yprime)**2 + (xprime/abratio)**2)**.5
			#print eqtheta, R[1], thetaconv, r, rprime, r/rprime
			r=rprime
			#
			coordtype='rect'
			if coordtype=='rect2':
				R=self.dispFrom(inloc=inloc, coordtype='rect')	# displacement...
				#
				#thistheta=math.atan(R[1]/R[0])
				#print thistheta, 40.0*deg2rad, thistheta-40.0*deg2rad, abs(math.cos(thistheta-40.0*deg2rad))
				# for rect. coords:
				#r=math.sqrt(R[0]*R[0] + R[1]*R[1])
				#
				# trying a convolution...
				#thistheta = math.acos(R[0]/r) - thetaconv
				xprime = R[0]*math.cos(thetaconv) - R[1]*math.sin(thetaconv)
				yprime = R[0]*math.sin(thetaconv) + R[1]*math.cos(thetaconv)
				#
				r = math.sqrt(abratio*yprime*yprime + xprime*xprime/abratio) # this gives an equal-area transformation.
			
		#
		#r/=abs(math.cos(thistheta + math.pi/4.0))	# for fun: this gives a cool dipole picture which is not (i assume) terribly physical.
		#
		# or use the proper spherical coords, directly:
		#
		#if r<10.0**(self.lruptlen):
		#	r=10.0**(self.lruptlen)
		if r<lrupture:
			r=lrupture		# this implies that the spatial-rate distribution saturates at (some factor of) lrupture
		#
		# the basic idea is to
		#1) calculate the absolute Omori or spatial rate/density
		#2) then, calculate the expected decay in the other.
		# aka, total omori rate at t=t', then the relative decay in the spatial dimension
		# 
	
		#orate=self.omoriRate(t=t, t0=self.gett0(), tau=self.gettau(), p=1.0)
		orate=self.omoriRate(t=t, t0=10**self.lt0, tau=10**self.ltau, p=self.p)
		#
		# note: the omori rate is being estimated directly from the magnitude
		# (by estimating the parameters tau, t0).
		#Nremaining=(10**(self.mag - 1.0 - self.mc)) - self.omoriN(t=t, t0=self.t0, tau=self.tau, p=self.p)
		#Nomori = 10.0**(self.mag - 1.0 - self.mc)
		#Nremaining = Nomori - self.omoriNln(t=t, t0=10.0**self.lt0, tau=self.tau)
		# and we use that as N in the spatial density formula.
		# the relative intensity factor will be
		#intensityFactor=Nremaining/Nomori
		#lrupture=10**self.lruptlen
		#spatialDecayFactor = (1.0 + r/lrupture)**(1.0-self.rho)
		#spatialDecayFactor = (1.0 + r/lrupture)**(-self.rho)		# this is the relative linear density (decay) at r. "1+" comes from using the CDF version.
																					# but it is probably bettter to just calc. the density at r, then integrate around the area.
																					# in fact, it is arguably the job of forecastlocation() to do that. earthquake just reports
																					# the local intensity at r
		#start with linear density (Felzer):
		# lambda = N/l_r = c * 1/(r0+r)**rho
		# integrate and set equal to Omori:
		# N_omori = integral [0 to inf] [ lambda dr ] = (c/(rho-1)) * r0**(1-rho)
		# solve for c:
		# c = N_omori * (rho-1) * r0**(rho-1)
		# r0 = l_rupt.
		# now, substitute omori-rate for N_omori to give a linear rate. to find the local rate at some radial distance
		# and some radial location theta, divide the total rate by 2pi*r
		#
		#spatialdensity = (orate/(2.0*math.pi*(lrupture**(-self.rho)))) * r**(-(1.0 + self.rho))
		#
		#spatialdensity = ((orate * (rho-1.0)*lrupture**(rho-1.0)))/( 2.0*math.pi*r**(1.0+rho))
		radialDens = self.radialDensity(r=r, lrupt=None, rho=None, dm=None)	# all events at distance r
		spatialdensity = radialDens/(2.0*math.pi*r)
		# noting that the 1+rho exponent comes from 1/r**rho * 1/2pi*r
		# also note that we are using r instead of (r0+r'). we correct for this 
		#
		#localIntensity = orate*spatialDecayFactor
		localIntensity = spatialdensity*orate
		
		#
		# so this is the total rate of aftershock production.
		# now, calculate the spatial distribuiton. the total number of remaining (expected) aftershocks is:
		#
		#lruptTime = self.lDeltat(m, alpha)
		
		#return 1.0+lruptTime			
		
		return localIntensity


def plotcolor(z, Z):
	# z is value, Z is max value
	#
	colorRange=int('0xFFFFFF', 0)
	#
	colorint=int(colorRange*z/Z)
	#colorhex=self.fillHexString(hex(colorint), 6)
	colorhex=fillHexString(hex(colorint), 6)
	colorhexstr='#' + colorhex.split('x')[1]
	#
	return colorhexstr

def fillHexString(hexstring='', strlen=6):
	strs=hexstring.split('x')
	while len(strs[1])<strlen:
		strs[1]='0'+strs[1]
	return strs[0] + 'x' + strs[1]

def getMFETAScatFromANSS(lons=[-121.0, -114.0], lats=[31.0, 37.0], dates=[None, None], mc=4.0):
	if dates==None: dates=[None, None]
	if dates[1]==None: dates[1]=dtm.datetime.now(pytz.timezone('UTC'))
	if dates[0]==None or dates[0]>=dates[1]: dates[0]=dates[1]-dtm.timedelta(days=840)
	#
	#print dates
	clist1=atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=dates, Nmax=999999, fout=None)
	catalog=[]
	#X,Y,M = [], [], []
	#
	for rw in clist1:
		catalog+=[[mpd.date2num(rw[0]), rw[1], rw[2], rw[3], rw[4]]]
		#catalog+=[[rw[0], rw[1], rw[2], rw[3], rw[4]]]
		#catalog+=[[mpd.datestr2num(rw[0]), rw[1], rw[2], rw[3], rw[4]]]
		#X+=[rw[2]]	# lon
		#Y+=[rw[1]]	# lat
		#M+=[rw[4]]
	return catalog

def checkPoly(inpoly, fignum=None):
	# polygons are getting mangled by the "find closest neighbor" algarithm.
	# try looking for dangling elements after a polygon has closed. we'll then either
	# discard them or insert them next to their closest neighbor.
	#
	# given local syntax, inpoly -> [[x], [y]]
	#
	x0=inpoly[0][0]
	y0=inpoly[1][0]
	outpoly=[[[inpoly[0][0]], [inpoly[1][0] ]]]
	#outpoly=[[[], []]]
	for i in xrange(1,len(inpoly[0])):
		thisx=inpoly[0][i]
		thisy=inpoly[1][i]
		#
		outpoly[-1][0]+=[thisx]
		outpoly[-1][1]+=[thisy]
		#if thisx==x0 and thisy==y0 and len(outpoly[-1][0])>1:
		if thisx==outpoly[-1][0][0] and thisy==outpoly[-1][1][0] and len(outpoly[-1][0])>1:
			print "(improperly) closed poly, %d/%d" % (i, len(inpoly[0])-1)
			outpoly+=[[[], []]]
	
	if len(outpoly[-1][0])<2: outpoly.pop()
	if outpoly[-1][0][0]!=outpoly[-1][0][-1] and outpoly[-1][1][0]!=outpoly[-1][1][-1]: outpoly.pop()
	#
	#plt.plot(inpoly[0], inpoly[1], 'o-', lw=3, alpha=.4)
	if fignum!=None:
		plt.figure(fignum)
		plt.clf()
	
		for ply in outpoly:
			plt.plot(ply[0], ply[1], '.--', alpha=.8)

	return outpoly

