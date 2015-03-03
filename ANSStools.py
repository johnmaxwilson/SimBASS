import datetime as dtm
import pytz
import calendar
import operator
import urllib
#
# note on datetimes:
# timezone awareness is confusing but not as bad as it looked a minute ago.
# datetime.utcnow() gives a UTC time, but without tz awareness.
# probably the better approach is to use datetime.now(pytz.timezone('UTC'))
tzutc=pytz.timezone('UTC')

def anssDateStr(x=dtm.datetime.now(pytz.timezone('UTC'))):
	yr=x.year
	mo=x.month
	dy=x.day
	hr=x.hour
	mn=x.minute
	sc=x.second
	ms=x.microsecond
	fsecs=float(sc) + float(ms)*(10**(-6.0))
	#
	# for now, assume the dates are properly aggregated -- no minutes=60, etc.
	# string-fixing is probably not necessary.
	#mo = ('00' + str(mo))[-2:]
	#dy = ('00' + str(dy))[-2:]
	#hr = ('00' + str(hr))[-2:]
	#mn = ('00' + str(mn))[-2:]
	#sc = ('00' + str(sc))[-2:]
	#
	#return '%d/%s/%s,%s:%s:%f' % (yr, mo, dy, hr, mn, secs)
	return '%d/%d/%d,%d:%d:%f' % (yr, mo, dy, hr, mn, fsecs)
	
#def getANSStoFilehandler(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.date(2001,01,01), dtm.date(2010, 12, 31)], Nmax=999999):
def getANSStoFilehandler(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.datetime(2001,01,01, tzinfo=tzutc), dtm.datetime(2010, 12, 31, tzinfo=tzutc)], Nmax=999999):
	# fetch data from ANSS; return a file handler.
	#
	# use urllib in "post" mode. an example from http://www.python.org/doc/current/library/urllib.html#urllib.FancyURLopener)
	# using "get" (aka, query-string method; note the ?%s string at the end of the URL, this is a single pram call to .urlopen):
	#
	#>>> import urllib
	#>>> params = urllib.urlencode({'spam': 1, 'eggs': 2, 'bacon': 0})
	#>>> f = urllib.urlopen("http://www.musi-cal.com/cgi-bin/query?%s" % params)
	#>>> print f.read()
	#
	# using "post" (note this is a 2 pram call):
	#>>> import urllib
	#>>> params = urllib.urlencode({'spam': 1, 'eggs': 2, 'bacon': 0})
	#>>> f = urllib.urlopen("http://www.musi-cal.com/cgi-bin/query", params)
	#>>> print f.read()
	#
	# make ANSS prams dictionary (thank james for the bash-template):
	# ANSSquery has day-resolution:
	# revision: ANSS has time resolution, but you have to replace "-" -> "/" and the " " (space) -> ","
	#dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	dates=dates0
	datestr1 = anssDateStr(dates[0])
	datestr2 = anssDateStr(dates[1])
	print datestr1, datestr2
	#
	#anssPrams={'format':'cnss', 'output':'readable', 'mintime':str(dates[0]).replace('-', '/'), 'maxtime':str(dates[1]).replace('-', '/'), 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':'E', 'searchlimit':Nmax}
	# so this is better, but i think it is still limited to 1 second resolution.
	anssPrams={'format':'cnss', 'output':'readable', 'mintime':datestr1, 'maxtime':datestr2, 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':'E', 'searchlimit':Nmax}
	f = urllib.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.urlencode(anssPrams))
	#
	# we might return f, a string of f, or maybe a list of lines from f. we'll work that out shortly...
	return f

#def catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.date(2005,01,01), None], Nmax=999999, fout='cats/mycat.cat'):
def catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.datetime(2005,01,01, tzinfo=tzutc), None], Nmax=999999, fout='cats/mycat.cat'):
	# get a basic catalog. then, we'll do a poly-subcat. we need a consistent catalog.
	# eventually, cut up "japancatfromANSS()", etc. to call this base function and move to yodapy.
	#
	# note: there may be a version inconsisency. older versions of this function may have returned catlist raw, in which
	# [..., depth, mag], where we regurn [..., mag, depth] here.
	#
	if dates0[1]==None:
		# i think this needs a "date" object, and datetime breaks.
		# so, make a Now() for date.
		#nowdtm=dtm.datetime.now()
		#dates0[1]=dtm.date(nowdtm.year, nowdtm.month, nowdtm.day)
		dates0[1]=dtm.datetime.now(tzutc)
	#	
	catlist=getANSSlist(lon, lat, minMag, dates0, Nmax, None)
	if fout==None: print " no file."
	
	if fout!=None:
		f=open(fout, 'w')
		f.write("#anss catalog\n")
		f.write("#lon=%s\tlat=%s\tm0=%f\tdates=%s\n" % (str(lon), str(lat), minMag, str(dates0)))
	
	rlist=[]
	for rw in catlist:
		# simple, right? except that ANSS has a habit of writing useless date-times like "2001/10/08 24:00:07.62" (hour=24), or
		# where minute=60. we could toss these. for now, assume 2001/10/8 24:00:00 -> 2001/10/9/00:00:00. change by proper time-arithmetic.
		# first, parse the date-string:
		strDt, strTm=rw[0].split()[0], rw[0].split()[1]
		if '/' in strDt: delim='/'
		if '-' in strDt: delim='-'
		strDts=strDt.split(delim)
		strTms=strTm.split(':')
		yr=int(strDts[0])
		mnth=int(strDts[1])
		dy=int(strDts[2])
		hr=int(strTms[0])
		mn=int(strTms[1])
		sc=float(strTms[2])
		microsecs=(10**6)*sc%1.
		# one approach is to start with year, month and add all the subsequent quantities using datetime.timedelta objects, which we have to
		# do once we get into callendar addition anyway...
		#so let's assume the date part is correct:
		myDt=dtm.datetime(yr, mnth, dy, tzinfo=tzutc)
		#mytimedelta=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(minutes=mn)
		myDt+=dtm.timedelta(seconds=sc)
		myDt+=dtm.timedelta(microseconds=microsecs)
		#
		rlist +=[[myDt, float(rw[1]), float(rw[2]), float(rw[4]), float(rw[3])]]
		if fout!=None:
			myDtStr='%d/%d/%d %d:%d:%d.%d' % (myDt.year, myDt.month, myDt.day, myDt.hour, myDt.minute, myDt.second, myDt.microsecond)	
			#
			#f.write('%s\t%s\t%s\t%s\n' % (rw[0], rw[1], rw[2], rw[4]))
			#f.write('%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4]))
			f.write('%s\t%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4], rw[3]))	# indlude depth...
	if fout!=None:
		f.close()
	
	 
	#return catlist
	return rlist

def getANSSlist(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.datetime(2001,01,01, tzinfo=tzutc), dtm.datetime(2010, 12, 31, tzinfo=tzutc)], Nmax=999999, fin=None):
	#
	# note: this appears to be a bad idea for global downloads. a full catalog is ~4GB, which kills my computer.
	#
	# note: this may be repeated exactly in ygmapbits.py
	# fetch new ANSS data; return a python list object of the data.
	# fin: data file handler. if this is None, then get one from ANSS.
	#dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	dates=dates0	# date/datetime issue is fixed.(ish)
	anssList=[]
	if fin==None:
		#print "get data from ANSS...(%s, %s, %s, %s, %s)" % (lon, lat, minMag, dates, Nmax)
		fin = getANSStoFilehandler(lon, lat, minMag, dates, Nmax)
		#fin = getANSStoFilehandler([-180, 180], [-90, 90], 0, [datetime.date(1910,01,01), datetime.date(2010, 01, 16)], 9999999)

		print "data handle fetched..."
		
	for rw in fin:
		if rw[0] in ["#", "<"] or rw[0:4] in ["Date", "date", "DATE", "----"]:
			#print "skip a row... %s " % rw[0:10]
			continue
		#anssList+=[rw[:-1]]
		# data are fixed width delimited
		# return date-time, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID (because those are all the available bits)
		#print "skip a row... %s " % rw
		rwEvdt=rw[0:22].strip()
		rwLat=rw[23:31].strip()
		if rwLat=='' or isnumeric(str(rwLat))==False or rwLat==None:
			continue
			#rwLat=0.0
		else:
			rwLat=float(rwLat)
		rwLon=rw[32:41].strip()
		if rwLon=='' or isnumeric(str(rwLon))==False or rwLon==None:
			#rwLon=0.0
			continue
		else:
			rwLon=float(rwLon)
		rwDepth=rw[42:48].strip()
		if rwDepth=='' or isnumeric(str(rwDepth))==False or rwDepth==None or str(rwDepth).upper() in ['NONE', 'NULL']:
			#rwDepth=0.0
			rwDepth=None
			continue
		else:
			rwDepth=float(rwDepth)
		rwMag=rw[49:54].strip()
		if rwMag=='' or isnumeric(str(rwMag))==False or rwMag==None:
			#rwMag=0.0
			continue
		else:
			rwMag=float(rwMag)
		rwMagType=rw[55:59].strip()
		rwNst=rw[60:64].strip()
		if rwNst=='':
			rwNst=0.0
		else:
			rwNst=float(rwNst)
		rwGap=rw[65:68].strip()
		rwClo=rw[69:73].strip()
		rwrms=rw[74:78].strip()
		if rwrms=='':
			rwrms=0.0
		else:
			rwrms=float(rwrms)		
		rwsrc=rw[79:83].strip()
		rwCatEventId=rw[84:96].strip()
		
		#anssList+=[[rw[0:22].strip(), float(rw[23:31].strip()), float(rw[32:41].strip()), float(rw[42:48].strip()), float(rw[49:54].strip()), rw[55:59].strip(), float(rw[60:64].strip()), rw[65:68].strip(), rw[69:73].strip(), float(rw[74:78].strip()), rw[79:83].strip(), rw[84:96].strip()]]
		anssList+=[[rwEvdt, rwLat, rwLon, rwDepth, rwMag, rwMagType, rwNst, rwGap, rwClo, rwrms, rwsrc, rwCatEventId]]
	return anssList

def isnumeric(value):
  return str(value).replace(".", "").replace("-", "").isdigit()
