# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:07:35 2015

@author: jmwilson
"""

#==============================================================================
# calipoly = [[-129.0, 42.75], [-121.5, 29.5], [-113.5, 29.5], [-114.0, 35.6], [-119.3, 39.4], [-119.5, 42.75], [-129.0, 42.75]]
#==============================================================================

#Schorlemmer 07 RELM testing area
arealats = [43.0, 43.0, 39.4, 35.7, 34.3, 32.9, 32.2, 31.7, 31.5, 31.9, 32.8, 33.7, 34.2, 37.7, 40.2, 40.5]
arealons = [-125.2, -119.0, -119.0, -114.0, -113.1, -113.5, -113.6, -114.5, -117.1, -117.9, -118.4, -121.0, -121.6, -123.8, -125.4, -125.4]

calipoly = zip(arealons, arealats)
#print calipoly

def isinPoly(pt, polyvecs=None):
    # be careful; this needs to be a polyvec, not just a polygon...
    if polyvecs==None: polyvecs=getPolyVecs(calipoly)
    if type(polyvecs[0][0])==type(0) or type(polyvecs[0][0]) == type(1.0):
        polyvecs=getPolyVecs(polyvecs)
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
    
        #we're in the y-domain of this vector.
        if x0<x1 and x0<x2:
            nlefts+=1
        else:
            #might still be left...
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
                     
def getPolyVecs(poly):
    polyvecs=[]
    for i in xrange(len(poly)):
        polyvecs+=[[poly[i], poly[(i+1)%(len(poly))]]]
        if polyvecs[-1][0]==polyvecs[-1][1]: polyvecs.pop()
    return polyvecs
    
    
points = [[-117,34]]

polyvectors = getPolyVecs(calipoly)

for point in points:
    print isinPoly(point, polyvecs = polyvectors)
    
    
plt.figure()  
m = Basemap(projection='cyl', llcrnrlat=31, urcrnrlat=44, llcrnrlon=-126, urcrnrlon=-113, resolution='i')
m.drawcoastlines()
m.drawmapboundary(fill_color='PaleTurquoise')
m.fillcontinents(color='lemonchiffon',lake_color='PaleTurquoise', zorder=0)
m.drawstates()
m.drawcountries()
m.drawparallels(np.arange(31, 45, 2),labels=[1,0,0,0])
m.drawmeridians(np.arange(-124, -113, 2),labels=[0,0,0,1])

for vect in polyvectors:
    m.plot([vect[0][0], vect[1][0]], [vect[0][1], vect[1][1]], 'm-')
m.plot(points[0][0], points[0][1], 'k.')

plt.show()