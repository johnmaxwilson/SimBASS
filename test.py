# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 14:55:53 2014

@author: jmwilson
"""
#==============================================================================
# import numpy as np
# import math
# import matplotlib.pyplot as plt
# import pickle
# import time
# 
# from mpl_toolkits.basemap import Basemap
# from scipy.stats import norm
# from scipy.spatial import KDTree
# 
# import ETAS_tools as tools
# x = np.arange(-125.4, -113.1, 0.1)
# y = np.arange(31.5, 43.0, 0.1)
# xx, yy = np.meshgrid(x, y)
# xy = np.c_[xx.ravel(), yy.ravel()]
# 
# binmask = []
# for coord in xy:
#     if tools.isinPoly(coord):
#         binmask.append([0,0])
#     else:
#         binmask.append([1,1])
# binmask = np.array(binmask)
# 
# 
# tree = KDTree(maskedxy)
# res = tree.query_ball_point([-125.5, 31.3], 0.4)
# for ind in res:
#     print xy[ind]
# 
# #==============================================================================
# # plt.pcolor(x, y, tree.query(xy)[1].reshape(11, 11))
# # plt.plot(points[:,0], points[:,1], 'ko')
# # plt.show()
# #==============================================================================
#==============================================================================

#==============================================================================
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# 
# fig = plt.figure()
# 
# def f(x, y):
#     return np.sin(x) + np.cos(y)
# 
# test = np.zeros((25,25))
# rad = 0
# ims = []
# for step in range(70):
#     for i in range(len(test)):
#         for j in range(len(test[0])):
#             cellrad = math.sqrt((i-12)**2+(j-12)**2)
#             if cellrad < rad and cellrad > rad-1:
#                 pass#test[i,j] += 1
#             #
#         #
#     #
#     im = plt.imshow(test)
#     ims.append([im])
#     rad += 1
# 
# ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
#     repeat_delay=0)
# 
# #ani.save('/home/jmwilson/Desktop/dynamic_images.mp4')
#==============================================================================


#==============================================================================
# plt.show()
#==============================================================================


x = np.linspace(0,1,100)
plt.plot(x, -x*np.log(x))
plt.show()