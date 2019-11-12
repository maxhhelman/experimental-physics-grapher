#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 14:18:56 2019

@author: maxhelman
"""

import matplotlib.pyplot as plt

heavyz = [.1, .4, .2, .9, .6, .5, .2, 0, 0, 0, -.2, -.4, -.3, -.2, -.5, -.4, \
-.7, -.2, -.4, -.7]

heavyx = [75.2, 74.8, 74.3, 73.7, 72.7, 72.2, 72.9, 73.7, 74.2, 74.1, 
  74.4, 74.4, 74.2, 74.2, 72.1, 72.6, 72.7, 73.7, 73.4, 73.6]

lightz = [2.1, 2.1, .5, .7, .9, .5, .4, .1, .3, .3, .4, .2, 
  0, .2, .2, -.3, -.7, 0, -.6, -.9]

lightx = [57, 59.7, 62.1, 58.7, 60.4, 61.4, 60.2, 61.6, 61, 60.9, 
  60.7, 60.4, 59.8, 59.4, 59.2, 59.0, 60.6, 60.8, 61.0, 61.7]

plt.xlabel("z (cm)")
plt.ylabel("x (cm)")
plt.title("x and z Values in cm for Heavy Ball")
plt.plot(0,78.7,'ro')
plt.errorbar(0,78.7,.31, fmt = 'ro')
plt.scatter(heavyz, heavyx)


plt.figure()


plt.xlabel("z (cm)")
plt.ylabel("x (cm)")
plt.title("x and z Values in cm for Light Ball")
plt.plot(0,71.4,'ro')
plt.errorbar(0,71.4,.3 , fmt = 'ro')
plt.scatter(lightz, lightx)

plt.figure()

plt.hist(heavyz, normed=False, bins=5)
plt.ylabel('Occurances')
plt.xlabel('z (cm)')
plt.title('Histogram of z values for heavy ball')

plt.figure()

plt.hist(heavyx, normed=False, bins=5)
plt.ylabel('Occurances')
plt.xlabel('x (cm)')
plt.title('Histogram of x values for heavy ball')

plt.figure()

plt.hist(lightz, normed=False, bins=5)
plt.ylabel('Occurances')
plt.xlabel('z (cm)')
plt.title('Histogram of z values for light ball')

plt.figure()

plt.hist(lightx, normed=False, bins=5)
plt.ylabel('Occurances')
plt.xlabel('x (cm)')
plt.title('Histogram of x values for light ball')