from os import walk
from numpy import loadtxt, zeros, random, mean, std, linspace, argmax, percentile

import matplotlib.pyplot as plt

from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap

from math import log10, sin, pi

from PyAstronomy.pyasl import foldAt
from PyAstronomy.funcFit import SinusFit1d

import numpy as np
from scipy.optimize import leastsq

import matplotlib.patches as mpatches





N = 100000
short_P = 0.1
long_P = 62.0
pd = linspace(short_P, long_P, float(N))
f = (2.0*pi)/pd
p0 = zeros(N)
for i in range (N):
    p0[i] = log10((2.0*pi)/f[i])
    
for root, dirs, files in walk("/Users/Leo/Desktop/test"):
    for file in files:
        if file.endswith("sysrem.txt"):
            periods = 31.5 * random.random_sample((500, 1)) + 0.1
            file00 = open(str(root)+'/'+str(file), 'r')
            time, mag, error, ra, dec, flag = loadtxt(file00, dtype='float',delimiter=' ', usecols=(0,1,2,3,4,7), unpack=True)
            inputLC00 = zip(time[0:], mag[0:], error[0:], ra[0:], dec[0:], flag[0:])
            inputLC01 = [[time,mag,error,ra,dec,flag] for (time,mag,error,ra,dec,flag) in inputLC00 if \
			((error != 0.) and (ra != 0.) and (dec != 0.) and (flag == 0.))]
	    std01 = std([item[1] for item in inputLC01])
	    amplitude = [std01*0.3, std01*0.6, std01*0.9, std01*1.2, std01*1.5]
	    time = [item[0] for item in inputLC01]
	    mag = [item[1] for item in inputLC01]
	    error = [item[2] for item in inputLC01]
	    
	    for i in range(0, len(amplitude)):
	        for j in range(0, len(periods)):
	            sinusoid = [amplitude[i] * sin(2*pi*(x-time[0])/periods[j]) for x in time]
	            newmag = [a+b for a,b in zip(mag, sinusoid)]
                        
                    pgram_g = lomb_scargle(time, newmag, error, f, generalized = True)
                    for k in range(0, len(pgram_g)):
                        if pgram_g[k] == max(pgram_g):
                            newperiod = pd[k]
                    
                    if max(pgram_g) > 0.2 and max(pgram_g) < 0.3:
                        plt.plot(periods[j], newperiod, marker='D', c='lightcoral', markersize=0.5, fillstyle='none')
                    
                    newpgram_g = [0] * len(pgram_g)
                    for m in range(0, len(pgram_g)):
                        if pgram_g[m] < max(newpgram_g):
                            newpgram_g[m] = pgram_g[m]
                            
                    if max(newpgram_g) <= 0.6 * max(pgram_g):
        
                        if max(pgram_g) > 0.3 and max(pgram_g) < 0.6:
                            plt.plot(periods[j], newperiod, marker='D', c='khaki', markersize=0.5, fillstyle='none')
                        
                        if max(pgram_g) > 0.6:
                            plt.plot(periods[j], newperiod, marker='D', c='black', markersize=0.5, fillstyle='none')
      
      
      

ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([0.1,32.0])
ax.set_ylim([0.1,32.0])
plt.xlabel('Injected P_rot (days)')
plt.ylabel('Recovered P_rot (days)')
red_patch = mpatches.Patch(color='lightcoral', label='0.2 < power < 0.3')
yellow_patch = mpatches.Patch(color='khaki', label='unique = 1 & 0.3 < power < 0.6')
black_patch = mpatches.Patch(color='black', label='unique = 1 & power > 0.6')
plt.legend(handles=[black_patch, yellow_patch, red_patch])    
plt.savefig("/Users/Leo/Desktop/Result 2.eps")