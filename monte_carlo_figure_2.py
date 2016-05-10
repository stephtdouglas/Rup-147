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
import astropy.io.ascii as at
      

for root, dirs, files in walk("/Users/Leo/Desktop/Result"):
    for file in files:
        if file.endswith(".csv"):
            
            tbl = at.read("/Users/Leo/Desktop/Result.csv")
            unique = tbl["unique"]
            injected_period = tbl["injected_period"]
            new_period = tbl["new_period"]
            max_power = tbl["max_power"]
            
            for i in range(0, len(unique)):
                if max_power[i] > 0.2 and max_power[i] < 0.3:
                    plt.plot(injected_period[i], new_period[i], marker='D', c='lightcoral', markersize=0.5, fillstyle='none')
                if unique[i] == 1:
                    if max_power[i] > 0.3 and max_power[i] < 0.6:
                        plt.plot(injected_period[i], new_period[i], marker='D', c='khaki', markersize=0.5, fillstyle='none')
                    if max_power[i] > 0.6:
                        plt.plot(injected_period[i], new_period[i], marker='D', c='black', markersize=0.5, fillstyle='none')        
                                  
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
plt.savefig("/Users/Leo/Desktop/Figure 2.eps")