import sys
import logging
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
      
logging.basicConfig(level=logging.DEBUG)  

def monte_carlo_figure2(files):
    
    for filename in files:
        logging.info(filename)
        if filename.endswith(".csv"):
            
            tbl = at.read(str(root)+'/'+str(file))
            unique = tbl["unique"]
            injected_period = tbl["injected_period"]
            new_period = tbl["new_period"]
            max_power = tbl["max_power"]
            
            a = np.where((max_power > 0.2) & (max_power < 0.3))[0]
            plt.scatter(injected_period[a], new_period[a], marker='D', color='lightcoral', s=0.5)
            
            b = np.where((unique == 1) & (max_power > 0.3) & (max_power < 0.6))[0]
            plt.scatter(injected_period[b], new_period[b], marker='D', color='khaki', s=0.5)
            
            c = np.where((unique == 1) & (max_power > 0.6))[0]
            plt.scatter(injected_period[c], new_period[c], marker='D', color='black', s=0.5)        
                                  
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
    plt.legend(handles=[black_patch, yellow_patch, red_patch], loc=2)    
    plt.savefig("/Users/Leo/Desktop/Figure 2.eps")
    
    
    
if __name__=="__main__":
    
    if len(sys.argv)==1:
        print("Please provide a list of injection files")
    else:
        listfile = at.read(sys.argv[1])
        file_list = listfile["filename"]
        monte_carlo_figure2(file_list)
