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
import astropy.io.ascii as at
import matplotlib.patches as mpatches

logging.basicConfig(level=logging.DEBUG)

def monte_carlo_injections(files):

    count = 0
    zone1 = 0
    zone2 = 0
    zone3 = 0
    zone4 = 0
    zone5 = 0
    zone6 = 0
    zone7 = 0
    zone8 = 0
    zone9 = 0
    zone10 = 0
    zone1count = 0
    zone2count = 0
    zone3count = 0
    zone4count = 0
    zone5count = 0
    zone6count = 0
    zone7count = 0
    zone8count = 0
    zone9count = 0
    zone10count = 0
    
    unique = 0
    uniquecount = 0
    unique1 = 0
    unique2 = 0
    unique3 = 0
    unique4 = 0
    unique5 = 0
    unique6 = 0
    unique7 = 0
    unique8 = 0
    unique9 = 0
    unique10 = 0
    unique1count = 0
    unique2count = 0
    unique3count = 0
    unique4count = 0
    unique5count = 0
    unique6count = 0
    unique7count = 0
    unique8count = 0
    unique9count = 0
    unique10count = 0

    N = 100000
    short_P = 0.1
    long_P = 62.0
    pd = np.linspace(short_P, long_P, N)
    f = (2.0*pi)/pd
    # 2pi/f = 2pi/(2pi/pd) = pd
    p0 = np.log10(pd)

    for filename in files:
        logging.info(filename)
        if filename.endswith("sysrem.txt"):
            periods = 31.5 * random.random_sample((500, 1)) + 0.1
            file00 = open(filename, 'r')
            time, mag, error, ra, dec, flag = loadtxt(file00, dtype='float',delimiter=' ',
                                                      usecols=(0,1,2,3,4,7), unpack=True)
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
                    logging.debug("%d %0.3f %d %0.3f",i,amplitude[i],j,periods[j])
	            sinusoid = [amplitude[i] * sin(2*pi*(x-time[0])/periods[j]) for x in time]
	            newmag = [a+b for a,b in zip(mag, sinusoid)]
                        
                    pgram_g = lomb_scargle(time, newmag, error, f, generalized = True)

                    # Find the peak of the new periodogram
                    max_loc = np.argmax(pgram_g)
                    newperiod = pd[max_loc]
                    max_power = pgram_g[max_loc]
                    if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                        count = count + 1

#part 1
                    if max_power > 0 and max_power <= 0.1:
                        zone1 = zone1 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone1count = zone1count + 1
                    elif max_power > 0.1 and max_power <= 0.2:
                        zone2 = zone2 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone2count = zone2count + 1
                    elif max_power > 0.2 and max_power <= 0.3:
                        zone3 = zone3 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone3count = zone3count + 1
                    elif max_power > 0.3 and max_power <= 0.4:
                        zone4 = zone4 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone4count = zone4count + 1
                    elif max_power > 0.4 and max_power <= 0.5:
                        zone5 = zone5 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone5count = zone5count + 1
                    elif max_power > 0.5 and max_power <= 0.6:
                        zone6 = zone6 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone6count = zone6count + 1
                    elif max_power > 0.6 and max_power <= 0.7:
                        zone7 = zone7 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone7count = zone7count + 1
                    elif max_power > 0.7 and max_power <= 0.8:
                        zone8 = zone8 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone8count = zone8count + 1
                    elif max_power > 0.8 and max_power <= 0.9:
                        zone9 = zone9 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone9count = zone9count + 1
                    elif max_power > 0.9 and max_power <= 1.0:
                        zone10 = zone10 + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            zone10count = zone10count + 1
    
                    
                    # I think this is deleting the main peak and looking for other peaks???
                    newpgram_g = np.zeros(N)
                    max_new = max(newpgram_g)
                    less_than = (pgram_g < max_new)
                    logging.debug(len(np.where(less_than)[0]))
                    newpgram_g[less_than] = pgram_g[less_than]
                    # There's something wrong with the next three (now commented out) lines
                    # the max of the new periodogram starts out at 0, so nothing in the old
                    # periodogram will ever be lower than that. So it's wasting a bunch of 
                    # code time to evaluate to the same array of 0s anyway
#                    for m in range(0, len(pgram_g)):
#                        if pgram_g[m] < max(newpgram_g):
#                            newpgram_g[m] = pgram_g[m]
                    if max(newpgram_g) <= 0.6 * max_power:
                        unique = unique + 1
                        if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                            uniquecount = uniquecount + 1
                            
                        if max_power > 0 and max_power <= 0.1:
                            unique1 = unique1 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique1count = unique1count + 1
                        elif max_power > 0.1 and max_power <= 0.2:
                            unique2 = unique2 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique2count = unique2count + 1
                        elif max_power > 0.2 and max_power <= 0.3:
                            unique3 = unique3 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique3count = unique3count + 1
                        elif max_power > 0.3 and max_power <= 0.4:
                            unique4 = unique4 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique4count = unique4count + 1
                        elif max_power > 0.4 and max_power <= 0.5:
                            unique5 = unique5 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique5count = unique5count + 1
                        elif max_power > 0.5 and max_power <= 0.6:
                            unique6 = unique6 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique6count = unique6count + 1
                        elif max_power > 0.6 and max_power <= 0.7:
                            unique7 = unique7 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique7count = unique7count + 1
                        elif max_power > 0.7 and max_power <= 0.8:
                            unique8 = unique8 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique8count = unique8count + 1
                        elif max_power > 0.8 and max_power <= 0.9:
                            unique9 = unique9 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique9count = unique9count + 1
                        elif max_power > 0.9 and max_power <= 1.0:
                            unique10 = unique10 + 1
                            if ((abs(newperiod - periods[j]))/periods[j]) <= 0.03:
                                unique10count = unique10count + 1
        break
                        
    print count, unique, uniquecount

    if zone1 == 0:
        percentage1 = 0
    else:
        percentage1 = zone1count/zone1

    if zone2 == 0:
        percentage2 = 0
    else:
        percentage2 = zone2count/zone2

    if zone3 == 0:
        percentage3 = 0
    else:
        percentage3 = zone3count/zone3

    if zone4 == 0:
        percentage4 = 0
    else:
        percentage4 = zone4count/zone4

    if zone5 == 0:
        percentage5 = 0
    else:
        percentage5 = zone5count/zone5

    if zone6 == 0:
        percentage6 = 0
    else:
        percentage6 = zone6count/zone6

    if zone7 == 0:
        percentage7 = 0
    else:
        percentage7 = zone7count/zone7

    if zone8 == 0:
        percentage8 = 0
    else:
        percentage8 = zone8count/zone8

    if zone9 == 0:
        percentage9 = 0
    else:
        percentage9 = zone9count/zone9

    if zone10 == 0:
        percentage10 = 0
    else:
        percentage10 = zone10count/zone10
    
    percentage = [percentage1, percentage2, percentage3, percentage4, percentage5, percentage6, percentage7, percentage8, percentage9, percentage10]

    if unique1 == 0:
        uniquepercentage1 = 0
    else:
        uniquepercentage1 = unique1count/unique1

    if unique2 == 0:
        uniquepercentage2 = 0
    else:
        uniquepercentage2 = unique2count/unique2

    if unique3 == 0:
        uniquepercentage3 = 0
    else:
        uniquepercentage3 = unique3count/unique3

    if unique4 == 0:
        uniquepercentage4 = 0
    else:
        uniquepercentage4 = unique4count/unique4

    if unique5 == 0:
        uniquepercentage5 = 0
    else:
        uniquepercentage5 = unique5count/unique5

    if unique6 == 0:
        uniquepercentage6 = 0
    else:
        uniquepercentage6 = unique6count/unique6

    if unique7 == 0:
        uniquepercentage7 = 0
    else:
        uniquepercentage7 = unique7count/unique7

    if unique8 == 0:
        uniquepercentage8 = 0
    else:
        uniquepercentage8 = unique8count/unique8

    if unique9 == 0:
        uniquepercentage9 = 0
    else:
        uniquepercentage9 = unique9count/unique9

    if unique10 == 0:
        uniquepercentage10 = 0
    else:
        uniquepercentage10 = unique10count/unique10

    uniquepercentage = [uniquepercentage1, uniquepercentage2, uniquepercentage3, uniquepercentage4, uniquepercentage5, uniquepercentage6, uniquepercentage7, uniquepercentage8, uniquepercentage9, uniquepercentage10]

    x1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    plt.bar(x1, percentage, 0.05, color="lightcoral")
    x2 = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
    plt.bar(x2, uniquepercentage, 0.05, color="khaki")
    plt.xlabel('Peak Periodogram Power')
    plt.ylabel('Fraction of P_rot Recovered')
    red_patch = mpatches.Patch(color='lightcoral', label='no cuts')
    yellow_patch = mpatches.Patch(color='khaki', label='unique = 1')
    plt.legend(handles=[red_patch, yellow_patch])    
    plt.savefig("result1.eps",bbox_inches="tight")


if __name__=="__main__":
    
    if len(sys.argv)==1:
        print("Please provide a list of light curve files")
    else:
        listfile = at.read(sys.argv[1])
        file_list = listfile["filename"]
        monte_carlo_injections(file_list)
