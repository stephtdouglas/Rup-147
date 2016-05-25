from __future__ import division
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

def monte_carlo_figure1(files):

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

    for filename in files:
        logging.info(filename)
        if filename.endswith(".csv"):

            tbl = at.read(filename)
            unique = tbl["unique"]
            injected_period = tbl["injected_period"]
            new_period = tbl["new_period"]
            max_power = tbl["max_power"]

            for i in range(0, len(unique)):
                if max_power[i] > 0 and max_power[i] <= 0.1:
                    zone1 = zone1 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone1count = zone1count + 1

                elif max_power[i] > 0.1 and max_power[i] <= 0.2:
                    zone2 = zone2 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone2count = zone2count + 1

                elif max_power[i] > 0.2 and max_power[i] <= 0.3:
                    zone3 = zone3 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone3count = zone3count + 1

                elif max_power[i] > 0.3 and max_power[i] <= 0.4:
                    zone4 = zone4 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone4count = zone4count + 1

                elif max_power[i] > 0.4 and max_power[i] <= 0.5:
                    zone5 = zone5 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone5count = zone5count + 1

                elif max_power[i] > 0.5 and max_power[i] <= 0.6:
                    zone6 = zone6 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone6count = zone6count + 1

                elif max_power[i] > 0.6 and max_power[i] <= 0.7:
                    zone7 = zone7 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone7count = zone7count + 1

                elif max_power[i] > 0.7 and max_power[i] <= 0.8:
                    zone8 = zone8 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone8count = zone8count + 1

                elif max_power[i] > 0.8 and max_power[i] <= 0.9:
                    zone9 = zone9 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone9count = zone9count + 1

                elif max_power[i] > 0.9 and max_power[i] <= 1.0:
                    zone10 = zone10 + 1
                    if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                        zone10count = zone10count + 1

                if unique[i] == 1:
                    if max_power[i] > 0 and max_power[i] <= 0.1:
                        unique1 = unique1 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique1count = unique1count + 1

                    elif max_power[i] > 0.1 and max_power[i] <= 0.2:
                        unique2 = unique2 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique2count = unique2count + 1

                    elif max_power[i] > 0.2 and max_power[i] <= 0.3:
                        unique3 = unique3 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique3count = unique3count + 1

                    elif max_power[i] > 0.3 and max_power[i] <= 0.4:
                        unique4 = unique4 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique4count = unique4count + 1

                    elif max_power[i] > 0.4 and max_power[i] <= 0.5:
                        unique5 = unique5 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique5count = unique5count + 1

                    elif max_power[i] > 0.5 and max_power[i] <= 0.6:
                        unique6 = unique6 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique6count = unique6count + 1

                    elif max_power[i] > 0.6 and max_power[i] <= 0.7:
                        unique7 = unique7 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique7count = unique7count + 1

                    elif max_power[i] > 0.7 and max_power[i] <= 0.8:
                        unique8 = unique8 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique8count = unique8count + 1

                    elif max_power[i] > 0.8 and max_power[i] <= 0.9:
                        unique9 = unique9 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique9count = unique9count + 1

                    elif max_power[i] > 0.9 and max_power[i] <= 1.0:
                        unique10 = unique10 + 1
                        if ((abs(new_period[i] - injected_period[i]))/injected_period[i]) <= 0.03:
                            unique10count = unique10count + 1


    if zone1 == 0:
        percentage1 = 0
    if zone1 != 0:
        percentage1 = zone1count/zone1
    if zone2 == 0:
        percentage2 = 0
    if zone2 != 0:
        percentage2 = zone2count/zone2
    if zone3 == 0:
        percentage3 = 0
    if zone3 != 0:
        percentage3 = zone3count/zone3
    if zone4 == 0:
        percentage4 = 0
    if zone4 != 0:
        percentage4 = zone4count/zone4
    if zone5 == 0:
        percentage5 = 0
    if zone5 != 0:
        percentage5 = zone5count/zone5
    if zone6 == 0:
        percentage6 = 0
    if zone6 != 0:
        percentage6 = zone6count/zone6
    if zone7 == 0:
        percentage7 = 0
    if zone7 != 0:
        percentage7 = zone7count/zone7
    if zone8 == 0:
        percentage8 = 0
    if zone8 != 0:
        percentage8 = zone8count/zone8
    if zone9 == 0:
        percentage9 = 0
    if zone9 != 0:
        percentage9 = zone9count/zone9
    if zone10 == 0:
        percentage10 = 0
    if zone10 != 0:
        percentage10 = zone10count/zone10

    percentage = [percentage1, percentage2, percentage3, percentage4, percentage5, percentage6, percentage7, percentage8, percentage9, percentage10]

    if unique1 == 0:
        uniquepercentage1 = 0
    if unique1 != 0:
        uniquepercentage1 = unique1count/unique1
    if unique2 == 0:
        uniquepercentage2 = 0
    if unique2 != 0:
        uniquepercentage2 = unique2count/unique2
    if unique3 == 0:
        uniquepercentage3 = 0
    if unique3 != 0:
        uniquepercentage3 = unique3count/unique3
    if unique4 == 0:
        uniquepercentage4 = 0
    if unique4 != 0:
        uniquepercentage4 = unique4count/unique4
    if unique5 == 0:
        uniquepercentage5 = 0
    if unique5 != 0:
        uniquepercentage5 = unique5count/unique5
    if unique6 == 0:
        uniquepercentage6 = 0
    if unique6 != 0:
        uniquepercentage6 = unique6count/unique6
    if unique7 == 0:
        uniquepercentage7 = 0
    if unique7 != 0:
        uniquepercentage7 = unique7count/unique7
    if unique8 == 0:
        uniquepercentage8 = 0
    if unique8 != 0:
        uniquepercentage8 = unique8count/unique8
    if unique9 == 0:
        uniquepercentage9 = 0
    if unique9 != 0:
        uniquepercentage9 = unique9count/unique9
    if unique10 == 0:
        uniquepercentage10 = 0
    if unique10 != 0:
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
    plt.savefig("Figure1.eps")  



if __name__=="__main__":

    if len(sys.argv)==1:
        print("Please provide a list of injection files")
    else:
        listfile = at.read(sys.argv[1])
        file_list = listfile["filename"]
        monte_carlo_figure1(file_list)
