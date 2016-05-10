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

            injected_period = []
            newperiod = []
            max_power = []
            unique = []

            periods = 31.5 * random.random_sample(500,) + 0.1
            new = np.array(periods).tolist()

            injected_period = injected_period + new * 5

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
                    newperiod = newperiod + [pd[max_loc]]
                    max_power = max_power + [pgram_g[max_loc]]

                    newpgram_g = np.zeros(N)
                    keep = np.where(pgram_g<max(pgram_g))[0]
                    newpgram_g[keep] = pgram_g[keep]
                    # for m in range(0, len(pgram_g)):
                    #     if pgram_g[m] < max(pgram_g):
                    #         newpgram_g[m] = pgram_g[m]

                    if max(newpgram_g) <= 0.6 * max(pgram_g):
                        unique = unique + [1]
                    elif max(newpgram_g) > 0.6 * max(pgram_g):
                        unique = unique + [0]

            data = {"injected_period": injected_period,
                    "new_period": newperiod,
                    "max_power": max_power,
                    "unique": unique}

            names = ["injected_period", "new_period", "max_power", "unique"]

            at.write(data,"/Users/Leo/Desktop/Result/"+str(filename)[0:19]+".csv",names=names,delimiter=" ")


if __name__=="__main__":

    if len(sys.argv)==1:
        print("Please provide a list of light curve files")
    else:
        listfile = at.read(sys.argv[1])
        file_list = listfile["filename"]
        monte_carlo_injections(file_list)
