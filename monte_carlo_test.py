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

def test_pgram(periods, powers, threshold, n_aliases, alias_with):
    """ID the most likely period and aliases.

    Adapted from k2spin by Stephanie Douglas and Kevin Covey

    Inputs
    ------
    periods, powers: array-like
    threshold: float
    n_aliases: int
    alias_with: float
        aliases to search for. Defaults to 0.25 days (6 hrs), the time
        between K2 thruster fires
    Outputs
    -------
    best_period: float
    best_power: float
    aliases: array
    is_clean: bool
    """
    logging.debug("Eval %d aliases with threshold %f",n_aliases, threshold)

    # Find the most likely period
    fund_loc = np.argmax(abs(powers))
    fund_period = periods[fund_loc]
    fund_power = powers[fund_loc]

    logging.debug("Fundamental %d Prot=%f Power=%f", fund_loc, fund_period,
                  fund_power)

    # and aliases
    for_aliases = np.arange(1.0, n_aliases+1.0) / alias_with
    inverse_fundamental = 1. / fund_period
    pos_aliases = 1. / (inverse_fundamental + for_aliases)
    neg_aliases = 1. / abs(inverse_fundamental - for_aliases)

    aliases = np.append(pos_aliases, neg_aliases)
    tot_aliases = len(aliases)

#    logging.debug("Aliases: {}".format(aliases))

    # Clip the best peak out of the periodogram
    to_clip = np.where(abs(periods - fund_period)<=fund_5percent)[0]
    #logging.debug(periods[to_clip])

    # Now clip out aliases
    for i, alias in enumerate(aliases):
        to_clip = np.union1d(to_clip,
                             np.where(abs(periods - alias)<=(0.02*alias))[0])

    clipped_periods = np.delete(periods, to_clip)
    clipped_powers = np.delete(powers, to_clip)


    # Find the maximum of the clipped periodogram
    clipped_max = np.argmax(clipped_powers)
    max_clip_power = clipped_powers[clipped_max]
    max_clip_period = clipped_periods[clipped_max]

    # Set clean = True if max of clipped periodogram < threshold*best_power
    if max_clip_power < (threshold * fund_power):
        is_clean = True
    else:
        is_clean = False
        logging.debug("Max clipped power = %f", max_clip_power)

    # Return best_period, best_power, aliases, is_clean
    return is_clean

def monte_carlo_injections(files):

    N = 100000
    short_P = 0.1
    long_P = 62.0
    pd = np.linspace(short_P, long_P, N)
    f = (2.0*pi)/pd
    # 2pi/f = 2pi/(2pi/pd) = pd
    p0 = np.log10(pd)


    for fcount,filename in enumerate(files):
        logging.info(filename)
        if filename.endswith("sysrem.txt"):

            injected_period = []
            newperiod = []
            max_power = []
            unique = []

            periods = 31.5 * random.random_sample(50) + 0.1
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

                    is_clean = test_pgram(pd, pgram_g, threshold=0.6, n_aliases=3, alias_with=1.0)

                    if is_clean is True:
                        unique = unique + [1]
                    else:
                        unique = unique + [0]

            data = {"injected_period": injected_period,
                    "new_period": newperiod,
                    "max_power": max_power,
                    "unique": unique}

            names = ["injected_period", "new_period", "max_power", "unique"]

            at.write(data,str(filename)[:-4]+".injections.csv",names=names,delimiter=",")
            if fcount>=10:
                break


if __name__=="__main__":

    if len(sys.argv)==1:
        print("Please provide a list of light curve files")
    else:
        listfile = at.read(sys.argv[1])
        file_list = listfile["filename"]
        monte_carlo_injections(file_list)
