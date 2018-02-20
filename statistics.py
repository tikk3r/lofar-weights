#!/usr/bin/env python
""" Provides acces to the weights of an MS file. """

import datetime
import math
import os
import sys

import astropy.time
import numpy as np

from casacore.tables import taql
from matplotlib import cm
from matplotlib import dates
from matplotlib.pyplot import figure, show
from matplotlib import ticker

__author__ = 'Frits Sweijen'
__version__ = '1.0.0'
__maintainer__ = 'Frits Sweijen'
__email__ = 'sweijen <at> strw.leidenuniv.nl'
__status__ = 'Development'

def plot_data_channel(msfile, pol=0):
    # Polarization indices are 0, 1, 2, 3 = XX, YY, XY, YX, respectively.
    print 'Plotting data vs. channel for %s' % (msfile,)
    imgname ='data_chan_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+'.png'
    # GMEANS (or MEANS(GAGGR(array), N) calculates the mean over the Nth axis.
    # Select only rows containing no flagged data.
    t1 = taql('SELECT TIME, ANTENNA1, ANTENNA2, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False')
    #t1 = taql('SELECT TIME, ANTENNA1, ANTENNA2, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE FLAG_ROW==False')
    #t1 = taql('SELECT TIME, ANTENNA1, ANTENNA2, DATA, WEIGHT_SPECTRUM FROM $msfile')
    # Average over all baselines (axis 0) and group the resulting data by antenna.
    t = taql('SELECT TIME, GSTDDEV(REAL(DATA)) AS DATA_REAL_STD, GSTDDEV(IMAG(DATA)) AS DATA_IMAG_STD FROM $t1')
    data = t.getcol('DATA_REAL_STD')
    print data.shape
    data_imag = data.imag
    data_real = data.real
    print len(data_imag[0,0,0])
    sys.exit(0)
    # Select all antennas, all channels and polarization pol.
    data = np.abs(t.getcol('DATA')[:, :, pol])
    antennas = t.getcol('ANTENNA')
    print len(antennas), ' antennas'
    antenna_names = taql('SELECT NAME FROM '+msfile+'/ANTENNA')
    # Obtain channel frequencies in Hz.
    #chan_freq = taql('SELECT CHAN_FREQ FROM '+msfile+'/SPECTRAL_WINDOW')
    # Select the first table, column CHAN_FREQ and convert to MHz.
    #freq = chan_freq[0]['CHAN_FREQ'] * 1e-6
    # Plot the results.
    print 'Plotting...'
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(data))))
    for d, c, a in zip(data, colors, antennas):
        #ax.scatter(freq, d, color=c, marker='.', label=antenna_names[a]['NAME'])
        ax.scatter(xrange(len(d)), d, color=c, marker='.', label=antenna_names[a]['NAME'])
    #ax.set_xlim(min(freq), max(freq))
    ax.set_xlim(0, len(d))
    ax.set_ylim(np.min(data), np.max(data))
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
    #major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
    #ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Data')
    ax.set_yscale('log')
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    #fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    return 

if __name__ == '__main__':
    # Get the MS filename.
    msfile = sys.argv[1]
    plot_data_channel(msfile)
