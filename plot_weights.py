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
__credits__= 'Francesco de Gasperin'
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
    t1 = taql('SELECT TIME, ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False')
    # Average over all baselines (axis 0) and group the resulting data by antenna.
    t = taql('SELECT TIME, ANTENNA, MEANS(GAGGR(DATA), 0) AS DATA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM $t1 GROUPBY ANTENNA')
    #t1 = taql('SELECT ANTENNA1, DATA, WEIGHT_SPECTRUM, TIME, FIELD_ID FROM $msfile WHERE ANY(FLAG)==False')
    #t = taql('select ANTENNA, MEANS(GAGGR(DATA), 0) as DATA, MEANS(GAGGR(WEIGHT_SPECTRUM), 0) as WEIGHT from [SELECT ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM from $t1] group by ANTENNA')
    d = t.getcol('DATA')
    print d.shape
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
    ax.set_xlabel('Channel')
    ax.set_ylabel('Data')
    ax.set_yscale('log')
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    return 

def plot_weight_channel(msfile, pol=0):
    # Polarization indices are 0, 1, 2, 3 = XX, YY, XY, YX, respectively.
    print 'Plotting weights vs. channels for %s' % (msfile,)
    imgname ='weight_chan_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+'.png'
    # GMEANS (or MEANS(GAGGR(array), N) calculates the mean over the Nth axis.
    # Select only rows containing no flagged data.
    #t1 = taql('SELECT TIME, ANTENNA1 AS ANTENNA, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False')
    # Average over all baselines (axis 0) and group the resulting data by antenna.
    #t = taql('SELECT TIME, ANTENNA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM $t1 GROUPBY ANTENNA')
    t1 = taql('SELECT TIME, ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False')
    #t1 = taql('SELECT TIME, ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE FLAG_ROW==False')
    # Average over all baselines (axis 0) and group the resulting data by antenna.
    t = taql('SELECT TIME, ANTENNA, DATA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM $t1 GROUPBY ANTENNA')
    w = t.getcol('WEIGHT')
    print w.shape
    # Select all antennas, all channels and polarization pol.
    weights = t.getcol('WEIGHT')[:, :, pol]
    antennas = t.getcol('ANTENNA')
    print len(antennas), ' antennas'
    antenna_names = taql('SELECT NAME FROM '+msfile+'/ANTENNA')
    # Obtain channel frequencies in Hz.
    chan_freq = taql('SELECT CHAN_FREQ FROM '+msfile+'/SPECTRAL_WINDOW')
    # Select the first table, column CHAN_FREQ and convert to MHz.
    freq = chan_freq[0]['CHAN_FREQ'] * 1e-6
    # Plot the results.
    print 'Plotting...'
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    for w, c, a in zip(weights, colors, antennas):
        ax.scatter(freq, w, color=c, marker='.', label=antenna_names[a]['NAME'])
        #ax.plot(freq, w, '-x', color=c, label=antenna_names[a]['NAME'])
    ax.set_xlim(min(freq), max(freq))
    ax.set_ylim(np.min(weights), np.max(weights))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
    major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
    ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Weights')
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    return 

def plot_weight_time(msfile, dt_elev=100, plot_time_unit='h'):
    print 'Plotting weights vs. time for %s' % (msfile,)
    imgname ='weight_time_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+'.png'
    # Select the time, weights and elevation of ANTENNA1 averaging over baselines/antennas.
    # Select only unflagged data.
    t1 = taql('select ANTENNA1, ANTENNA2, DATA, WEIGHT_SPECTRUM, TIME, FIELD_ID from $msfile where any(FLAG)==False')
    w = t1.getcol('WEIGHT_SPECTRUM')
    print w.shape
    print len(np.unique(t1.getcol('ANTENNA1')))
    print len(np.unique(t1.getcol('ANTENNA2')))
    # Select time, weights and elevation.
    # This gets the average elevation w.r.t. to antenna 1, where the average is taken over all baselines.
    #t = taql('SELECT TIME, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT, MSCAL.AZEL()[1] AS ELEV FROM $t1 GROUPBY TIME')
    t = taql('SELECT TIME, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT, MEANS(GAGGR(MSCAL.AZEL1()[1]), 0) AS ELEV FROM $t1 GROUPBY TIME')
    t2 = taql('SELECT TIME, REAL(DATA) AS DATA_REAL FROM $t1 GROUPBY TIME')
    weights = t.getcol('WEIGHT')
    time = t.getcol('TIME')
    elevation = t.getcol('ELEV')

    datar = t2.getcol('DATA_REAL')
    # Calculate variance over 10 timestamp intervals. The last interval may be shorter.
    # datar shape is (timestamps, channels, polarizations)
    delta = 10
    variance = np.ones(shape=(len(time)//delta, weights.shape[1], weights.shape[2]))
    for i in xrange(len(time)//delta):
        datar_shifted = np.roll(datar, -1)
        datar = datar - datar_shifted
        v = np.var(datar[delta*i: delta*i+delta,:,:], axis=0)
        variance[i] = 1. / v
    print variance.shape
    print len(variance[:, 0, 0])
    print variance[0, :, 0]
    print len(time), len(time[::delta])
    
    # Select the weights for all timestamps one channel and one polarization (in that order).
    weights = t.getcol('WEIGHT')[:, 0, 0]
    # Plot weights for elevation and save the image.
    print 'Plotting...'
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)

    if plot_time_unit.lower() == 'h':
        time_h = time / 3600
        time_h = (time_h - math.floor(time_h[0]/24) * 24)
        ax.scatter(time_h, weights, marker='.', color='k', label='Weights')
        #print len(time_h[::delta])
        #print time_h[::delta].shape, variance[:, 0, 0].shape
        ax.scatter(time_h[::delta][:-1], variance[:, 0, 0], marker='x', label='Boxed variance $\\Delta=%d$'%(delta,))
        ax.set_yscale('log')
        # Plot the elevation as a function of time.
        ax_elev = ax.twinx()
        ax_elev.plot(time_h, elevation * 180/np.pi, 'b', label='Elevation')

        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%.2d:%.2d:%.2d' % (int(x), (x%1)*60, (((x%1)*60)%1 * 60))))
        ax.set_xlim(min(time_h), max(time_h))
        ax.set_xlabel('Time [h]')
    elif plot_time_unit.lower() == 's':
        ax.scatter(time, weights, marker='.', color='k', label='Weights')
        # Plot the elevation as a function of time.
        ax_elev = ax.twinx()
        ax_elev.plot(time, elevation * 180/np.pi, 'b', label='Elevation')
        ax.set_xlim(min(time), max(time))
        ax.set_xlabel('Time [s]')
    #ax.set_ylim(np.min(weights), np.max(weights))
    #ax.set_ylim(np.min(variance), np.max(variance))
    ax.set_ylim(min(np.min(weights), np.min(variance)), 10 * max(np.max(weights), np.max(variance)))
    ax.set_ylabel('Weights')
    ax_elev.set_ylabel('Elevation [deg]')

    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax_elev.get_legend_handles_labels()
    leg = ax.legend(handles+handles2, labels+labels2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    return

if __name__ == '__main__':
    # Get the MS filename.
    msfile = sys.argv[1]
    #plot_weight_channel(msfile)
    plot_weight_time(msfile, plot_time_unit='h')
    #plot_data_channel(msfile, pol=3)
