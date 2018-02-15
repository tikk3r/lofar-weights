#!/usr/bin/env python
# Self-explanatory imports.
import os
import sys

import numpy as np

from casacore.tables import taql
from matplotlib import cm
from matplotlib.pyplot import figure, show
from matplotlib import ticker

@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%.2f" % x

@ticker.FuncFormatter
def minor_formatter(x, pos):
    return "%.2f" % x

def plot_weight_channel(msfile):
    print 'Plotting weights vs. channels for %s' % (msfile,)
    imgname ='weight_chan_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+'.png'
    # GMEANS (or MEANS(GAGGR(array), N) calculates the mean over the Nth axis.
    #t = taql('SELECT TIME, ANTENNA, MEANS(GAGGR(DATA),0) AS DATA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM [SELECT ANTENNA1 AS ANTENNA, TIME, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False], [SELECT ANTENNA2 AS ANTENNA, TIME, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False] GROUPBY ANTENNA')
    #t = taql('SELECT TIME, ANTENNA, MEANS(GAGGR(DATA),0) AS DATA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM [SELECT ANTENNA1 AS ANTENNA, TIME, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False] GROUPBY ANTENNA')
    #t = taql('SELECT TIME, ANTENNA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM [SELECT ANTENNA1 AS ANTENNA, TIME, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False] GROUPBY ANTENNA')
    #t = taql('SELECT TIME, ANTENNA1 AS ANTENNA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM $msfile WHERE ANY(FLAG)==False GROUPBY ANTENNA1')
    # Select unflagged data.
    t1 = taql('SELECT TIME, ANTENNA1 AS ANTENNA, WEIGHT_SPECTRUM FROM $msfile WHERE ANY(FLAG)==False')
    # Average over all baselines (axis 0) and group the resulting data by antenna.
    t = taql('SELECT TIME, ANTENNA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM $t1 GROUPBY ANTENNA')
    w = t.getcol('WEIGHT')
    print w.shape
    # Select all stations(?), all channels and polarization 3.
    weights = t.getcol('WEIGHT')[:, :, 3]
    antennas = t.getcol('ANTENNA')
    print len(antennas), ' antennas'
    antenna_names = taql('SELECT NAME FROM '+msfile+'/ANTENNA')
    # Obtain channel frequencies in Hz.
    chan_freq = taql('SELECT CHAN_FREQ FROM '+msfile+'/SPECTRAL_WINDOW')
    print 'Channel frequencies: '
    # Select the first table, column CHAN_FREQ and convert to MHz.
    freq = chan_freq[0]['CHAN_FREQ'] * 1e-6
    print freq
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
    ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Weights')
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    return 

def plot_weight_time(msfile):
    print 'Plotting weights vs. time for %s' % (msfile,)
    imgname ='weight_time_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+'.png'
    # Select the time, weights and elevation of ANTENNA1 averaging over baselines/antennas.
    #t = taql('SELECT TIME, ANTENNA1 AS ANTENNA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT, MSCAL.AZEL1()[1] AS ELEV FROM $msfile WHERE ANY(FLAG)==False GROUPBY TIME')
    # Select only unflagged data.
    t1 = taql('select ANTENNA1, DATA, WEIGHT_SPECTRUM, TIME, FIELD_ID from $msfile where any(FLAG)==False')
    w = t1.getcol('WEIGHT_SPECTRUM')
    print w.shape
    # Select time, weights and elevation after averaging over all baselines (axis 0).
    t = taql('SELECT TIME, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT, MSCAL.AZEL1()[1] AS ELEV FROM $t1 GROUPBY TIME')
    #antennas = t.getcol('ANTENNA')
    #print len(antennas), ' antennas'
    #antenna_names = taql('SELECT NAME FROM '+msfile+'/ANTENNA')

    weights = t.getcol('WEIGHT')
    print weights.shape
    time = t.getcol('TIME')
    elevation = t.getcol('ELEV')
    # Select the weights for all timestamps one channel and one polarization (in that order).
    weights = t.getcol('WEIGHT')[:, 0, 0]
    print weights.shape
    
    # Plot weights for elevation and save the image.
    print 'Plotting...'
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.scatter(time, weights, marker='.', color='k', label='Weights')
    # Plot the elevation as a function of time.
    ax_elev = ax.twinx()
    ax_elev.plot(time, elevation * 180/np.pi, 'b', label='Elevation')

    ax.set_xlim(min(time), max(time))
    ax.set_ylim(np.min(weights), np.max(weights))
    ax.set_xlabel('Time [s]')
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
    plot_weight_channel(msfile)
    #plot_weight_time(msfile)
