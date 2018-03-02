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

def normalize(x, nmin, nmax):
    normed = (x - nmin) / (nmax - nmin)
    return normed

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

def plot_weight_channel(msfile, pol=0, delta=10, threshold=1e5):
    # Polarization indices are 0, 1, 2, 3 = XX, YY, XY, YX, respectively.
    print 'Plotting weights vs. channels for %s' % (msfile,)
    # GMEANS (or MEANS(GAGGR(array), N) calculates the mean over the Nth axis.
    # Select only rows where not all data is flagged.
    t1 = taql('SELECT TIME, ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM FROM $msfile WHERE ALL(FLAG)==False')
    # Average over all baselines (axis 0) and group the resulting data by antenna.
    t = taql('SELECT TIME, ANTENNA, MEANS(GAGGR(REAL(DATA)), 0) AS DATA_REAL, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT FROM $t1 GROUPBY ANTENNA')
    w = t.getcol('WEIGHT')
    datar = t.getcol('DATA_REAL')
    # Calculate the variance in the visibilities over channels.
    print w.shape
    print np.any(np.isfinite(datar))
    print datar.shape
    #print datar[0,:,0]
    #print datar[0,:,1]
    #print datar[0,:,2]
    #print datar[0,:,3]
    
    # Select all antennas, all channels and polarization pol.
    print 'Calculating visibility variance.'
    variance = np.ones(shape=(w.shape[0], w.shape[1]//delta, w.shape[2]))
    # Subtract adjacent channels to eliminate physical signal.
    datar_shifted = np.roll(datar, -1, axis=1)
    datar -= datar_shifted
    print datar.shape
    print datar[:, :, pol]
    for i in xrange(datar.shape[1]//delta):
        # Take a frequency bin of delta channels.
        v = np.nanvar(datar[:,delta*i: delta*i+delta,:], axis=1, keepdims=True)
        if not np.any(v):
            variance[:,delta*i: delta*i+delta,:] = -np.inf
        else:
            variance[:,delta*i: delta*i+delta,:] = 1. / v
    weights = t.getcol('WEIGHT')[:, :, pol]
    #weights = np.where(weights < threshold, weights, 0)
    print weights.shape
    print weights[:, :]
    #print variance[:, :, pol]
    #print variance[:, :, pol].shape
    antennas = t.getcol('ANTENNA')
    print len(antennas), ' antennas'
    antenna_names = taql('SELECT NAME FROM '+msfile+'/ANTENNA')
    # Obtain channel frequencies in Hz.
    chan_freq = taql('SELECT CHAN_FREQ FROM '+msfile+'/SPECTRAL_WINDOW')
    # Select the first table, column CHAN_FREQ and convert to MHz.
    freq = chan_freq[0]['CHAN_FREQ'] * 1e-6
    #print 'Frequencies [MHz]: '
    #print freq
    #print freq[::delta]
    #print len(freq[::delta])
    # Plot the results.
    print 'Plotting weights...'
    imgname ='weight_chan_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+str(delta)+'.png'
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    weights = normalize(weights, np.min(weights), np.max(weights))
    print weights
    f = freq[::delta]
    print len(freq), weights.shape[1]
    for w, c, a in zip(weights, colors, antennas):
        #w = normalize(w, np.nanmin(w), np.nanmax(w))
        if len(freq) > weights.shape[1]:
            ax.scatter(freq[:-1], w, color=c, marker='.', label=antenna_names[a]['NAME'])
        elif len(freq) == weights.shape[1]:
            ax.scatter(freq, w, color=c, marker='.', label=antenna_names[a]['NAME'])
    ax.set_xlim(min(freq), max(freq))
    ax.set_ylim(np.min(weights), np.max(weights))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
    ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Weights')
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)

    print 'Plotting variance weights...'
    imgname ='var_weight_chan_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+'.png'
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    variance = variance[:, :, pol]
    variance = normalize(variance, np.nanmin(variance), np.nanmax(variance))
    print variance
    indices = ((np.asarray(range(0, variance.shape[1])) + 0.5) * delta).astype(int)
    f = freq[indices]
    for v, c, a in zip(variance, colors, antennas):
        if len(f) > variance.shape[1]:
            ax.scatter(f[:-1], v, color=c, marker='.', label=antenna_names[a]['NAME'])
        elif len(f) == variance.shape[1]:
            ax.scatter(f, v, color=c, marker='.', label=antenna_names[a]['NAME'])
    ax.set_xlim(min(freq), max(freq))
    ax.set_ylim(np.nanmin(variance), np.nanmax(variance))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
    ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Variance Normalized w.r.t. XX')
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)

    return 

def plot_weight_time(msfile, dt_elev=100, plot_time_unit='h', delta=10):
    print 'Plotting weights vs. time for %s' % (msfile,)
    imgname ='weight_time_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+'_'+str(delta)+'.png'
    # Select the time, weights and elevation of ANTENNA1 averaging over baselines/antennas.
    # Select only unflagged data.
    t1 = taql('select ANTENNA1, ANTENNA2, DATA, WEIGHT_SPECTRUM, TIME, FIELD_ID from $msfile where ALL(FLAG)==False')
    w = t1.getcol('WEIGHT_SPECTRUM')
    # Select time, weights and elevation.
    # This gets the average elevation w.r.t. to antenna 1, where the average is taken over all baselines.
    t = taql('SELECT TIME, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT, MEANS(GAGGR(MSCAL.AZEL1()[1]), 0) AS ELEV FROM $t1 GROUPBY TIME')
    t2 = taql('SELECT TIME, MEANS(GAGGR(REAL(DATA)),0) AS DATA_REAL FROM $t1 GROUPBY TIME')
    weights = t.getcol('WEIGHT')
    time = t.getcol('TIME')
    elevation = t.getcol('ELEV')

    datar = t2.getcol('DATA_REAL')
    # Calculate variance over 10 timestamp intervals. The last interval may be shorter.
    # datar shape is (timestamps, channels, polarizations)
    variance = np.ones(shape=(len(time)//delta, weights.shape[1], weights.shape[2]))
    datar_shifted = np.roll(datar, -1, axis=1)
    datar = datar - datar_shifted
    for i in xrange(len(time)//delta):
        v = np.var(datar[delta*i: delta*i+delta,:,:], axis=0)
        if not np.any(v):
            variance[i] = -1. / v
        else:
            variance[i] = 1. / v

    # Select the weights for all timestamps one channel and one polarization (in that order).
    weights = t.getcol('WEIGHT')
    # Plot weights for elevation and save the image.
    print 'Plotting...'
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)

    if plot_time_unit.lower() == 'h':
        time_h = time / 3600
        time_h = (time_h - math.floor(time_h[0]/24) * 24)
        
        # Deal with outlier weights if necessary.        
        wexponents = np.floor(np.log10(weights)).astype(int)
        print np.unique(wexponents)
        weights = np.where(wexponents < 0, weights, 1e-10)
        
        # Normalize the weights w.r.t. the XX polarization.
        nmin = np.min(weights[:, 5, 0])
        nmax = np.max(weights[:, 5, 0])
        
        ax.scatter(time_h, normalize(weights[:, 5, 0], nmin, nmax), marker='.', color='C1', alpha=0.25, label='XX Weights')
        ax.scatter(time_h, normalize(weights[:, 5, 1], nmin, nmax), marker='.', color='C2', alpha=0.25, label='XY Weights')
        ax.scatter(time_h, normalize(weights[:, 5, 2], nmin, nmax), marker='.', color='C3', alpha=0.25, label='YX Weights')
        ax.scatter(time_h, normalize(weights[:, 5, 3], nmin, nmax), marker='.', color='C4', alpha=0.25, label='YY Weights')
        del nmin, nmax
        
        # Normalize the statistic w.r.t. the XX polarization.
        nmin = np.min(variance[:, 5, 0])
        nmax = np.max(variance[:, 5, 0])
        indices = ((np.asarray(range(0, len(variance[:, 5, 0]))) + 0.5) * delta).astype(int)
        ax.plot(time_h[indices], normalize(variance[:, 5, 0], nmin, nmax), '--d', color='C1', label='XX Boxed variance $\\Delta=%d$'%(delta,))
        ax.plot(time_h[indices], normalize(variance[:, 5, 1], nmin, nmax), '--d', color='C2', label='XY Boxed variance $\\Delta=%d$'%(delta,))
        ax.plot(time_h[indices], normalize(variance[:, 5, 2], nmin, nmax), '--d', color='C3', label='YX Boxed variance $\\Delta=%d$'%(delta,))
        ax.plot(time_h[indices], normalize(variance[:, 5, 3], nmin, nmax), '--d', color='C4', label='YY Boxed variance $\\Delta=%d$'%(delta,))
        # Plot the elevation as a function of time.
        ax_elev = ax.twinx()
        ax_elev.plot(time_h, elevation * 180/np.pi, 'k', linewidth=2, label='Elevation')

        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%.2d:%.2d:%.2d' % (int(x), (x%1)*60, (((x%1)*60)%1 * 60))))
        ax.set_xlim(min(time_h), max(time_h))
        ax.set_ylim(0, 1.5)
        ax.set_xlabel('Time [h]')
    elif plot_time_unit.lower() == 's':
        # To do.
        pass
    ax.set_ylabel('Weights Normalized w.r.t. XX')
    ax_elev.set_ylabel('Elevation [deg]')

    # Deal with the legend.
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax_elev.get_legend_handles_labels()
    leg = ax.legend(handles+handles2, labels+labels2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    #show()
    return

if __name__ == '__main__':
    # Get the MS filename.
    msfile = sys.argv[1]
    #plot_weight_channel(msfile, delta=16)
    plot_weight_time(msfile, delta=100)
    #plot_data_channel(msfile, pol=3)
