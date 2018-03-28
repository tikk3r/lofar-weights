#!/usr/bin/env python
""" Provides acces to the weights of an MS file. """
from __future__ import division

import matplotlib.colors as colors
import matplotlib.patheffects as pe
import numpy as np

from casacore.tables import taql
from matplotlib import cm
from matplotlib.pyplot import figure, subplots
from matplotlib import ticker
from matplotlib.pyplot import show

__author__ = 'Frits Sweijen'
__credits__= 'Francesco de Gasperin'
__version__ = '1.0.0'
__maintainer__ = 'Frits Sweijen'
__email__ = 'sweijen <at> strw.leidenuniv.nl'
__status__ = 'Development'

def filter_IQR(x, threshold=3):
    Q1 = np.percentile(x, 25)
    Q3 = np.percentile(x, 75)
    IQR = Q3 - Q1
    if type(x) is np.ma.core.MaskedArray:
        med = np.ma.median(x)
    else:
        med = np.median(x)
    bound_upper = med + threshold * IQR
    bound_lower = med - threshold * IQR
    y = np.ma.masked_where(x > bound_upper, x)
    z = np.ma.masked_where(y < bound_lower, y)
    return z


def normalize(x, nmin, nmax):
    normed = (x - nmin) / (nmax - nmin)
    return normed

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

class WeightPlotter:
    def __init__(self, msfile):
        if msfile.endswith('/'):
            self.msfile = msfile[:-1]
        else:
            self.msfile = msfile

        self.colors = ['C1', 'C2', 'C3', 'C4']
        # Open the table and ignore rows that are completely flagged.
        print 'Loading '+msfile+'...'
        mstable_init = taql('SELECT TIME, ANTENNA1, FIELD_ID, CORRECTED_DATA, WEIGHT_SPECTRUM, FLAG FROM $msfile WHERE !ALL(FLAG)')
        self.mstable_init = mstable_init
        self.mstable = taql('SELECT TIME, MEANS(GAGGR(CORRECTED_DATA), 0) AS DATA, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHTS, MEANS(GAGGR(MSCAL.AZEL1()[1]), 0) AS ELEV, FLAG FROM $mstable_init GROUPBY TIME')

        # Additional processing before we can use the data.
        self.flags = self.mstable.getcol('FLAG')
        self.time = self.mstable.getcol('TIME')
        self.elevation = self.mstable.getcol('ELEV')

        data_um = self.mstable.getcol('DATA')
        self.data = np.ma.MaskedArray(data=data_um, mask=self.flags)

        weights_um = self.mstable.getcol('WEIGHTS')
        self.weights = np.ma.MaskedArray(data=weights_um, mask=self.flags)
        print 'FLAG table applied to DATA and WEIGHT_SPECTRUM.'
        # Obtain polarization setup.
        temp = taql('SELECT CORR_TYPE from '+msfile+'/POLARIZATION')
        if temp.getcol('CORR_TYPE')[0] in np.asarray([5, 6, 7, 8]):
            # Circular polarization.
            self.polarization = ['RR', 'LL', 'RL', 'LR']
        elif temp.getcol('CORR_TYPE')[0] in np.asarray([9, 10, 11, 12]):
            self.polarization = ['XX', 'YY', 'XY', 'YX']
        print 'Polarzation setup:', ', '.join(self.polarization)

        # Obtain channel frequencies.
        chan_freq = taql('SELECT CHAN_FREQ FROM '+msfile+'/SPECTRAL_WINDOW')
        # Select the first table, column CHAN_FREQ and convert to MHz.
        self.freq = chan_freq[0]['CHAN_FREQ'] * 1e-6

        print self.msfile+' loaded.'

    def plot_data_2D(self):
        print 'Plotting data matrix...'
        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Normalized log10(1 + |DATA|) '+self.msfile, fontweight='bold')
        axes = axes.ravel()

        for (pol, pol_label) in enumerate(self.polarization):
            axes[pol].set_title(pol_label)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel='Time')
            axes[pol].label_outer()
            #print np.min(np.mean(weights[:, :, pol], axis=0)), np.max(np.mean(weights[:, :, pol], axis=0))
            #weights = np.ma.masked_where(self.weights > 1e5, self.weights)
            data = np.log10(1 + np.abs(self.data))
            if pol == 0:
                nmin = np.nanmin(data[:, :, pol])
                nmax = np.nanmax(data[:, :, pol])
            nweights = normalize(data[:, :, pol], nmin, nmax)
            #print nmin, nmax
            #print nweights
            #im = axes[pol].imshow(nweights, extent=[freq[0], freq[-1], time[0], time[-1]])
            X, Y = np.meshgrid(self.freq, self.time - self.time[0])
            im = axes[pol].pcolor(X, Y, nweights, vmin=0.0, vmax=1.0)
            axes[pol].set_aspect('auto')
            axes[pol].label_outer()
            fig.colorbar(im, ax=axes[pol])
        imgname = 'data_2D.png'
        fig.savefig(imgname, bbox_inches='tight', dpi=250)
        print 'Saved plot '+imgname

    def plot_variance_2D(self, delta=3):
        print 'Plotting variance matrix...'
        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Inverse Variance '+self.msfile+', box =%dx%d'%(delta, delta), fontweight='bold')
        axes = axes.ravel()
        for (pol, pol_label) in enumerate(self.polarization):
            print 'Processing polarization '+pol_label
            axes[pol].set_title(pol_label)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel='Time')
            axes[pol].label_outer()
            # Calculate variance in the visibilities.
            var = np.abs(self.visibility_variance2D(self.data, self.weights, pol, delta=delta))
            var = np.ma.MaskedArray(data=var, mask=self.flags[:, :, pol])
            if pol == 0:
                nmin = np.nanmin(var)#; print nmin
                nmax = np.nanmax(var)#; print nmax
                #print np.median(var)
                #print np.percentile(var, 10)
                nmax = np.percentile(var, 99)
            nvar = normalize(var, nmin, nmax)
            #im = axes[pol].imshow(nvar, extent=[freq[0], freq[-1], time[0], time[-1]])
            #im = axes[pol].imshow(nvar, extent=[freq[0], freq[-1], time[0], time[-1]], norm=colors.LogNorm(vmin=1e-4, vmax=1.0), origin='lower')
            X, Y = np.meshgrid(self.freq, self.time - self.time[0])
            im = axes[pol].pcolor(X, Y, nvar, norm=colors.LogNorm(vmin=1e-4, vmax=1.0))
            axes[pol].set_aspect('auto')
            fig.colorbar(im, ax=axes[pol])
            #break
        imgname = 'variance_2D_box_%dx%d.png' % (delta, delta)
        fig.savefig(imgname, bbox_inches='tight', dpi=250)
        print 'Saved plot '+imgname

    def plot_weight_2D(self):
        print 'Plotting weight matrix...'
        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Normalized Weights '+self.msfile, fontweight='bold')
        axes = axes.ravel()

        for (pol, pol_label) in enumerate(self.polarization):
            axes[pol].set_title(pol_label)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel='Time')
            axes[pol].label_outer()
            weights = self.weights[:, :, pol]
            # Deal with outliers in the weights to avoid plotting issues/biases.
            Q1 = np.percentile(weights, 25)
            Q3 = np.percentile(weights, 75)
            IQR = Q3 - Q1
            wmedian = np.median(weights)
            # Treat weights that are further away than 3 times the IQR from the median as outliers.
            outlier_thresh_max = (wmedian + 5*IQR)
            outlier_thresh_min = (wmedian - 5*IQR)
            if np.max(weights) > outlier_thresh_max or np.min(weights) < outlier_thresh_min:
                fig.suptitle('Normalized Weights '+self.msfile+' (IQR filtered)', fontweight='bold')
                weights = np.ma.masked_where(weights > outlier_thresh_max, weights)
                weights = np.ma.masked_where(weights < outlier_thresh_min, weights)
            if pol == 0:
                nmin = np.nanmin(weights[:, :])
                nmax = np.nanmax(weights[:, :])
            nweights = normalize(weights[:, :], nmin, nmax)
            #print nmin, nmax
            #print nweights
            #im = axes[pol].imshow(nweights, extent=[freq[0], freq[-1], time[0], time[-1]])
            X, Y = np.meshgrid(self.freq, self.time - self.time[0])
            im = axes[pol].pcolor(X, Y, nweights, vmin=0.0, vmax=1.0)
            axes[pol].set_aspect('auto')
            axes[pol].label_outer()
            fig.colorbar(im, ax=axes[pol])
        imgname = 'weight_2D.png'
        fig.savefig(imgname, bbox_inches='tight', dpi=250)
        print 'Saved plot '+imgname

    def plot_weight_frequency(self, delta=10):
        print 'Plotting weights vs. frequency...'
        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Time-mean Weights ' + self.msfile + ' $\\Delta=%d$'%(delta,), fontweight='bold')
        axes = axes.ravel()
        handles = []
        labels = []

        for (pol, pol_label) in enumerate(self.polarization):
            data_shift = np.roll(self.data, -1, axis=0)
            data_sub = self.data - data_shift
            # Calculate the variance in frequency space.
            fdata = np.nanmean(data_sub[:, :, pol], axis=0)
            fvar = np.zeros(shape=fdata.shape[0] // delta)
            f_axis = np.empty(shape=fdata.shape[0] // delta)
            for i in xrange(fdata.shape[0] // delta):
                #vr = np.nanvar(fdata.real[delta * i:delta * (i+1)])
                #vi = np.nanvar(fdata.imag[delta * i:delta * (i+1)])
                #v = (vr + vi) / 2.
                v = np.nanvar(fdata[delta * i:delta * (i+1)])
                if v != 0 and np.isfinite(v):
                    fvar[i] = 1. / v
                f_axis[i] = np.mean(self.freq[delta * i:delta * (i+1)])
            nmin = np.nanmin(fvar)
            nmax = np.nanmax(fvar)
            nvar = normalize(fvar, nmin, nmax)
            nvar = filter_IQR(nvar, 3)
            # Deal with outliers in the weights to avoid plotting issues/biases.
            Q1 = np.percentile(self.weights, 25)
            Q3 = np.percentile(self.weights, 75)
            IQR = Q3 - Q1
            wmedian = np.median(self.weights)
            # Treat weights that are further away than 3 times the IQR from the median as outliers.
            outlier_thresh_high = (wmedian + 5*IQR)
            outlier_thresh_low = (wmedian - 5*IQR)
            if np.max(self.weights) > outlier_thresh_high or np.min(self.weights) < outlier_thresh_low:
                weights = np.ma.masked_where(self.weights > outlier_thresh_high, self.weights)
                weights = np.ma.masked_where(weights < outlier_thresh_low, weights)
                axes[pol].set_title(pol_label+' (IQR filtered)')
            else:
                weights = self.weights
                axes[pol].set_title(pol_label)
            fweights = np.mean(weights, axis=0)
            # Normalize w.r.t. XX or RR.
            nmin = np.nanmin(fweights[:, 0])
            nmax = np.nanmax(fweights[:, 0])
            nweights = normalize(fweights[:, pol], nmin, nmax)
            # Take mean in time to plot weights as a function of frequency.
            #tweights = np.mean(weights[:, :, :], axis=0)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel=self.polarization[0]+' Normalized Weight')
            axes[pol].label_outer()
            axes[pol].plot(self.freq, nweights, color=self.colors[pol], label=pol_label+' MS Weights')
            axes[pol].plot(f_axis, nvar, '--d', color=self.colors[pol], label=pol_label+' Inv. Variance', linewidth=2, path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])

            lhandle, llabel = axes[pol].get_legend_handles_labels()
            handles.append(lhandle)
            labels.append(llabel)
        leg = axes[3].legend(np.ravel(handles), np.ravel(labels), loc='upper center', bbox_to_anchor=(-0.1, -0.3), ncol=2, borderaxespad=0.0)
        imgname = 'weight_frequency_mean_%s.png' % (self.msfile,)
        fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
        print 'Saved plot '+imgname

        fig, axes = subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Time-median Weights ' + self.msfile + ', $\\Delta=%d$'%(delta,), fontweight='bold')
        axes = axes.ravel()
        handles = []
        labels = []
        for (pol, pol_label) in enumerate(self.polarization):
            # Deal with outliers.
            Q1 = np.percentile(self.weights, 25)
            Q3 = np.percentile(self.weights, 75)
            IQR = Q3 - Q1
            wmedian = np.median(self.weights)
            # Treat weights that are further away than 3 times the IQR from the median as outliers.
            outlier_thresh_high = (wmedian + 5*IQR)
            outlier_thresh_low = (wmedian - 5*IQR)
            if np.max(self.weights) > outlier_thresh_high or np.min(self.weights) < outlier_thresh_low:
                weights = np.ma.masked_where(self.weights > outlier_thresh_high, self.weights)
                weights = np.ma.masked_where(weights < outlier_thresh_low, weights)
                fweights = np.median(weights, axis=0)
                axes[pol].set_title(pol_label+' (IQR filtered)')
            else:
                fweights = np.median(self.weights, axis=0)
                axes[pol].set_title(pol_label)

            # Normalize w.r.t. XX or RR.
            nmin = np.nanmin(fweights[:, 0])
            nmax = np.nanmax(fweights[:, 0])
            nweights = normalize(fweights[:, pol], nmin, nmax)
            axes[pol].set(xlabel='Frequency [MHz]', ylabel=self.polarization[0]+' Normalized Weight')
            axes[pol].label_outer()
            axes[pol].plot(self.freq, nweights, color=self.colors[pol], label=pol_label+' MS Weights')
            axes[pol].plot(f_axis, nvar, '--d', color=self.colors[pol], label=pol_label+' Inv. Variance', linewidth=2, path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
            lhandle, llabel = axes[pol].get_legend_handles_labels()
            handles.append(lhandle)
            labels.append(llabel)
        leg = axes[3].legend(np.ravel(handles), np.ravel(labels), loc='upper center', bbox_to_anchor=(-0.1, -0.3), ncol=2, borderaxespad=0.0)
        imgname = 'weight_frequency_median_%s.png' % (self.msfile,)
        fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
        print 'Saved plot '+imgname

    def plot_weight_frequency_antenna(self, filtered=True):
        ''' Requires resorting the init table, so it probably takes a time of the order of initial run.
        '''
        print 'Plotting WEIGHT vs frequency per antenna.'
        tab = self.mstable_init
        t = taql('SELECT ANTENNA1 AS ANTENNA, CORRECTED_DATA as DATA, WEIGHT_SPECTRUM AS WEIGHTS, FLAG FROM $tab GROUPBY ANTENNA1')

        # Additional processing before we can use the data.
        flags = t.getcol('FLAG')
        data_um = t.getcol('DATA')
        data = np.ma.MaskedArray(data=data_um, mask=flags)

        weights_um = t.getcol('WEIGHTS')
        weights = np.ma.MaskedArray(data=weights_um, mask=flags)
        weights = filter_IQR(weights, threshold=5)
        print 'FLAG table applied to DATA and WEIGHT_SPECTRUM.'
        antennas = t.getcol('ANTENNA')
        print len(antennas), ' antennas'
        antenna_names = taql('SELECT NAME FROM '+msfile+'/ANTENNA')
        antenna_names = antenna_names.getcol('NAME')
        print data.shape
        for pol, pol_label in enumerate(self.polarization):
            colors = iter(cm.rainbow(np.linspace(0, 1, weights.shape[0])))
            fig = figure()
            fig.suptitle('WEIGHT vs frequency (IQR filtered) for %s'%(self.msfile,))
            ax = fig.add_subplot(111)
            for w, c, a in zip(weights[:, :, pol], colors, antennas):
                if filtered:
                    y = savitzky_golay(w, window_size=7, order=1)
                    ax.plot(self.freq, y, color=c, label=antenna_names[a])
                else:
                    ax.plot(self.freq, w, color=c, label=antenna_names[a])
            ax.set_xlim(min(self.freq), max(self.freq))
            ax.set_ylim(np.min(weights), np.max(weights))
            ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
            major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
            ax.xaxis.set_major_formatter(major_formatter)
            ax.set_xlabel('Frequency [MHz]')
            ax.set_ylabel('Normalized Weights')
            leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
            imgname = 'weight_frequency_antenna_%s_%s.png' % (self.msfile, pol_label)
            print 'Saving plot as %s' % (imgname,)
            fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)

    def plot_weight_time(self, mode='mean', delta=100):
        print 'Plotting weights vs. time...'
        fig, axes = subplots(nrows=1, ncols=1)
        fig.suptitle('Weights ('+mode+') for '+self.msfile+', $\\Delta=%d$'%(delta,), fontweight='bold')
        axes.set(xlabel='Time', ylabel=self.polarization[0]+' Normalized Weight')
        axes.label_outer()
        axes_elev = axes.twinx()
        axes_elev.set(ylabel='Elevation [deg]')

        # Subtract adjacent channels to remove any signal. This should leave us with noise.
        data_shift = np.roll(self.data, -1, axis=1)
        data_sub = self.data - data_shift

        # Take the mean or median over frequency to marginalize into the time domain.
        data_sub = np.ma.masked_where(~np.isfinite(data_sub), data_sub)
        if mode == 'median':
            tdata = np.ma.median(data_sub, axis=1)
        elif mode == 'mean':
            tdata = np.ma.mean(data_sub, axis=1)
        else:
            print 'Unknown mode, using median.'
            tdata = np.ma.median(data_sub, axis=1)
        variance = np.zeros(shape=(tdata.shape[0] // delta, tdata.shape[1]))

        tdatar = tdata.real
        tdatai = tdata.imag
        tdataf = np.abs(tdata)
        t_axis = []
        for i in xrange(tdata.shape[0] // delta):
            t = np.mean(self.time[delta * i: delta * (i+1)])
            t_axis.append(t)
            #vr = np.nanvar(tdatar[delta * i: delta * (i+1), :], axis=0)
            #vi = np.nanvar(tdatai[delta * i: delta * (i+1), :], axis=0)
            #v = (vr + vi) / 2.
            v = np.nanvar(tdataf[delta * i: delta * (i+1), :], axis=0)
            variance[i] = np.where(np.isfinite(1. / v), 1. / v, np.nan)

        for (pol, pol_label) in enumerate(self.polarization):
            Q1 = np.percentile(self.weights[:, :, pol], 25)
            Q3 = np.percentile(self.weights[:, :, pol], 75)
            IQR = Q3 - Q1
            wmedian = np.median(self.weights[:, :, pol])
            outlier_thresh_high = (wmedian + 5*IQR)
            outlier_thresh_low = (wmedian - 5*IQR)
            if np.max(self.weights[:, :, pol]) > outlier_thresh_high or np.min(self.weights[:, :, pol]) < outlier_thresh_low:
                weights = np.ma.masked_where(self.weights[:, :, pol] > outlier_thresh_high, self.weights[:, :, pol])
                weights = np.ma.masked_where(weights < outlier_thresh_low, weights)
                fig.suptitle('Weights ('+mode+') for '+self.msfile+', $\\Delta=%d$ (IQR filtered)'%(delta,), fontweight='bold')
            else:
                weights = self.weights[:, :, pol]

            # Take median or mean in frequency to plot weights as a function of time.
            if mode == 'median':
                tweights = np.median(weights, axis=1)
            elif mode == 'mean':
                tweights = np.mean(weights, axis=1)
            else:
                print 'Unknown mode, using median.'
                tweights = np.median(weights, axis=1)
            #if mask_threshold is not None:
                # Deal with some extreme outliers that screw with the median in this case.
            tweights = np.ma.masked_where(tweights > 1e3, tweights)
            # Normalize w.r.t. XX or RR.
            if pol == 0:
                nmin = np.nanmin(tweights); nmax = np.nanmax(tweights)
            nweights = normalize(tweights, nmin, nmax)
            var = filter_IQR(variance[:, pol])
            nvar = normalize(var, np.nanmin(var), np.nanmax(var))


            axes.plot(self.time, nweights, label=pol_label, alpha=0.5, color=self.colors[pol])
            axes.plot(t_axis, nvar, '--d', label='Variance '+pol_label, color=self.colors[pol])
        axes_elev.plot(self.time, self.elevation * 180/np.pi, '-k', linewidth=2, label='Elevation')
        handles, labels = axes.get_legend_handles_labels()
        handles2, labels2 = axes_elev.get_legend_handles_labels()
        leg = axes.legend(handles+handles2, labels+labels2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, borderaxespad=0.0)
        imgname = 'weight_'+mode+'_time_%s.png' % (self.msfile,)
        fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
        print 'Saved plot '+imgname

    def visibility_variance2D(self, vis, msweights, pol=None, delta=3):
        if delta % 2 == 0:
            raise Exception('Box size must be odd!')
        s = delta // 2
        if pol is not None:
            vis = vis[:, :, pol]
            weight = np.copy(msweights[:, :, pol])
            msweights = msweights[:, :, pol]
        max_x = vis.shape[1] - 1
        max_y = vis.shape[0] - 1

        for (y, x), _ in np.ndenumerate(vis):
            #if x == 0:
            #    print 'At row', y
            # Horizontal bounds.
            if (x - s) <= 0:
                bound_left = 0
            else:
                bound_left = (x - s)
            if (x + s) >= max_x:
                bound_right = max_x
            else:
                bound_right = (x + s) + 1
            # Vertical bounds.
            if (y - s) <= 0:
                bound_low = 0
            else:
                bound_low = (y - s)
            if (y + s) >= max_y:
                bound_high = max_y
            else:
                bound_high = (y + s) + 1
            box_vis = vis[bound_low:bound_high, bound_left:bound_right]
            if np.any(box_vis) and not np.all(box_vis.mask):
                v = np.nanvar(box_vis)
                if v != 0.0:
                    w = 1. / v
                else:
                    w = 0.
            else:
                w = 0.
            if np.isfinite(w):
                weight[y, x] = w
        return weight

if __name__ == '__main__':
    import sys
    # Get the MS filename.
    msfile = sys.argv[1]
    wp = WeightPlotter(msfile)
    #wp.plot_data_2D()
<<<<<<< HEAD
    #wp.plot_weight_time(mode='mean', delta=100)
    #wp.plot_weight_time(mode='median', delta=100)
    wp.plot_weight_frequency(delta=20)
    #wp.plot_weight_frequency_antenna()
    #wp.plot_weight_2D()
    #wp.plot_variance_2D(delta=3)
=======
    wp.plot_weight_time(mode='mean', delta=100)
    wp.plot_weight_time(mode='median', delta=100)
    wp.plot_weight_frequency(delta=20)
    wp.plot_weight_frequency_antenna()
    wp.plot_weight_2D()
    wp.plot_variance_2D(delta=3)
>>>>>>> f26b3e8e9cd0d08a0395e13ba419a282ebb7950c
    #wp.plot_variance_2D(delta=7)
    #wp.plot_variance_2D(delta=9)

    #from matplotlib.pyplot import show
    #show()
