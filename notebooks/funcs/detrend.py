import copy
import numpy as np
import pandas as pd

from altaipony.altai import find_iterative_median
from altaipony.utils import sigma_clip

from collections import defaultdict

import astropy.units as u

from scipy.interpolate import UnivariateSpline
from scipy import optimize
from scipy.fftpack import fft


def estimate_detrended_noise(flc, mask_pos_outliers_sigma=2.5, 
                             std_window=100, padleft=3, padright=10):

    flcc = copy.deepcopy(flc)
    flcc = flcc.find_gaps()

    for (le, ri) in flc.gaps:

        flcd = copy.deepcopy(flcc[le:ri])
        mask = sigma_clip(flcd.detrended_flux, max_sigma=mask_pos_outliers_sigma, longdecay=2)

        flcd.detrended_flux[~mask] = np.nan
        # apply rolling window std and interpolate the masked values
        flcd.detrended_flux_err[:] = pd.Series(flcd.detrended_flux).rolling(std_window,
                                                                 center=True,
                                                                 min_periods=1).std().interpolate()
        
        # remove nan
       # flcd = flcd[np.isfinite(flcd.detrended_flux)]
        
        # and refine it:
        flcd = find_iterative_median(flcd)
        
        
        # make a copy first
        filtered = copy.deepcopy(flcd.detrended_flux)
        
        # get right bound of flux array
        tf = filtered.shape[0]

        # pick outliers
        mask = sigma_clip(filtered, max_sigma=mask_pos_outliers_sigma, longdecay=2)
       # excl = np.where(flcd.detrended_flux - flcd.it_med > 
       #                 mask_pos_outliers_sigma * flcd.detrended_flux_err)[0]

        # mask outliers with a padding to make sure flare flux does factor in
        #exclpad = []
        #ran = np.arange(-padleft,padright)

        # add padded indices to list
        #[exclpad.append(x + i) for x in excl for i in ran ]
        #exclpad = np.array(exclpad)

        # remove out of bounds indices
        #exclpad = exclpad[(exclpad > -1) & (exclpad < tf)]

        # mask strong positive outliers so that they don't add to std
        #if len(exclpad) > 0:
        #    filtered[exclpad] = np.nan
        filtered[~mask] = np.nan    

        # apply rolling window std and interpolate the masked values
        flcc.detrended_flux_err[le:ri]= pd.Series(filtered).rolling(std_window,
                                                                 center=True,
                                                                 min_periods=1).std().interpolate()
    return flcc



def custom_detrending(flc, pad=25):
    """Wrapper"""
    f = flc.flux[np.isfinite(flc.flux)]
   
    if np.abs(f[0]-f[-1])/np.median(f) > .2:
        print("Do a coarse spline interpolation to remove trends.")
        flc = fit_spline(flc, spline_coarseness=12)
        flc.flux[:] = flc.detrended_flux[:]
    
    # Iteratively remove fast sines with Periods of 0.1 to 2 day periods (the very fast rotators)
    flc = iteratively_remove_sines(flc.remove_nans()) # nans come from too short LC chunks in spline fit
    flc.flux[:] = flc.detrended_flux[:]
    
    
    # remove some rolling medians on a 10 hours time scale
  #  dt = np.nanmin(np.diff(flc.time))
  #  roll = int(np.rint(0.5 / dt))
  #  flc.flux[:] = flc.flux - pd.Series(flc.flux).rolling(roll, center=True, min_periods=1).median() + np.nanmedian(flc.flux)
    
    # Determine the window length for the SavGol filter for each continuous observation gap
    flc = find_iterative_median(flc)
    flc = flc.remove_nans()
    
    w, flc = search_gaps_for_window_length(flc)
    
    # Use lightkurve's SavGol filter while padding outliers with some data points around the outliers/flare candidates
    flc = flc.detrend("savgol", window_length=w, pad=pad)
    flc.flux[:] = flc.detrended_flux[:]

    # After filtering, always use a 2.5 hour window to remove the remaining 
    flcd = flc.detrend("savgol", window_length=75, pad=pad)

    return flcd

def search_gaps_for_window_length(flc):
    """Search continuous light curve chunks for
    appropriate window_length to apply to 
    SavGol filter.
    
    Parameters:
    ------------
    flc : FlareLightCurve
    
    Return:
    -------
    list of odd ints
    """
    flc = flc.remove_nans()
    flc = flc.find_gaps()
    wls = []
    for le,ri in flc.gaps:
    
        wls.append(select_window_length(flc.flux[le:ri]))
      
    
    return wls, flc


def select_window_length(flux):
    """Performs an FFT and defines a window
    length that is smaller than the most prominent
    frequency.
    
    Parameters:
    -----------
    flux : array
        
    Return:
    --------
    odd int
    """
    #normalize flux and FFT it:
    yf = fft(flux/np.nanmean(flux)-1.)
    
    maxfreq = len(yf) // 5
    minfreq = 1

    # choose window length
    w = np.rint(len(yf) / (minfreq + np.argmax(yf[minfreq:maxfreq])) / 3)

    # w must be odd
    if w%2==0:
        w += 1
        
    # if w is too large don't do it at all
    if w > len(yf) // 2:
        return None
    else:
        return int(max(w, 75))


def fit_spline(flc, spline_coarseness=36, spline_order=3):
    """Do a spline fit on a coarse sampling of data points.
    
    Parameters:
    ------------
    flc : FlareLightCurve
    
    spline_coarseness : int
 
    spline_order : int
        order of spline fit
        
    Return:
    --------
    FlareLightCurve with new flux attribute
    """
    flc = flc[np.where(np.isfinite(flc.flux))]
    flcp = copy.deepcopy(flc)

    flcp = flcp.find_gaps()
    flux_med = np.nanmedian(flcp.flux)
    n = int(np.rint(spline_coarseness/ 24 / (flcp.time[1] - flcp.time[0])))
    k = spline_order
    #do a first round
    model = np.full_like(flcp.flux, np.nan)
    for le, ri in flcp.gaps:

        rip = flcp.flux[le:ri].shape[0] + le
        t, f = np.zeros((rip - le)//n+2), np.zeros((rip - le)//n+2)
    
        t[1:-1] = np.mean(flcp.time[le:rip - (rip - le)%n].reshape((rip - le)//n, n), axis=1)
        f[1:-1] =  np.median(flcp.flux[le:rip - (rip - le)%n].reshape((rip - le)//n, n), axis=1)
        t[0], t[-1] = flcp.time[le], flcp.time[rip-1]
        f[0], f[-1] = np.nanmedian(flcp.flux[le:le+n]), np.nanmedian(flcp.flux[rip-1-n:rip-1])#flcp.flux[le], flcp.flux[rip-1]
        
        # if the LC chunk is too short, drop it
        if t.shape[0] <= k:
            p1 = UnivariateSpline(t[1:-1], f[1:-1], k=1)
            flcp.detrended_flux[le:ri] = flcp.flux[le:ri] - p1(flcp.time[le:ri]) + flux_med
        else:
            #cond11 = flcp.time[le:ri] < t[2]
            #cond3 = (flcp.time[le:ri] > t[2]) & (flcp.time[le:ri] < t[-2]) 
            #cond12 = flcp.time[le:ri] > t[-2]
            #p11 = UnivariateSpline(t[:2], f[:2], k=1)
            p3 = UnivariateSpline(t, f, k=1)
            #p12 = UnivariateSpline(t[-2:], f[-2:], k=1)
            #cond11 = flcp.time[le:ri] < t[3]
            #cond3 = (flcp.time[le:ri] > t[3]) & (flcp.time[le:ri] < t[-3]) 
            #cond12 = flcp.time[le:ri] > t[-3]
            #flcp.detrended_flux[le:ri][cond3] = flcp.flux[le:ri][cond3] - p3(flcp.time[le:ri][cond3]) + flux_med
            #flcp.detrended_flux[le:ri][cond11] = flcp.flux[le:ri][cond11] - p11(flcp.time[le:ri][cond11]) + flux_med
            #flcp.detrended_flux[le:ri][cond12] = flcp.flux[le:ri][cond12] - p12(flcp.time[le:ri][cond12]) + flux_med
            #model[le:ri][cond3] = p3(flcp.time[le:ri][cond3])
            #model[le:ri][cond11] = p11(flcp.time[le:ri][cond11])
            #model[le:ri][cond12] = p12(flcp.time[le:ri][cond12])
            flcp.detrended_flux[le:ri] = flcp.flux[le:ri] - p3(flcp.time[le:ri]) + flux_med
            model[le:ri] = p3(flcp.time[le:ri])
    
    return flcp, model


def iteratively_remove_sines(flcd, freq_unit=1/u.day, 
                             maximum_frequency=10, 
                             minimum_frequency=0.5):
    def cosine(x, a, b, c, d):
        return a * np.cos(b * x + c) + d

    snr = 3
    flct = copy.deepcopy(flcd)
    for le, ri in flct.find_gaps().gaps:
        print(le,ri)
        flc = copy.deepcopy(flct[le:ri])
        flc = find_iterative_median(flc)

        #mask flares
        mask = sigma_clip(flc.flux, max_sigma=3.5, longdecay=2)
        flc = flc[mask]
        # only remove sines if LC chunk is larger than one full period of the fastest frequency
        if flc.flux.shape[0] > 1 / maximum_frequency / np.nanmin(np.diff(flc.time)):
            pg = flc.remove_nans().to_periodogram(freq_unit=freq_unit,
                                              maximum_frequency=maximum_frequency,
                                              minimum_frequency=minimum_frequency)
            snr = pg.flatten().max_power
        else:
            snr = 0.
    #    print("Found peak in periodogram at ", pg.frequency_at_max_power)
    #    print("SNR at ", snr)
        n = 0
        while ((snr > 1.) & (n < 10)):
            print(pg.frequency_at_max_power.value)
            pg = flc.remove_nans().to_periodogram(freq_unit=freq_unit,
                                                  maximum_frequency=maximum_frequency,
                                                  minimum_frequency=minimum_frequency)
            
            cond = np.invert(np.isnan(flc.time)) & np.invert(np.isnan(flc.flux)) 
            p, p_cov = optimize.curve_fit(cosine, flc.time[cond], flc.flux[cond],
                                          p0=[np.nanstd(flc.flux),
                                          2*np.pi*pg.frequency_at_max_power.value,
                                          0, np.nanmean(flc.flux)])
            flc.flux = np.nanmean(flc.flux) + flc.flux-cosine(flc.time, p[0], p[1], p[2], p[3])
           # print(snr)
            snr = pg.flatten().max_power
            n += 1
          #  print(snr)

        flcd.detrended_flux[le:ri] = flc.flux
    return flcd
