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




def fit_spline(flc, spline_coarseness=30, spline_order=3):
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
        f[0], f[-1] = flcp.flux[le], flcp.flux[rip-1]#np.nanmedian(flcp.flux[le:le+n]), np.nanmedian(flcp.flux[rip-1-n:rip-1])#
        
        # if the LC chunk is too short, interpolate linearly
        if t.shape[0] <= k:
            p1 = UnivariateSpline(t[1:-1], f[1:-1], k=1)
            flcp.detrended_flux[le:ri] = flcp.flux[le:ri] - p1(flcp.time[le:ri]) + flux_med
            
        # otherwise fit a spline
        else:
            p3 = UnivariateSpline(t, f, k=k)
            flcp.detrended_flux[le:ri] = flcp.flux[le:ri] - p3(flcp.time[le:ri]) + flux_med
            model[le:ri] = p3(flcp.time[le:ri])
    
    return flcp, model

