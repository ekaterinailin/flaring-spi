import copy
import numpy as np
import pandas as pd

from altaipony.altai import find_iterative_median


def estimate_detrended_noise(flcd, mask_pos_outliers_sigma=2.5, 
                             std_window=100, padleft=3, padright=10):


    # start with a first approximation to std
    flcd.detrended_flux_err[:] =  np.nanstd(flcd.detrended_flux)
    
    # remove nan
    flcd = flcd[np.isfinite(flcd.detrended_flux)]
    
    # and refine it:
    flcd = find_iterative_median(flcd)
    
    
    # make a copy first
    filtered = copy.deepcopy(flcd.detrended_flux)
    
    # get right bound of flux array
    tf = filtered.shape[0]

    # pick outliers
    excl = np.where(flcd.detrended_flux - flcd.it_med > 
                    mask_pos_outliers_sigma * flcd.detrended_flux_err)[0]

    # mask outliers with a padding to make sure flare flux does factor in
    exclpad = []
    ran = np.arange(-padleft,padright)

    # add padded indices to list
    [exclpad.append(x + i) for x in excl for i in ran ]
    exclpad = np.array(exclpad)

    # remove out of bounds indices
    print(exclpad)
    exclpad = exclpad[(exclpad > -1) & (exclpad < tf)]

    # mask strong positive outliers so that they don't add to std
    if len(exclpad) > 0:
        filtered[exclpad] = np.nan

    # apply rolling window std and interpolate the masked values
    flcd.detrended_flux_err[:] = pd.Series(filtered).rolling(std_window,
                                                             center=True,
                                                             min_periods=1).std().interpolate()
    return flcd
