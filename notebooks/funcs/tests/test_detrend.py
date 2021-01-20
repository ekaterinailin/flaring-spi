import pytest

import numpy as np

from altaipony.flarelc import FlareLightCurve

from ..detrend import estimate_detrended_noise


def test_estimate_detrended_noise():
    
    # setup light curve
    time = np.linspace(10,30,200)
    
    # seed numpy to get the same error array
    np.random.seed(30)
    
    # define flux with gaussian noise and baseline flux
    flux = np.random.normal(0,40, time.shape[0]) + 200.
    
    # define light curve
    flc = FlareLightCurve(time=time, detrended_flux=flux)

    # this should work
    flces = estimate_detrended_noise(flc, mask_pos_outliers_sigma=2.5, 
                                 std_window=100, padleft=3, padright=10)
    
    # error should be similar to input error of 40
    np.median(flces.detrended_flux_err) == pytest.approx(41.38048677022836)

    # re-seed and add a flare
    np.random.seed(30)
    flux = np.random.normal(0,40, time.shape[0]) + 200.
    flux[120:124] = [500,380,300,270]
    flc = FlareLightCurve(time=time, detrended_flux=flux)

    # should mask flare, error should not grow
    flces = estimate_detrended_noise(flc, mask_pos_outliers_sigma=2.5, 
                                 std_window=100, padleft=3, padright=10)

    np.median(flces.detrended_flux_err) == pytest.approx(41.24232394552432)

    # re-seed and add some NaNs
    np.random.seed(30)
    flux = np.random.normal(0,40, time.shape[0]) + 200.
    flux[120:124] = [500,380,300,270]
    flux[30:40] = np.nan
    flc = FlareLightCurve(time=time, detrended_flux=flux)

    # should work regardless
    flces = estimate_detrended_noise(flc, mask_pos_outliers_sigma=2.5, 
                                     std_window=100, padleft=3, padright=10)

    # error should not change too much
    np.median(flces.detrended_flux_err) == pytest.approx(41.23144256208637)