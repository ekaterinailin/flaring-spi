import numpy as np
import pandas as pd
import pytest

import batman
import transitleastsquares as tls

from altaipony.flarelc import FlareLightCurve

from ..transitmask import (get_full_transit_mask,
                           get_transit_mid_epochs,
                           fit_and_remove_transits
                          )



def test_get_full_transit_mask():
    flc = FlareLightCurve(time=np.linspace(10,30,101))

    system = pd.DataFrame({"pl_trandur":[12.,20.,np.nan],
                           "pl_tranmidepoch" : [9., 37., np.nan], 
                           "pl_orbper": [5.,10.,19.],
                           "pl_tranmidepocherr" : [0.001, 0.001, .1], 
                           "pl_orbpererr": [.001,.001,.1]})

    # if the table contains NaN, demand cleanup of system parameters
    error_message = ("Some planets in your system are not transiting"
                             ", or lack necessary info about transit epoch"
                             ", period, and duration.")
    with pytest.raises(ValueError, match=error_message):
        get_full_transit_mask(system, flc)

    # drop last row everything should work again    
    system = system.iloc[:2]

    # the results are manually verified here:
    mask = get_full_transit_mask(system, flc)
   
    assert (flc.time.value[mask] == 
            pytest.approx(np.array([13.8, 14. , 14.2, 16.6, 16.8, 17. , 
                                    17.2, 17.4, 18.8, 19. , 19.2, 23.8, 
                                    24. , 24.2, 26.6, 26.8, 27. , 27.2, 
                                    27.4, 28.8, 29. , 29.2])))
    

    assert mask.dtype == "bool"

    # Check overlapping transits are treated okay
    system = pd.DataFrame({"pl_trandur":[12., 20.],
                           "pl_tranmidepoch" : [9., 9.], 
                           "pl_orbper": [5.,10.],
                           "pl_tranmidepocherr" : [0.001, 0.001], 
                           "pl_orbpererr": [.001,.001]})

    mask = get_full_transit_mask(system, flc)
    

    # These are also manually confirmed
    assert (flc.time.value[mask] == 
            pytest.approx(np.array([13.8, 14. , 14.2, 18.6, 18.8,
                                    19. , 19.2, 19.4, 23.8, 24. ,
                                    24.2, 28.6, 28.8, 29. , 29.2, 
                                    29.4])))
    
    # Check that single transits are treated okay
    system = pd.DataFrame({"pl_trandur":[12., 20.],
                           "pl_tranmidepoch" : [9., 9.], 
                           "pl_orbper": [5.,np.nan],
                           "pl_tranmidepocherr" : [0.001, 0.001], 
                           "pl_orbpererr": [.001,.001]})

    mask = get_full_transit_mask(system, flc)


    # These are also manually confirmed
    assert (flc.time.value[mask] == 
            pytest.approx(np.array([13.8, 14. , 14.2, 18.8,
                                    19. , 19.2, 23.8, 24. ,
                                    24.2, 28.8, 29. , 29.2,])))

def test_get_transit_mid_epochs(): 
    assert (get_transit_mid_epochs(10., 0., 3., 0., 20., 40.) ==
            np.array([(19., 0.), (22., 0.), (25., 0.), (28., 0.), (31., 0.), (34., 0.), (37., 0.), (40., 0.)])).all()

    assert (get_transit_mid_epochs(46., 0., 3., 0., 20., 40.) ==
            np.array([(19., 0.), (22., 0.), (25., 0.), (28., 0.), (31., 0.), (34., 0.), (37., 0.), (40., 0.)])).all()

    with pytest.raises(ValueError) as e:
        get_transit_mid_epochs(46., .1, np.nan, .1, 20., 40.)

    with pytest.raises(ValueError) as e:
        get_transit_mid_epochs(np.nan, .1, 3., .1, 20., 40.)

    with pytest.raises(ValueError) as e:
        get_transit_mid_epochs(np.nan, .1, 3., .1, 20., np.nan)

    with pytest.raises(ValueError) as e:
        get_transit_mid_epochs(np.nan, .1, 3., .1, np.nan, 40.)
        
    # test a very long orbital period for single transit case
    assert (get_transit_mid_epochs(30., .0, 300., .0, 20., 40.) ==
        np.array([(-270., 0.), (30., 0.),  (330., 0.)])).all()
    
    
def test_fit_and_remove_transits():

    # --------------------------------------------------------
    # test case: transitleastsuqares ORIGINAL TEST reproduced
    # --------------------------------------------------------

    np.random.seed(seed=0)  # reproducibility 

    # Create test data
    time_start = 3.14
    data_duration = 100
    samples_per_day = 48
    samples = int(data_duration * samples_per_day)
    time = np.linspace(time_start, time_start + data_duration, samples)

    # Use batman to create transits
    ma = batman.TransitParams()
    ma.t0 = time_start  # time of inferior conjunction; first transit is X days after start
    ma.per = 10.123  # orbital period
    ma.rp = 6371 / 696342  # 6371 planet radius (in units of stellar radii)
    ma.a = 19  # semi-major axis (in units of stellar radii)
    ma.inc = 90  # orbital inclination (in degrees)
    ma.ecc = 0  # eccentricity
    ma.w = 90  # longitude of periastron (in degrees)
    ma.u = [0.4, 0.4]  # limb darkening coefficients
    ma.limb_dark = "quadratic"  # limb darkening model
    m = batman.TransitModel(ma, time)  # initializes model
    synthetic_signal = m.light_curve(ma)  # calculates light curve

    # Create noise and merge with flux
    ppm = 50  # Noise level in parts per million
    noise = np.random.normal(0, 10**-6 * ppm, int(samples))
    flux = synthetic_signal + noise


    # call the fit function
    trflux, results = fit_and_remove_transits(time, flux)

    # test transit parameters
    assert ma.rp == pytest.approx(results.rp_rs, rel=.05)
    assert results.period == pytest.approx(ma.per, rel=1e-3)
    assert results.duration == pytest.approx(0.171, rel=1e-2)
    assert results.T0 == pytest.approx(ma.per + time_start, rel=1e-2)

    in_transit = tls.transit_mask(time, results.period, results.duration, results.T0)

    assert (np.where(in_transit)[0] == 
            np.array([   0,    1,    2,    3,  482,  483,  484,  485,  486,  487,  488,
                       489,  968,  969,  970,  971,  972,  973,  974,  975, 1453, 1454,
                      1455, 1456, 1457, 1458, 1459, 1460, 1461, 1939, 1940, 1941, 1942,
                      1943, 1944, 1945, 1946, 1947, 2425, 2426, 2427, 2428, 2429, 2430,
                      2431, 2432, 2911, 2912, 2913, 2914, 2915, 2916, 2917, 2918, 3397,
                      3398, 3399, 3400, 3401, 3402, 3403, 3404, 3883, 3884, 3885, 3886,
                      3887, 3888, 3889, 3890, 4369, 4370, 4371, 4372, 4373, 4374, 4375,
                      4376])).all()

    # out of transit the flux is conserved
    assert ((trflux-flux)[~in_transit] < 1e-16).all()

    # -------------------------------
    # test case: SINGLE KNOWN TRANSIT
    # -------------------------------

    np.random.seed(seed=0) 
    ma.per = 70.123  # orbital period upped to get a single transit
    ma.t0 = time_start - 10 # shift the first transit such there is really only one transit left
    m = batman.TransitModel(ma, time)  # initializes model
    synthetic_signal = m.light_curve(ma)  # calculates light curve

    # Create noise and merge with flux
    ppm = 50  # Noise level in parts per million
    noise = np.random.normal(0, 10**-6 * ppm, int(samples))
    flux = synthetic_signal + noise

    # call the fit function
    trflux, results = fit_and_remove_transits(time, flux)

    # the period recovered is too short when we have a single transit
    assert results.period == pytest.approx(15.48, rel=1e-2)

    # call the fit function with constrained orbital period 
    # and uncertainty improves the result drastically
    trflux, results = fit_and_remove_transits(time, flux, pl_orbper=70., pl_orbpererr=2.)

    # test transit parameters
    assert ma.rp == pytest.approx(results.rp_rs, rel=.1)

    # single transits have trouble finding a correct period obviously for a single transit
    assert results.period == pytest.approx(71.75, rel=1e-2)

    # but the centering on the transit is fine
    assert results.T0 == pytest.approx(ma.per + ma.t0, rel=1e-2)

    in_transit = tls.transit_mask(time, results.period, results.duration, results.T0)

    # single transit recovered, yes
    assert (np.where(in_transit)[0] == np.arange(2865,2915)).all()

    # out of transit the flux is conserved
    assert ((trflux-flux)[~in_transit] < 1e-16).all()
