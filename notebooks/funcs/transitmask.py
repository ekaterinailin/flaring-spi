"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2021, MIT License

Functions that use information about transiting planets
to mask transit signatures in FlareLightCurves.
"""

from .notebook import *

import transitleastsquares as tls

def get_transit_mid_epochs(tranmidepoch, tranmidepocherr, orbper, orbpererr, t0, tf):
    """Get transit mid epochs between t0 and tf for a
    planet with known orbital period and one recorded
    transit mid epoch.
    
    Parameters:
    -----------
    trandmidepoch : float
        a transit epoch in BKJD or BTJD
    orbper : float
        orbital period in days
    t0, tf : float, float
        start and end time of light curve
    
    Return:
    -------
    array of transit mid epochs between t0 and tf
    """
    # get number of periods between t0 and recorded transit
    dper = (t0 - tranmidepoch) // orbper #works for both >0 and <0
    
    # get the epoch of a transit right before t0
    tran0 = dper * orbper + tranmidepoch
    tran0err_sqr = (dper * orbpererr)**2 + tranmidepocherr**2
    
    tranmidtimes = np.arange(tran0, tf + orbper, orbper)
    tranmidtimeserr = np.sqrt((np.arange(len(tranmidtimes)) * orbpererr) **2 + tran0err_sqr)
    
    return list(zip(tranmidtimes, tranmidtimeserr))

def get_full_transit_mask(system, flc, pad=0):
    """Get full transit mask for a light curve
    with one or more transiting planets.
    
    The padding is NOT TESTED! -- but should be easy to test.
    
    Parameters:
    ------------
    system : DataFrame
        table of planets with columns
        "pl_trandur", "pl_tranmidepoch", "pl_orbper"
    flc :  FlareLightCurve
        light curve for which to get the mask
    pad : float
        add this buffer time before and after the
        transit to the mask
    
    Return:
    --------
    boolean array of light curve length with 1=transit, and
    0=no transit.
    """

    time = flc.time.value
    
    # start and end time of light curve
    t0, tf = time[0], time[-1]
    
    # mask single transits
    for i, row in system.iterrows():
        if (np.isnan(row.pl_orbper)) & (~np.isnan(row.pl_tranmidepoch)):
            system.loc[i,"pl_orbper"] = (tf - t0) * 1e6 
    
    # half transit duration in days
    system["durhalf"] = system["pl_trandur"] / 48. # convert to days and cut in half

    print(system.durhalf, "<<< durhalf")
    print(system.pl_tranmidepocherr, "<<< pl_tranmidepoch")
    print(system.pl_orbpererr, "<<< pl_orbpererr")
    # flag all time stamps withing transit duration time
    try:
        res = system.apply(lambda x:
                           np.sum(np.array([((time > (mid - x.durhalf - pad - x.pl_orbpererr - miderr)) &
                                             (time < (mid + x.durhalf + pad + x.pl_orbpererr + miderr)))
                                            for mid, miderr in get_transit_mid_epochs(x.pl_tranmidepoch,
                                                                              x.pl_tranmidepocherr,
                                                                              x.pl_orbper, 
                                                                              x.pl_orbpererr,
                                                                              t0, tf)]),
                                  axis=0), # <<< merge masks all transits in one planet into a single mask
                           axis=1).values
   
    # throw error if the planets are not properly characterized:
    except ValueError as e:
        print(e)
        raise ValueError("Some planets in your system are not transiting"
                         ", or lack necessary info about transit epoch"
                         ", period, and duration.")
    
    plt.plot(time, np.sum(res, axis=0))
    # merge masks for all planets into a single mask:
    return np.sum(res, axis=0).astype(bool)


def fit_and_remove_transits(time, flux, ID=None, mission=None, pl_orbper=None,  pl_orbpererr=None, **kwargs):
    """Fit and subtract transits of a planet from light curves
    using prior information about the star-planet system.
    
    Parameters:
    -----------
    time : n-array
        time array in units of days
    flux : n-array
        flux array
    mission : str
        "TESS" or "Kepler"
    pl_orbper : float
        orbital period of planet in units of days
    pl_orbpererr : float
        uncertainty on orbital period
        
    Return:
    -------
    n-array - flux with transit signal removed
    """

    # transit durations in TESS vary from 0.2 hours to over 30 h
    # at 1 min cadence this would be at least 12 datapoints per transit
    # which should be fine given that we have 
    # more constraints on the transits given anyways

    # So: Downsample to 1 min cadence if needed:
    dt = np.diff(time).mean() 
    one_min = 1. / 60. / 24.
    if dt < one_min:
        new_time, new_flux = tls.resample(time, flux / np.median(flux), factor=one_min/dt)
    else:
        new_time, new_flux = time, flux / np.median(flux)

    # Fetch stellar parameters from TESS and Kepler catalogs
    if mission is not None:
        
        if mission == "TESS":
            kwarg = {"TIC_ID":ID}
            
        elif mission == "Kepler":
            kwarg = {"KIC_ID":ID}

        ab, mass, mass_min, mass_max, radius, radius_min, radius_max = tls.catalog_info(**kwarg)
        
        r = radius
        rmin = radius - 3 * radius_min
        rmax = radius + 3 * radius_max
        m = mass
        mmin=mass - 3 * mass_min
        mmax=mass + 3 * mass_max
        
    else: # pick defaults
        ab = [0.4804, 0.1867] # (a G2V star in the Kepler bandpass), 
        m, mmin, mmax = 1, 0.08, 1.3
        r, rmin, rmax = 1, 0.13, 1.3 # roughly everything with a convection zone
        
    if (pl_orbper is None) | (np.isnan(pl_orbper)):
        print("No transit duration given.")
        n_transits_min = 2
        period_min, period_max =  0., np.inf
    else:
        print(new_time[-1], new_time[0], pl_orbper)
        n_transits_min = int((new_time[-1] - new_time[0]) // pl_orbper)
        if pl_orbpererr is None:
            print("No transit duration uncertainty given.")
            period_min, period_max = 0.75 * pl_orbper, 1.5 * pl_orbper
        else:
            print("Transit duration and uncertainty given.")
            period_min = pl_orbper -  3 * pl_orbpererr
            period_max = pl_orbper +  3 * pl_orbpererr 

    # Initialize the LST model
    model = tls.transitleastsquares(new_time, new_flux)
    print(new_time, new_flux, new_time.shape, new_flux.shape)
    print(period_min, period_max, n_transits_min, mmin, mmax, rmin, rmax, n_transits_min)
    # Run the search with system contraints from the catalogs
    results = model.power(u=ab, 
                          R_star=r,
                          R_star_min=rmin,
                          R_star_max=rmax,
                          M_star=m,
                          M_star_min=mmin,
                          M_star_max=mmax,
                          period_min=period_min,
                          period_max=period_max, 
                          n_transits_min = max(n_transits_min, 1),# this could be 1, 0 would throw an error
                          **kwargs) 
    
    # resample model light curve to original cadence
    if type(results.model_lightcurve_model) == float:
        print("Transit not recovered.")
        return flux, results
    else:

        model_flux = np.interp(time, results.model_lightcurve_time, results.model_lightcurve_model)
       
        # get the transit mask
        in_transit = tls.transit_mask(time, results.period, results.duration, results.T0)

        # get new flux array with transits subtracted
        notransit_flux = np.copy(flux)
        notransit_flux[in_transit] = notransit_flux[in_transit] /  model_flux[in_transit]

        return notransit_flux, results
