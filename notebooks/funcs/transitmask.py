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