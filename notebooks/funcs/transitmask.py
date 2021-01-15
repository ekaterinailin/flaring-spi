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

def get_transit_mid_epochs(tranmidepoch, orbper, t0, tf):
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
    return np.arange(tran0, tf + orbper, orbper)

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
    # start and end time of light curve
    t0, tf = flc.time[0], flc.time[-1]
    
    # half transit duration in days
    system["durhalf"] = system["pl_trandur"] / 48. # convert to days and cut in half
    
    # flag all time stamps withing transit duration time
    try:
        res = system.apply(lambda x:
                           np.sum(np.array([((flc.time > (mid - x.durhalf - pad)) &
                                             (flc.time < (mid + x.durhalf + pad)))
                                            for mid in get_transit_mid_epochs(x.pl_tranmidepoch, 
                                                                              x.pl_orbper, 
                                                                              t0, tf)]),
                                  axis=0), # <<< merge masks all transits in one planet into a single mask
                           axis=1).values
        
    # throw error if the planets are not properly characterized:
    except ValueError:
        raise ValueError("Some planets in your system are not transiting"
                         ", or lack necessary info about transit epoch"
                         ", period, and duration.")
    
    # merge masks for all planets into a single mask:
    return np.sum(res, axis=0).astype(bool)