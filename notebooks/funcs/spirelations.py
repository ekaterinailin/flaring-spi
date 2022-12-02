import numpy as np
import astropy.units as u
from astropy.constants import R_jup, R_sun


def convert_datapoints_to_obstime(n, mission):
    """Calculate the observing time in days.
    based on the number of datapoints in the LC.

    Parameters
    ----------
    n : int
        Number of datapoints.
    mission : str
        Mission name, either "Kepler" or "TESS". 

    Returns
    -------
    obstime : float
        Observing time in days.
    """
    if n > 50000:
        if mission == "TESS":
            return n / 3. / 60. / 24.  # 20 sec cadence
        elif mission == "Kepler":
            return n / 60. / 24.  # 1 min cadence
        else:
            raise ValueError("Mission name not recognized for a large LC.")
    else:
        if mission == "TESS":
            return n * 2. / 60. / 24.  # 2 min cadence
        elif mission == "Kepler":
            return n / 60. / 24.  # 1 min cadence
        else:
            raise ValueError("Mission name not recognized for a small LC.")


def wrap_obstimes(TIC, flares):
    """Wrap the obstime function for use with multiple
    LCs for a given TIC ID. NOT TESTED!

    Parameters
    ----------
    TIC : str
        TIC ID.
    flares : pandas.DataFrame
        Flare table with total_n_valid_data_points and mission columns.

    Returns
    -------
    obstime_d: float
        Total observing time in days.
    """
    # pick up all the rows with LCs searched for flares
    sysvalues = flares[flares.TIC.astype(str) == TIC]

    # drop the duplicates with multiple flares per LC
    vals = sysvalues[["total_n_valid_data_points",
                      "mission"]].drop_duplicates()
  
    # calculate the observing time based on the 
    # number of data points in the LC
    # and summing up the observing time for all LCs
    obstime_d = vals.apply(lambda x: convert_datapoints_to_obstime(
        x.total_n_valid_data_points, x.mission), axis=1).sum()

    return obstime_d


def p_spi_lanza12(v_rel, B, pl_rad, Bp=1.):
    """Lanza 2012 scaling relation and its adaptation to absence of
    planet magnetic field (Bp=0).
    
    Parameters
    ----------
    v_rel : float
        Relative velocity between stellar rotation and orbital velocity in km/s.
    B : float
        X-ray luminosity in erg/s.
    pl_rad : float
        Planet radius in R_jup.
    Bp : float
        Planet magnetic field in Gauss.

    Returns
    -------
    p_spi : float
        Power of the star-planet interaction in erg/s
    """
    # convert from Rjup to cm
    pl_rad = pl_rad * R_jup.to(u.cm).value

    # convert from km/s to cm/s
    v_rel = v_rel * 1e6
    if Bp > 0.:
        return B**(4./3.) * v_rel * pl_rad**2 * Bp**(2./3.)
    elif Bp == 0.:
        return B**2 * v_rel * pl_rad**2



def b_from_lx_reiners(lx, r):
    """Reiners et al. 2022 Eq. 4 inverted to get magnetic field in Gauss.
    
    Parameters
    ----------
    lx : float
        X-ray luminosity in erg/s.
    r : float
        Stellar radius in R_sun.
    """
    rcm = (r * R_sun).to(u.cm).value
  
    B =  np.power(lx / 3.28e-12, 1. / 1.58) / (4. * np.pi * rcm**2) 
    return B


def calculate_relative_velocity(a_au, orbper, rotper):
    """Calculate the relative velocity between stellar rotation at the planetary orbit
    and the orbital distance of the planet in km/s.

    Parameters
    ----------
    au : float
        Semi-major axis in AU.
    orbper : float
        Orbital period in days.
    rotper : float
        Stellar rotation period in days.

    Returns
    -------
    v_rel_km_s : float
        Relative velocity between stellar rotation at the planetary orbit
        and the orbital velocity of the planet in km/s.
    """
    v_rel = 2 * np.pi * a_au * (1/orbper - 1/rotper) 
    
    # return in km/s
    # 1AU = 149597870.700 km
    # 1day = 24h * 3600s
    return v_rel / 24. / 3600. * 149597870.700