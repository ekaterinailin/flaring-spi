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
    else:
        raise ValueError("Planet magnetic field must be >= 0.")


def b_from_lx_reiners(lx, r, error=False, lx_err=None, r_err=None):
    """Reiners et al. 2022 Eq. 4 inverted to get magnetic field in Gauss.
    
    Parameters
    ----------
    lx : float
        X-ray luminosity in erg/s.
    r : float
        Stellar radius in R_sun.
    error : bool
        If True, return the error in the magnetic field by
        estimating it from the scatter in the relation, and the
        error in the stellar radius and X-ray luminosity.
    lx_err : float
        Error in X-ray luminosity in erg/s.
    r_err : float
        Error in stellar radius in R_sun.

    Returns
    -------
    B : float
        Magnetic field in Gauss.
    B_err : float
        Error in magnetic field in Gauss.
    """
    if np.isnan(lx) | np.isnan(r) | (lx == 0) | (r == 0):
        B = np.nan
        highB, lowB = np.nan, np.nan
        if error:
            return B, highB, lowB
        else:
            return B
    else:
        # convert stellar radius to cm
        rcm = (r * R_sun).to(u.cm).value

        # constants from Reiners et al. 2022 Eq. 4    
        a, b = 1 / 3.28e-12, 1.58

        # error on b from Reiners et al. 2022 Eq. 4
        b_err = 0.06

        # calculate magnetic field from inverse of Reiners et al. 2022 Eq. 4
        get_B = lambda lx, rcm, a, b: (a * lx)**(1/b) / (4. * np.pi * rcm**2)

        B =  get_B(lx, rcm, a, b) 

        if error:

            # convert error in radius
            rcm_err = (r_err * R_sun).to(u.cm).value

            # calculate upper and lower 1-sigma range on magnetic field
            highB = get_B(lx+lx_err, rcm+rcm_err, a, b - b_err)
            lowB = get_B(lx-lx_err, rcm-rcm_err, a, b + b_err)

            return B, highB, lowB
        else:
            return B

    

def calculate_relative_velocity(a_au, orbper, rotper, error=False,
                                a_au_err=None, orbper_err=None, rotper_err=None):
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
    error : bool
        If True, return the error in the relative velocity from
        error propagation.
    a_au_err : float
        Error in semi-major axis in AU.
    orbper_err : float
        Error in orbital period in days.
    rotper_err : float
        Error in stellar rotation period in days.

    Returns
    -------
    v_rel_km_s : float
        Relative velocity between stellar rotation at the planetary orbit
        and the orbital velocity of the planet in km/s.
    """
    
    # return in km/s
    # 1AU = 149597870.700 km
    # 1day = 24h * 3600s

    # convet au to km
    a_au = a_au * 149597870.700

    if (np.isnan(a_au) | np.isnan(orbper) | np.isnan(rotper) | 
        (a_au == 0) | (orbper == 0) | (rotper == 0)):
        v_rel = np.nan
        v_rel_err = np.nan

    else:
        # convert orbper to s
        orbper = orbper * 24. * 3600.

        # convert rotper to s
        rotper = rotper * 24. * 3600.
        
        v_rel = 2 * np.pi * a_au * (1/orbper - 1/rotper)

        if error:
            # get error in km/s
            dv_dau = 2 * np.pi * (1/orbper - 1/rotper)
            dv_dorbper = -2 * np.pi * a_au * (1/orbper**2)
            dv_drotper = 2 * np.pi * a_au * (1/rotper**2)

            # convert a_au_err to km
            a_au_err = a_au_err * 149597870.700

            # convert orbper_err to s
            orbper_err = orbper_err * 24. * 3600.

            # convert rotper_err to s
            rotper_err = rotper_err * 24. * 3600.

            # quadratic error propagation
            v_rel_err = np.sqrt((dv_dau * a_au_err)**2 +
                                (dv_dorbper * orbper_err)**2 +
                                (dv_drotper * rotper_err)**2)

    if error:

        return v_rel, v_rel_err

    else:
        return v_rel

