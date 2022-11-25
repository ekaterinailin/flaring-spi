import numpy as np
import astropy.units as u
from astropy.constants import R_jup, R_sun

def p_spi_lanza12(v_rel, B, pl_rad, Bp=1.):
    """Lanza 2012 scaling relation
    
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
    return B**(4./3.) * v_rel * pl_rad**2 * Bp**(2./3.)



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
