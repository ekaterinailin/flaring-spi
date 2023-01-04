"""
UTF-8, Python 3

------------------
FlaringSPI
------------------

Ekaterina Ilin, 2022, MIT License


This module contains functions used to 
calculate 

- semi-major axis of a planet from given stellar mass,
  planet mass, and orbital period

  
"""

import numpy as np

from astropy import units as u
from astropy.constants import G, M_sun, M_jup



def get_semimajor_axis(period, mass_star, mass_planet):
    """
    Calculate the semimajor axis of a planet orbiting a star
    
    $G (M_star + M_planet) / (4 pi^2) = a^3 / T^2$.
    
    Parameters
    ----------
    period : float
        The orbital period of the planet in days.
    mass_star : float
        The mass of the star in solar masses.
    mass_planet : float
        The mass of the planet in Jupiter masses.
    
    Returns
    -------
    semimajor_axis : float
        The semimajor axis of the planet in AU.
    """
    # calculate the total mass
    total_mass = mass_star * M_sun + mass_planet * M_jup


    return (G * total_mass * (period * u.d)**2 / (4 * np.pi**2))**(1/3)




def get_distance_from_planet_range(semimajor_axis, eccentricity, semimajor_axis_err1):
    """If eccentricity is known, then the range of the distance is either 
    
    - the error on the semi-major axis or 
    - half of the difference between the semi-major and semi-minor axis, 
        as calculated from the eccentricity 
    
    -- whichever is larger. 
    
    If eccentricity is NOT known, then the rnage of the distance is either
    
    - the 25% error on the semi-major axis, that is, assuming e=0.5 or 
    - the uncertainty on the semimajor axis 
    
    -- whichever is larger.
    
    Parameters
    ----------
    semimajor_axis : float
        Semi-major axis in AU
    eccentricity : float
        Eccentricity
    semimajor_axis_err1 : float
        Upper error on the semi-major axis in AU

    Returns
    -------
    float
        Range of the semi-major axis in AU
    """

    if np.isnan(eccentricity):

        if ~np.isnan(semimajor_axis_err1):
            return np.max([semimajor_axis_err1, 0.25 * semimajor_axis])
        else:    
            return 0.25 * semimajor_axis

    else:
        # calculate the semi-minor axis
        semiminor_axis = semimajor_axis * np.sqrt(1 - eccentricity**2)

        # calculate the distance from the planet
        distance_range = semimajor_axis - semiminor_axis

        # return the maximum of distance range and error
        # if error is not nan
        if ~np.isnan(semimajor_axis_err1):
            return np.max([distance_range, semimajor_axis_err1])
        else:
            return distance_range
