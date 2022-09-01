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