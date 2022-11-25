"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script uses Kepler's third law to fill in missing semi-major axes for
several planets, where the NASA Exoplanet Archive has the orbital period,
stellar mass and planet mass given. It then also adds some literature values
for systems where any of these parameters are missing, so that Kepler III
cannot be applied, namely TIC 67646988 and TIC 236387002, so far.
"""

import numpy as np
import pandas as pd

from astropy import units as u
from astropy.constants import M_jup, M_sun

from funcs.keplerian import get_semimajor_axis

def wrap_get_semimajor_axis(x):
    """Wrap get_semimajor_axis to work with pandas.apply
    and take available period frm either kepler or TESS."""

    # If Kepler orbital period is available, prefer that
    if x.pl_orbper_kepler is not np.nan:
        p = x.pl_orbper_kepler
    else:
        p = x.pl_orbper_tess
        
    # Return semi-major axis
    return get_semimajor_axis(p, x.st_mass, x.pl_bmassj).to(u.AU).value
    
if __name__ == "__main__":

    # read in SPS table
    path = "../results/params_of_star_planet_systems_with_AD_tests.csv"
    sps_w_ad = pd.read_csv(path)
    print(f"\n\nReading SPS parameters from {path}.\n")

    # get stellar mass, and planet mass and orbital period
    mass_star = sps_w_ad.st_mass # in solar masses
    mass_planet = sps_w_ad.pl_bmassj # in Jupiter masses, best mass estimate either mass, or mass*sin(i)
    period = sps_w_ad.pl_orbper_kepler.fillna(sps_w_ad.pl_orbper_tess) # in days

    # calculate semimajor axis using stellar mass, and planet mass and orbital period
    semimajor_axis = get_semimajor_axis(period, mass_star, mass_planet)

    # select rows where semi-major axis is not known
    a_is_none = sps_w_ad[sps_w_ad.pl_orbsmax.isna()]

    print(f"Semi-major axis is missing for {a_is_none.shape[0]} SPSs.\n")

    # calculate semi-major axis for these rows
    sps_w_ad.loc[sps_w_ad.pl_orbsmax.isna(), 
                "pl_orbsmax"] = a_is_none.apply(wrap_get_semimajor_axis, axis=1)

    # For the remaining unknown semi-major axes, supplement literature values for
    # stellar mass, and planet mass

    # ------------------------------------------------------------------------------
    # TIC 67646988

    # from Irwin+2018 Table 8
    st_mass_67646988 = 0.30
    pl_bmassj_67646988 = 0.0650 * M_sun / M_jup

    pl_orbsmax_67646988 = get_semimajor_axis(1.881719, st_mass_67646988,
                                            pl_bmassj_67646988).to(u.AU).value

    sps_w_ad.loc[sps_w_ad.tic_id == 67646988,"pl_orbsmax"] = pl_orbsmax_67646988

    # ------------------------------------------------------------------------------
    # TIC 236387002

    # from Canas+2022 Table 3
    pl_orbsmax_236387002 = 0.064

    sps_w_ad.loc[sps_w_ad.tic_id == 236387002,"pl_orbsmax"] = pl_orbsmax_236387002

    # ------------------------------------------------------------------------------

    remaining_missing = sps_w_ad[sps_w_ad.pl_orbsmax.isna()]
    print(f"{remaining_missing.shape[0]} SPSs still miss semi-major axis info.\n")
    if remaining_missing.shape[0] != 0:
        print(f"These are:\n")
        print(remaining_missing.TIC)

    # write the result to the original file
    sps_w_ad.to_csv("../results/params_of_star_planet_systems_with_AD_tests.csv", index=False)