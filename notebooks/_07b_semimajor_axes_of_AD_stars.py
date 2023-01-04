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

from funcs.keplerian import get_semimajor_axis, get_distance_from_planet_range

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

    # make tic_id a string
    sps_w_ad.tic_id = sps_w_ad.tic_id.astype(str)
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
    # sps_w_ad.loc[sps_w_ad.pl_orbsmax.isna(), 
    #             "pl_orbsmax"] = a_is_none.apply(wrap_get_semimajor_axis, axis=1)


    # -------------------------------------------------------------------------
    # Add some values manually

    # HIP 67522 from rizzuto+2020 11.7 stellar radii
    # with 1.38 R sun, and 0.00465047 AU per solar radius
    au_per_rad = 0.00465047
    strad = sps_w_ad.loc[sps_w_ad.tic_id == "166527623", "st_rad_kepler"]
    sps_w_ad.loc[sps_w_ad.tic_id == "166527623", "pl_orbsmax"] = 11.7 * au_per_rad * strad
    sps_w_ad.loc[sps_w_ad.tic_id == "166527623", "pl_orbsmaxerr2"] = -0.27 * au_per_rad * strad
    sps_w_ad.loc[sps_w_ad.tic_id == "166527623", "pl_orbsmaxerr1"] = 0.24 * au_per_rad * strad

    reflink = ("<a refstr=RIZZUTO_ET_AL__2020 "
               "href=https://ui.adsabs.harvard.edu/abs/2020AJ....160...33R/abstract "
               "target=ref>Rizzuto et al. 2020</a>")
    sps_w_ad.loc[sps_w_ad.tic_id == "166527623", "pl_orbsmax_reflink"] = reflink


    # TOI 837 from Bouma+2020 17.26 stellar radii
    # with 1.04 R sun, and 0.00465047 AU per solar radius
    strad = sps_w_ad.loc[sps_w_ad.tic_id == "460205581", "st_rad_kepler"]
    sps_w_ad.loc[sps_w_ad.tic_id == "460205581", "pl_orbsmax"] = 17.26 * au_per_rad * strad
    sps_w_ad.loc[sps_w_ad.tic_id == "460205581", "pl_orbsmaxerr2"] = -0.6 * au_per_rad * strad
    sps_w_ad.loc[sps_w_ad.tic_id == "460205581", "pl_orbsmaxerr1"] = 0.6 * au_per_rad * strad

    reflink = ("<a refstr=BOUMA_ET_AL__2020 "
               "href=https://ui.adsabs.harvard.edu/abs/2020AJ....160..239B/abstract "
               "target=ref>Bouma et al. 2020</a>")

    # TOI 837 eccentricity from Bouma+2020 is not constrained
    # so put NaN
    sps_w_ad.loc[sps_w_ad.tic_id == "460205581", "pl_orbeccen"] = np.nan


    # K2-354 from DeLeon+2021 19.05 stellar radii
    # with 0.41 R sun, and 0.00465047 AU per solar radius
    strad = sps_w_ad.loc[sps_w_ad.tic_id == "468989066", "st_rad_kepler"]
    sps_w_ad.loc[sps_w_ad.tic_id == "468989066", "pl_orbsmax"] = 19.05 * au_per_rad * strad
    sps_w_ad.loc[sps_w_ad.tic_id == "468989066", "pl_orbsmaxerr2"] = -0.21 * au_per_rad * strad
    sps_w_ad.loc[sps_w_ad.tic_id == "468989066", "pl_orbsmaxerr1"] = 0.21 * au_per_rad * strad

    reflink = ("<a refstr=DE_LEON_ET_AL__2021 "
               "href=https://ui.adsabs.harvard.edu/abs/2021MNRAS.tmp.2120D/abstract "
               "target=ref>de Leon et al. 2021</a>")
    sps_w_ad.loc[sps_w_ad.tic_id == "468989066", "pl_orbsmax_reflink"] = reflink


    # TAP 26 eccentricity from Yu+2017 0.16 =/- 0.15
    sps_w_ad.loc[sps_w_ad.tic_id == "435907158", "pl_orbeccen"] = 0.16
    sps_w_ad.loc[sps_w_ad.tic_id == "435907158", "pl_orbeccenerr1"] = 0.15
    sps_w_ad.loc[sps_w_ad.tic_id == "435907158", "pl_orbeccenerr2"] = -0.15

    reflink = ("<a refstr=YU_ET_AL__2017 "
               "href=https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.1342Y/abstract "
               "target=ref>Yu et al. 2017</a>")
    sps_w_ad.loc[sps_w_ad.tic_id == "435907158", "pl_orbeccen_reflink"] = reflink

    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # CALCULATE THE BROWN DWARF SEMI-MAJOR AXIS

    # For the remaining unknown semi-major axes, supplement literature values for
    # stellar mass, and planet mass

    # ------------------------------------------------------------------------------
    # TIC 67646988

    # from Irwin+2018 Table 8
    st_mass_67646988 = 0.30
    pl_bmassj_67646988 = 0.0650 * M_sun / M_jup

    pl_orbsmax_67646988 = get_semimajor_axis(1.881719, st_mass_67646988,
                                            pl_bmassj_67646988).to(u.AU).value

    sps_w_ad.loc[sps_w_ad.tic_id == "67646988","pl_orbsmax"] = pl_orbsmax_67646988

    # ------------------------------------------------------------------------------
    # TIC 236387002

    # from Canas+2022 Table 3
    pl_orbsmax_236387002 = 0.064

    sps_w_ad.loc[sps_w_ad.tic_id == "236387002","pl_orbsmax"] = pl_orbsmax_236387002

    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------

    remaining_missing = sps_w_ad[sps_w_ad.pl_orbsmax.isna()]
    print(f"{remaining_missing.shape[0]} SPSs still miss semi-major axis info.\n")
    if remaining_missing.shape[0] != 0:
        print(f"These are:\n")
        print(remaining_missing.tic_id)

    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # NOW COMPUTE THE RANGE OF THE SEMI-MAJOR AXIS USING ECCENTRICITY AND
    # THE ERROR ON THE SEMI-MAJOR AXIS

    # if eccentricity is given as 0 with no uncertainties, fill in np.nan
    sps_w_ad.loc[sps_w_ad.pl_orbeccen == 0, "pl_orbeccen"] = np.nan

    # if eccentricity is given as 0 with no a value >0 for the upper error,
    # fill in that value
    sps_w_ad.loc[(sps_w_ad.pl_orbeccen == 0) & 
                 (sps_w_ad.pl_orbeccenerr1 > 0.), "pl_orbeccen"] = sps_w_ad.loc[(sps_w_ad.pl_orbeccen == 0) & 
                               (sps_w_ad.pl_orbeccenerr1 > 0.),
                                "pl_orbeccenerr1"]

    sps_w_ad["a_au_err"] = sps_w_ad.apply(lambda x: 
                                            get_distance_from_planet_range(x.pl_orbsmax, x.pl_orbeccen, x.pl_orbsmaxerr1), axis=1)

    # write the result to the original file
    sps_w_ad.to_csv("../results/params_of_star_planet_systems_with_AD_tests.csv", index=False)