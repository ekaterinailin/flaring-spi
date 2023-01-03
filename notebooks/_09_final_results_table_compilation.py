"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script calculates a bunch of values about the AD tests and stellar systems,
 including:

- the magnetic field of the star based on X-ray luminosity (Reiners+2022)
- the total observing time per star
- number of covered orbits per star
- relative velocity between star and planet at the planets orbit
- SPI power (Lanya 2009 and Lanza 2012)

Final table is called: results.csv

* https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PSCompPars
** https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
"""


import numpy as np
import pandas as pd

from funcs.ad import aggregate_pvalues
from funcs.spirelations import (b_from_lx_reiners,
                                wrap_obstimes,
                                p_spi_lanza12,
                                calculate_relative_velocity)


def map_bibkey(reflink, bibkeys):
    """Map the bibkey to the reflink.
    
    Parameters
    ----------
    reflink : str
        The reflink to the source.
    bibkeys : pd.DataFrame 
        The table with the reflink to bibkey mapping.

    Returns
    -------
    str
        The bibkey.
    """

    try:
        return bibkeys[bibkeys.pl_orbper_reflink == reflink]["pl_orbper_bibkey"].values[0]
    except:
        return np.nan


if __name__ == "__main__":

    # -------------------------------------------------------------------------
    # First, read in all the data you want to combine, and will need for the 
    # calculations

    # [READ] 
    # in AD test results
    adtests = pd.read_csv("../results/multistar_adtest.csv")

    # [READ] 
    # in the table of known parameters from the NASA Exoplanet Archive and TESS
    sps_w_ad = pd.read_csv(
        "../results/params_of_star_planet_systems_with_AD_tests.csv")

    # rename column tic_id to TIC
    sps_w_ad.rename(columns={"tic_id": "TIC"}, inplace=True)

    # [READ] 
    # in age and rotation period from the literature search
    literature_params = pd.read_csv(
        "../results/2022_08_stellar_params.csv")

    # [READ] 
    # in the flare table
    flares = pd.read_csv("../results/2022_07_flares_vetted.csv")

    # -------------------------------------------------------------------------
    # Next, initialize the final table by aggregating the p-values for the AD tests

    # aggregate the p-values
    mean_std = aggregate_pvalues(adtests, subsample="ED>1s", period="orbit")

    # Then merge in the NASA Exoplanet Archive table and the literature 
    # search table

    # merge the AD test results with the stellar parameters
    mean_std = mean_std.merge(literature_params, on=["TIC"], how="left")

    # add semi-major axis in AU, stellar radius in Rsun, 
    # planet radius in Rjup, and orbital period to the table:

    # merge the relevant part of the table
    mean_std = mean_std.merge(sps_w_ad[["TIC", "pl_orbsmax","st_rad_kepler",
                                        "st_rad_tess", "pl_radj", "pl_orbper_kepler",
                                        "pl_orbper_tess","pl_orbper_reflink"]],
                                on="TIC", how="left")

    # rename columns and fill NaNs will TESS values
    mean_std = mean_std.rename(columns = {"pl_orbsmax": "a_au",})
    mean_std["st_rad"] =  mean_std.st_rad_kepler.fillna(mean_std.st_rad_tess)
    mean_std["orbper_d"] =  mean_std.pl_orbper_kepler.fillna(mean_std.pl_orbper_tess)

    # -------------------------------------------------------------------------
    # Add some values manually

    # HIP 67522 from rizzuto+2020 56 stellar radii
    # with 1.38 R sun, and 0.00465047 AU per solar radius
    mean_std.loc[mean_std.TIC == 166527623, "a_au"] = 1.38 * 56 * 0.00465047


    # -------------------------------------------------------------------------
    # Now it's time to calculate some properties of the systems with
    # AD tests:

    # calculate B field from X-ray luminosity and stellar radius
    mean_std["B_G"] = mean_std.apply(lambda x: b_from_lx_reiners(x.xray_flux_erg_s, x.st_rad ), axis=1)

    # get the observing time for each star using the flare tables
    mean_std["obstime_d"] = mean_std.apply(lambda x: wrap_obstimes(str(x.TIC), 
                                                    flares), axis=1)

    # calculate the number of covered orbits
    mean_std["orbits_covered"] = mean_std["obstime_d"] / mean_std["orbper_d"]

    # calculate relative velocity between stellar rotation at the planetary orbit
    # and the orbital velocity of the planet in km/s
    mean_std["v_rel_km_s"] = mean_std.apply(lambda x: 
                                            calculate_relative_velocity(x.a_au,
                                            x.orbper_d, x.st_rotp), axis=1)


    # calculate the SPI power from the Lanza 2012 scaling relation
    mean_std["p_spi_erg_s"] = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                            x.B_G, x.pl_radj), axis=1)   

    # calculate the SPI power from the Lanza 2012 scaling relation with Bp=0
    mean_std["p_spi_erg_s_bp0"] = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                            x.B_G, x.pl_radj, Bp=0.), axis=1) 

    # -------------------------------------------------------------------------
    # For transparency, add bibkeys to the table for the literature values

    # 1. orbital period

    # read in the reflink to bibkey mapping table
    bibkeys = pd.read_csv("../results/pl_orbper_source_to_bibkey.csv")

    # add the bibkey to the table
    mean_std["pl_orbper_bibkey"] = mean_std.apply(lambda x: map_bibkey(x.pl_orbper_reflink, bibkeys),
                                                    axis=1)

    # delete the reflink column
    mean_std.drop(columns=["pl_orbper_reflink"], inplace=True)





    # -------------------------------------------------------------------------
    # Finally, save the table

    # [WRITE] 
    # to file
    mean_std.to_csv("../results/results.csv", index=False)

    # [WRITE]
    # to paper folder
    mean_std.to_csv("../../../002_writing/flaring-spi-paper/src/data/results.csv",
                    index=False)