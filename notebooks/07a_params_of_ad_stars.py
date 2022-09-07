"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script takes and combines the information 
from the NASA Exoplanet Archive
and the TESS table on the star-planet systems
for which AD tests were performed."""

import pandas as pd
import numpy as np

if __name__ == "__main__":

    # Read in known parameters for Kepler and TESS SPSs
    path_kepler = "../data/2022_07_27_confirmed_uncontroversial_innermost_transiting_exoplanet.csv"
    kepler_params = pd.read_csv(path_kepler)

    path_tess = "../data/2022_07_27_tess_confirmed_and_known_planets.csv"
    tess_params = pd.read_csv(path_tess)

    # read in AD test results
    adtests = pd.read_csv("../results/multistar_adtest.csv")

    # get the unique TIC IDs
    unique_ids = adtests.TIC.unique()

    # find the TIC IDs that are in the AD test results in the stellar parameters table
    idx_kepler = np.where(kepler_params.tic_id.str[4:].isin(unique_ids))
    idx_tess = np.where(tess_params.tic_id.astype(str).isin(unique_ids))

    # select the rows with the stars for which we have AD test results
    sps_w_ad_kepler = kepler_params.iloc[idx_kepler]
    sps_w_ad_tess = tess_params.iloc[idx_tess]

    # rename the the TIC prefix from the TIC ID column
    sps_w_ad_kepler.tic_id = sps_w_ad_kepler.tic_id.str[4:]

    # convert the TIC ID column to string
    sps_w_ad_tess.tic_id = sps_w_ad_tess.tic_id.astype(str)

    # merge the two tables
    sps_w_ad = sps_w_ad_kepler.merge(sps_w_ad_tess, on="tic_id", how="outer",
                                    suffixes=("_kepler", "_tess"))

    # write the table to a CSV file
    path_to_params = "../results/params_of_star_planet_systems_with_AD_tests.csv"
    sps_w_ad.to_csv(path_to_params, index=False)

    print(f"Wrote {path_to_params}")
    print(f"{len(sps_w_ad)} rows")
