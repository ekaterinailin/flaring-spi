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

    # Read in known parameters for transiting Kepler SPSs
    path_kepler = "../data/2022_07_27_confirmed_uncontroversial_innermost_transiting_exoplanet.csv"
    kepler_params = pd.read_csv(path_kepler)

    # Read in known parameters for transiting TESS SPSs
    path_tess = "../data/2022_07_27_tess_confirmed_and_known_planets.csv"
    tess_params = pd.read_csv(path_tess)

    # Read in known parameters for non-transiting SPSs
    path_nontrans = "../data/2022_11_15_input_catalog_NONtransit_star_planet_systems.csv"
    rv_params = pd.read_csv(path_nontrans)

    # read in AD test results
    adtests = pd.read_csv("../results/multistar_adtest.csv")

    # get the unique TIC IDs
    unique_ids = adtests.TIC.astype(str).unique()

    # rename the the TIC prefix from the TIC ID column
    kepler_params.tic_id = kepler_params.tic_id.str[4:]

    # convert the TIC ID column to string
    tess_params.tic_id = tess_params.tic_id.astype(str)
    rv_params.tic_id = rv_params.tic_id.astype(str)

    # get the indices of AD SPSs in all the tables
    idx_kepler = kepler_params.tic_id.isin(unique_ids)
    idx_tess = tess_params.tic_id.isin(unique_ids)
    idx_rv = rv_params.tic_id.isin(unique_ids)

    # get the AD SPSs from all the tables
    sps_w_ad_kepler = kepler_params[idx_kepler]
    sps_w_ad_tess = tess_params[idx_tess]
    sps_w_ad_rv = rv_params[idx_rv]
    
    # concatenate the transiting and non-transiting SPSs in NASA tables
    sps_w_ad_transnontrans = pd.concat([sps_w_ad_kepler, sps_w_ad_rv])

    # merge in the TESS table
    sps_w_ad = sps_w_ad_transnontrans.merge(sps_w_ad_tess, on="tic_id", how="outer",
                                    suffixes=("_kepler", "_tess"))

    # add for AU Mic the values uncertainty in the orbital period 
    # from Martioli et al. 2021
    # and remove the Kepler enty
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", "pl_orbper_tess"] = 8.463
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", "pl_orbper_kepler"] = np.nan
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", "pl_orbpererr2_tess"] = 2e-6
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", "pl_orbpererr1_tess"] = 2e-6
    
    # add for Kepler-396 the values uncertainty in the orbital period
    # from Battley et al. 2022
    # the results are from TESS, but fill in Kepler anyways
    sps_w_ad.loc[sps_w_ad.tic_id == "27769688", "pl_orbper_kepler"] = 42.99292140
    sps_w_ad.loc[sps_w_ad.tic_id == "27769688", "pl_orbpererr2_kepler"] = 0.00002072
    sps_w_ad.loc[sps_w_ad.tic_id == "27769688", "pl_orbpererr1_kepler"] = 0.00002072



    # write the table to a CSV file
    path_to_params = "../results/params_of_star_planet_systems_with_AD_tests.csv"
    sps_w_ad.to_csv(path_to_params, index=False)

    print(f"Wrote {path_to_params}")
    print(f"{len(sps_w_ad)} rows")

    # make sure all SPS are covered
    assert len(sps_w_ad) == len(unique_ids), "Not all SPS are covered!"
