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

import forecaster

from astropy.constants import R_earth, R_jup

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
    # from Zicher et al. 2022
    # and remove the Kepler enty
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", "pl_orbper_tess"] = 8.463
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", "pl_orbper_kepler"] = np.nan
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", "pl_orbpererr2_tess"] = -2e-6
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", "pl_orbpererr1_tess"] = 2e-6
    
    # fill in Zicher et al. 2022 for AU Mic instead of Cale et al. 2021
    sps_w_ad.loc[sps_w_ad.tic_id == "441420236", 
                 "pl_orbper_reflink"] = ("<a refstr=ZICHER_ET_AL__2022 "
                                         "href=https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.3060Z/abstract "
                                         "target=ref>Zicher et al. 2022</a>")


    # add for Kepler-396 the values uncertainty in the orbital period
    # from Battley et al. 2021
    # the results are from TESS, but fill in Kepler anyways
    sps_w_ad.loc[sps_w_ad.tic_id == "27769688", "pl_orbper_kepler"] = 42.99292140
    sps_w_ad.loc[sps_w_ad.tic_id == "27769688", "pl_orbpererr2_kepler"] = -0.00002072
    sps_w_ad.loc[sps_w_ad.tic_id == "27769688", "pl_orbpererr1_kepler"] = 0.00002072
    sps_w_ad.loc[sps_w_ad.tic_id == "27769688", 
                 "pl_orbper_reflink"] = ("<a refstr=BATTLEY_ET_AL__2021 "
                                         "href=https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.4092B/abstract "
                                         "target=ref>Battley et al. 2021</a>")


    # for HAT-P-11, fill in the orbital period error from
    # Bakos et al. 2010
    sps_w_ad.loc[sps_w_ad.tic_id == "28230919", "pl_orbpererr2_kepler"] = -7.1e-6
    sps_w_ad.loc[sps_w_ad.tic_id == "28230919", "pl_orbpererr1_kepler"] = 7.1e-6
    sps_w_ad.loc[sps_w_ad.tic_id == "28230919", 
                 "pl_orbper_reflink"] = ("<a refstr=BAKOS_ET_AL__2010 "
                                         "href=https://ui.adsabs.harvard.edu/abs/2010ApJ...710.1724B/abstract "
                                         "target=ref>Bakos et al. 2010</a>")
    

    # for GJ 674, tic_id 218263393, add uncertainty on pl_bmassj as 0.3 Mearth,
    # or equivalently 0.00314558 * 0.3 MJup
    # estimated from Bonfils et al. 2007 and Boisse et al. 2011 discrepancy
    # see Boisse et al. 2011 Table 2
    err = 0.00314558 * 0.3
    sps_w_ad.loc[sps_w_ad.tic_id == "218263393", "pl_bmassjerr2"] = -err
    sps_w_ad.loc[sps_w_ad.tic_id == "218263393", "pl_bmassjerr1"] = err

    # for planets with no radius error, use Chen and Kipping 2017 forecaster
    # to get the radius from M sin i
    radfrommass = lambda x: forecaster.Mstat2R(mean=x.pl_bmassj, 
                                              onesig_neg=x.pl_bmassjerr2,
                                              onesig_pos=x.pl_bmassjerr1, 
                                              unit='Jupiter',
                                              n_mass_samples=int(1e3),
                                               classify=False)
    noradius_but_masserr = (sps_w_ad.pl_radjerr1.isna()) & (~sps_w_ad.pl_bmassjerr1.isna())
    rrr = sps_w_ad[noradius_but_masserr].apply(lambda x: radfrommass(x), axis=1)

    rrr = np.array(rrr.to_list())
    rrr = pd.DataFrame({"pl_radj":rrr[:,0],
                        "pl_radjerr1":rrr[:,1],
                        "pl_radjerr2":rrr[:,2]})

    sps_w_ad.loc[noradius_but_masserr, "pl_radj"] = rrr.pl_radj.values
    sps_w_ad.loc[noradius_but_masserr, "pl_radjerr1"] = rrr.pl_radjerr1.values
    sps_w_ad.loc[noradius_but_masserr, "pl_radjerr2"] = rrr.pl_radjerr2.values

    reflink = ("<a refstr=CHEN_AND_KIPPING__2017 "
               "href=https://ui.adsabs.harvard.edu/abs/2017ApJ...834...17C "
               "target=ref>Chen and Kipping 2017</a>")
    sps_w_ad.loc[noradius_but_masserr, "pl_radj_reflink"] = reflink


    # NOW FOR THE BROWN DWARF SYSTEMS

    # LP 261-75
    lp261 = sps_w_ad.tic_id == "67646988"

    # fill in the source for the orbital period
    # Irwiin+2018 is consistent with the table, but we got this from TESS
    sps_w_ad.loc[lp261, "pl_orbper_reflink"]

    # write the table to a CSV file
    path_to_params = "../results/params_of_star_planet_systems_with_AD_tests.csv"
    sps_w_ad.to_csv(path_to_params, index=False)

    print(f"Wrote {path_to_params}")
    print(f"{len(sps_w_ad)} rows")

    # make sure all SPS are covered
    assert len(sps_w_ad) == len(unique_ids), "Not all SPS are covered!"
