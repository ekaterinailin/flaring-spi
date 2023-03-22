"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script calculates a bunch of values about the AD tests and stellar systems,
 including:

- the magnetic field of the star based on Rossby number (Reiners+2014/2022)
- the total observing time per star
- number of covered orbits per star
- relative velocity between star and planet at the planets orbit
- SPI power (Lanza 2012 w/ and w/o planetary field, Saur+2013/Kavanagh+2022)

Final table is called: results.csv
"""


import numpy as np
import pandas as pd

from funcs.ad import aggregate_pvalues
from funcs.spirelations import (wrap_obstimes,
                                p_spi_lanza12,
                                calculate_relative_velocity,
                                rossby_reiners2014,
                                b_from_ro_reiners2022,
                                pspi_kavanagh2022)


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
    path = "../results/multistar_adtest.csv"
    print(f"[UP] Reading in AD test results from {path}")
    adtests = pd.read_csv(path)

    # make TIC a string
    adtests["TIC"] = adtests["TIC"].astype(str)

    # [READ] 
    # in the table of known parameters from the NASA Exoplanet Archive and TESS
    path = "../results/params_of_star_planet_systems_with_AD_tests.csv"
    print(f"[UP] Reading in the table of known parameters from "
          f"the NASA Exoplanet Archive {path}")
    sps_w_ad = pd.read_csv(path)

    # rename column tic_id to TIC
    sps_w_ad.rename(columns={"tic_id": "TIC"}, inplace=True)

    # make TIC a string
    sps_w_ad["TIC"] = sps_w_ad["TIC"].astype(str)

    # [READ] 
    # in age and rotation period from the literature search
    print("[UP] Reading in age and rotation period "
          "from the literature search from")
    path = "../results/2022_08_stellar_params.csv"

    print(f"[UP] Reading in age and rotation period "
          f"from the literature search from {path}")
    literature_params = pd.read_csv(path)
    # make TIC a string
    literature_params["TIC"] = literature_params["TIC"].astype(str)

    # [READ] 
    # in the flare table
    path = "../results/2022_07_flares_vetted.csv"
    print(f"[UP] Reading in the flare table from {path}")
    flares = pd.read_csv(path)

    # -------------------------------------------------------------------------
    # Next, initialize the final table by aggregating the p-values for the 
    # AD tests

    print("[CALC] Aggregating the p-values of the AD tests.")

    # aggregate the p-values
    mean_std = aggregate_pvalues(adtests, subsample="ED>1s", period="orbit")

    # Then merge in the NASA Exoplanet Archive table and the literature 
    # search table

    # merge the AD test results with the stellar parameters
    mean_std = mean_std.merge(literature_params, on=["TIC"], how="left")

    # add semi-major axis in AU, stellar radius in Rsun, 
    # planet radius in Rjup, and orbital period to the table:

    print("[MERGE] Merging in the NASA Exoplanet Archive table to the AD tests.")

    # merge the relevant part of the table
    mean_std = mean_std.merge(sps_w_ad[["TIC", "pl_orbsmax","st_rad_kepler",
                                        "st_raderr1_kepler", "st_raderr2_kepler",
                                        "st_rad_reflink", "st_lum", "st_lumerr1",
                                        "pl_orbeccen", "pl_orbeccenerr1",
                                        "pl_orbeccenerr2", "pl_orbeccen_reflink",
                                        "st_lumerr2","st_lum_reflink", "pl_orbsmax_reflink",
                                        "st_rad_tess", "pl_radj", "pl_orbper_kepler",
                                        "pl_orbpererr1_kepler", "pl_orbpererr2_kepler",
                                        "pl_orbpererr1_tess", "pl_orbpererr2_tess",
                                        "pl_orbper_tess","pl_orbper_reflink",
                                        "a_au_err", "pl_radjerr1","pl_radjerr2",
                                        "pl_radj_reflink",
                                        "sy_dist", "sy_snum"]],
                                on="TIC", how="left")

    # rename columns and fill NaNs will TESS values
    mean_std = mean_std.rename(columns = {"pl_orbsmax": "a_au",
                                          "st_rad_kepler": "st_rad",
                                          "st_raderr1_kepler": "st_rad_err1",
                                          "st_raderr2_kepler": "st_rad_err2",})

    print("[CALC] Calculating stellar radius uncertainty from the mean"
          "of the upper and lower uncertainties.")
    mean_std["st_rad_err"] =  (mean_std.st_rad_err1 + mean_std.st_rad_err2) / 2

    print("[MERGE] Fill in TESS orbital period values where Kepler values are NaN.")
    mean_std["orbper_d"] =  mean_std.pl_orbper_kepler.fillna(mean_std.pl_orbper_tess)

    print("[MERGE] Fill in TESS orbital period uncertainty values where Kepler "
          "values are NaN.")
    filledin1 = mean_std.pl_orbpererr1_kepler.fillna(mean_std.pl_orbpererr1_tess)
    filledin2 = mean_std.pl_orbpererr2_kepler.fillna(mean_std.pl_orbpererr2_tess)
    mean_std["orbper_d_err"] =  (filledin1 - filledin2) / 2

    # -------------------------------------------------------------------------
    # Fill in missing values for the rotation period uncertainty with 10% of the
    # rotation period

    print("[MERGE] Fill in missing values for the rotation period uncertainty with "
          "10 per cent of the rotation period.")
    rot_no_err = np.isnan(mean_std.st_rotp_err) & ~np.isnan(mean_std.st_rotp)
    mean_std.loc[rot_no_err, "st_rotp_err"] = mean_std[rot_no_err].st_rotp * 0.1
    
    # add [*] footnote as the st_rotp_source2
    mean_std.loc[rot_no_err, "st_rotp_source2"] = "[*]"

    # -------------------------------------------------------------------------
    # Now it's time to calculate some properties of the systems with
    # AD tests:

   
    # where still nan, calculate Ro and B from Ro with Reiners et al. (2022)
    # luminosity is in log10(L/Lsun) in Exoplanet Archive    
    print("[CALC] Calculating Rossby number from the rotation period and luminosity.")
    res = mean_std.apply(lambda x: rossby_reiners2014(10**x.st_lum, 
                                                      x.st_rotp,
                                                      error=True,
                                                      Lbol_high=10**(x.st_lum+x.st_lumerr1),
                                                      Lbol_low=10**(x.st_lum+x.st_lumerr2),
                                                      Prot_high=x.st_rotp+x.st_rotp_err, 
                                                      Prot_low=x.st_rotp-x.st_rotp_err), 
                        axis=1)
    # convert res into a 2d array
    res = np.array(res.tolist())    

    # write to columns            
    mean_std["Ro"] = res[:,0]
    mean_std["Ro_high"] = res[:,1]
    mean_std["Ro_low"] = res[:,2]


    # cond = np.isnan(mean_std.B_G) [cond]
    print("[CALC] Calculating magnetic field strength from the Rossby number.")
    res = mean_std.apply(lambda x: b_from_ro_reiners2022(x.Ro, error=True,
                                                        Ro_high=x.Ro_high,
                                                        Ro_low=x.Ro_low),
                        axis=1)

    # convert res into a 2d array
    res = np.array(res.tolist())

    # write to columns
    mean_std["B_G"] = res[:,0]
    mean_std["B_G_high"] = res[:,1]
    mean_std["B_G_low"] = res[:,2]

   
    # OBSERVING BASELINE and ORBITS COVERED

    # get the observing time for each star using the flare tables
    print("[CALC] Calculating the observing baseline for each star in days.")
    mean_std["obstime_d"] = mean_std.apply(lambda x: wrap_obstimes(str(x.TIC), 
                                                    flares), axis=1)

    # calculate the number of covered orbits
    print("[CALC] Calculating the number of orbits covered by the observing baseline.")
    mean_std["orbits_covered"] = mean_std["obstime_d"] / mean_std["orbper_d"]


    # RELATIVE VELOCITY

    # calculate relative velocity between stellar rotation at the planetary orbit
    # and the orbital velocity of the planet in km/s and
    # calculate relative velocity errors using the mean uncertainty in the
    # orbital period
    print("[CALC] Calculating the relative velocity between the stellar rotation "
            "and the orbital velocity of the planet in km/s.")
    res = mean_std.apply(lambda x: calculate_relative_velocity(x.a_au,
                                                               x.orbper_d,
                                                               x.st_rotp,
                                                               a_au_err=x.a_au_err,
                                                               orbper_err=x.orbper_d_err,
                                                               rotper_err=x.st_rotp_err,
                                                               error=True),
                         axis=1)

    # convert to 2d array
    res = np.array(res.tolist())

    # write to columns
    mean_std["v_rel_km_s"] = res[:,0]
    mean_std["v_rel_err_km_s"] = res[:,1]


    # -------------------------------------------------------------------------
    # SPI POWER FROM LANZA 2012
    # -------------------------------------------------------------------------

    # -----------
    # B_P = 1 G
    # -----------

    # calculate the SPI power from the Lanza 2012 scaling relation
    print("[CALC] Calculating the SPI power from the Lanza 2012 scaling relation.")
    res = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                 x.B_G, 
                                                 x.pl_radj,
                                                 x.a_au, 
                                                 x.st_rad,
                                                 error=True,
                                                 Blow=x.B_G_low,
                                                 Bhigh=x.B_G_high,
                                                 pl_radhigh=x.pl_radj + x.pl_radjerr1,
                                                 pl_radlow=x.pl_radj + x.pl_radjerr2,
                                                 v_rel_err=np.abs(x.v_rel_err_km_s),
                                                 Bp_err=0.,
                                                 a_err=x.a_au_err,
                                                 rstarhigh=x.st_rad + x.st_rad_err1,
                                                 rstarlow=x.st_rad + x.st_rad_err2),
                        axis=1)

    # convert res into an a 2d array
    res = np.array(res.values.tolist()).T

    # write to columns
    mean_std["p_spi_erg_s"] = res[0]
    mean_std["p_spi_erg_s_high"] = res[1]
    mean_std["p_spi_erg_s_low"] = res[2]

     # get normalization value from AU Mic p_spi_erg_s value
    print("[CALC] Calculating the normalization value for the SPI power.")
    norm = mean_std.loc[mean_std.TIC == "441420236", "p_spi_erg_s"].values[0]

    # normalize the SPI power
    print("[CALC] Normalizing the SPI power.")
    mean_std["p_spi_sb_bp1_norm"] = mean_std["p_spi_erg_s"] / norm
    mean_std["p_spi_sb_bp1_norm_high"] = mean_std["p_spi_erg_s_high"] / norm
    mean_std["p_spi_sb_bp1_norm_low"] = mean_std["p_spi_erg_s_low"] / norm

    # -------------------------------------------------------------------------

    # -----------
    # B_P = 0 G
    # -----------

    # calculate the SPI power from the Lanza 2012 scaling relation with Bp=0
    print("[CALC] Calculating the SPI power from the Lanza 2012 scaling relation "
            "with Bp=0.")
    res = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                 x.B_G,
                                                 x.pl_radj,
                                                 x.a_au, 
                                                 x.st_rad,
                                                 error=True,
                                                 Blow=x.B_G_low,
                                                 Bhigh=x.B_G_high, 
                                                 Bp=0., # THIS IS THE ONLY DIFFERENCE
                                                 pl_radhigh=x.pl_radj + x.pl_radjerr1,
                                                 pl_radlow=x.pl_radj + x.pl_radjerr2,
                                                 v_rel_err=np.abs(x.v_rel_err_km_s),
                                                 Bp_err=0.,
                                                 a_err=x.a_au_err,
                                                 rstarhigh=x.st_rad + x.st_rad_err1,
                                                 rstarlow=x.st_rad + x.st_rad_err2),
                         axis=1)
    
    
    # convert res into an a 2d array
    res = np.array(res.values.tolist()).T

    # write to columns
    mean_std["p_spi_erg_s_bp0"] = res[0]
    mean_std["p_spi_erg_s_bp0_high"] = res[1]
    mean_std["p_spi_erg_s_bp0_low"] = res[2]

    # get normalization value from AU Mic p_spi_erg_s_bp0 value
    norm = mean_std.loc[mean_std.TIC == "441420236", "p_spi_erg_s_bp0"].values[0]
    print(norm)

    # normalize the SPI power to AU Mic
    print("[CALC] Normalizing the SPI power to AU Mic.")
    mean_std["p_spi_sb_bp0_norm"] = mean_std["p_spi_erg_s_bp0"] / norm
    mean_std["p_spi_sb_bp0_norm_high"] = mean_std["p_spi_erg_s_bp0_high"] / norm
    mean_std["p_spi_sb_bp0_norm_low"] = mean_std["p_spi_erg_s_bp0_low"] / norm


    # -----------------------------------------------------------------------
    # KAVANAGH 2022 / SAUR 2013
    # -----------------------------------------------------------------------

    # -----------
    # B_P = 1 G
    # -----------

    # calculate the SPI power from the Saur et al. 2013 / Kavanagh et al. (2022)
    # scaling relation with Bp=1
    print("[CALC] Calculating the SPI power from the Saur et al. 2013 / Kavanagh et al. (2022) "
            "scaling relation with Bp=1.")
    res = mean_std.apply(lambda x: pspi_kavanagh2022(x.pl_radj,
                                                     x.B_G,
                                                     np.abs(x.v_rel_km_s), 
                                                     x.a_au,
                                                     error=True, 
                                                     Rphigh=x.pl_radj + x.pl_radjerr1, 
                                                     Rplow=x.pl_radj + x.pl_radjerr2,
                                                     Bhigh=x.B_G_high,
                                                     Blow=x.B_G_low,
                                                     vrelhigh=np.abs(x.v_rel_km_s) + x.v_rel_err_km_s,
                                                     vrellow=np.abs(x.v_rel_km_s) - x.v_rel_err_km_s,
                                                     ahigh=x.a_au + x.a_au_err, 
                                                     alow=x.a_au - x.a_au_err,
                                                     Bphigh=1.,
                                                     Bplow=1.),
                        axis=1)

    # convert res into an a 2d array
    res = np.array(res.values.tolist()).T

    # write to columns
    mean_std["p_spi_kav22"] = res[0]
    mean_std["p_spi_kav22_high"] = res[1]
    mean_std["p_spi_kav22_low"] = res[2]
    
    # get normalization value from AU Mic p_spi_kav22 value
    norm = mean_std.loc[mean_std.TIC == "441420236", "p_spi_kav22"].values[0]
    
    # normalize the SPI power to AU Mic
    print("[CALC] Normalizing the SPI power to AU Mic.")
    mean_std["p_spi_aw_bp1_norm"] = mean_std["p_spi_kav22"] / norm
    mean_std["p_spi_aw_bp1_norm_high"] = mean_std["p_spi_kav22_high"] / norm
    mean_std["p_spi_aw_bp1_norm_low"] = mean_std["p_spi_kav22_low"] / norm

    # -----------
    # B_P = 0 G
    # -----------

    # calculate the SPI power from the Saur et al. 2013 / Kavanagh et al. (2022)
    # scaling relation with Bp=0

    print("[CALC] Calculating the SPI power from the Saur et al. 2013 / Kavanagh et al. (2022) "
            "scaling relation with Bp=0.")

    res = mean_std.apply(lambda x: pspi_kavanagh2022(x.pl_radj,
                                                        x.B_G,
                                                        np.abs(x.v_rel_km_s),
                                                        x.a_au,
                                                        error=True,
                                                        Bp=0.,
                                                        Rphigh=x.pl_radj + x.pl_radjerr1,
                                                        Rplow=x.pl_radj + x.pl_radjerr2,
                                                        Bhigh=x.B_G_high,
                                                        Blow=x.B_G_low,
                                                        vrelhigh=np.abs(x.v_rel_km_s) + x.v_rel_err_km_s,
                                                        vrellow=np.abs(x.v_rel_km_s) - x.v_rel_err_km_s,
                                                        ahigh=x.a_au + x.a_au_err,
                                                        alow=x.a_au - x.a_au_err,
                                                        Bphigh=0.,
                                                        Bplow=0.),
                            axis=1)

    # convert res into an a 2d array
    res = np.array(res.values.tolist()).T

    # write to columns
    mean_std["p_spi_kav22_bp0"] = res[0]
    mean_std["p_spi_kav22_bp0_high"] = res[1]
    mean_std["p_spi_kav22_bp0_low"] = res[2]

    # get normalization value from AU Mic p_spi_kav22 value
    norm = mean_std.loc[mean_std.TIC == "441420236", "p_spi_kav22_bp0"].values[0]

    # normalize the SPI power to AU Mic
    print("[CALC] Normalizing the SPI power to AU Mic.")
    mean_std["p_spi_aw_bp0_norm"] = mean_std["p_spi_kav22_bp0"] / norm
    mean_std["p_spi_aw_bp0_norm_high"] = mean_std["p_spi_kav22_bp0_high"] / norm
    mean_std["p_spi_aw_bp0_norm_low"] = mean_std["p_spi_kav22_bp0_low"] / norm

    # -------------------------------------------------------------------------
    # Reverse the 10% error on the rotation period
    print("[REVERSE] rotation period uncertainty with "
          "10 per cent of the rotation period.")
    
    mean_std.loc[mean_std.st_rotp_source2 == "[*]", "st_rotp_err"] = np.nan
    
    del mean_std["st_rotp_source2"]


    # -------------------------------------------------------------------------
    # BIBKEYS
    # -------------------------------------------------------------------------

    # For transparency, ADD BIBKEYS to the table for the literature values


    # read in the reflink to bibkey mapping table
    print("[UP] Reading in the reflink to bibkey mapping table.")
    bibkeys = pd.read_csv("../results/reflink_to_bibkey.csv")

    # 1. orbital period
    print("[MERGE] Adding the bibkeys for orbital period to the results table.")
    mean_std["pl_orbper_bibkey"] = mean_std.apply(lambda x: 
                                                  map_bibkey(x.pl_orbper_reflink, 
                                                             bibkeys),
                                                  axis=1)

    # 2. stellar radius
    print("[MERGE] Adding the bibkeys for stellar radius to the results table.")
    mean_std["st_rad_bibkey"] = mean_std.apply(lambda x: 
                                                map_bibkey(x.st_rad_reflink,
                                                           bibkeys),
                                                axis=1)

    # 3. orbital semi-major axis
    print("[MERGE] Adding the bibkeys for orbital semi-major axis to the results table.")
    mean_std["pl_orbsmax_bibkey"] = mean_std.apply(lambda x:
                                                    map_bibkey(x.pl_orbsmax_reflink,
                                                                bibkeys),
                                                axis=1)

    # 4. stellar luminosity
    print("[MERGE] Adding the bibkeys for stellar luminosity to the results table.")
    mean_std["st_lum_bibkey"] = mean_std.apply(lambda x:
                                                map_bibkey(x.st_lum_reflink,
                                                            bibkeys),
                                                axis=1) 

    # 5. planet radius
    print("[MERGE] Adding the bibkeys for planet radius to the results table.")
    mean_std["pl_radj_bibkey"] = mean_std.apply(lambda x:
                                                map_bibkey(x.pl_radj_reflink,
                                                            bibkeys),
                                                axis=1) 

    # 6. orbital eccentricity
    print("[MERGE] Adding the bibkeys for orbital eccentricity to the results table.")       
    mean_std["pl_orbeccen_bibkey"] = mean_std.apply(lambda x:
                                                      map_bibkey(x.pl_orbeccen_reflink,
                                                                  bibkeys),
                                                      axis=1)

                                                      



    # delete the reflink columns
    mean_std.drop(columns=["pl_orbper_reflink"], inplace=True)
    mean_std.drop(columns=["pl_orbsmax_reflink"], inplace=True)
    mean_std.drop(columns=["st_rad_reflink"], inplace=True)
    mean_std.drop(columns=["st_lum_reflink"], inplace=True)
    mean_std.drop(columns=["pl_radj_reflink"], inplace=True)
    mean_std.drop(columns=["pl_orbeccen_reflink"], inplace=True)



    # -------------------------------------------------------------------------
    # SAVE THE RESULTS TABLE
    # -------------------------------------------------------------------------

    # [WRITE] 
    # to file
    print("[DOWN] writing final results table to")
    path = "../results/results.csv"
    print(path)
    mean_std.to_csv(path, index=False)

    # [WRITE]
    # to paper folder
    print("[DOWN] writing final results table to")
    path = "../../../002_writing/flaring-spi-paper/src/data/results.csv"
    print(path)
    mean_std.to_csv(path, index=False)

                    