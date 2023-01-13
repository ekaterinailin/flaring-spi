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
                                        "st_raderr1_kepler", "st_raderr2_kepler",
                                        "st_rad_reflink", "st_lum", "st_lumerr1",
                                        "st_lumerr2","st_lum_reflink",
                                        "st_rad_tess", "pl_radj", "pl_orbper_kepler",
                                        "pl_orbpererr1_kepler", "pl_orbpererr2_kepler",
                                        "pl_orbpererr1_tess", "pl_orbpererr2_tess",
                                        "pl_orbper_tess","pl_orbper_reflink",
                                        "a_au_err", "pl_radjerr1","pl_radjerr2",
                                        "sy_dist"]],
                                on="TIC", how="left")

    # rename columns and fill NaNs will TESS values
    mean_std = mean_std.rename(columns = {"pl_orbsmax": "a_au",
                                          "st_rad_kepler": "st_rad",
                                          "st_raderr1_kepler": "st_rad_err1",
                                          "st_raderr2_kepler": "st_rad_err2",})
    mean_std["st_rad_err"] =  (mean_std.st_rad_err1 + mean_std.st_rad_err2) / 2
    mean_std["orbper_d"] =  mean_std.pl_orbper_kepler.fillna(mean_std.pl_orbper_tess)
    mean_std["orbper_d_err"] =  (mean_std.pl_orbpererr1_kepler.fillna(mean_std.pl_orbpererr1_tess) +
                                 mean_std.pl_orbpererr2_kepler.fillna(mean_std.pl_orbpererr2_tess)) / 2

    # -------------------------------------------------------------------------
    # Fill in missing values for the rotation period uncertainty with 10% of the
    # rotation period

    rot_but_no_err = np.isnan(mean_std.st_rotp_err) & ~np.isnan(mean_std.st_rotp)
    mean_std.loc[rot_but_no_err, "st_rotp_err"] = mean_std[rot_but_no_err].st_rotp * 0.1
    
    # add [*] footnote as the st_rotp_source
    mean_std.loc[rot_but_no_err, "st_rotp_source"] = "[*]"

    # -------------------------------------------------------------------------
    # Now it's time to calculate some properties of the systems with
    # AD tests:

   
    # X-RAY LUMINOSITY

    # assume uncertainty on Lx is 50% of the value based on Wright et al 2011/2018
    # for the values from Foster and Poppenhäger 2022
    mean_std["xray_flux_err_erg_s"] = mean_std.xray_flux_erg_s * 0.3

   
    # MAGNETIC FIELD

    # calculate B field from X-ray luminosity and stellar radius
    # mean_std["B_G"] = mean_std.apply(lambda x: b_from_lx_reiners(x.xray_flux_erg_s, x.st_rad), axis=1)

    # # get uncertainty in B field from X-ray luminosity and stellar radius uncertainty
    # # and use the intrinsic scatter, too
    # mean_std["high_B_G"] = mean_std.apply(lambda x: b_from_lx_reiners(x.xray_flux_erg_s,
    #                                                 x.st_rad, error=True, 
    #                                                 r_err=x.st_rad_err,
    #                                                 lx_err=x.xray_flux_err_erg_s)[1], axis=1)

    # mean_std["low_B_G"] = mean_std.apply(lambda x: b_from_lx_reiners(x.xray_flux_erg_s,
    #                                                 x.st_rad, error=True,
    #                                                 r_err=x.st_rad_err,
    #                                                 lx_err=x.xray_flux_err_erg_s)[2], axis=1)
   
    # where still nan, calculate Ro and B from Ro with Reiners et al. (2022)
    # luminosity is in log10(L/Lsun) in Exoplanet Archive
    mean_std["Ro"] = mean_std.apply(lambda x: rossby_reiners2014(10**x.st_lum, x.st_rotp), axis=1)
    mean_std["Ro_high"] = mean_std.apply(lambda x: rossby_reiners2014(10**x.st_lum, x.st_rotp, error=True,
                                                                    Lbol_high=10**(x.st_lum+x.st_lumerr1),
                                                                    Lbol_low=10**(x.st_lum+x.st_lumerr2),
                                                                    Prot_high=x.st_rotp+x.st_rotp_err, 
                                                                    Prot_low=x.st_rotp-x.st_rotp_err)[1], axis=1)        
    mean_std["Ro_low"] = mean_std.apply(lambda x: rossby_reiners2014(10**x.st_lum, x.st_rotp, error=True,
                                                                    Lbol_high=10**(x.st_lum+x.st_lumerr1),
                                                                    Lbol_low=10**(x.st_lum+x.st_lumerr2),
                                                                    Prot_high=x.st_rotp+x.st_rotp_err, 
                                                                    Prot_low=x.st_rotp-x.st_rotp_err)[2], axis=1)        


    # cond = np.isnan(mean_std.B_G) [cond]
    mean_std["B_G"] = mean_std.apply(lambda x: b_from_ro_reiners2022(x.Ro), axis=1)
    mean_std["high_B_G"] = mean_std.apply(lambda x: b_from_ro_reiners2022(x.Ro, error=True, 
                                                        Ro_high=x.Ro_high, Ro_low=x.Ro_low)[1], axis=1)
    mean_std["low_B_G"] = mean_std.apply(lambda x: b_from_ro_reiners2022(x.Ro, error=True,
                                                        Ro_high=x.Ro_high, Ro_low=x.Ro_low)[2], axis=1)

   
    # OBSERVING BASELINE and ORBITS COVERED

    # get the observing time for each star using the flare tables
    mean_std["obstime_d"] = mean_std.apply(lambda x: wrap_obstimes(str(x.TIC), 
                                                    flares), axis=1)

    # calculate the number of covered orbits
    mean_std["orbits_covered"] = mean_std["obstime_d"] / mean_std["orbper_d"]


    # RELATIVE VELOCITY

    # calculate relative velocity between stellar rotation at the planetary orbit
    # and the orbital velocity of the planet in km/s
    mean_std["v_rel_km_s"] = mean_std.apply(lambda x: 
                                            calculate_relative_velocity(x.a_au,
                                            x.orbper_d, x.st_rotp), axis=1)
    # cqalculate relative velocity errors using the mean uncertainty in the
    # orbital period
    mean_std["v_rel_err_km_s"] = mean_std.apply(lambda x: 
                                            calculate_relative_velocity(x.a_au,
                                            x.orbper_d, x.st_rotp, a_au_err=x.a_au_err,
                                            orbper_err=x.orbper_d_err,
                                            rotper_err=x.st_rotp_err, error=True)[1], axis=1)


    # SPI POWER 

    # calculate the SPI power from the Lanza 2012 scaling relation
    mean_std["p_spi_erg_s"] = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                            x.B_G, x.pl_radj, x.a_au), axis=1)   

    # calculate the uncertainty in the SPI power from the Lanza 2012 scaling relation
    mean_std["p_spi_erg_s_high"] = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                            x.B_G, x.pl_radj,x.a_au, error=True,
                                                            Blow=x.low_B_G, Bhigh=x.high_B_G,
                                                            pl_radhigh=x.pl_radj + x.pl_radjerr1,
                                                            pl_radlow=x.pl_radj + x.pl_radjerr2,
                                                            v_rel_err=np.abs(x.v_rel_err_km_s),
                                                            Bp_err=0.,a_err=x.a_au_err)[1], axis=1)
  
    mean_std["p_spi_erg_s_low"] = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                            x.B_G, x.pl_radj, x.a_au, error=True,
                                                            Blow=x.low_B_G, Bhigh=x.high_B_G,
                                                            pl_radhigh=x.pl_radj + x.pl_radjerr1,
                                                            pl_radlow=x.pl_radj + x.pl_radjerr2,
                                                            v_rel_err=np.abs(x.v_rel_err_km_s),
                                                            Bp_err=0., a_err=x.a_au_err)[2], axis=1)



    # calculate the SPI power from the Lanza 2012 scaling relation with Bp=0
    mean_std["p_spi_erg_s_bp0"] = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                            x.B_G, x.pl_radj,x.a_au, Bp=0.), axis=1) 

    # calculate the uncertainty in the SPI power from the Lanza 2012 scaling relation
    mean_std["p_spi_erg_s_bp0_high"] = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                            x.B_G, x.pl_radj,x.a_au, error=True,
                                                            Blow=x.low_B_G, Bhigh=x.high_B_G,
                                                            pl_radhigh=x.pl_radj + x.pl_radjerr1,
                                                            pl_radlow=x.pl_radj + x.pl_radjerr2,
                                                            v_rel_err=np.abs(x.v_rel_err_km_s), Bp=0.,
                                                            Bp_err=0.,a_err=x.a_au_err)[1], axis=1)
  
    mean_std["p_spi_erg_s_bp0_low"] = mean_std.apply(lambda x: p_spi_lanza12(np.abs(x.v_rel_km_s),
                                                            x.B_G, x.pl_radj, x.a_au, error=True,
                                                            Blow=x.low_B_G, Bhigh=x.high_B_G,
                                                            pl_radhigh=x.pl_radj + x.pl_radjerr1,
                                                            pl_radlow=x.pl_radj + x.pl_radjerr2,
                                                            v_rel_err=np.abs(x.v_rel_err_km_s), Bp=0.,
                                                            Bp_err=0., a_err=x.a_au_err)[2], axis=1)

    # calculate the SPI power from the Saur et al. 2013 / Kavanagh et al. (2022)
    # scaling relation with Bp=0
    res = mean_std.apply(lambda x: pspi_kavanagh2022(x.pl_radj, x.B_G, np.abs(x.v_rel_km_s), 
                                                     x.a_au, error=True, 
                                                     Rphigh=x.pl_radj + x.pl_radjerr1, 
                                                     Rplow=x.pl_radj + x.pl_radjerr2,
                                                     Bhigh=x.high_B_G, Blow=x.low_B_G,
                                                     vrelhigh=np.abs(x.v_rel_km_s) + x.v_rel_err_km_s,
                                                     vrellow=np.abs(x.v_rel_km_s) - x.v_rel_err_km_s,
                                                     ahigh=x.a_au + x.a_au_err, 
                                                     alow=x.a_au - x.a_au_err,
                                                     Bphigh=1., Bplow=1.), axis=1)
    # convert res into an a 2d array
    res = np.array(res.values.tolist()).T

    # write to columns
    mean_std["p_spi_kav22"] = res[0]
    mean_std["p_spi_kav22_high"] = res[1]
    mean_std["p_spi_kav22_low"] = res[2]
    
    # -------------------------------------------------------------------------
    # For transparency, add bibkeys to the table for the literature values

    # read in the reflink to bibkey mapping table
    bibkeys = pd.read_csv("../results/reflink_to_bibkey.csv")

    # 1. orbital period
    mean_std["pl_orbper_bibkey"] = mean_std.apply(lambda x: map_bibkey(x.pl_orbper_reflink, bibkeys),
                                                    axis=1)

    # 2. stellar radius
    mean_std["st_rad_bibkey"] = mean_std.apply(lambda x: map_bibkey(x.st_rad_reflink, bibkeys),
                                                    axis=1)


    # delete the reflink columns
    mean_std.drop(columns=["pl_orbper_reflink"], inplace=True)
    mean_std.drop(columns=["st_rad_reflink"], inplace=True)





    # -------------------------------------------------------------------------
    # Finally, save the table

    # [WRITE] 
    # to file
    mean_std.to_csv("../results/results.csv", index=False)

    # [WRITE]
    # to paper folder
    mean_std.to_csv("../../../002_writing/flaring-spi-paper/src/data/results.csv",
                    index=False)