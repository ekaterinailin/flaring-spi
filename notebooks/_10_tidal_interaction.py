"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script uses the procedures from Ilic et al. (2022)  to calculate the
tidal interaction between the star and the planet in the systems, and combines
the results with the results of the AD tests."""

import numpy as np
import pandas as pd
from funcs.ad import aggregate_pvalues
from funcs.masses_and_radii import calculate_abs_Ks, mann_mass_from_abs_Ks

import subprocess
import forecaster

from astroquery.gaia import Gaia

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

    # ------------------------------------------------------------------------
    # read in stellar parameters
    params = pd.read_csv("../results/params_of_star_planet_systems_with_AD_tests.csv")

    # rename columns
    p = params[['st_mass','st_masserr1','st_masserr2',
                'pl_bmassj','pl_bmassjerr1','pl_bmassjerr2','tic_id',
                'st_mass_reflink', 'pl_bmassj_reflink',
                'pl_radj','pl_radjerr1','pl_radjerr2',
                'sy_kmag','sy_kmagerr1', "hostname"]]


    # make TIC a string
    p.tic_id = p.tic_id.astype(str)

    # read in the reflink to bibkey mapping table
    print("[UP] Reading in the reflink to bibkey mapping table.")
    bibkeys = pd.read_csv("../results/reflink_to_bibkey.csv")

    # 1. stellar mass
    print("[MERGE] Adding the bibkeys for stellar mass to the results table.")
    p["st_mass_bibkey"] = p.apply(lambda x: map_bibkey(x.st_mass_reflink, 
                                                             bibkeys),
                                                  axis=1)
    
    # 2. planet mass
    print("[MERGE] Adding the bibkeys for planetary mass to the results table.")
    p["pl_bmassj_bibkey"] = p.apply(lambda x: map_bibkey(x.pl_bmassj_reflink,
                                                                bibkeys),
                                                                axis=1)
    
    # ------------------------------------------------------------------------


    # calculate the planetary masses in the missing entries
    massfromrad = lambda x: forecaster.Rstat2M(mean=x.pl_radj, 
                                               onesig_neg=x.pl_radjerr2,
                                               onesig_pos=x.pl_radjerr1, 
                                               unit='Jupiter',
                                               n_radii_samples=int(1e3),
                                               classify=False,
                                               )

    nomasserr = ((p.pl_bmassjerr1.isna()) | p.pl_bmassj_bibkey.isna()) & (~p.pl_radjerr2.isna())

    rrr = p[nomasserr].iloc[::-1].apply(lambda x: massfromrad(x), axis=1)

    rrr = np.array(rrr.to_list())

    rrr = pd.DataFrame({"M_pl": rrr[:,0],
                        "M_pl_err1": rrr[:,1],
                        "M_pl_err2": rrr[:,2]})

    p.loc[nomasserr, "pl_bmassj"] = rrr.M_pl.values
    p.loc[nomasserr, "pl_bmassjerr1"] = rrr.M_pl_err1.values
    p.loc[nomasserr, "pl_bmassjerr2"] = -rrr.M_pl_err2.values
    p.loc[nomasserr, "pl_bmassj_bibkey"] = "chen2017probabilistic"

    # -----------------------------------------------------------------------    
    # add the stellar masses in the missing entries

    # add missing distance from Gaia DR3 Bailer-Jones to GJ 3323, GJ 674 and GJ 3082
    p.loc[p.hostname == "GJ 3323", "dist_pc"] = 5.373613
    p.loc[p.hostname == "GJ 674", "dist_pc"] = 4.552019
    p.loc[p.hostname == "GJ 3082", "dist_pc"] = 16.632895

    # add errors to the distances
    p.loc[p.hostname == "GJ 3323", "dist_pc_err"] =  (5.374399 - 5.372851) / 2.
    p.loc[p.hostname == "GJ 674", "dist_pc_err"] = (4.55265 - 4.5514603) / 2.
    p.loc[p.hostname == "GJ 3082", "dist_pc_err"] = (16.63665 - 16.62884) / 2.

    ids = ["GJ 3323", "GJ 674", "GJ 3082"]

    # calculate distance modulus and err
    p.loc[p.hostname.isin(ids), "dist_mod"] = 5 * np.log10(p.loc[p.hostname.isin(ids), "dist_pc"]) - 5
    p.loc[p.hostname.isin(ids), "dist_mod_err"] = (5 * p.loc[p.hostname.isin(ids), "dist_pc_err"] / 
                                             (p.loc[p.hostname.isin(ids), "dist_pc"] * np.log(10)))
    p.loc[p.hostname.isin(ids), "st_mass_bibkey"] = "mann2015how"

    absk = p[p.hostname.isin(ids)].apply(lambda x: calculate_abs_Ks(x.dist_mod,
                                                                          x.dist_mod_err,
                                                                            x.sy_kmag,
                                                                            x.sy_kmagerr1),
                                                                            axis=1)
    abskmag, abskmagerr = np.array(absk.to_list()).T   


    p.loc[p.hostname.isin(ids), "st_mass"], p.loc[p.hostname.isin(ids), "st_masserr1"] = mann_mass_from_abs_Ks(abskmag,
                                                                                                 abskmagerr)
    
    # set second error to first error
    p.loc[p.hostname.isin(ids), "st_masserr2"] = - p.loc[p.hostname.isin(ids), "st_masserr1"].values


    ocs = ['st_mass', 'st_masserr1', 'st_masserr2',
           'pl_bmassj', 'pl_bmassjerr1', 'pl_bmassjerr2',
           'tic_id', "pl_bmassj_bibkey", "st_mass_bibkey"]
    newcs = ['M_star', 'M_star_up_err','M_star_low_err',
             'M_pl','M_pl_up_err','M_pl_low_err',
             'TIC',"pl_bmassj_bibkey", "st_mass_bibkey"]

    p = p.rename(index=str, columns=dict(zip(ocs, newcs)))

    p = p[newcs]

    print(p)


    # ------------------------------------------------------------------------
    # read in the other results with more stellar and planetary parameters
    nr = pd.read_csv("../results/results.csv")

    #select columns
    nr = nr[["multiple_star","TIC", "ID",  "st_rotp","st_rotp_err", "orbper_d", 
             "orbper_d_err", "a_au", "a_au_err", "st_rad","st_rad_err1", "st_rad_err2"] ]

    # rename columns
    cols =  ["multiple_star", "TIC", "ID",   "P_rot", "P_rot_err", "P_orb", 
             "P_orb_err", "a", "a_err", "R_star", "R_star_up_err","R_star_low_err"]
    
    nr.columns = cols

    # placeholder just to make sure we are convective envelope stars
    nr["Teff"] = 4000. 

    nr.TIC = nr.TIC.astype(str)
    new = pd.merge(p, nr, on="TIC", how="inner")

    new['P_rot_err'] = new['P_rot_err'].fillna(0.1*new['P_rot'])


    # add Kepler-42 c upper limit of 2.06 Earth mass from muirhead2012characterizing
    new.loc[new.ID == "Kepler-42", "M_pl_up_err"] = 0
    new.loc[new.ID == "Kepler-42", "M_pl_low_err"] = 2.06 /  317.82838
    new.loc[new.ID == "Kepler-42", "M_pl"] = 2.06 / 317.82838

    # add bibkey for Kepler-42 c
    new.loc[new.ID == "Kepler-42", "pl_bmassj_bibkey"] = "muirhead2012characterizing"


    new = new.dropna(subset=cols[1:], how="any")
   
    new.to_csv("tidal/params.csv", index=False)

    # call Nikoleta's code to calculate expected power of interaction
    subprocess.run(["python", "tidal/tidal_interaction_strength.py"])


    # read in the results
    res = pd.read_csv("tidal/TIS_all.csv")
    res.TIC = res.TIC.astype(str)

    # read in the p-values
    pvals = pd.read_csv("../results/multistar_adtest.csv")

    # aggregate the p-values
    mean_std = aggregate_pvalues(pvals, subsample="ED>1s", period="half_orbit")

    # make TIC a string
    mean_std.TIC = mean_std.TIC.astype(str)

    # merge the results
    df = pd.merge(res, mean_std, on="TIC", how="inner")

    # remove GJ 1061 because it does not have a rotation period
    df = df[df["ID"] != "GJ 1061"]

    # drop rows with TIC 67646988 and 236387002, the brown dwarfs
    df = df[df.TIC != "67646988" ]
    df = df[df.TIC != "236387002" ]

    # the old Kepler-411 instance
    df = df[df.TIC != "399954349(c)" ]

    # remove multiple stars
    res = df[df["multiple_star"].isnull()]

    # write to file
    res.to_csv("tidal/TIS_with_ADtests.csv", index=False)
    res.to_csv("/home/ekaterina/Documents/002_writing/flaring-spi-paper/src/data/TIS_with_ADtests.csv", index=False)
