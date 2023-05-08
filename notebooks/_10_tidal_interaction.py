"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script uses the procedures from Ilic et al. (2022)  to calculate the
tidal interaction between the star and the planet in the systems, and combines
the results with the results of the AD tests."""

import pandas as pd
from funcs.ad import aggregate_pvalues

import subprocess

if __name__ == "__main__":

    # read in stellar parameters
    params = pd.read_csv("../results/params_of_star_planet_systems_with_AD_tests.csv")

    # rename columns
    p = params[['st_mass','st_masserr1',
                'pl_bmassj','pl_bmassjerr1','tic_id',]]
    p.columns = ['M_star','M_star_err','M_pl','M_pl_err','TIC']

    # add in a very generous estimate of the masses of the star and planet
    p["M_pl_err"] = p["M_pl_err"].fillna(.5*p["M_pl"])
    p["M_star_err"] = p["M_star_err"].fillna(.2*p["M_star"])

    # make TIC a string
    p.TIC = p.TIC.astype(str)


    # read in the other results with more stellar and planetary parameters
    nr = pd.read_csv("../results/results.csv")

    #select columns
    nr = nr[["multiple_star","TIC", "ID",  "st_rotp","st_rotp_err", "orbper_d", "orbper_d_err", "a_au", "a_au_err", "st_rad","st_rad_err1"] ]

    # rename columns
    cols =  ["multiple_star", "TIC", "ID",   "P_rot", "P_rot_err", "P_orb", "P_orb_err", "a", "a_err", "R_star", "R_star_err"]
    nr.columns = cols

    # placeholder just to make sure we are convective envelope stars
    nr["Teff"] = 4000. 
    nr.TIC = nr.TIC.astype(str)
    new = pd.merge(p, nr, on="TIC", how="inner")
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
