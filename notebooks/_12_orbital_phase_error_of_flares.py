"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2023, MIT License

This script calculates the phase error of the orbital period of the innermost 
planet. The error is calculated from the uncertainty on the transit time or the
time of the first flare (for RV planets), and the uncertainty in the orbital period.

The formula is:

phase error = abs((transit time or first flare time) - flare time) / (orbital period^2) * (orbital period error)
"""

import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt


if __name__ == "__main__":

    # read in orbital period, error, first absolute flare start time
    params = pd.read_csv("../results/results.csv")

    s2 = params[["ID","abs_tstart_min"]]
    # s2.TIC = s2.TIC.astype(str)

    # read in transit midtime and error
    spisys = pd.read_csv("../data/2022_07_27_input_catalog_star_planet_systems.csv")

    # take TESS unless it's NaN, then take Kepler
    spisys["tranmid"] = spisys["pl_tranmid"]
    # take the average of the two errors and fill in NaNs if necessary
    spisys["tranmid_err"] = (spisys["pl_tranmiderr1"] - spisys["pl_tranmiderr2"]) / 2
    spisys["pl_orbper_err"] = (spisys["pl_orbpererr1"] - spisys["pl_orbpererr2"]) / 2

    # rename tic_ic to TIC  
    spisys.hostname = spisys.hostname.fillna(spisys.TIC)
    spisys = spisys.rename(columns={"hostname":"ID"})


    # select only the columns we need
    s1 = spisys[["ID","tranmid","tranmid_err",'pl_orbper','pl_orbper_err']]
    # s1.TIC = s1.TIC.astype(str)

    # merge the two dataframes
    s = pd.merge(s1,s2,on="ID",how="left", validate="1:1")

    # if tranmid is NaN, use the absolute flare start time and its error with 2 min uncertainty
    s["tranmid"] = s["tranmid"].fillna(s["abs_tstart_min"])
    s["tranmid_err"] = s["tranmid_err"].fillna(2 / 60 / 24) # 2 min uncertainty for the flare start time

    # -------------------------------------------------------------------------

    # Read in the flare table
    df = pd.read_csv("../results/PAPER_flare_table.csv")
    if "orbital_phase_err" in df.columns:   
        df = df.drop(columns=["orbital_phase_err"])
   
    # merge s on flare table
    df2 = pd.merge(df,s,on="ID",how="left", validate="m:1")

    # calculate the orbital phase error
    df2["orbital_phase_err"] = np.abs(df2["tranmid"] - df2["abs_tstart"]) / (df2["pl_orbper"]**2) * df2["pl_orbper_err"]

    # # plot the orbital phase error vs. absolute flare start time for each star
    # for tic, g in df.groupby("TIC"):
    #     g = g.sort_values("abs_tstart",ascending=True)
    #     plt.plot(g["abs_tstart"],g["orbital_phase_err"],linewidth=2)

    # plt.yscale("log")
    # plt.xscale("log")

    print(df2[df2.TIC == "299096355"])

    # drop the helper columns
    df2 = df2.drop(columns=["tranmid","tranmid_err","abs_tstart_min","pl_orbper","pl_orbper_err"])


    assert df.shape[0] == df2.shape[0]
    assert df.shape[1] == df2.shape[1] -1


    # save the dataframe
    df2.to_csv("../results/PAPER_flare_table.csv",index=False)

    # save to flaring spi paper directory
    path = "/home/ekaterina/Documents/002_writing/flaring-spi-paper/src/data/"
    df2.to_csv(path + "PAPER_flare_table.csv", index=False)