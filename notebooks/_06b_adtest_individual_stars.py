"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script runs the A-D test on individual stars.
Requires the paper flare table, and the de-trended light curves with phases.

Output goes to multistar_adtest.csv
"""


import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


from funcs.phaseanalysismulti import (get_observed_phases,
                                      get_cumulative_distribution,
                                      get_null_hypothesis_distribution,
                                     ) 

from funcs.ad import (sample_AD_for_custom_distribution,
                      get_pvalue_from_AD_statistic,
                      )

import time

if __name__ == '__main__':
        

    # time stamp
    tstamp = time.strftime("%Y_%m_%d")

    # Define the number of step for the A2 sampling
    N = 10000

    # phaseshift zero for now
    phaseshift = 0.75
    print('phaseshift = ', phaseshift)

    # Get flare table with final flares
    flare_table = pd.read_csv("../results/PAPER_flare_table.csv")

    # Select only flare with ED :
    flare_table = flare_table[flare_table['ED'] > 0.]
    # flare_table = flare_table[flare_table.TIC == 267749737]
    # pick a star to test
    for TIC, flare_table_single_star in flare_table.groupby("TIC"):

    # get entries for the TIC
    # flare_table_single_star = flare_table[flare_table.TIC == TIC]
        if (flare_table_single_star.shape[0] >= 1):

            ID = flare_table_single_star.ID.iloc[0]
            if ID == "AU Mic":
                print("We exclude AU Mic for now.")
                continue
            
            print(ID)
            # Sort the real flares by their phases in ascending order
            real = (flare_table_single_star.orbital_phase != -1) 
            p = flare_table_single_star[real].orbital_phase
            p = p.sort_values(ascending=True).values

            # Select the LCs that were searched for flares for this star
            lcs = flare_table_single_star[["timestamp","TIC","ID","mission",
                                        "quarter_or_sector"]].drop_duplicates()

            # Check that LCs were not searched multiple times with 
            # different timestamps
            unique_cols = ["mission","quarter_or_sector","TIC","ID"]

            # make sure that the LCs are unique, and none are searched
            # multiple times
            lccounts = lcs.groupby(by=unique_cols).count().timestamp.values
            assert (lccounts == 1).all(), \
                print(lcs)

            # Get the total observing time in each phase bin
            location = "/media/ekaterina/USB DISK/lcs_w_phases/"
            observed_phases, binmids  = get_observed_phases(p, lcs, location)

            # Get the (cumulative) flare phase distributions j
            n_exp, cum_n_exp = get_cumulative_distribution(flare_table_single_star,
                                                        observed_phases, lcs)
            # Get the null hypothesis distribution of flare phases
            # assuming flares are distributed uniformly
            f = get_null_hypothesis_distribution(p, cum_n_exp)
            print(len(p))
            print(p)
            p = np.sort(np.random.rand(len(p)))
            print(p)

            if len(p) > 2:

                # Make a diagnostic plot
                plt.figure(figsize=(8,6))
                p = np.insert(p,0,0)
                p = np.append(p,1)
                plt.plot(p,f(p))
                cumsum =  np.cumsum(np.ones_like(p)) / len(p)
                plt.scatter(p, cumsum, c="r")
                plt.title(f"TIC {TIC}, {ID}")
                plt.xlim(0,1)
                plt.ylim(0,1)
                plt.savefig(f"../results/plots/{tstamp}_TIC_{TIC}_cumhist.png")
                plt.close()

                # Finally, run the A-D test
                A2 = sample_AD_for_custom_distribution(f, p.shape[0], N)

                # This should go into the function above
                # select only the finite values
                A2 = A2[np.isfinite(A2)]

                # Calculate the p-value and A2 value using the distribution 
                # of A2 values
                pval, atest = get_pvalue_from_AD_statistic(p, f, A2)
                print(pval, atest)

                with open("../results/multistar_adtest.csv", "a") as f:
                    f.write(f"{tstamp},{TIC},{ID},"
                            f"{len(p)-2},"# account for the added 0 and 1
                            f"{observed_phases.sum().sum()},{phaseshift},"
                            f"{N},{pval},{atest},all,randomized_orbit\n")
