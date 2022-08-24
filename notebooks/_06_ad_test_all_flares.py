"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script calculates the cumulative phase distribution of all
flares in stars with known transiting exoplanets. 
"""

import time
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




if __name__ == "main":

    # time stamp
    tstamp = time.strftime("%Y_%m_%d")

    # Define the number of step for the A2 sampling
    N = 10000

    # phaseshift zero for now
    phaseshift = 0.

    # -------------------------------------------------------------------------
    # Get the flare phases

    # Get flare table with final flares
    flare_table = pd.read_csv("../results/PAPER_flare_table.csv")

    # Sort the real flares by their phases in ascending order
    real = (flare_table.orbital_phase != -1)
    p = flare_table[real].orbital_phase.sort_values(ascending=True).values

    # Select the LCs that were searched for flares for this star
    lcs = flare_table[["timestamp", "TIC", "ID", "mission",
                        "quarter_or_sector"]].drop_duplicates()

    # Check that LCs were not searched multiple times with different timestamps
    unique_cols = ["mission", "quarter_or_sector", "TIC", "ID"]
    assert (lcs.groupby(by=unique_cols).count().timestamp.values == 1).all(), \
        print(lcs.groupby(by=unique_cols).count().timestamp.sort_values(ascending=False))

    # add phaseshift to the phases
    p = (p + phaseshift) % 1.

    # -------------------------------------------------------------------------
    # Get the observed phases
    
    # Get the total observing time in each phase bin
    location = "/media/ekaterina/USB DISK/lcs_w_phases/"
    observed_phases, binmids = get_observed_phases(p, lcs, location,
                                                phaseshift=phaseshift)

    # -------------------------------------------------------------------------
    # Get the cumulative flare phase distributions

    # Get the (cumulative) flare phase distributions
    n_exp, cum_n_exp = get_cumulative_distribution(flare_table,
                                                    observed_phases, lcs)

    # Get the null hypothesis distribution of flare phases
    # assuming flares are distributed uniformly
    f = get_null_hypothesis_distribution(p, cum_n_exp)

    # calculate the observed cumulative frequency of flares
    p = np.insert(p, 0, 0)
    p = np.append(p, 1)
    cumsum = np.cumsum(np.ones_like(p)) / len(p)

    # -------------------------------------------------------------------------
    # Make a diagnostic plot

    plt.figure(figsize=(8, 6))
    plt.plot(p, f(p), label="expected")
    plt.scatter(p, cumsum, c="r", label="observed")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('orbital phase')
    plt.ylabel('cumulative distribution of flares')
    plt.tight_layout()
    plt.savefig(f'../results/plots/{tstamp}_cum_dist_all_stars.png')

    # -------------------------------------------------------------------------
    # Save null hypothesis distribution with flare phase distribution in a file
    # for PAPER plots

    fname = f"../results/PAPER_flare_phase_distribution_all_flares_{tstamp}.csv"
    df = pd.DataFrame({"flare_phase": p, 
                    "null_hypothesis_cum_dist": cum_n_exp,
                    "observed_cum_dist": cumsum})
    df.to_csv(fname, index=False, header=True)

    # -------------------------------------------------------------------------
    # Finally, run the A-D test

    # Sample A-D statistics from the null hypothesis distribution
    A2 = sample_AD_for_custom_distribution(f, p.shape[0], N)

    # This should go into the function above
    # select only the finite values
    A2 = A2[np.isfinite(A2)]

    # Calculate the p-value and A2 value using the distribution of A2 values
    pval, atest = get_pvalue_from_AD_statistic(p, f, A2)
    print(pval, atest)

    with open("../results/multistar_adtest.csv", "a") as f:
        f.write(f"{tstamp},all,all,"
                f"{len(p)-2},"# account for the added 0 and 1
                f"{observed_phases.sum().sum()},{phaseshift},{N},{pval},{atest}\n")

    # save A2 values for all stars in a file for PAPER plots
    fname = f"../results/{tstamp}_PAPER_A2_all_stars_{phaseshift}.txt"
    np.savetxt(fname, A2, fmt="%f")