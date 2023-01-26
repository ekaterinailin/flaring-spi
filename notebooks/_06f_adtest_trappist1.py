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
    for  phaseshift in [0., 0.5,0.25,0.75]:
        print('phaseshift = ', phaseshift)

        # Get flare table with final flares
        flare_table = pd.read_csv("../results/PAPER_flare_table.csv")

        # Select only flare with ED :
        flare_table = flare_table[flare_table['ED'] > 0.]

        #select TRAPPIST-1
        TIC = 278892590
        flare_table_single_star = flare_table[flare_table.TIC == TIC]
        
        ID = flare_table_single_star.ID.iloc[0]
                
        print(ID)
        # Sort the real flares by their phases in ascending order
        real = (flare_table_single_star.orbital_phase != -1) 
        p = flare_table_single_star[real].orbital_phase
        p = p.sort_values(ascending=True).values

        # read in the observed phases for trappist 1 with numpy
        k2_phases = (np.loadtxt(f"../results/obsphases_trappist1_b.txt") + phaseshift) % 1.

        # bin the observed phases
        bins = np.insert(p, 0, 0)
        bins = np.append(bins, 1.)
        observedphases, bins = np.histogram(k2_phases, bins=bins)

        # total observed time in days
        tot_obstimes_k = k2_phases.shape[0] / 60. / 24. 

        # Get the number of flares in this LC
        F_k = len(p)

        # calculate the expected frequency of flare per time bin
        n_exp = (observedphases * F_k / tot_obstimes_k)

        # Calculate cumulative frequency distribution
        _ =  n_exp.cumsum()

        # CDF maximum should be 1 
        cum_n_exp = _ / np.nanmax(_)

        # CDF minimum should be 0
        cum_n_exp = np.insert(cum_n_exp, 0, 0)

        # Get the null hypothesis distribution of flare phases
        # assuming flares are distributed uniformly
        print(len(p), len(cum_n_exp))
        f = get_null_hypothesis_distribution(p, cum_n_exp)
        print(len(p))

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
                        f"{tot_obstimes_k},{phaseshift},"
                        f"{N},{pval},{atest},all,orbit\n")
