"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script uses runs AD tests on individual RV stars using orbital periods
from the NASA Exoplanet Archive* (column description**) 

* https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PSCompPars
** https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
"""


import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


from funcs.phaseanalysismulti import (get_rotational_phases,
                                      get_observed_phases,
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
    # phaseshift = 0.

    # Get flare table with final flares
    flare_table = pd.read_csv("../results/PAPER_flare_table.csv")

    # Get table with orbital periods from NASA Exoplanet Archive
    orbital_table = pd.read_csv("../data/2022_11_15_input_catalog_NONtransit_star_planet_systems.csv")

    # Select only flare with ED > 1.:
    flare_table = flare_table[flare_table['ED'] > 0.]

    for phaseshift in [0., 0.25, 0.5, 0.75]:
        print('---------------\nphaseshift = ', phaseshift)


        # pick a star to test
        for TIC, flare_table_single_star in flare_table.groupby("TIC"):

        # get entries for the TIC
        # flare_table_single_star = flare_table[flare_table.TIC == TIC]
            if (flare_table_single_star.shape[0] >= 1):

                ID = flare_table_single_star.ID.iloc[0]
                
                # Get rotation period
                # print(rotation_table.TIC, TIC)
                try:
                    orbital_period = orbital_table[orbital_table.TIC == int(TIC)].pl_orbper.iloc[0]
                except IndexError:
                    print('No rotation period for TIC', TIC)
                    continue
                print(ID, orbital_period)

                # real flares are those where an ED is given
                real = (~flare_table_single_star.ED.isnull()) 

                # -----------------------------------------------------------------
                # Get the pseudo-orbitak phase instead of orbital
                real_flares = flare_table_single_star[real]
                p = get_rotational_phases(real_flares.tstart, orbital_period,
                                        real_flares.mission)

                # sort values and apply phaseshift
                p = np.sort((p + phaseshift) % 1)

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
                observed_rotational_phases, binmids = get_observed_phases(p, lcs, location, phaseshift=phaseshift,
                                                    period="rotation", rotper=orbital_period)

                print(observed_rotational_phases)
                print(p)

                # -----------------------------------------------------------------
                # Get the (cumulative) flare phase distributions j
                # CHECK HERE IF THE PHASES ARE CORRECTLY CALCULATED
                n_exp, cum_n_exp = get_cumulative_distribution(flare_table_single_star,
                                                            observed_rotational_phases, lcs)

                # Get the null hypothesis distribution of flare phases
                # assuming flares are distributed uniformly
                f = get_null_hypothesis_distribution(p, cum_n_exp)
            
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
                                f"{observed_rotational_phases.sum().sum()},{phaseshift},"
                                f"{N},{pval},{atest},all,orbit\n")
