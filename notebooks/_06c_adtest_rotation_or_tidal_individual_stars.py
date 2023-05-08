"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script performs AD tests with rotational variability, or whatever
period is chosen. Here, we also use P_orb/2 to calculate the deviation
for potential tidal interaction."""



import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


import warnings
warnings.filterwarnings("ignore")

from funcs.phaseanalysismulti import (get_rotational_phases,
                                      get_observed_phases,
                                      get_cumulative_distribution,
                                      get_null_hypothesis_distribution,
                                     ) 

from funcs.ad import (sample_AD_for_custom_distribution,
                      get_pvalue_from_AD_statistic,
                      )

def shuffle_column(df, column_name):
    # Store the column values in a variable
    column_values = df[column_name].values
    
    # Use numpy's random.shuffle function to shuffle the values
    np.random.shuffle(column_values)
    
    # Update the column in the DataFrame with the shuffled values
    df[column_name] = column_values
    
    return df

import time

if __name__ == '__main__':
        

    # time stamp
    tstamp = time.strftime("%Y_%m_%d")

    # Define the number of step for the A2 sampling
    N = 10000

    # Get flare table with final flares
    flare_table = pd.read_csv("../results/PAPER_flare_table.csv")

    # Get table with rotation periods
    rotation_table = pd.read_csv("../results/params_of_star_planet_systems_with_AD_tests.csv")

    # Fill in missing pl_orbper_kepler values with pl_orbper_tess
    rotation_table.pl_orbper_kepler.fillna(rotation_table.pl_orbper_tess, inplace=True)
   
    #pd.read_csv("../results/2022_08_stellar_params.csv")

    # Select only flare with ED > 1.:
    flare_table = flare_table[flare_table['ED'] >= 1.]

    # remove EPIC 200164267
    flare_table = flare_table[flare_table.ID != "EPIC 200164267"]

    readres = pd.read_csv("../results/multistar_adtest.csv")

    # get all IDs that have been tested with ED>1s and phasshift 0.25
    tested_IDs = readres[(readres.subsample=="ED>1s") &
                         (readres.phaseshift == 0.25) &
                         (readres.period == "half_orbit")].ID.values

    print(flare_table.ID.unique().shape)

    # remove all stars that have been tested
    flare_table = flare_table[~flare_table.ID.isin(tested_IDs)]

    # remove  TIC 399954349 (c)
    flare_table = flare_table[flare_table.TIC != "399954349(c)"]

    print(flare_table.ID.unique().shape)
    
    for phaseshift in [0, 0.25, 0.5, .75]:
        
        print('phaseshift = ', phaseshift)

        # pick a star to test
        for TIC, flare_table_single_star in flare_table.groupby("TIC"):

            # get entries for the TIC
            # flare_table_single_star = flare_table[flare_table.TIC == TIC]
            if (flare_table_single_star.shape[0] >= 1):

                ID = flare_table_single_star.ID.iloc[0]
                
                # Get rotation period
                # print(rotation_table.TIC, TIC)
                try:
                    _ = rotation_table[rotation_table.tic_id == int(TIC)]
                    # rotation_period = _.st_rotp.iloc[0]
                    # use half the orbital period instead
                    rotation_period = _.pl_orbper_kepler.iloc[0] / 2.
                except IndexError:
                    print('No rotation period for TIC', TIC)
                    continue
                print(ID, rotation_period)

                # Sort the real flares by their phases in ascending order
                # real = (flare_table_single_star.orbital_phase != -1) 

                # -----------------------------------------------------------------
                # Get the rotational phase instead of orbital
                # replace p = flare_table_single_star[real].orbital_phase
                real_flares = flare_table_single_star#[real]
                p = get_rotational_phases(real_flares.tstart, rotation_period,
                                        real_flares.mission)

                p = np.sort(p)
                # print(p)

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
                observed_rotational_phases, binmids = get_observed_phases(p, lcs,
                                                                        location, period="rotation",
                                                                        rotper=rotation_period,
                                                                       phaseshift=phaseshift)
             
                # Get the (cumulative) flare phase distributions j
                # CHECK HERE IF THE PHASES ARE CORRECTLY CALCULATED
                n_exp, cum_n_exp = get_cumulative_distribution(flare_table_single_star,
                                                            observed_rotational_phases, lcs)
                
              
                # Get the null hypothesis distribution of flare phases
                # assuming flares are distributed uniformly
                f = get_null_hypothesis_distribution(p, cum_n_exp)
                print(len(p))
                # print(p)

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
                                f"{N},{pval},{atest},ED>1s,half_orbit\n")
