"""
Python 3.8 -- UTF-8

Ekaterina Ilin MIT License (2023)

This script downloads the K2 light curve of TRAPPIST-1, to get the observed 
phases of TRAPPIST-1 b, and the phases at which the flares occured. The latter 
is strored in the flare table, the former in a text file. Both results are used 
in script _06f where the test for star-planet interaction is performed.
"""

import lightkurve as lk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":

    # download the light curve
    lcs = lk.search_lightcurve("EPIC 200164267", cadence="short", mission="K2")
    k2lcshort = lcs[0].download()   

    # from Gillon et al. 2017, get the tranit time in BJD_TDB
    bjdtt = 7322.51736 + 2450000    

    # from Gillon et al. 2017, get the orbital period of TRAPPIST-1 b
    orbper = 1.51087081

    # take the light curve and convert time to BJD
    bjdtime = k2lcshort.time.value + 2454833. # adding the K2 offset

    # only use finite flux values to calculate observed phases
    finite_flux = ~np.isnan(np.array(k2lcshort.flux.value))

    # calculate the observed phases of the light curve, subtracting the transit
    # midtime and dividing by the orbital period mod 1, so that zero phase is the
    # transit midtime, and the phases are between 0 and 1
    obsphases = (bjdtime[finite_flux] - bjdtt) / orbper % 1.

    # save the observed phases
    print("[DOWN] Saving observed phases of TRAPPIST-1 to")
    path = "../results/obsphases_trappist1_b.txt"
    print(path)
    np.savetxt(path, obsphases)

    # read in flare table
    print("[UP] Reading in flare table from")
    path = "../results/2022_07_flares_vetted.csv"
    print(path)
    flares = pd.read_csv(path)

    # selection criterion for TRAPPIST-1 flares
    trappist = (flares.ID=="EPIC 200164267")

    # calculate the observed phases of the flares
    # by adding the K2 offset to the flare times, subtracting the transit midtime
    # and dividing by the orbital period mod 1, so that zero phase is the transit
    # midtime, and the phases are between 0 and 1
    flare_phases = (flares.loc[trappist, "tstart"]  + 2454833. - bjdtt) / orbper % 1.
    
    # add flare phases to flare table
    flares.loc[trappist, "phase"] = flare_phases

    # save the flare table
    print("[DOWN] Saving flare table with TRAPPIST-1 flares to")
    path = "../results/2022_07_flares_vetted.csv"
    print(path)
    flares.to_csv(path, index=False)

    # -------------------------------------------------------------------------
    # MAKE A VALIDATION PLOT

    # bin the light curve data
    bins = np.linspace(0,1,3000)

    # make a dataframe with finite flux values
    df = pd.DataFrame({'phase': obsphases, 
                    'flux': k2lcshort.flux.value[finite_flux]})

    # group data by phase bins and use the mean
    df = df.groupby(pd.cut(df.phase, bins=bins)).mean()

    # plot the light curve
    plt.figure(figsize=(10,5))
    plt.plot((df.phase + .5)%1., df.flux/np.mean(df.flux), 'k.')
    plt.plot([.5 - 36.4/60/24, .5 + 36.4/60/24], [1.001, 1.001], 'r-')
    plt.ylim(.98,1.02)
    plt.title("Phase-folded K2 light curve of TRAPPIST-1 b")
    plt.tight_layout()

    # save the plot
    print("[DOWN] Saving phase-folded light curve of TRAPPIST-1 b to")
    path = "../results/plots/phasefolded_trappist1_b.png"
    print(path)
    plt.savefig(path, dpi=300)




