import pandas as pd
import numpy as np

from funcs.phaseanalysismulti import OFFSET

import matplotlib.pyplot as plt

def check_output(tminmax, TIC, missions, per_and_errs):
    """Check the output of the coherence timescale calculation.

    Parameters
    ----------
    tminmax : list
        List of minimum and maximum times to consider.
    TIC : int
        TIC ID of the target.
    separate : bool
        Whether to calculate the coherence timescale for each mission
        separately.
    missions : list
        List of missions to consider.
    per_and_errs : list
        List of periods and their errors.
    
    Returns
    -------
    None
    """

    system_string = fr"System {TIC}:\n" \
                     fr"  tminmax = {tminmax}\n" \
                     fr"  missions = {missions}\n" \
                     fr"  periods = {per_and_errs}\n" 
    
    # print a message if the orbital period is not finite
    if ~np.isfinite(tminmax["orbper"]):
        print(f"Orbital period not finite for TIC {row.tic_id}\n")
        print(system_string)
    
    # print a message if the orbital period error is not finite
    if ~np.isfinite(tminmax["orbper_err"]):
        print(f"Orbital period error not finite for TIC {row.tic_id}\n")
        print(system_string)


def coherence_timescale(period, error):
    """Calculate the coherence timescale from the period and the error on the
    period.

    Explanation for dummies like myself:

    If the error is the same as the period, after 1 period the undertainty
    the orbital phase is the same as the period. If the error is 0.5 times
    the period, after 2 periods the uncertainty is the same as the period.
    In these cases, the coherence timescale is 1 or 2 times the period.

    Parameters
    ----------
    period : float
        The period in days.
    error : float
        The error on the period in days.

    Returns
    -------
    float
        The coherence timescale in days.

    """
    if error == 0:
        return np.inf
    elif ~np.isfinite(error):
        return np.inf

    return period**2 / error

if __name__ == "__main__":
    
    # ------------------------------------------------------------------------------
    # 1. ROTATION PERIOD COHERENCE

    # Read in stellar parameters table
    path = "../results/2022_08_stellar_params.csv"
    df = pd.read_csv(path)
    print(f"Get stellar parameters from\n{path}\n")
    if "time_span_d" in df.columns:
        for l in ["abs_tstart_min",	"abs_tstart_max", "time_span_d", 
                  "coherence_timescale_rotation_d", "coherence_ratio_rotation",
                  "tstart_min",	"tstart_max", "orbper", "orbper_err",
                  "timespan_d",	"coherence_timescale_orbit_d", "coherence_ratio_orbit"]:

            del df[l]

    # Read in flare table
    path = "../results/PAPER_flare_table.csv"
    flares = pd.read_csv(path)
    print(f"Get flare table from\n{path}\n")

    # shift each tstart by mission specific offset using mission column and OFFSET
    # column
    flares["abs_tstart"] = flares.tstart + flares["mission"].apply(lambda x: OFFSET[x])

    # Get first and last absolute flare times for each star
    tmin = flares.groupby("TIC").abs_tstart.min()
    tmax = flares.groupby("TIC").abs_tstart.max()

    # merge tmin and tmax into a single DataFrame
    tminmax = pd.merge(tmin, tmax, on="TIC", suffixes=("_min", "_max"))

    # Calculate the time span between first and last flare
    tminmax["time_span_d"] = tminmax["abs_tstart_max"] - tminmax["abs_tstart_min"]

    # merge tminmax with the table of stellar parameters
    df_timespan = pd.merge(df, tminmax, on="TIC")

    g = lambda x: coherence_timescale(x.st_rotp, x.st_rotp_err)
    df_timespan["coherence_timescale_rotation_d"] = df_timespan.apply(g, axis=1)

    # Calculate the coherence ratio 
    # as the ratio between the total time between first and last flare
    # and the time it takes to lose the information on rotational phase
    # Ideally, the number is << 1
    df_timespan["coherence_ratio_rotation"] = (df_timespan["time_span_d"] / 
                                    df_timespan["coherence_timescale_rotation_d"])

    # ------------------------------------------------------------------------------
    # plot the result and save the file
    plt.figure(figsize=(7,5.5))
    plt.hist(df_timespan["coherence_ratio_rotation"], bins=np.logspace(-3,1.,30),
            facecolor="orange", edgecolor="maroon")
    plt.xscale("log")

    plt.xlabel("time span covered by flares / coherence time of rotation period",
                size=12)
    plt.ylabel("number of star-planet systems")
    plt.tight_layout()
    plt.savefig("../results/plots/coherence_time_of_rotation_period.png", dpi=300)
    # --------------------------------------------------------------------------

    # write to file
    df_timespan.to_csv("../results/2022_08_stellar_params.csv", 
                    index=False)

    # --------------------------------------------------------------------------

    # 2. ORBITAL PERIOD COHERENCE TIMES
    # This is a little more complicated because we sometimes have Kepler orbital
    # periods only but TESS observations without TESS orbital periods

    columns = [ "tstart_min",	"tstart_max", "TIC","orbper", "orbper_err",
                  "timespan_d",	"coherence_timescale_orbit_d", "coherence_ratio_orbit"]

    # Read in stellar parameters table
    path = "../results/2022_08_stellar_params.csv"
    df = pd.read_csv(path)
    print(f"Get stellar parameters from\n{path}\n")

    # Get orbital period of SPSs
    path = "../results/params_of_star_planet_systems_with_AD_tests.csv"
    orbits = pd.read_csv(path)
    # cast tic_id to str
    orbits["tic_id"] = orbits["tic_id"].astype(str)
    print(f"Get get orbital periods and uncertainties from\n{path}\n")
    
    # Read in flare table
    path = "../results/PAPER_flare_table.csv"
    flares = pd.read_csv(path)
    print(f"Get flare table from\n{path}\n")

    # shift each tstart by mission specific offset using mission column and OFFSET
    # column
    flares["abs_tstart"] = flares.tstart + flares["mission"].apply(lambda x: OFFSET[x])

    # Init a coherence timescales table
    ct = pd.DataFrame(columns=columns)

    # go through the SPSs table
    for i, row in orbits.iterrows():

        # take the real flares for this TIC
        f = flares[(flares["TIC"].astype(str)==str(row.tic_id)) & (np.isfinite(flares["ED"]))]
    
        # If both transit midtimes for Kepler and TESS are present
        if (np.isfinite(row.pl_orbper_kepler)) & (np.isfinite(row.pl_orbper_tess)):
            separate = True
        # Otherwise use the one available
        elif (np.isfinite(row.pl_orbper_kepler)) | (np.isfinite(row.pl_orbper_tess)):
            separate = False
        else:
            # If no transit midtimes are available, skip this SPS
            continue

        # in which missions were flares observed?
        missions = f.mission.unique()

        # in how many missions were flares observed?
        nmissions = missions.size
    
        # If only one mission, pick the min and max flare times without 
        # grouping by mission
        if nmissions == 1:

            # no need to get absolute flare times, just get min and max
            tmin = f.tstart.min()
            tmax = f.tstart.max()
            

            # merge the two values into one Series
            tminmax = pd.Series([tmin, tmax, row.tic_id], 
                                index=["tstart_min", "tstart_max", "TIC"])

            # select the mission and get the orbital period and the error on the
            # orbital period

            # define coloumns with the orbital period and uncertainty to use
            cols = ["pl_orbper_kepler","pl_orbpererr1_kepler", "pl_orbpererr2_kepler"]

            # if TESS ...
            if missions[0] == "TESS":
                for i, col in enumerate(cols):
                    # if the column is not empty, use it
                    if np.isfinite(row[col]):
                        continue
                    # otherwise, use the TESS column
                    else:
                        cols[i] = cols[i].replace("kepler", "tess")
                

            # add them to the Series
            tminmax["orbper"] = row[cols[0]]
            tminmax["orbper_err"] = (row[cols[1]] - row[cols[2]]) / 2.


        # if there are multiple missions, pick the first and last flare times while
        # grouping by mission if separate is True
        elif nmissions == 2:
            # This case does not occur!
            if separate:
                raise ValueError("This case does not occur! And if it would,"
                                "I would have to rewrite the code!")
            # if separate is False, group by TIC only    
            # because you have both missions present
            else:

                tmin = f.abs_tstart.min()
                tmax = f.abs_tstart.max()
                tminmax = pd.Series([tmin, tmax, row.tic_id],
                                    index=["tstart_min", "tstart_max", "TIC"])


                # select the mission where tranmid time is from and get the orbital 
                # period and the error on the orbital period
                if np.isfinite(row.pl_orbper_kepler):
                    # define coloumns with the orbital period and uncertainty to use
                    cols = ["pl_orbper_kepler","pl_orbpererr1_kepler", "pl_orbpererr2_kepler"]
                    
                elif np.isfinite(row.pl_orbper_tess):
                    # define coloumns with the orbital period and uncertainty to use
                    cols = ["pl_orbper_tess","pl_orbpererr1_tess", "pl_orbpererr2_tess"]
                else:
                    raise ValueError("No orbital period available")
            tminmax["orbper"] = row[cols[0]]
            tminmax["orbper_err"] = (row[cols[1]] - row[cols[2]]) / 2. 
            

        elif nmissions==0:

            continue

        # Check if the output is ok
        # tminmax, TIC, missions, per_and_errs

        check_output(tminmax, row.tic_id, missions, row[cols])

        # Calculate the time span covered by the flares
        tminmax["timespan_d"] = tminmax.tstart_max - tminmax.tstart_min


        # Calculate the coherence time of the orbital period
        tminmax["coherence_timescale_orbit_d"] = coherence_timescale(tminmax.orbper, 
                                                                     tminmax.orbper_err)        
        # Calculate the coherence ratio for the orbital period
        tminmax["coherence_ratio_orbit"] = (tminmax.timespan_d / 
                                                    tminmax.coherence_timescale_orbit_d)

        ct = pd.concat([ct, pd.DataFrame(tminmax).T], ignore_index=True)

    ct.TIC = ct.TIC.astype(str)
    df.TIC = df.TIC.astype(str)
    df_timespan = pd.merge(df, ct, on="TIC")

    # ------------------------------------------------------------------------------
    # check if the table has all SPSs in it that have ad tests
    
    assert df_timespan.shape[0] == 40


    print("All 40 SPSs with AD tests are in the table")

    # ------------------------------------------------------------------------------
    # plot the result and save the file
    plt.figure(figsize=(7,5.5))
    plt.hist(df_timespan["coherence_ratio_orbit"], 
            bins=np.logspace(-6,-2.,30),
            facecolor="orange", edgecolor="maroon")
    plt.xscale("log")

    plt.xlabel("time span covered by flares / coherence time of orbital period",
                size=12)
    plt.ylabel("number of star-planet systems")
    plt.tight_layout()
    plt.savefig("../results/plots/coherence_time_of_orbital_period.png", dpi=300)


    # write to file
    df_timespan.to_csv("../results/2022_08_stellar_params.csv", 
                    index=False)
    # ------------------------------------------------------------------------------