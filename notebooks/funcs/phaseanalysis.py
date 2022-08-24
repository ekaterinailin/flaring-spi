import numpy as np
import pandas as pd


import matplotlib.pyplot as plt

def get_cumulative_distributions(df, dfphases, get_secs_cadences):
    """Calculate the cumulative distribution function of observed
    flares, while accounting for the number of observed flares in
    a light curve with a given detection threshold.

    Parameters:
    -----------
    df : DataFrame
        detected flares table, with a flag qcs that denotes the LC
    dfphases : DataFrame
        table that lists the observed phases, with a column for each LC
    get_secs_cadences : list
        format: [(Sector, cadence in minutes), (.., ..), (.., ..)]
        Sector is an int
        cadence is a float
    """
    # Measured number of flares in each bin
    F = dfphases.shape[0]
    n_i = np.full(F, 1)

    # Calculate cumulative frequency distribution
    cum_n_i = np.cumsum(n_i) / F

    # Expected number of flares in each bin

    # setup array for each bin
    n_exp = np.zeros_like(n_i).astype(float)

    # sum over different sectors to account for different detection thresholds
    for sec, cad in get_secs_cadences:
        obstimes_k = dfphases[sec]
        tot_obstimes_k = obstimes_k.sum()
        F_k = df[df.qcs==sec].shape[0]
        n_exp_k = (obstimes_k * F_k / tot_obstimes_k).values
        n_exp += n_exp_k
        
    # for debugging: 
    # calculate the old cumulative distribution that 
    # ignores the different detection thresholds
    # cum_n_exp_alt = dfphases.sum(axis=1).values.cumsum() / dfphases.sum(axis=1).values.sum()
    
    # Calculate cumulative frequency distribution
    cum_n_exp = n_exp.cumsum() / F
    
    return n_i, n_exp, cum_n_exp, cum_n_i#, cum_n_exp_alt return alternative only if necessary

def get_flare_phases(df, mode, rotper=None, starname="AU Mic", lcn=0, mission="TESS",
                     path="../results/observedtimes/"):
    phases = []
    
    if mode=="Orbit":
        for j, row in df.iterrows():
            try:
                lc = pd.read_csv(f"{path}{starname}_{row.qcs}_{lcn}_{mission}.csv")
                phases.append(lc.phase[np.argmin(np.abs(lc.time - row.tstart))])
            except Exception as e:
                print(e)
                phases.append(np.nan)
        df["phases"] = phases
                
    elif mode=="Rotation":
        
        df["phases"] = (df.tstart % rotper) / rotper
        
    return df

def get_observed_phases(mode, p, get_secs_cadences, rotper=None, phaseshift=0., 
                        starname="AU Mic", lcn=0, mission="TESS", test="AD"):
    """Takes 
    
    
    Parameters:
    -----------
    mode : str
        either "Rotation" or "Orbit"
    p : array
        array of phases, for the KS-test it should be the 
        measured flare phases, for the AD test, any well 
        sampled grid of phases works
    get_secs_cadences : list
        format: [(Sector, cadence in minutes), (.., ..), (.., ..)]
        Sector is an int
        cadence is a float
    rotper : float
        Rotation period
    phaseshift : float
        default = 0, the starting phase is arbitrary, 
        so adding a phaseshift should not change the result for
        the AD-test, but for the KS-test
    starname : str
        name of the object, needed to identify dataset
    lcn : int
        number of light curve, needed to identify dataset
    mission : str
        mission name, needed to identify dataset
    test : str
        "KS" or "AD"
        
    Return:
    -------
    pandas.DataFrame - one column per Sector or Quarter,
                       and observing time per phase in each row
    list - mids of bins for each phase row in the DataFrame
    """
    # bin array is two elements longer than the number of bins
    # to include 0 and 1
    bins = np.zeros(len(p) +2) 
    
    # add zero and one
    bins[0] = 0 
    bins[-1] = 1 
    
    # and the others are kept as defined by observations
    bins[1:-1] = p 
#     print(bins)
    aumicphases = pd.DataFrame()
#     plt.figure(figsize=(8, 3))
    for qcs, cadence in get_secs_cadences:
        lc = pd.read_csv(f"../results/observedtimes/{starname}_{qcs}_{lcn}_{mission}.csv")
        
        if mode=="Orbit":
            # add phaseshift as modulo 1 to wrap phases > 1 back
            counts, bins = np.histogram((lc.phase.values + phaseshift) % 1, bins=bins)
            
        elif mode=="Rotation":
            # add phaseshift as modulo 1 to wrap phases > 1 back
            counts, bins = np.histogram(((lc.time % rotper) / rotper + phaseshift) % 1, bins=bins)
            
        if qcs==27:
            linestyle="solid"
        else:
            linestyle = "dashed"
#         plt.plot((bins[1:] + bins[:-1])/2, counts * cadence, c="k", linestyle=linestyle, label=qcs)
        
        if test == "KS": 
            # circular boundary condition
            counts[0] = counts[0] + counts[-1]

            # remove last bin to avoid double counting
            counts = counts[:-1]

        # get observing times for each Sector
        aumicphases[qcs] = counts * cadence

        #     plt.legend()
#     plt.savefig(f"../results/plots/{mode}_histogram_{phaseshift}.png", dpi=300)
    return aumicphases, (bins[1:] + bins[:-1])/2
