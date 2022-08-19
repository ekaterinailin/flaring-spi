import pandas as pd 
import numpy as np

def get_cumulative_distribution(flare_table_single_star, observedphases, lcs):
    """Calculate the cumulative distribution function of observed
    flares, while accounting for the number of observed flares in
    a light curve with a given detection threshold.

    Parameters:
    -----------
    flare_table_single_star : pd.DataFrame
        detected flares table, quarter/sector and mission columns 
        that denote the LC
    observedphases : DataFrame
        table that lists the observed time in each phase bin, 
        with a column for each LC
    lcs : pd.DataFrame
        table that lists the target's LCs that were searched for flares

    Returns:
    --------
    n_exp : np.array
        df of expected frequency of flares in each bin
    cum_n_exp : np.array
        cdf of expected frequency of flares in each bin
    """
    # Measured number of flares in each bin
    F = observedphases.shape[0]

    # Expected number of flares in each bin

    # setup array for each bin
    n_exp = np.zeros(F).astype(float)

    # sum over different sectors to account for different detection thresholds
    for i, row in lcs.iterrows():

        # Get the observed phases for this LC, the values indicate the
        # observed time at a given phase???!!!
        lcid = f"{row.mission}_{row.quarter_or_sector}"

        # Get the observed phases for this LC, or if it does not exist,
        # throw an error
        try:
            obstimes_k = observedphases[lcid]
        except KeyError:
            raise KeyError(f"{lcid} column not found in observedphases table.")

        print(lcid)

        # Check that flare phases are correct:
        assert (obstimes_k <= 1).all(), \
               f"One or more observedphases[{lcid}] > 1!"
        assert (obstimes_k >= 0).all(), \
               f"One or more observedphases[{lcid}] < 0!"

        # Get the total number of obstimes in this LC
        tot_obstimes_k = obstimes_k.sum()
        pick_lc = ((flare_table_single_star.quarter_or_sector == row.quarter_or_sector) &
                   (flare_table_single_star.mission == row.mission))
        F_k = flare_table_single_star[pick_lc].shape[0]
        print(F_k, tot_obstimes_k)
        n_exp_k = (obstimes_k * F_k / tot_obstimes_k).values
        n_exp += n_exp_k
        
    # for debugging: 
    # calculate the old cumulative distribution that 
    # ignores the different detection thresholds
    # cum_n_exp_alt = dfphases.sum(axis=1).values.cumsum() / dfphases.sum(axis=1).values.sum()
    
    # Calculate cumulative frequency distribution
    _ =  n_exp.cumsum()
    
    # CDF maximum should be 1 
    cum_n_exp = _ / np.nanmax(_)

    # CDF minimum should be 0
    cum_n_exp = np.insert(cum_n_exp, 0, 0)
    
    return n_exp, cum_n_exp#, cum_n_exp_alt return alternative only if necessary
