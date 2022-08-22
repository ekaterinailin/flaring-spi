import glob


import pandas as pd 
import numpy as np

from astropy.io import fits


import glob
def read_lightcurves_from_series(series, location):
    """
    Read lightcurves using the information in the pandas.Series.

    Parameters
    ----------
    series : pandas.Series
        Dataframe containing the lightcurve information.
    location : str
        Location of the lightcurves ending with a slash.

    Returns
    -------
    list of lightcurves
    """

    # read in the LC: get ID to look for in path
    if series.mission == "Kepler":
        ID = series.ID
    elif series.mission == "TESS":
        ID = series.TIC

    # read in the LCs
    try:
        # find all paths that contain the following string
        string = f"{series.timestamp}_{ID}_{series.quarter_or_sector:02}_altai_"
      
        paths = glob.glob(f"{location}{string}*")
       
        lcs = [fits.open(path)[1].data for path in paths]
        if len(lcs) > 0:
            return lcs
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        string = f"{ID}_{series.quarter_or_sector:02}_altai_0"
       
        return [fits.open(f"{location}{series.timestamp}_{string}.fits")[1].data]

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
            allcols = observedphases.columns.values
            cols = allcols[np.isin(allcols, lcid)]
            obstimes_k = observedphases[cols]

            if len(cols) == 0:
                raise ValueError(f"{lcid} column not found in observedphases table.")
        except KeyError:
            print("no LC")
            raise KeyError(f"{lcid} column not found in observedphases table.")

        # Check that flare phases are correct:
        assert (obstimes_k.values <= 1).all(), \
               f"One or more observedphases[{lcid}] > 1!"
        assert (obstimes_k.values >= 0).all(), \
               f"One or more observedphases[{lcid}] < 0!"

        # Get the total number of obstimes in this LC
        tot_obstimes_k = obstimes_k.sum().sum()
        pick_lc = ((flare_table_single_star.quarter_or_sector == row.quarter_or_sector) &
                   (flare_table_single_star.mission == row.mission))
        F_k = flare_table_single_star[pick_lc].shape[0]
      
        n_exp_k = (obstimes_k * F_k / tot_obstimes_k).sum(axis=1)
        n_exp += n_exp_k
        
    # for debugging: 
    # calculate the old cumulative distribution that 
    # ignores the different detection thresholds
    # cum_n_exp_alt = dfphases.sum(axis=1).values.cumsum() / dfphases.sum(axis=1).values.sum()
    
    # Calculate cumulative frequency distribution
    _ =  n_exp.cumsum().values
    
    # CDF maximum should be 1 
    cum_n_exp = _ / np.nanmax(_)

    # CDF minimum should be 0
    print("cumnexp", cum_n_exp)
    cum_n_exp = np.insert(cum_n_exp, 0, 0)
    
    return n_exp, cum_n_exp#, cum_n_exp_alt return alternative only if necessary



def get_observed_phases(p, lcs, location):
    """Calculate the observed phases for a given set of flares.
    
    Parameters:
    -----------
    p : array
        array of phases for the flares
    lcs : pd.DataFrame
        table that lists the target's LCs that were searched for flares
    location : str
        path to the folder with the LCs

    Returns:
    --------
    observedphases : pd.DataFrame
        table that lists the observed phases, with a column for each LC
    binmids : array
        array of bin midpoints
    """
    # bin array is two elements longer than the number of bins
    # to include 0 and 1
    bins = np.zeros(len(p) + 2) 

    # add zero and one
    bins[0] = 0 
    bins[-1] = 1 

    # and the others are kept as defined by observations
    bins[1:-1] = p 

    # define the bin mids
    binmids = (bins[1:] + bins[:-1]) / 2

    # initialize the observed phases table
    observedphases = pd.DataFrame()

    # loop over the LCs
    for i, row in lcs.iterrows():
   
        # read in the lightcurve
        lightcurves = read_lightcurves_from_series(row, location)
       
        for i, lightcurve in enumerate(lightcurves):
         
            # get the phases and bin them
            counts, bins = np.histogram(lightcurve["phase"], bins=bins)

            # circular boundary condition
            counts[0] = counts[0] + counts[-1]

            # remove last bin to avoid double counting
            counts = counts[:-1]

            # get observing times for each lightcurve using cadence
            cadence = np.nanmin(np.diff(lightcurve["time"]))
      
            observedphases[f"{row.mission}_{row.quarter_or_sector}_{i}"] = counts * cadence

    return observedphases, binmids

