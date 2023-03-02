import os

from lightkurve import search_lightcurve

import numpy as np
import pandas as pd

# offsets for BTJD and BKJD
offset = {"K2":2454833.,
          "Kepler":2454833.,
          "TESS":2457000.,
          'Transiting Exoplanet Survey Satellite (TESS)':2457000.}    


def get_flare_phases(flare_cadenceno, observed_phases, observed_cadenceno):
    """Get flare phases for flare candidates.

    Parameters
    ----------
    flare_cadenceno : array_like
        Cadence numbers of flare candidates.
    observed_phases : array_like
        Observed phases in light curve.
    observed_cadenceno : array_like
        Cadence numbers of observed phases in light curve.

    Returns
    -------
    flare_phases : array_like
        Flare phases.
    """
    # convert to integers for quick and precise search
    observed_cadenceno = observed_cadenceno.astype(int)
    flare_cadenceno = flare_cadenceno.astype(int)

    # get indices into observed phases array
    indices = np.where(np.isin(observed_cadenceno, flare_cadenceno, 
                               assume_unique=True))[0]
    # get flare phases
    flare_phases = observed_phases[indices]
    
    # output should be same length as flare_cadenceno
    message = ("Input and output number of flares do not match."
               "Any NaN? Any out of range cadence numbers?")
    assert len(flare_phases) == len(flare_cadenceno), message

    # check if flare phases are finite
    if np.isfinite(flare_phases).all():
        return flare_phases
    else:
        raise ValueError("Flare phases are not finite.")
    
   

def get_observed_phases(time, transit_midtime, orbital_period):
    """Get observed phases of a light curve.
    
    Parameters:
    -----------
    time: numpy.array
        time of light curve
    transit_midtime: float
        transit midtime of the planet
    orbital_period: float
        orbital period of the planet
    
    Return:
    -------
    phases: numpy.array
        observed phases of the light curve
    """
    # check if inputs are all finite
    assert np.isfinite(time).all(), "Time is not finite."
    assert np.isfinite(transit_midtime), "Transit midtime is not finite."
    assert np.isfinite(orbital_period), "Orbital period is not finite."
    assert orbital_period > 0, "Orbital period is not positive."

    # return the observed phases
    phases = ((time - transit_midtime) % orbital_period) / orbital_period
    
    if np.isfinite(phases).all():
        return phases
    else:
        raise ValueError("Observed phases are not finite.")


def get_midtime(input_target, mission):
    """Get transit midtime of innermost planet of the target in BKJD or BTJD.
    
    Parameters:
    -----------
    input_target: pandas.Series
        target info including pl_tranmid and/or pl_tranmid_tess
    mission: string
        mission name, either 'Kepler' or 'TESS'

    Return:
    -------
    midtime: float

    """
    # calculate midtime when there is no transit

    # TBD

    # calculate midtime of transit in TESS or Kepler time
    if mission == "TESS":
        # if transit has been measured in TESS
        if np.isfinite(input_target.pl_tranmid_tess):
            return input_target.pl_tranmid_tess - offset[mission]
        # else transit has not been measured in TESS but only in Kepler
        # still use TESS offset because the light curve is from TESS
        elif np.isfinite(input_target.pl_tranmid):
            return input_target.pl_tranmid - offset[mission]
        else:
            raise ValueError("No transit midtime found. Neither"
                             "pl_tranmid_tess nor pl_tranmid found.")
    # if transit has been measured in Kepler
    elif mission == "Kepler":
        if np.isfinite(input_target.pl_tranmid):
            return input_target.pl_tranmid - offset[mission]
        else:
            raise ValueError("No transit midtime found. No pl_tranmid found.")



def add_meta_data_and_write(ff, dflcn, ID, TIC, sector, mission,
                            lc_n, w, tstamp, mask_pos_outliers_sigma,
                            path="../results/2022_07_flares.csv",
                            header=False):
    """Write out flare table to file.
    
    Parameters:
    -----------

    ff: pd.DataFrame
        flare table
    dflcn: FlareLightCurve
        detrended light curve
    ID: int or string
        star ID
    TIC: int or string
        TIC ID without the leading 'TIC '
    sector: int or string
        Sector or Quarter number
    mission: string
        Kepler or TESS
    lc_n: int
        number of light curves for this star
    w: float
        window size of rolling std
    tstamp: string
        timestamp of this run
    mask_pos_outliers_sigma: float
        number of sigma to mask positive outliers
    path: string
        path to save flare table

    Return:
    -------
    None
    """
    
    
    if ff.shape[0]==0:
        ff["phase"]=-1
        ff["total_n_valid_data_points"] = dflcn.detrended_flux.shape[0]
        ff["ID"] = ID
        ff["TIC"] = TIC
        ff["qcs"] = sector
        ff["mission"] = mission
        ff["tstamp"] = tstamp
        ff["lc_n"] = lc_n
        ff["w"] = w
        ff["mask_pos_outliers_sigma"] = mask_pos_outliers_sigma
        ff["real"]=-1
        ff = ff.append({"phase":-1,
                        "total_n_valid_data_points":dflcn.detrended_flux.shape[0],
                        "ID":ID,
                        "TIC":TIC,
                        "qcs" : sector,
                        "mission":mission,
                        "tstamp":tstamp,
                        "lc_n":lc_n,
                        "w":w,
                        "mask_pos_outliers_sigma":mask_pos_outliers_sigma,
                        "real":-1},
                         ignore_index=True)

    # otherwise add ID, QCS and mission
    else:
        ff["total_n_valid_data_points"] = dflcn.detrended_flux.shape[0]
        ff["ID"] = ID
        ff["TIC"] = TIC
        ff["qcs"] = sector
        ff["mission"] = mission
        ff["tstamp"] = tstamp
        ff["lc_n"] = lc_n
        ff["w"] = w
        ff["mask_pos_outliers_sigma"] = mask_pos_outliers_sigma

    # add results to file
    with open(path, "a") as file:
        ff.to_csv(file, index=False, header=header)
            

def write_flc_to_file(dflcn, flc, path_dflcn, overwrite=True):
    """Write detrended light curve to fits.
    Checks if required columns exist, otherwise
    throws error.
    
    Parameters:
    -----------
    dflcn: FlareLightCurve
        detrended light curve
    flc: lightkurve.LightCurve
        un-detrended light curve
    path_dflcn: string
        path to save detrended light curve
    overwrite: bool
        overwrite existing file, default=True

    Return:
    -------
    None

    """
    
    # check if required columns exist
    if "flux" not in flc.columns:
        raise KeyError("flux not in columns")
    if "detrended_flux" not in dflcn.columns:
        raise KeyError("detrended_flux not in columns")
    if "detrended_flux_err" not in dflcn.columns:
        raise KeyError("detrended_flux_err not in columns")
    if "time" not in dflcn.columns:
        raise KeyError("time not in columns")
    if "it_med" not in dflcn.columns:
        raise KeyError("it_med not in columns")
    if "phase" not in dflcn.columns:
        raise KeyError("phase not in columns")
    if "flux_model" not in dflcn.columns:
        raise KeyError("flux_model not in columns")    

    dflcn.to_fits(path_dflcn, 
                  FLUX=flc.flux.value.astype(float),
                  DETRENDED_FLUX=dflcn.detrended_flux.value.astype(float),
                  DETRENDED_FLUX_ERR=dflcn.detrended_flux_err.value.astype(float),
                  IT_MED=dflcn.it_med.value.astype(float),
                  FLUX_MODEL=dflcn.flux_model.value.astype(float),
                  PHASE = dflcn.phase.astype(float),
                  overwrite=True)
                  

def write_no_lc(input_target, path="../results/2022_07_nolc.txt"):
    """Write TIC to file if no light curve is found.
    
    Parameters:
    -----------
    input_target: pd.Series
        input target with TIC column
    path: string
        path to save file

    Return:
    -------
    None
    """
    # string to write
    s = f"{input_target.TIC}\n"

    # check if file at path exists
    # if not write header "TIC" to path
    if os.path.exists(path):
        with open(path, "a") as f:
            f.write(s)  
    else:
        with open(path, "w") as f:
            f.write("TIC\n")
            f.write(s)
        

def get_table_of_light_curves(input_target, path="../results/2022_07_nolc.txt"):
    """Get table of light curves for a given target using lightkurve query.
    
    Parameters:
    -----------
    input_target: pd.Series
        input catalog row
    path: string
        path to save file
    
    Return:
    -------
    lc_table: pd.DataFrame
    """
    try:
        # search for light curves
        lcs  = search_lightcurve(input_target.hostname) 
        # only keep short and fast cadence
        # no K2 data, no TASOC asteroseismic light curves
        conditions = (lcs.exptime.value < 130)  & (lcs.author != "TASOC") & (lcs.author != "K2")
        lc_table = lcs[conditions]
        
    except KeyError:
     
        try:
            # search for light curves using TESS TIC if ID fails
            lcs  = search_lightcurve(f"TIC {input_target.TIC}")
            
            # only keep short and fast cadence
            # no K2 data, no TASOC asteroseismic light curves
            conditions = (lcs.exptime.value < 130) & (lcs.author != "TASOC") & (lcs.author != "K2")
            lc_table = lcs[conditions]
        
        except KeyError:
            # if no light curve is found, write TIC to file
            write_no_lc(input_target, path=path)
            return None

    # if lc_table is empty, write TIC to file and return
    if len(lc_table)==0:
        write_no_lc(input_target, path=path)
        return None
   
    # if lc_table is not empty, convert to pandas dataframe 
    # and sort by increasing exposure time
    lc_table = lc_table.table.to_pandas().sort_values(by="t_exptime",
                                                      ascending=True)
    
    # select only the first light curve in the TESS subset, i.e.
    # the one with the shortest cadence
    lcs_sel_tess = lc_table.loc[lc_table.mission.str[:4]=="TESS",:]
    lcs_sel_tess = lcs_sel_tess.drop_duplicates(subset=["mission"], 
                                                keep="first")
    # select only the first light curve in the Kepler subset, i.e.
    # the one with the shortest cadence
    lcs_sel_kepler = lc_table[lc_table.mission.str[:6]=="Kepler"]
    lcs_sel_kepler = lcs_sel_kepler.drop_duplicates(subset=["mission"])

    # combine the two tables
    lc_table = pd.concat([lcs_sel_kepler,lcs_sel_tess])
   
    return lc_table
