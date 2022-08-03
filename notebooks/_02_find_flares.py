"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2021, MIT License

De-trending Kepler and TESS

- get table of Kepler exoplanet system light curves
- fetch FLC
- get system info from table
- mask transits
- apply custom detrending
- search flares
- save results
"""

import copy
import time

from funcs.notebook import *
from funcs.detrend import (custom_detrending,
                           estimate_detrended_noise)


from altaipony.lcio import from_mast, from_path
from altaipony.flarelc import FlareLightCurve
from altaipony.altai import find_iterative_median

from lightkurve import search_lightcurvefile

from astropy.io import fits


def add_meta_data_and_write(ff, dflcn, ID, TIC, sector, mission,
                            lc_n, w, tstamp, mask_pos_outliers_sigma):
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
    with open("../results/2022_07_flares.csv", "a") as file:
        ff.to_csv(file, index=False, header=False)
            

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
                  FLUX=flc.flux.value,
                  DETRENDED_FLUX=dflcn.detrended_flux.value,
                  DETRENDED_FLUX_ERR=dflcn.detrended_flux_err.value,
                  IT_MED=dflcn.it_med.value,
                  FLUX_MODEL=dflcn.flux_model.value,
                  PHASE = dflcn.phase,
                  overwrite=True)

def write_no_lc(input_target):
    """Write TIC to file if no light curve is found.
    
    Parameters:
    -----------
    input_target: pd.Series
        input target with TIC column

    Return:
    -------
    None
    """
    with open("../results/2022_07_nolc.txt","a") as f:
        s = f"TIC {input_target.TIC}\n"
        f.write(s)
    
def run_analysis(flc, input_target, sector, mission, lc_n, download_dir,
                 i=0, mask_pos_outliers_sigma = 2.5, addtail = True):

    """Run de-trending and flare finding on a single light curve. 
    Get orbital phases for all observed times and flare times as well.
    Write out resulting light curve and flare table to file.

    Parameters:
    -----------
    flc: FlareLightCurve
        un-detrended light curve
    input_target: pd.Series
        target dataframe row with orbital period,
        transit midtime, ID
    sector: int or string
        Sector or Quarter number
    mission: string
        Kepler or TESS
    lc_n: int
        number of light curve for this star
    download_dir: string
        path to download directory, which is also where the
        detrended light curve is saved
    i: int
        index of light curve, if multiple light curves are 
        found per Quarter
    mask_pos_outliers_sigma: float
        number of sigma to mask positive outliers
    addtail: bool
        add tail of less than 3-sigma outliers to flare
    
    Return:
    -------
    number of flares found


    
    """    
    # get timestamp for result
    tstamp = time.strftime("%Y_%m_%d", time.localtime())
    print(f"date: {tstamp}")

    dflc = custom_detrending(flc)
    print("LC successfully detrended.")

    # define two hour window for rolling std
    w = np.floor(1. / 12. / np.nanmin(np.diff(dflc.time.value)))
    if w%2==0: 
        w+=1

    # use window to estimate the noise in the LC
    dflcn = estimate_detrended_noise(dflc, std_window=int(w), 
                                  mask_pos_outliers_sigma=mask_pos_outliers_sigma)

    # search the residual for flares
    ff = dflcn.find_flares(addtail=addtail).flares


    # calculate the observed phases
    # calculate midtime of transit in TESS or Kepler time
    if mission == "TESS":
        if np.isfinite(input_target.pl_tranmid_tess):
            midtime = input_target.pl_tranmid_tess - offset[mission]
        else:
            midtime = input_target.pl_tranmid - offset[mission]
    elif mission == "Kepler":
        midtime = input_target.pl_tranmid - offset[mission]
    print(f"Transit midtime in {mission} time: {midtime}")

    # calculate phases for the light curve
    dflcn['phase'] = ((dflcn.time.value - midtime) % input_target.pl_orbper) / input_target.pl_orbper

    # calculate the phase at which the flare was observed
    ff["phase"] = ff.cstart.apply(lambda x: dflcn["phase"][np.where(x==dflcn.cadenceno)][0])
    

    # this is just to get the order of columns right, will be added later again
    if ff.shape[0]>0:
        del ff["total_n_valid_data_points"]

    # chop out all phases where we have no data points to look 
    # for flares in:
    dflcn["phase"][~np.isfinite(dflcn["detrended_flux"])] = np.nan

    # Write result to screen
    fshow = ff[["tstart",'tstop',"phase","ampl_rec","dur"]]
    if fshow.shape[0]>0:
        print(f"Flares found:\n{fshow}")
    else:
        print(f'No flares found in LC.')

    # add meta info to flare table
    # if no flares found, add empty row and write to file
    add_meta_data_and_write(ff, dflcn, input_target.hostname, 
                            input_target.TIC, sector,
                            mission, lc_n, w, tstamp,
                            mask_pos_outliers_sigma)


    #write out detrended light curve
    if mission=="TESS":
        path_dflcn = f"{download_dir}/{tstamp}_{input_target.TIC}_{sector}_altai_{i}.fits"
    elif mission=="Kepler":
        name = input_target.hostname.replace(" ","_").replace("-","_")
        path_dflcn = f"{download_dir}/{tstamp}_{input_target.hostname}_{sector}_altai_{i}.fits"
        
    # write light curve to file and notify user
    write_flc_to_file(dflcn, flc, path_dflcn)
    print(f"Wrote out LC to {path_dflcn}.")

    return ff.shape[0]

def get_table_of_light_curves(input_target):
    """Get table of light curves for a given target using lightkurve query.
    
    Parameters:
    -----------
    input_target: pd.Series
        input catalog row
    
    Return:
    -------
    lc_table: pd.DataFrame
    """
    try:
        # search for light curves
        lcs  = search_lightcurvefile(input_target.hostname) 
        # only keep short and fast cadence
        # no K2 data, no TASOC asteroseismic light curves
        conditions = (lcs.exptime.value < 130)  & (lcs.author != "TASOC") & (lcs.author != "K2")
        lc_table = lcs[conditions]
    
    except KeyError:
     
        try:
            # search for light curves using TESS TIC if ID fails
            lcs  = search_lightcurvefile(f"TIC {input_target.TIC}")

            # only keep short and fast cadence
            # no K2 data, no TASOC asteroseismic light curves
            conditions = (lcs.exptime.value < 130) & (lcs.author != "TASOC") & (lcs.author != "K2")
            lc_table = lcs[conditions]
        
        except KeyError:
            # if no light curve is found, write TIC to file
            write_no_lc(input_target)
            return

    # if lc_table is empty, write TIC to file and return
    if len(lc_table)==0:
        write_no_lc(input_target)
        return
   
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

# separator for logging    
sep = "-----------------------------------------"

def mprint(message):
    """Pretty print message."""
    print(sep)
    print(message)
    print(sep)
    
# offsets for BTJD and BKJD
offset = {"K2":2454833.,
          "Kepler":2454833.,
          "TESS":2457000.,
          'Transiting Exoplanet Survey Satellite (TESS)':2457000.}    


if __name__=="__main__":

    # Composite Table of confirmed exoplanets
    path = "../data/2022_07_27_input_catalog_star_planet_systems.csv"

    mprint(f"[UP] Using compiled input catalog from {path}")

    input_catalog = pd.read_csv(path) 

    # init counter for target
    counter = 0

    # init counter for light curve
    Nflares = 0 

    # define download directory
    download_dir = "/home/ekaterina/Documents/001_science/lcs"

    # stop if too many flares are found, bc it's sus
    while Nflares < 100:

        # init analysis
        print(f"\nCOUNT: {count}\n")
        lcs_sel=pd.DataFrame()
        
        # move while no light curve is found
        while lcs_sel.shape[0]==0:

            # get next target
            input_target = input_catalog.iloc[count]
            
            # get table of light curves
            lcs_sel = get_table_of_light_curves(input_target)

            # if no light curve is found, move to next target
            if lcs_sel is None:
                lcs_sel=pd.DataFrame()
            
            # advance counter to move to next target
            count+=1

        # format TIC
        TIC = "TIC " + str(input_target.TIC)
        ID = input_target.hostname

        # loop over the table of light curves
        n = 0

        # reset number of flares
        Nflares = 0

        # move until all light curves are analyzed
        while n < lcs_sel.shape[0]:
            
            # get sector and mission
            sector = lcs_sel.iloc[n].mission[-2:]
            mission = lcs_sel.iloc[n].mission.split(" ")[0]

            # light curve number starts with 1
            lc_n = n + 1

            # disambiguate short and fast cadence
            if lcs_sel.iloc[n].exptime < 30:
                cadence = "fast"
            else: 
                cadence = "short"

            print(f"Get {mission} Sector/Quarter {sector}, {TIC}, {ID}, {cadence} cadence.")

            # fetch light curve from MAST and analyze

            if mission=="TESS":

                flc = from_mast(TIC, mission=mission, c=sector,
                            cadence=cadence, author="SPOC",
                            download_dir=download_dir)

                # handle case when no light curve is found
                if flc is None:
                    print(f"No LC found for {mission}, {ID}, Quarter {sector}.")
                    with open("../results/2022_07_listed_but_nothing_found.txt", "a") as f:
                        string = f"{mission},{ID},{TIC},{sector},{cadence}\n"
                        f.write(string)
                    # advance to next light curve
                    n += 1
                # otherwise, analyze light curve
                else:
                    Nflares += run_analysis(flc, input_target, sector, mission, lc_n, download_dir, i=0)
                    # advance to next light curve
                    n += 1

            elif mission=="Kepler":

                flcl = from_mast(ID, mission=mission, c=sector,
                            cadence=cadence,
                            download_dir=download_dir)

                # handle case when no light curve is found
                if flcl is None:
                    print(f"No LC found for {mission}, {ID}, Quarter {sector}.")
                    with open("../results/2022_07_listed_but_nothing_found.txt", "a") as f:
                        string = f"{mission},{ID},{TIC},{sector},{cadence}\n"
                        f.write(string)
                        
                    # advance to next light curve
                    n += 1

                # handle case when only one light curve is found
                elif type(flcl) != list:

                    print(f"1 LC found for {mission}, {ID}, Quarter {sector}.")
                    Nflares += run_analysis(flcl, input_target, sector, mission, lc_n, download_dir, i=0)
                    
                    # advance to next light curve
                    n += 1

                # handle case when multiple light curves are found
                else:
                    print(f"{len(flcl)} LCs found for {mission}, {ID}, Quarter {sector}.")
                    # analyze all light curves
                    for i, flc in enumerate(flcl):
                        Nflares += run_analysis(flc, input_target, sector, mission, lc_n, download_dir, i=i)
                    
                    # advance to next light curve
                    n += 1

        print(f"\n---------------------\n{Nflares} flares found!\n-------------------\n")
    print(f"\nNext count is {count}.\n")



