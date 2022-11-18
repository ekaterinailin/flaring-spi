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

import time

from funcs.notebook import *
from funcs.detrend import (custom_detrending,
                           estimate_detrended_noise)

from funcs.findflares import (get_observed_phases, get_flare_phases,
                              get_midtime,
                              add_meta_data_and_write,
                              write_flc_to_file,
                              get_table_of_light_curves,
                            )

from altaipony.lcio import from_mast
 


    
def run_analysis(flc, input_target, sector, mission, lc_n, download_dir,
                 i=0, mask_pos_outliers_sigma=2.5, addtail = True):

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
    # -------------------------------------------------------------------------
    # get timestamp for result
    tstamp = time.strftime("%Y_%m_%d", time.localtime())
    print(f"date: {tstamp}")

    # -------------------------------------------------------------------------
    # detrend light curve
    dflc = custom_detrending(flc)
    print("LC successfully detrended.")

    # -------------------------------------------------------------------------
    # get an estimate for flux uncertainty

    # define two hour window for rolling std
    w = np.floor(1. / 12. / np.nanmin(np.diff(dflc.time.value)))
    if w%2==0: 
        w+=1

    # use window to estimate the noise in the LC
    dflcn = estimate_detrended_noise(dflc, std_window=int(w), 
                                  mask_pos_outliers_sigma=mask_pos_outliers_sigma)

    # -------------------------------------------------------------------------
    # search the residual for flare candidates
    ff = dflcn.find_flares(addtail=addtail).flares

    # -------------------------------------------------------------------------
    # get orbital phases for all observed times and the flare candidate times
    # get midtime in BTJD or BKJD, depending on mission

    # THIS IS ONLY SET TO KEPLER because not TESS DATA EXIST
    # midtime = get_midtime(input_target, mission)
    # print(f"Transit midtime in {mission} time: {midtime}")

    # # calculate phases for the light curve
    # # we always take the Kepler one, as it uses the longer baseline, and is
    # # otherwise filled in from the TESS one in input_catalog
    # dflcn['phase'] = get_observed_phases(dflcn.time.value, midtime, 
    #                                      input_target.pl_orbper)

    # # calculate the phase at which the flare was observed
    # ff['phase'] = get_flare_phases(ff.cstart, dflcn["phase"], 
    #                                dflcn["cadenceno"])

    dflcn["phase"] = -1.
    ff["phase"] = -1.

    # -------------------------------------------------------------------------

    # The next line is just to get the order of columns right, 
    # will be added later again
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

    # -------------------------------------------------------------------------
    # add meta info to flare table
    # if no flares found, add empty row and write to file
    add_meta_data_and_write(ff, dflcn, input_target.hostname, 
                            input_target.TIC, sector,
                            mission, lc_n, w, tstamp,
                            mask_pos_outliers_sigma)

    # -------------------------------------------------------------------------
    # write out detrended light curve with observed phases
    if mission=="TESS":
        path_dflcn = f"{download_dir}/{tstamp}_{input_target.TIC}_{sector}_altai_{i}.fits"
    elif mission=="Kepler":
        name = input_target.hostname.replace(" ","_").replace("-","_")
        path_dflcn = f"{download_dir}/{tstamp}_{name}_{sector}_altai_{i}.fits"
        
    # write light curve to file and notify user
    write_flc_to_file(dflcn, flc, path_dflcn)
    print(f"Wrote out LC to {path_dflcn}.")

    return ff.shape[0]


def add_helpid(x):
    """Add help ID to light curve table.
    
    Parameters:
    -----------
    x: pd.Series
        input table row        

    Return:
    -------
    string with help ID
    """
    if x.mission == "Kepler":
        return f"{x.ID}_{x.mission}_{int(x.qcs)}"
    elif x.mission == "TESS":
        return f"{x.TIC}_{x.mission}_{int(x.qcs)}"


if __name__=="__main__":

    mask_pos_outliers_sigma = 2.5

    # Composite Table of confirmed exoplanets
    # path = "../data/2022_07_27_input_catalog_star_planet_systems.csv"
    # path = "../data/2022_08_16_input_catalog_systems_with_KOI_IDs.csv"
    path = "../data/2022_11_15_input_catalog_NONtransit_star_planet_systems.csv"

    mprint(f"[UP] Using compiled input catalog from {path}")

    input_catalog = pd.read_csv(path) 

    # define download directory
    download_dir = "/home/ekaterina/Documents/001_science/lcs"

    # read in table of already found flares
    found_flares = pd.read_csv("../results/2022_07_flares.csv")

    # make an ID list of all stars that have already been processed
    found_flares["helpid"] = found_flares.apply(add_helpid, axis=1)

    # stop if too many flares are found, bc it's sus
    for n, input_target in input_catalog.iterrows():

        # init analysis
        print(f"\nCOUNT: {n}\n")

        print(f"Target: {input_target.hostname}")
        
        lcs_sel = get_table_of_light_curves(input_target)
        print(lcs_sel)
        # if no light curve is found, move to next target
        if lcs_sel is None:
            continue

        # format TIC
        TIC = "TIC " + str(input_target.TIC)
        ID = input_target.hostname
        Nflares = 0

        # move until all light curves are analyzed
        for lc_idx, lc_row in lcs_sel.iterrows():
            
            # get sector and mission
            sector = lc_row.mission[-2:]
            mission = lc_row.mission.split(" ")[0]
            print(f"\nSECTOR: {sector}\n")

            # light curve number starts with 1
            lc_n = lc_idx + 1

            # disambiguate short and fast cadence
            if lc_row.exptime < 30:
                cadence = "fast"
            else: 
                cadence = "short"

            # Check whether this light curve has already been searched
            if mission == "Kepler":
                id_str = f"{ID}_{mission}_{int(sector)}"
            elif mission == "TESS":
                id_str = f"{TIC[4:]}_{mission}_{int(sector)}"

            print(f"\nID string: {id_str}")

            # If not, start analysis:
            if len(np.where(found_flares["helpid"].str.contains(id_str))[0]) == 0:

                print(f"Get {mission} Sector/Quarter {sector}, {TIC}, {ID}, "
                f"{cadence} cadence.")
          
                # fetch light curve from MAST and analyze
                parameters = {"c":int(sector), "mission":mission,
                            "cadence":cadence,"download_dir":download_dir}
                print(f"SECTOR/QUARTER {sector}")
                if mission=="TESS":
                    flc = from_mast(TIC, author="SPOC", **parameters)
                elif mission=="Kepler":
                    flc = from_mast(ID, **parameters)

                
                # handle case when no light curve is found
                if flc is None:
                    print(f"No LC found for {mission}, {ID}, {TIC} "
                        f"Quarter/Sector {sector}.")
                    with open("../results/2022_07_listed_but_nothing_found.txt", "a") as f:
                        string = f"{mission},{ID},{TIC},{sector},{cadence}\n"
                        f.write(string)
                    

                # otherwise, analyze light curve(s)
                elif type(flc) != list:

                        print(f"1 LC found for {mission}, {ID}, Quarter {sector}.")
                        Nflares += run_analysis(flc, input_target, sector, mission,
                                                lc_n, download_dir, i=0, 
                                                mask_pos_outliers_sigma=mask_pos_outliers_sigma)

                # handle case when multiple light curves are found
                else:
                    print(f"{len(flc)} LCs found for {mission}, {ID}, Quarter {sector}.")
                    # analyze all light curves
                    for i, lc in enumerate(flc):
                        Nflares += run_analysis(lc, input_target, sector, mission,
                                                lc_n, download_dir, i=i, 
                                                mask_pos_outliers_sigma=mask_pos_outliers_sigma)
            else:
                print(f"\nAlready searched {mission}, {ID}, Quarter/Sector {sector}.\n")
                continue
        print(f"\n---------------------\n{Nflares} flares found!\n-------------------\n")

    print(f"\nNext input target is row index {n+1}.\n")
