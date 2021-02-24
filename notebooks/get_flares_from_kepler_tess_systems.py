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
from funcs.detrend import estimate_detrended_noise, custom_detrending
from funcs.transitmask import get_full_transit_mask

from altaipony.lcio import from_mast
from altaipony.flarelc import FlareLightCurve


sep = "-----------------------------------------"

def mprint(message):
    print(sep)
    print(message)
    print(sep)
    
offset = {"K2":2454833.,"Kepler":2454833.,"TESS":2457000., "Transiting Exoplanet Survey Satellite (TESS)" : 2457000.}    


if __name__ == "__main__":
    
    # Composite Table of confirmed exoplanets
    path = "20_01_2021_confirmed_uncontroversial_exoplanet_systems.csv"
  
    mprint(f"[UP] Using confirmed and uncontroversial "
          "entries in NASA Composite Table from {path}")
    
    exokepler = pd.read_csv(f"../data/{path}") # composite table

    # read in TESS-TOI sample 
    path = "../data/2021_01_13_TESS_TOI_CATALOG.csv"
    
    mprint(f"[UP] Using TESS-TOI Table from {path}")
    
    exotess = pd.read_csv(path, skiprows=4)

    # rename the relevant columns for transit masking
    exotess = exotess.rename(index=str, 
                             columns={'Transit Duration Value':"pl_trandur",
                                      'Orbital Period Value': "pl_orbper", 
                                      'Epoch Value':"pl_tranmidepoch"})


    # read in list of LCs to search
    es = pd.read_csv("../data/20_01_2021_full_kepler_k2_tess_exoplanet_lcs_some_excluded.csv")

    # select only Kepler and TESS, ignore K2 for now
    eskeptess = es[(es.mission=="TESS") | (es.mission=="Kepler")]

    # read in searched LCs with and without flares
    fla = pd.read_csv("../results/flares.csv")

    # pick only LC that were not yet searched
    eskeptess = eskeptess[~eskeptess.ID.isin(fla.ID.unique())]

    mprint(f"LC left to search: {eskeptess.shape[0]}")

    #work through a subset first
    eskeptess = eskeptess.iloc[160:360]
    
    #track progress
    N, n = eskeptess.shape[0], 0
    
    # take how long it takes to process
    TSTART = time.time()
    
    for i, row in eskeptess.iterrows():

        # TIC is unique ID for star
        system_tess = exotess[(exotess.TIC == row.TIC)]
        
        # ID is unique, also ignore entries that have no transits
        # because there is nothing to mask (they are still searched for flares)
        system_kepler = exokepler[(exokepler.hostname == row.ID) &
                           (exokepler.discoverymethod == "Transit")]
        

        if system_kepler.shape[0] > 0:
            system_kepler["pl_tranmidepoch"] = (system_kepler.pl_tranmid -
                                                offset[system_kepler.iloc[0].disc_facility])

        system = pd.concat([system_kepler, system_tess],ignore_index=True)

        # fetch light curve from MAST
        flc = from_mast(row.ID, mission=row.mission, c=row.qcs, cadence="short",
                        download_dir="/home/ekaterina/Documents/001_science/lcs")

        # make it a list of LCs even if only one LC is returned
        if type(flc) == FlareLightCurve:

            flc = [flc]

        elif type(flc) == list:

            flc = flc

        # info
        mprint(f"{len(flc)} light curves available for {row.ID} in {row.mission}.")

        # loop over all LCs for the system    
        for i, f in enumerate(flc):

            # If any planet transiting
            if system.shape[0] > 0:

                # mask transits
                tranmask = get_full_transit_mask(system, f, pad=0)
                f.flux[tranmask] = np.nan

            # get timestamp for result
            tstamp = time.strftime("%d_%m_%Y_%H_%M_%S", time.localtime())
            # apply custom detrending
            try:
                ts = time.clock()
                fd = custom_detrending(f)
                tf = time.clock()


                # define two hour window for rolling std
                w = np.floor(1. / 12. / np.nanmin(np.diff(fd.time)))
                if w%2==0: 
                    w+=1

                # use window to estimate the noise in the LC
                df = estimate_detrended_noise(fd, std_window=int(w), mask_pos_outliers_sigma=2.5)

                # search the residual for flares
                ff = df.find_flares(addtail=True).flares

                # add meta info to flare table
                # if no flares found, add empty row
                if ff.shape[0]==0:
                    ff["total_n_valid_data_points"] = df.detrended_flux.shape[0]
                    ff["ID"] = row.ID
                    ff["qcs"] = row.qcs
                    ff["mission"] = row.mission
                    ff["tstamp"] = tstamp
                    ff["dur_detrend"] = tf - ts
                    ff["lc_n"] = i
                    ff = ff.append({"total_n_valid_data_points":df.detrended_flux.shape[0],
                                    "ID":row.ID,
                                    "qcs" : row.qcs,
                                    "mission":row.mission,
                                    "tstamp":tstamp,
                                    "dur_detrend":tf-ts,
                                    "lc_n":i},
                                     ignore_index=True)

                # otherwise add ID, QCS and mission
                else:
                    ff["ID"] = row.ID
                    ff["qcs"] = row.qcs
                    ff["mission"] = row.mission
                    ff["tstamp"] = tstamp
                    ff["dur_detrend"] = tf - ts
                    ff["lc_n"] = i

                # add results to file
                with open("../results/flares.csv", "a") as file:
                    ff.to_csv(file, index=False, header=False)
            except Exception as err:
                mprint(f"{row.ID}, QCS={row.qcs} ({row.mission}) not de-trended!")
                with open("../results/nodetrend.txt", "a") as file:
                    s = f"{row.ID},{row.qcs},{row.mission},{tstamp},{i},{err}\n"
                    file.write(s)

        # info
        n += 1
        print(f"{n / N * 100.:.1f}%, [{n}/{N}]")

        # breathe
        time.sleep(10)
            
    TSTOP = time.time()
    mprint(f"Analysis of {N} light curves took {(TSTOP - TSTART) / 60. / 60.:.1f} hours.")
