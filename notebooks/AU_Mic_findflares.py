import copy
import time
import sys

from funcs.notebook import *
from funcs.detrend import estimate_detrended_noise, custom_detrending

from altaipony.lcio import from_path
from altaipony.flarelc import FlareLightCurve

def write_log(s):
    print(s)
    with open("results/log.txt", "a") as file:        
        file.write(s + "\n")

if __name__ == "__main__":
    
    
    
    path = sys.argv[1]

    print(path)

    # get timestamp for result
    tstamp = time.strftime("%d_%m_%Y_%H_%M_%S", time.localtime())

    TIC = int(path.split("-0000")[1].split("-")[0])
    sector = int(path.split("-s00")[1].split("-00000")[0])
    
    try:

        # fetch light curves from vanir
        write_log(f"{tstamp} Fetching TIC {TIC} (Sector {sector}).")
        lc = from_path(path, mode="LC", mission="TESS")

        # apply custom detrending

        # try:

        # start timer
        ts = time.process_time()

        # run de-trending
        fd, ax = custom_detrending(lc, spline_coarseness=30)


        # stop timer
        tf = time.process_time()

        # log that de-trending worked
        write_log(s=f"{tstamp}:: TIC {TIC}, Sector {sector} was "
                    f"successfully de-trended.")

        # define two hour window for rolling std
        w = np.floor(1. / 12. / np.nanmin(np.diff(fd.time.value)))
        if w%2==0: 
            w+=1

        # use window to estimate the noise in the LC
        df = estimate_detrended_noise(fd, std_window=int(w), 
                                      mask_pos_outliers_sigma=2.5)

        # search the residual for flares
        ff = df.find_flares(addtail=True).flares

        # log that flare finding worked
        write_log(s=f"{tstamp}:: TIC {TIC}, Sector {sector} has "
                    f"{ff.shape[0]} flare candidates.")

        # add meta info to flare table
        # if no flares found, add empty row
        if ff.shape[0]==0:
            ff["total_n_valid_data_points"] = df.detrended_flux.value.shape[0]
            ff["ID"] = TIC
            ff["qcs"] = sector
            ff["tstamp"] = tstamp
            ff["dur_detrend"] = tf - ts

            ff = ff.append({"total_n_valid_data_points":df.detrended_flux.value.shape[0],
                            "ID":TIC,
                            "qcs" : sector,
                            "tstamp":tstamp,
                            "dur_detrend":tf-ts},
                             ignore_index=True)

        # otherwise add ID, QCS and mission
        else:
            ff["ID"] = TIC
            ff["qcs"] = sector
            ff["tstamp"] = tstamp
            ff["dur_detrend"] = tf - ts

        # add results to file
        with open("results/flares.csv", "a") as file:
            ff.to_csv(file, index=False, header=False)
            write_log(s = f"{tstamp}:: TIC {TIC}, Sector {sector} saved with {ff.shape[0]} lines.")

        # plot diagnostic
        ax.plot(df.time.value,df.it_med+df.flux_err.value)
        ax.set_xlim(lc.time.value[0], lc.time.value[-1])
        ax.set_xlabel("time [BTJD]", fontsize=14)
        ax.set_ylabel(r"flux e$^-$s$^{-1}$", fontsize=14)
        ax.legend(frameon=False, fontsize=14)    
        plt.savefig(f"diag_plots/{tstamp}_{lc.targetid}_{lc.sector}.png", dpi=300)
        write_log(s = f"{tstamp}:: TIC {TIC}, Sector {sector} de-trended LC diagnostic plot saved.")
        
    except Exception as e: 
        write_log(s = f"{tstamp}:: TIC {TIC}, Sector {sector} failed with error {e}.")