"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script uses the Composite Table of confirmed exoplanets
from the NASA Exoplanet Archive* (column description**) AND 
the TESS-TOI table of confirmed planets and known planets***
(column description****) to compile a sample of all 
Kepler/K2/TESS short cadence light curves available 
for all currently known exoplanet systems.

* https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PSCompPars
** https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
*** https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=TOI
**** https://exoplanetarchive.ipac.caltech.edu/docs/API_toi_columns.html
"""

import pandas as pd
import time

if __name__ == '__main__':

    tstamp = time.strftime("%Y_%m_%d", time.localtime())
    print(f"Timestamp: {tstamp}")



    # --------------------------------------------------
    # READ IN NASA EXOPLANET ARCHIVE TABLE
    # Composite Table of confirmed exoplanets
    path = "../data/2022_07_26_NASA_COMPOSITE.csv"
    print(f"[UP] Using NASA Composite Table from {path}")
    df_nasa_full = pd.read_csv(path, skiprows=320) # composite table

    # select only uncontroversial detections of planets known to transit
    conditions = ((df_nasa_full.pl_controv_flag==0) &
                  (df_nasa_full.tran_flag==1))
    sel = df_nasa_full[conditions]

    # select the indices of the innermost planets
    sel = sel.sort_values("pl_orbper",ascending=True)
    df_nasa = sel.groupby("tic_id").first().reset_index()

    path =(f"../data/{tstamp}_confirmed_uncontroversial_"  
          f"innermost_transiting_exoplanet.csv")
    print(f"[DOWN] Saving {df_nasa.shape[0]} uncontroversial "
           f"transiting exosystems from NASA Composite Table to {path}")
    df_nasa.to_csv(path, index=False)

    # --------------------------------------------------

    # --------------------------------------------------
    # READ in TESS-TOI CATALOG

    path = "../data/2022_07_26_TESS_TOI_CATALOG.csv"
    print(f"[UP] Using TESS-TOI Table from {path}")
    df_tess_full = pd.read_csv(path, skiprows=90)

    # select only known candidates and confirmed planets
    # KP = known planet, CP = confirmed planet
    df1 = df_tess_full[df_tess_full["tfopwg_disp"].isin(["KP","CP"])]

    # sort by orbital period
    df1 = df1.sort_values("pl_orbper",ascending=True)

    # select innermost planet
    df_tess  = df1.groupby("tid").first().reset_index()
    df_tess = df_tess.rename(index=str, columns={"tid":"tic_id"})

    # save to file
    path = f"../data/{tstamp}_tess_confirmed_and_known_planets.csv"
    print(f"[DOWN] Saving CP and KP sample of {df_tess.shape[0]} "
        f"planet hosts from TESS TOI to {path}")
    df_tess.to_csv(path, index=False)

    # --------------------------------------------------


    # --------------------------------------------------
    # MERGE BOTH TABLES on TIC ID

    # Note: 
    # All NASA ARCHIVE entries have a TIC ID expect 
    # for some RV, direct imaging, and microlensing targets
    # and two transiting planets in the galactic bulge.

    # the TIC ID should look the same in both data frames
    df_nasa["TIC"] = df_nasa.tic_id.str[4:]
    df_tess["TIC"] = df_tess.tic_id.astype(str)

    # select final columns for the merged table
    columns = ["pl_orbper","pl_orbpererr1","pl_orbpererr2",
            "pl_tranmid","pl_tranmiderr1", "pl_tranmiderr2"]

    df_merged = df_nasa.merge(df_tess, how="outer",on=["TIC"],
                              suffixes=["","_tess"])

    # --------------------------------------------------

    # --------------------------------------------------
    # SELECT FINAL COLUMNS and FILL IN MISSING VALUES IN
    # NASA table from TESS TABLE

    # select columns
    df_final = df_merged[["TIC"] +
                        ["hostname","pl_name","sy_pnum","sy_snum"] +
                        columns +
                        [c+"_tess" for c in columns]]

    # fill in TESS orbital parameter if NASA table is missing them
    condition = ((df_final["pl_orbper"].isnull()) &
                 (~df_final["pl_orbper_tess"].isnull()))

    # set flag if values are filled in from TESS table
    df_final["obrparams_tess"] = 0
    df_final.loc[condition,"obrparams_tess"] = 1

    # fill in each column
    for col in columns:
        df_final[col].fillna(df_final[col+"_tess"], inplace=True)
        
    # --------------------------------------------------

    # --------------------------------------------------
    # WRITE RESULT TO FILE

    df_final.to_csv(f"../data/{tstamp}_input_catalog_star_planet_systems.csv",
                    index=False)

    # --------------------------------------------------
