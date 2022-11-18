"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script uses the Composite Table of confirmed exoplanets
from the NASA Exoplanet Archive* (column description**)
to compile a sample of all Kepler/K2/TESS short cadence light curves available 
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

    # select only uncontroversial detections of planets 
    # that have orbital periods shorter than 50 days
    conditions = ((df_nasa_full.pl_controv_flag == 0) &
                  (df_nasa_full.pl_orbper < 50.))
    sel = df_nasa_full[conditions]

    # select the indices of the innermost planets
    sel = sel.sort_values("pl_orbper",ascending=True)
    df_nasa = sel.groupby("tic_id").first().reset_index()

    # now select those that DO NOT TRANSIT
    df_nasa = df_nasa[df_nasa.tran_flag == 0]

    # rename tic id column to TIC
    df_nasa = df_nasa.rename(columns={"tic_id":"TIC"})
    df_nasa["TIC"] = df_nasa["TIC"].str[4:].astype(int)
    print(df_nasa.TIC.values)

    path =f"../data/{tstamp}_input_catalog_NONtransit_star_planet_systems.csv"
    print(f"[DOWN] Saving {df_nasa.shape[0]} uncontroversial "
           f"NON-transiting exosystems from NASA Composite Table to {path}")
    df_nasa.to_csv(path, index=False)
