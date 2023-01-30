"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2022, MIT License

This script fills in missing luminosities and their errors for stars with AD
tests, using Gaia DR3 data. If a luminosity is given, and the error is not,
we still replace the NASA Exoplanet archive value with the Gaia DR3 value.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astroquery.gaia import Gaia



if __name__ == "__main__":

    # log in to Gaia
    Gaia.login()

    # read in SPS table
    path = "../results/params_of_star_planet_systems_with_AD_tests.csv"
    print(f"[UP] Loading stellar parameters of stars with AD tests from {path}")
    df = pd.read_csv(path)

    print(len(df))
    for col in ["lum_flame", "lum_flame_lower", "lum_flame_upper", "gaiadr3_id"]:
        if col in df.columns:
            del df[col]

    # fill in NaN in ra_kepler and dec_kepler with ra_tess and dec_tess
    df.ra_kepler = df.ra_kepler.fillna(df.ra_tess)
    df.dec_kepler = df.dec_kepler.fillna(df.dec_tess)

    # delete the old file in Gaia with the same name
    job = Gaia.delete_user_table("spistars")    

    # upload the table to Gaia
    job = Gaia.upload_table(upload_resource=path, table_name="spistars", format="csv")

    # first match the ra and dec columns with gaiadr3.gaia_source within
    # 0.5 arcsec, then match the source_id with gaiadr3.astrophysical_parameters

    query = """
    SELECT
        spistars.hostname,
        gaiadr3.astrophysical_parameters.source_id,
        gaiadr3.astrophysical_parameters.lum_flame,
        gaiadr3.astrophysical_parameters.lum_flame_lower,
        gaiadr3.astrophysical_parameters.lum_flame_upper

    FROM
        spistars
    INNER JOIN  gaiadr3.gaia_source
    ON
        1=CONTAINS(
            POINT('ICRS', spistars.ra_kepler, spistars.dec_kepler),
            CIRCLE('ICRS', gaiadr3.gaia_source.ra, gaiadr3.gaia_source.dec, 1./3600.0)
        )
    INNER JOIN  gaiadr3.astrophysical_parameters
    ON
        gaiadr3.gaia_source.source_id = gaiadr3.astrophysical_parameters.source_id
    """

    # launch the query and get the results
    job = Gaia.launch_job(query)
    results = job.get_results()

    # convert the results to a pandas dataframe
    gaialum = results.to_pandas()

    # rename source_id to gaiadr3_id
    gaialum = gaialum.rename(columns={"source_id": "gaiadr3_id"})

    # drop duplicates
    gaialum = gaialum.drop_duplicates(subset="hostname", keep="first")

    # merge the gaialum table with the df table on hostname
    df = df.merge(gaialum, on="hostname", how="left")

    # convert lum_flame from linear to log
    df["lum_flame"] = np.log10(df["lum_flame"])

    # convert lum flape upper and lower from to error in log
    df["lum_flame_upper"] = np.log10(df["lum_flame_upper"] - df["lum_flame"])
    df["lum_flame_lower"] = -np.log10(df["lum_flame"] - df["lum_flame_lower"])

    # define reflink for the values you want to replace
    reflink = ("<a refstr=FOUESNEAU__ET_AL__2022 "
               "href=https://ui.adsabs.harvard.edu/abs/2022arXiv220605992F "
               "target=ref>Fousneau et al. 2022</a>")

    # pick where no error on luminosity
    no_err_st_lum = df["st_lumerr1"].isna()

    # replace with FLAME luminosity and errors
    df.loc[no_err_st_lum, "st_lum"] = df.loc[no_err_st_lum, "lum_flame"]
    df.loc[no_err_st_lum, "st_lumerr1"] = df.loc[no_err_st_lum, "lum_flame_upper"]
    df.loc[no_err_st_lum, "st_lumerr2"] = df.loc[no_err_st_lum, "lum_flame_lower"]

    # replace the reflink in st_lum_reflink
    df.loc[no_err_st_lum, "st_lum_reflink"] = reflink
 
    # assert table is 41 entries long
    assert len(df) == 41, len(df)

    # write back to file
    print(f"[DOWN] Writing the supplemented table to {path}")
    df.to_csv(path, index=False)



