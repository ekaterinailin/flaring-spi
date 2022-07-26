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

1. Preselect TOI table: exclude FPs, binaries, etc.
1b. Search the preselected table for LCs.
2. Preselect NASA table: uncontroversial confirmed systems with transits.
2b. Search the preselected table for LCs.
3. Match TIC IDs to star identifiers into 1b. result
4. Merge TESS entries from 3. and 2b.
5. Add Kepler/K2 entries from 3 to 4.
6. Add Kepler/K2 entries from 2b. to 5.
7. Some diagnostic output and file saving
"""

from funcs.notebook import *
from lightkurve import search_lightcurvefile
from lightkurve.search import _query_mast

import time

def get_TIC(ID):
    """Get TIC ID by searching MAST and
    selecting TESS mission ID.
    
    Parameter:
    ----------
    ID : str
        stellar identifier
        
    Return:
    --------
    TIC ID
    """
    f = _query_mast(ID).to_pandas()
    res = f[(f["dataproduct_type"]=="timeseries") & (f["obs_collection"]=="TESS")].target_name.iloc[0]
    return res

    
def preselect_toi_catalog(tstamp):

    # read in TESS-TOI sample
    path = "../data/2022_07_26_TESS_TOI_CATALOG.csv"
    mprint(f"[UP] Using TESS-TOI Table from {path}")
    df = pd.read_csv(path, skiprows=90)

    # remove EB, V, O flagged entries
    df1 = df[df["tfopwg_disp"].isin(["KP","CP"])]
     
    # one entry per system
    dfs = df1.groupby("tid").first().reset_index()
    dfs = dfs.rename(index=str, columns={"tid":"TIC"})

    # save to file
    path = f"../data/{tstamp}_tess_toi_candidates_known_planets.csv"
    mprint(f"[DOWN] Saving CP and KP sample of {dfs.shape[0]} planet hosts from TESS TOI to {path}")
    dfs.to_csv(path, index=False)
    
    return dfs  


def preselect_nasa_catalog(tstamp):

    # Composite Table of confirmed exoplanets
    path = "../data/2022_07_26_NASA_COMPOSITE.csv"
    mprint(f"[UP] Using NASA Composite Table from {path}")
    df = pd.read_csv(path, skiprows=320) # composite table
    # df.columns.values

    # select only uncontroversial detections
    sel = df[(df.pl_controv_flag==0) & (df.tran_flag==1)]
#     sel.shape, sel.groupby("hostname").first().shape

    dfs = sel.groupby("hostname").first()
    path =f"../data/{tstamp}_confirmed_uncontroversial_transiting_exoplanet_systems.csv"
    mprint(f"[DOWN] Saving {dfs.shape[0]} uncontroversial transiting exosystems from NASA Composite Table to {path}")
    dfs.to_csv(path, index=False)
    
    return dfs


if __name__ == "__main__":

    tstamp = time.strftime("%Y_%m_%d", time.localtime())

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 1. SELECTION in TESS-TOI CATALOG
    # -----------------------------------------------------------------------

    dfs = preselect_toi_catalog(tstamp)

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 1b. MAST SEARCH in TESS-TOI CATALOG
    # -----------------------------------------------------------------------

#     mprint("Beginning search for SC LCs in TOI catalog.")
#     x = 0
#     N = dfs.shape[0]
#     # save first row:
#     with open(f"../data/{tstamp}_toi_gotlc.csv", "w") as f:
#                f.write("sector,mission,TIC\n")
#     for i, row in dfs.iterrows():
#         TIC = row.TIC
#         x += 1
#         try:
#             lst = search_lightcurvefile(f"TIC {TIC}", cadence="short")
#             print("TIC")
#             print(lst.table.to_pandas()[["mission","project"]])
#             print(lst.table.to_pandas().head().T)
#             _ = {"sector":lst.table.to_pandas()["mission"].str[-2:].values,
#                  "mission":lst.table.to_pandas()["project"].values,
#                  "TIC":TIC}
#             r = pd.DataFrame(_)
#             with open(f"../data/{tstamp}_toi_gotlc.csv", "a") as f:
#                 r.to_csv(f,header=False, index=False)
#             print(f"[{x / N * 100:.0f} %]", "TIC ", i, r.shape[0], x)
#         except:
#             with open(f"../data/{tstamp}_toi_nolc.txt", "a") as f:
#                 f.write(f"TIC {TIC}\n")
#     mprint("Finished search for SC LCs in TOI catalog.")


    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 2. SELECTION in NASA CATALOG
    # -----------------------------------------------------------------------

    # # search for light curves
    # # last update: 2022/07/26

    # # -----------------------------------------------------------------------
    # # read in confirmed and unconroversial NASA exoplanets and find LCs
    # # -----------------------------------------------------------------------

    dfs = preselect_nasa_catalog(tstamp)

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 2b. MAST SEARCH in NASA CATALOG
    # -----------------------------------------------------------------------

    mprint("Begin search for SC LCs in NASA catalog.")
    x = 0
    N = dfs.shape[0]
    print(N)
    # save first row:
    with open(f"../data/{tstamp}_nasa_gotlc.csv", "w") as f:
               f.write("qcs,mission,ID\n")
    for i, row in dfs.iterrows()[1000:]:
        hostname = row.name
        x += 1
        try:
            lst = search_lightcurvefile(hostname, cadence="short")
            print(hostname)
            print(lst.table.to_pandas()[["mission","project"]])
            _ = {"qcs":lst.table.to_pandas()["mission"].str[-2:].values,
                 "mission":lst.table.to_pandas()["project"].values,
                 "ID":hostname}
            r = pd.DataFrame(_)
            mprint(r)
            with open(f"../data/{tstamp}_nasa_gotlc.csv", "a") as f:
                r.to_csv(f,header=False, index=False)
            print(f"[{x / N * 100:.0f} %]", hostname, r.shape[0], x)
        except:
            with open(f"../data/{tstamp}_nasa_nolc.txt", "a") as f:
                f.write(f"{hostname}\n")
            print("exception")
    mprint("Finished search for SC LCs in NASA catalog.")

    # # -----------------------------------------------------------------------
    # # last updated 2022/07/26

    # # To find the overlap between the NASA and the TOI catalogs 
    # # we need the TIC IDs for all stars in the NASA catalog 
    # # that were also found in TESS

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 3. Match TIC IDs to star identifiers in 1b. result
    # -----------------------------------------------------------------------

    # read in results from the NASA LC search
    path = f"../data/{tstamp}_nasa_gotlc.csv"
    mprint(f"[UP] Using SC LC matches in NASA catalog from {path}")
    dd = pd.read_csv(path, names=["qcs","mission","ID"])

    # select the light curves found in the TESS database
    dd = dd.loc[dd.mission=="TESS"]

    # get only the system names
    dd = dd.drop_duplicates(subset="ID")

    # find TIC ID for every star in NASA table that has TESS light curves
    mprint("Beginning NASA ID to TESS TIC matching via MAST.")
    dd["TIC"] = dd.apply(lambda x: get_TIC(x.ID), axis=1)
    mprint("Finished NASA ID to TESS TIC matching via MAST.")

    # save to file
    path = f"../data/{tstamp}_tic_name_match.csv"
    mprint(f"[DOWN] NASA ID to TESS TIC matches to {path}")
    dd.to_csv(path, index=False)

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 4. MERGE TESS ENTRIES from NASA and TOI CATALOGS
    # -----------------------------------------------------------------------

    # read in TIC ID matches
    path = f"../data/{tstamp}_tic_name_match.csv"
    mprint(f"[UP] Using TIC ID matches from {path}")
    dd = pd.read_csv(path)

    dd["TIC"] = dd["TIC"].astype(int)

    # read in results from the NASA LC search
    path = f"../data/{tstamp}_nasa_gotlc.csv"
    mprint(f"[UP] Using NASA exoplanet LCs table from {path}")
    nasa = pd.read_csv(path)

    # merge TIC IDs into the main LC list
    nasa2 = nasa.merge(dd[["TIC","ID"]], how="left", 
                       on=["ID"]).drop_duplicates(keep="first")

    # set NASA catalog flag
    nasa2["catalog_nasa"] = 1

    # select NASA catalog entries that appear in TESS
    nasa2_tic = nasa2.loc[nasa2.mission == "TESS",
                          ["qcs","TIC","ID","catalog_nasa"]]

    # read in TOI catalog LCs
    path = f"../data/{tstamp}_toi_gotlc.csv"
    mprint(f"[UP] Using TESS-TOI exoplanet LCs table from {path}")
    toi = pd.read_csv(path)

    # set TOI catalog ID
    toi["catalog_toi"] = 1

    # pick TESS mission entries from TOI catalog
    toi_tic = toi.loc[toi.mission=="TESS", 
                      ["qcs","TIC", "catalog_toi"]]
    toi_tic = toi_tic.astype(int).drop_duplicates(keep="first")

    # find out what TESS mission entries are found in both catalogs
    inboth = pd.merge(nasa2_tic, toi_tic, how ='inner',
                      on =['qcs', 'TIC']).shape[0]

    # merge TESS mission entries from both catalogs
    together = pd.merge(toi_tic, nasa2_tic, how ='outer', on =['qcs', 'TIC'])

    # check that the merging went right
    first, second = toi_tic.shape[0], nasa2_tic.shape[0]
    assert together.shape[0] == first + second - inboth

    # fill catalog flag NaNs will zeros
    together[["catalog_nasa", "catalog_toi"]] = together[["catalog_nasa", "catalog_toi"]].fillna(0.)

    # fill missing IDs with empty strings
    together["ID"] = together["ID"].fillna("")

    # the mission is TESS for all entries in this table 
    # anyways, but we dropped the column before for convenience, 
    # so let's re-introduce it
    together["mission"] = "TESS"

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 5. ADD KEPLER/K2 ENTRIES from NASA catalog
    # -----------------------------------------------------------------------

    # select NASA LCs table entries from Kepler and K2 missions
    keplerk2 = nasa2[(nasa2.mission=="Kepl") | (nasa2.mission=="K2 C")]

    # fill in missing TIC IDs with zeros
    keplerk2["TIC"] = keplerk2["TIC"].fillna(0).astype(int) 


    # add NASA LCs table entries from Kepler and K2 missions 
    # to the full TESS mission table
    tesskeplerk2 = pd.merge(keplerk2, together, how="outer")

    # there should be no overlap, so check that:
    tesskeplerk2.shape[0], keplerk2.shape[0] + together.shape[0]

    # set the TOI catalog flag to zero in the NASA table Kepler/K2 entries
    # that we just added:
    tesskeplerk2.loc[tesskeplerk2.catalog_toi.isnull(),"catalog_toi"] = 0

    # rename the mission flags
    tesskeplerk2.loc[tesskeplerk2.mission=="Kepl", "mission"] = "Kepler"
    tesskeplerk2.loc[tesskeplerk2.mission=="K2 C", "mission"] = "K2"


    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 6. ADD KEPLER/K2 ENTRIES FROM TOI CATALOG
    # -----------------------------------------------------------------------

    # Now we have treated all TESS entries from the TOI catalog,
    # all entries from the NASA catalog, but not yet
    # the Kepler/K2 entries from the TOI catalog
    # (Note: Are these the planets that were detected or confirmed in TESS
    # but had a Kepler or K2 light curve already?)

    # pick Kepler/K2 mission entries from TOI catalog
    keplerk2extra = toi.loc[toi.mission!="TESS",
                            ["qcs", "TIC", "catalog_toi", "mission"]]
    keplerk2extra = keplerk2extra.drop_duplicates(keep="first")

    # rename mission entries 
    keplerk2extra.loc[keplerk2extra.mission=="Kepl", "mission"] = "Kepler"
    keplerk2extra.loc[keplerk2extra.mission=="K2 C", "mission"] = "K2"

    # merge the entries Kepler/K2/TESS-NASA+TESS-TOI catalog onto
    # the Kepler/K2-TOI entries
    notfoundinnasa = pd.merge(keplerk2extra, tesskeplerk2, how="left", 
                              suffixes=("","_tkk"), on=["mission","qcs","TIC"])

    # In the merged table pick the entries come from the Kepler/K2-TOI catalog ALONE
    notfoundinnasa = notfoundinnasa[notfoundinnasa.catalog_nasa.isnull()]

    # Find the overlap between 
    # Kepler/K2/TESS-NASA+TESS-TOI catalog and Kepler/K2-TOI catalog
    alreadyfoundinnasa = pd.merge(keplerk2extra, tesskeplerk2, 
                                  how="inner", 
                                  on = ["mission","qcs","TIC"]).shape[0]

    # already found in Kepler/K2/TESS-NASA+TESS-TOI and not found in it
    # should add up to the full Kepler/K2-TOI list
    assert notfoundinnasa.shape[0] + alreadyfoundinnasa, keplerk2extra.shape[0]


    # Add the entries that are new from the  Kepler/K2-TOI catalog to get the
    # Kepler/K2/TESS-TOI + Kepler/K2/TESS-NASA full catalog 
    full = pd.concat([tesskeplerk2, notfoundinnasa[["qcs","TIC","catalog_toi",
                                                    "ID","catalog_nasa","mission"]]], 
                     ignore_index=True)

    # fill missing values
    full["catalog_nasa"] = full["catalog_nasa"].fillna(0)
    full["ID"] = full["ID"].fillna("")

    # sanity check that concat worked as expected
    assert notfoundinnasa.shape[0] + tesskeplerk2.shape[0] == full.shape[0]

    # fill in ID column with TIC IDs if no other ID is available
    full.loc[full.ID=="","ID"] = full.loc[full.ID=="","TIC"].apply(lambda x: f"TIC {x}")


    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # 7. GET SOME DIAGNOSTIC NUMBERS FROM THE FINAL SAMPLE and STDOUT them
    # -----------------------------------------------------------------------

    # STDOUT
    mprint("Number of light curves from each mission:")
    mprint(full.groupby("mission").mission.count())

    # Get an estimate of expeted flare number
    m = {"K2":80.,"Kepler":90.,"TESS":25.}
    estimate_obs_time = lambda x: m[x]
    full["est_obs_time_days"] = full.mission.apply(estimate_obs_time)
    minflares = full.est_obs_time_days.sum() / 365 * 0.01 * 12

    # STDOUT
    mprint(f"Expected number of flares assuming \n"
          f"1% flaring stars detected per month of observation:\n"
          f">>>> {minflares:.0f}-{minflares*10:.0f} flares <<<<")

    # STDOUT
    path = f"../data/{tstamp}_full_kepler_k2_tess_exoplanet_lcs.csv"
    mprint(f"[DOWN] Saving full sample to {path}")
    # SAVE TO FILE
    full.to_csv(path, index=False)

    # STDOUT
    nstars = full.groupby("ID").first().shape[0]
    mprint(f"{nstars} individual systems in sample with short cadence LC.")
    
    # STDOUT
    path = f"../data/{tstamp}_nasa_nolc.txt"
    mprint(f"[UP] Using IDs of targets without SC LCs in NASA database from {path}")
    d = pd.read_csv(path, names=["ID"])
    nosclcnasa = d.drop_duplicates().shape[0]
    mprint(f"{nosclcnasa} systems in the NASA database have no short cadence LC.")
    butkepler = d[d.ID.str.contains("Kepler")].shape[0]
    mprint(f"Out of these, {butkepler} systems have a Kepler ID nonetheless.")
    
    # STDOUT
    path = f"../data/{tstamp}_toi_nolc.txt"
    mprint(f"[UP] Using IDs of targets without SC LCs in TOI database from {path}")
    d = pd.read_csv(path, names=["ID"])
    nosclctoi = d.drop_duplicates().shape[0]
    mprint(f"{nosclctoi} systems in the TOI database have no short cadence LC.") 
    
