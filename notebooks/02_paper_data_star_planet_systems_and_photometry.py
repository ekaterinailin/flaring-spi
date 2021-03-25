"""
UTF-8, Python 3

------------------
Flaring SPI
------------------

Ekaterina Ilin, 2021, MIT License

This script generates the numbers quoted in the paper
in the Data section on Star-Planet Systems and Kepler
and TESS photometry.
"""


from funcs.notebook import *


if __name__ == "__main__":

    path = "/home/ekaterina/Documents/002_writing/flaring-spi-draft/flaring-spi-draft/values/"

    # Get all sample tables

    nasa = pd.read_csv("../data/20_01_2021_confirmed_uncontroversial_exoplanet_systems.csv")
    tess = pd.read_csv("../data/20_01_2021_tess_toi_candidates_known_planets.csv")

    # All uncontroversial systems in NASA Composite Table

    with open(f"{path}nasa_all_systems.txt","w") as f:
        f.write(str(nasa.shape[0]))

    mprint(f"All uncontroversial systems in NASA Composite Table: {nasa.shape[0]}")

    # -----------------------------------------------------------

    # All unique TESS TOIs

    _ = tess.TIC.unique().shape[0]

    with open(f"{path}tess_unique_tois.txt","w") as f:
        f.write(str(_))

    mprint(f"All unique TESS TOIs: {_}")

    # -----------------------------------------------------------

    mprint("Now look at the light curves we found:")

    # -----------------------------------------------------------

    # Get light curve tables:

    lcs = pd.read_csv("../data/20_01_2021_full_kepler_k2_tess_exoplanet_lcs_complete.csv")
    lcsselect = pd.read_csv("../data/20_01_2021_full_kepler_k2_tess_exoplanet_lcs_some_excluded.csv") #Kepler-451 excluded, dubbed 2MASS J19383260+4603591


    # Light curves of systems in both missions

    for mission in ["Kepler", "TESS"]:
        n = lcs.groupby("mission").mission.count()[mission]
        # light curves of systems
        with open(f"{path}lcs_{mission}.txt","w") as f:
            f.write(str(n))
            mprint(f"Light curves from the {mission} mission: {n}")

    # -----------------------------------------------------------

    # IGNORE K2!

    lcs = lcs[lcs.mission != "K2"]

    # -----------------------------------------------------------

    # How many systems were observed with both missions?
    in_kepler_and_tess = (np.where(lcs.groupby(["ID","mission"]) # select ID mission pairs
                                   .first() #remove all but one light curve from each mission
                                   .reset_index() 
                                   .groupby("ID") #if there is more than one entry for an ID 
                                   .count() # it's because both Kepler and TESS observed it
                                   .mission.values == 2)[0].shape[0]) # how many were observed by both missions?

    with open(f"{path}systems_observed_by_both_kepler_and_tess.txt","w") as f:
        f.write(str(in_kepler_and_tess))


    mprint(f"Systems observed with both missions: {in_kepler_and_tess}")

    # -----------------------------------------------------------

    # How many systems were observed either TESS or Kepler, but not both?

    # get all ID+mission pairs
    _ = lcs.groupby(["ID","mission"]).first().reset_index()

    # select the Kepler subsample 
    # and subtract systems that were also observed by the other mission
    observed_by_kepler_only = _[_.mission == "Kepler"].shape[0] - in_kepler_and_tess//2
    observed_by_tess_only = _[_.mission == "TESS"].shape[0] - in_kepler_and_tess//2

    # sanity check
    total_systems_observed = lcs.groupby("ID").first().shape[0]
    assert observed_by_kepler_only + observed_by_tess_only == total_systems_observed

    with open(f"{path}systems_observed_by_kepler_only.txt","w") as f:
        f.write(str(observed_by_kepler_only))

    mprint(f"Systems observed by Kepler only: {observed_by_kepler_only}")

    with open(f"{path}systems_observed_by_tess_only.txt","w") as f:
        f.write(str(observed_by_tess_only))

    mprint(f"Systems observed by TESS only: {observed_by_tess_only}")
   
    mprint(f"Systems observed in total by Kepler and TESS: {total_systems_observed}")

    # -----------------------------------------------------------