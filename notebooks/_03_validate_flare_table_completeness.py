import pandas as pd
import numpy as np

def unique_and_count(df, col="TIC"):
    unique = df[col].unique().astype("<U12")
    return unique, len(unique)


if __name__=="__main__":

    # Composite Table of confirmed exoplanets
    path = "../data/2022_07_27_input_catalog_star_planet_systems.csv"
    input_catalog = pd.read_csv(path) 
    print(f"[UP] Using compiled input catalog from {path}")

    # Table of RV systems that we added later
    path = "../data/2022_11_15_input_catalog_NONtransit_star_planet_systems.csv"
    extra_input_catalog = pd.read_csv(path) 
    print(f"[UP] Adding compiled input catalog of RV systems from {path}")

    # combine the two catalogs
    input_catalog = pd.concat([input_catalog, extra_input_catalog],
                               ignore_index=True)

    # Table of vetted flares
    path_vetted = "../results/2022_07_flares_vetted.csv"
    vetted_table = pd.read_csv(path_vetted)
    print(f"[UP] Using vetted flares table from {path_vetted}")

    # Table of unvetted flares
    path_unvetted = "../results/2022_07_flares.csv"
    unvetted_table = pd.read_csv(path_unvetted)
    print(f"[UP] Using unvetted flares table from {path_unvetted}")

    # Table of TICs without short or fast cadence LCs
    path_no_lc = "../results/2022_07_nolc.txt"
    nolc = pd.read_csv(path_no_lc)
    print(f"[UP] Using TICs without LCs from {path_no_lc}")

    # Get unique TICs in input catalog and number of TICs
    unique_input_catalog, n_unique_input_catalog = unique_and_count(input_catalog,
                                                                    col="TIC")
    # Assert number did not change for some reason
    assert n_unique_input_catalog == 2993 + 191, \
           (f"{n_unique_input_catalog} != 2993 transiting + "
           f"191 non-transiting systems")

    # Get unique TICs in flare tables and nolc table with number of TICs
    unique_vt, n_unique_in_vetted_table = unique_and_count(vetted_table)
    unique_uvt, n_unique_in_unvetted_table = unique_and_count(unvetted_table)
    unique_nolc, n_unique_nolc = unique_and_count(nolc)

    # Assert all unvetted flares ended up in the vetted table
    assert n_unique_in_vetted_table == n_unique_in_unvetted_table, \
           f"{n_unique_in_vetted_table} != {n_unique_in_unvetted_table}"

    # Assert number of unique TICs in vetted table is the same as in input catalog
    # minues the number of TICs without LCs
    assert n_unique_input_catalog == n_unique_in_vetted_table + n_unique_nolc, \
           f"{n_unique_input_catalog} != {n_unique_in_vetted_table + n_unique_nolc}"

    print(f"[FIN] Flare table completeness check finished successfully")