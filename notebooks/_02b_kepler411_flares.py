import pandas as pd
from astropy.io import fits

if __name__ == "__main__":

    # read in Jackman et al. data
    path = '../data/jackman2021.fits'
    hdu = fits.open(path)
    df = pd.DataFrame(hdu[1].data)
    print(f"[UP] Jackman et al. 2021 flare table fits file from {path}")

    # big endian to little endian in each column
    df = df.apply(lambda x: x.values.byteswap().newbyteorder())

    # pick Kepler-411
    kep411 = df[df["Source ID"] == "Gaia DR2 2132768952604988672"]

    # rename columns Flare Start, Flare End, Detected Flare Amplitude to
    # tstart, tstop, ampl_rec
    kep411 = kep411.rename(columns={"Flare Start": "tstart", 
                                    "Flare End": "tstop", 
                                    "Detected Flare Amplitude": "ampl_rec"})

    # read in flare table from results
    path = '../results/2022_07_flares_vetted.csv'
    flare_table = pd.read_csv(path)
    print(f"[UP] Our flare table file from {path}")

    # rename rows where ID is Kepler-411 to Kepler-411(c) like contaminated
    flare_table.loc[flare_table["ID"] == "Kepler-411", "ID"] = "Kepler-411(c)"

    # rename rows where TIC is 399954349 to 399954349(c) like contaminated
    flare_table.loc[flare_table["TIC"] == 399954349, "TIC"] = "399954349(c)"

    # pick the columns from kep411 from above plus add the TIC, ID, qcs, mission, and tstamp
    qcs = 13 #
    mission = 'Kepler'
    tstamp = '2023_02_13'
    TIC = 399954349
    note = "Jackman2021"
    ID = "Kepler-411"

    # add the columns to the kep411 table
    kep411['qcs'] = qcs
    # except for the 1366.59396 flare, which is a qcs 14 one
    kep411.loc[kep411["tstart"] > 1366.5939, "qcs"] = 14
    kep411['mission'] = mission
    kep411['tstamp'] = tstamp
    kep411['TIC'] = TIC
    kep411['ID'] = ID
    kep411['note'] = note
    kep411["real"] = 1
    kep411["helpid"] = f"{ID}_{mission}_{qcs}"


    # read in the systems parameters from data
    path = "../results/params_of_star_planet_systems_with_AD_tests.csv"
    params = pd.read_csv(path)
    print(f"[UP] System parameters table from {path}")

    # get transit midtime
    bjdtt = params[params.ID == "Kepler-411", "pl_tranmid_kepler"].iloc[0]

    # get orbital period
    orbper = params[params.ID == "Kepler-411", "pl_orbper_kepler"].iloc[0]


    # calculate the flare phases
    kep411["phase"] = (kep411["tstart"]  + 2454833. - bjdtt) / orbper % 1.

    # old number of rows
    oldsize = flare_table.shape[0]

    # pick these columns, and the above columns, and add them to flare_table
    flare_table = flare_table.append(kep411[['tstart', 'tstop', 'ampl_rec', 
                                            'qcs', 'mission', 'tstamp', 'TIC', 
                                            'ID', 'note', 'real', 'helpid',
                                            'phase']])

    # assert that the number of rows has increased by 7
    assert flare_table.shape[0] == oldsize + 7

    # assert that the old part of the table is unaffected
    assert flare_table.iloc[:oldsize].equals(flare_table.iloc[:oldsize])

    # write the flare table to a csv
    flare_table.to_csv('../results/2022_07_flares_vetted.csv', index=False)
    print(f"[DOWN] Our flare table file with new Kepler-411 flares to {path}")