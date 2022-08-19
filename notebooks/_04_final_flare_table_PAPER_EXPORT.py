import pandas as pd

if __name__ == "__main__":

    # Read the vetted flare table
    flares = pd.read_csv("../results/2022_07_flares_vetted.csv")

    # Select only the confirmed flares
    final_table = flares[(flares.real == 1)]

    # Add lines with LCs without flares for star with flares found in 
    # other LCs
    selection_criteria = ((flares.real == -1) & # no flares found in this LC
                          (flares.TIC.isin(final_table.TIC))) # flares found in other LCs
    final_table = final_table.append(flares[selection_criteria])

    # For the final csv table, sort by TIC ID
    final_table = final_table.sort_values(by="TIC")

    # Select columns to keep
    cols = ['TIC', 'ID', 'mission', 'qcs', 'tstart', 'tstop',
            'ampl_rec','ed_rec','ed_rec_err','phase','tstamp']
    final_table = final_table[cols]

    # Rename columns to be more informative
    final_table = final_table.rename(index=str,
                                    columns={"qcs": "quarter_or_sector",
                                            "phase": "orbital_phase",
                                            "ampl_rec": "rel_amplitude",
                                            "ed_rec": "ED",
                                            "ed_rec_err": "ED_err",
                                            "tstamp": "timestamp",
                                            })

    # Write table to the results folder in the analysis repo
    final_table.to_csv("../results/PAPER_flare_table.csv",
                    index=False, header=True)
                    
    # Write table to the input folder in the paper repo
    paperpath = "../../../002_writing/flaring-spi-paper/src/data/PAPER_flare_table.csv"
    final_table.to_csv(paperpath, index=False, header=True)