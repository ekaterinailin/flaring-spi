# Flaring Star-Planet Interactions in Kepler and TESS Data

This repository contains the statistical analysis modules, and scripts to produce the figures and tables in the two following works in progress:

1. Ilin, E. and Poppenhäger K., (2022, in prep.). *Searching for flaring star-planet interactions in AU Mic TESS observations.*
2. Ilin, E., Poppenhäger K. et al. (2022, in prep.) *Searching for flaring star-planet interactions in Kepler and TESS observations.*

#### How to navigate this repository

All file names than contain `AU_Mic` belong to project 1. The other files mostly belong to project 2., but can also be modules or scripts that are used in both.

- `notebooks/` contains all the notebooks and scripts as well as all the modules stored in `funcs/`
- `data/` contains ancillary data, such as the TESS transmission curve
- `results/` contain the flare table, the results of the A-D tests

## Project 1: Searching for flaring star-planet interactions in AU Mic TESS observations

Notebooks that produce the figures and tables in the paper found in `notebooks/`:

- `AU_Mic_flare_catalog_with_phases_TABLE1.ipynb`
  - electronic table of all flares: `results/2021_12_AU_Mic_2021_12_AU_Mic_final_flare_table.csv`
- `AU_Mic_AD_test_analysis_TABLE2.ipynb`
  - all A-D tests: `results/adtests.csv` (time stamped within the table) 
- `AU_Mic_illustrate_flares_FIGURE1.ipynb`
- `AU_Mic_FFD_FIGURE2.ipynb`
- `AU_Mic_cumdist_plot_FIGURE3.ipynb`
- `AU_Mic_illustrate_observability_FIGURE4.ipynb`

Look into A-D tests with doubled and tripled sample: `AU_Mic_AD_test_doubled_sample_DISCUSSION.ipynb`

Calculation of quiescent flux in the TESS band: `AU_Mic_quiescent_flux_in_TESS_band_RESULTS.ipynb`



## Required Python packages

- numpy
- scipy
- astropy
- matplotlib
- lightkurve
- altaipony
- transitleastsquares
- emcee
- corner`
