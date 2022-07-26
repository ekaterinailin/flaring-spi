# Flaring Star-Planet Interactions in Kepler and TESS Data

This repository contains the statistical analysis modules, and scripts to produce the figures and tables in the two following works in progress:

1. [Ilin, E. and Poppenhäger K., (2022). *Searching for flaring star-planet interactions in AU Mic TESS observations.*](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.4579I/abstract)
2. Ilin, E., Poppenhäger K. et al. (2023, in prep.) *Searching for flaring star-planet interactions in Kepler and TESS observations.*

#### How to navigate this repository

**Nomenclature:**

- File names that contain `AU_Mic` or `AUMic` belong to project 1. 
- File names that contain `ALL` belong to project 2.
- File names with neither belong to both, like `notebooks/findflares.py`.
- File names with a `_` prefix are work in progress.

**Folders:**

- `notebooks/` contains all the notebooks and scripts 
  - also contains the modules used in both, stored in `funcs/`
- `data/` contains ancillary data, such as the TESS transmission curve
- `results/` contain the flare table, the results of the A-D tests

## Project 1: Searching for flaring star-planet interactions in AU Mic TESS observations

Notebooks that produce the figures and tables in the paper found in `notebooks/`:

- `AU_Mic_flare_catalog_with_phases_TABLE1.ipynb`
  - example script of how to use the de-trending and flare finding described in Methods: `notebooks/findflares.py`
  - flaring finding method described in 2.1: `notebooks/funs/detrend.py`
  - electronic table of all flares: `results/2021_12_AU_Mic_2021_12_AU_Mic_final_flare_table.csv`
- `AU_Mic_AD_test_analysis_TABLE2.ipynb`
  - all A-D tests: `results/adtests.csv` (time stamped within the table) 
- `AU_Mic_illustrate_flares_FIGURE1.ipynb`
- `AU_Mic_FFD_FIGURE2.ipynb`
- `AU_Mic_cumdist_plot_FIGURE3.ipynb`
- `AU_Mic_illustrate_observability_FIGURE4.ipynb`

Look into A-D tests with doubled and tripled sample: `AU_Mic_AD_test_doubled_sample_DISCUSSION.ipynb`

Calculation of quiescent flux in the TESS band: `AU_Mic_quiescent_flux_in_TESS_band_RESULTS.ipynb`

## Project 2: Searching for flaring star-planet interactions in Kepler and TESS observations

*work in progress*

`notebooks/`





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
