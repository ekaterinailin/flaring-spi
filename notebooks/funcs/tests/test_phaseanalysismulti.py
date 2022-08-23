import os

import pytest
import pandas as pd
import numpy as np

from astropy.io import fits
from astropy.table import Table

from ..phaseanalysismulti import (get_cumulative_distribution,
                                  get_observed_phases,
                                  get_null_hypothesis_distribution,
                                  )


def _make_lightcurve(length, duration, period, phaseshift=0., 
                     mission="Kepler", ID="123",TIC="456",qcs=1, i=0,
                     tstamp = "2019-01-01", save=False, location="./"):
    """TEST HELPER FUNCTION. Make synthetic light curve and save to file.

    Parameters:
    -----------
    length : int
        length of the light curve
    duration : float
        duration of the light curve in some unit
    period : float
        orbital period of the target in some unit
    phaseshift : float
        shift the phases by this amount if you want
    save : bool
        save the light curve to file
    path : str
        path to the folder where the light curve is saved

    Returns:
    --------
    path to the file containing the light curve
    """
    # make a light curve phases
    time = np.linspace(0, duration, length)
    print(time[-3:], time[:3])
    phases = (time % period) / period

    # make a light curve with a phase shift
    phases_shift = (phases + phaseshift) % 1
   
    # create pd.Series with information about the LC
    lc = pd.Series({"mission": mission, 
                    "ID": ID, 
                    "TIC": TIC,
                    "quarter_or_sector": qcs,
                    "timestamp": tstamp, 
                    "duration": duration})

    if mission == "Kepler":
        hostname = ID
    elif mission == "TESS":
        hostname = TIC

    # create a path to the file using the information about the LC
    path = f"{location}{tstamp}_{hostname}_{qcs:02}_altai_{i}.fits"

    # save the light curve to a fits file
    if save:
        # make a fits file
        hdu = fits.PrimaryHDU(np.array([1.,2.,3.]))
        data = fits.BinTableHDU(data=Table.from_pandas(pd.DataFrame({"time":time, 
                                   "phase":phases_shift})))
        hdulist = fits.HDUList([hdu, data])
        hdulist.writeto(path, overwrite=True)

    return path, lc
    

def test_get_cumulative_distribution():
    """Test the get_cumulative_distribution function.
    
    
    """
    # -------------------------------------------------------------------------
    # PREPARE DATA
    # 1. Define table with LCs
    timestamps = ["2022_08_01","2022_08_02","2022_08_01"]
    TICs = [1,2,3]
    IDs = ["a","b","c"]
    missions = ["TESS","TESS","Kepler"]
    qcs = ["1","2","3"]
    lcs = pd.DataFrame({"timestamp": timestamps,
                        "TIC": TICs,
                        "ID": IDs,
                        "mission":missions,
                        "quarter_or_sector":qcs})

    # 2. Define table with observed phases
    
    # Define column names by combining qcs and missions
    colnames = [f"{mission}_{qc}" for qc, mission in zip(qcs,missions)]

    # Define table with observed phases
    observedphases = pd.DataFrame(np.random.rand(300,3), columns=colnames)

    # 3. Define table with flares
    timestamps = ["2022_08_01","2022_08_02","2022_08_01",
                  "2022_08_01","2022_08_02","2022_08_01"]
    TICs = [1,2,3,1,2,3]
    IDs = ["a","b","c","a","b","c"]
    missions = ["TESS","TESS","Kepler","TESS","TESS","Kepler"]
    qcs = ["1","2","3","1","2","3"]
    orbital_phase = [0.1,0.2,-1,0.4,0.5,-1]
    flares = pd.DataFrame({"timestamp": timestamps,
                            "TIC": TICs,
                            "ID": IDs, 
                            "mission":missions,
                            "quarter_or_sector":qcs,
                            "orbital_phase":orbital_phase})



    # -------------------------------------------------------------------------
    # TEST FUNCTION

    # Test 1: full table

    # Calculate cumulative distribution function
    n_exp, cum_n_exp = get_cumulative_distribution(flares, observedphases, lcs)

    # Assert that the cumulative distribution function is correct
    assert len(cum_n_exp) == len(n_exp) + 1
    
    # Assert that the length of the cumulative distribution function is
    # the same as the length of the observed phases table
    assert len(cum_n_exp) == 301

    # Assert that the maximum value of the cumulative distribution function
    # is 1, and the minimum 0
    assert cum_n_exp.max() == 1
    assert cum_n_exp.min() == 0

    # 2 flare per 150 days in each of 3 LCs, and 0.5 days in each bin
    assert np.mean(n_exp) == pytest.approx(2 / 150. * 3 *.5) 

    
    # Test 2: partial table

    # Define partial table with only two flares
    flares_partial = flares.iloc[:2]

    # Calculate cumulative distribution function
    n_exp, cum_n_exp = get_cumulative_distribution(flares_partial, observedphases, lcs)

    # Assert that the cumulative distribution function is correct
    assert len(cum_n_exp) == len(n_exp) + 1

    # Assert that the length of the cumulative distribution function is
    # the same as the length of the observed phases table
    assert len(cum_n_exp) == 301

    # Assert that the maximum value of the cumulative distribution function
    # is 1, and the minimum 0
    assert cum_n_exp.max() == 1
    assert cum_n_exp.min() == 0

    # 2 flare per 150 days in one of 3 LCs, and 0.5 days in each bin
    assert np.mean(n_exp) == pytest.approx(2 / 150. * 1 *.5) 



    # -------------------------------------------------------------------------
    # TEST FAILURE

    # Test 1: phases are not within the 0-1 range
    with pytest.raises(AssertionError) as e:
        observedphases_alt = pd.DataFrame(np.random.rand(300,3)*2., columns=colnames)
        get_cumulative_distribution(flares, observedphases_alt, lcs)
        assert "One or more observedphases" in str(e.value)

    # Test 2: LC is missing in the table of observed phases
    with pytest.raises(ValueError) as e:
        lcs["mission"] = ["TESS"] * 3
        get_cumulative_distribution(flares, observedphases, lcs)
        assert "column not found in observedphases" in str(e.value)


def test_get_observed_phases():
    """Test the get_observed_phases function. Integration test only so far"""

    # create three light curves
    path1, lc1 = _make_lightcurve(length=10000, duration=10., period=1.,
                                phaseshift=0, save=True, location="./")
    path2, lc2 = _make_lightcurve(length=20000, duration=10., period=1., i=1,
                                phaseshift=0., save=True, location="./")
    path3, lc3 = _make_lightcurve(length=20000, duration=10., period=1., i=0,
                                mission="TESS",
                                phaseshift=0., save=True, location="./")

    # concatenate lc1 and lc2 to a DataFrame
    lcs = pd.concat([lc1, lc2, lc3], axis=1).T

    # create an array of flare phases
    phases = [0.5, 0.6, 0.7, .8]

    # get the observed phases for the flares
    observedphases, binmids = get_observed_phases(phases, lcs, "./")

    # check that the observed phases are correct
    assert observedphases.iloc[0,:].values == pytest.approx(5., rel=1e-3)
    assert observedphases.iloc[1:-1,:].values == pytest.approx(1., rel=1e-3)
    assert observedphases.iloc[-1,:].values == pytest.approx(2., rel=1e-3)

    # total observing time should be correct
    assert (observedphases.sum(axis=0).values == pytest.approx(10.,rel=1e-2))

    # assert all columns are correctly read in
    assert (observedphases.columns.values == 
            ["Kepler_1_0", "Kepler_1_1", "TESS_1_0"]).all()

    # assert the number of bins is correct
    assert observedphases.shape[0] == len(phases) + 1

    # delete the file created
    for p in [path1, path2, path3]:
        os.remove(p)


def test_get_null_hypothesis_distribution():
    """Test the get_null_hypothesis_distribution function."""

    p = [.1, .5, .6, .9]
    cum_n_exp = [0.,.1, .3, .6, .9, 1.]

    f = get_null_hypothesis_distribution(p, cum_n_exp)

    # make sure the interpolation function is hitting the endpoints
    assert f(0.) == 0.
    assert f(1.) == 1.

    # check a random value inbetween the endpoints
    assert f(.5) == .3

    #plt.plot(p, cum_n_exp[1:-1], "o")
    #plt.plot(p, f(p))
