import os
import pytest

import pandas as pd
import numpy as np

from altaipony.flarelc import FlareLightCurve

from ..findflares import (add_meta_data_and_write,
                              write_flc_to_file,
                              write_no_lc,
                              get_table_of_light_curves,
                              get_midtime,
                              get_observed_phases,
                              get_flare_phases,
                              )

def test_get_flare_phases():
    """Test getting the orbital phases of flare candidates."""

    # Test valid input
    flare_cadenceno = np.array([0.,31.,100.])
    observed_phases = np.linspace(0,100,num=1001) % 1.
    observed_cadenceno = np.arange(1001).astype(float)

    phases = get_flare_phases(flare_cadenceno, observed_phases, 
                              observed_cadenceno)

    # phases should be between 0 and 1
    assert np.all(phases >= 0.)
    assert np.all(phases <= 1.)

    # first phase should be 0.
    assert phases[0] == 0.

    # second phase should be 0. as well
    assert phases[1] == pytest.approx(0.1, rel=1e-9)

    # third phase should be 0.3
    assert phases[2] == 0.

    # -------------------------------------------------------------------------
    # Test invalid input that should raise an error

    # Test with invalid NaN cadenceno
    flare_cadenceno = np.array([np.nan, 1., 100.])

    with pytest.raises(AssertionError) as e:
        get_flare_phases(flare_cadenceno, observed_phases, observed_cadenceno)
        assert "Input and output number of flares do not match." in str(e.value)

    # Test with out of range cadenceno
    flare_cadenceno = np.array([0., 1., 10000.])

    with pytest.raises(AssertionError) as e:
        get_flare_phases(flare_cadenceno, observed_phases, observed_cadenceno)
        assert "Input and output number of flares do not match." in str(e.value)


def test_get_observed_phases():
    """Test getting the observed phases in a light curve with a planet
    with given transit midtime and orbital period."""

    # -------------------------------------------------------------------------
    # Test valid input 
    time = np.linspace(10,20,100)
    transit_midtime = 15.
    orbital_period = 1.

    phases = get_observed_phases(time, transit_midtime, orbital_period)

    # phases should be between 0 and 1
    assert np.all(phases >= 0.)
    assert np.all(phases <= 1.)

    # first phase should be 0.
    assert phases[0] == 0.

    # -------------------------------------------------------------------------
    # Test invalid input

    # Test with invalid NaN time
    with pytest.raises(AssertionError) as e:
        time = np.array([np.nan, 1., 100.])
        get_observed_phases(time, transit_midtime, orbital_period)
        assert "Time is not finite." in str(e.value)
        
    # Test with invalid NaN transit_midtime
    with pytest.raises(AssertionError) as e:
        transit_midtime = np.nan
        get_observed_phases(time, transit_midtime, orbital_period)
        assert "Transit midtime is not finite." in str(e.value)

    # Test with invalid NaN orbital_period
    with pytest.raises(AssertionError) as e:
        orbital_period = np.nan
        get_observed_phases(time, transit_midtime, orbital_period)
        assert "Orbital period is not finite." in str(e.value)

    # Test with invalid negative orbital_period
    with pytest.raises(AssertionError) as e:
        orbital_period = -1.
        get_observed_phases(time, transit_midtime, orbital_period)
        assert "Orbital period is not positive." in str(e.value)


def test_add_meta_data_and_write():
    """
    Test the add_meta_data_and_write function.
    """
    # define default columns
    default_columns = ["istart", "istop", "cstart", "cstop", "tstart", "tstop",
                       "ed_rec", "ed_rec_err","ampl_rec", "dur", "phase"] 

    # Create a FlareLightCurve object.
    flc = FlareLightCurve(time=[1,1.1,1.2])
    flc["detrended_flux"] = [1, 2, 3]

    ff = pd.DataFrame(columns=default_columns)  # empty dataframe
    ID, TIC, sector, mission = "AB", 12345, 1, "TESS"
    lc_n, w, tstamp, mask_pos_outliers_sigma = 3, 1, "2019-01-01", 2.5

    if os.path.exists("test.csv"):
        os.remove("test.csv")
   
    print(ff)
    print(flc)

    add_meta_data_and_write(ff, flc, ID, TIC, sector, mission, lc_n, w, tstamp,
                            mask_pos_outliers_sigma, path="test.csv", 
                            header=True)

    # Check that the file was created.
    assert os.path.exists("test.csv")

    # Check that the file has the correct data.
    df = pd.read_csv("test.csv")
    assert df.shape == (1, 21)

    # Check that the file has the correct header.
    assert (df.columns.values == np.array(default_columns +
                                         ["total_n_valid_data_points", 
                                          "ID", "TIC", "qcs","mission",	"tstamp",
                                          "lc_n", "w", "mask_pos_outliers_sigma",
                                          "real"
                                          ])).all()

    # Check that the file has the correct data.
    assert df.total_n_valid_data_points.values[0] == 3
    assert df.ID.values[0] == "AB"
    assert df.TIC.values[0] == 12345
    assert df.mission.values[0] == "TESS" 
    assert df.qcs.values[0] == 1
    assert df.tstamp.values[0] == "2019-01-01"
    assert df.lc_n.values[0] == 3
    assert df.w.values[0] == 1
    assert df.mask_pos_outliers_sigma.values[0] == 2.5
    # real and phase should be -1 because the DataFrame is empty.
    assert df.real.values[0] == -1
    assert df.phase.values[0] == -1
    for col in default_columns[:-1]:
        assert np.isnan(df[col].values[0])

    # delete the file.
    os.remove("test.csv")


def test_write_flc_to_file():
    """
    Test the write_flc_to_file function.
    """
    # Create a FlareLightCurve object.
    flc = FlareLightCurve(time=[1,1.1,1.2])

    # Create a detrended FlareLightCurve object.
    dflcn = FlareLightCurve(time=[1,1.1,1.2], 
                            cadenceno=[1.,2.,3.],
                            centroid_col=[1.,2.,3.],
                            centroid_row=[1.,2.,3.],
                            mission="TESS", targetid="123")

    # Try to write the FlareLightCurve object to a file.
    # Will fail because asserts are not satisfied.
    with pytest.raises(KeyError) as e:
        
        write_flc_to_file(dflcn, flc, "test.fits")
        # check if error message is correct.
        assert "flux not in columns" in str(e.value)
    
    # Add flux column to original LC
    flc["flux"] = [1.,2.,3.]

    # Try to write the FlareLightCurve object to a file.
    # Will fail because asserts are not satisfied.
    columns = ["time", "flux", "detrended_flux", "detrended_flux_err", "it_med", 
               "phase", "flux_model"]
    for col in columns:
        with pytest.raises(KeyError) as e:
            write_flc_to_file(dflcn, flc, "test.fits")
            # check if error message is correct.
            assert f"{col} not in columns" in str(e.value) 
        # Add column
        dflcn[col] = [1.,1.,1.]
   
    # Now all necessary columns have been created, 
    # function should run without error.
    write_flc_to_file(dflcn, flc, "test.fits")

    # Check that the file was created.
    assert os.path.exists("test.fits")
    
    # delete the file.
    os.remove("test.fits")

def test_write_no_lc():
    """Test that writing to no_lc file works."""

    # create input target as pandas Series with TIC column
    target = pd.Series({"TIC": 12345})

    # write to file
    write_no_lc(target, path="test.csv")

    # check that file was created
    assert os.path.exists("test.csv")

    # read file
    df = pd.read_csv("test.csv")

    # check if file has correct data
    assert df.TIC.iloc[0] == 12345

    # delete the file.
    os.remove("test.csv")

    # create faulty input target as pandas Series without TIC column
    target = pd.Series({"not_TIC": 12345})

    with pytest.raises(AttributeError) as e:
        write_no_lc(target, path="test.csv")
        # check if error message is correct.
        assert "TIC" in str(e.value)
        # delete the file.
        os.remove("test.csv")   



def test_get_table_of_light_curves():
    """Test if light curve search handles all cases properly."""

    # -------------------------------------------------------------------------
    # Check fo for invalid input.
    input_target = pd.Series({"TIC": 12345, "hostname": "AB"})

    lc_table = get_table_of_light_curves(input_target, path="test.csv")

    # check if file was created
    assert os.path.exists("test.csv")

    # read file
    df = pd.read_csv("test.csv")

    # check if file has correct data
    assert df.TIC.iloc[0] == 12345

    # remove file
    os.remove("test.csv")

    # -------------------------------------------------------------------------
    # Check for valid input from TESS

    input_target = pd.Series({"TIC": 410214986, "hostname": "DS Tuc A"})

    lc_table = get_table_of_light_curves(input_target, path="test.csv")
    
    # check that no file was created
    assert ~os.path.exists("test.csv")

    # check if table is not empty
    assert lc_table.shape[0] > 0

    # check if table has correct columns
    assert "mission" in lc_table.columns.values
    assert "exptime" in lc_table.columns.values
   
    # check if there are no duplicate light curves as expected
    assert lc_table.shape[0] == lc_table.drop_duplicates("mission").shape[0]

    # check if only short cadence is retained
    assert (lc_table.exptime < 130.).all()

    # check that author column is neither TASOC nor K2
    assert (lc_table.author != "TASOC").all()
    assert (lc_table.author != "K2").all()

    # -------------------------------------------------------------------------
    # Check for valid input from Kepler

    input_target = pd.Series({"TIC": 377909730, "hostname": "Kepler-304"})

    lc_table = get_table_of_light_curves(input_target, path="test.csv")

    # check that no file was created
    assert ~os.path.exists("test.csv")

    # check if table is not empty
    assert lc_table.shape[0] > 0

    # check if table has correct columns
    assert "mission" in lc_table.columns.values
    assert "exptime" in lc_table.columns.values

    # check if there are no duplicate light curves as expected
    assert lc_table.shape[0] == lc_table.drop_duplicates("mission").shape[0]

    # check if only short cadence is retained
    assert (lc_table.exptime < 130.).all()

    # check that author column is neither TASOC nor K2
    assert (lc_table.author != "TASOC").all()



def test_get_midtime():
    """Test the get_midtime function."""

    # -------------------------------------------------------------------------
    # Test valid input with Kepler parameters

    # create input target as pandas Series with pl_tranmid column
    target = pd.Series({"pl_tranmid": 2454833. + 20.})

    # get midtime
    midtime = get_midtime(target, "Kepler")

    # check if midtime is correct
    assert midtime == 20.

    # -------------------------------------------------------------------------
    # Test valid input with TESS parameters

    # create input target as pandas Series with pl_tranmid_tess column
    target = pd.Series({"pl_tranmid_tess": 2457000. + 20.})

    # get midtime
    midtime = get_midtime(target, "TESS")

    # check if midtime is correct
    assert midtime == 20.

    # -------------------------------------------------------------------------
    # Test invalid input with Kepler parameters

    # create input target as pandas Series with empty pl_tranmid column
    target = pd.Series({"pl_tranmid": np.nan})

    with pytest.raises(ValueError) as e:
        get_midtime(target, "Kepler")
        # check if error message is correct.
        assert "pl_tranmid" in str(e.value)

    # -------------------------------------------------------------------------
    # Test invalid input with TESS parameters

    # create input target as pandas Series with empty pl_tranmid_tess column
    # and empty pl_tranmid column
    target = pd.Series({"pl_tranmid_tess": np.nan,
                        "pl_tranmid": np.nan})

    with pytest.raises(ValueError) as e:
        get_midtime(target, "TESS")
        # check if error message is correct.
        assert "pl_tranmid_tess" in str(e.value)
        assert "pl_tranmid" in str(e.value)

    # -------------------------------------------------------------------------
    # Test input with no TESS parameters but mission being TESS 
    # and Kepler parameters being present

    # create input target as pandas Series without pl_tranmid_tess column
    # but with pl_tranmid column
    target = pd.Series({"pl_tranmid_tess": np.nan,
                        "pl_tranmid": 2457000 + 20.})

    # get midtime
    tranmid = get_midtime(target, "TESS")

    # check if midtime is correct
    assert tranmid == 20.










