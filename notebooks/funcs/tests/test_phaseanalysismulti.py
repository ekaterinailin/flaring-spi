import pytest
import pandas as pd
import numpy as np

from ..phaseanalysismulti import get_cumulative_distribution

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
    with pytest.raises(KeyError) as e:
        lcs["mission"] = ["TESS"] * 3
        get_cumulative_distribution(flares, observedphases, lcs)
        assert "column not found in observedphases" in str(e.value)