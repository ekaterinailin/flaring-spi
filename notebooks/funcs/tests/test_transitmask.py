import numpy as np
import pandas as pd
import pytest

from altaipony.flarelc import FlareLightCurve

from ..transitmask import get_full_transit_mask, get_transit_mid_epochs



def test_get_full_transit_mask():
    flc = FlareLightCurve(time=np.linspace(10,30,101))

    system = pd.DataFrame({"pl_trandur":[12.,20.,np.nan],
                           "pl_tranmidepoch" : [9., 37., np.nan], 
                           "pl_orbper": [5.,10.,19.]})

    # if the table contains NaN, demand cleanup of system parameters
    error_message = ("Some planets in your system are not transiting"
                             ", or lack necessary info about transit epoch"
                             ", period, and duration.")
    with pytest.raises(ValueError, match=error_message):
        get_full_transit_mask(system, flc)

    # drop last row everything should work again    
    system = system.iloc[:2]

    # the results are manually verified here:
    mask = get_full_transit_mask(system, flc)

    assert (flc.time[mask] == 
            pytest.approx(np.array([13.8, 14. , 14.2, 16.6, 16.8, 17. , 
                                    17.2, 17.4, 18.8, 19. , 19.2, 23.8, 
                                    24. , 24.2, 26.6, 26.8, 27. , 27.2, 
                                    27.4, 28.8, 29. , 29.2])))
    assert mask.dtype == "bool"

    # Check overlapping transits are treated okay
    system = pd.DataFrame({"pl_trandur":[12., 20.],
                           "pl_tranmidepoch" : [9., 9.], 
                           "pl_orbper": [5.,10.]})

    mask = get_full_transit_mask(system, flc)
    flc.time[mask]

    # These are also manually confirmed
    assert (flc.time[mask] == 
            pytest.approx(np.array([13.8, 14. , 14.2, 18.6, 18.8,
                                    19. , 19.2, 19.4, 23.8, 24. ,
                                    24.2, 28.6, 28.8, 29. , 29.2, 
                                    29.4])))
    
    # Check that single transits are treated okay
    system = pd.DataFrame({"pl_trandur":[12., 20.],
                           "pl_tranmidepoch" : [9., 9.], 
                           "pl_orbper": [5.,np.nan]})

    mask = get_full_transit_mask(system, flc)
    flc.time[mask]

    # These are also manually confirmed
    assert (flc.time[mask] == 
            pytest.approx(np.array([13.8, 14. , 14.2, 18.8,
                                    19. , 19.2, 23.8, 24. ,
                                    24.2, 28.8, 29. , 29.2,])))

def test_get_transit_mid_epochs(): 
    assert (get_transit_mid_epochs(10., 3., 20., 40.) ==
            np.array([19., 22., 25., 28., 31., 34., 37., 40.])).all()

    assert (get_transit_mid_epochs(46., 3., 20., 40.) ==
            np.array([19., 22., 25., 28., 31., 34., 37., 40.])).all()

    with pytest.raises(ValueError) as e:
        get_transit_mid_epochs(46., np.nan, 20., 40.)

    with pytest.raises(ValueError) as e:
        get_transit_mid_epochs(np.nan, 3., 20., 40.)

    with pytest.raises(ValueError) as e:
        get_transit_mid_epochs(np.nan, 3., 20., np.nan)

    with pytest.raises(ValueError) as e:
        get_transit_mid_epochs(np.nan, 3., np.nan, 40.)
        
    # test a very long orbital period for single transit case
    assert (get_transit_mid_epochs(30., 300., 20., 40.) ==
        np.array([-270.,   30.,  330.])).all()
