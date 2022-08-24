import os
import numpy as np
import pandas as pd
import pytest

from ..phaseanalysis import (get_observed_phases,
                             get_flare_phases,
                             get_cumulative_distributions,)


def test_get_cumulative_distributions():
    """Test two setup of two-Sector datasets.
    
    """
    # 20 flares in each of two sectors
    # In this setup we could ignore havin two Sectors involved 
    # because the flare frequency in both is the same 

    qcs = ["99"] * 20 +  ["9"] * 20
    get_secs_cadences = [("99",1),("9",1)]

    # define data with flares marked by Sector
    df = pd.DataFrame({"qcs": qcs})

    # define random observing time in each of the 40 phase bins
    np.random.seed(42)
    dfphases = pd.DataFrame({"99": np.random.rand(len(qcs)) * 10,
                             "9": np.random.rand(len(qcs)) * 10,})
    # calculate distributions
    n_i, n_exp, cum_n_exp, cum_n_i = get_cumulative_distributions(df, dfphases, get_secs_cadences)

    # measured distribution is all ones
    assert (n_i == np.ones(40)).all()

    # cumulative distribution should be incrementing from first phase to the last
    assert cum_n_i == pytest.approx(np.linspace(0.025,1,40))

    # expected number of flares should be roughly 1 flare per bin, since we distributed the observing times randomly
    assert n_exp.mean() == pytest.approx(1)

    # difference between the measured and expected distribution should be small 
    # as they are drawn form the uniform distribution
    assert (np.abs(cum_n_exp - cum_n_i) < .05).all()

    # 100 + 300 flares 
    # In this setup we cannot ignore having two Sectors involved 
    # because the flare frequency in Sector 9 is mcuh higher than in sector 99, 
    # 6 times as high to be precise 

    qcs = ["99"] * 100 +  ["9"] * 300
    get_secs_cadences = [("99",1),("9",1)]

    # define data with flares marked by Sector
    df = pd.DataFrame({"qcs": qcs})

    # define random observing time in each of the 40 phase bins
    np.random.seed(42)
    dfphases = pd.DataFrame({"99": np.random.rand(len(qcs)) * 10,
                             "9": np.concatenate([np.zeros(len(qcs)//2),(np.random.rand(len(qcs)//2)) * 10])})
    # calculate distributions
    n_i, n_exp, cum_n_exp, cum_n_i = get_cumulative_distributions(df, dfphases, get_secs_cadences)

    # measured distribution is all ones
    assert (n_i == np.ones(len(qcs))).all()

    # cumulative distribution should be incrementing from first phase to the last
    assert cum_n_i == pytest.approx(np.linspace(1/len(qcs),1,len(qcs)))

    # expected number of flares should be roughly 1 flare per bin, since we distributed the observing times randomly
    assert n_exp.mean() == pytest.approx(1)

    # the kink in the distribution is at phase 0.5
    assert np.argmax((np.abs(cum_n_exp - cum_n_i))) == 200


def test_get_observed_phases():
    
    # set up synthetic light curve 
    time = np.linspace(10,20,10000)
    PER = 2.
    lc = pd.DataFrame({"time" : time, 
                       "phase" : (time % PER) / PER})
    lc.to_csv(f"../results/observedtimes/AU Mic_999_0_TESS.csv")

    # 
    mode = "Rotation"
    get_secs_cadences = [(999, time[1]-time[0])]

    p = np.random.rand(300)

    # fails if phases are not ordered
    with pytest.raises(ValueError):
        get_observed_phases("Rotation", p, get_secs_cadences, rotper=PER)
        
    # fails if phases are not ordered
    with pytest.raises(ValueError):
        get_observed_phases("Orbit", p, get_secs_cadences, rotper=PER)

    # we defined the light curve phases such that 
    # both modes give the same results
    for mode in ["Rotation", "Orbit"]:
        
        # TEST NR 1
        
        p = np.random.rand(300)

        # fails if phases are not ordered
        with pytest.raises(ValueError):
            get_observed_phases(mode, p, get_secs_cadences, rotper=PER)
            
        # TEST NR 2
        
        p = np.linspace(0.0198,.998,50)
        phasesobs, binmids = get_observed_phases(mode, p, get_secs_cadences, 
                                                 rotper=PER)
        # check that observing bins are filled with the right amount of data
        # the rhs value is the time resolution of the lc 
        # which should give the precision limit
        print((time[-1] - time[0]) / len(p))
        print(phasesobs.iloc[0],len(p))
        print(get_secs_cadences[0][1])
        assert (((time[-1] - time[0]) / len(p) - phasesobs)[1:-1] 
                < get_secs_cadences[0][1]).T.values.all()

        # we get one more bins than we get flares
        assert binmids.shape[0] == len(p)+1



def test_get_flare_phases():
    # set up synthetic light curve 
    time = np.linspace(10,20,10001)
    PER = 2.
    lc = pd.DataFrame({"time" : time, 
                       "phase" : (time % PER) / PER})
    lc.to_csv(f"AU Mic_999_0_TESS.csv")

    # define flare table with start times, and giving the Sector to pick it up
    df = pd.DataFrame({"tstart":np.linspace(10,20,21), "qcs":"999"})

    # apply function to both modes
    for mode in ["Orbit", "Rotation"]:
        ph = get_flare_phases(df, mode, rotper=PER, path="./")

        # test exact results
        assert (ph.phases.values == np.array([0., .25, .5, 0.75] * 5 + [0])).all()
    
    # remove the created file in the end
    os.remove("./AU Mic_999_0_TESS.csv")
