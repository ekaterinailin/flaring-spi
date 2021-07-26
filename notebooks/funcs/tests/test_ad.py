import pytest

import numpy as np

from ..ad import (anderson_darling_statistic,
                  get_pvalue_from_AD_statistic,
                  sample_AD_for_custom_distribution,
                 )

def test_anderson_darling_statistic():

    # These should fail because the AD statistic does not work on [0,1]
    c = np.array([[0.,.2, .3]])
    with pytest.raises(ValueError):
        anderson_darling_statistic(c)

    c = np.array([[1.,.2, .3]])
    with pytest.raises(ValueError):
        anderson_darling_statistic(c)

    # This is the ideal case of a perfect uniform distribution
    # in this case the AD statistic becomes AD = - len(c) + 4 sum( c * log (c) )
    # see Jaentschi and Bolboaca 2018
    c = np.array([[.05, .15,.25, .35, .45, .55, .65, .75, .85, .95]])
    assert (anderson_darling_statistic(c)[0] == 
            pytest.approx(-len(c[0]) - 4 * np.sum(c * np.log(c))))

    # This one just shows that the function works on 2D-samples, too
    c = np.random.rand(100,10)
    assert anderson_darling_statistic(c).shape[0] == 100
    
    

def test_get_pvalue_from_AD_statistic():
    A2 = np.random.normal(2., .4, 10000)
    cphases = np.linspace(0.05,0.95,100)

    p, a = get_pvalue_from_AD_statistic(cphases, A2)

    # this value should not depend on the randomness of A2
    assert a == pytest.approx(.8329200828193137)

    # high significance because .83 is about 3 sigma outlier in A2
    assert p < .01

    A2 = np.random.normal(.4, .1, 10000)

    p, a = get_pvalue_from_AD_statistic(cphases, A2)
    # high significance because .83 is a >4 sigma outlier in A2
    assert p < .01
    
    

def test_sample_AD_for_custom_distribution():
    # define observed flare phases
    cphases = np.linspace(1e-2,1 - 1e-2, 100)

    # define observing duration for each phase
    cobs = np.random.rand(100)/10 + 100
    cobs[70:80] -= 3

    # sample from distribution
    A2  = sample_AD_for_custom_distribution(cphases, cobs, 1000)

    # number of samples is preserved
    assert A2.shape[0] == 1000

    # AD statistic is always positive
    assert (A2 > 0).all()
