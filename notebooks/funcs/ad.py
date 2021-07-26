"""
UTF-8, Python 3

------------
Flaring SPI
------------

Ekaterina Ilin, 2021, MIT License


This module contains functions used to 
calculate and evaluate the Anderson-Darling 
statistic for custom distributions and arbitrary
sample sizes using emcee to sample from
the custom flare phase distribution.
"""

import numpy as np
import emcee

from scipy.interpolate import interp1d
from scipy.stats import percentileofscore


def anderson_darling_statistic(c):
    """Compute the AD statistic 
    from a sample of flare phases.
    
    Parameters:
    ------------
    c : n, m-array
        sample of flare phases
    
    Return:
    -------
    m-array - AD statistic of c
    """
    if (c == 0.).any() | (c == 1.).any():
        raise ValueError("AD statistic is defined on (0,1), not on [0,1]!")
    
    # To calculate the AD statistic the sample must be sorted
    # in ascendin order
    c = np.sort(c)
    
    # calculate and return AD statistic
    n = len(c[0])
    i = np.arange(n)
    lis = (2. * (i + 1) - 1.) * np.log(c) + (2 * n + 1 - 2 * (i + 1)) * np.log(1 - c)
    
    return -n - (1 / n) * np.sum(lis, axis=1)



def get_pvalue_from_AD_statistic(cphases, A2):
    """Calculate the p-value of the deviation of the
    observed phases using the Ensemble sampling for 
    the expected distribution of the AD-statistic.
    
    Parameters:
    ------------
    cphases : n-array
        n-long array of observed flare phases
    A2 : m-array
        distribution of the AD statistic (m values)
        
    Returns:
    ---------
    pval, atest - p-value and measured AD statistic
    """
    # AD statistic of the observed flare phases
    atest = anderson_darling_statistic(cphases.reshape(1, len(cphases)))[0]
    
    # percentile of the AD distribution the statistic falls into
    perc = percentileofscore(A2, atest)
    
    # calculate p-value for a two sided test
    pval = 2 * np.min([100 - perc, perc]) / 100
    
    return pval, atest


    
def sample_AD_for_custom_distribution(cphases, cobs, N):
    """
    
    Parameters:
    ------------
    cphases : n-array
        phases at which a flare was observed
    cobs : n-array
        number of observing data point 
        between the last two flares
    N : int
        number of samples to draw from the expected 
        distribution of flare phases
    """

    # Interpolate!
    # for values near 0 and near 1  - extrapolate! 
    f = interp1d(cphases, cobs, bounds_error=False, 
                 fill_value=(cobs[0],cobs[-1]))

    # plot diagnostic with the interpolated callable cdf for the KS test
#     plt.figure(figsize=(8, 3))
#     x = np.sort(np.random.rand(100))
#     y = f(x)
#     plt.plot(x,y)
#     plt.scatter(cphases, cobs)
#     plt.xlim(0,1)

    # define the of the distribution function to draw samples from
    def func(x):
        if (x > 1) | (x < 0):
            return -np.inf
        else:
            return np.log(f(x))
        
    # number of flares in the sample 
    # crucial to get the correct critical values    
    nobs = len(cphases) 
    
    # Apply ensemble sampler from emcee
    
    # Sample from a 1-D distribution
    # Sample nobs values from the distribution
    # because the distrbution of the AD statistic
    # crucially depends on the sample size
    ndim, nwalkers = 1, nobs
    
    # initial state of the sampler is random values between 0 and 1
    p0 = np.random.rand(nwalkers, ndim)
    
    # Define the sampler with func as the distribution to sample from
    sampler = emcee.EnsembleSampler(nwalkers, ndim, func)
    
    # Run MCMC for N steps
    sampler.run_mcmc(p0, N, progress=True)
    
    # Get the samples
    samples = sampler.get_chain()

    # remove third dimension in samples: (N,nobs,1) --> (N, nobs)
    c = samples.reshape((N,nobs))
    
    # calculate the AD statistic for all samples
    A2 = anderson_darling_statistic(c)
    
    # make sure that converting shapes and calculating the statistic
    # preserved the number of samples correctly
    assert len(A2) == N
    
    return A2
