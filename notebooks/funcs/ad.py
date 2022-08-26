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
import corner

from scipy.interpolate import interp1d
from scipy.stats import percentileofscore
from scipy.misc import derivative


import matplotlib.pyplot as plt

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



def get_pvalue_from_AD_statistic(x, dist, A2):
    """Calculate the p-value of the deviation of the
    observed phases using the Ensemble sampling for 
    the expected distribution of the AD-statistic.
    
    Parameters:
    ------------
    x : n-array
        n-long array of observed flare phases
    dist : func
        expected cum. distribution func (EDF)
    A2 : m-array
        distribution of the AD statistic (m values)
        
    Returns:
    ---------
    pval, atest - p-value and measured AD statistic
    """
#    print(x)
    # AD statistic of the observed flare phases
    atest = anderson_custom(x, dist)
    
    # percentile of the AD distribution the statistic falls into
    perc = percentileofscore(A2, atest)
    
    # calculate p-value for a two sided test
    pval = 2 * np.min([100 - perc, perc]) / 100

    
    return pval, atest


    
def sample_AD_for_custom_distribution(f, nobs, N, savefig=False,
                                      fig_name="AD_statistic.png"):
    """
    
    Parameters:
    ------------
    f : func
        expected cum. dist. function (EDF)
    nobs : int
        size of data sample
    N : int
        number of samples to draw from the expected 
        distribution of flare phases
    savefig : bool
        whether to save the figure
    fig_name : str
        name of the figure to save

    Returns:
    ---------
    A2 : N-array
    """

   
    def func(x):
        if (x >= 1) | (x <= 0):
            return -np.inf
        else:
            lp =  np.log(derivative(f, x, dx=np.min([1e-4, (1-x)/3, x/3])))
            if np.isnan(lp):
                return -np.inf
            else:
                return lp
        

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

    # replace infs
    missing = np.where(~np.isfinite(samples))[0]
    if len(missing)>0:
        samples[missing] = f(np.random.rand(*missing.shape))
    c = samples.reshape((N,nobs))

    # Make a figure of the sampled distribution
    fig = corner.corner(sampler.get_chain(discard=100, thin=15, flat=True))
    yt = list(plt.yticks()[0])
    plt.xlabel("phase")
    
    if savefig == True:
        fig.savefig(fig_name, dpi=300)
 
    # calculate the AD statistic for all samples
    A2 = np.array([anderson_custom(c[i,:], f) for i in range(N)])
    
    # make sure that converting shapes and calculating the statistic
    # preserved the number of samples correctly
    assert len(A2) == N
    
    return A2



def anderson_custom(x, dist):
    """
    Anderson-Darling test for data coming from a particular distribution
    The Anderson-Darling test is a modification of the Kolmogorov-
    Smirnov test `kstest` for the null hypothesis that a sample is
    drawn from a population that follows a particular distribution.
    For the Anderson-Darling test, the critical values depend on
    which distribution is being tested against.  
    Parameters
    ----------
    x : array_like
        array of sample data
    dist : func
        epected cum. distribution func (EDF)
    Returns
    -------
    A2 : float
        The Anderson-Darling test statistic
    """

    y = np.sort(x)
    z = dist(y)

    # A2 statistic is undefined for 1 and 0
    z = z[(z<1.) & (z>0.)]
  
    N = len(z)
    i = np.arange(1, N + 1)
    S = np.sum((2 * i - 1.0) / N * (np.log(z) + np.log(1 - z[::-1])), axis=0)
    A2 = - N - S

    return A2
