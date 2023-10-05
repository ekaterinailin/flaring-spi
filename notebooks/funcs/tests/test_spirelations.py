import pytest
import numpy as np
import pandas as pd
from astropy.constants import R_sun
from astropy import units as u


from ..spirelations import (convert_datapoints_to_obstime,
                            eq45_lanza2012,
                            p_spi_lanza12,
                            b_from_lx_reiners,
                            calculate_relative_velocity)



def test_convert_datapoints_to_obstime():
    """Test that conversion works as expected,
    and that error are thrown at the right places."""


    # 1 TESS data point is two minute
    assert convert_datapoints_to_obstime(1, "TESS") * 24. * 60. == 2. 

    # 1 Kepler data point is one minute
    assert convert_datapoints_to_obstime(1, "Kepler") * 24. * 60. == 1.

    # 1 TESS data point is 20 seconds
    assert convert_datapoints_to_obstime(1e5, "TESS") * 24. * 60. * 60. == pytest.approx(20. * 1e5)


    with pytest.raises(ValueError):
        convert_datapoints_to_obstime(1, "Chandra")
        
    with pytest.raises(ValueError):
        convert_datapoints_to_obstime(1e5, "Chandra")


def test_eq45_lanza2012():
    """Tables 1 and 3 from Lanza et al. (2012), also notes from the text to fill
    in the missing values. Results should be consistent with 2%.

    """

    # from Table 1
    sname = ["HD 179949", "HD 189733", "ups And", "tau Boo"]
    B = np.array([10, 10, 10, 10]) # G
    Bp = np.array([100, 100, 100, 100]) # G
    a = [0.045, 0.031, 0.059, 0.046] # AU
    a_over_R = [7.72, 8.56, 10.16, 7.38] 
    vrel = [88.4, 125.4, 71.4, 22.8] # km/s

    # from Table 3
    res_tab3 = [1.23e26, 1.04e26, 0.35e26, 0.38e26] # erg/s

    # planet radii from text in 3.2
    Rp = [1.2, 1.127, 1.2, 1.2] # R_jup

    # make DataFrame from lists
    df = pd.DataFrame({"sname": sname, "Bs": B, "Bp": Bp, "a": a, "Rp": Rp,
                    "a_over_R": a_over_R, "vrel": vrel, "res_tab3": res_tab3})

    # calculate Rstar from a/R and R in units of R_sun
    df["Rs"] = df["a"]*u.AU.to("m") / df["a_over_R"] / R_sun

    # use the eigenvalue and n from 3.3
    lambda_05 = .82343 # lambda squared
    n = 0.5

    # calculate SPI power
    df["power_SPI_erg_s"] = df.apply(lambda x: eq45_lanza2012(x.Bs, x.Bp, x.Rs, 
                                                            x.Rp, x.a, x.vrel, n, 
                                                            lambda_05), axis=1)


    # calculate ratio with Table 3 results
    df["power_SPI_erg_s_over_res_tab3"] = df["power_SPI_erg_s"] / df["res_tab3"]

    # check that the results are within 2% of the values in Table 3
    assert (np.round(df["power_SPI_erg_s_over_res_tab3"].values, 2) ==
            [0.99, 0.99, 0.98, 0.98]).all()

    # check that Bp=0 gives a positive value
    assert (df.apply(lambda x: eq45_lanza2012(x.Bs, 0, x.Rs, x.Rp, x.a, x.vrel, n, 
                                              lambda_05,), axis=1).values > 0.).all()  

def test_p_spi_lanza12():
    """Test SPI power function with and without errors"""

    # check that no error reduces to the original Lanza equation.
    assert p_spi_lanza12(1, 1, 1, 1, 1, Bp=1.) == eq45_lanza2012(1, 1, 1, 1, 1, 1, 0.25, 1.01203)

    res = p_spi_lanza12(1, 1, 1, 1, 1, Bp=1., error=True, 
                        v_rel_err=0, Bhigh=1, Blow=0, 
                        pl_radhigh=1, pl_radlow=0, a_err=0,
                        Bp_err=0, rstarhigh=1, rstarlow=0, 
                        n=0.25, lambda_squared=1.01203)
    
    assert res[1] == res[0]
    assert res[2] == 0

    res = p_spi_lanza12(1, 1, 1, 1, 1, Bp=1., error=True, 
                        v_rel_err=.1, Bhigh=2, Blow=0, 
                        pl_radhigh=1., pl_radlow=0, a_err=0.,
                        Bp_err=.1, rstarhigh=1., rstarlow=0, 
                        n=0.25, lambda_squared=1.01203)
    
    assert res[1] > 0
    assert res[1] > res[0]
    assert res[2] == 0

    res = p_spi_lanza12(1, 1, 1, 1, 1, Bp=1., error=True,
                        v_rel_err=0, Bhigh=1, Blow=0.9,
                        pl_radhigh=1, pl_radlow=0.9, a_err=0,
                        Bp_err=0, rstarhigh=1, rstarlow=0.9, 
                        n=0.25, lambda_squared=1.01203)

    assert res[1] == res[0]
    assert res[2] > 0

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., Bhigh=0.1, error=True)

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., Blow=1.1, error=True)

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., pl_radhigh=0.1, error=True)

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., pl_radlow=1.1, error=True)

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., a_err=-1.1, error=True)

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., Bp_err=-1.1, error=True)

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., v_rel_err=-1.1, error=True)

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., rstarhigh=0.1, error=True)

    with pytest.raises(AssertionError):
        p_spi_lanza12(1, 1., 1., 1., 1., rstarlow=1.1, error=True)

    with pytest.raises(ValueError):
        p_spi_lanza12(1, 1., 1., 1., 1., Bp=-1.)

    res = p_spi_lanza12(1, 1, 1, 1, 1, Bp=0., error=True, 
                        v_rel_err=0, Bhigh=1, Blow=0, 
                        pl_radhigh=1, pl_radlow=0, a_err=0,
                        Bp_err=0, rstarhigh=1, rstarlow=0, 
                        n=0.25, lambda_squared=1.01203)
    
    assert res[0] > 0
    assert res[1] == res[0]
    assert res[2] == 0


def test_b_from_lx_reiners():
    """Tiny test for Reiners et al. 2022 Eq. 4"""

    # input unit conversion check
    assert b_from_lx_reiners(3.28e-12, 1. / (R_sun.to("cm").value)) * 4. * np.pi == 1.

    # input NaN returns NaN
    assert np.isnan(b_from_lx_reiners(3.28e-12, np.nan))

    # input 0 returns nan
    assert np.isnan(b_from_lx_reiners(3.28e-12, 0))

    # error propagation
    lx, r = 3.28e-12, 1. / (R_sun.to("cm").value)
    lx_err, r_err = 0. * lx, 0. * r
    B, highB, lowB = b_from_lx_reiners(lx, r, error=True, lx_err=lx_err, r_err=r_err)

    # as before
    assert B == 1. / (4. * np.pi)

    # if input value is 1 output is 1 always
    assert highB == B
    assert lowB == B



def test_calculate_relative_velocity():
    "Test if the formula converts properly."
    # synchronized orbit
    assert calculate_relative_velocity(1, 1, 1) == 0

    # no orbital distance, planet is in the star
    assert np.isnan(calculate_relative_velocity(0, 1, 3))

    # unit check
    assert (calculate_relative_velocity(1 /2 / np.pi, 1, 2) / 
            149597870.700 * 24. * 3600. == 0.5)

    # error propagation
    assert (calculate_relative_velocity(1, 1, 1, error=True,
                                        a_au_err=1, orbper_err=1, rotper_err=1) ==
            (0, np.sqrt(8) * np.pi / 24. / 3600. * 149597870.700))

    
    # nans
    assert np.isnan(calculate_relative_velocity(1, 1, 0))
    assert np.isnan(calculate_relative_velocity(1, 0, 1))
    assert np.isnan(calculate_relative_velocity(np.nan, 1, 1))
    assert np.isnan(calculate_relative_velocity(np.nan, 1, 1, error=True)).all()

