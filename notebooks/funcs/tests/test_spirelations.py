import pytest
import numpy as np
from astropy.constants import R_jup, R_sun


from ..spirelations import (convert_datapoints_to_obstime,
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


def test_p_spi_lanza():
    """Tiny test for SPI power function"""
    # input unit conversion check
    assert p_spi_lanza12(1e-6, 1., 1./ R_jup.to("cm").value) == 1.

    # input NaN returns NaN
    with pytest.raises(ValueError):
        np.isnan(p_spi_lanza12(1e-6, 1., 1./ R_jup.to("cm").value, Bp=np.nan))


def test_b_from_lx_reiners():
    """Tiny test for Reiners et al. 2022 Eq. 4"""

    # input unit conversion check
    assert b_from_lx_reiners(3.28e-12, 1. / (R_sun.to("cm").value)) * 4. * np.pi == 1.

    # input NaN returns NaN
    np.isnan(b_from_lx_reiners(3.28e-12, np.nan))


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

