import numpy as np
from astropy.constants import R_jup, R_sun

from ..spirelations import (p_spi_lanza12,
                            b_from_lx_reiners)

def test_p_spi_lanza():
    """Tiny test for SPI power function"""
    # input unit conversion check
    assert p_spi_lanza12(1e-6, 1., 1./ R_jup.to("cm").value) == 1.

    # input NaN returns NaN
    assert np.isnan(p_spi_lanza12(1e-6, 1., 1./ R_jup.to("cm").value, Bp=np.nan))


def test_b_from_lx_reiners():
    """Tiny test for Reiners et al. 2022 Eq. 4"""

    # input unit conversion check
    (3.28e-12, 1. / (R_sun.to("cm").value)) * 4. * np.pi == 1.

    # input NaN returns NaN
    np.isnan(b_from_lx_reiners(3.28e-12, np.nan))