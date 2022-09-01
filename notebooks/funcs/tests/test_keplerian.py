import pytest

from astropy.constants import M_jup, M_earth
from astropy import units as u

from ..keplerian import get_semimajor_axis

def test_get_semimajor_axis():
    """Two integration tests to see if semimajor axis
    is computed sanely in the case of the Sun-Earth system.
    """

    # check that the Earth is placed at 1 AU
    a = get_semimajor_axis(365.25, 1, M_earth/M_jup).to(u.AU).value
    assert a == pytest.approx(1., rel=1e-4)

    # check that the Earth's mass is negligible compared to the Sun's
    a = get_semimajor_axis(365.25, 1, 0).to(u.AU).value
    assert a == pytest.approx(1., rel=1e-4)