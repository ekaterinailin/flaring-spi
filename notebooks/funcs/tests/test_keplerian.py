import pytest

import numpy as np
from astropy.constants import M_jup, M_earth
from astropy import units as u

from ..keplerian import get_semimajor_axis, get_distance_from_planet_range

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


def test_get_distance_from_planet_range():

    # test the case where eccentricity is not known

    # ... and the error on the semimajor axis is small
    assert get_distance_from_planet_range(1, np.nan, .1) == 0.25

    # ... and the error on the semimajor axis is large
    assert get_distance_from_planet_range(1, np.nan, .5) == 0.5

    # ... and the error on the semimajor axis is nan
    assert get_distance_from_planet_range(1, np.nan, np.nan) == 0.25

    # test the case where eccentricity is known to be 0.1

    # ... and the error on the semimajor axis is larger than eccentricity
    assert get_distance_from_planet_range(1, 0.1, .1) == 0.1

    # ... and the error on the semimajor axis is smaller than eccentricity
    assert get_distance_from_planet_range(1, 0.1, .001) == 1 - np.sqrt(0.99)

    # ... and the error on the semimajor axis is nan
    assert get_distance_from_planet_range(1, 0.1, np.nan) == 1 - np.sqrt(0.99)