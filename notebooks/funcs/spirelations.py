import numpy as np
import astropy.units as u
from astropy.constants import R_jup, R_sun


def convert_datapoints_to_obstime(n, mission):
    """Calculate the observing time in days.
    based on the number of datapoints in the LC.

    Parameters
    ----------
    n : int
        Number of datapoints.
    mission : str
        Mission name, either "Kepler" or "TESS". 

    Returns
    -------
    obstime : float
        Observing time in days.
    """
    if n > 50000:
        if mission == "TESS":
            return n / 3. / 60. / 24.  # 20 sec cadence
        elif mission == "Kepler":
            return n / 60. / 24.  # 1 min cadence
        else:
            raise ValueError("Mission name not recognized for a large LC.")
    else:
        if mission == "TESS":
            return n * 2. / 60. / 24.  # 2 min cadence
        elif mission == "Kepler":
            return n / 60. / 24.  # 1 min cadence
        elif mission == "K2":
            return n / 60. / 24.  # 1 min cadence
        else:
            raise ValueError("Mission name not recognized for a small LC.")


def wrap_obstimes(TIC, flares):
    """Wrap the obstime function for use with multiple
    LCs for a given TIC ID. NOT TESTED!

    Parameters
    ----------
    TIC : str
        TIC ID.
    flares : pandas.DataFrame
        Flare table with total_n_valid_data_points and mission columns.

    Returns
    -------
    obstime_d: float
        Total observing time in days.
    """
    # pick up all the rows with LCs searched for flares
    sysvalues = flares[flares.TIC.astype(str) == TIC]

    # drop the duplicates with multiple flares per LC
    vals = sysvalues[["total_n_valid_data_points",
                      "mission"]].drop_duplicates()
  
    # calculate the observing time based on the 
    # number of data points in the LC
    # and summing up the observing time for all LCs
    obstime_d = vals.apply(lambda x: convert_datapoints_to_obstime(
        x.total_n_valid_data_points, x.mission), axis=1).sum()

    return obstime_d



def p_spi_lanza12(v_rel, B, pl_rad, a, rstar, Bp=1., error=False, 
                  v_rel_err=None, Bhigh=None, Blow=None, 
                  pl_radhigh=None, pl_radlow=None, a_err=None,
                  Bp_err=None, rstarhigh=None, rstarlow=None):
    """Lanza 2012 scaling relation and its adaptation to absence of
    planet magnetic field (Bp=0).

    Parameters
    ----------
    v_rel : float
        Relative velocity between stellar rotation and orbital velocity in km/s.
    B : float
        X-ray luminosity in erg/s.
    pl_rad : float
        Planet radius in R_jup.
    a : float
        Semi-major axis in AU.
    rstar : float
        Stellar radius in R_sun.
    Bp : float
        Planet magnetic field in Gauss.
    error : bool
        If True, return the error in the power of the star-planet interaction
        through estimating max and min values of the relation.
    v_rel_err : float
        Error in relative velocity between stellar rotation and orbital velocity
        in km/s.
    Bhigh : float
        Upper limit of magnetic field in Gauss.
    Blow : float
        Lower limit of magnetic field in Gauss.       
    pl_radhigh : float
        Upper limit of planet radius in R_jup.
    pl_radlow : float
        Lower limit of planet radius in R_jup.
    a_err : float
        Error in semi-major axis in AU.
    Bp_err : float
        Error in planet magnetic field in Gauss.
    rstarhigh : float
        Upper limit of stellar radius in R_sun.
    rstarlow : float
        Lower limit of stellar radius in R_sun.


    Returns
    -------
    p_spi, p_spi_err : float, float
        Power of the star-planet interaction in erg/s, error in erg/s.
    """
    # --------------------
    # Unit conversion:

    # convert from Rjup to cm
    pl_rad = pl_rad * R_jup.to(u.km).value

    # convert from rsun to cm
    rstar = rstar * R_sun.to(u.km).value

    # convert a from AU to km
    a = a * u.AU.to(u.km)


    if error:
        # convert from Rjup to cm
        pl_radhigh = pl_radhigh * R_jup.to(u.km).value
        pl_radlow = pl_radlow * R_jup.to(u.km).value

        # convert rstar, rstarhigh, rstarlow from rsun to cm
      
        rstarhigh = rstarhigh * R_sun.to(u.km).value
        rstarlow = rstarlow * R_sun.to(u.km).value

        #convert a_err from AU to km
        a_err = a_err * u.AU.to(u.km)


    # --------------------

    # --------------------
    # Three cases: Bp > 0 and Bp = 0 and Bp < 0

    if Bp > 0.:

        pspi = B**(4./3.) * v_rel * pl_rad**2 * Bp**(2./3.) / (a**4) * (rstar**4)
        
        if error:
            pspi_high = Bhigh**(4./3.) * (v_rel + v_rel_err) * pl_radhigh**2 * (Bp + Bp_err)**(2./3.) * rstarhigh**4 /  (a - a_err)**4 
            pspi_low = Blow**(4./3.) * (v_rel - v_rel_err) * pl_radlow**2 * (Bp - Bp_err)**(2./3.)  * rstarlow**4 /  (a + a_err)**4
            return pspi, pspi_high, pspi_low
        
        else:
            return pspi
    
    elif Bp == 0.:

        pspi = B**2 * v_rel * pl_rad**2  / (a**4) * (rstar**4)

        if error:
            pspi_high = Bhigh**2 * (v_rel + v_rel_err) * pl_radhigh**2  * rstarhigh**4 / (a-a_err**4) 
            pspi_low = Blow**2 * (v_rel - v_rel_err) * pl_radlow**2 * rstarlow**4 / (a+a_err**4) 
            return pspi, pspi_high, pspi_low
    
        else:
            return pspi
    
    else:
        raise ValueError("Planet magnetic field must be >= 0.")

def b_from_lx_reiners(lx, r, error=False, lx_err=None, r_err=None):
    """Reiners et al. 2022 Eq. 4 inverted to get magnetic field in Gauss.
    
    Parameters
    ----------
    lx : float
        X-ray luminosity in erg/s.
    r : float
        Stellar radius in R_sun.
    error : bool
        If True, return the error in the magnetic field by
        estimating it from the scatter in the relation, and the
        error in the stellar radius and X-ray luminosity.
    lx_err : float
        Error in X-ray luminosity in erg/s.
    r_err : float
        Error in stellar radius in R_sun.

    Returns
    -------
    B : float
        Magnetic field in Gauss.
    B_err : float
        Error in magnetic field in Gauss.
    """
    if np.isnan(lx) | np.isnan(r) | (lx == 0) | (r == 0):
        B = np.nan
        highB, lowB = np.nan, np.nan
        if error:
            return B, highB, lowB
        else:
            return B
    else:
        # convert stellar radius to cm
        rcm = (r * R_sun).to(u.cm).value

        # constants from Reiners et al. 2022 Eq. 4    
        a, b = 1 / 3.28e-12, 1.58

        # error on b from Reiners et al. 2022 Eq. 4
        b_err = 0.06

        # calculate magnetic field from inverse of Reiners et al. 2022 Eq. 4
        get_B = lambda lx, rcm, a, b: (a * lx)**(1/b) / (4. * np.pi * rcm**2)

        B =  get_B(lx, rcm, a, b) 

        if error:

            # convert error in radius
            rcm_err = (r_err * R_sun).to(u.cm).value

            # calculate upper and lower 1-sigma range on magnetic field
            highB = get_B(lx+lx_err, rcm+rcm_err, a, b - b_err)
            lowB = get_B(lx-lx_err, rcm-rcm_err, a, b + b_err)

            return B, highB, lowB
        else:
            return B

    

def calculate_relative_velocity(a_au, orbper, rotper, error=False,
                                a_au_err=None, orbper_err=None, rotper_err=None):
    """Calculate the relative velocity between stellar rotation at the planetary orbit
    and the orbital distance of the planet in km/s.

    Parameters
    ----------
    au : float
        Semi-major axis in AU.
    orbper : float
        Orbital period in days.
    rotper : float
        Stellar rotation period in days.
    error : bool
        If True, return the error in the relative velocity from
        error propagation.
    a_au_err : float
        Error in semi-major axis in AU.
    orbper_err : float
        Error in orbital period in days.
    rotper_err : float
        Error in stellar rotation period in days.

    Returns
    -------
    v_rel_km_s : float
        Relative velocity between stellar rotation at the planetary orbit
        and the orbital velocity of the planet in km/s.
    """
    
    # return in km/s
    # 1AU = 149597870.700 km
    # 1day = 24h * 3600s

    # convet au to km
    a_au = a_au * 149597870.700

    if (np.isnan(a_au) | np.isnan(orbper) | np.isnan(rotper) | 
        (a_au == 0) | (orbper == 0) | (rotper == 0)):
        v_rel = np.nan
        v_rel_err = np.nan

    else:
        # convert orbper to s
        orbper = orbper * 24. * 3600.

        # convert rotper to s
        rotper = rotper * 24. * 3600.
        
        v_rel = 2 * np.pi * a_au * (1/orbper - 1/rotper)

        if error:
            # get error in km/s
            dv_dau = 2 * np.pi * (1/orbper - 1/rotper)
            dv_dorbper = -2 * np.pi * a_au * (1/orbper**2)
            dv_drotper = 2 * np.pi * a_au * (1/rotper**2)

            # convert a_au_err to km
            a_au_err = a_au_err * 149597870.700

            # convert orbper_err to s
            orbper_err = orbper_err * 24. * 3600.

            # convert rotper_err to s
            rotper_err = rotper_err * 24. * 3600.

            # quadratic error propagation
            v_rel_err = np.sqrt((dv_dau * a_au_err)**2 +
                                (dv_dorbper * orbper_err)**2 +
                                (dv_drotper * rotper_err)**2)

    if error:

        return v_rel, v_rel_err

    else:
        return v_rel



def rossby_reiners2014(Lbol, Prot, error=False, Lbol_high=None, 
                       Lbol_low=None, Prot_high=None, Prot_low=None):
    """Calculate the Rossby number for a given bolometric luminosity.
    Based on Reiners et al. (2014), as noted in Reiners et al. (2022).
    
    Parameters
    ----------
    Lbol : float
        Bolometric luminosity in solar units.
    Prot : float
        Rotation period in days.
    error : bool
        If True, return the error in the Rossby number from
        error propagation.
    Lbol_high : float
        Upper 1-sigma error in bolometric luminosity in solar units.
    Lbol_low : float
        Lower 1-sigma error in bolometric luminosity in solar units.
    Prot_high : float
        Upper 1-sigma error in rotation period in days.
    Prot_low : float    
        Lower 1-sigma error in rotation period in days.


    Returns
    -------
    Ro : float
        Rossby number.
    """
    # convective turnover time
    tau = 12.3 / (Lbol**0.5)

    if error:
        tau_high = 12.3 / (Lbol_low**0.5)
        tau_low = 12.3 / (Lbol_high**0.5)

    # Rossby number
    Ro = Prot / tau

    if error:
        Ro_high = Prot_high / tau_low
        Ro_low = Prot_low / tau_high

    if error:
        return Ro, Ro_high, Ro_low
    else:
        return Ro


def b_from_ro_reiners2022(Ro, error=False, Ro_high=None, Ro_low=None):
    """Calculate the manetic field from the Rossby number.
    Based on Reiners et al. (2022), Table 2.
    
    Parameters
    ----------
    Ro : float
        Rossby number.
    error : bool
        If True, return the error in the magnetic field from
        scatter in relation.
    Ro_high : float
        Upper 1-sigma limit in Rossby number.
    Ro_low : float
        Lower 1-sigma limit in Rossby number.

    Returns
    -------
    B : float
        Magnetic field in Gauss.
    """

    # slow rotator
    if Ro > 0.13:
        B = 199 * Ro**(-1.26)#pm .1
        if error:
            if Ro>=1:
                B_high = 199 * Ro_low**(-1.26 + 0.1)
                B_low = 199 * Ro_high**(-1.26 - 0.1)
            else:
                B_low = 199 * Ro_high**(-1.26 + 0.1)
                B_high = 199 * Ro_low**(-1.26 - 0.1)
            
    # fast rotator
    elif Ro < 0.13:
        B = 2050 * Ro**(-0.11) #pm 0.03
        if error:
            B_low = 2050 * Ro_low**(-0.11 + 0.03)
            B_high = 2050 * Ro_high**(-0.11 - 0.03)
    else:
        B, B_high, B_low = np.nan, np.nan, np.nan

    if error:
        return B, B_high, B_low
    else:
        return B



def pspi_kavanagh2022(Rp, B, vrel, a, Bp=1., error=False, Rphigh=None, Bphigh=1.,
                      Bhigh=None, vrelhigh=None, alow=None, Rplow=None, Bplow=1., 
                      Blow=None, vrellow=None, ahigh=None):
    """Power of star-plaet interactions following the
    Saur et al. 2013 model, put in a scaling law by
    Kavanagh et al. 2022.
    
    Parameters
    ----------
    Rp : float
        Planet radius in Jupiter radii
    B : float
        Stellar magnetic field in G
    vrel : float
        Relative velocity in km/s
    a : float
        Orbital separation in AU
    Bp : float
        Planet magnetic field in G
    error : bool
        If True, return the error on the pspi
    Rphigh : float
        Upper limit on planet radius in Jupiter radii
    Bphigh : float
        Upper limit on planet magnetic field in G
    Bhigh : float   
        Upper limit on stellar magnetic field in G
    vrelhigh : float
        Upper limit on relative velocity in km/s
    alow : float    
        Lower limit on orbital separation in AU
    Rplow : float
        Lower limit on planet radius in Jupiter radii
    Bplow : float
        Lower limit on planet magnetic field in G
    Blow : float
        Lower limit on stellar magnetic field in G
    vrellow : float
        Lower limit on relative velocity in km/s
    ahigh : float
        Upper limit on orbital separation in AU

    Returns
    -------
    pspi : float
        prop. to power of star-planet interaction
    """
    # convert a from AU to km
    a = a * u.AU.to(u.km)

    # convert Rp from Jupiter radii to km
    Rp = Rp * u.R_jup.to(u.km)


    pspi = Rp**2 * Bp**(2/3) * B**(1/3) * vrel**2 * a**(-2)

    if error:
        # convert alow and ahigh from AU to km
        alow = alow * u.AU.to(u.km)
        ahigh = ahigh * u.AU.to(u.km)

        # convert Rplow and Rphigh from Jupiter radii to km
        Rplow = Rplow * u.R_jup.to(u.km)
        Rphigh = Rphigh * u.R_jup.to(u.km)

        pspi_high = Rphigh**2 * Bphigh**(2/3) * Bhigh**(1/3) * vrelhigh**2 * alow**(-2)
        pspi_low = Rplow**2 * Bplow**(2/3) * Blow**(1/3) * vrellow**2 * ahigh**(-2)
    
        return pspi, pspi_high, pspi_low

    else:
        return pspi