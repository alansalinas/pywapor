import numpy as np
from pywapor.et_look_dev import constants as con


def longitude_rad(lon_deg):
    r"""
    Converts longitude from degrees to radians.

    Parameters
    ----------
    lon_deg : float
        longitude in degrees
        :math:`\phi`
        [deg]

    Returns
    -------
    lon : float
        longitude
        :math:`\phi`
        [rad]

    """
    return lon_deg * np.pi/180.0


def latitude_rad(lat_deg):
    r"""
    Converts latitude from degrees to radians.

    Parameters
    ----------
    lat_deg : float
        latitude in degrees
        :math:`\lambda`
        [deg]

    Returns
    -------
    lat : float
        latitude
        :math:`\lambda`
        [rad]

    """
    return lat_deg * np.pi/180.0


def slope_rad(slope_deg):
    r"""
    Converts slope from degrees to radians.

    Parameters
    ----------
    slope_deg : float
        slope in degrees
        :math:`s`
        [deg]

    Returns
    -------
    slope : float
        slope
        :math:`\Delta`
        [rad]

    """
    return slope_deg * np.pi/180.0


def aspect_rad(aspect_deg):
    r"""
    Converts aspect from degrees to radians.

    Parameters
    ----------
    aspect_deg : float
        aspect in degrees
        :math:`s`
        [deg]

    Returns
    -------
    aspect : float
        aspect (0 is north; pi is south)
        :math:`\alpha`
        [rad]
    """
    return aspect_deg * np.pi/180.0


def declination(doy):
    r"""
    Computes the solar declination which is the angular height of the sun
    above the astronomical equatorial plane in radians

    .. math ::
        \delta=0.409\sin\left(\frac{2\pi J}{365}-1.39\right)

    Parameters
    ----------
    doy : float
        julian day of the year
        :math:`J`
        [-]

    Returns
    -------
    decl : float
        declination
        :math:`\delta`
        [rad]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> solrad.declination(180)
    0.40512512455439242
    """
    B = 360./365 * (doy-81) 
    decl = np.arcsin(np.sin(np.deg2rad(23.45))*np.sin(np.deg2rad(B)))
    
    return decl

def inverse_earth_sun_distance(doy):
    r"""
    Computes the inverse earth sun distance (iesd) in Angstrom Unit where 1 AU is 1.496e8 km

    .. math ::
        d_{r}=1+0.033\cos\left(\frac{2\pi J}{365}\right)

    Parameters
    ----------
    doy : float
        julian day of the year
        :math:`J`
        [-]

    Returns
    -------
    iesd : float
        inverse earth sun distance
        :math:`d_{r}`
        [AU]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> solrad.inverse_earth_sun_distance(180)
    0.96703055420162642
    """

    return 1 + 0.033 * np.cos(doy * 2 * np.pi / 365.0)


def actual_earth_sun_distance(iesd):
    r"""
    Computes the earth sun distance (esd) in Angstrom Unit where 1 AU is 1.496e8 km

    .. math ::
        d_{r}=1+0.033\cos\left(\frac{2\pi J}{365}\right)

    Parameters
    ----------
    iesd : float
        inverse earth sun distance
        :math:`d_{i}_{r}`
        [AU]

    Returns
    -------
    esd : float
        earth sun distance
        :math:`d_{r}`
        [AU]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> solrad.actual_earth_sun_distance(180)
    1.034093489244084
    """

    return 1 / iesd


def seasonal_correction(doy):
    r"""
    Computes the seasonal correction for solar time  in hours

    .. math ::
        b=\frac{2\pi\left(J-81\right)}{364}

    .. math ::
        s_{c}=0.1645sin\left(2b\right)-0.1255cos\
              left(b\right)-0.025\left(b\right)

    Parameters
    ----------
    doy : float
        julian day of the year
        :math:`J`
        [-]

    Returns
    -------
    sc : float
        seasonal correction
        :math:`s_{c}`
        [hours]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> solrad.seasonal_correction(180)
    -0.052343379605521212
    """
    b = 2 * np.pi * (doy - 81) / 365.0
    return 0.1645 * np.sin(2 * b) - 0.1255 * np.cos(b) - 0.025 * np.sin(b)


def sunset_hour_angle(lat, decl):
    r"""
    Computes the sunset hour angle

    .. math ::
        w_{s}=\arccos(-\tan(\lambda)\tan(\delta))

    Parameters
    ----------
    decl : float
        solar declination
        :math:`\delta`
        [rad]
    lat : float
        latitude
        :math:`\lambda`
        [rad]

    Returns
    -------
    ws : float
        sunset hour angle
        :math:`w_{s}`
        [rad]

    """
    return np.arccos(-(np.tan(lat) * np.tan(decl)))


def hour_angle(sc, dtime, lon=0):
    r"""
    Computes the hour angle which is zero at noon and -pi at 0:00 am and
    pi at 12:00 pm

    .. math ::
        \omega=\left(\frac{\pi}{12}\right)\left(t+s_{c}-12\right)

    Parameters
    ----------
    sc : float
        seasonal correction
        :math:`s_{c}`
        [hours]
    dtime : float
        decimal time
        :math:`t`
        [hours]
    lon : float
        longitude
        :math:`\phi`
        [rad]

    Returns
    -------
    ha : float
        hour_angle
        :math:`\omega`
        [rad]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> solrad.hour_angle(sc=solrad.seasonal_correction(75), dtime=11.4)
    -0.19793970172084141
    """
    dtime = dtime #+ (lon / (15*np.pi/180.0))
    return (np.pi / 12.0) * (dtime + sc - 12.0)


def inst_solar_radiation_toa(csza, iesd):
    r"""
    Computes the instantaneous solar radiation at the top of
    the atmosphere [Wm-2]

    .. math ::
        S_{toa}^{i} = S_{sun}d_{r}\phi

    Parameters
    ----------
    csza : float
        cosine solar zenith angle
        :math:`\phi`
        [-]
    iesd : float
        inverse earth sun distance
        :math:`d_{r}`
        [AU]

    Returns
    -------
    ra_i_toa : float
        instantaneous solar radiation at top of atmosphere
        :math:`S_{toa}^{i}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> doy = 1
    >>> sc = solrad.seasonal_correction(doy)
    >>> ha = solrad.hour_angle(sc, dtime=12)
    >>> decl = solrad.declination(doy)
    >>> csza = solrad.cosine_solar_zenith_angle(ha, decl, 0)
    >>> iesd = solrad.inverse_earth_sun_distance(doy)
    >>> solrad.inst_solar_radiation_toa(csza, iesd)
    1299.9181944414036
    """
    return csza * con.sol * iesd


def daily_solar_radiation_toa(sc, decl, iesd, lat, slope=0, aspect=0):
    r"""
    Computes the daily solar radiation at the top of the atmosphere.

    .. math ::
        S_{toa}=S_{sun}d_{r}\int_{i=-\pi}^{i=\pi}S_{toa}^{i}

    Parameters
    ----------
    iesd : float
        inverse earth sun distance
        :math:`d_{r}`
        [AU]
    decl : float
        solar declination
        :math:`\delta`
        [rad]
    sc : float
        seasonal correction
        :math:`s_{c}`
        [hours]
    lat : float
        latitude
        :math:`\lambda`
        [rad]
    slope : float
        slope
        :math:`\Delta`
        [rad]
    aspect : float
        aspect (0 is north; pi is south)
        :math:`\alpha`
        [rad]

    Returns
    -------
    ra_24_toa : float
        daily solar radiation at the top of atmosphere
        :math:`S_{toa}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> from math import pi
    >>> doy = 1
    >>> sc = solrad.seasonal_correction(doy)
    >>> decl = solrad.declination(doy)
    >>> iesd = solrad.inverse_earth_sun_distance(doy)
    >>> solrad.daily_solar_radiation_toa(sc, decl, iesd, lat=25*pi/180.0)
    265.74072308978026
    """

    # hour angle for the whole day in half-hourly intervals (0:15-23:45)
    t_start = 0.25
    t_end = 24.00
    interval = 0.5
    times = [t_start+i*interval for i in range(0, 48)]
    hours = [hour_angle(sc, t) for t in times]

    ra24 = 0

    for t in hours:
        csza = cosine_solar_zenith_angle(t, decl, lat, slope, aspect)
        ra24 += inst_solar_radiation_toa(csza, iesd) / len(hours)

    # return the average daily radiation
    return ra24


def cosine_solar_zenith_angle(ha, decl, lat, slope=0, aspect=0):
    r"""
    computes the cosine of the solar zenith angle [-]

    .. math ::
        \phi = \begin{array}{c}
        \sin\delta\sin\lambda\cos\Delta-\\
        \sin\delta\cos\lambda\sin\Delta+\\
        \cos\delta\cos\lambda\cos\Delta\cos\left(\omega\right)+\\
        \cos\delta\sin\lambda\sin\Delta\sin\alpha\cos\left(\omega\right)+\\
        \cos\delta\sin\Delta\sin\alpha\sin\left(\omega\right)
        \end{array}

    Parameters
    ----------
    ha : float
        hour angle
        :math:`\omega`
        [rad]
    decl : float
        declination
        :math:`\delta`
        [rad]
    lat : float
        latitude
        :math:`\lambda`
        [rad]
    slope : float
        slope
        :math:`\Delta`
        [rad]
    aspect : float
        aspect (0 is north; pi is south)
        :math:`\alpha`
        [rad]

    Returns
    -------
    csza : float
        cosine solar zenith angle
        :math:`\phi`
        [-]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> sc = solrad.seasonal_correction(1)
    >>> ha = solrad.hour_angle(sc, dtime=12)
    >>> solrad.cosine_solar_zenith_angle(ha, decl=solrad.declination(1), lat=0)
    0.92055394167363314
    """
    t1 = np.sin(decl) * np.sin(lat) * np.cos(slope)
    t2 = np.sin(decl) * np.cos(lat) * np.sin(slope) * np.cos(aspect - np.pi)
    t3 = np.cos(decl) * np.cos(lat) * np.cos(slope)
    t4 = np.cos(decl) * np.sin(lat) * np.sin(slope) * np.cos(aspect - np.pi)
    t5 = np.cos(decl) * np.sin(slope) * np.sin(aspect - np.pi)
    csza = t1 - t2 + t3 * np.cos(ha) + t4 * np.cos(ha) + t5 * np.sin(ha)

    # check if the sun is above the horizon
    check = np.sin(decl) * np.sin(lat) + np.cos(decl) * \
        np.cos(lat) * np.cos(ha)

    nans = np.logical_or(np.isnan(csza), np.isnan(check))

    res = np.where(np.logical_and(csza > 0, check >= 0), csza, 0)
    res[nans] = np.nan

    return res


def daily_solar_radiation_toa_flat(decl, iesd, lat, ws):
    r"""
    Computes the daily solar radiation at the top of the atmosphere for a flat
    surface.

    .. math ::
        S_{toa,f}=\frac{S_{sun}}{\pi}d_{inv,r}*(w_{s}\sin(\lambda))\sin(\delta) +
                  \cos(\lambda))\cos(\delta)\sin(w_{s})

    Parameters
    ----------
    decl : float
        solar declination
        :math:`\delta`
        [rad]
    iesd : float
        inverse earth sun distance
        :math:`d_{inv,r}`
        [AU]
    lat : float
        latitude
        :math:`\lambda`
        [rad]
    ws : float
        sunset hour angle
        :math:`w_{s}`
        [rad]

    Returns
    -------
    ra_24_toa_flat : float
        daily solar radiation at the top of atmosphere for a flat surface
        :math:`S_{toa,f}`
        [Wm-2]

    """
    ra_flat = (con.sol / np.pi) * iesd * (ws * np.sin(lat) * np.sin(decl) +
                                        np.cos(lat) * np.cos(decl) * np.sin(ws))

    return ra_flat


def daily_solar_radiation_flat(ra_24_toa_flat, trans_24):
    r"""
    Computes the daily solar radiation at the earth's surface

    .. math ::
        S^{\downarrow} = \tau S_{toa}

    Parameters
    ----------
    ra_24_toa_flat : float
        daily solar radiation at the top of atmosphere for a flat surface
        :math:`S_{toa}`
        [Wm-2]
    trans_24 : float
        daily atmospheric transmissivity
        :math:`\tau`
        [-]

    Returns
    -------
    ra_24 : float
        daily solar radiation for a flat surface
        :math:`S^{\downarrow}`
        [Wm-2]

    """
    return ra_24_toa_flat * trans_24


def diffusion_index(trans_24, diffusion_slope=-1.33, diffusion_intercept=1.15):
    r"""
    Computes the diffusion index, the ratio between diffuse and direct
    solar radiation. The results are clipped between 0 and 1.

    .. math ::
        I_{diff} = a_{diff}+b_{diff}\tau

    Parameters
    ----------
    trans_24 : float
        daily atmospheric transmissivity
        :math:`\tau`
        [-]
    diffusion_slope : float
        slope of diffusion index vs transmissivity relationship
        :math:`b_{diff}`
        [-]
    diffusion_intercept : float
        intercept of diffusion index vs transmissivity relationship
        :math:`a_{diff}`
        [-]

    Returns
    -------
    diffusion_index : float
        diffusion_index
        :math:`I_{diff}`
        [-]

    """
    res = diffusion_intercept + trans_24 * diffusion_slope

    res = np.clip(res, 0, 1)

    return res


def daily_total_solar_radiation(ra_24_toa, ra_24_toa_flat, diffusion_index, trans_24):
    r"""
    Computes the daily solar radiation at the earth's surface taken
    diffuse and direct solar radiation into account

    .. math ::
        S^{\downarrow} = I_{diff} \tau S_{toa,f} +(1-I_{diff}) \tau S_{toa}

    Parameters
    ----------
    ra_24_toa : float
        daily solar radiation at the top of atmosphere
        :math:`S_{toa}`
        [Wm-2]
    ra_24_toa_flat : float
        daily solar radiation at the top of atmosphere for a flat surface
        :math:`S_{toa,f}`
        [Wm-2]
    diffusion_index : float
        diffusion_index
        :math:`I_{diff}`
        [-]
    trans_24 : float
        daily atmospheric transmissivity
        :math:`\tau`
        [-]

    Returns
    -------
    ra_24 : float
        daily solar radiation
        :math:`S^{\downarrow}`
        [Wm-2]

    """
    diffuse = trans_24 * ra_24_toa_flat * diffusion_index
    direct = trans_24 * ra_24_toa * (1 - diffusion_index)
    return diffuse + direct

def daily_solar_radiation_toa_new(sc, decl, iesd, lat, doy, slope=0, aspect=0):
    r"""
    Computes the daily solar radiation at the top of the atmosphere.

    .. math ::
        S_{toa}=S_{sun}d_{r}\int_{i=-\pi}^{i=\pi}S_{toa}^{i}

    Parameters
    ----------
    iesd : float
        inverse earth sun distance
        :math:`d_{r}`
        [AU]
    decl : float
        solar declination
        :math:`\delta`
        [rad]
    sc : float
        seasonal correction
        :math:`s_{c}`
        [hours]
    lat : float
        latitude
        :math:`\lambda`
        [rad]
    slope : float
        slope
        :math:`\Delta`
        [rad]
    aspect : float
        aspect (0 is north; pi is south)
        :math:`\alpha`
        [rad]

    Returns
    -------
    ra_24_toa : float
        daily solar radiation at the top of atmosphere
        :math:`S_{toa}`
        [Wm-2]

    Examples
    --------
    >>> import ETLook.solar_radiation as solrad
    >>> from math import pi
    >>> doy = 1
    >>> sc = solrad.seasonal_correction(doy)
    >>> decl = solrad.declination(doy)
    >>> iesd = solrad.inverse_earth_sun_distance(doy)
    >>> solrad.daily_solar_radiation_toa(sc, decl, iesd, lat=25*pi/180.0)
    265.74072308978026
    """
    #print(type(slope))
    if type(slope)==int:
        slope = np.zeros(lat.shape)
    #print(type(aspect))    
    if type(aspect)==int:
        aspect = np.zeros(lat.shape)
        
    gamma = np.deg2rad(np.rad2deg(aspect)-180)                               # Surface aspect angle (radians)
    a,b,c = Constants(decl,slope,gamma,lat)

    ra24 = np.zeros(np.shape(lat))*np.nan
    dr = 1 + 0.033 * np.cos(doy*2*np.pi/365)  # Inverse relative distance Earth-Sun
    constant=con.sol*dr/(2*np.pi)
    TwoPeriod= TwoPeriods(decl,slope,lat)  # all input in radians

    #2.) calculate the 24-hours extraterrestrial radiation (2 periods)
    ID = np.where(np.ravel(TwoPeriod==True))
    ra24.flat[ID]=TwoPeriodSun(constant, decl, slope.flat[ID], gamma.flat[ID], lat.flat[ID])

    #3.) calculate the 24-hours extraterrestrial radiation (1 period)
    ID = np.where(np.ravel(TwoPeriod==False))
    ra24.flat[ID]=OnePeriodSun(constant, decl, slope.flat[ID], gamma.flat[ID], lat.flat[ID])

    # Horizontal surface
    ws = np.arccos(-np.tan(decl) * np.tan(lat))  # Sunrise/sunset time angle

    # Extraterrestial radiation for a horizontal surface for 24-h period:
    Ra_hor_24 = (con.sol * dr / np.pi * (np.sin(decl) * np.sin(lat) * ws + np.cos(decl) * np.cos(lat) * np.sin(ws)))
    # cos_theta_flat = (np.sin(delta) * np.sin(phi) + np.cos(delta) * np.cos(phi) * np.cos(w))

    # Mountain radiation
    ra24 = np.where(ra24 > 0.1 * Ra_hor_24, ra24 / np.cos(slope),
                           Ra_hor_24)
    ra24[ra24 > 600.0] = 600.0
    
    return ra24

def OnePeriodSun(constant,delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006
    Calculate the 24-hours extraterrestrial radiation when there is one sun period
    '''
    sunrise,sunset = SunHours(delta,s,gamma,phi)
    Vals=IntegrateSlope(constant,sunrise,sunset,delta,s,gamma,phi)

    return(Vals)

def TwoPeriodSun(constant,delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006
    Calculate the 24-hours extraterrestrial radiation when there are two sun period
    '''
    A1, A2 = SunHours(delta,s,gamma,phi)
    a,b,c = Constants(delta,s,gamma,phi)
    riseSlope, setSlope = BoundsSlope(a,b,c)
    B1 = np.maximum(riseSlope,setSlope)
    B2 = np.minimum(riseSlope,setSlope)
    Angle_B1 = AngleSlope(a,b,c,B1)
    Angle_B2 = AngleSlope(a,b,c,B2)

    B1[abs(Angle_B1) > 0.001] = np.pi - B1[abs(Angle_B1) > 0.001]
    B2[abs(Angle_B2) > 0.001] = -np.pi - B2[abs(Angle_B2) > 0.001]

    # Check if two periods really exist
    ID = np.ravel_multi_index(np.where(np.logical_and(B2 >= A1, B1 >= A2) == True),a.shape)
    Val = IntegrateSlope(constant,B2.flat[ID],B1.flat[ID],delta,s.flat[ID],gamma.flat[ID],phi.flat[ID])
    ID = ID[Val < 0]

    # Finally calculate resulting values
    Vals = np.zeros(B1.shape)

    Vals.flat[ID] = (IntegrateSlope(constant,A1.flat[ID],B2.flat[ID],delta,s.flat[ID],gamma.flat[ID],phi.flat[ID])  +
                   IntegrateSlope(constant,B1.flat[ID],A2.flat[ID],delta,s.flat[ID],gamma.flat[ID],phi.flat[ID]))
    ID = np.ravel_multi_index(np.where(Vals == 0),a.shape)
    Vals.flat[ID] = IntegrateSlope(constant,A1.flat[ID],A2.flat[ID],delta,s.flat[ID],gamma.flat[ID],phi.flat[ID])

    return(Vals)

#------------------------------------------------------------------------------
def IntegrateSlope(constant,sunrise,sunset,delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006 equation 5
    Calculate the 24 hours extraterrestrial radiation
    '''
    # correct the sunset and sunrise angels for days that have no sunset or no sunrise
    SunOrNoSun = np.logical_or(((np.abs(delta + phi)) > (np.pi/2)),((np.abs(delta - phi)) > (np.pi/2)))
    integral=np.zeros(s.shape)
    ID = np.where(np.ravel(SunOrNoSun==True))

    # No sunset
    IDNoSunset = np.where(np.ravel(abs(delta+phi.flat[ID])>(np.pi/2)))
    if np.any(IDNoSunset) == True:
        sunset1=np.pi
        sunrise1=-np.pi
        integral.flat[IDNoSunset] = constant * (np.sin(delta)*np.sin(phi)*np.cos(s)*(sunset1-sunrise1)
            - np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma)*(sunset1-sunrise1)
            + np.cos(delta)*np.cos(phi)*np.cos(s)*(np.sin(sunset1)-np.sin(sunrise1))
            + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)*(np.sin(sunset1)-np.sin(sunrise1))
            - np.cos(delta)*np.sin(s)*np.sin(gamma)*(np.cos(sunset1)-np.cos(sunrise1)))

    # No sunrise
    elif np.any(IDNoSunset) == False:
        integral.flat[IDNoSunset==False]=constant * (np.sin(delta)*np.sin(phi)*np.cos(s)*(0)
            - np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma)*(0)
            + np.cos(delta)*np.cos(phi)*np.cos(s)*(np.sin(0)-np.sin(0))
            + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)*(np.sin(0)-np.sin(0))
            - np.cos(delta)*np.sin(s)*np.sin(gamma)*(np.cos(0)-np.cos(0)))

    ID = np.where(np.ravel(SunOrNoSun==False))
    integral.flat[ID] = constant * (np.sin(delta)*np.sin(phi)*np.cos(s)*(sunset-sunrise)
            - np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma)*(sunset-sunrise)
            + np.cos(delta)*np.cos(phi)*np.cos(s)*(np.sin(sunset)-np.sin(sunrise))
            + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)*(np.sin(sunset)-np.sin(sunrise))
            - np.cos(delta)*np.sin(s)*np.sin(gamma)*(np.cos(sunset)-np.cos(sunrise)))

    return(integral)

def BoundsHorizontal(delta,phi):
    ''''
    Based on Richard G. Allen 2006
    This function calculates sunrise hours based on earth inclination and latitude
    If there is no sunset or sunrise hours the values are either set to 0 (polar night) or pi (polar day)
    '''
    bound = np.arccos(-np.tan(delta)*np.tan(phi))
    bound[abs(delta+phi) > np.pi/2] = np.pi
    bound[abs(delta-phi) > np.pi/2] = 0

    return(bound)

def TwoPeriods(delta,s,phi):
    '''
    Based on Richard G. Allen 2006
    Create a boolean map with True values for places with two sunsets
    '''
    TwoPeriods = (np.sin(s) > np.ones(phi.shape)*np.sin(phi)*np.cos(delta)+np.cos(phi)*np.sin(delta))

    return(TwoPeriods)

def SunHours(delta,slope,slopedir,lat):
    # Define sun hours in case of one sunlight period

    a,b,c = Constants(delta,slope,slopedir,lat)
    riseSlope, setSlope = BoundsSlope(a,b,c)
    bound = BoundsHorizontal(delta,lat)

    Calculated = np.zeros(slope.shape, dtype = bool)
    RiseFinal = np.zeros(slope.shape)
    SetFinal = np.zeros(slope.shape)

    # First check sunrise is not nan
    # This means that their is either no sunrise (whole day night) or no sunset (whole day light)
    # For whole day light, use the horizontal sunrise and whole day night a zero..
    Angle4 = AngleSlope(a,b,c,-bound)
    RiseFinal[np.logical_and(np.isnan(riseSlope),Angle4 >= 0.0)] = -bound[np.logical_and(np.isnan(riseSlope),Angle4 >= 0.0)]
    Calculated[np.isnan(riseSlope)] = True

    # Step 1 > 4
    Angle1 = AngleSlope(a,b,c,riseSlope)
    Angle2 = AngleSlope(a,b,c,-bound)

    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(Angle2 < Angle1+0.001 ,Angle1 < 0.001),Calculated == False) == True),a.shape)
    RiseFinal.flat[ID] = riseSlope.flat[ID]
    Calculated.flat[ID] = True
    # step 5 > 7
    Angle3 = AngleSlope(a,b,c,-np.pi - riseSlope)

    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(-bound<(-np.pi-riseSlope),Angle3 <= 0.001),Calculated == False) == True),a.shape)
    RiseFinal.flat[ID] = -np.pi -riseSlope.flat[ID]
    Calculated.flat[ID] = True

    # For all other values we use the horizontal sunset if it is positive, otherwise keep a zero
    RiseFinal[Calculated == False] = -bound[Calculated == False]

    # Then check sunset is not nan or < 0
    Calculated = np.zeros(slope.shape, dtype = bool)

    Angle4 = AngleSlope(a,b,c,bound)
    SetFinal[np.logical_and(np.isnan(setSlope),Angle4 >= 0.0)] = bound[np.logical_and(np.isnan(setSlope),Angle4 >= 0.0)]
    Calculated[np.isnan(setSlope)] = True

    # Step 1 > 4
    Angle1 = AngleSlope(a,b,c,setSlope)
    Angle2 = AngleSlope(a,b,c,bound)

    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(Angle2 < Angle1+0.001,Angle1 < 0.001),Calculated == False) == True),a.shape)
    SetFinal.flat[ID] = setSlope.flat[ID]
    Calculated.flat[ID] = True
    # step 5 > 7
    Angle3 = AngleSlope(a,b,c,np.pi - setSlope)

    ID = np.ravel_multi_index(np.where(np.logical_and(np.logical_and(bound>(np.pi-setSlope),Angle3 <= 0.001),Calculated == False) == True),a.shape)
    SetFinal.flat[ID] = np.pi - setSlope.flat[ID]
    Calculated.flat[ID] = True

    # For all other values we use the horizontal sunset if it is positive, otherwise keep a zero
    SetFinal[Calculated == False] = bound[Calculated == False]

    #    Angle4 = AngleSlope(a,b,c,bound)
    #    SetFinal[np.logical_and(Calculated == False,Angle4 >= 0)] = bound[np.logical_and(Calculated == False,Angle4 >= 0)]

    # If Sunrise is after Sunset there is no sunlight during the day
    SetFinal[SetFinal <= RiseFinal] = 0.0
    RiseFinal[SetFinal <= RiseFinal] = 0.0

    return(RiseFinal,SetFinal)

def AngleSlope(a,b,c,w):
    '''
    Based on Richard G. Allen 2006
    Calculate the cos zenith angle by using the hour angle and constants
    '''
    angle = -a + b*np.cos(w) + c*np.sin(w)

    return(angle)

def BoundsSlope(a,b,c):
    '''
    Based on Richard G. Allen 2006 equation 13
    This function calculates candidate values for sunrise and sunset hour angles
    '''
    Div = (b**2+c**2)
    Div[Div <= 0] = 0.00001
    sinB = (a*c + b*np.sqrt(b**2+c**2-a**2)) / Div
    sinA = (a*c - b*np.sqrt(b**2+c**2-a**2)) / Div

    sinB[sinB < -1] = -1; sinB[sinB > 1] = 1    # Limits see appendix A.2.i
    sinA[sinA < -1] = -1; sinA[sinA > 1] = 1    # Limits see appendix A.2.i

    sunrise = np.arcsin(sinA)
    sunset = np.arcsin(sinB)

    return(sunrise,sunset)

def Constants(delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006 equation 11
    determines constants for calculating the exterrestial solar radiation
    '''
    a = np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma) - np.sin(delta)*np.sin(phi)*np.cos(s)
    b = np.cos(delta)*np.cos(phi)*np.cos(s) + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)
    c = np.cos(delta)*np.sin(s)*np.sin(gamma)

    return(a,b,c)    
    
