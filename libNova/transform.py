
from math import cos, sin, asin, acos, atan2
from libNova.ln_types import *
from libNova.utility import *
from libNova.sidereal_time import *
from libNova.precession import *

# def ln_get_rect_from_helio( object ) -> ln_rect_posn
# param object Object heliocentric coordinates
# param position Pointer to store new position

# Transform an objects heliocentric ecliptical coordinates
# into heliocentric rectangular coordinates.


def ln_get_rect_from_helio( object : ln_helio_posn ) -> ln_rect_posn:
    posn = ln_rect_posn()

    # ecliptic J2000
    sin_e = 0.397777156
    cos_e = 0.917482062

    # calc common values
    cos_B = cos(ln_deg_to_rad(object.B))
    cos_L = cos(ln_deg_to_rad(object.L))
    sin_B = sin(ln_deg_to_rad(object.B))
    sin_L = sin(ln_deg_to_rad(object.L))

    # equ 37.1
    posn.X = object.R * cos_L * cos_B
    posn.Y = object.R * (sin_L * cos_B * cos_e - sin_B * sin_e)
    posn.Z = object.R * (sin_L * cos_B * sin_e + sin_B * cos_e)
    return posn


# def ln_get_hrz_from_equ( object, observer, JD) -> ln_hrz_posn
# param object Object coordinates.
# param observer Observer cordinates.
# param JD Julian day
# returns position Pointer to store new position.

# Transform an objects equatorial coordinates into horizontal coordinates for the given julian day and observers position.

# 0 deg azimuth = south, 90 deg = west.

def ln_get_hrz_from_equ( object : ln_equ_posn, observer : ln_lnlat_posn, JD : float ) -> ln_hrz_posn:
    posn = ln_hrz_posn()

    # get mean sidereal time in hours*
    sidereal = ln_get_mean_sidereal_time(JD)
    posn = ln_get_hrz_from_equ_sidereal_time (object, observer, sidereal)
    return posn



def ln_get_hrz_from_equ_sidereal_time( object : ln_equ_posn, observer : ln_lnlat_posn, sidereal : float ) -> ln_hrz_posn:
    posn = ln_hrz_posn()

    # change sidereal_time from hours to radians
    sidereal *= 2.0 * M_PI / 24.0

    # calculate hour angle of object at observers position
    ra = ln_deg_to_rad(object.ra)
    H = sidereal + ln_deg_to_rad(observer.lng) - ra

    # hence formula 12.5 and 12.6 give
    # convert to radians - hour angle, observers latitude, object declination
    latitude = ln_deg_to_rad(observer.lat)
    declination = ln_deg_to_rad(object.dec)

    # formula 12.6 * missuse of A (you have been warned)
    A = sin(latitude) * sin(declination) + cos(latitude) * cos(declination) * cos(H)
    h = asin(A)

    # convert back to degrees
    posn.alt = ln_rad_to_deg(h)

    # zenith distance, Telescope Control 6.8a
    z = acos(A)
    zs = sin(z)

    # sane check for zenith distance don't try to divide by 0
    if fabs(zs) < 1e-5:
        if object.dec > 0.0:
            posn.az = 180.0
        else:
            posn.az = 0.0
        if (object.dec > 0.0 and observer.lat > 0.0) or (object.dec < 0.0 and observer.lat < 0.0):
            posn.alt = 90.0
        else:
            posn.alt = -90.0
        return posn

    As = (cos(declination) * sin(H)) / zs
    Ac = (sin(latitude) * cos(declination) * cos(H) - cos(latitude) * sin(declination)) / zs

    # don't blom at atan2
    if Ac == 0.0 and As == 0.0:
        if object.dec > 0:
            posn.az = 180.0
        else:
            posn.az = 0.0
        return posn

    A = atan2(As, Ac)

    # convert back to degrees
    posn.az = ln_range_degrees(ln_rad_to_deg(A))
    return posn


# def ln_get_equ_from_hrz( object, observer, JD) -> ln_equ_posn
# param object Object coordinates.
# param observer Observer cordinates.
# param JD Julian day
# returns position Pointer to store new position.

# Transform an objects horizontal coordinates into equatorial coordinates for the given julian day and observers position.

def ln_get_equ_from_hrz( object : ln_hrz_posn, observer : ln_lnlat_posn, JD : float ) -> ln_equ_posn:
    position = ln_equ_posn()

    # change observer/object position into radians

    # object alt/az
    A = ln_deg_to_rad(object.az)
    h = ln_deg_to_rad(object.alt)

    # observer long / lat
    longitude = ln_deg_to_rad(observer.lng)
    latitude = ln_deg_to_rad(observer.lat)

    # equ on pg89
    H = atan2(sin(A), (cos(A) * sin(latitude) + tan(h) * cos(latitude)))
    declination = sin(latitude) * sin(h) - cos(latitude) * cos(h) * cos(A)
    declination = asin(declination)

    # get ra = sidereal - longitude + H and change sidereal to radians
    sidereal = ln_get_apparent_sidereal_time(JD)
    sidereal *= 2.0 * M_PI / 24.0

    position.ra = ln_range_degrees(ln_rad_to_deg(sidereal - H + longitude))
    position.dec = ln_rad_to_deg(declination)
    return position


# def ln_get_equ_from_ecl( object, JD ) -> ln_equ_posn
# param object Object coordinates.
# param JD Julian day
# returns position Pointer to store new position.

# Transform an objects ecliptical coordinates into equatorial coordinates for the given julian day.

def ln_get_equ_from_ecl( object : ln_lnlat_posn, JD : float ) -> ln_equ_posn:
    position = ln_equ_posn()

    # get obliquity of ecliptic and change it to rads
    nutation = ln_get_nutation(JD)
    nutation.ecliptic = ln_deg_to_rad(nutation.ecliptic)

    # change object's position into radians

    # object
    longitude = ln_deg_to_rad(object.lng)
    latitude = ln_deg_to_rad(object.lat)

    # Equ 12.3, 12.4
    ra = atan2((sin(longitude) * cos(nutation.ecliptic) - tan(latitude) * sin(nutation.ecliptic)), cos(longitude))
    declination = sin(latitude) * cos(nutation.ecliptic) + cos(latitude) * sin(nutation.ecliptic) * sin(longitude)
    declination = asin(declination)

    # store in position
    position.ra = ln_range_degrees(ln_rad_to_deg(ra))
    position.dec = ln_rad_to_deg(declination)
    return position

# def ln_get_ecl_from_equ( object, JD) -> ln_lnlat_posn
# param object Object coordinates in B1950. Use ln_get_equ_prec2 to transform from J2000.
# param JD Julian day
# returns position Pointer to store new position.

# Transform an objects equatorial cordinates into ecliptical coordinates for the given julian day.


def ln_get_ecl_from_equ( object : ln_equ_posn, JD : float ) -> ln_lnlat_posn:
    posn = ln_lnlat_posn()

    # object position
    ra = ln_deg_to_rad(object.ra)
    declination = ln_deg_to_rad(object.dec)
    nutation = ln_get_nutation(JD)
    nutation.ecliptic = ln_deg_to_rad(nutation.ecliptic)

    # Equ 12.1, 12.2
    longitude = atan2((sin(ra) * cos(nutation.ecliptic) + tan(declination) * sin(nutation.ecliptic)), cos(ra))
    latitude = sin(declination) * cos(nutation.ecliptic) - cos(declination) * sin(nutation.ecliptic) * sin(ra)
    latitude = asin(latitude)

    # store in position
    posn.lat = ln_rad_to_deg(latitude)
    posn.lng = ln_range_degrees(ln_rad_to_deg(longitude))
    return posn


# def ln_get_ecl_from_rect( rect ) -> ln_lnlat_posn
# param rect Rectangular coordinates.
# returns posn Pointer to store new position.

# Transform an objects rectangular coordinates into ecliptical coordinates.

def ln_get_ecl_from_rect( rect : ln_rect_posn ) -> ln_lnlat_posn:
    posn = ln_lnlat_posn()
    t = sqrt(rect.X * rect.X + rect.Y * rect.Y)
    posn.lng = ln_range_degrees(ln_rad_to_deg(atan2(rect.X, rect.Y)))
    posn.lat = ln_rad_to_deg(atan2(t, rect.Z))
    return posn

# def ln_get_equ_from_gal( gal ) -> ln_equ_posn
# param gal Galactic coordinates.
# param equ B1950 equatorial coordinates. Use ln_get_equ_prec2 to transform to J2000.

# Transform an object galactic coordinates into B1950 equatorial coordinate.

def ln_get_equ_from_gal( gal : ln_gal_posn ) -> ln_equ_posn:
    equ = ln_equ_posn()

    RAD_27_4 = ln_deg_to_rad(27.4)
    SIN_27_4 = sin(RAD_27_4)
    COS_27_4 = cos(RAD_27_4)

    l_123 = ln_deg_to_rad(gal.l - 123)
    cos_l_123 = cos(l_123)

    rad_gal_b = ln_deg_to_rad(gal.b)

    sin_b = sin(rad_gal_b)
    cos_b = cos(rad_gal_b)

    y = atan2(sin(l_123), cos_l_123 * SIN_27_4 - (sin_b / cos_b) * COS_27_4)
    equ.ra  = ln_range_degrees(ln_rad_to_deg(y) + 12.25)
    equ.dec = ln_rad_to_deg(asin(sin_b * SIN_27_4 + cos_b * COS_27_4 * cos_l_123))
    return equ


# def ln_get_equ2000_from_gal( gal ) -> ln_equ_posn
# param gal Galactic coordinates.
# returns: equ J2000 equatorial coordinates.

# Transform an object galactic coordinates into equatorial coordinate.

def ln_get_equ2000_from_gal( gal : ln_gal_posn ) -> ln_equ_posn:
    equ = ln_equ_posn()
    equ = ln_get_equ_from_gal(gal)
    equ = ln_get_equ_prec2(equ, B1950, JD2000)
    return equ


# def ln_get_gal_from_equ( equ ) -> ln_gal_posn
# param equ B1950 equatorial coordinates.
# returns gal Galactic coordinates.

# Transform an object B1950 equatorial coordinate into galactic coordinates.

def ln_get_gal_from_equ( equ : ln_equ_posn ) -> ln_gal_posn:
    gal = ln_gal_posn()

    RAD_27_4 = ln_deg_to_rad(27.4)
    SIN_27_4 = sin(RAD_27_4)
    COS_27_4 = cos(RAD_27_4)

    ra_192_25 = ln_deg_to_rad(192.25 - equ.ra)
    cos_ra_192_25 = cos(ra_192_25)

    rad_equ_dec = ln_deg_to_rad(equ.dec)

    sin_dec = sin(rad_equ_dec)
    cos_dec = cos(rad_equ_dec)

    x = atan2(sin(ra_192_25), cos_ra_192_25 * SIN_27_4 - (sin_dec / cos_dec) * COS_27_4)
    gal.l = ln_range_degrees(303 - ln_rad_to_deg(x))
    gal.b = ln_rad_to_deg(asin(sin_dec * SIN_27_4 + cos_dec * COS_27_4 * cos_ra_192_25))
    return gal


# def ln_get_gal_from_equ2000( equ ) -> ln_gal_posn
# param equ J2000 equatorial coordinates.
# returns gal Galactic coordinates.

# Transform an object J2000 equatorial coordinate into galactic coordinates.

def ln_get_gal_from_equ2000( equ : ln_equ_posn ) -> ln_gal_posn:
    gal = ln_gal_posn()

    equ_1950 = ln_get_equ_prec2(equ, JD2000, B1950)
    gal = ln_get_gal_from_equ(equ_1950)
    return gal

