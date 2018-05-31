from libNova.ln_types import *
from libNova.utility import *
from libNova.nutation import *
from libNova.earth import *
from libNova.transform import *
from libNova.rise_set import ln_get_body_rst_horizon
from libNova.julian_day import ln_get_date

LN_SOLAR_STANDART_HORIZON = -0.8333
LN_SOLAR_CIVIL_HORIZON = -6.0
LN_SOLAR_NAUTIC_HORIZON = -12.0
LN_SOLAR_ASTRONOMICAL_HORIZON = -18.0

# def ln_get_solar_geom_coords( JD ) -> ln_helio_posn
# param JD Julian day
# param position Pointer to store calculated solar position.

# Calculate geometric coordinates and radius vector accuracy 0.01 arc second error - uses VSOP87 solution.

# Latitude and Longitude returned are in degrees, whilst radius vector returned is in AU.


def ln_get_solar_geom_coords( JD : float ) -> ln_helio_posn:
    position = ln_helio_posn()

    # get earths heliocentric position
    position = ln_get_earth_helio_coords(JD)

    position.L += 180.0
    position.L = ln_range_degrees(position.L)
    position.B *= -1.0
    return position


# def ln_get_solar_equ_coords( JD ) -> ln_equ_posn
# param JD Julian day
# param position Pointer to store calculated solar position.

# Calculate apparent equatorial solar coordinates for given julian day.
# This function includes the effects of aberration and nutation.

def ln_get_solar_equ_coords( JD : float ) -> ln_equ_posn:
    position = ln_equ_posn()
    LB = ln_lnlat_posn()

    # get geometric coords */
    sol = ln_get_solar_geom_coords(JD)

    # add nutation
    nutation = ln_get_nutation(JD)
    sol.L += nutation.longitude

    # aberration
    aberration = (20.4898 / (360.0 * 60.0 * 60.0)) / sol.R
    sol.L -= aberration

    # transform to equatorial
    LB.lat = sol.B
    LB.lng = sol.L
    position = ln_get_equ_from_ecl(LB, JD)
    return position


# def ln_get_solar_ecl_coords( JD ) -> ln_lnlat_posn
# param JD Julian day
# param position Pointer to store calculated solar position.

# Calculate apparent ecliptical solar coordinates for given julian day.
# This function includes the effects of aberration and nutation.

def ln_get_solar_ecl_coords( JD : float ) -> ln_lnlat_posn:
    position = ln_lnlat_posn()

    # get geometric coords
    sol = ln_get_solar_geom_coords(JD)

    # add nutation
    nutation = ln_get_nutation(JD)
    sol.L += nutation.longitude

    # aberration
    aberration = (20.4898 / (360.0 * 60.0 * 60.0)) / sol.R
    sol.L -= aberration

    position.lng = sol.L
    position.lat = sol.B
    return position


# def ln_get_solar_geo_coords( JD ) -> ln_rect_posn:
# param JD Julian day
# returns position Pointer to store calculated solar position.

# Calculate geocentric coordinates (rectangular) for given julian day.
# Accuracy 0.01 arc second error - uses VSOP87 solution. Position returned is in units of AU.

def ln_get_solar_geo_coords( JD : float ) -> ln_rect_posn:
    position = ln_rect_posn()

    # get earths's heliocentric position
    sol = ln_get_earth_helio_coords(JD)

    # now get rectangular coords
    position = ln_get_rect_from_helio(sol)
    position.X *=-1.0
    position.Y *=-1.0
    position.Z *=-1.0
    return position


def ln_get_solar_rst_horizon( JD : float, observer : ln_lnlat_posn, horizon : float ) -> (int, ln_rst_time):
    return ln_get_body_rst_horizon(JD, observer, ln_get_solar_equ_coords, horizon)



# def ln_get_solar_rst( JD, observer) -> (int, ln_rst_time)
# Calls get_solar_rst_horizon with horizon set to LN_SOLAR_STANDART_HORIZON.


def ln_get_solar_rst( JD : float, observer : ln_lnlat_posn ) -> (int, ln_rst_time):
    return ln_get_solar_rst_horizon(JD, observer, LN_SOLAR_STANDART_HORIZON)


# def ln_get_solar_sdiam( JD ) -> double
# param JD Julian day
# return Semidiameter in arc seconds

# Calculate the semidiameter of the Sun in arc seconds for the given julian day.

def ln_get_solar_sdiam( JD : float ) -> float:
    So = 959.63                                 # at 1 AU
    dist = ln_get_earth_solar_dist(JD)
    return So / dist

def ln_get_day_length( JD : float, observer : ln_lnlat_posn ) -> float:
    (cir, rst) = ln_get_solar_rst(JD, observer)
    deltaJ = rst.set - rst.rise
    return deltaJ * 24