
from libNova.ln_types import *
from libNova.utility import *
from libNova.nutation import *
from libNova.earth import *
from math import sin, cos

# def ln_get_heliocentric_time_diff(JD, object) -> float
# param JD: float Julian day
# returns object Pointer to object (RA, DEC) for which heliocentric correction will be caculated
# return Heliocentric correction in fraction of day
# Calculate heliocentric corection for object at given coordinates and on given date.

def ln_get_heliocentric_time_diff( JD : float, object : ln_equ_posn ) -> float:
    nutation = ln_get_nutation(JD)
    earth = ln_get_earth_helio_coords(JD)

    theta = ln_deg_to_rad(ln_range_degrees(earth.L + 180))
    ra = ln_deg_to_rad(object.ra)
    dec = ln_deg_to_rad(object.dec)
    c_dec = cos(dec)
    obliq = ln_deg_to_rad(nutation.ecliptic)
    return -0.0057755 * earth.R * (cos(theta) * cos(ra) * c_dec + sin(theta) * (sin(obliq) * sin(dec) + cos(obliq) * c_dec * sin(ra)))

