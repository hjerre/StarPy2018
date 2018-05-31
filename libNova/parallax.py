from libNova.ln_types import *
from libNova.utility import *
from libNova.sidereal_time import *
from libNova.earth import *
from math import atan2


# def ln_get_parallax( object, au_distance, observer, height, JD) -> ln_equ_posn
# param object Object geocentric coordinates
# param au_distance Distance of object from Earth in AU
# param observer Geographics observer positions
# param height Observer height in m
# param JD  Julian day of observation
# param parallax RA and DEC parallax

# Calculate body parallax, which is need to calculate topocentric position of the body.


def ln_get_parallax( object : ln_equ_posn, au_distance : float, observer : ln_lnlat_posn, height : float, JD : float ) -> ln_equ_posn:
	parallax = ln_equ_posn()

	H = ln_get_apparent_sidereal_time(JD) + (observer.lng - object.ra) / 15.0
	parallax = ln_get_parallax_ha(object, au_distance, observer, height, H)
	return parallax



# def ln_get_parallax_ha( object, au_distance, observer, height, H) -> ln_equ_posn
# param object Object geocentric coordinates
# param au_distance Distance of object from Earth in AU
# param observer Geographics observer positions
# param height Observer height in m
# param H Hour angle of object in hours
# returns parallax RA and DEC parallax

# Calculate body parallax, which is need to calculate topocentric position of the body.
# Uses hour angle as time reference (handy in case we already compute it).


def ln_get_parallax_ha( object : ln_equ_posn, au_distance : float, observer : ln_lnlat_posn, height : float, H : float ) -> ln_equ_posn:
    parallax = ln_equ_posn()

    (ro_sin, ro_cos) = ln_get_earth_centre_dist (height, observer.lat)
    sin_pi = sin(ln_deg_to_rad((8.794 / au_distance) / 3600.0))

    # change hour angle from hours to radians
    H *= M_PI / 12.0

    sin_H = sin(H)
    cos_H = cos(H)

    dec_rad = ln_deg_to_rad(object.dec)
    cos_dec = cos(dec_rad)

    parallax.ra = atan2(-ro_cos * sin_pi * sin_H, cos_dec  - ro_cos * sin_pi * cos_H)
    parallax.dec = atan2((sin(dec_rad) - ro_sin * sin_pi) * cos(parallax.ra), cos_dec - ro_cos * sin_pi * cos_H)

    parallax.ra = ln_rad_to_deg(parallax.ra)
    parallax.dec = ln_rad_to_deg(parallax.dec) - object.dec
    return parallax

