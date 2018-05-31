
from libNova.ln_types import *
from libNova.utility import *
from math import cos,sin,atan2


# def ln_get_angular_separation( posn1, posn2 ) -> float
# param posn1:ln_equ_posn		 Equatorial position of body 1
# param posn2:ln_equ_posn		 Equatorial position of body 2
# returns: Angular separation in degrees

# Calculates the angular separation of 2 bodies.
# This method was devised by Mr Thierry Pauwels of the Royal Observatory Belgium.

def ln_get_angular_separation( posn1 : ln_equ_posn, posn2 : ln_equ_posn) -> float:
	
	# covert to radians
	a1 = ln_deg_to_rad(posn1.ra)
	d1 = ln_deg_to_rad(posn1.dec)
	a2 = ln_deg_to_rad(posn2.ra)
	d2 = ln_deg_to_rad(posn2.dec)
	
	x = (cos(d1) * sin(d2)) - (sin(d1) * cos(d2) * cos(a2 - a1))
	y = cos(d2) * sin(a2 - a1)
	z = (sin(d1) * sin(d2)) + (cos(d1) * cos(d2) * cos(a2 - a1))

	x = x * x
	y = y * y
	d = atan2(sqrt(x + y), z)
	
	return ln_rad_to_deg(d)

# def ln_get_rel_posn_angle( posn1, posn2 : ln_equ_posn) -> float
# param posn1 Equatorial position of body 1
# param posn2 Equatorial position of body 2
# returns Position angle in degrees

# Calculates the position angle of a body with respect to another body.

def ln_get_rel_posn_angle(posn1 : ln_equ_posn, posn2 : ln_equ_posn) -> float:
	
	# covert to radians
	a1 = ln_deg_to_rad(posn1.ra)
	d1 = ln_deg_to_rad(posn1.dec)
	a2 = ln_deg_to_rad(posn2.ra)
	d2 = ln_deg_to_rad(posn2.dec)
	
	y = sin(a1 - a2)
	x = (cos(d2) * tan(d1)) - (sin(d2) * cos(a1 - a2))
	
	P = atan2(y, x)
	return ln_rad_to_deg(P)
