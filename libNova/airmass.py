

from math import sin, sqrt,asin
from libNova.ln_types import *
from libNova.utility import *

# def ln_get_airmass( alt, airmass_scale ) -> float
# param alt: float			            Altitude in degrees
# param airmass_scale: float		    Airmass scale - usually 750.
# returns:  Airmass:float for give altitude.

def ln_get_airmass(alt : float, airmass_scale : float) -> float:
	a = airmass_scale * sin( ln_deg_to_rad(alt) )
	return sqrt(a * a + 2 * airmass_scale + 1) - a


# def ln_get_alt_from_airmass( X : float, airmass_scale : float) -> float
# param X: float		              Airmass
# param airmass_scale: float		  Airmass scale - usually 750.
# returns  Altitude:float for give airmass.

def ln_get_alt_from_airmass( X : float, airmass_scale : float) -> float:
	return ln_rad_to_deg(asin((2 * airmass_scale + 1 - X * X) / (2 * X * airmass_scale)))
