
from libNova.ln_types import *
from libNova.aberration import *
from libNova.proper_motion import *
from libNova.precession import *
from libNova.nutation import *

#  Apparent place of an Object


# def ln_get_apparent_posn( mean_position, proper_motion, JD ) -> ln_equ_posn
# param mean_position Mean position of object
# param proper_motion Proper motion of object
# param JD Julian Day
# returns: position Pointer to store new object position

# Calculate the apparent equatorial position of a star from its mean equatorial position.
# This function takes into account the effects of proper motion, precession, nutation,
# annual aberration when calculating the stars apparent position. The effects of annual
# parallax and the gravitational deflection of light (Einstein effect) are NOT used in this calculation.

def ln_get_apparent_posn( mean_position : ln_equ_posn, proper_motion : ln_equ_posn, JD : float ) -> ln_equ_posn:
    proper_position = ln_get_equ_pm( mean_position, proper_motion, JD )
    aberration_position = ln_get_equ_aber( proper_position, JD )
    precession_position = ln_get_equ_prec( aberration_position, JD )
    return ln_get_equ_nut( precession_position, JD )

