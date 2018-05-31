
from libNova.ln_types import *
from libNova.utility import *

# def ln_get_equ_pm( mean_position, proper_motion, JD ) -> ln_equ_posn
# param mean_position Mean position of object.
# param proper_motion Annual Proper motion of object.
# param JD Julian Day.
# returns position Pointer to store new object position.

# Calculate a stars equatorial coordinates from it's mean coordinates (J2000.0) with the effects of proper motion for a given Julian Day.

def ln_get_equ_pm( mean_position : ln_equ_posn, proper_motion : ln_equ_posn, JD : float ) -> ln_equ_posn:
    posn = ln_equ_posn()
    posn = ln_get_equ_pm_epoch (mean_position, proper_motion, JD, JD2000)
    return posn


# def ln_get_equ_pm_epoch( mean_position, proper_motion, JD, epoch_JD) -> ln_equ_posn
# param mean_position Mean position of object.
# param proper_motion Annual Proper motion of object.
# param JD Julian Day.
# param JD_epoch Mean position epoch in JD
# returns position Pointer to store new object position.

# Calculate a stars equatorial coordinates from it's mean coordinates and epoch with the effects of proper motion for a given Julian Day.

def  ln_get_equ_pm_epoch( mean_position : ln_equ_posn, proper_motion : ln_equ_posn, JD : float, epoch_JD : float ) -> ln_equ_posn:
    posn = ln_equ_posn()
    T = (JD - epoch_JD) / 365.25

    # calc proper motion
    posn.ra = mean_position.ra + T * proper_motion.ra
    posn.dec = mean_position.dec + T * proper_motion.dec

    # change to degrees
    posn.ra = ln_range_degrees(posn.ra)
    return posn

