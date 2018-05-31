from libNova.ln_types import *
from libNova.utility import *
from math import cos, sin, atan2, acos, asin

# def ln_get_equ_prec( mean_position, JD ) -> ln_equ_posn
# param mean_position Mean object position
# param JD Julian day
# returns position Pointer to store new object position.

# Calculate equatorial coordinates with the effects of precession for a given Julian Day.
# Uses mean equatorial coordinates and is only for initial epoch J2000.0

def ln_get_equ_prec( mean_position : ln_equ_posn, JD : float ) -> ln_equ_posn:
	posn = ln_equ_posn()

	# change original ra and dec to radians
	mean_ra = ln_deg_to_rad(mean_position.ra)
	mean_dec = ln_deg_to_rad(mean_position.dec)

	t = (JD - JD2000) / 36525.0
	t *= 1.0 / 3600.0
	t2 = t * t
	t3 = t2 *t
	zeta = 2306.2181 * t + 0.30188 * t2 + 0.017998 * t3
	eta = 2306.2181 * t + 1.09468 * t2 + 0.041833 * t3
	theta = 2004.3109 * t - 0.42665 * t2 - 0.041833 * t3
	zeta = ln_deg_to_rad(zeta)
	eta = ln_deg_to_rad(eta)
	theta = ln_deg_to_rad(theta)

	A = cos(mean_dec) * sin(mean_ra + zeta)
	B = cos(theta) * cos(mean_dec) * cos(mean_ra + zeta) - sin(theta) * sin(mean_dec)
	C = sin (theta) * cos (mean_dec) * cos(mean_ra + zeta) + cos(theta) * sin(mean_dec)
	
	ra = atan2(A, B) + eta
	
	# check for object near celestial pole
	if mean_dec > (0.4 * M_PI) or mean_dec < (-0.4 * M_PI):
		# close to pole
		dec = acos(sqrt(A * A + B * B))
		if mean_dec < 0.0:
			dec *= -1
	else:
		# not close to pole
		dec = asin(C)


	# change to degrees
	posn.ra = ln_range_degrees(ln_rad_to_deg(ra))
	posn.dec = ln_rad_to_deg(dec)
	return posn

# def ln_get_equ_prec2( mean_position, fromJD, toJD) -> ln_equ_posn
# param mean_position Mean object position
# param fromJD Julian day (start)
# param toJD Julian day (end)
# returns: position Pointer to store new object position.

# Calculate the effects of precession on equatorial coordinates, between arbitary Jxxxx epochs.
# Use fromJD and toJD parameters to specify required Jxxxx epochs.

def ln_get_equ_prec2( mean_position : ln_equ_posn, fromJD : float, toJD : float ) -> ln_equ_posn:
    posn = ln_equ_posn()

    # change original ra and dec to radians */
    mean_ra = ln_deg_to_rad(mean_position.ra)
    mean_dec = ln_deg_to_rad(mean_position.dec)

    # calc t, T, zeta, eta and theta Equ 20.2
    T = float((fromJD - JD2000)) / 36525.0
    T *= 1.0 / 3600.0
    t = (float(toJD - fromJD)) / 36525.0
    t *= 1.0 / 3600.0
    T2 = T * T
    t2 = t * t
    t3 = t2 *t
    zeta = (2306.2181 + 1.39656 * T - 0.000139 * T2) * t + (0.30188 - 0.000344 * T) * t2 + 0.017998 * t3
    eta = (2306.2181 + 1.39656 * T - 0.000139 * T2) * t + (1.09468 + 0.000066 * T) * t2 + 0.018203 * t3
    theta = (2004.3109 - 0.85330 * T - 0.000217 * T2) * t - (0.42665 + 0.000217 * T) * t2 - 0.041833 * t3
    zeta = ln_deg_to_rad(zeta)
    eta = ln_deg_to_rad(eta)
    theta = ln_deg_to_rad(theta)

    A = cos(mean_dec) * sin(mean_ra + zeta);
    B = cos(theta) * cos(mean_dec) * cos(mean_ra + zeta) - sin(theta) * sin(mean_dec)
    C = sin(theta) * cos(mean_dec) * cos(mean_ra + zeta) + cos(theta) * sin(mean_dec)

    ra = atan2(A, B) + eta

    # check for object near celestial pole
    if mean_dec > (0.4 * M_PI) or  mean_dec < (-0.4 * M_PI):
        # close to pole
        dec = acos(sqrt(A * A + B * B))
        if mean_dec < 0.0:
            dec *= -1
    else:
        # not close to pole
        dec = asin(C)


    # change to degrees
    posn.ra = ln_range_degrees(ln_rad_to_deg(ra))
    posn.dec = ln_rad_to_deg(dec)
    return posn

# def ln_get_ecl_prec( mean_position, JD ) -> ln_lnlat_posn
# param mean_position Mean object position
# param JD Julian day
# returns position Pointer to store new object position.

# Calculate ecliptical coordinates with the effects of precession for a given Julian Day.
# Uses mean ecliptical coordinates and is only for initial epoch J2000.0

def ln_get_ecl_prec( mean_position : ln_lnlat_posn, JD : float ) -> ln_lnlat_posn:
    pass
