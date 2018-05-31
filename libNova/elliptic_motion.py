from libNova.ln_types import *
from libNova.rise_set import ln_get_motion_body_rst_horizon
from libNova.utility import *
from libNova.solar import ln_get_solar_geo_coords
from libNova.transform import ln_get_rect_from_helio
from libNova.earth import ln_get_earth_solar_dist, ln_get_earth_helio_coords
from math import cos, sin, atan, atan2, acos, asin


# number of steps in calculation, 3.32 steps for each significant digit required
KEPLER_STEPS = 53

# the BASIC SGN() function  for doubles
def sgn( x : float ) -> float:
	if x == 0.0:
		return x
	else:
		if x < 0.0:
			return -1.0
		else:
			return 1.0


# def ln_solve_kepler( E, M ) -> float
# param E Orbital eccentricity
# param M Mean anomaly
# returns: Eccentric anomaly

# Calculate the eccentric anomaly.
# This method was devised by Roger Sinnott. (Sky and Telescope, Vol 70, pg 159)

def ln_solve_kepler(e : float, M : float) -> float:
	Eo = M_PI_2
	D = M_PI_4

	# covert to radians
	M = ln_deg_to_rad(M)
	
	F = sgn(M)
	M = fabs(M) / (2.0 * M_PI)
	M = (M - int(M)) * 2.0 * M_PI * F
	
	if M < 0:
		M = M + 2.0 * M_PI
	F = 1.0
	
	if M > M_PI:
		F = -1.0
	
	if M > M_PI:
		M = 2.0 * M_PI - M
	
	for i in range(0, KEPLER_STEPS ):
		M1 = Eo - e * sin(Eo)
		Eo = Eo + D * sgn(M - M1)
		D /= 2.0
	Eo *= F
	
	# back to degrees
	Eo = ln_rad_to_deg(Eo)
	return Eo


# def ln_get_ell_mean_anomaly( n, delta_JD ) -> double
# param n Mean motion (degrees/day)
# param delta_JD Time since perihelion
# returns Mean anomaly (degrees)

# Calculate the mean anomaly.

def ln_get_ell_mean_anomaly(n : float, delta_JD : float) -> float:
	return delta_JD * n


# def ln_get_ell_true_anomaly( e, E ) -> float
# param e Orbital eccentricity
# param E Eccentric anomaly
# returns True anomaly (degrees)

# Calculate the true anomaly.

def ln_get_ell_true_anomaly( e : float, E : float) -> float:
	E = ln_deg_to_rad(E)
	v = sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0)
	v = 2.0 * atan(v)
	v = ln_range_degrees(ln_rad_to_deg(v))
	return v


# double ln_get_ell_radius_vector( a, e, E ) -> float
# param a Semi-Major axis in AU
# param e Orbital eccentricity
# param E Eccentric anomaly
# returns Radius vector AU

# Calculate the radius vector.

def ln_get_ell_radius_vector(a : float, e : float, E : float) -> float:
	return a * (1.0 - e * cos(ln_deg_to_rad(E)))



# def ln_get_ell_smajor_diam( e, q ) > float
# param e Eccentricity
# param q Perihelion distance in AU
# returns Semi-major diameter in AU

# Calculate the semi major diameter.

def ln_get_ell_smajor_diam(e : float, q : float) -> float:
	return q / (1.0 - e)


# double ln_get_ell_sminor_diam( e, a ) -> float
# param e: float		 Eccentricity
# param a: float		 Semi-Major diameter in AU
# return Semi-minor diameter in AU

# Calculate the semi minor diameter.

def ln_get_ell_sminor_diam(e : float, a : float) -> float:
	return a * sqrt(1 - e * e)


# def ln_get_ell_mean_motion( a ) -> float
# param a: float		 Semi major diameter in AU
# returns Mean daily motion (degrees/day)

# Calculate the mean daily motion (degrees/day).

def ln_get_ell_mean_motion(a : float) -> float:
	q = 0.9856076686                # Gaussian gravitational constant (degrees)
	return q / (a * sqrt(a))


# def ln_get_ell_helio_rect_posn( orbit, JD ) -> ln_rect_posn
# param orbit: ln_ell_orbit		 Orbital parameters of object.
# param JD: float		 		Julian day
# returns posn Position pointer to store objects position

# Calculate the objects rectangular heliocentric position given it's orbital elements for the given julian day.

def ln_get_ell_helio_rect_posn(orbit : ln_ell_orbit, JD : float) -> ln_rect_posn:
	posn = ln_equ_posn()

	# J2000 obliquity of the ecliptic
	sin_e = 0.397777156
	cos_e = 0.917482062

	# equ 33.7
	sin_omega = sin(ln_deg_to_rad(orbit.omega))
	cos_omega = cos(ln_deg_to_rad(orbit.omega))
	sin_i = sin(ln_deg_to_rad(orbit.i))
	cos_i = cos(ln_deg_to_rad(orbit.i))
	F = cos_omega
	G = sin_omega * cos_e
	H = sin_omega * sin_e
	P = -sin_omega * cos_i
	Q = cos_omega * cos_i * cos_e - sin_i * sin_e
	R = cos_omega * cos_i * sin_e + sin_i * cos_e

	# equ 33.8
	A = atan2(F,P)
	B = atan2(G,Q)
	C = atan2(H,R)
	a = sqrt(F ** 2 + P ** 2)
	b = sqrt(G ** 2 + Q ** 2)
	c = sqrt(H ** 2 + R ** 2)

	# get daily motion
	if orbit.n == 0:
		orbit.n = ln_get_ell_mean_motion(orbit.a)

	# get mean anomaly
	M = ln_get_ell_mean_anomaly(orbit.n, JD - orbit.JD)

	# get eccentric anomaly
	E = ln_solve_kepler(orbit.e, M)

	# get true anomaly
	v = ln_get_ell_true_anomaly(orbit.e, E)

	# get radius vector
	r = ln_get_ell_radius_vector(orbit.a, orbit.e, E)

	# equ 33.9
	posn.X = r * a * sin(A + ln_deg_to_rad(orbit.w + v))
	posn.Y = r * b * sin(B + ln_deg_to_rad(orbit.w + v))
	posn.Z = r * c * sin(C + ln_deg_to_rad(orbit.w + v))
	return posn

# def ln_get_ell_geo_rect_posn( orbit, JD) -> ln_rect_posn
# param orbit: ln_ell_orbit		 Orbital parameters of object.
# param JD:float				 Julian day
# returns ln_rect_posn

# Calculate the objects rectangular geocentric position given it's orbital elements for the given julian day.

def ln_get_ell_geo_rect_posn( orbit : ln_ell_orbit, JD : float ) -> ln_rect_posn:
	posn = ln_rect_posn()

	# elliptic helio rect coords
	p_posn = ln_get_ell_helio_rect_posn(orbit, JD)

	# earth rect coords
	earth  =  ln_get_earth_helio_coords(JD)
	e_posn =  ln_get_rect_from_helio(earth)

	posn.X = e_posn.X - p_posn.X
	posn.Y = e_posn.Y - p_posn.Y
	posn.Z = e_posn.Z - p_posn.Z
	return posn


# def ln_get_ell_body_equ_coords(JD, orbit) -> ln_equ_posn
# param JD: float			 	Julian Day.
# param orbit: ln_ell_orbit		 Orbital parameters.
# returns						equatorial position
# Calculate a bodies equatorial coordinates for the given julian day.

def ln_get_ell_body_equ_coords( JD : float, orbit : ln_ell_orbit ) -> ln_equ_posn:
	posn = ln_equ_posn()

	# get solar and body rect coords
	body_rect_posn = ln_get_ell_helio_rect_posn(orbit, JD)
	sol_rect_posn = ln_get_solar_geo_coords(JD)

	# calc distance and light time
	dist = ln_get_rect_distance(body_rect_posn, sol_rect_posn)
	t = ln_get_light_time(dist)

	# repeat calculation with new time (i.e. JD - t)
	body_rect_posn = ln_get_ell_helio_rect_posn(orbit, JD - t)

	# calc equ coords equ 33.10
	x = sol_rect_posn.X + body_rect_posn.X
	y = sol_rect_posn.Y + body_rect_posn.Y
	z = sol_rect_posn.Z + body_rect_posn.Z

	posn.ra = ln_range_degrees(ln_rad_to_deg(atan2(y,x)))
	posn.dec = ln_rad_to_deg(asin(z / sqrt(x ** 2 + y ** 2 + z **2)))
	return posn


# def ln_get_ell_orbit_len( orbit) -> float
# param orbit : ln_ell_orbit	 Orbital parameters
# return Orbital length in AU(float)

# Calculate the orbital length in AU.
#
# Accuracy:
# - 0.001% for e < 0.88
# - 0.01% for e < 0.95
# - 1% for e = 0.9997
# - 3% for e = 1

def ln_get_ell_orbit_len( orbit : ln_ell_orbit ) -> float:
	b = ln_get_ell_sminor_diam(orbit.e, orbit.a)

	A = (orbit.a + b) / 2.0
	G = sqrt(orbit.a * b)
	H = (2.0 * orbit.a * b) / (orbit.a + b)
	return M_PI * ((21.0 * A - 2.0 * G - 3.0 * H) / 8.0)


# def ln_get_ell_orbit_vel( JD,  orbit ) -> float
# param JD: float			 	Julian day.
# param orbit: ln_ell_orbit		 Orbital parameters
# return Orbital velocity in km/s.

# Calculate orbital velocity in km/s for the given julian day.

def ln_get_ell_orbit_vel( JD : float, orbit : ln_ell_orbit ) -> float:

	r = ln_get_ell_body_solar_dist(JD, orbit)
	V = 1.0 / r - 1.0 / (2.0 * orbit.a)
	V = 42.1219 * sqrt(V)
	return V


# def ln_get_ell_orbit_pvel( orbit ) -> float
# param orbit : ln_ell_orbit	 Orbital parameters
# return Orbital velocity in km/s.

# Calculate orbital velocity at perihelion in km/s.

def ln_get_ell_orbit_pvel( orbit : ln_ell_orbit ) -> float:
	V = 29.7847 / sqrt(orbit.a)
	V *= sqrt((1.0 + orbit.e) / (1.0 - orbit.e))
	return V


# def ln_get_ell_orbit_avel( orbit ) -> float
# param orbit : ln_ell_orbit	 Orbital parameters
# return Orbital velocity in km/s.

# Calculate the orbital velocity at aphelion in km/s.

def ln_get_ell_orbit_avel(orbit : ln_ell_orbit) -> float:
	V = 29.7847 / sqrt(orbit.a)
	V *= sqrt((1.0 - orbit.e) / (1.0 + orbit.e))
	return V



# def ln_get_ell_body_solar_dist( JD, orbit ) -> float
# param JD: float				 Julian Day.
# param orbit: ln_ell_orbit		 Orbital parameters
# return The distance in AU between the Sun and the body.

# Calculate the distance between a body and the Sun.

def ln_get_ell_body_solar_dist( JD : float, orbit : ln_ell_orbit ) -> float:
	# get solar and body rect coords
	sol_rect_posn = ln_rect_posn()

	body_rect_posn = ln_get_ell_helio_rect_posn(orbit, JD)
	sol_rect_posn.X = 0
	sol_rect_posn.Y = 0
	sol_rect_posn.Z = 0

	# calc distance
	return ln_get_rect_distance(body_rect_posn, sol_rect_posn)



# def ln_get_ell_body_earth_dist( JD, orbit ) -> float
# param JD: float				 Julian Day.
# param orbit: ln_ell_orbit		 Orbital parameters
# returns Distance in AU

# Calculate the distance between an body and the Earth for the given julian day.

def ln_get_ell_body_earth_dist( JD : float,  orbit : ln_ell_orbit ) -> float:
	earth_rect_posn = ln_rect_posn()

	# get solar and body rect coords
	body_rect_posn = ln_get_ell_geo_rect_posn(orbit, JD)
	earth_rect_posn.X = 0
	earth_rect_posn.Y = 0
	earth_rect_posn.Z = 0
	return ln_get_rect_distance(body_rect_posn, earth_rect_posn)



# def ln_get_ell_body_phase_angle( JD, orbit) -> float
# param JD: float				 Julian Day.
# param orbit: ln_ell_orbit		 Orbital parameters
# return Phase angle.

#  Calculate the phase angle of the body. The angle Sun - body - Earth.

def ln_get_ell_body_phase_angle( JD : float,  orbit : ln_ell_orbit) -> float:
	if orbit.n == 0.0:
		orbit.n = ln_get_ell_mean_motion(orbit.a)
	M = ln_get_ell_mean_anomaly(orbit.n, JD - orbit.JD)
	
	# get eccentric anomaly
	E = ln_solve_kepler(orbit.e, M)
	
	# get radius vector
	r = ln_get_ell_radius_vector(orbit.a, orbit.e, E)
	
	# get solar and Earth distances
	R = ln_get_ell_body_earth_dist(JD, orbit)
	d = ln_get_ell_body_solar_dist(JD, orbit)
	
	phase = (r * r + d * d - R * R) / ( 2.0 * r * d)
	return ln_range_degrees(acos(ln_deg_to_rad(phase)))



# def ln_get_ell_body_elong( JD, orbit ) -> float
# param JD: float				 Julian Day.
# param orbit: ln_ell_orbit		 Orbital parameters
# return Elongation to the Sun.

# Calculate the bodies elongation to the Sun..

def ln_get_ell_body_elong( JD : float, orbit : ln_ell_orbit ) -> float:

	# time since perihelion
	t = JD - orbit.JD
	
	# get mean anomaly
	if orbit.n == 0.0:
		orbit.n = ln_get_ell_mean_motion(orbit.a)
	M = ln_get_ell_mean_anomaly(orbit.n, t)
	
	# get eccentric anomaly
	E = ln_solve_kepler(orbit.e, M)
	
	# get radius vector
	r = ln_get_ell_radius_vector(orbit.a, orbit.e, E)
	
	# get solar and Earth-Sun distances
	R = ln_get_earth_solar_dist(JD)
	d = ln_get_ell_body_solar_dist(JD, orbit)

	elong = (R * R + d * d - r * r) / ( 2.0 * R * d)
	return ln_range_degrees(ln_rad_to_deg(acos(elong)))


# def ln_get_ell_body_rst( JD, observer, orbit) -> (int, rst_time)
# param JD: float					 Julian day
# param observer: ln_lnlat_posn		 Observers position
# param orbit: ln_ell_orbit			 Orbital parameters
# return 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
# returns rst Pointer to store Rise, Set and Transit time in JD


# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an elliptic orbit for the given Julian day.

# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_ell_body_rst( JD : float, observer : ln_lnlat_posn, orbit : ln_ell_orbit) -> (int, ln_rst_time):
	return ln_get_ell_body_rst_horizon(JD, observer, orbit, LN_STAR_STANDART_HORIZON)

# def ln_get_ell_body_rst_horizon(JD : float, observer, orbit, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param rst Pointer to store Rise, Set and Transit time in JD
#return 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an elliptic orbit for the given Julian day.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_ell_body_rst_horizon(JD : float, observer : ln_lnlat_posn, orbit : ln_ell_orbit, horizon : float) -> (int, ln_rst_time):
	return ln_get_motion_body_rst_horizon(JD, observer, ln_get_ell_body_equ_coords, orbit, horizon)


# def ln_get_ell_body_next_rst(JD, observer, orbit) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param rst Pointer to store Rise, Set and Transit time in JD
# return 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an elliptic orbit for the given Julian day.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_ell_body_next_rst(JD : float, observer : ln_lnlat_posn, orbit : ln_ell_orbit) -> (int, ln_rst_time):
	return ln_get_ell_body_next_rst_horizon(JD, observer, orbit, LN_STAR_STANDART_HORIZON)


# def ln_get_ell_body_next_rst_horizon(JD, observer, orbit, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param rst Pointer to store Rise, Set and Transit time in JD
# return 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an elliptic orbit for the given Julian day.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_ell_body_next_rst_horizon(JD : float, observer : ln_lnlat_posn, orbit : ln_ell_orbit,  horizon : float) -> (int, ln_rst_time):
	return ln_get_motion_body_next_rst_horizon(JD, observer, ln_get_ell_body_equ_coords, orbit, horizon)


# def ln_get_ell_body_next_rst_horizon( JD, observer, orbit, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param day_limit Maximal number of days that will be searched for next rise and set
# param rst Pointer to store Rise, Set and Transit time in JD
# return 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an elliptic orbit for the given Julian day.
#
# This function guarantee, that rise, set and transit will be in <JD, JD + day_limit> range.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_ell_body_next_rst_horizon_future(JD : float, observer : ln_lnlat_posn, orbit : ln_ell_orbit, horizon : float, day_limit : int) -> (int, ln_rst_time):
	return ln_get_motion_body_next_rst_horizon_future(JD, observer, ln_get_ell_body_equ_coords, orbit, horizon, day_limit)


# def ln_get_ell_last_perihelion( epoch_JD, M, n ) -> float
# param epoch_JD: float			 Julian day of epoch
# param M: float				 Mean anomaly
# param n: float				 daily motion in degrees
#
# Calculate the julian day of the last perihelion.


def ln_get_ell_last_perihelion(epoch_JD : float, M : float, n : float) -> float:
	return epoch_JD - (M / n)

