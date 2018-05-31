from math import atan, fabs, tan, sqrt, cos, atan2, sin
from libNova.utility import cbrt
from libNova.ln_types import *
from libNova.transform import *
from libNova.earth import ln_get_earth_helio_coords, ln_get_earth_solar_dist
from libNova.solar import ln_get_solar_geo_coords

LN_STAR_STANDART_HORIZON = -0.5667

GAUS_GRAV =	0.01720209895	# Gaussian gravitational constant k
PREC =  1e-10

# double ln_solve_hyp_barker (double Q1, double G, double t)
# param Q1 See 35.0
# param G See 35.0
# param t Time since perihelion in days
# return Solution of Barkers equation

# Solve Barkers equation. LIAM add more

def ln_solve_hyp_barker(Q1 : float, G : float, t : float) -> float:
	hasmore1 = True
	Q2 = Q1 * t
	S = 2 / (3 * fabs(Q2))
	S = 2 / tan(2 * atan(cbrt (tan(atan(S) / 2))))

	if t< 0:
		S = -S

	L = 0

	# we have initial s, so now do the iteration
	while hasmore1:
		S0 = S
		Z = 1
		Y = S * S
		G1 = -Y * S
		Q3 = Q2 + 2.0 * G * S * Y / 3.0

#		next_z:
		Z+=1
		G1 = -G1 * G * Y
		Z1 = (Z - (Z + 1) * G) / (2.0 * Z + 1.0)
		F = Z1 * G1
		Q3 = Q3 + F

		if Z > 100 or fabs(F) > 10000:
			return nan("0")

		#if fabs(F) > PREC:
		#		goto next_z

		L+=1
		if L > 100:
			return nan("0")

		while hasmore2:
			S1 = S
			S = (2 * S * S * S / 3 + Q3) / (S * S + 1)
			hasmore2 = (fabs(S - S1) > PREC)
	hasmore1 = (fabs(S - S0) > PREC)

	return S

# def ln_get_hyp_true_anomaly( q, e, t)
# param q Perihelion distance in AU
# param e Orbit eccentricity
# param t Time since perihelion
# return True anomaly (degrees)
#
# Calculate the true anomaly.

def ln_get_hyp_true_anomaly( q : float, e : float, t : float ) -> float:
	Q = (GAUS_GRAV / (2.0 * q)) * sqrt((1.0 + e) / q)
	gama = (1.0 - e) / (1.0 + e)
	
	s = ln_solve_hyp_barker(Q, gama, t)
	v = 2.0 * atan(s)

	return ln_range_degrees(ln_rad_to_deg(v))


# def ln_get_hyp_radius_vector( q, e, t) -> float
# param q Perihelion distance in AU
# param e Orbit eccentricity
# param t Time since perihelion in days
# return Radius vector AU
#
# Calculate the radius vector.

def ln_get_hyp_radius_vector( q : float, e : float, t : float ) -> float:
	return q * (1.0 + e) / (1.0 + e * cos(ln_deg_to_rad(ln_get_hyp_true_anomaly(q, e, t))))


# def ln_get_hyp_helio_rect_posn( orbit, JD) -> ln_rect_posn
# param orbit Orbital parameters of object.
# param JD Julian day
# param posn Position pointer to store objects position
#
# Calculate the objects rectangular heliocentric position given it's orbital elements for the given julian day. 

def ln_get_hyp_helio_rect_posn (orbit : ln_hyp_orbit, JD : float ) -> ln_rect_posn:
	posn = ln_rect_posn()
	
	# time since perihelion 
	t = JD - orbit.JD

	# J2000 obliquity of the ecliptic 
	sin_e = 0.397777156
	cos_e = 0.917482062
	
	# equ 33.7 */
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
	
	# equ 33.8 */
	A = atan2(F, P)
	B = atan2(G, Q)
	C = atan2(H, R)
	a = sqrt(F * F + P * P)
	b = sqrt(G * G + Q * Q)
	c = sqrt(H * H + R * R)

	# get true anomaly */
	v = ln_get_hyp_true_anomaly(orbit.q, orbit.e, t)
	
	# get radius vector */
	r = ln_get_hyp_radius_vector(orbit.q, orbit.e, t)

	# equ 33.9 */
	posn.X = r * a * sin(A + ln_deg_to_rad(orbit.w + v))
	posn.Y = r * b * sin(B + ln_deg_to_rad(orbit.w + v))
	posn.Z = r * c * sin(C + ln_deg_to_rad(orbit.w + v))
	return posn

# def ln_get_hyp_geo_rect_posn( orbit, JD) -> ln_rect_posn
# param orbit Orbital parameters of object.
# param JD Julian day
# param posn Position pointer to store objects position
#
# Calculate the objects rectangular geocentric position given it's orbital
# elements for the given julian day.

def  ln_get_hyp_geo_rect_posn( orbit : ln_hyp_orbit, JD : float ) -> ln_rect_posn:
	posn = ln_rect_posn()

	# parabolic helio rect coords
	p_posn = ln_get_hyp_helio_rect_posn(orbit, JD)

	# earth rect coords */
	earth = ln_get_earth_helio_coords(JD)

	e_posn = ln_get_rect_from_helio(earth)
	posn.X = p_posn.X - e_posn.X
	posn.Y = p_posn.Y - e_posn.Y
	posn.Z = p_posn.Z - e_posn.Z
	return posn


# def ln_get_hyp_body_equ_coords( JD, orbit) -> ln_equ_posn
# param JD Julian Day.
# param orbit Orbital parameters.
# param posn Pointer to hold asteroid position.
#
# Calculate a bodies equatorial coordinates for the given julian day.

def ln_get_hyp_body_equ_coords( JD : float, orbit : ln_hyp_orbit ) -> ln_equ_posn:
	posn = ln_equ_posn()

	# get solar and body rect coords */
	body_rect_posn = ln_get_hyp_helio_rect_posn(orbit, JD)
	sol_rect_posn = ln_get_solar_geo_coords(JD)

	# calc distance and light time */
	dist = ln_get_rect_distance(body_rect_posn, sol_rect_posn)
	t = ln_get_light_time(dist)
	
	# repeat calculation with new time (i.e. JD - t)
	body_rect_posn = ln_get_hyp_helio_rect_posn(orbit, JD - t)
	
	# calc equ coords equ 33.10
	x = sol_rect_posn.X + body_rect_posn.X
	y = sol_rect_posn.Y + body_rect_posn.Y
	z = sol_rect_posn.Z + body_rect_posn.Z

	posn.ra = ln_range_degrees(ln_rad_to_deg(atan2(y, x)))
	posn.dec = ln_rad_to_deg(atan2(z,sqrt(x * x + y * y)))
	return posn


# def ln_get_hyp_body_earth_dist( JD, orbit) -> float
# param JD Julian day.
# param orbit Orbital parameters
# returns Distance in AU
#
# Calculate the distance between a body and the Earth
# for the given julian day.

def ln_get_hyp_body_earth_dist( JD : float, orbit : ln_hyp_orbit ) -> float:
	earth_rect_posn = ln_rect_posn()

	# get solar and body rect coords
	body_rect_posn = ln_get_hyp_geo_rect_posn(orbit, JD)
	earth_rect_posn.X = 0.0
	earth_rect_posn.Y = 0.0
	earth_rect_posn.Z = 0.0
	
	# calc distance */
	return ln_get_rect_distance(body_rect_posn, earth_rect_posn)



# def ln_get_hyp_body_solar_dist( JD, ln_hyp_orbit) -> float
# param JD Julian Day.
# param orbit Orbital parameters
#return The distance in AU between the Sun and the body.
#
# Calculate the distance between a body and the Sun.

def ln_get_hyp_body_solar_dist( JD : float, orbit : ln_hyp_orbit ) -> float:
	sol_rect_posn = ln_rect_posn()
	
	# get solar and body rect coords */
	body_rect_posn = ln_get_hyp_helio_rect_posn (orbit, JD)
	sol_rect_posn.X = 0.0
	sol_rect_posn.Y = 0.0
	sol_rect_posn.Z = 0.0
	return ln_get_rect_distance(body_rect_posn, sol_rect_posn)


# def ln_get_hyp_body_phase_angle( JD, ln_hyp_orbit) -> float
# param JD Julian day
# param orbit Orbital parameters
# return Phase angle.
#
# Calculate the phase angle of the body. The angle Sun - body - Earth.

def ln_get_hyp_body_phase_angle( JD : float, orbit : ln_hyp_orbit ) -> float:
	
	# time since perihelion
	t = JD - orbit.JD
	
	# get radius vector
	r = ln_get_hyp_radius_vector(orbit.q, orbit.e, t)
	
	# get solar and Earth-Sun distances
	R = ln_get_earth_solar_dist(JD)
	d = ln_get_hyp_body_solar_dist(JD, orbit)

	phase = (r * r + d * d - R * R) / (2.0 * r * d)
	return ln_range_degrees(ln_rad_to_deg(acos(phase)))


# def ln_get_hyp_body_elong( JD, orbit) -> float
# param JD Julian day
# param orbit Orbital parameters
# return Elongation to the Sun.
#
# Calculate the bodies elongation to the Sun..

def ln_get_hyp_body_elong( JD : float, orbit : ln_hyp_orbit ) -> float:

	# time since perihelion
	t = JD - orbit.JD
	
	# get radius vector */
	r = ln_get_hyp_radius_vector(orbit.q, orbit.e, t)
	
	# get solar and Earth-Sun distances */
	R = ln_get_earth_solar_dist(JD)
	d = ln_get_hyp_body_solar_dist(JD, orbit)

	elong = (R * R + d * d - r * r) / ( 2.0 * R * d )
	return ln_range_degrees(ln_rad_to_deg(acos(elong)))


# def ln_get_hyp_body_rst( JD, observer, orbit) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param rst Pointer to store Rise, Set and Transit time in JD
# return 0 for success, else 1 for circumpolar.
#
# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with a parabolic orbit for the given Julian day.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above or below the horizon. Returns -1 when it remains whole day below the horizon.

def ln_get_hyp_body_rst( JD : float, observer : ln_lnlat_posn, orbit : ln_hyp_orbit ) -> (int, ln_rst_time):
	return ln_get_hyp_body_rst_horizon(JD, observer, orbit, LN_STAR_STANDART_HORIZON)


# def ln_get_hyp_body_rst_horizon( JD, observer, orbit, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param rst Pointer to store Rise, Set and Transit time in JD
# return 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with a parabolic orbit for the given Julian day.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above or below the horizon. Returns -1 when it remains whole day below the horizon.

def ln_get_hyp_body_rst_horizon( JD : float, observer : ln_lnlat_posn, orbit : ln_hyp_orbit, horizon : float ) -> (int, ln_rst_time):
	return ln_get_motion_body_rst_horizon(JD, observer, ln_get_hyp_body_equ_coords, orbit, horizon)


# def ln_get_hyp_body_next_rst(JD, observer, orbit) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param rst Pointer to store Rise, Set and Transit time in JD
# return 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an hyperbolic orbit for the given Julian day.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_hyp_body_next_rst( JD : float, observer : ln_lnlat_posn, orbit : ln_hyp_orbit ) -> (int, ln_rst_time):
	return ln_get_hyp_body_next_rst_horizon(JD, observer, orbit, LN_STAR_STANDART_HORIZON)


# def ln_get_hyp_body_next_rst_horizon( JD, observer, orbit, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param rst Pointer to store Rise, Set and Transit time in JD
#return 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an hyperbolic orbit for the given Julian day.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_hyp_body_next_rst_horizon( JD : float, observer : ln_lnlat_posn, orbit : ln_hyp_orbit, horizon : float ) -> (int, ln_rst_time):
	return ln_get_motion_body_next_rst_horizon(JD, observer, ln_get_hyp_body_equ_coords, orbit, horizon)


# def ln_get_hyp_body_next_rst_horizon_future( JD, observer, orbit, horizon, day_limit) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param day_limit Maximal number of days that will be searched for next rise and set
# param rst Pointer to store Rise, Set and Transit time in JD
# return 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an hyperbolic orbit for the given Julian day.
#
# This function guarantee, that rise, set and transit will be in <JD, JD + day_limit> range.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_hyp_body_next_rst_horizon_future( JD : float, observer : ln_lnlat_posn, orbit : ln_hyp_orbit, horizon : float, day_limit : int ) -> (int, ln_rst_time):
	return ln_get_motion_body_next_rst_horizon_future(JD, observer, ln_get_hyp_body_equ_coords, orbit, horizon, day_limit)

