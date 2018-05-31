
from libNova.ln_types import *
from libNova.utility import *
from libNova.solar import ln_get_solar_geo_coords
from libNova.earth import ln_get_earth_helio_coords, ln_get_earth_solar_dist

from math import cos, sin, acos, asin, atan2, sqrt, pi, atan

# def ln_solve_barker( q, t) -> float
# param q Perihelion distance in AU
# param t Time since perihelion in days
# returns Solution of Barkers equation

# Solve Barkers equation. LIAM add more 

def ln_solve_barker( q : float, t : float ) -> float:
	W = ((0.03649116245) / (q * sqrt(q))) * t

	G = W / 2.0
	Y = cbrt(G + sqrt(G * G + 1.0))
	return Y - 1.0 / Y

# def ln_get_par_true_anomaly( q, t) -> float
# param q Perihelion distance in AU
# param t Time since perihelion
# returns True anomaly (degrees)

# Calculate the true anomaly. 


def ln_get_par_true_anomaly( q : float, t : float ) -> float:
	s = ln_solve_barker(q, t)
	v = 2.0 * atan(s)
	
	return ln_range_degrees(ln_rad_to_deg(v))

# def ln_get_par_radius_vector( q, t) -> float
# param q Perihelion distance in AU
# param t Time since perihelion in days
# returns Radius vector AU
#
# Calculate the radius vector. 

def ln_get_par_radius_vector( q : float, t : float ) -> float:
	s = ln_solve_barker(q, t)
	return q * (1.0 + s * s)



# def ln_get_par_helio_rect_posn( orbit, JD) -> ln_rect_posn
# param orbit Orbital parameters of object.
# param JD Julian day
# param posn Position pointer to store objects position
#
# Calculate the objects rectangular heliocentric position given it's orbital elements for the given julian day.

def ln_get_par_helio_rect_posn( orbit : ln_par_orbit, JD : float ) -> ln_rect_posn:
	posn = ln_rect_posn()

	# time since perihelion
	t = JD - orbit.JD
	
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
	A = atan2(F, P)
	B = atan2(G, Q)
	C = atan2(H, R)
	a = sqrt(F * F + P * P)
	b = sqrt(G * G + Q * Q)
	c = sqrt(H * H + R * R)
	
	# get true anomaly 
	v = ln_get_par_true_anomaly(orbit.q, t)
	
	# get radius vector 
	r = ln_get_par_radius_vector(orbit.q, t)

	# equ 33.9 
	posn.X = r * a * sin(A + ln_deg_to_rad(orbit.w + v))
	posn.Y = r * b * sin(B + ln_deg_to_rad(orbit.w + v))
	posn.Z = r * c * sin(C + ln_deg_to_rad(orbit.w + v))
	return posn


# def ln_get_par_geo_rect_posn( orbit,  JD) -> ln_rect_posn
# param orbit Orbital parameters of object.
# param JD Julian day
# param posn Position pointer to store objects position
#
# Calculate the objects rectangular geocentric position given it's orbital elements for the given julian day.

def ln_get_par_geo_rect_posn( orbit : ln_par_orbit, JD : float ) -> ln_rect_posn:
	posn = ln_rect_posn()
	
	# parabolic helio rect coords 
	p_posn = ln_get_par_helio_rect_posn(orbit, JD)
	
	# earth rect coords 
	earth = ln_get_earth_helio_coords(JD)
	
	e_posn = ln_get_rect_from_helio(earth)
	posn.X = p_posn.X - e_posn.X
	posn.Y = p_posn.Y - e_posn.Y
	posn.Z = p_posn.Z - e_posn.Z
	return posn



# def ln_get_par_body_equ_coords( JD, orbit) -> ln_equ_posn
# param JD Julian Day.
# param orbit Orbital parameters.
# param posn Pointer to hold asteroid position.
#
# Calculate a bodies equatorial coordinates for the given julian day.

def ln_get_par_body_equ_coords( JD : float, orbit : ln_par_orbit ) -> ln_equ_posn:
	posn = ln_equ_posn()

	# get solar and body rect coords 
	body_rect_posn = ln_get_par_helio_rect_posn (orbit, JD)
	sol_rect_posn = ln_get_solar_geo_coords(JD)

	# Calc distance and light time
	dist = ln_get_rect_distance (body_rect_posn, sol_rect_posn)
	t = ln_get_light_time (dist)
	
	# repeat calculation with new time (i.e. JD - t) 
	body_rect_posn = ln_get_par_helio_rect_posn (orbit, JD - t)
	
	# Calc equ coords equ 33.10
	x = sol_rect_posn.X + body_rect_posn.X
	y = sol_rect_posn.Y + body_rect_posn.Y
	z = sol_rect_posn.Z + body_rect_posn.Z

	posn.ra = ln_range_degrees(ln_rad_to_deg(atan2(y, x)))
	posn.dec = ln_rad_to_deg(atan2(z,sqrt(x * x + y * y)))
	return posn


# def ln_get_par_body_earth_dist( JD, ln_par_orbit) -> float
# param JD Julian day.
# param orbit Orbital parameters
# returnss Distance in AU

# Calculate the distance between a body and the Earth for the given julian day.

def ln_get_par_body_earth_dist( JD : float, orbit : ln_par_orbit ) -> float:
	earth_rect_posn = ln_rect_posn()

	# get solar and body rect coords
	body_rect_posn = ln_get_par_geo_rect_posn(orbit, JD)
	earth_rect_posn.X = 0
	earth_rect_posn.Y = 0
	earth_rect_posn.Z = 0
	return ln_get_rect_distance(body_rect_posn, earth_rect_posn)


# def ln_get_par_body_solar_dist( JD, orbit) -> float
# param JD Julian Day.
# param orbit Orbital parameters
# returns The distance in AU between the Sun and the body. 

# Calculate the distance between a body and the Sun.

def ln_get_par_body_solar_dist( JD : float, orbit : ln_par_orbit ) -> float:
	sol_rect_posn = ln_rect_posn()
	
	# get solar and body rect coords 
	body_rect_posn = ln_get_par_helio_rect_posn(orbit, JD)
	sol_rect_posn.X = 0
	sol_rect_posn.Y = 0
	sol_rect_posn.Z = 0
	
	# Calc distance
	return ln_get_rect_distance(body_rect_posn, sol_rect_posn)


# def ln_get_par_body_phase_angle( JD, orbit ) -> float
# param JD Julian day
# param orbit Orbital parameters
# returns Phase angle.

# Calculate the phase angle of the body. The angle Sun - body - Earth. 

def ln_get_par_body_phase_angle( JD : float, orbit : ln_par_orbit ) -> float:

	# time since perihelion 
	t = JD - orbit.JD
	
	# get radius vector 
	r = ln_get_par_radius_vector(orbit.q, t)
	
	# get solar and Earth-Sun distances 
	R = ln_get_earth_solar_dist(JD)
	d = ln_get_par_body_solar_dist(JD, orbit)

	phase = (r * r + d * d - R * R) / (2.0 * r * d )
	return ln_range_degrees(ln_rad_to_deg(acos(phase)))


# def ln_get_par_body_elong( JD, orbit ) -> float
# param JD Julian day
# param orbit Orbital parameters
# returns Elongation to the Sun.

# Calculate the bodies elongation to the Sun.. 

def ln_get_par_body_elong( JD : float, orbit : ln_par_orbit ) -> float:

	# time since perihelion 
	t = JD - orbit.JD
	
	# get radius vector 
	r = ln_get_par_radius_vector (orbit.q, t)
	
	# get solar and Earth-Sun distances 
	R = ln_get_earth_solar_dist(JD)
	d = ln_get_par_body_solar_dist(JD, orbit)

	elong = (R * R + d * d - r * r) / ( 2.0 * R * d )
	return ln_range_degrees(ln_rad_to_deg(acos(elong)))

# def ln_get_par_body_rst( JD, observer, orbit) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)

# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with a parabolic orbit for the given Julian day.

# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day either above the horizon. Returns -1 when it remains whole day below the horizon.

def ln_get_par_body_rst( JD : float, observer : ln_lnlat_posn, orbit : ln_par_orbit ) -> (int, ln_rst_time):
	return ln_get_par_body_rst_horizon(JD, observer, orbit, LN_STAR_STANDART_HORIZON)

# def ln_get_par_body_rst_horizon( JD, observer, orbit, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, else 1 for circumpolar.
#
# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with a parabolic orbit for the given Julian day.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day either above the horizon. Returns -1 when it remains whole day below the horizon.

def ln_get_par_body_rst_horizon( JD : float, observer : ln_lnlat_posn, orbit : ln_par_orbit, horizon : float ) -> (int, ln_rst_time):
	return ln_get_motion_body_rst_horizon(JD, observer, ln_get_par_body_equ_coords, orbit, horizon)


# def ln_get_par_body_next_rst( JD, observer, orbit) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)

# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an parabolic orbit for the given Julian day.

# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.

# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_par_body_next_rst( JD : float, observer : ln_lnlat_posn, orbit : ln_par_orbit ) -> (int, ln_rst_time):
	return ln_get_par_body_next_rst_horizon(JD, observer, orbit, LN_STAR_STANDART_HORIZON)


# def ln_get_par_body_next_rst_horizon( JD, observer, orbit, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an parabolic orbit for the given Julian day.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_par_body_next_rst_horizon( JD : float, observer : ln_lnlat_posn, orbit : ln_par_orbit, horizon : float ) -> (int, ln_rst_time):
	return ln_get_motion_body_next_rst_horizon(JD, observer, ln_get_par_body_equ_coords, orbit, horizon)


# def ln_get_par_body_next_rst_horizon_future( JD, observer, orbit, horizon, day_limit) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param orbit Orbital parameters
# param horizon Horizon height
# param day_limit Maximal number of days that will be searched for next rise and set
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, else 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)

# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of a body with an parabolic orbit for the given Julian day.

# This function guarantee, that rise, set and transit will be in <JD, JD + day_limit> range.

# Note: this functions returns 1 if the body is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day below the horizon.

def ln_get_par_body_next_rst_horizon_future( JD : float, observer : ln_lnlat_posn, orbit : ln_par_orbit, horizon : float,  day_limit : int ) -> (int, ln_rst_time):
	return ln_get_motion_body_next_rst_horizon_future(JD, observer, ln_get_par_body_equ_coords, orbit, horizon, day_limit)

