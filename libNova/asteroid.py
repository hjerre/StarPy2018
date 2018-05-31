from libNova.ln_types import *
from libNova.elliptic_motion import *
from libNova.solar import *
from libNova.earth import *

from math import tan, exp, log10

# def ln_get_asteroid_mag( JD, orbit, H, G ) -> float
# param JD Julian day.
# param orbit Orbital parameters
# param H Mean absolute visual magnitude
# param G Slope parameter
# return The visual magnitude.

# Calculate the visual magnitude of an asteroid.

def ln_get_asteroid_mag( JD : float, orbit : ln_ell_orbit, H : float, G : float ) -> float:
	
	# get phase angle
	b = ln_get_ell_body_phase_angle(JD, orbit)
	b = ln_deg_to_rad(b)
	
	# get mean anomaly
	if orbit.n == 0:
		orbit.n = ln_get_ell_mean_motion(orbit.a)
	M = ln_get_ell_mean_anomaly(orbit.n, JD - orbit.JD)
	
	# get eccentric anomaly
	E = ln_solve_kepler(orbit.e, M)
	
	# get radius vector
	r = ln_get_ell_radius_vector(orbit.a, orbit.e, E)
	d = ln_get_ell_body_solar_dist(JD, orbit)
	
	t1 = exp(-3.33 * pow(tan(b / 2.0), 0.63))
	t2 = exp(-0.187 * pow(tan(b / 2.0), 1.22))
	
	return H + 5.0 * log10(r * d) - 2.5 * log10((1.0 - G) * t1 + G * t2)


# def ln_get_asteroid_sdiam_km( H, A ) -> float
# param H Absolute magnitude of asteroid
# param A Albedo of asteroid
# return Semidiameter in km

# Calculate the semidiameter of an asteroid in km.

# Note: Many asteroids have an irregular shape and therefore this function returns an approximate value of the diameter.

def ln_get_asteroid_sdiam_km (H : float, A : float) -> float:
	return 3.13 - 0.2 * H - (0.5 * log10(A))


# def ln_get_asteroid_sdiam_arc( JD, orbit, H, A ) -> float
# param JD Julian day
# param orbit Orbital parameters
# param H Absolute magnitude of asteroid
# param A Albedo of asteroid
# return Semidiameter in seconds of arc

# Calculate the semidiameter of an asteroid in arc seconds.

# Note: Many asteroids have an irregular shape and therefore this function returns an approximate value of the diameter.

def ln_get_asteroid_sdiam_arc(JD : float, orbit : ln_ell_orbit, H : float, A : float) -> float:
	# calc distance to Earth in AU
	dist = ln_get_ell_body_earth_dist(JD, orbit)

	d = 3.13 - 0.2 * H - (0.5 * log10(A))
	return 0.0013788 * d / dist

