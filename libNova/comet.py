
from libNova.ln_types import *
from libNova.elliptic_motion import *
from libNova.solar import *
from libNova.parabolic_motion import *
from math import log10

# def ln_get_ell_comet_mag( JD, orbit, g, k ) -> float
# param JD Julian day.
# param orbit Orbital parameters
# param g Absolute magnitude
# param k Comet constant
# return The visual magnitude.

# Calculate the visual magnitude of a comet in an elliptical orbit.

def ln_get_ell_comet_mag( JD : float, orbit : ln_ell_orbit, g : float, k : float) -> float:
	
	# get mean anomaly
	if orbit.n == 0:
		orbit.n = ln_get_ell_mean_motion (orbit.a)
	M = ln_get_ell_mean_anomaly(orbit.n, JD - orbit.JD)
	
	# get eccentric anomaly
	E = ln_solve_kepler(orbit.e, M)
	
	# get radius vector
	r = ln_get_ell_radius_vector(orbit.a, orbit.e, E)
	d = ln_get_ell_body_solar_dist(JD, orbit)
	
	return g + 5.0 * log10(d) + k * log10(r)


# def ln_get_par_comet_mag( JD, orbit, g, k ) -> float
# param JD 			Julian day.
# param orbit 	 	Orbital parameters
# param absmag	 	Absolute magnitude
# param comconst	 Comet constant
# return The visual magnitude.

# Calculate the visual magnitude of a comet in a parabolic orbit.

def ln_get_par_comet_mag(JD : float, orbit : ln_par_orbit, absmag : float, comconst : float) -> float:
	# time since perihelion
	t = JD - orbit.JD

	# get radius vector
	r = ln_get_par_radius_vector(orbit.absmag, t)
	d = ln_get_par_body_solar_dist(JD, orbit)

	return absmag + 5.0 * log10(d) + comconst * log10(r)

