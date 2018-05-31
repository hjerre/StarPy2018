from libNova.ln_types import *
from libNova.utility import *
from math import cos,sin, atan2,acos

TERMS = 36

# data structures to hold arguments and coefficients of Ron-Vondrak theory
class arg:
	def __init__(self):
		self.a_L2 = 0.0
		self.a_L3 = 0.0
		self.a_L4 = 0.0
		self.a_L5 = 0.0
		self.a_L6 = 0.0
		self.a_L7 = 0.0
		self.a_L8 = 0.0
		self.a_LL = 0.0
		self.a_D = 0.0
		self.a_MM = 0.0
		self.a_F = 0.0


class XYZ:
	def __init__(self, sin1 = 0.0, sin2 = 0.0, cos1 = 0.0, cos2 = 0.0):
		self.sin1 = sin1
		self.sin2 = sin2
		self.cos1 = cos1
		self.cos2 = cos2


arguments = [
	[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
	[0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
	[0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
	[0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
	[0, 2, 0, -1, 0, 0, 0, 0, 0, 0, 0],
	[0, 3, -8, 3, 0, 0, 0, 0, 0, 0, 0],
	[0, 5, -8, 3, 0, 0, 0, 0, 0, 0, 0],
	[2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
	[0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
	[0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
	[2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0],
	[0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 3, 0, -2, 0, 0, 0, 0, 0, 0, 0],
	[1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0],
	[2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 3, -2, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 1, 2, -1, 0],
	[8, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[8, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
	[3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0],
	[3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0]
]
	   
x_coefficients = [
  [-1719914, -2, -25, 0],
  [6434, 141, 28007, -107],
  [715, 0, 0, 0],
  [715, 0, 0, 0],
  [486, -5, -236, -4],
  [159, 0, 0, 0],
  [0, 0, 0, 0],
  [39, 0, 0, 0],
  [33, 0, -10, 0],
  [31, 0, 1, 0],
  [8, 0, -28, 0],
  [8, 0, -28, 0],
  [21, 0, 0, 0],
  [-19, 0, 0, 0],
  [17, 0, 0, 0],
  [16, 0, 0, 0],
  [16, 0, 0, 0],
  [11, 0, -1, 0],
  [0, 0, -11, 0],
  [-11, 0, -2, 0],
  [-7, 0, -8, 0],
  [-10, 0, 0, 0],
  [-9, 0, 0, 0],
  [-9, 0, 0, 0],
  [0, 0, -9, 0],
  [0, 0, -9, 0],
  [8, 0, 0, 0],
  [8, 0, 0, 0],
  [-4, 0, -7, 0],
  [-4, 0, -7, 0],
  [-6, 0, -5, 0],
  [-1, 0, -1, 0],
  [4, 0, -6, 0],
  [0, 0, -7, 0],
  [5, 0, -5, 0],
  [5, 0, 0, 0]
]

y_coefficients = [
    [25, -13, 1578089, 156],
    [25697, -95, -5904, -130],
    [6, 0, -657, 0],
	[0, 0, -656, 0],
    [-216, -4, -446, 5],
    [2, 0, -147, 0],
	[0, 0, 26, 0],
	[0, 0, -36, 0],
    [-9, 0, -30, 0],
    [1, 0, -28, 0],
    [25, 0, 8, 0],
    [-25, 0, -8, 0],
	[0, 0, -19, 0],
	[0, 0, 17, 0],
	[0, 0, -16, 0],
	[0, 0, 15, 0],
    [1, 0, -15, 0],
    [-1, 0, -10, 0],
    [-10, 0, 0, 0],
    [-2, 0, 9, 0],
    [-8, 0, 6, 0],
	[0, 0, 9, 0], 
	[0, 0, -9, 0], 
	[0, 0, -8, 0],
    [-8, 0, 0, 0],
    [8, 0, 0, 0],
	[0, 0, -8, 0],
	[0, 0, -7, 0],
    [-6, 0, -4, 0],
    [6, 0, -4, 0],
    [-4, 0, 5, 0],
    [-2, 0, -7, 0],
    [-5, 0, -4, 0],
    [-6, 0, 0, 0],
    [-4, 0, -5, 0],
    [0, 0, -5, 0]
]

z_coefficients = [
    [10, 32, 684185, -358],
    [11141, -48, -2559, -55],
    [-15, 0, -282, 0],
	[0, 0, -285, 0],
    [-94, 0, -193, 0],
    [-6, 0, -61, 0],
	[0, 0, 59, 0],
	[0, 0, 16, 0],
    [-5, 0, -13, 0],
	[0, 0, -12, 0],
    [11, 0, 3, 0],
    [-11, 0, -3, 0],
	[0, 0, -8, 0],
	[0, 0, 8, 0],
	[0, 0, -7, 0],
    [1, 0, 7, 0],
    [-3, 0, -6, 0],
    [-1, 0, 5, 0],
    [-4, 0, 0, 0],
    [-1, 0, 4, 0],
    [-3, 0, 3, 0],
	[0, 0, 4, 0],
	[0, 0, -4, 0],
	[0, 0, -4, 0],
    [-3, 0, 0, 0],
    [3, 0, 0, 0],
	[0, 0, -3, 0],
	[0, 0, -3, 0],
    [-3, 0, 2, 0],
    [3, 0, -2, 0],
    [-2, 0, 2, 0],
    [1, 0, -4, 0],
    [-2, 0, -2, 0],
    [-3, 0, 0, 0],
    [-2, 0, -2, 0],
	[0, 0, -2, 0]
]

# def ln_get_equ_aber( mean_position :  ln_equ_posn, JD : float, position : ln_equ_posn)
#  param: mean_position Mean position of object
#  param: JD Julian Day
#  param: position Pointer to store new object position.

# Calculate a stars equatorial coordinates from it's mean equatorial coordinates
# with the effects of aberration for a given Julian Day.

# Equ 22.3, 22.4

def ln_get_equ_aber( mean_position : ln_equ_posn, JD : float) -> ln_equ_posn:
	position = ln_equ_posn()
	X = 0.0
	Y = 0.0
	Z = 0.0

	# speed of light in 10-8 au per day
	c = 17314463350.0

	# calc T
	T = (JD - 2451545.0) / 36525.0

	# calc planetary perturbutions
	L2 = 3.1761467 + 1021.3285546 * T
	L3 = 1.7534703 + 628.3075849 * T
	L4 = 6.2034809 + 334.0612431 * T
	L5 = 0.5995464 + 52.9690965 * T
	L6 = 0.8740168 + 21.329909095 * T
	L7 = 5.4812939 + 7.4781599 * T
	L8 = 5.3118863 + 3.8133036 * T
	LL = 3.8103444 + 8399.6847337 * T
	D = 5.1984667 + 7771.3771486 * T
	MM = 2.3555559 + 8328.6914289 * T
	F = 1.6279052 + 8433.4661601 * T

	X = 0
	Y = 0
	Z = 0

	# sum the terms
	for i in range(0, TERMS + 1):
		A = arguments[i][0] * L2 + arguments[i].a_L3 * L3 + arguments[i].a_L4 * L4 + arguments[i].a_L5 * L5 + arguments[i].a_L6 * L6 + arguments[i].a_L7 * L7 + \
			arguments[i].a_L8 * L8 + arguments[i].a_LL * LL + arguments[i].a_D * D + arguments[i].a_MM * MM + arguments[i].a_F * F

		X += (x_coefficients[i].sin1 + x_coefficients[i].sin2 * T) * sin(A) + (x_coefficients[i].cos1 + x_coefficients[i].cos2 * T) * cos(A)
		Y += (y_coefficients[i].sin1 + y_coefficients[i].sin2 * T) * sin(A) + (y_coefficients[i].cos1 + y_coefficients[i].cos2 * T) * cos(A)
		Z += (z_coefficients[i].sin1 + z_coefficients[i].sin2 * T) * sin(A) + (z_coefficients[i].cos1 + z_coefficients[i].cos2 * T) * cos(A)


	# Equ 22.4
	mean_ra = ln_deg_to_rad(mean_position.ra)
	mean_dec = ln_deg_to_rad(mean_position.dec)

	if mean_dec < M_PI * 0.4999:
		delta_ra = (Y * cos(mean_ra) - X * sin(mean_ra)) / cos(mean_dec)
		delta_ra /= c
		delta_dec = (X * cos(mean_ra) + Y * sin(mean_ra)) * sin(mean_dec) - Z * cos(mean_dec)
		delta_dec /= -c
	
		position.ra = ln_rad_to_deg(mean_ra + delta_ra)
		position.dec = ln_rad_to_deg(mean_dec + delta_dec)
	else:
		X /= c
		Y /= c
		Z /= c
		cos_dec = cos(mean_dec)

		px = cos_dec * cos(mean_ra)
		py = cos_dec * sin(mean_ra)

		px += X
		py += Y

		ra = atan2(py, px)
		dec = acos(sqrt(px * px + py * py))

		dec += cos_dec * Z

		position.ra = ln_rad_to_deg(ra)
		position.dec = ln_rad_to_deg(dec)
	return position


# def ln_get_ecl_aber( mean_position, JD ) -> ln_lnlat_posn
#  param mean_position: ln_lnlat_posn		 Mean position of object
#  param JD: float							 Julian Day
#  param ecliptical position: ln_lnlat_posn.

# Calculate a stars ecliptical coordinates from it's mean ecliptical coordinates with the effects of aberration for a given Julian Day.

def ln_get_ecl_aber( mean_position : ln_lnlat_posn,  JD : float) -> ln_lnlat_posn:
	position = ln_lnlat_posn()
	sol_position = ln_helio_posn()

	# constant of aberration
	k = ln_deg_to_rad(20.49552 *  (1.0 / 3600.0))

	# Equ 21.1
	T = (JD - 2451545) / 36525
	T2 = T * T

	# suns longitude in radians
	sol_position = ln_get_solar_geom_coords(JD)
	true_longitude = ln_deg_to_rad(sol_position.B)

	# Earth orbit ecentricity
	e = 0.016708617 - 0.000042037 * T - 0.0000001236 * T2
	e = ln_deg_to_rad(e)

	# longitude of perihelion Earths orbit
	t = 102.93735 + 1.71953 * T + 0.000046 * T2
	t = ln_deg_to_rad(t)

	# change object long/lat to radians
	mean_lng = ln_deg_to_rad(mean_position.lng)
	mean_lat = ln_deg_to_rad(mean_position.lat)

	# equ 22.2 */
	delta_lng = (-k * cos(true_longitude - mean_lng) + e * k * cos(t - mean_lng)) / cos(mean_lat)
	delta_lat = -k * sin(mean_lat) * (sin(true_longitude - mean_lng) - e * sin(t - mean_lng))

	mean_lng += delta_lng
	mean_lat += delta_lat

	position.lng = ln_rad_to_deg(mean_lng)
	position.lat = ln_rad_to_deg(mean_lat)
	return position

