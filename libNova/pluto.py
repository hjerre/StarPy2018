from libNova.ln_types import *
from libNova.utility import *
from libNova.earth import *
from libNova.transform import *
from math import cos, sin, asin, atan2, acos, sqrt, log10

PLUTO_COEFFS = 43

argument = [
	[0, 0, 1],
	[0, 0, 2],
	[0, 0, 3],
	[0, 0, 4],
	[0, 0, 5],
	[0, 0, 6],
	[0, 1, -1],
	[0, 1, 0],
	[0, 1, 1],
	[0, 1, 2],
	[0, 1, 3],
	[0, 2, -2],
	[0, 2, -1],
	[0, 2, 0],
	[1, -1, 0],
	[1, -1, 1],
	[1, 0, -3],
	[1, 0, -2],
	[1, 0, -1],
	[1, 0, 0],
	[1, 0, 1],
	[1, 0, 2],
	[1, 0, 3],
	[1, 0, 4],
	[1, 1, -3],
	[1, 1, -2],
	[1, 1, -1],
	[1, 1, 0],
	[1, 1, 1],
	[1, 1, 3],
	[2, 0, -6],
	[2, 0, -5],
	[2, 0, -4],
	[2, 0, -3],
	[2, 0, -2],
	[2, 0, -1],
	[2, 0, 0],
	[2, 0, 1],
	[2, 0, 2],
	[2, 0, 3],
	[3, 0, -2],
	[3, 0, -1],
	[3, 0, 0],
]


longitude = [
	[-19799805, 19850055],
	[897144, -4954829],
	[611149, 1211027],
	[-341243, -189585],
	[129287, -34992],
	[-38164, 30893],
	[20442, -9987],
	[-4063, -5071],
	[-6016, -3336],
	[-3956, 3039],
	[-667, 3572],
	[1276, 501],
	[1152, -917],
	[630, -1277],
	[2571, -459],
	[899, -1449],
	[-1016, 1043],
	[-2343, -1012],
	[7042, 788],
	[1199, -338],
	[418, -67],
	[120, -274],
	[-60, -159],
	[-82, -29],
	[-36, -20],
	[-40, 7],
	[-14, 22],
	[4, 13],
	[5,2],
	[-1,0],
	[2,0],
	[-4, 5],
	[4, -7],
	[14, 24],
	[-49, -34],
	[163, -48],
	[9, 24],
	[-4, 1],
	[-3,1],
	[1,3],
	[-3, -1],
	[5, -3],
	[0,0],
]

latitude = [
	[-5452852, -14974862],
	[3527812, 1672790],
	[-1050748, 327647],
	[178690, -292153],
	[18650, 100340],
	[-30697, -25823],
	[4878, 11248],
	[226, -64],
	[2030, -836],
	[69, -604],
	[-247, -567],
	[-57, 1],
	[-122, 175],
	[-49, -164],
	[-197, 199],
	[-25, 217],
	[589, -248],
	[-269, 711],
	[185, 193],
	[315, 807],
	[-130, -43],
	[5, 3],
	[2, 17],
	[2, 5],
	[2, 3],
	[3, 1],
	[2, -1],
	[1, -1],
	[0, -1],
	[0, 0],
	[0, -2],
	[2, 2],
	[-7, 0],
	[10, -8],
	[-3, 20],
	[6, 5],
	[14, 17],
	[-2, 0],
	[0, 0],
	[0, 0],
	[0, 1],
	[0, 0],
	[1, 0],
] 	

radius = [
	[66865439, 68951812],
	[-11827535, -332538],
	[1593179, -1438890],
	[-18444, 483220],
	[-65977, -85431],
	[31174, -6032],
	[-5794, 22161],
	[4601, 4032],
	[-1729, 234],
	[-415, 702],
	[239, 723],
	[67, -67],
	[1034, -451],
	[-129, 504],
	[480, -231],
	[2, -441],
	[-3359, 265],
	[7856, -7832],
	[36, 45763],
	[8663, 8547],
	[-809, -769],
	[263, -144],
	[-126, 32],
	[-35, -16],
	[-19, -4],
	[-15, 8],
	[-4, 12],
	[5, 6],
	[3, 1],
	[6, -2],
	[2, 2],
	[-2, -2],
	[14, 13],
	[-63, 13],
	[136, -236],
	[273, 1065],
	[251, 149],
	[-25, -9],
	[9, -2],
	[-8, 7],
	[2, -10],
	[19, 35],
	[10, 2],
]


# def ln_get_pluto_equ_coords( JD ) -> ln_equ_posn
# param  JD julian Day
# param  position Pointer to store position

# Calculates Pluto's equatorial position for the given julian day.

def ln_get_pluto_equ_coords( JD : float ) -> ln_equ_posn:
	position = ln_equ_posn()
	hasmore = True

	# need typdef for solar heliocentric coords
	h_sol = ln_get_solar_geom_coords(JD)
	g_sol = ln_get_rect_from_helio(h_sol)
	
	while hasmore:
		last = t
		h_pluto = ln_get_pluto_helio_coords(JD - t)
		g_pluto = ln_get_rect_from_helio(h_pluto)

		# equ 33.10 pg 229 
		a = g_sol.X + g_pluto.X
		b = g_sol.Y + g_pluto.Y
		c = g_sol.Z + g_pluto.Z
	
		delta = a * a + b * b + c * c
		delta = sqrt(delta)
		t = delta * 0.0057755183
		diff = t - last
		hasmore = (diff > 0.0001 or diff < -0.0001)
		
	ra = atan2(b, a)
	dec = c / delta
	dec = asin(dec)

	# back to hours, degrees 
	position.ra = ln_range_degrees(ln_rad_to_deg(ra))
	position.dec = ln_rad_to_deg(dec)
	return position


# def ln_get_pluto_helio_coords( JD) -> ln_helio_posn
# param  JD Julian Day
# param  position Pointer to store new heliocentric position

# Calculate Pluto's heliocentric coordinates for the given julian day.
# This function is accurate to within 0.07" in longitude, 0.02" in latitude
# and 0.000006 AU in radius vector.
#
# Note: This function is not valid outside the period of 1885-2099.


def ln_get_pluto_helio_coords( JD : float ) -> ln_helio_posn:
	position = ln_helio_posn()


	# get julian centuries since J2000
	t =(JD - 2451545.0) / 36525.0

	# calculate mean longitudes for jupiter, saturn and pluto
	J =  34.35 + 3034.9057 * t
	S =  50.08 + 1222.1138 * t
	P = 238.96 +  144.9600 * t
	sum_longitude = 0.0
	sum_latitude = 0.0
	sum_radius = 0.0

	for i in range(PLUTO_COEFFS):
		a = argument[i][0] * J + argument[i][1] * S + argument[i][2] * P
		sin_a = sin(ln_deg_to_rad(a))
		cos_a = cos(ln_deg_to_rad(a))

		# longitude
		sum_longitude += longitude[i][0] * sin_a + longitude[i][1] * cos_a

		# latitude
		sum_latitude += latitude[i][0] * sin_a + latitude[i][1] * cos_a

		# radius
		sum_radius += radius[i][0] * sin_a + radius[i][1] * cos_a


	# calc L, B, R
	position.L = 238.958116 + 144.96 * t + sum_longitude * 0.000001
	position.B = -3.908239 + sum_latitude * 0.000001
	position.R = 40.7241346 + sum_radius * 0.0000001
	return position

# def ln_get_pluto_earth_dist( JD ) -> float
# param  JD Julian day
# returns Distance in AU

# Calculates the distance in AU between the Earth and Pluto for the given julian day.

def ln_get_pluto_earth_dist (JD : float ) -> float:
	# get heliocentric positions
	h_pluto = ln_get_pluto_helio_coords(JD)
	h_earth = ln_get_earth_helio_coords(JD)

	# get geocentric coords
	g_pluto = ln_get_rect_from_helio(h_pluto)
	g_earth = ln_get_rect_from_helio(h_earth)

	# use pythag
	x = g_pluto.X - g_earth.X
	y = g_pluto.Y - g_earth.Y
	z = g_pluto.Z - g_earth.Z
	x = x * x
	y = y * y
	z = z * z

	return sqrt(x + y + z)


# def ln_get_pluto_solar_dist( JD ) -> float
# param  JD Julian day
# returns Distance in AU
#
# Calculates the distance in AU between the Sun and Pluto for the given julian day.
 
def ln_get_pluto_solar_dist( JD : float ) -> float:
	# get heliocentric position
	h_pluto = ln_get_pluto_helio_coords(JD)
	return h_pluto.R

	
# def ln_get_pluto_magnitude( JD ) -> float
# param  JD Julian day
# returns Visible magnitude of Pluto
#
# Calculate the visible magnitude of Pluto for the given julian day.
 
def ln_get_pluto_magnitude( JD : float ) -> float:
	# get distances
	r = ln_get_pluto_solar_dist(JD)
	delta = ln_get_pluto_earth_dist(JD)

	return -1.0 + 5.0 * log10(r * delta)


# def ln_get_pluto_disk( JD ) -> float
# param  JD Julian day
# returns Illuminated fraction of Plutos disk
#
# Calculate the illuminated fraction of Pluto's disk for the given julian day.

def ln_get_pluto_disk( JD : float ) -> float:
	# get distances
	R = ln_get_earth_solar_dist(JD)
	r = ln_get_pluto_solar_dist(JD)
	delta = ln_get_pluto_earth_dist(JD)

	# calc fraction angle
	return (((r + delta) * (r + delta)) - R * R) / (4.0 * r * delta)


# def ln_get_pluto_phase( JD ) -> float
# param  JD Julian day
# returns Phase angle of Pluto (degrees)
#
# Calculate the phase angle of Pluto (Sun - Pluto - Earth) for the given julian day.
 
def ln_get_pluto_phase( JD : float ) -> float:
	# get distances
	R = ln_get_earth_solar_dist(JD)
	r = ln_get_pluto_solar_dist(JD)
	delta = ln_get_pluto_earth_dist(JD)

	# calc phase
	i = (r * r + delta * delta - R * R) / (2.0 * r * delta)
	i = acos(i)
	return ln_rad_to_deg(i)



# def ln_get_pluto_rst( JD, observer) -> (int, ln_rst_time)
# param  JD Julian day
# param  observer Observers position
# param  rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, else 1 for circumpolar.
#
# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination)
# time of Pluto for the given Julian day.
#
# Note: this functions returns 1 if Pluto is circumpolar, that is it remains the whole
# day either above or below the horizon.

def ln_get_pluto_rst( JD : float, observer : ln_lnlat_posn ) -> (int, ln_rst_time):
	return ln_get_body_rst_horizon(JD, observer, ln_get_pluto_equ_coords, LN_STAR_STANDART_HORIZON)



# def ln_get_pluto_sdiam( JD ) -> float
# param  JD Julian day
# returns Semidiameter in arc seconds

# Calculate the semidiameter of Pluto in arc seconds for the given julian day.

def ln_get_pluto_sdiam( JD : float ) -> float:
	So = 2.07        # at 1 AU

	dist = ln_get_pluto_earth_dist(JD)
	return So / dist


# def ln_get_pluto_rect_helio( JD )-> ln_rect_posn
# param  JD Julian day.
# param  position pointer to return position
#
# Calculate Plutos rectangular heliocentric coordinates for the given Julian day. Coordinates are in AU.

def ln_get_pluto_rect_helio( JD : float ) -> ln_rect_posn:
	pluto = ln_get_pluto_helio_coords(JD)
	return ln_get_rect_from_helio(pluto)

