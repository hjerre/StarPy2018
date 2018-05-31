from libNova.utility import *
from libNova.nutation import *

# def ln_get_mean_sidereal_time( JD ) -> float
# param JD Julian Day
# return Mean sidereal time.

# Calculate the mean sidereal time at the meridian of Greenwich of a given date.

def ln_get_mean_sidereal_time( JD : float ) -> float:
	T =(JD - 2451545.0) / 36525.0

	# calc mean angle
	sidereal = 280.46061837 + (360.98564736629 *(JD - 2451545.0)) + (0.000387933 * T * T) - (T * T * T / 38710000.0)

	# add a convenient multiple of 360 degrees
	sidereal = ln_range_degrees(sidereal)

	# change to hours *
	sidereal *= 24.0 / 360.0

	return sidereal


# def ln_get_apparent_sidereal_time( JD ) -> float
# param JD Julian Day
# return Apparent sidereal time (hours).

# Calculate the apparent sidereal time at the meridian of Greenwich of a given date.

def ln_get_apparent_sidereal_time( JD : float ) -> float:
    # get the mean sidereal time
    sidereal = ln_get_mean_sidereal_time(JD)

    # add corrections for nutation in longitude and for the true obliquity of the ecliptic
    nutation = ln_get_nutation(JD)
    correction = (nutation.longitude / 15.0 * cos(ln_deg_to_rad(nutation.obliquity)))
    sidereal += correction
    return sidereal
