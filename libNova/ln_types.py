

# define some useful constants if they are not already defined
M_PI   = 3.1415926535897932384626433832795
M_PI_2 = 1.5707963267948966192313216916398
M_PI_4 = 0.78539816339744830961566084581988

# sidereal day length in seconds and days(for JD)
LN_SIDEREAL_DAY_SEC = 86164.09
LN_SIDEREAL_DAY_DAY = (LN_SIDEREAL_DAY_SEC / 86400.0)

# 1.1.2000 Julian Day & others
JD2000 = 2451545.0
JD2050 = 2469807.50

B1900 = 2415020.3135
B1950 = 2433282.4235

LN_SOLAR_STANDART_HORIZON = -0.8333
LN_SOLAR_CIVIL_HORIZON = -6.0
LN_SOLAR_NAUTIC_HORIZON = -12.0
LN_SOLAR_ASTRONOMICAL_HORIZON = -18.0
LN_STAR_STANDART_HORIZON = -0.5667

solsystem = [ "Mercury", "Venus", "Earth", "Mars" , "Jupiter", "Saturn", "Uranus", "Neptune",  "Pluto" ]




# Date
# class ln_date
#   brief Human readable Date and time used by libnova

# This is the Human readable (easy printf) date format used
# by libnova.It's always in UTC. For local time, please use ln_zonedate.

class ln_date:
    def __init__(self, yy = 0, mm = 0, dd = 0, hh = 0, min = 0, sec = 0.0):
        self.years = yy     # Years.All values are valid
        self.months = mm    # Months.Valid values: 1(January) - 12(December)
        self.days = dd      # Days.Valid values 1 - 28, 29, 30, 31 Depends on month.
        self.hours = hh     # Hours.Valid values 0 - 23.
        self.minutes = min  # Minutes.Valid values 0 - 59.
        self.seconds = sec  # Seconds.Valid values 0 - 59.99999....

    def __str__(self):
        return str(self.days)+"/"+str(self.days)+"/"+str(self.years)+" "+str(self.hours)+":"+str(self.minutes)+":"+str(self.seconds)


# Zonedate
#  class ln_zonedate
#   brief Human readable Date and time with timezone information used by libnova

# This is the Human readable (easy printf) date with timezone format used by libnova.

class ln_zonedate:
    def __init__(self, yy = 0, mm = 0, dd = 0, hh = 0, min = 0, sec = 0.0, goff = 0):
        self.years = yy            # Years.All values are valid
        self.months = mm           # Months.Valid values: 1(January) - 12(December)
        self.days = dd             # Days.Valid values 1 - 28, 29, 30, 31 Depends on month
        self.hours = hh            # Hours.Valid values 0 - 23
        self.minutes = min         # Minutes.Valid values 0 - 59
        self.seconds = sec         # Seconds.Valid values 0 - 59.99999....
        self.gmtoff = goff         # Timezone offset.Seconds east of UTC.Valid values 0..86400


# class ln_dms
#  brief: Degrees, minutes and seconds.

# Human Readable Angle in degrees, minutes and seconds


class ln_dms:
    def __init__(self, n = -1, deg = 0, mn = 0, sec = 0.0):
        self.neg =n         # Non zero if negative
        self.degrees = deg  # Degrees.Valid 0 - 360
        self.minutes = mn   # Minutes.Valid 0 - 59
        self.seconds = sec  # Seconds.Valid 0 - 59.9999...

    def __str__(self):
        return self.neg+" "+str(self.degrees)+":"+str(self.minutes)+":"+str(self.seconds)


# class ln_hms
#  brief: Hours, minutes and seconds.

# Human readable Angle in hours, minutes and seconds


class ln_hms:
    def __init__(self, hh = 0, mn = 0, sec = 0.0):
        self.hours = hh      # Hours.Valid 0 - 23
        self.minutes = mn    # Minutes.Valid 0 - 59
        self.seconds = sec   # Seconds.Valid 0 - 59.9999...

    def __str__(self):
        return str(self.hours)+":"+str(self.minutes)+":"+str(self.seconds)


# Class lnh_equ_posn
#  brief: Right Ascension and Declination.

# Human readable Equatorial Coordinates.


class lnh_equ_posn:
    def __init__(self):
        self.ra = ln_hms()    # Object right ascension.
        self.dec = ln_dms()   # Object declination

    def __str__(self):
        return str(self.ra)+" : " + str(self.dec)


# class lnh_hrz_posn
#  brief: Azimuth and Altitude.

# Human readable Horizontal Coordinates.


class lnh_hrz_posn:
    def __init(self):
        self.az = ln_dms()    # Object azimut
        self.alt = ln_dms()   # Object altitude

    def __str__(self):
        return str(self.az)+" : "+str(self.alt)


# class lnh_lnlat_posn
#  brief: Ecliptical( or celestial) Latitude and Longitude.
#  Human readable Ecliptical( or celestial) Longitude and Latitude.


class lnh_lnlat_posn:
    def __init__(self):
        self.lng = ln_dms()     # Object Longitude
        self.lat = ln_dms()     # Object Latitude

    def __str__(self):
        return str(self.lng)+" : " + str(self.lat)

# class ln_equ_posn
#  brief: Equatorial Coordinates.

# The Right Ascension and Declination of an object.
# Angles are expressed in degrees.

class ln_equ_posn:
    def __init__(self, ra = 0.0, dec = 0.0):
        self.ra = ra     # right ascension in degrees.
        self.dec = dec   # declination in degrees

    def __str__(self):
        return str(self.ra)+" : "+str(self.dec)

# class ln_hrz_posn
# brief: Horizontal Coordinates.
# The Azimuth and Altitude of an object.
# Angles are expressed in degrees.


class ln_hrz_posn:
    def __init__(self, az = 0.0, alt = 0.0):
        self.az = az        #  AZ.Object azimuth. 0 deg = South, 90 deg = West, 180 deg = Nord, 270 deg = East
        self.alt = alt        #  ALT.Object altitude. 0 deg = horizon, 90 deg = zenit, -90 deg = nadir

    def __str__(self):
        return str(self.az)+" : "+str(self.alt)

# class ln_lnlat_posn
#  brief: Ecliptical( or celestial) Longitude and Latitude.
#  The Ecliptical( or celestial) Latitude and Longitude of and object.
# Angles are expressed in degrees.East is positive, West negative.


class ln_lnlat_posn:
    def __init__(self, lo = 0.0, la = 0.0):
        self.lng = lo     #  Object longitude.
        self.lat = la     #  Object latitude

    def __str__(self):
        return str(self.lng)+" : " + str(self.lat)

# class ln_helio_posn
# brief: Heliocentric position
# A heliocentric position is an objects position relative to the centre of the Sun.

# Angles are expressed in degrees.
#  Radius vector is in AU.

class ln_helio_posn:
    def __init__(self, l = 0.0, b = 0.0, r = 0.0):
        self.L = l      #  Heliocentric longitude
        self.B = b      #  Heliocentric latitude
        self.R = r      #  Heliocentric radius vector

    def __str__(self):
        return str(self.L)+" : " + str(self.B)+":"+str(self.R)

# class ln_rect_posn
# brief: Rectangular coordinates
# Rectangular Coordinates of a body. These coordinates can be either geocentric or heliocentric.

# A heliocentric position is an objects position relative to the centre of the Sun.
# A geocentric position is an objects position relative to the centre of the Earth.

# Position is in units of AU for planets and in units of km for the Moon.

class ln_rect_posn:
    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.X = x         # Rectangular X coordinate
        self.Y = y         # Rectangular Y coordinate
        self.Z = z         # Rectangular Z coordinate

    def __str__(self):
        return str(self.X)+" : " + str(self.Y)+":"+str(self.Z)

# Class ln_gal_posn
#  brief: Galactic coordinates
#   The Galactic Latitude and Longitude of and object.

# Angles are expressed in degrees.

class ln_gal_posn:
    def __init__(self, l = 0.0, b = 0.0):
        self.l = l               # Galactic longitude (degrees)
        self.b = b               # Galactic latitude (degrees)

    def __str__(self):
        return str(self.l)+" : " + str(self.b)


# class ln_ell_orbit
#  brief: Elliptic Orbital elements
#  Angles are expressed in degrees.

class ln_ell_orbit:
    def __init__(self, a = 0.0, e = 0.0, i = 0.0, w = 0.0, omega = 0.0, n = 0.0, JD = 0.0):
        self.a = a                #  Semi major axis, in AU
        self.e = e                #  Eccentricity
        self.i = 0                #  Inclination in degrees
        self.w = w                #  Argument of perihelion in degrees
        self.omega = omega        #  Longitude of ascending node in degrees
        self.n = n                #  Mean motion, in degrees / day
        self.JD = JD              #  Time of last passage in Perihelion, in julian day



# class ln_par_orbit
#  brief: Parabolic Orbital elements
#  Angles are expressed in degrees.

class ln_par_orbit:
    def __init__(self, q = 0.0, i = 0.0, w = 0.0, omega = 0.0, JD = 0.0):
        self.q = q          # Perihelion distance in AU
        self.i = i          # Inclination in degrees
        self.w = w          # Argument of perihelion in degrees
        self.omega = omega  # Longitude of ascending node in degrees
        self.JD = JD        # Time of last passage in Perihelion, in julian day


# Class ln_hyp_orbit
#  brief: Hyperbolic Orbital elements
#  Angles are expressed in degrees.


class ln_hyp_orbit:
    def __init__(self, q = 0.0, e = 0.0, i = 0.0, w = 0.0, omega = 0.0, JD = 0.0):
        self.q = q                     #  Perihelion distance in AU
        self.e = e                     #  Eccentricity
        self.i = i                     #  Inclination in degrees
        self.w = w                     #  Argument of perihelion in degrees
        self.omega = omega             #  Longitude of ascending node in degrees
        self.JD = JD                   #  Time of last passage in Perihelion, in julian day

class ln_sao:
    def __init__(self, ra = 0.0, dec = 0.0, ep = 0, rapm = 0.0, decpm = 0.0, magv = 0.0, fk5no = 0, HDno = 0):
        self.ra = ra
        self.dec = dec
        self.ep = ep
        self.rapm = rapm
        self.decpm = decpm
        self.magv = magv
        self.fk5no = fk5no
        self.HDno = HDno



# Class ln_rst_time
# brief: Rise, Set and Transit times.
#  Contains the Rise, Set and transit times for a body.
# Angles are expressed in degrees.

class ln_rst_time:
    def __init__(self, rise = 0.0, set = 0.0, transit = 0.0):
        self.rise = rise         # Rise time in JD
        self.set = set           # Set time in JD
        self.transit = transit   # Transit time in JD



# class ln_nutation
# brief: Nutation in longitude, ecliptic and obliquity.
# Contains Nutation in longitude, obliquity and ecliptic obliquity.
# Angles are expressed in degrees.

class ln_nutation:
    def __init__(self, long = 0.0, ob = 0.0, ecl = 0.0):
        self.longitude = long     #  Nutation in longitude, in degrees
        self.obliquity = ob       #  Nutation in obliquity, in degrees
        self.ecliptic = ecl       #  Mean obliquity of the ecliptic, in degrees


