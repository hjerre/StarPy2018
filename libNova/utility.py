# $Id: utility.c,v 1.18 2009-04-20 07:17:00 pkubanek Exp $
#
# Copyright (C) 1999, 2000 Juan Carlos Remis
# Copyright (C) 2002 Liam Girdwood <lgirdwood@gmail.com>
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#

import sqlite3
from math import cos, sin, tan, fabs, sqrt, pow
from libNova.ln_types import *
import json

# Conversion factors between degrees and radians
D2R = (1.7453292519943295769e-2)  # deg.radian
R2D = (5.7295779513082320877e1)  # radian.deg
R2S = (2.0626480624709635516e5)  # arc seconds per radian
S2R = (4.8481368110953599359e-6)  # radians per arc second

# Golden ratio
GOLDEN = 1.61803398875

DM_PI = (2 * M_PI)
RADIAN = (180.0 / M_PI                                                                                                                                                                                                                                  )


# convert radians to degrees

def ln_rad_to_deg( radians: float ) -> float:
    return (radians * R2D)


# convert degrees to radians
def ln_deg_to_rad( degrees: float ) -> float:
    return (degrees * D2R)


# convert hours:mins:secs to degrees
def ln_hms_to_deg( hms: ln_hms ) -> float:
    degrees = float(hms.hours) / 24.0 * 360.0
    degrees += float(hms.minutes) / 60.0 * 15.0
    degrees += float(hms.seconds) / 60.0 * 0.25
    return degrees


# convert hours:mins:secs to radians
def ln_hms_to_rad( hms: ln_hms ) -> float:
    radians = float(hms.hours) / 24.0 * 2.0 * M_PI
    radians += float(hms.minutes) / 60.0 * 2.0 * M_PI / 24.0
    radians += float(hms.seconds) / 60.0 * 2.0 * M_PI / 1440.0
    return radians


# convert degrees to hh:mm:ss
def ln_deg_to_hms( degrees: float ) -> ln_hms:
    hms = ln_hms()

    degrees = ln_range_degrees(degrees)

    # divide degrees by 15 to get the hours
    dtemp = degrees / 15.0
    hms.hours = int(dtemp)

    # multiply remainder by 60 to get minutes
    dtemp = 60.0 * (dtemp - hms.hours)
    hms.minutes = int(dtemp)

    # multiply remainder by 60 to get seconds
    hms.seconds = 60.0 * (dtemp - hms.minutes)

    # catch any overflows
    if hms.seconds > 59:
        hms.seconds = 0.0
        hms.minutes += 1

    if hms.minutes > 59:
        hms.minutes = 0
        hms.hours += 1
    return hms


# convert radians to hh:mm:ss
def ln_rad_to_hms( radians: float ) -> ln_hms:
    dd = ln_hms()

    radians = ln_range_radians(radians)
    degrees = ln_rad_to_deg(radians)

    dd = ln_deg_to_hms(degrees)
    return dd


# convert dms to degrees
def ln_dms_to_deg( dms: ln_dms ) -> float:
    degrees = fabs(float(dms.degrees))
    degrees += fabs(float(dms.minutes) / 60.0)
    degrees += fabs(float(dms.seconds / 3600.0))

    # negative ?
    if dms.neg == 1:
        degrees *= -1.0

    return degrees


# convert dms to radians
def ln_dms_to_rad( dms: ln_dms ) -> float:
    radians = fabs(float(dms.degrees) / 360.0 * 2.0 * M_PI)
    radians += fabs(float(dms.minutes) / 21600.0 * 2.0 * M_PI)
    radians += fabs(float(dms.seconds) / 1296000.0 * 2.0 * M_PI)

    # negative ?
    if dms.neg:
        radians *= -1.0

    return radians

# puts a large angle in the correct range -2PI - 2PI radians preserve sign
def ln_range_radians2( angle : float ) -> float:
    if angle > (-2.0 * M_PI) and angle < (2.0 * M_PI):
        return angle

        temp = (int)(angle / (M_PI * 2.0))
        temp *= (M_PI * 2.0)
        return angle - temp

# puts a large angle in the correct range 0 - 360 degrees
def ln_range_degrees( angle : float) -> float:
    if angle >= 0.0 and angle < 360.0:
        return angle

        temp = (int)(angle / 360)
        if angle < 0.0:
                temp -= 1
    temp *= 360
    return angle - temp



# convert radians to degrees
def ln_rad_to_deg( radians : float ) -> float:
        return (radians * R2D)

# convert degrees to radians
def ln_deg_to_rad( degrees : float ) -> float:
        return (degrees * D2R)

# convert degrees to dms
def ln_deg_to_dms( degrees: float ) -> ln_dms:
    dms = ln_dms()


    if degrees >= 0.0:
        dms.neg = 0
    else:
        dms.neg = 1

    degrees = fabs(degrees)
    dms.degrees = int(degrees)

    # multiply remainder by 60 to get minutes
    dtemp = 60.0 * (degrees - dms.degrees)
    dms.minutes = int(dtemp)

    # multiply remainder by 60 to get seconds
    dms.seconds = 60.0 * (dtemp - dms.minutes)

    # catch any overflows
    if dms.seconds > 59:
        dms.seconds = 0.0
        dms.minutes += 1

    if dms.minutes > 59:
        dms.minutes = 0
        dms.degrees += 1
    return dms


# convert radians to dms
def ln_rad_to_dms( radians: float ) -> ln_dms:
    degrees = ln_rad_to_deg(radians)
    return ln_deg_to_dms(degrees)


# puts a large angle in the correct range 0 - 360 degrees
def ln_range_degrees( angle: float ) -> float:
    if angle >= 0.0 and angle < 360.0:
        return angle

    temp = int(angle / 360)
    if angle < 0.0:    temp -= 1
    temp *= 360.0
    return angle - temp


# puts a large angle in the correct range 0 - 2PI radians
def ln_range_radians( angle: float ) -> float:
    if angle >= 0.0 and angle < (2.0 * M_PI):
        return angle

    temp = int((angle / (M_PI * 2.0)))

    if angle < 0.0:
        temp -= 1
    temp *= (M_PI * 2.0)
    return angle - temp


# puts a large angle in the correct range -2PI - 2PI radians preserve sign
def ln_range_radians2( angle: float ) -> float:
    if angle > (-2.0 * M_PI) and angle < (2.0 * M_PI):
        return angle
    temp = int((angle / (M_PI * 2.0)))
    temp *= (M_PI * 2.0)
    return angle - temp


# add seconds to hms
def ln_add_secs_hms( hms: ln_hms, seconds: float ):
    source_hms = ln_hms()

    # breaks double seconds int hms
    source_hms.hours = int(seconds / 3600)
    seconds -= source_hms.hours * 3600
    source_hms.minutes = int(seconds / 60)
    seconds -= source_hms.minutes * 60
    source_hms.seconds = seconds

    # add hms to hms
    ln_add_hms(source_hms, hms)


# add hms to hms
def ln_add_hms( source: ln_hms, dest: ln_hms ):
    dest.seconds += source.seconds
    if dest.seconds >= 60.0:
        # carry
        source.minutes += 1
        dest.seconds -= 60.0
    else:
        if dest.seconds < 0.0:
            # carry
            source.minutes -= 1
            dest.seconds += 60.0
    dest.minutes += source.minutes
    if dest.minutes >= 60:
        # carry
        source.hours += 1
        dest.minutes -= 60
    else:
        if dest.seconds < 0.0:
            # carry
            source.hours -= 1
            dest.minutes += 60

    dest.hours += source.hours


# def ln_hequ_to_equ( hpos ) -> ln_equ_posn
#  brief human readable equatorial position to double equatorial position

def ln_hequ_to_equ( hpos: lnh_equ_posn ) -> ln_equ_posn:
    pos = ln_equ_posn()
    pos.ra = ln_hms_to_deg(hpos.ra)
    pos.dec = ln_dms_to_deg(hpos.dec)
    return pos


# def ln_equ_to_hequ( pos ) -> lnh_equ_posn
#  brief human double equatorial position to human readable equatorial position
#  ingroup conversion

def ln_equ_to_hequ( pos: ln_equ_posn ) -> ln_equ_posn:
    hpos = lnh_equ_posn()
    hpos.ra = ln_deg_to_hms(pos.ra)
    hpos.dec = ln_deg_to_dms(pos.dec)
    return hpos


# def ln_hhrz_to_hrz( hpos ) -> ln_hrz_posn
#  brief human readable horizontal position to double horizontal position

def ln_hhrz_to_hrz( hpos: ln_hrz_posn ) -> ln_hrz_posn:
    pos = ln_hrz_posn()
    pos.alt = ln_dms_to_deg(hpos.alt)
    pos.az = ln_dms_to_deg(hpos.az)
    return pos


# def ln_hrz_to_hhrz( pos ) -> lnh_hrz_posn
#  brief double horizontal position to human readable horizontal position

def ln_hrz_to_hhrz( pos: ln_hrz_posn ) -> lnh_hrz_posn:
    hpos = lnh_hrz_posn()
    hpos.alt = ln_deg_to_dms(pos.alt)
    hpos.az = ln_deg_to_dms(pos.az)
    return hpos


# def ln_hrz_to_nswe( pos ) -> str
#  brief returns direction of given azimuth - like N,S,W,E,NSW,...

def ln_hrz_to_nswe( pos: ln_hrz_posn ) -> str:
    directions = ["S", "SSW", "SW", "SWW", "W", "NWW", "NW", "NNW", "N", "NNE", "NE", "NEE", "E", "SEE", "SE", "SSE"]
    return directions[int(pos.az / 22.5)]


# def ln_hlnlat_to_lnlat( hpos ) -> ln_lnlat_posn
#  brief human readable long/lat position to double long/lat position

def ln_hlnlat_to_lnlat( hpos: lnh_lnlat_posn ) -> ln_lnlat_posn:
    pos = ln_lnlat_posn()
    pos.lng = ln_dms_to_deg(hpos.lng)
    pos.lat = ln_dms_to_deg(hpos.lat)
    return pos


# def ln_lnlat_to_hlnlat( pos ) -> lnh_lnlat_posn
#  brief double long/lat position to human readable long/lat position

def ln_lnlat_to_hlnlat( pos: ln_lnlat_posn ) -> lnh_lnlat_posn:
    hpos = lnh_lnlat_posn()
    hpos.lng = ln_deg_to_dms(pos.lng)
    hpos.lat = ln_deg_to_dms(pos.lat)
    return hpos


# def ln_get_rect_distance( a, b ) -> float
#  param a First rectangular coordinate
#  param b Second rectangular coordinate
#  return Distance between a and b.

# Calculate the distance between rectangular points a and b.

def ln_get_rect_distance( a: ln_rect_posn, b: ln_rect_posn ) -> float:
    x = a.X - b.X
    y = a.Y - b.Y
    z = a.Z - b.Z

    x *= x
    y *= y
    z *= z

    return sqrt(x + y + z)


# def ln_get_light_time( dist ) -> double
#  param dist Distance in AU
#  return Distance in light days.

# Convert units of AU into light days.

def ln_get_light_time( dist: float ) -> float:
    return dist * 0.005775183


"""
# local types and macros 
typedef int BOOL;
#define TRUE 1
#define FALSE 0
#define iswhite(c)  ((c)== ' ' || (c)=='\t')

/*
[]------------------------------------------------------------------------[]
|  trim() & strip()                                                        |
|                                                                          |
|  strips trailing whitespaces from buf.                                   |
|                                                                          |
[]------------------------------------------------------------------------[]
*/
static char *trim(char *x)
{
	char *y;

	if(!x)
		return(x);
	y = x + strlen(x)-1;
	while (y >= x && isspace(*y)) 
		*y-- = 0; /* skip white space */
	return x;
}


/*
[]------------------------------------------------------------------------[]
|                                                                          |
|   skipwhite()                                                            |
|   salta espacios en blanco                                               |
|                                                                          |
[]------------------------------------------------------------------------[]
*/
static void skipwhite(char **s)
{
   while (iswhite(**s))
		(*s)++;
}


/*! \fn double ln_get_dec_location(char * s)
* \param s Location string
* \return angle in degrees
*
* Obtains Latitude, Longitude, RA or Declination from a string.
*
*  If the last char is N/S doesn't accept more than 90 degrees.            
*  If it is E/W doesn't accept more than 180 degrees.                      
*  If they are hours don't accept more than 24:00                          
*                                                                          
*  Any position can be expressed as follows:                               
*  (please use a 8 bits charset if you want                                
*  to view the degrees separator char '0xba')                              
*
*  42.30.35,53                                                             
*  90º0'0,01 W                                                             
*  42º30'35.53 N                                                           
*  42º30'35.53S                                                            
*  42º30'N                                                                 
*  - 42.30.35.53                                                           
*   42:30:35.53 S                                                          
*  + 42.30.35.53                                                           
*  +42º30 35,53                                                            
*   23h36'45,0                                                             
*                                                                          
*                                                                          
*  42:30:35.53 S = -42º30'35.53"                                           
*  + 42 30.35.53 S the same previous position, the plus (+) sign is        
*  considered like an error, the last 'S' has precedence over the sign     
*                                                                          
*  90º0'0,01 N ERROR: +- 90º0'00.00" latitude limit                        
*
*/
double ln_get_dec_location(char *s)
{
	char *ptr, *dec, *hh, *ame, *tok_ptr;
	BOOL negative = FALSE;
	char delim1[] = " :.,;DdHhMm'\n\t";
	char delim2[] = " NSEWnsew\"\n\t";
	int dghh = 0, minutes = 0;
	double seconds = 0.0, pos;
	short count;
	enum {
	HOURS, DEGREES, LAT, LONG
	} type;

	if (s == NULL || !*s)
	return -0.0;

	count = strlen(s) + 1;
	if ((ptr = (char *) alloca(count)) == NULL)
	return -0.0;

	memcpy(ptr, s, count);
	trim(ptr);
	skipwhite(&ptr);
	if (*ptr == '+' || *ptr == '-')
	negative = (char) (*ptr++ == '-' ? TRUE : negative);

	/* the last letter has precedence over the sign */
	if (strpbrk(ptr,"SsWw") != NULL)
	negative = TRUE;

	skipwhite(&ptr);
	if ((hh = strpbrk(ptr,"Hh")) != NULL && hh < ptr + 3) {
	type = HOURS;
	if (negative) /* if RA no negative numbers */
		negative = FALSE;
	} else if ((ame = strpbrk(ptr,"SsNn")) != NULL) {
		type = LAT;
		if (ame == ptr) /* the North/South found before data */
			ptr++;
		} else
			type = DEGREES; /* unspecified, the caller must control it */
	if ((ptr = strtok_r(ptr,delim1, &tok_ptr)) != NULL)
		dghh = atoi (ptr);
	else
		return (-0.0);
	if ((ptr = strtok_r(NULL,delim1, &tok_ptr)) != NULL) {
		minutes = atoi (ptr);
		if (minutes > 59)
			return -0.0;
	} else
		return -0.0;

	if ((ptr = strtok_r(NULL,delim2,&tok_ptr)) != NULL) {
		if ((dec = strchr(ptr,',')) != NULL)
			*dec = '.';
		seconds = strtod (ptr, NULL);
		if (seconds >= 60.0)
			return -0.0;
	}

	if ((ptr = strtok(NULL," \n\t")) != NULL) {
		skipwhite(&ptr);
		if (*ptr == 'S' || *ptr == 'W' || *ptr == 's' || *ptr == 'w')
			negative = TRUE;
	}

	pos = dghh + minutes /60.0 + seconds / 3600.0;
	if (type == HOURS && pos > 24.0)
		return -0.0;
	if (type == LAT && pos > 90.0)
		return -0.0;
	if (negative == TRUE)
		pos = 0.0 - pos;

	return pos;
}
"""

# def  ln_get_humanr_location(double location)    -> str
# param location Location angle in degress
# return Angle string
#
# Obtains a human readable location in the form: ddºmm'ss.ss"
# String must be freed after use.

def ln_get_humanr_location( location : float) -> str:
    pass
    #sec = fabs(60.0 * (modf(location, &deg)))
    #sec = 60.0 * (modf(sec, &min))
    #buf = str(deg) + "º" + str(int(min)) + "'" + str(sec)"

    #return buf




# def ln_interpolate3 ( n : double, y1 : float, y2 : float, y3 : float) -> float
# return interpolation value
#  param n Interpolation factor
#  param y1 Argument 1
#  param y2 Argument 2
#  param y3 Argument 3

# Calculate an intermediate value of the 3 arguments for the given interpolation factor.

def ln_interpolate3( n: float, y1: float, y2: float, y3: float ) -> float:
    # equ 3.2
    a = y2 - y1
    b = y3 - y2
    c = b - a

    # equ 3.3
    y = y2 + n / 2.0 * (a + b + n * c)

    return y


# def ln_interpolate5 (double n, double y1, double y2, double y3, double y4, double y5) -> float
#  return interpolation value
#  param n Interpolation factor
#  param y1 Argument 1
#  param y2 Argument 2
#  param y3 Argument 3
#  param y4 Argument 4
#  param y5 Argument 5

# Calculate an intermediate value of the 5 arguments for the given interpolation factor.

def ln_interpolate5( n: float, y1: float, y2: float, y3: float, y4: float, y5: float ) -> float:
    # equ 3.8
    A = y2 - y1
    B = y3 - y2
    C = y4 - y3
    D = y5 - y4
    E = B - A
    F = C - B
    G = D - C
    H = F - E
    J = G - F
    K = J - H

    y = 0.0
    n2 = n * n
    n3 = n2 * n
    n4 = n3 * n

    y += y3
    y += n * ((B + C) / 2.0 - (H + J) / 12.0)
    y += n2 * (F / 2.0 - K / 24.0)
    y += n3 * ((H + J) / 12.0)
    y += n4 * (K / 24.0)

    return y


def frange(start, stop, step = 1.0):
    i = start
    while i < stop:
        yield i
        i += step



# double ln_find_zero(func, from:float, to:float, double *arg) -> double
#  param f Function to find zero (root place)
#  param from Lower bound of search interval
#  param to Upper bound of search interval
#  param arg Pointer to the other parameters of the function f

# Find zero of function f() at given interval by Newton method.

def  ln_find_zero( func, ffrom:float, to:float, arg) -> float:
    i = 0
    cont = True

    x1 = to
    x = ffrom

    while cont:
        f = func(x1, arg)
        x2 = x1 - f * (x1 - x) / (f - func(x, arg))
        x = x1
        x1 = x2
        i += 1
        cont = ((fabs(x - x1) > 1e-6) and (i < 1000))

    return x2

# double ln_find_max(double (*f) (double, double *), double from, double to, double *arg)
# param f Function to find maximum
# param from Lower bound of search interval
# param to Upper bound of search interval
# param arg Pointer to the other parameters of the function f
#
# Find local maximum of function f() at given interval by Golden Section method.

"""
double ln_find_max( func, from : float, to : float, double *arg)
{
	double	a, b, xl, xu, eps;

	a = from;
	b = to;
	xl = b - (b - a) / GOLDEN;
	xu = a + (b - a) / GOLDEN;
	eps = fabs(b - a);

	do {
		if (func(xl, arg) > func(xu, arg)) {
			b = xu;
			xu = xl;
			xl = b - (b - a) / GOLDEN;
		} else {
			a = xl;
			xl = xu;
			xu = a + (b - a) / GOLDEN;
		}
		eps = fabs(b - a);

	} while (eps > 1e-6);

	return (xu + xl) * 0.5;
}

/* Catches calls to the POSIX strtok_r and converts them to a related WIN32 version. */
char *strtok_r(char *str, const char *sep, char **last)
{
#ifndef __MINGW__
	return strtok_s(str, sep, last);
#else
	return strtok(str, sep);
#endif // !__MINGW__
}

#endif /* __WIN32__ */

/* C89 substitutions for C99 functions. */
#ifdef __C89_SUB__

"""

# Simple cube root
def cbrt(x : float) -> float:
    return pow(x, 1.0 / 3.0)


def find_comet( name: str ) -> ln_ell_orbit:
    with open("/usr/local/lib/libnova/NEWMPCORB.DAT", "r") as file:
       for line in file:
           if name in line:
                orbit = ln_ell_orbit()
                arr = line.split("#")
                #print("a : ", arr[10])
                #print("e : ", arr[8])
                #print("n : ", arr[9])
                #print("i : ", arr[7])
                #print("perdist : ", arr[5])
                #print("Node : ", arr[6])
                orbit.a = float(arr[10])
                orbit.e = float(arr[8])
                orbit.i = float(arr[7])
                orbit.n = float(arr[9])
                #orbit.w = float(line["w"])
                orbit.omega = float(arr[6])
                return orbit

    return None


def find_asteroid( name: str ) -> ln_ell_orbit:
    with open("/usr/local/lib/libnova/NEWNEA.DAT", "r") as file:
        for line in file:
            if name in line:
                orbit = ln_ell_orbit()
                arr = line.split("#")
                #print("a : ", arr[10])
                #print("e : ", arr[8])
                #print("n : ", arr[9])
                #print("i : ", arr[7])
                #print("perdist : ", arr[5])
                #print("Node : ", arr[6])
                orbit.a = float(arr[10])
                orbit.e = float(arr[8])
                orbit.i = float(arr[7])
                orbit.n = float(arr[9])
                # orbit.w = float(line["w"])
                orbit.omega = float(arr[6])
                return orbit

    return None


def find_sao( name: str ):
    with open("/usr/local/lib/libnova/NEWSAO.DAT", "r") as file:
       for line in file:
            if name in line:
                arr = line.split("$")
                data = ln_sao()
                data.ra = float(arr[1])
                data.dec = float(arr[2])
                data.ep = float(arr[3])
                data.rapm = float(arr[4])
                data.decpm = float(arr[5])
                data.magv = float(arr[6])
                data.fk5no = float(arr[7])
                data.spec = float(arr[9])
                data.HDno = float(arr[10])
                return data

    return None

def print_date( title : str, date : ln_zonedate):
    print(title)
    print(" Year    : ",  date.years)
    print(" Month   : " ,  date.months)
    print(" Day     : ", date.days)
    print(" Hours   : " ,  date.hours)
    print(" Minutes : ", date.minutes)
    print(" Seconds : " ,  date.seconds)