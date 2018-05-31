import pandas as pd
from libNova.ln_types import *
import time

# def ln_get_julian_day( date ) -> float
# param: date Date required.
# returns: Julian day

# Calculate the julian day from a calendar day. 
# Valid for positive and negative years but not for negative JD.


def ln_get_julian_day_local( date : ln_date ) -> float:
	local_date = ln_date()
	local_date = date
	if local_date.months < 3:
		local_date.years -= 1
		local_date.months += 12

	
	a = local_date.years / 100
	if local_date.years > 1582 or (local_date.years == 1582 and (local_date.months > 10 or (local_date.months == 10 and local_date.days >= 4))):
		b = int(2 - a + (a / 4))
	else:
		b = 0
	
	days = local_date.days + float(local_date.hours / 24.0) + float(local_date.minutes / 1440.0) + float(local_date.seconds /  86400.0)
	JD = int(365.25 * (local_date.years + 4716)) + int(30.6001 * (local_date.months + 1)) + days + b - 1524.5
	return JD

def ln_get_julian_day( date ):
	if type(date) is ln_date:
		return ln_get_julian_day(date)
	elif type(date) is str:
		ldt = pd.Timestamp(date)
		return ln_get_julian_day_local(
			ln_date( ldt.year, ldt.month, ldt.day, ldt.hour, ldt.minute, ldt.second))
	else:
		return None

# ln_julian_range
# return array of JD from d2 to either d2 OR d1 periods=d2

def ln_julian_range(d1 : str , d2  ) -> list:
	if type(d2) is int:
		dater = pd.date_range(d1, periods=d2)
	else:
		dater = pd.date_range(d1, d2)
	dlist = [ ln_get_julian_day_local(ln_date( dat.year, dat.month, dat.day, dat.hour, dat.minute, dat.second)) for dat in dater ]
	return dlist


# def ln_get_day_of_week( date ) -> int
# param: date:ln_date		 Date required
# returns: Day of the week
#
# Calculate the day of the week.
# Returns 0 = Sunday .. 6 = Saturday

def ln_get_day_of_week( date : ln_date ) -> int:
	# get julian day
	JD = ln_get_julian_day (date)
	JD += 1.5
	day = int(JD) % 7

	return day


# def ln_get_date( JD ) -> date
# param: JD Julian day
# param: date Pointer to new calendar date.
#
# Calculate the date from the Julian day

def ln_get_date( JD : float ) -> ln_date:
	date = ln_date()

	JD += 0.5
	Z = int(JD)
	F = JD - Z
   
	if Z < 2299161:
		A = int(Z)
	else:
		a = int((Z - 1867216.25) / 36524.25)
		A = int(Z + 1 + a - (int)(a / 4))

	B = A + 1524
	C = (int) ((B - 122.1) / 365.25)
	D = (int) (365.25 * C)
	E = (int) ((B - D) / 30.6001)

	# get the hms
	date.hours = int(F * 24)
	F -= date.hours / 24.0
	date.minutes = int(F * 1440)
	F -= date.minutes / 1440.0
	date.seconds = F * 86400

	# get the day
	date.days = B - D - int(30.6001 * E)

	# get the month
	if E < 14:
		date.months = E - 1
	else:
		date.months = E - 13
   
	# get the year
	if date.months > 2:
		date.years = C - 4716
	else:
		date.years = C - 4715
	return date

"""
#  void ln_get_date_from_timet (time_t *t, struct ln_date *date)
# param: t system time
# param: date Pointer to new calendar date.
*
* Set date from system time

void ln_get_date_from_timet (time_t *t, struct ln_date *date)
{
	struct tm gmt;

	# convert to UTC time representation 
	gmtime_r(t, &gmt);
    	
	ln_get_date_from_tm(&gmt, date);
}

#  void ln_get_date_from_tm(struct tm *t, struct ln_date *date)
# param: tm system tm structure
# param: date Pointer to new calendar date.
*
* Set date from system tm structure

void ln_get_date_from_tm(struct tm *t, struct ln_date *date)
{
	# fill in date struct 
	date->seconds = t->tm_sec;
	date->minutes = t->tm_min;
	date->hours = t->tm_hour;
	date->days = t->tm_mday;
	date->months = t->tm_mon + 1;
	date->years = t->tm_year + 1900;
}

#  void ln_get_date_from_sys(struct ln_date *date)
# param: date Pointer to store date.
*
* Calculate local date from system date.

void ln_get_date_from_sys(struct ln_date *date)
{
	struct tm * gmt;
	struct timeval tv;
	struct timezone tz;

	# get current time with microseconds precission
	gettimeofday(&tv, &tz);

	# convert to UTC time representation 
	gmt = gmtime(&tv.tv_sec);
    	
	# fill in date struct 
	date->seconds = gmt->tm_sec + ((float)tv.tv_usec / 1000000);
	date->minutes = gmt->tm_min;
	date->hours = gmt->tm_hour;
	date->days = gmt->tm_mday;
	date->months = gmt->tm_mon + 1;
	date->years = gmt->tm_year + 1900;
}


#  float ln_get_julian_from_timet (time_t * in_time)
# param: time The time_t.
# returns: Julian day.
*
* Calculate Julian day from time_t.

float ln_get_julian_from_timet(time_t * in_time)
{
	# 1.1.1970 = JD 2440587.5 
	return (float)(2440587.5 + (float)(*in_time / (float) 86400.0));
}

#ifndef HAVE_ROUND

# Simple round to nearest 
static inline float round(float x)
{
	return floor(x + 0.5);
}

#endif # ! HAVE_ROUND 

#  void ln_get_timet_from_julian(float JD, time_t * in_time)
# param: JD Julian day
# param: in_time Pointer to store time_t
*
* Calculate time_t from julian day

void ln_get_timet_from_julian(float JD, time_t * in_time)
{	
	*in_time = (time_t)round((JD - (float) 2440587.5) * (float) 86400.0);
}

#  float ln_get_julian_from_sys()
# returns: Julian day (UT)
*
* Calculate the julian day (UT) from the local system time

float ln_get_julian_from_sys()
{
	float JD;
	struct ln_date date;
		
	# get sys date 
	ln_get_date_from_sys(&date);
	JD = ln_get_julian_day(&date);

	return JD;
}
"""
# def ln_get_julian_local_date( zonedate ) -> float
# param: zonedate Local date
# returns: Julian day (UT)
#
# Calculate Julian day (UT) from zone date

def ln_get_julian_local_date( zonedate : ln_zonedate ) -> float:
	date = ln_zonedate_to_date(zonedate)
	return ln_get_julian_day(date)


# def ln_get_local_date( JD) -> ln_zonedate
# param: JD Julian day
# param: zonedate Pointer to new calendar date.
#
# Calculate the zone date from the Julian day (UT). Gets zone info from
# system using either _timezone or tm_gmtoff fields.


def ln_get_local_date( JD : float ) -> ln_zonedate:
	date = ln_get_date(JD)

	curtime = time.time()
	loctime = time.localtime(curtime)
	gmtoff = loctime.tm_gmtoff

	zonedate = ln_date_to_zonedate( date, gmtoff )
	return zonedate


# def ln_get_date_from_mpc( date, mpc_date ) -> int
# param: date Pointer to new calendar date.
# param: mpc_date Pointer to string MPC date
# returns: 0 for valid date
#
# Calculate the local date from the a MPC packed date.
# See http://cfa-www.harvard.edu/iau/info/PackedDates.html for info.

def translate_mpc_segment( l : str ) -> int:
	if l >= "0" and l <= "9":
		return int(l)
	else:
		return 75 - ord(l[0])

def ln_get_date_from_mpc( date : ln_date, mpc_date : str) -> int:
	if len(mpc_date) != 5:
		return -1

	# get the century
	if mpc_date[0] == "I":
			date.years = 1800
	elif mpc_date[0] ==  'J':
			date.years = 1900
	elif mpc_date[0] == 'K':
			date.years = 2000
	else:
			return -1
	
	# get the year 
	year = mpc_date[1] + mpc_date[2]
	date.years += int(year)

	# month 
	month = translate_mpc_segment(mpc_date[3])
	date.months = month
	
	# day
	day = translate_mpc_segment(mpc_date[4])
	date.days = day
	
	# reset hours,min,secs to 0 
	date.hours = 0
	date.minutes = 0
	date.seconds = 0
	return 0


# def ln_get_julian_from_mpc( mpc_date ) -> float
# param: mpc_date Pointer to string MPC date
# returns: Julian day.
#
# Calculate the julian day from the a MPC packed date.
# See http://cfa-www.harvard.edu/iau/info/PackedDates.html for info.

def ln_get_julian_from_mpc( mpc_date : str ) -> float:
	date = ln_date()

	ln_get_date_from_mpc(date, mpc_date)
	JD = ln_get_julian_day(date)
	return JD


# def ln_date_to_zonedate( date, zonedate, gmtoff) -> ln_zonedata
# param: zonedate Ptr to zonedate
# param: gmtoff Offset in seconds from UT
# param: date Ptr to date
#
# Converts a ln_date (UT) to a ln_zonedate (local time).


def ln_date_to_zonedate( date : ln_date, gmtoff : int ) -> ln_zonedate:
	zonedate = ln_zonedate()

	jd = ln_get_julian_day(date)
	jd += gmtoff / 86400.0
	dat = ln_get_date(jd)

	zonedate.years   = dat.years
	zonedate.months  = dat.months
	zonedate.days    = dat.days
	zonedate.hours   = dat.hours
	zonedate.minutes = dat.minutes
	zonedate.seconds = dat.seconds

	zonedate.gmtoff  = gmtoff
	return zonedate


# def ln_zonedate_to_date( zonedate) -> ln_date
# param: zonedate Ptr to zonedate
# param: date Ptr to date
#
# Converts a ln_zonedate (local time) to a ln_date (UT).

def ln_zonedate_to_date( zonedate : ln_zonedate ) -> ln_date:
	dat = ln_date()

	dat.years   = zonedate.years
	dat.months  = zonedate.months
	dat.days    = zonedate.days
	dat.hours   = zonedate.hours
	dat.minutes = zonedate.minutes
	dat.seconds = zonedate.seconds

	jd = ln_get_julian_day(dat)
	jd -= zonedate.gmtoff / 86400.0
	return ln_get_date(jd)
