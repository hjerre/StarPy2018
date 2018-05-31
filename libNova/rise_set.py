from libNova.ln_types import *
from libNova.utility import *
from libNova.sidereal_time import ln_get_apparent_sidereal_time
from libNova.dynamical_time import *

from math import fabs, cos, acos, asin, isnan


# helper function to check if object can be visible
def check_coords( observer : ln_lnlat_posn, H1 : float, horizon : float, object : ln_equ_posn ) ->int:

	# check if body is circumpolar
	if fabs(H1) > 1.0:
		# check if maximal height < horizon
		# h = asin(cos(ln_deg_to_rad(observer.lat - object.dec)))
		h = 90.0 + object.dec - observer.lat
		# normalize to <-90+90>
		if h > 90.0:
			h = 180.0 - h
		if h < -90.0:
			h = -180.0 - h
		if h < horizon:
			return -1
		# else it must be above horizon
		return 1
	return 0


# def ln_get_object_rst( JD, observer, object) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param object Object position
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination) time of the object for the given Julian day.

# Note: this functions returns 1 if the object is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains the whole day bellow the horizon.

def ln_get_object_rst( JD : float, observer : ln_lnlat_posn, object : ln_equ_posn ) -> (int, ln_rst_time):
	return ln_get_object_rst_horizon(JD, observer, object, LN_STAR_STANDART_HORIZON)	# standard altitude of stars


# def ln_get_object_rst_horizon( JD, observer, object, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param object Object position
# param horizon Horizon height
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)

# Calculate the time the rise, set and transit (crosses the local meridian at upper culmination) time of the object for the given Julian day and horizon.

# Note: this functions returns 1 if the object is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains whole day bellow the horizon.



def ln_get_object_rst_horizon_offset( JD : float, observer : ln_lnlat_posn, object : ln_equ_posn, horizon : float, ut_offset : float ) -> ln_rst_time:
	rst = ln_rst_time()

	#if (isnan(ut_offset)) {
	#	JD_UT = JD
	#} else {
		# convert local sidereal time into degrees for 0h of UT on day JD
	jd = int(JD)
	JD_UT = jd + ut_offset


	O = ln_get_apparent_sidereal_time(JD_UT)
	O *= 15.0

	# equ 15.1
	H0 = (sin(ln_deg_to_rad(horizon)) -
		 sin(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(object.dec)))
	H1 = (cos(ln_deg_to_rad(observer.lat)) * cos(ln_deg_to_rad(object.dec)))

	H1 = H0 / H1

	ret = check_coords(observer, H1, horizon, object)
	if ret:
		return 1, ret

	H0 = acos(H1)
	H0 = ln_rad_to_deg(H0)

	# equ 15.2
	mt = (object.ra - observer.lng - O) / 360.0
	mr = mt - H0 / 360.0
	ms = mt + H0 / 360.0

	for i in range(0, 4):
		# put in correct range
		if mt > 1.0:
			mt-= 1
		elif mt < 0:
			mt+=1
		if mr > 1.0:
			mr-=1
		elif mr < 0:
			mr+=1
		if ms > 1.0:
			ms-=1
		elif ms < 0:
			ms+=1

		# find sidereal time at Greenwich, in degrees, for each m
		mst = O + 360.985647 * mt
		msr = O + 360.985647 * mr
		mss = O + 360.985647 * ms

		# find local hour angle
		Hat = mst + observer.lng - object.ra
		Har = msr + observer.lng - object.ra
		Has = mss + observer.lng - object.ra

		# find altitude for rise and set
		altr = sin(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(object.dec)) + cos(ln_deg_to_rad(observer.lat)) * cos(ln_deg_to_rad(object.dec)) * cos(ln_deg_to_rad(Har))
		alts = sin(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(object.dec)) + cos(ln_deg_to_rad(observer.lat)) * cos(ln_deg_to_rad(object.dec)) * cos(ln_deg_to_rad(Has))

		# must be in degrees
		altr = ln_rad_to_deg(altr)
		alts = ln_rad_to_deg(alts)

		#corrections for m
		ln_range_degrees(Hat)
		if Hat > 180.0:
			Hat -= 360

		dmt = -(Hat / 360.0)
		dmr = (altr - horizon) / (360 * cos(ln_deg_to_rad(object.dec)) * cos(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(Har)))
		dms = (alts - horizon) / (360 * cos(ln_deg_to_rad(object.dec)) * cos(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(Has)))

		# add corrections and change to JD
		mt += dmt
		mr += dmr
		ms += dms

		if mt <= 1.0 and mt >= 0.0 and mr <= 1.0 and mr >= 0.0 and ms <= 1.0 and ms >= 0.0:
			break


	rst.rise = JD_UT + mr
	rst.transit = JD_UT + mt
	rst.set = JD_UT + ms

	# not circumpolar
	return 0, rst



def ln_get_object_rst_horizon( JD : float, observer : ln_lnlat_posn, object : ln_equ_posn, horizon : float ) -> (int, ln_rst_time):
	return ln_get_object_rst_horizon_offset(JD, observer, object, horizon, 0.5)


# def ln_get_object_next_rst( JD, observer, object) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param object Object position
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of the object for the given Julian day and horizon.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note: this functions returns 1 if the object is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains whole day bellow the horizon.

def ln_get_object_next_rst( JD : float, observer : ln_lnlat_posn, object : ln_equ_posn ) -> (int, ln_rst_time):
	return ln_get_object_next_rst_horizon(JD, observer, object, LN_STAR_STANDART_HORIZON)


# helper functions for ln_get_object_next_rst_horizon
def set_next_rst( rst : ln_rst_time, diff : float ):
	out = ln_rst_time()
	out.rise = rst.rise + diff
	out.transit = rst.transit + diff
	out.set = rst.set + diff
	return out

def find_next( JD : float, jd1 : float, jd2 : float, jd3 : float) -> float:
	if isnan(jd1) and isnan(jd2):
		return jd3

	if JD < jd1:
		return jd1

	if JD < jd2:
		return jd2

	return jd3


# def ln_get_object_next_rst_horizon(  JD, observer, object, horizon) -> (int, ln_rst_time)
# param JD Julian day
# param observer Observers position
# param object Object position
# param horizon Horizon height
# param rst Pointer to store Rise, Set and Transit time in JD
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination)
# time of the object for the given Julian day and horizon.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note: this functions returns 1 if the object is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains whole day bellow the horizon.

def ln_get_object_next_rst_horizon( JD : float, observer : ln_lnlat_posn, object : ln_equ_posn, horizon : float ) -> (int, ln_rst_time):
	ret, rst = ln_get_object_rst_horizon_offset(JD, observer, object, horizon, nan("0"))
	if ret:
		# circumpolar
		return ret, None

	if rst.rise > (JD + 0.5) or rst.transit > (JD + 0.5) or rst.set > (JD + 0.5):
		x, rst1 = ln_get_object_rst_horizon_offset(JD - 1.0, observer, object, horizon, nan("0"))
	else:
		rst_1 = set_next_rst(rst, -1.0)

	if rst.rise < JD or rst.transit < JD or rst.set < JD:
		x, rst_2 = ln_get_object_rst_horizon_offset(JD + 1.0, observer, object, horizon, nan("0"))
	else:
		rst_2 = set_next_rst (rst, 1.0)

	rst.rise = find_next(JD, rst_1.rise, rst.rise, rst_2.rise)
	rst.transit = find_next(JD, rst_1.transit, rst.transit, rst_2.transit)
	rst.set = find_next(JD, rst_1.set, rst.set, rst_2.set)

	if isnan (rst.rise):
		return ret

	return 0


# def ln_get_body_rst_horizon( JD, observer, func, horizon) -> (int, ln_rst_time)
# param JD Julian day 
# param observer Observers position 
# param get_equ_body_coords Pointer to get_equ_body_coords() function
# param horizon Horizon, see LN_XXX_HORIZON constants
# param rst Pointer to store Rise, Set and Transit time in JD 
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (below the horizon)
#
# Calculate the time the rise, set and transit (crosses the local meridian at
# upper culmination) time of the body for the given Julian day and given horizon.
# Note 1: this functions returns 1 if the object is circumpolar, that is it remains the whole
# day above the horizon. Returns -1 when it remains whole day bellow the horizon.
#
# Note 2: this function will not work for body, which ra changes more
# then 180 deg in one day (get_equ_body_coords changes so much). But
# you should't use that function for any body which moves to fast..use
# some special function for such things.

# Only method called by other python scripts



def ln_get_body_rst_horizon_offset( JD : float, observer : ln_lnlat_posn, get_equ_body_coords, horizon : float, ut_offset : float ) -> (int, ln_rst_time):
	rst = ln_rst_time()

	# dynamical time diff
	T = ln_get_dynamical_time_diff(JD)
	if isnan(ut_offset):
		JD_UT = JD
	else:
		jd = int(JD)
		JD_UT = jd + ut_offset

	# convert local sidereal time into degrees for 0h of UT on day JD
	JD_UT = JD
	O = ln_get_apparent_sidereal_time(JD_UT)
	O *= 15.0
	# get body coords for JD_UT -1, JD_UT and JD_UT + 1
	sol1 = get_equ_body_coords(JD_UT - 1.0)
	sol2 = get_equ_body_coords(JD_UT)
	sol3 = get_equ_body_coords(JD_UT + 1.0)

	# equ 15.1
	H0 = (sin(ln_deg_to_rad(horizon)) - sin(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(sol2.dec)))
	H1 = (cos(ln_deg_to_rad(observer.lat)) * cos(ln_deg_to_rad(sol2.dec)))
	H1 = H0 / H1

	ret = check_coords (observer, H1, horizon, sol2)
	if ret:
		return ret, None

	H0 = acos(H1)
	H0 = ln_rad_to_deg(H0)

	# correct ra values for interpolation	- put them to the same side of circle 
	if (sol1.ra - sol2.ra) > 180.0:
		sol2.ra += 360.0

	if (sol2.ra - sol3.ra) > 180.0:
		sol3.ra += 360.0

	if (sol3.ra - sol2.ra) > 180.0:
		sol3.ra -= 360.0

	if (sol2.ra - sol1.ra) > 180.0:
		sol3.ra -= 360.0

	mt = (sol2.ra - observer.lng - O) / 360.0
	mr = mt - H0 / 360.0
	ms = mt + H0 / 360.0

	for i in range(0, 4):
		# put in correct range
		if mt > 1.0:
			mt -= 1
		elif mt < 0:
			mt += 1
		if mr > 1.0:
			mr-= 1
		elif mr < 0:
			mr+=1
		if ms > 1.0:
			ms-= 1
		elif ms < 0:
			ms+=1

		# find sidereal time at Greenwich, in degrees, for each m 
		mst = O + 360.985647 * mt
		msr = O + 360.985647 * mr
		mss = O + 360.985647 * ms
	
		nt = mt + T / 86400.0
		nr = mr + T / 86400.0
		ns = ms + T / 86400.0

		# interpolate ra and dec for each m, except for transit dec (dec2)
		posr = ln_equ_posn()
		post = ln_equ_posn()
		poss = ln_equ_posn()
		posr.ra = ln_interpolate3(nr, sol1.ra, sol2.ra, sol3.ra)
		posr.dec = ln_interpolate3(nr, sol1.dec, sol2.dec, sol3.dec)
		post.ra = ln_interpolate3(nt, sol1.ra, sol2.ra, sol3.ra)
		poss.ra = ln_interpolate3(ns, sol1.ra, sol2.ra, sol3.ra)
		poss.dec = ln_interpolate3(ns, sol1.dec, sol2.dec, sol3.dec)

		# find local hour angle 
		Hat = mst + observer.lng - post.ra
		Har = msr + observer.lng - posr.ra
		Has = mss + observer.lng - poss.ra

		# find altitude for rise and set 
		altr = sin(ln_deg_to_rad(observer.lat)) * \
			sin(ln_deg_to_rad(posr.dec)) + \
				cos(ln_deg_to_rad(observer.lat)) * \
				cos(ln_deg_to_rad(posr.dec)) * \
				cos(ln_deg_to_rad(Har))
		alts = sin(ln_deg_to_rad(observer.lat)) * \
				sin(ln_deg_to_rad(poss.dec)) + \
				cos(ln_deg_to_rad(observer.lat)) * \
				cos(ln_deg_to_rad(poss.dec)) * \
				cos(ln_deg_to_rad(Has))

		# must be in degrees 
		altr = ln_rad_to_deg(altr)
		alts = ln_rad_to_deg(alts)

		# corrections for m 
		ln_range_degrees(Hat)
		if  Hat > 180.0:
			Hat -= 360

		dmt = -(Hat / 360.0)
		dmr = (altr - horizon) / (360.0 * cos(ln_deg_to_rad(posr.dec)) * \
			cos(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(Har)))
		dms = (alts - horizon) / (360.0 * cos(ln_deg_to_rad(poss.dec)) * \
			cos(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(Has)))

		# add corrections and change to JD 
		mt += dmt
		mr += dmr
		ms += dms

		if mt <= 1.0 and mt >= 0.0 and \
				1.0 >= mr >= 0.0 and \
				1.0 >= ms >= 0.0:
			break


	rst.rise = JD_UT + mr
	rst.transit = JD_UT + mt
	rst.set = JD_UT + ms

	# not circumpolar
	print("rst.rise = ", rst.rise)
	return 0, rst


def ln_get_body_rst_horizon( JD : float, observer :  ln_lnlat_posn, get_equ_body_coords, horizon : float ) -> (int, ln_rst_time):
	return ln_get_body_rst_horizon_offset(JD, observer, get_equ_body_coords, horizon, 0.5)


# def ln_get_body_next_rst_horizon( JD, observer, object, horizon) -> (int, ln_rst_time)
# param JD Julian day 
# param observer Observers position 
# param get_equ_body_coords Pointer to get_equ_body_coords() function
# param horizon Horizon, see LN_XXX_HORIZON constants
# param rst Pointer to store Rise, Set and Transit time in JD 
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at
# upper culmination) time of the body for the given Julian day and given horizon.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note 1: this functions returns 1 if the body is circumpolar, that is it remains
# the whole day either above or below the horizon.
#
# Note 2: This function will not work for body, which ra changes more
# then 180 deg in one day (get_equ_body_coords changes so much). But you should't use that function for any body which moves to fast..use
# some special function for such things.

def ln_get_body_next_rst_horizon( JD : float, observer : ln_lnlat_posn, get_equ_body_coords, horizon : float ) -> (int, ln_rst_time):
	return ln_get_body_next_rst_horizon_future(JD, observer, get_equ_body_coords, horizon, 1)


# def ln_get_body_next_rst_horizon_future( JD, observer, func, horizon, day_limit) -> (int, ln_rst_time)
# param JD Julian day 
# param observer Observers position 
# param get_equ_body_coords Pointer to get_equ_body_coords() function
# param horizon Horizon, see LN_XXX_HORIZON constants
# param day_limit Maximal number of days that will be searched for next rise and set
# param rst Pointer to store Rise, Set and Transit time in JD 
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at
# upper culmination) time of the body for the given Julian day and given horizon.
#
# This function guarantee, that rise, set and transit will be in <JD, JD + day_limit> range.
#
# Note 1: this functions returns 1 if the body is circumpolar, that is it remains
# the whole day either above or below the horizon.
#
# Note 2: This function will not work for body, which ra changes more
# than 180 deg in one day (get_equ_body_coords changes so much). But
# you should't use that function for any body which moves to fast..use
# some special function for such things.

def ln_get_body_next_rst_horizon_future( JD : float, observer :  ln_lnlat_posn, get_equ_body_coords, horizon : float, day_limit : int ) -> (int, ln_rst_time):
	rst_1 = ln_rst_time()
	rst_2 = ln_rst_time()
	rst = None

	ret = ln_get_body_rst_horizon_offset(JD, observer, get_equ_body_coords, horizon, 0.0)

	if ret and day_limit == 1:
		# circumpolar
		return ret, None

	if not ret and (rst.rise > (JD + 0.5) or rst.transit > (JD + 0.5) or rst.set >(JD + 0.5)):
		ret, rst_1 = ln_get_body_rst_horizon_offset(JD - 1, observer, get_equ_body_coords, horizon, 0.0)
		if ret:
			rst = set_next_rst (rst_1, -1)
		else:
			rst.rise = 0.0
			rst.transit = 0.0
			rst.set = 0.0

		rst_1 = set_next_rst(rst, -1)


	if ret or (rst.rise < JD or rst.transit < JD or rst.set < JD):
		# find next day when it will rise, up to day_limit days
		day = 1

		while day <= day_limit:
			ret, rst_2 = ln_get_body_rst_horizon_offset(JD + day, observer, get_equ_body_coords, horizon, nan ("0"))

			if not ret:
				day = day_limit + 2
				break
			day += 1

		if day == day_limit + 1:
			# it's then really circumpolar in searched period
			return ret, None
	#else:
	#	set_next_rst(rst, +1, &rst_2)

	rst.rise = find_next(JD, rst_1.rise, rst.rise, rst_2.rise)
	rst.transit = find_next(JD, rst_1.transit, rst.transit, rst_2.transit)
	rst.set = find_next(JD, rst_1.set, rst.set, rst_2.set)
	if isnan (rst.rise):
		return ret, None

	return 0


# def ln_get_body_rst_horizon( JD, observer, func, horizon) -> (int, ln_rst_time)
# param JD Julian day 
# param observer Observers position 
# param get_motion_body_coords Pointer to ln_get_ell_body_equ_coords. ln_get_para_body_equ_coords or ln_get_hyp_body_equ_coords function
# param horizon Horizon, see LN_XXX_HORIZON constants
# param rst Pointer to store Rise, Set and Transit time in JD 
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time the rise, set and transit (crosses the local meridian at
# upper culmination) time of the body for the given Julian day and given horizon.
#
# Note 1: this functions returns 1 if the body is circumpolar, that is it remains the whole day either above or below the horizon.

def ln_get_motion_body_rst_horizon( JD : float, observer : ln_lnlat_posn, get_motion_body_coords, orbit, horizon : float ) -> (int, ln_rst_time):
	return ln_get_motion_body_rst_horizon_offset(JD, observer, get_motion_body_coords, orbit, horizon, 0.5)


def ln_get_motion_body_rst_horizon_offset( JD : float, observer : ln_lnlat_posn, get_motion_body_coords, orbit, horizon : float, ut_offset : float ) -> (int, ln_rst_time):
	sol1, sol2, sol3 = (None, None, None)
	rst = ln_rst_time()
	posr = ln_equ_posn()
	poss = ln_equ_posn()
	post = ln_equ_posn()
		
	# dynamical time diff 
	T = ln_get_dynamical_time_diff(JD)

	if isnan(ut_offset):
		JD_UT = JD
	else:
		jd = int(JD)
		JD_UT = jd + ut_offset
	O = ln_get_apparent_sidereal_time(JD_UT)
	O *= 15.0
	
	# get body coords for JD_UT -1, JD_UT and JD_UT + 1 
	sol1 = get_motion_body_coords(JD_UT - 1.0, orbit)
	sol2 = get_motion_body_coords(JD_UT, orbit)
	sol3 = get_motion_body_coords(JD_UT + 1.0, orbit)

	# equ 15.1 
	H0 = (sin(ln_deg_to_rad(horizon)) - sin(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(sol2.dec)))
	H1 = (cos(ln_deg_to_rad(observer.lat)) * cos(ln_deg_to_rad(sol2.dec)))

	H1 = H0 / H1

	ret = check_coords(observer, H1, horizon, orbit)
	if ret:
		return ret, None

	H0 = acos(H1)
	H0 = ln_rad_to_deg(H0)

	# correct ra values for interpolation	- put them to the same side of circle 
	if (sol1.ra - sol2.ra) > 180.0:
		sol2.ra += 360.0

	if (sol2.ra - sol3.ra) > 180.0:
		sol3.ra += 360.0

	if (sol3.ra - sol2.ra) > 180.0:
		sol3.ra -= 360.0

	if (sol2.ra - sol1.ra) > 180.0:
		sol3.ra -= 360.0

	for i in range(3):
		# equ 15.2
		mt = (sol2.ra - observer.lng - O) / 360.0
		mr = mt - H0 / 360.0
		ms = mt + H0 / 360.0

		# put in correct range 
		if mt > 1.0:
			mt -=1
		elif mt < 0.0:
			mt += 1
		if mr > 1.0:
			mr -= 1
		elif mr < 0.0:
			mr += 1
		if ms > 1.0:
			ms -= 1
		elif ms < 0.0:
			ms += 1

		# find sidereal time at Greenwich, in degrees, for each m
		mst = O + 360.985647 * mt
		msr = O + 360.985647 * mr
		mss = O + 360.985647 * ms

		nt = mt + T / 86400.0
		nr = mr + T / 86400.0
		ns = ms + T / 86400.0

		# interpolate ra and dec for each m, except for transit dec (dec2) 
		posr.ra = ln_interpolate3(nr, sol1.ra, sol2.ra, sol3.ra)
		posr.dec = ln_interpolate3(nr, sol1.dec, sol2.dec, sol3.dec)
		post.ra = ln_interpolate3(nt, sol1.ra, sol2.ra, sol3.ra)
		poss.ra = ln_interpolate3(ns, sol1.ra, sol2.ra, sol3.ra)
		poss.dec = ln_interpolate3(ns, sol1.dec, sol2.dec, sol3.dec)

		# find local hour angle 
		Hat = mst + observer.lng - post.ra
		Har = msr + observer.lng - posr.ra
		Has = mss + observer.lng - poss.ra

		# find altitude for rise and set 
		altr = sin(ln_deg_to_rad(observer.lat)) * \
				sin(ln_deg_to_rad(posr.dec)) + \
				cos(ln_deg_to_rad(observer.lat)) * \
				cos(ln_deg_to_rad(posr.dec)) * \
				cos(ln_deg_to_rad(Har))
		alts = sin(ln_deg_to_rad(observer.lat)) * \
				sin(ln_deg_to_rad(poss.dec)) + \
				cos(ln_deg_to_rad(observer.lat)) * \
				cos(ln_deg_to_rad(poss.dec)) * \
				cos(ln_deg_to_rad(Has))

		# corrections for m 
		dmt = - (Hat / 360.0)
		dmr = (altr - horizon) / (360.0 * cos(ln_deg_to_rad(posr.dec)) * cos(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(Har)))
		dms = (alts - horizon) / (360.0 * cos(ln_deg_to_rad(poss.dec)) * cos(ln_deg_to_rad(observer.lat)) * sin(ln_deg_to_rad(Has)))

		# add corrections and change to JD 
		mt += dmt
		mr += dmr
		ms += dms

		if mt <= 1.0 and mt >= 0.0 and mr <= 1.0 and mr >= 0.0 and ms <= 1.0 and ms >= 0.0:
			break

	rst.rise = JD_UT + mr
	rst.transit = JD_UT + mt
	rst.set = JD_UT + ms
	
	# not circumpolar 
	return 0, rst


# def ln_get_body_next_rst_horizon( JD, observer, func, horizon) -> (int, ln_rst_time)
# param JD Julian day 
# param observer Observers position 
# param get_motion_body_coords Pointer to ln_get_ell_body_equ_coords. ln_get_para_body_equ_coords or ln_get_hyp_body_equ_coords function
# param horizon Horizon, see LN_XXX_HORIZON constants
# param rst Pointer to store Rise, Set and Transit time in JD 
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)
#
# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination) time of the body for the given Julian day and given horizon.
#
# This function guarantee, that rise, set and transit will be in <JD, JD+1> range.
#
# Note 1: this functions returns 1 if the body is circumpolar, that is it remains the whole day either above or below the horizon.

def ln_get_motion_body_next_rst_horizon( JD : float, observer : ln_lnlat_posn, get_motion_body_coords, orbit,  horizon : float ) -> (int, ln_rst_time):
	return ln_get_motion_body_next_rst_horizon_future(JD, observer, get_motion_body_coords, orbit, horizon, 1)


# def ln_get_motion_body_next_rst_horizon_future( JD, observer, func, horizon, day_limit) -> (int, ln_rst_time)
# param JD Julian day 
# param observer Observers position 
# param get_motion_body_coords Pointer to ln_get_ell_body_equ_coords. ln_get_para_body_equ_coords or ln_get_hyp_body_equ_coords function
# param horizon Horizon, see LN_XXX_HORIZON constants
# param day_limit Maximal number of days that will be searched for next rise and set
# param rst Pointer to store Rise, Set and Transit time in JD 
# returns 0 for success, 1 for circumpolar (above the horizon), -1 for circumpolar (bellow the horizon)

# Calculate the time of next rise, set and transit (crosses the local meridian at upper culmination) time of the body for the given Julian day and given
# horizon.

# This function guarantee, that rise, set and transit will be in <JD, JD + day_limit> range.

# Note 1: this functions returns 1 if the body is circumpolar, that is it remains
# the whole day either above or below the horizon.

def ln_get_motion_body_next_rst_horizon_future( JD : float, observer : ln_lnlat_posn,
											   get_motion_body_coords, orbit, horizon : float, day_limit : int ) -> (int, ln_rst_time):
	rst = ln_rst_time()

	ret, rst = ln_get_motion_body_rst_horizon_offset(JD, observer, get_motion_body_coords, orbit, horizon, 0)
	if ret and day_limit == 1:
		# circumpolar
		return ret, None

	if not ret and (rst.rise >(JD + 0.5) or rst.transit >(JD + 0.5) or rst.set >(JD + 0.5)):
		ret, rst_1 = ln_get_motion_body_rst_horizon_offset(JD - 1.0, observer, get_motion_body_coords, orbit, horizon, 0)
		if ret:
			rst_1 = set_next_rst(rst, -1.0)
	else:
		rst.rise = 0.0
		rst.transit = 0.0
		rst.set = 0.0

		rst_1 = set_next_rst(rst, -1.0)

	if ret or (rst.rise < JD or rst.transit < JD or rst.set < JD):
		# find next day when it will rise, up to day_limit days
		day = 1

		while day <= day_limit:

			ret, rst_2 = ln_get_motion_body_rst_horizon_offset(JD + day, observer, get_motion_body_coords, orbit, horizon, 0)

			if not ret:
				day = day_limit + 2
				break
			day += 1

		if day == day_limit + 1:
			# it's then really circumpolar in searched period
			return ret, None
	else:
		rst_2 = set_next_rst(rst, +1.0)

	rst.rise = find_next(JD, rst_1.rise, rst.rise, rst_2.rise)
	rst.transit = find_next(JD, rst_1.transit, rst.transit, rst_2.transit)
	rst.set = find_next(JD, rst_1.set, rst.set, rst_2.set)

	if isnan (rst.rise):
		return ret, None

	return 0, rst

