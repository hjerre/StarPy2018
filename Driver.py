from libNova.ln_types import *
from libNova.ln_classes import *
from libNova.mars import *
from libNova.julian_day import *
from libNova.lunar import *

observer = ln_lnlat_posn(lo=12.50, la=57.5)



def solar1():
    jd = ln_get_julian_day( ln_date(2018,1,8))
    print("jd = ", jd)
    pos = ln_get_solar_equ_coords(jd)
    print("Lng.RA   : ", pos.ra)
    print("Lat.Decl : ", pos.dec)

def solar2():
    jd = ln_get_julian_day(ln_date(2018,1, 8))
    (cir, rst) = ln_get_solar_rst(jd, observer)
    dtr = ln_get_date(rst.rise)
    dtt = ln_get_date(rst.transit)
    dts = ln_get_date(rst.set)
    print("Rise    : ", dtr.hours, ":", dtr.minutes, "    Set : ", dts.days, " ", dts.hours, ":", dts.minutes)



def moon1():
    jd = ln_get_julian_day(ln_date(2018,1, 8))
    pos = ln_get_lunar_equ_coords(jd)
    print("JD = ", jd)
    print("lunar RA : %f    decl : %f " % (pos.ra, pos.dec))
    print("Lunar-Earth distance : %f " % ln_get_lunar_earth_dist(jd))
    print("lunar disk          :  %f " % ln_get_lunar_disk(jd))
    print("lunar phase       :  %f " % ln_get_lunar_phase(jd))
    print("lunar bright limb :  %f " % ln_get_lunar_bright_limb(jd))

    (cp, rst) = ln_get_lunar_rst(jd, 0)
    if cp:
        print("Circumpolar")
    else:
        rise    = ln_get_local_date(rst.rise)
        transit = ln_get_local_date(rst.transit)
        set     = ln_get_local_date(rst.set)
        print_date("Rise", rise)
        print_date("Transit", transit)
        print_date("Set", set)



def moon2():
    for day in range(1, 31):
        jd = ln_get_julian_day(ln_date(2018,1, 8))
        cir, ph = ln_get_lunar_rst(jd, observer)
        dt = ln_get_date(ph.rise)
        print("Day : ", day, "   Rise : ", dt.days, " ", dt.hours, ":",dt.minutes)


def database():
    orbit1 = find_comet("Halley")
    print("ORBIT(Halley) : ", orbit1)
    orbit2 = find_asteroid("Ceres")
    print("ORBIT(Ceres)  : ", orbit2)



def print_date( title : str, date : ln_zonedate):
    print(title)
    print(" Year    : ", date.years)
    print(" Month   : ", date.months)
    print(" Day     : ",  date.days)
    print(" Hours   : " ,  date.hours)
    print(" Minutes : " ,  date.minutes)
    print(" Seconds : %6.2f" % date.seconds)


def planet():
    plan = lookup("Jupiter")
    jd = 2458126.3725
    print("Earth Distance : ", plan.earth_distance(jd))
    print("Sun   Distance : ", plan.solar_distance(jd))


def solarsystem():
    jd = 2458126.3725
    for name in solsystem:
        planet = lookup(name)
        if planet is not None:
            try:
                print("Name : %s \t  solar distance : %f A.U  \t earth distance : %f A.U " % (name, planet.solar_distance(jd), planet.earth_distance(jd)))
            except UnimplementedMethod as e:
                print("Name : %s  \t  is unimplemented" % name)


def comet():
    orbit = find_comet("Halley")
    if orbit is not None:
        print("Orbit(Halley) : ", orbit)
    else:
        print("No orbital elements found")

def asteroid():
    orbit = find_comet("Eros")
    if orbit is not None:
        print("Orbit(Eros) : ", orbit)
        print("inc         : ", orbit.i)
    else:
        print("No orbital elements found")





#mars()
#moon1()
#solar2()
solarsystem()
#venus()
#planet()
#database()
#jsonreader()
asteroid()
#comet()
