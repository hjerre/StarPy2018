from libNova.julian_day import *
from libNova.lunar import *
from libNova.ln_classes import *
from libNova.asteroid import *
import pandas as pd
import matplotlib.pyplot as plt

#JD = 2458126.37250    #ln_get_julian_day(ln_date(2018, 1, 11))
JD = ln_get_julian_day('1/2/2018')

env = Environment("København")
print("Latitude : ", env.observer.lat)
print("Longitude: ", env.observer.lng)


def solar():
    print("----------------  SOLAR --------------\n")

    print("JD : ", JD)

    # geometric coordinates
    pos = ln_get_solar_geom_coords(JD)
    print("Solar Coords longitude (deg)    : %10.6f" % pos.L)
    print("             latitude (deg)     : ", pos.B)
    print("             radius vector (AU) : ", pos.R)

    # ra, dec
    equ = ln_get_solar_equ_coords(JD)
    print("Solar Position RA  : ", equ.ra)
    print("               DEC :" , equ.dec)

    # rise, set and transit
    (cp, rst) = ln_get_solar_rst(JD, env.observer)
    if cp == 1:
        print("Sun is circumpolar")
    else:
        print("Rise    : ", rst.rise)
        print("Transit : ", rst.transit)
        print("Set     : ", rst.set)
        print("Length  : %2.2f" % ln_get_day_length(JD, env.observer))
        rise = ln_get_local_date(rst.rise)
        transit = ln_get_local_date(rst.transit)
        set = ln_get_local_date(rst.set)
        print_date("Rise", rise)
        print_date("Transit", transit)
        print_date("Set", set)


def lunar():
    print("--------------- Lunar ------------\n")

    print("JD : %6.2f " % JD)

    # get the lunar geopcentric position in km, earth is at 0, 0, 0
    moon = ln_get_lunar_geo_posn(JD, 0)
    print("lunar x : %f    y : %f    z : %f " % (moon.X, moon.Y, moon.Z))

    # Long Lat
    ecl = ln_get_lunar_ecl_coords(JD, 0)
    print("lunar long : %f   lat : %f " % ( ecl.lng, ecl.lat) )

    # RA, DEC
    equ = ln_get_lunar_equ_coords(JD)
    print("lunar RA : %f     Dec : %f " % (equ.ra, equ.dec))

    # moon earth distance
    print("lunar distance km : %f" % ln_get_lunar_earth_dist(JD))

    # lunar disk, phase and bright limb
    print( "lunar disk        :  %f "  % ln_get_lunar_disk(JD) )
    print( "lunar phase       :  %f "  % ln_get_lunar_phase(JD) )
    print( "lunar bright limb :  %f "  % ln_get_lunar_bright_limb(JD) )

    # rise, set and transit time
    (cp, rst) = ln_get_lunar_rst(JD, env.observer)
    if cp == 1:
        print( "Moon is circumpolar")
    else:
        #print("Rise    : ", rst.rise)
        #print("Transit : ", rst.transit)
        #print("Set     : ", rst.set)
        rise = ln_get_local_date(rst.rise)
        transit = ln_get_local_date(rst.transit)
        set = ln_get_local_date(rst.set)
        print_date("Rise",  rise)
        print_date("Transit",  transit)
        print_date("Set",  set)



def mars():
    print("-------------------- Mars -----------------\n")

    print( "JD : " , JD )

    # longitude, latitude and radius vector
    pos = ln_get_mars_helio_coords(JD)
    print( "Mars   L : %f    B : %f     R : %f " % (pos.L, pos.B, pos.R))

    # RA, DEC
    equ = ln_get_mars_equ_coords(JD)
    hequ = ln_equ_to_hequ( equ)
    print("Mars RA %d:%d:%f Dec %d:%d:%f" % (hequ.ra.hours, hequ.ra.minutes, hequ.ra.seconds, hequ.dec.degrees,
       hequ.dec.minutes, hequ.dec.seconds))

    # Earth - Mars dist AU
    au = ln_get_mars_earth_dist(JD)
    print( "mars -> Earth dist (AU) : %f "  % au )

    # Sun - Mars Dist AU
    au = ln_get_mars_solar_dist(JD)
    print( "mars -> Sun dist (AU) : %f"  % au )

    # Mars disk, magnitude and phase
    au = ln_get_mars_disk(JD)
    print( "mars -> illuminated disk : %f"  % au )
    au = ln_get_mars_magnitude(JD)
    print( "mars -> magnitude : %f"  % au )
    au = ln_get_mars_phase(JD)
    print( "mars -> phase : %f" % au )

    # rise, set and transit time
    (cp, rst) = ln_get_mars_rst(JD, env.observer)
    if cp == 1:
        print( "Moon is circumpolar")
    else:
        print("Rise    : ", rst.rise)
        print("Transit : ", rst.transit)
        print("Set     : ", rst.set)
        rise = ln_get_local_date(rst.rise)
        transit = ln_get_local_date(rst.transit)
        set = ln_get_local_date(rst.set)
        print_date("Rise",  rise)
        print_date("Transit",  transit)
        print_date("Set",  set)


def comet():
    print("---------------- COMET ------------\n")
    epoch_date = ln_date()
    orbit = ln_ell_orbit()

    print( "JD : " , JD )

    # calc epoch JD
    epoch_date.years = 1990
    epoch_date.months = 10
    epoch_date.days = 28
    epoch_date.hours = 12
    epoch_date.minutes = 30
    epoch_date.seconds = 0
    e_JD = ln_get_julian_day(epoch_date)

    # Enckle orbital elements
    orbit.JD = e_JD
    orbit.a = 2.2091404
    orbit.e = 0.8502196
    orbit.i = 11.94525
    orbit.omega = 334.75006
    orbit.w = 186.23352
    orbit.n = 0

    # solve kepler for orbit
    E = ln_solve_kepler(0.1, 5.0)
    print( "(Equation of kepler) E when e is 0.1 and M is 5.0  : %f " % E )

    # true anomaly
    v = ln_get_ell_true_anomaly(0.1, E)
    print( "(True Anomaly) v when e is 0.1 and E is 5.5545  : %f " % v )

    # mean anomaly
    v = ln_get_ell_mean_anomaly(orbit.n, JD - orbit.JD)
    print("(Mean Anomaly) v when e is 0.1 and E is 5.5545  : %f " % v)

    # radius vector
    r = ln_get_ell_radius_vector(0.5, 0.1, E)
    print( "(Radius Vector) r when v is , e is 0.1 and E is 5.5545  : %f"  % r )

    # geocentric rect coords
    posn = ln_get_ell_geo_rect_posn( orbit, JD)
    print( "(Geocentric Rect Coords X) for comet Enckle  : %f " % posn.X )
    print( "(Geocentric Rect Coords Y) for comet Enckle  : %f " % posn.Y )
    print( "(Geocentric Rect Coords Z) for comet Enckle  : %f " % posn.Z )

    # rectangular coords
    posn = ln_get_ell_helio_rect_posn( orbit, JD)
    print( "(Heliocentric Rect Coords X) for comet Enckle  : %f " % posn.X )
    print( "(Heliocentric Rect Coords Y) for comet Enckle   : %f " %  posn.Y )
    print( "(Heliocentric Rect Coords Z) for comet Enckle   : %f " %  posn.Z )

    # ra, dec
    equ = ln_get_ell_body_equ_coords(JD,  orbit)
    print( "(RA) for comet Enckle   : %f " %  equ.ra )
    print( "(Dec) for comet Enckle  : %f " % equ.dec )

    # orbit length
    l = ln_get_ell_orbit_len( orbit)
    print( "(Orbit Length) for comet Enckle in AU  : %f " %  l )

    # orbital velocity at perihelion
    V = ln_get_ell_orbit_pvel( orbit)
    print( "(Orbit Perihelion Vel) for comet Enckle in kms  : %f " %  V )

    # orbital velocity at aphelion
    V = ln_get_ell_orbit_avel( orbit)
    print( "(Orbit Aphelion Vel) for comet Enckle in kms  : %f " %  V )

    # average orbital velocity
    V = ln_get_ell_orbit_vel(JD, orbit)
    print( "(Orbit Vel JD) for comet Enckle in kms  : %f " %  V )

    # comet sun distance
    dist = ln_get_ell_body_solar_dist(JD, orbit)
    print( "(Body Solar Dist) for comet Enckle in AU  : %f " %  dist )

    # comet earth distance
    dist = ln_get_ell_body_earth_dist(JD, orbit)
    print( "(Body Earth Dist) for comet Enckle in AU  : %f " %  dist )

    # rise, set and transit
    (cp, rst)  = ln_get_ell_body_rst(JD, env.observer, orbit)
    if cp == 1:
        print("Comet is circumpolar")
    else:
        print("Rise    : ", rst.rise)
        print("Transit : ", rst.transit)
        print("Set     : ", rst.set)
        rise = ln_get_local_date(rst.rise)
        transit = ln_get_local_date(rst.transit)
        set = ln_get_local_date(rst.set)
        print_date("Rise",  rise)
        print_date("Transit",  transit)
        print_date("Set",  set)

def asteroid():
    M_epoch = "K036A"
    print("---------------- Asteroid ------------\n")
    epoch_date = ln_date()
    orbit = ln_ell_orbit()

    print("JD : ", JD)

    # Pallas orbital parameters Taken from MPCORB.DAT

    orbit.a = 2.7730346
    orbit.e = 0.2299839
    orbit.i = 34.84989
    orbit.omega = 173.16479
    orbit.w = 310.45917
    orbit.n = 0.21343771
    H = 4.13
    G = 0.11

    # calc last passage in Perihelion, in julian day

    M_JD = ln_get_julian_from_mpc(M_epoch)

    orbit.JD = ln_get_ell_last_perihelion(M_JD, 260.69458, orbit.n)
    print("JD (Perihelion) %f : " % orbit.JD)

    # calc the earth centered position

    posn = ln_get_ell_geo_rect_posn(orbit, JD)

    print("(Geocentric Rect Coords X) for Pallas   %f : " % posn.X)
    print("(Geocentric Rect Coords Y) for Pallas   %f : " % posn.Y)
    print("(Geocentric Rect Coords Z) for Pallas   %f : " % posn.Z)

    # calc the sun centered position

    posn = ln_get_ell_helio_rect_posn( orbit, JD)
    print("(Heliocentric Rect Coords X) for Pallas   %f : " % posn.X)
    print("(Heliocentric Rect Coords Y) for Pallas    %f : " % posn.Y)
    print("(Heliocentric Rect Coords Z) for Pallas    %f : " % posn.Z)

    # get the RA and Dec

    equ_posn = ln_get_ell_body_equ_coords(JD, orbit)
    print("(RA) for Pallas    %f : " % equ_posn.ra)
    print("(Dec) for Pallas   %f : " % equ_posn.dec)

    # get Alt, Az

    hrz = ln_get_hrz_from_equ( equ_posn, env.observer, JD)
    print("Az  %f : " % hrz.az)
    print("Alt %f : " % hrz.alt)

    # orbit length
    l = ln_get_ell_orbit_len( orbit )
    print("(Orbit Length) for Pallas in AU   %f : " % l)

    # orbit velocities *
    V = ln_get_ell_orbit_pvel(orbit)
    print("(Orbit Perihelion Vel) for Pallas in kms   %f : " % V)
    V = ln_get_ell_orbit_avel( orbit)
    print("(Orbit Aphelion Vel) for Pallas in kms   %f : " % V)
    V = ln_get_ell_orbit_vel(JD, orbit)
    print("(Orbit Vel JD) for Pallas in kms   %f : " % V)

    # earth and solar distance
    dist = ln_get_ell_body_solar_dist(JD, orbit)
    print("Solar Dist (AU)  : %f : " % dist)
    dist = ln_get_ell_body_earth_dist(JD, orbit)
    print("Earth Dist (AU)  : %f : " % dist)

    # phase angle, elongation *
    ph = ln_get_ell_body_phase_angle(JD, orbit)
    print("Phase angle      : %f : " % ph)
    elong = ln_get_ell_body_elong(JD, orbit)
    print("Elongation       : %f : " % elong)

    # magnitude
    mag = ln_get_asteroid_mag(JD, orbit, H, G)
    print("Magnitude        : %f : " % mag)

    # rise, set and transit time
    (cp, rst) = ln_get_ell_body_rst(JD, env.observer, orbit)
    if cp != 0:
        print("Pallas is circumpolar")
    else:
        print("Rise    : ", rst.rise)
        print("Transit : ", rst.transit)
        print("Set     : ", rst.set)
        rise = ln_get_local_date(rst.rise)
        transit = ln_get_local_date(rst.transit)
        set = ln_get_local_date(rst.set)
        print_date("Rise", rise)
        print_date("Transit", transit)
        print_date("Set", set)

def series():
    dates = ln_julian_range('1/1/2018', 365)

    data1 = [ ln_get_mars_earth_dist(jd) for jd in dates]
    data2 = [ ln_get_venus_earth_dist(jd) for jd in dates]
    data3 = [ ln_get_mercury_earth_dist(jd) for jd in dates]

    df = pd.DataFrame(
        { "Mercury" : data3, "Venus" : data2, "Mars" : data1 }, index=dates
    )
    print("Mars smallest distance : ", df["Mercury"].min())
    print("Mars average  distance : ", df["Mercury"].mean())
    print("Mars Biggest  distance : ", df["Mercury"].max())



    # plotter = ln_plotter("distances", "JD", "AU")
    # plotter.addseries("mercury-earth", data3, dates)
    # plotter.addseries("venus-earth", data2, dates)
    # plotter.addseries("mars-earth", data1, dates)
    # plotter.show()


def laenge():
    leng = ln_get_day_length( JD, env.observer )
    print("Day length : ", leng)

def main():
    #solar()
    #lunar()
    #mars()
    #comet()
    #asteroid()
    series()
    #laenge()

main()

