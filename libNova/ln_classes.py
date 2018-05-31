import abc
import os
import numpy as np
import pandas as pd
import platform
import sqlite3
from matplotlib import pyplot as plt
from libNova.ln_types import *
from libNova.venus import *
from libNova.mercury import *
from libNova.mars import *
from libNova.jupiter import *
from libNova.saturn import *
from libNova.uranus import *
from libNova.neptune import *
from libNova.pluto import *

cities = "/usr/local/lib/libnova/worldcities-basic.csv"



def find_city( name : str) -> ln_lnlat_posn:
    if not os.path.isfile(cities):
        return None
    with open(cities) as file:
        for line in file:
            if name in line:
                att = line.split(",")
                return ln_lnlat_posn( float(att[3]), float(att[2]))
        return None

def find_comet( name : str) -> ln_ell_orbit:
    pass


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Environment(metaclass=Singleton):
    def __init__(self, city = ""):
        self.system = platform.system()
        self.home = "/usr/local/lib/libnova"
        self.platform = platform.platform()
        self.altitude = 0.0
        self.airpressure = 960.0
        self.location = None   #ep_lnlat_pos()
        if city != "":
            self.observer = find_city(city)
            if self.observer is None:
                self.observer = ln_lnlat_posn(0.0, 0.0)



class CelObject(metaclass=abc.ABCMeta):

    def get_distance(self):
        for i in range(1, 9):
            self.distance()
        self.phase()

    @abc.abstractmethod
    def earth_distance(self):
       pass

    @abc.abstractmethod
    def solar_distance(self):
        pass

    @abc.abstractmethod
    def magnitude(self):
        pass

    @abc.abstractmethod
    def phase(self):
        pass

    @abc.abstractmethod
    def illuminatedfraction(self):
        pass

    @abc.abstractmethod
    def riseset(self):
        pass

    @abc.abstractmethod
    def equ_semidiameter(self):
        pass

    @abc.abstractmethod
    def pol_semidiameter(self):
        pass

    @abc.abstractmethod
    def equ_position(self):
        pass

    @abc.abstractmethod
    def helio_position(self):
        pass

class Object(CelObject):
    def __init__(self):
        self.orbitalType="Elliptic"

    def orbit_is_hyperbolic(self):
        self.orbitalType="Hyperbolic"

    def get_earth_distance(self):
        pass

    def get_solar_distance(self):
        pass

    def get_equ_pos(self):
        pass



class Comet:
    def __init__(self, name):
        if type(name) is str:
            self.name = name
            d = CometDatabase()
            self.orbit = d.find(self.name)
            print("Orbit : ", self.orbit)
        elif name is ln_ell_orbit or name is ln_par_orbit:
            self.orbit = name

        def load(self):
            d = CometDatabase()
            self.orbit = d.find(self.name)


        def get_magnitude(self, jd):
            if self.name is ln_ell_orbit:
                return ln_get_ell_comet_mag(jd, self.orbit)
            elif name is ln_par_orbit:
                return ln_get_par_orbit(jd, self.orbit )

    def get_earth_distance(self):
        print("Comet::distance")

    def get_solar_distance(self):
        pass

    def phase(self):
        print("Comet::phase")

    def equ_coord(self):
        print("EquCoord")

class Asteroid(CelObject):
    def get_earth_distance(self):
        pass

    def get_solar_distance(self):
        pass

    def phase(self):
        print("Comet::phase")

    def equ_coord(self):
        print("EquCoord")


class Luna(CelObject):
    def __phase__(self):
        pass

    def rise_time(self):
        pass

    def set_time(self):
        pass


class Sun(CelObject):
    def equ_pos(self, jd):
        return ln_get_solar_equ_coords(jd)

    def ecl_pos(self, jd):
        return ln_get_solar_ecl_coords(jd)

    def geo_pos(self, jd):
        return ln_get_solar_geo_coords(jd)

    def solar_diameter(self, jd):
        return ln_get_solar_sdiam(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_solar_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst

class Mercury(CelObject):
    def earth_distance(self, jd):
        return ln_get_mercury_earth_dist(jd)

    def solar_distance(self, jd):
        return ln_get_mercury_solar_dist(jd)

    def phase(self, jd):
        return ln_get_mercury_phase(jd)

    def magnitude(self, jd):
        return ln_get_mercury_magnitude(jd)

    def equ_semidiameter(self, jd):
        return ln_get_mercury_sdiam(jd)

    def pol_semidiameter(self, jd):
        return ln_get_mercury_sdiam(jd)

    def illuminatedfraction(self, jd):
        return ln_get_mercury_disk(jd)

    def equ_position(self, jd):
        return ln_get_mercury_equ_coords(jd)

    def helio_position(self, jd):
        return ln_get_mercury_helio_coords(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_mercury_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst

class Venus(CelObject):
    def earth_distance(self, jd):
        return ln_get_venus_earth_dist(jd)

    def phase(self, jd):
        return ln_get_venus_phase(jd)

    def magnitude(self, jd):
        return ln_get_mercury_magnitude(jd)

    def solar_distance(self, jd):
        return ln_get_venus_solar_dist(jd)

    def equ_semidiameter(self, jd):
        return ln_get_venus_sdiam(jd)

    def pol_semidiameter(self, jd):
        return ln_get_venus_sdiam(jd)

    def illuminatedfraction(self, jd):
        return ln_get_venus_disk(jd)

    def equ_position(self, jd):
        return ln_get_venus_equ_coords(jd)

    def helio_position(self, jd):
        return ln_get_venus_helio_coords(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_venus_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst



class Earth(CelObject):
    def earth_distance(self, jd):
        return 0.0

    def solar_distance(self, jd):
        return ln_get_earth_solar_dist(jd)

    def phase(self, jd):
        return 0.0

    def magnitude(self, jd):
        return 0.0

    def riseset(self, jd):
        return None

    def equ_semidiameter(self, jd):
        return 0.0

    def pol_semidiameter(self, jd):
        return 0.0

    def illuminatedfraction(self, jd):
        return 0.0

    def equ_position(self, jd):
        return None

    def helio_position(self, jd):
        return None




class Mars(CelObject):
    def earth_distance(self, jd):
        return ln_get_mars_earth_dist(jd)

    def solar_distance(self, jd):
        return ln_get_mars_solar_dist(jd)

    def phase(self, jd):
        return ln_get_mars_phase(jd)

    def magnitude(self, jd):
        return ln_get_mars_magnitude(jd)

    def equ_semidiameter(self, jd):
        return ln_get_mars_sdiam(jd)

    def pol_semidiameter(self, jd):
        return ln_get_mars_sdiam(jd)

    def illuminatedfraction(self, jd):
        return ln_get_mars_disk(jd)

    def equ_position(self, jd):
        return ln_get_mars_equ_coords(jd)

    def helio_position(self, jd):
        return ln_get_mars_helio_coords(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_mars_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst


class Jupiter(CelObject):
    def earth_distance(self, jd):
        return ln_get_jupiter_earth_dist(jd)

    def solar_distance(self, jd):
        return ln_get_jupiter_solar_dist(jd)

    def phase(self, jd):
        return ln_get_jupiter_phase(jd)

    def magnitude(self, jd):
        return ln_get_jupiter_magnitude(jd)

    def equ_semidiameter(self, jd):
        return ln_get_jupiter_equ_sdiam(jd)

    def pol_semidiameter(self, jd):
        return ln_get_jupiter_pol_sdiam(jd)

    def illuminatedfraction(self, jd):
        return ln_get_jupiter_disk(jd)

    def equ_position(self, jd):
        return ln_get_jupiter_equ_coords(jd)

    def helio_position(self, jd):
        return ln_get_jupiter_helio_coords(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_jupiter_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst


class Saturn(CelObject):
    def earth_distance(self, jd):
        return ln_get_saturn_earth_dist(jd)

    def solar_distance(self, jd):
        return ln_get_saturn_solar_dist(jd)

    def phase(self, jd):
        return ln_get_saturn_phase(jd)

    def magnitude(self, jd):
        return ln_get_saturn_magnitude(jd)

    def equ_semidiameter(self, jd):
        return ln_get_saturn_equ_sdiam(jd)

    def pol_semidiameter(self, jd):
        return ln_get_saturn_pol_sdiam(jd)

    def illuminatedfraction(self, jd):
        return ln_get_saturn_disk(jd)

    def equ_position(self, jd):
        return ln_get_saturn_equ_coords(jd)

    def helio_position(self, jd):
        return ln_get_saturn_helio_coords(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_saturn_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst


class Uranus(CelObject):
    def earth_distance(self, jd):
        return ln_get_uranus_earth_dist(jd)

    def solar_distance(self, jd):
        return ln_get_uranus_solar_dist(jd)

    def phase(self, jd):
        return ln_get_uranus_phase(jd)

    def magnitude(self, jd):
        return ln_get_uranus_magnitude(jd)

    def equ_semidiameter(self, jd):
        return ln_get_uranus_sdiam(jd)

    def pol_semidiameter(self, jd):
        return ln_get_uranus_sdiam(jd)

    def illuminatedfraction(self, jd):
        return ln_get_uranus_disk(jd)

    def equ_position(self, jd):
        return ln_get_uranus_equ_coords(jd)

    def helio_position(self, jd):
        return ln_get_uranus_helio_coords(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_uranus_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst


class Neptune(CelObject):
    def earth_distance(self, jd):
        return ln_get_neptune_earth_dist(jd)

    def solar_distance(self, jd):
        return ln_get_neptune_solar_dist(jd)

    def phase(self, jd):
        return ln_get_neptune_phase(jd)

    def magnitude(self, jd):
        return ln_get_neptune_magnitude(jd)

    def equ_semidiameter(self, jd):
        return ln_get_neptune_sdiam(jd)

    def pol_semidiameter(self, jd):
        return ln_get_neptune_sdiam(jd)

    def illuminatedfraction(self, jd):
        return ln_get_neptune_disk(jd)

    def equ_position(self, jd):
        return ln_get_neptune_equ_coords(jd)

    def helio_position(self, jd):
        return ln_get_neptune_helio_coords(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_neptune_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst


class Pluto(CelObject):
    def earth_distance(self, jd):
        return ln_get_pluto_earth_dist(jd)

    def solar_distance(self, jd):
        return ln_get_pluto_solar_dist(jd)

    def phase(self, jd):
        return ln_get_pluto_phase(jd)

    def magnitude(self, jd):
        return ln_get_pluto_magnitude(jd)

    def equ_semidiameter(self, jd):
        return ln_get_pluto_sdiam(jd)

    def pol_semidiameter(self, jd):
        return ln_get_pluto_sdiam(jd)

    def illuminatedfraction(self, jd):
        return ln_get_pluto_disk(jd)

    def equ_position(self, jd):
        return ln_get_pluto_equ_coords(jd)

    def helio_position(self, jd):
        return ln_get_pluto_helio_coords(jd)

    def riseset(self, jd):
        rst = ln_lnlat_pos()
        env = Environment()
        (circ, rst) = ln_get_pluto_rst(jd, env.location)
        if circ:
            raise IsCircumpolar()
        return rst

solsystem = [ "Mercury", "Venus", "Earth", "Mars" , "Jupiter", "Saturn", "Uranus", "Neptune",  "Pluto" ]

solsystem_creator = {
 "Mercury" : Mercury(),
 "Venus" : Venus(),
 "Earth" : Earth(),
 "Mars"  : Mars(),
 "Jupiter" : Jupiter(),
 "Saturn"  : Saturn(),
 "Uranus"  : Uranus(),
 "Neptune" : Neptune(),
 "Pluto"   : Pluto(),
}

def lookup(obj : str) -> str:
    tmp = solsystem_creator[obj]
    return tmp

class ln_plotter:
    def __init__(self, headline = "", xlab = "", ylab = ""):
        self.data = {}
        self.df = None
        plt.xlabel(xlab)
        plt.ylabel(ylab)

    def addseries(self, key : str, data : list, indices : list):
       self.data[key] = pd.Series(data, index=indices, dtype=float)

    def show(self):
        if self.df is None:
            self.df = pd.DataFrame(self.data)
        self.df.plot()
        plt.show()


class AstException(Exception):
    pass

class UnimplementedMethod(AstException):
    pass

class IsCircumpolar(AstException):
    pass
