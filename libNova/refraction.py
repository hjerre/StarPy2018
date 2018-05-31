from math import tan, sin
from libNova.utility import *

# def ln_get_refraction_adj( altitude, atm_pres, temp) -> float
# param altitude The altitude of the object above the horizon in degrees
# param atm_pres Atmospheric pressure in milibars
# param temp Temperature in degrees C.
# return Adjustment in objects altitude in degrees.

# Calculate the adjustment in altitude of a body due to atmosphric refraction. This value varies over altitude, pressure and temperature.
# Note: Default values for pressure and teperature are 1010 mBar and 10C respectively.

def ln_get_refraction_adj( altitude : float, atm_pres : float, temp : float ) -> float:
    # equ 16.3
    R = 1.0 / tan(ln_deg_to_rad(altitude + (7.31 / (altitude + 4.4))))
    R -= 0.06 * sin(ln_deg_to_rad(14.7 * (R / 60.0) + 13.0))

    # take into account of atm press and temp
    R *= ((atm_pres / 1010.0) * (283.0 / (273.0 + temp)))

    # convert from arcminutes to degrees
    R /= 60.0

    return R
