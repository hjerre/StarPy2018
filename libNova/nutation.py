from libNova.ln_types import *
from libNova.utility import *
from math import cos,sin, fabs

TERMS = 63
LN_NUTATION_EPOCH_THRESHOLD = 0.1

# arguments and coefficients taken from table 21A on page 133 */

arguments = [
    [0.0,	0.0,	0.0,	0.0,	1.0],
    [-2.0,	0.0,	0.0,	2.0,	2.0],
    [0.0,	0.0,	0.0,	2.0,	2.0],
    [0.0,	0.0,	0.0,	0.0,	2.0],
    [0.0,	1.0,	0.0,	0.0,	0.0],
    [0.0,	0.0,	1.0,	0.0,	0.0],
    [-2.0,	1.0,	0.0,	2.0,	2.0],
    [0.0,	0.0,	0.0,	2.0,	1.0],
    [0.0,	0.0,	1.0,	2.0,	2.0],
    [-2.0,	-1.0,	0.0,	2.0,	2.0],
    [-2.0,	0.0,	1.0,	0.0,	0.0],
    [-2.0,	0.0,	0.0,	2.0,	1.0],
    [0.0,	0.0,	-1.0,	2.0,	2.0],
    [2.0,	0.0,	0.0,	0.0,	0.0],
    [0.0,	0.0,	1.0,	0.0,	1.0],
    [2.0,	0.0,	-1.0,	2.0,	2.0],
    [0.0,	0.0,	-1.0,	0.0,	1.0],
    [0.0,	0.0,	1.0,	2.0,	1.0],
    [-2.0,	0.0,	2.0,	0.0,	0.0],
    [0.0,	0.0,	-2.0,	2.0,	1.0],
    [2.0,	0.0,	0.0,	2.0,	2.0],
    [0.0,	0.0,	2.0,	2.0,	2.0],
    [0.0,	0.0,	2.0,	0.0,	0.0],
    [-2.0,	0.0,	1.0,	2.0,	2.0],
    [0.0,	0.0,	0.0,	2.0,	0.0],
    [-2.0,	0.0,	0.0,	2.0,	0.0],
    [0.0,	0.0,	-1.0,	2.0,	1.0],
    [0.0,	2.0,	0.0,	0.0,	0.0],
    [2.0,	0.0,	-1.0,	0.0,	1.0],
    [-2.0,	2.0,	0.0,	2.0,	2.0],
    [0.0,	1.0,	0.0,	0.0,	1.0],
    [-2.0,	0.0,	1.0,	0.0,	1.0],
    [0.0,	-1.0,	0.0,	0.0,	1.0],
    [0.0,	0.0,	2.0,	-2.0,	0.0],
    [2.0,	0.0,	-1.0,	2.0,	1.0],
    [2.0,	0.0,	1.0,	2.0,	2.0],
    [0.0,	1.0,	0.0,	2.0,	2.0],
    [-2.0,	1.0,	1.0,	0.0,	0.0],
    [0.0,	-1.0,	0.0,	2.0,	2.0],
    [2.0,	0.0,	0.0,	2.0,	1.0],
    [2.0,	0.0,	1.0,	0.0,	0.0],
    [-2.0,	0.0,	2.0,	2.0,	2.0],
    [-2.0,	0.0,	1.0,	2.0,	1.0],
    [2.0,	0.0,	-2.0,	0.0,	1.0],
    [2.0,	0.0,	0.0,	0.0,	1.0],
    [0.0,	-1.0,	1.0,	0.0,	0.0],
    [-2.0,	-1.0,	0.0,	2.0,	1.0],
    [-2.0,	0.0,	0.0,	0.0,	1.0],
    [0.0,	0.0,	2.0,	2.0,	1.0],
    [-2.0,	0.0,	2.0,	0.0,	1.0],
    [-2.0,	1.0,	0.0,	2.0,	1.0],
    [0.0,	0.0,	1.0,	-2.0,	0.0],
    [-1.0,	0.0,	1.0,	0.0,	0.0],
    [-2.0,	1.0,	0.0,	0.0,	0.0],
    [1.0,	0.0,	0.0,	0.0,	0.0],
    [0.0,	0.0,	1.0,	2.0,	0.0],
    [0.0,	0.0,	-2.0,	2.0,	2.0],
    [-1.0,	-1.0,	1.0,	0.0,	0.0],
    [0.0,	1.0,	1.0,	0.0,	0.0],
    [0.0,	-1.0,	1.0,	2.0,	2.0],
    [2.0,	-1.0,	-1.0,	2.0,	2.0],
    [0.0,	0.0,	3.0,	2.0,	2.0],
    [2.0,	-1.0,	0.0,	2.0,	2.0]
]

coefficients = [
    [-171996.0,	-174.2,	92025.0,8.9],
    [-13187.0,	-1.6,  	5736.0,	-3.1],
    [-2274.0, 	-0.2,  	977.0,	-0.5],
    [2062.0,   	0.2,    -895.0,    0.5],
    [1426.0,    -3.4,    54.0,    -0.1],
    [712.0,    0.1,    -7.0,    0.0],
    [-517.0,    1.2,    224.0,    -0.6],
    [-386.0,    -0.4,    200.0,    0.0],
    [-301.0,    0.0,    129.0,    -0.1],
    [217.0,    -0.5,    -95.0,    0.3],
    [-158.0,    0.0,    0.0,    0.0],
    [129.0,	0.1,	-70.0,	0.0],
    [123.0,	0.0,	-53.0,	0.0],
    [63.0,	0.0,	0.0,	0.0],
    [63.0,	0.1,	-33.0,	0.0],
    [-59.0,	0.0,	26.0,	0.0],
    [-58.0,	-0.1,	32.0,	0.0],
    [-51.0,	0.0,	27.0,	0.0],
    [48.0,	0.0,	0.0,	0.0],
    [46.0,	0.0,	-24.0,	0.0],
    [-38.0,	0.0,	16.0,	0.0],
    [-31.0,	0.0,	13.0,	0.0],
    [29.0,	0.0,	0.0,	0.0],
    [29.0,	0.0,	-12.0,	0.0],
    [26.0,	0.0,	0.0,	0.0],
    [-22.0,	0.0,	0.0,	0.0],
    [21.0,	0.0,	-10.0,	0.0],
    [17.0,	-0.1,	0.0,	0.0],
    [16.0,	0.0,	-8.0,	0.0],
    [-16.0,	0.1,	7.0,	0.0],
    [-15.0,	0.0,	9.0,	0.0],
    [-13.0,	0.0,	7.0,	0.0],
    [-12.0,	0.0,	6.0,	0.0],
    [11.0,	0.0,	0.0,	0.0],
    [-10.0,	0.0,	5.0,	0.0],
    [-8.0,	0.0,	3.0,	0.0],
    [7.0,	0.0,	-3.0,	0.0],
    [-7.0,	0.0,	0.0,	0.0],
    [-7.0,	0.0,	3.0,	0.0],
    [-7.0,	0.0,	3.0,	0.0],
    [6.0,	0.0,	0.0,	0.0],
    [6.0,	0.0,	-3.0,	0.0],
    [6.0,	0.0,	-3.0,	0.0],
    [-6.0,	0.0,	3.0,	0.0],
    [-6.0,	0.0,	3.0,	0.0],
    [5.0,	0.0,	0.0,	0.0],
    [-5.0,	0.0,	3.0,	0.0],
    [-5.0,	0.0,	3.0,	0.0],
    [-5.0,	0.0,	3.0,	0.0],
    [4.0,	0.0,	0.0,	0.0],
    [4.0,	0.0,	0.0,	0.0],
    [4.0,	0.0,	0.0,	0.0],
    [-4.0,	0.0,	0.0,	0.0],
    [-4.0,	0.0,	0.0,	0.0],
    [-4.0,	0.0,	0.0,	0.0],
    [3.0,	0.0,	0.0,	0.0],
    [-3.0,	0.0,	0.0,	0.0],
    [-3.0,	0.0,	0.0,	0.0],
    [-3.0,	0.0,	0.0,	0.0],
    [-3.0,	0.0,	0.0,	0.0],
    [-3.0,	0.0,	0.0,	0.0],
    [-3.0,	0.0,	0.0,	0.0],
    [-3.0,	0.0,	0.0,	0.0]
]


# def ln_get_nutation( JD) -> ln_nutation
# param JD Julian Day.
# returns nutation Pointer to store nutation
# Calculate nutation of longitude and obliquity in degrees from Julian Ephemeris Day


def ln_get_nutation( JD : float ) -> ln_nutation:
    nutation = ln_nutation()
    c_longitude = 0.0
    c_obliquity = 0.0
    c_ecliptic = 0.0


    # get julian ephemeris day
    JDE = float( JD )

    # calc T
    T = (JDE - 2451545.0) / 36525.0
    T2 = T * T
    T3 = T2 * T

    # calculate D,M,M',F and Omega
    D = 297.85036 + 445267.111480 * T - 0.0019142 * T2 + T3 / 189474.0
    M = 357.52772 + 35999.050340 * T - 0.0001603 * T2 - T3 / 300000.0
    MM = 134.96298 + 477198.867398 * T + 0.0086972 * T2 + T3 / 56250.0
    F = 93.2719100 + 483202.017538 * T - 0.0036825 * T2 + T3 / 327270.0
    O = 125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000.0

    # convert to radians
    D = ln_deg_to_rad(D)
    M = ln_deg_to_rad(M)
    MM = ln_deg_to_rad(MM)
    F = ln_deg_to_rad(F)
    O = ln_deg_to_rad(O)


    for i in range(0, len(coefficients)):
        coeff_sine = (coefficients[i][0] + (coefficients[i][1] * T))
        coeff_cos = (coefficients[i][2] + (coefficients[i][3] * T))

        argument = arguments[i][0] * D + arguments[i][1] * M + arguments[i][2] * MM + arguments[i][3] * F + arguments[i][4] * O

        c_longitude += coeff_sine * sin(argument)
        c_obliquity += coeff_cos * cos(argument)


    # change to arcsecs
    c_longitude /= 10000.0
    c_obliquity /= 10000.0

    # change to degrees
    c_longitude /= (60.0 * 60.0)
    c_obliquity /= (60.0 * 60.0)

    # calculate mean ecliptic - Meeus 2nd edition, eq. 22.2
    c_ecliptic = 23.0 + 26.0 / 60.0 + 21.448 / 3600.0 - 46.8150 / 3600.0 * T - 0.00059 / 3600.0 * T2 + 0.001813 / 3600.0 * T3

    # c_ecliptic += c_obliquity;

    # return results
    nutation.longitude = c_longitude
    nutation.obliquity = c_obliquity
    nutation.ecliptic = c_ecliptic
    return nutation


# def ln_get_equ_nut( mean_position, JD) -> ln_equ_posn
# param mean_position Mean position of object
# param JD Julian Day.
# returns position Pointer to store new object position.

# Calculate a stars equatorial coordinates from it's mean equatorial coordinates with the effects of nutation for a given Julian Day.

def ln_get_equ_nut( mean_position : ln_equ_posn, JD : float ) -> ln_equ_posn:
    posn = ln_equ_posn()
    nut = ln_get_nutation(JD)

    mean_ra = ln_deg_to_rad(mean_position.ra)
    mean_dec = ln_deg_to_rad(mean_position.dec)

    nut_ecliptic = ln_deg_to_rad(nut.ecliptic + nut.obliquity)
    sin_ecliptic = sin(nut_ecliptic)

    sin_ra = sin(mean_ra)
    cos_ra = cos(mean_ra)

    tan_dec = tan(mean_dec)

    delta_ra = (cos (nut_ecliptic) + sin_ecliptic * sin_ra * tan_dec) * nut.longitude - cos_ra * tan_dec * nut.obliquity;
    delta_dec = (sin_ecliptic * cos_ra) * nut.longitude + sin_ra * nut.obliquity;

    posn.ra = mean_position.ra + delta_ra
    posn.dec = mean_position.dec + delta_dec
    return posn
