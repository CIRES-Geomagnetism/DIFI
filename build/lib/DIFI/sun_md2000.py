import math
import numpy as np
from typing import Union

# Celestial Computing with MATLAB but with Python
# input
# mjd2000 = modified julian day 2000
# output
# rasc = right ascension of the sun (radians)
#   (0 <= rasc <= 2 pi)
# decl = declination of the sun (radians)
#   (-pi/2 <= decl <= pi/2)
# rsun = eci position vector of the sun (km)
# note
#  coordinates are inertial, geocentric,
#  equatorial and true-of-date

# Modified by Nils Oslen, DSRI
# Translated to python by Adam Woods, NCEI


def sun_md2000(mjd2000: Union[float, list]) -> list[np.ndarray, np.ndarray]:
    atr = math.pi/648000
    # time arguments
    djd = mjd2000 - 0.5
    t = djd/36525.0 + 1

    # fundamental arguments (converted from revolutions to radians with r2r)
    gs = r2r(0.993126 + 0.0027377785 * djd)
    lm = r2r(0.606434 + 0.03660110129 * djd)
    ls = r2r(0.779072 + 0.00273790931 * djd)
    g2 = r2r(0.140023 + 0.00445036173 * djd)
    g4 = r2r(0.053856 + 0.00145561327 * djd)
    g5 = r2r(0.056531 + 0.00023080893 * djd)
    rm = r2r(0.347343 - 0.00014709391 * djd)

    # geocentric, ecliptic longitude of the sun (radians)
    plon = 6910 * np.sin(gs) + 72 * np.sin(2 * gs) - 17 * t * np.sin(gs)
    plon = plon - 7 * np.cos(gs - g5) + 6 * np.sin(lm - ls) \
        + 5 * np.sin(4 * gs - 8 * g4 + 3 * g5)
    plon = plon - 5 * np.cos(2 * (gs - g2)) - 4 * (
        np.sin(gs - g2)
        - np.cos(4 * gs - 8 * g4 + 3 * g5)
    )
    plon = plon + 3 * (
        np.sin(2 * (gs - g2)) - np.sin(g5) - np.sin(2 * (gs - g5))
    )
    plon = ls + atr * (plon - 17 * np.sin(rm))
    # plon%2pi*180/pi

    # geocentric distance of the sun (kilometers)
    # rsm = 149597870.691 * (
    #     1.00014 - 0.01675 * np.cos(gs) - 0.00014 * np.cos(2 * gs)
    # )

    # obliquity of the ecliptic (radians)
    obliq = atr * (84428 - 47 * t + 9 * np.cos(rm))

    # geocentric, equatorial right ascension and declination (radians)
    a = np.sin(plon)*np.cos(obliq)
    b = np.cos(plon)

    rasc = np.arctan2(a, b)
    decl = np.arcsin(np.sin(obliq)*np.sin(plon))

    return [rasc, decl]


def r2r(x: float) -> float:
    # revolutions to radians function
    # input
    #  x = argument(revolutions; 0 <= x <= 1)
    # returns
    #  equivalent x (radians; 0 <= y <= 2 pi)

    return 2 * math.pi * (x % 1)
