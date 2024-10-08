# Autogenerated with SMOP 0.29
# Modified by Adam
import numpy as np


def geod2geoc(alpha, h, B=None):
    # [r, theta] = geod2geoc(alpha, h);
    # [r, theta, B_r, B_theta] = geod2geoc(alpha, h, X, Z);
    # conversion from geodetic X,Z components to geocentric B_r, B_theta
    # Input:   geodetic latitude alpha (rad)
    #          altitude h [km]
    #          X, Z
    # Output:  theta (rad)
    #          r (km)
    #          B_r, B_theta

    # Nils Olsen, DSRI Copenhagen, September 2001.
    # After Langel (1987), eq. (52), (53), (56), (57)
    # Converted to python by Adam (2017)

    # Latest changes: NIO, 170903 optional conversion of magnetic components
    #                 Update to WGS-84 ellipsoid

    # WGS-84 Ellipsoid parameters
    a = 6378.137
    b = 6356.7523142
    sin_alpha_2 = np.sin(alpha) ** 2
    cos_alpha_2 = np.cos(alpha) ** 2
    tmp = np.multiply(
        h,
        np.sqrt(np.dot(a ** 2, cos_alpha_2) + np.dot(b ** 2, sin_alpha_2))
    )
    beta = np.arctan(
        np.multiply((tmp + b ** 2) / (tmp + a ** 2), np.tan(alpha))
    )
    theta = np.pi / 2 - beta
    r = np.sqrt(
        h ** 2
        + np.dot(2, tmp)
        + np.dot(
            a ** 2, (1 - np.dot((1 - (b / a) ** 4), sin_alpha_2))
        ) / (1 - np.dot((1 - (b / a) ** 2), sin_alpha_2))
    )

    # convert also magnetic components from geodetic to B_r, B_theta
    if B is not None:
        psi = np.multiply(
            np.sin(alpha),
            np.sin(theta)
        ) - np.multiply(
            np.cos(alpha),
            np.cos(theta)
        )
        X = np.copy(B[0])
        Z = np.copy(B[1])
        B[0] = np.multiply(- np.sin(psi), X) - np.multiply(np.cos(psi), Z)
        B[1] = np.multiply(- np.cos(psi), X) + np.multiply(np.sin(psi), Z)
        return r, theta, B
    else:
        return r, theta
