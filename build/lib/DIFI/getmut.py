import math
from typing import Union

import numpy as np

from DIFI import sun_md2000
def getmut(mjd2000: list, theta_NGP: Union[float, list], phi_NGP: float, *args, **kwargs) -> np.ndarray:
    # mut = getmut(mjd2000,theta_NGP,phi_NGP)

    # Calculate the Magnetic Universal Time for a given MJD2000 and a given
    # position of the North Geomagnetic Pole.

    # Input:    mjd2000         modified julian day 2000 (days)
    #           theta_NGP       geographic co-latitude of North Geomagnetic
    #                           Pole (deg)
    #           phi_NGP         geographic longitude of North Geomagnetic
    #                           Pole (deg)
    # Output:   mut             magnetic universal time (hr)

    # A. Chulliat, October 2008
    # Converted to python by A. Woods

    ######################################################################

    rad = math.pi / 180
    vec = np.ones(mjd2000.size)
    theta_0_rad = np.deg2rad(theta_NGP*vec)
    phi_0 = phi_NGP*vec
    ut = (mjd2000 - np.floor(mjd2000))*24
    phi_s = 15*(12*vec - ut)

    rasc, decl = sun_md2000.sun_md2000(mjd2000)
    theta_s = 90*vec - decl / rad

    theta_s_rad = np.deg2rad(theta_s)
    phi_diff_rad = np.deg2rad(phi_s - phi_0)

    eopp = np.sin(theta_s_rad)*np.sin(phi_diff_rad)
    eadj = np.cos(theta_0_rad)*np.sin(theta_s_rad)*np.cos(phi_diff_rad) \
        - np.sin(theta_0_rad)*np.cos(theta_s_rad)
    mut = 12*vec - np.arctan2(eopp, eadj) / rad / 15

    return mut
