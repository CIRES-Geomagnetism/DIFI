import warnings

import numpy as np
import os

from DIFI import jd2000_dt
from DIFI import forward_Sq_d_Re
# from DIFI import get_f107_index
from geomaglib import util, magmath
from typing import Optional, Union


def theta_to_geod_lat(theta: float) -> float:

    f = 1 / 298.257223563
    inv_f = 1 / (1 - f)
    lat = np.arctan(inv_f * inv_f * np.tan(np.radians(theta)))
    lat = np.degrees(lat)

    return lat

if __name__ == '__main__':
    
    start_time = 0#5114.0
    year = [2000,2001,2025]
    month = [1,1,12]
    day = [1,1,31]
    hour = [1,12,23]
    minutes = [1,0,59]
    sq_t = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)
    # print(sq_t)
    warn_year_2001 = jd2000_dt.jd2000_dt(2001, 1, 1, 0, 0)
    warn_year_2014= jd2000_dt.jd2000_dt(2014, 1, 1, 0, 0)    
    warn_year_2024= jd2000_dt.jd2000_dt(2024, 1, 1, 0, 0)
    end_DIFI_valid = 9496.5

    if np.any(sq_t) > end_DIFI_valid:
        raise Exception(
            "DIFI is not valid after noon 12/31/2026 . Input time data contains a date outside DIFI's validity range."
        )
    start_f107_time = jd2000_dt.jd2000_dt(2000, 1, 1, 0, 0)
    end_f107_time =  jd2000_dt.jd2000_dt(2026, 1, 1, 12, 30)


    """("importing coeff from xdifi")"""
    # from DIFI import get_f107_index_xDIFI as get_f107_index
    from DIFI.get_f107_index_all import get_f107_index, load_coefs, load_swarm_xDIFI
    difi_t_f107, difi_f107 = load_coefs()
    swarm_data = load_swarm_xDIFI()
    end_time = float(difi_t_f107[-1])

    f107_1 = get_f107_index(sq_t, start_f107_time, end_f107_time, difi_f107, difi_t_f107)

    # print(f107_1)
    # sq_t = np.linspace(3900, 4000, 10000)
    # f107_1 = get_f107_index(sq_t, start_f107_time, end_f107_time, difi_f107, difi_t_f107)
    # import matplotlib.pyplot as plt
    # plt.plot(sq_t, f107_1)
    # plt.xlabel('days/jd2000 time')
    # plt.ylabel('f107 magnitude')
    # plt.grid()
    # plt.show()
    assert(np.abs(np.sum(np.array(f107_1) - np.array([129.9 ,171 ,222.24])) < 1e-15))
    print("f107 was read properly")