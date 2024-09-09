import numpy as np

from DIFI import difi_t_f107, difi_f107
from DIFI import SwarmL2_F107_Read
from DIFI import jd2000_dt
from DIFI import difi_f107, difi_t_f107, swarm_data
from DIFI import geod2geoc
from DIFI import forward_Sq_d_Re

def getSQfield(lat, lon, year, month, day, hour=0, minutes=0, h=0):
    """
    Input:
        Latitude, lat (in WGS-84 coordinates)
        Longtitude, lon
        An array of year, year (Only good between 2014.0 and 2025.0
        An array of months, month
        An array of days, day

    Optional Input:
        An array of hours, hour
        An array of minutes, minutes
        Height above WGS84 ellipsoid, h

    Output:
        B, the magnetic field due to the SQ in WGS-84 coordinates
    """
    a = 6371.2
    start_time = 5114.0
    end_time = float(difi_t_f107[-1])
    sq_t = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)
    frac_arr = sq_t - np.floor(sq_t)
    f107_1 = np.array([])
    for i in range(np.size(sq_t)):
        if sq_t[i] < start_time:
            raise Exception(
                "Request before 2014.0 reached difi calculation improper"
            )
        elif sq_t[i] > end_time:
            raise Exception(
                "Request after 2024.0 reached difi calculation improperly"
            )
        while sq_t[i] < 5114.0:
            sq_t[i] += 365
        j = 0
        while difi_t_f107[j] < sq_t[i]:
            j += 1
        f107_1 = np.append(
            f107_1,
            difi_f107[j] * frac_arr[i] + difi_f107[j - 1] * (1 - frac_arr[i])
        )

    B_XYZ = {}
    if h < 20:  # LIMIT ALTITUDE REQUESTS TO 20 KM.  Model invalid above!
        [r_gc, theta_gc] = geod2geoc.geod2geoc(np.radians(lat), h)
        r_gc = a + h
        theta_gc = np.degrees(theta_gc)
        # print "Difi input", r_gc, theta_gc, RV['lon'], sq_t, f107_1
        [B_1, B_2] = forward_Sq_d_Re.forward_Sq_d_Re(
            r_gc,
            theta_gc,
            lon,
            sq_t,
            f107_1,
            swarm_data,
        )
        # print "Difi output", B_1, B_2
        B_C = B_1 + B_2
        # B_C = B_1
        B_XYZ['Z'] = -1 * B_C[0]
        B_XYZ['Y'] = B_C[2]
        B_XYZ['X'] = -1 * B_C[1]
    else:
        raise ValueError(
            "Requested altitude {h} km is higher than 20 km".format(h=repr(h))
        )
        B_XYZ['Z'] = 0
        B_XYZ['Y'] = 0
        B_XYZ['X'] = 0
    B = RotateMagneticVector(B_XYZ, 90 - theta_gc, lat)
    # B = B_XYZ
    return B, f107_1