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

def getSQfield(lat: Union[float, list], lon: Union[float, list], year: Union[int, list], month: Union[int, list], day: Union[int, list], hour: Union[int, list]=0, minutes: Union[int, list]=0, h: Union[float, list]=0,r: Union[float, list]=0, f107_1: Optional[Union[float, list]]=None, model_name: Optional[Union[str]]="xdifi2", geoc:Optional[bool] = False, return_geoc:Optional[bool] = False) -> dict:
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

    earth_radius_km = 6371.2
    
    #forward_Sq_d_Re implicitly expects flat numpy arrays as inputs
    #some inputs are converted to equivalent numpy arrays by the
    #coordinate or time transformations, but there a few which 
    #need to be explicitly cast:
    lat = np.array(lat).flatten()
    lon = np.array(lon).flatten()
    h = np.array(h).flatten()
    if not geoc:
        r_gc, theta_gc = util.geod_to_geoc_lat(lat, h)
    else:
        theta_gc = lat
        r_gc = r

    
    # r_gc = earth_radius_km + h
    cotheta_gc = 90 - theta_gc

    start_time = 0#5114.0
    # print("trying to use model_name", model_name, model_name.lower())
    # end_time = float(get_f107_index.difi_t_f107[-1])
    sq_t = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)

    if (model_name.lower() == 'xdifi2'):
        """("importing coeff from xdifi")"""
        # from DIFI import get_f107_index_xDIFI as get_f107_index
        from DIFI.get_f107_index_all import get_f107_index, load_coefs, load_swarm_xDIFI
        difi_t_f107, difi_f107 = load_coefs()
        swarm_data = load_swarm_xDIFI()
        end_time = float(difi_t_f107[-1])
    elif (model_name.lower() == 'difi8'):
        """("importing coeff from difi8")"""
        from DIFI.get_f107_index_all import get_f107_index, load_coefs, load_swarm_DIFI8
        difi_t_f107, difi_f107 = load_coefs()
        swarm_data = load_swarm_DIFI8()
        end_time = float(difi_t_f107[-1])
    elif (model_name.lower() == 'difi7'):
        """("importing coeff from difi8")"""
        from DIFI.get_f107_index_all import get_f107_index, load_coefs, load_swarm_DIFI7
        difi_t_f107, difi_f107 = load_coefs()
        swarm_data = load_swarm_DIFI7()
        end_time = float(difi_t_f107[-1])
    else:
        raise ValueError("Input model_name = difi8, xdifi2, or difi7. Input model name didn't match any models.")
    
    if f107_1 is None:
        f107_1 = get_f107_index(sq_t, start_time, end_time, difi_f107, difi_t_f107)
    else:
        f107_1 = np.array(f107_1).flatten()
        
    B_XYZ = {}
    # print "Difi input", r_gc, theta_gc, RV['lon'], sq_t, f107_1
    
    [B_1, B_2] = forward_Sq_d_Re.forward_Sq_d_Re(
        r_gc,
        cotheta_gc,
        lon,
        sq_t,
        f107_1,
        swarm_data
    )
    B_C = B_1 + B_2
    B_XYZ['Z'] = -1 * B_C[0]
    B_XYZ['Y'] = B_C[2]
    B_XYZ['X'] = -1 * B_C[1]

    if not return_geoc:
        if geoc:
            lat = theta_to_geod_lat(lat)
        Bx, By, Bz = magmath.rotate_magvec(B_XYZ['X'], B_XYZ['Y'], B_XYZ['Z'], theta_gc, lat)

        B_XYZ["X"] = Bx
        B_XYZ["Y"] = By
        B_XYZ["Z"] = Bz
    else:
        print("B_XYZ values are returning in geocentric (NEC) coord")
    return B_XYZ