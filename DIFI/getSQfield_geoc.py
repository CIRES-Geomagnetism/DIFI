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
def getSQfield_geoc(theta: Union[float, list], r: Union[float, list], lon: Union[float, list], year: Union[int, list], month: Union[int, list], day: Union[int, list], hour: Union[int, list]=0, minutes: Union[int, list]=0, h: Union[float, list]=0, f107_1: Optional[Union[float, list]]=None) -> dict:
    """
    Input:
        Latitude, lat (-90,90)
        Longtitude, lon 
        An array of year, year (Only good between 2014.0 and 2025.0
        An array of months, month
        An array of days, day

    Optional Input:
        An array of hours, hour
        An array of minutes, minutes
        Height above sphere, h

    Output:
        B, the magnetic field due to the SQ in WGS-84 coordinates
    """

    earth_radius_km = 6371.2
    
    #forward_Sq_d_Re implicitly expects flat numpy arrays as inputs
    #some inputs are converted to equivalent numpy arrays by the
    #coordinate or time transformations, but there a few which 
    #need to be explicitly cast:
    theta = np.array(theta).flatten()
    r = np.array(r).flatten()
    h = np.array(h).flatten()
    
    theta_gc = theta
    #r_gc = earth_radius_km + h
    r_gc = r
    cotheta_gc = 90 - theta
    start_time = 5114.0
    
    # end_time = float(get_f107_index.difi_t_f107[-1])
    sq_t = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)
    if (sq_t < 5114.0):
        """("importing coeff from xdifi")"""
        from DIFI import get_f107_index_xDIFI as get_f107_index
        end_time = float(get_f107_index.difi_t_f107[-1])
    else:
        """("importing coeff from difi8")"""
        from DIFI import get_f107_index_DIFI8 as get_f107_index
        end_time = float(get_f107_index.difi_t_f107[-1])
    if f107_1 is None:
        f107_1 = get_f107_index.get_f107_index(sq_t, start_time, end_time)
    else:
        f107_1 = np.array(f107_1).flatten()
        
    B_XYZ = {}

    swarm_data_difi = get_f107_index.load_swarm_DIFI8()
    
    [B_1, B_2] = forward_Sq_d_Re.forward_Sq_d_Re(
        r_gc,
        cotheta_gc,
        lon,
        sq_t,
        f107_1,
        swarm_data_difi
    )
    B_C = B_1 + B_2
    B_XYZ['Z'] = -1 * B_C[0]
    B_XYZ['Y'] = B_C[2]
    B_XYZ['X'] = -1 * B_C[1]

    lat = theta_to_geod_lat(theta_gc)


    Bx, By, Bz = magmath.rotate_magvec(B_XYZ['X'], B_XYZ['Y'], B_XYZ['Z'], theta_gc, lat)

    B_XYZ["X"] = Bx
    B_XYZ["Y"] = By
    B_XYZ["Z"] = Bz

    return B_XYZ

def RotateMagneticVector(B_geoc, lat_geoc, latitude):
    # Rotate the Magnetic Vectors to Geodetic Coordinates
    # Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
    # Equation 16, WMM Technical report
    # INPUT : lat_geoc: geocentric latitude
    #                latitude:  geodetic latitude
    #
    #                B_geoc: Dictionary with the following elements
    #                        X      North
    #                        Y      East
    #                        Z      Down
    #
    # OUTPUT: MagneticResultsGeo Dictionary with the following elements
    #                        X      North
    #                        Y      East
    #                        Z      Down
    #
    # Difference between the spherical and Geodetic latitudes
    Psi = np.radians(lat_geoc - latitude)
    B = {}

    # Rotate spherical field components to the Geodetic system
    B['Z'] = B_geoc['X'] * np.sin(Psi) + B_geoc['Z'] * np.cos(Psi)
    B['X'] = B_geoc['X'] * np.cos(Psi) - B_geoc['Z'] * np.sin(Psi)
    B['Y'] = B_geoc['Y']
    return B