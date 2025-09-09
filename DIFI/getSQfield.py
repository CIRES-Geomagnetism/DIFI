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

def getSQfield(lat: Union[float, list], lon: Union[float, list], year: Union[int, list], month: Union[int, list],
                day: Union[int, list], hour: Union[int, list]=0, minutes: Union[int, list]=0, 
                h: Union[float, list]=0,r: Union[float, list]=0, f107_1: Optional[Union[float, list]]=None, 
                model_name: Optional[Union[str]]="xdifi2", geoc:Optional[bool] = False, 
                return_geoc:Optional[bool] = False) -> dict:
    """
    Input:
        Latitude, lat (in WGS-84 coordinates)
        Longtitude, lon
        An array of year, year (Only good between 2000.0 and 2025.0
        An array of months, month
        An array of days, day

    Optional Input:
        An array of hours, hour
        An array of minutes, minutes
        Height above WGS84 ellipsoid, h
        Radius in geocentric coordinates, r
        Name of DIFI model, DIFI8,xDIFI2, DIFI7, model_name
        Treat inputs as geocentric (must use r input, not h), geoc
        Return geocentric B, return_geoc

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
    
    if not geoc: #input `lat` is assumed geo*detic*, input `h` is used, input `r` is ignored
        r_gc, theta_gc = util.geod_to_geoc_lat(lat, h) #convert geodetic latitude and ellipsoidal height to geocentric lat and radius (distance from center of earth to your location)
    else: #input `lat` is assumed be geo*centric* latitude, input `r` is used as geocentric radius (distance from center of earth to your location), input `h` is ignored
        theta_gc = lat
        r_gc = r

    if np.any(np.abs((r_gc - 6371.2) - 110) < 30):
        warnings.warn("Warning: The altitude is within 30 km of 110 km, where the model switches between internal and external ionospheric current sources. The model is not expected to accurately represent the field within this transition zone.")

    cotheta_gc = 90 - theta_gc

    start_time = 0#5114.0

    sq_t = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)
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

    if (model_name.lower() == 'xdifi2'):
        """("importing coeff from xdifi")"""
        # from DIFI import get_f107_index_xDIFI as get_f107_index
        from DIFI.get_f107_index_all import get_f107_index, load_coefs, load_swarm_xDIFI
        difi_t_f107, difi_f107 = load_coefs()
        swarm_data = load_swarm_xDIFI()
        end_time = float(difi_t_f107[-1])
        if np.any(sq_t < warn_year_2001):
            warnings.warn("Dataset contains date before 2001.0, outside xDIFI2's recommended range", UserWarning)
        if np.any(sq_t >= warn_year_2024):

            warnings.warn("Dataset contains date after 2024.0, outside xDIFI2's recommended range", UserWarning)
    elif (model_name.lower() == 'difi8'):
        """("importing coeff from difi8")"""
        from DIFI.get_f107_index_all import get_f107_index, load_coefs, load_swarm_DIFI8
        difi_t_f107, difi_f107 = load_coefs()
        swarm_data = load_swarm_DIFI8()
        end_time = float(difi_t_f107[-1])
        if np.any(sq_t < warn_year_2014):
            warnings.warn("Dataset contains date before 2014.0, outside DIFI8's recommended range")
        if np.any(sq_t >= warn_year_2024):
            warnings.warn("Dataset contains date after 2024.0, outside DIFI8's recommended range")

    elif (model_name.lower() == 'difi7'):
        raise ValueError("Difi 7 is not supported since the release of difi8")
        """("importing coeff from difi8")"""
        from DIFI.get_f107_index_all import get_f107_index, load_coefs, load_swarm_DIFI7
        difi_t_f107, difi_f107 = load_coefs()
        swarm_data = load_swarm_DIFI7()
        end_time = float(difi_t_f107[-1])
    else:
        raise ValueError("Input model_name = difi8, xdifi2, or difi7. Input model name didn't match any models.")
    
    if f107_1 is None:
        f107_1 = get_f107_index(sq_t, start_f107_time, end_f107_time, difi_f107, difi_t_f107)
    else:
        f107_1 = np.array(f107_1).flatten()
        
    B_XYZ = {}
    # print "Difi input", r_gc, theta_gc, RV['lon'], sq_t, f107_1
    
    # calculate radii in units of reference radius
    a = 6371.2
    rho_Sq = (a + swarm_data['h']) / a
    rho = np.array(r_gc / a)
    if (np.min(rho) < rho_Sq) and (np.max(rho) >= rho_Sq):#If min below SQ and max above SQ
        #Split dataset into above and below SQ
        above_SQ_mask = rho < rho_Sq
        below_SQ_mask = rho >= rho_Sq
        # print('rioiarsnt\n\n\n',r_gc[above_SQ_mask], 
        # cotheta_gc[above_SQ_mask], lon[above_SQ_mask], sq_t[above_SQ_mask], f107_1[above_SQ_mask], 
        # swarm_data)
        [above_output1, above_output2] = forward_Sq_d_Re.forward_Sq_d_Re(r_gc[above_SQ_mask], 
        cotheta_gc[above_SQ_mask], lon[above_SQ_mask], sq_t[above_SQ_mask], f107_1[above_SQ_mask], 
        swarm_data)
        [below_output1, below_output2] = forward_Sq_d_Re.forward_Sq_d_Re(r_gc[below_SQ_mask], 
        cotheta_gc[below_SQ_mask], lon[below_SQ_mask], sq_t[below_SQ_mask], f107_1[below_SQ_mask], 
        swarm_data)
        B_1,B_2 = np.zeros((3, len(r_gc))), np.zeros((3, len(r_gc)))
        B_1[:, above_SQ_mask] = above_output1
        B_1[:, below_SQ_mask] = below_output1

        B_2[:, above_SQ_mask] = above_output2
        B_2[:, below_SQ_mask] = below_output2
    else:
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