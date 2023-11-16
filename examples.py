import math

from DIFI import forward_Sq_d_Re
from DIFI import geod2geoc
from DIFI import jd2000_dt
from DIFI import SwarmL2_F107_Read
from DIFI import SwarmL2_MIO_SHA_Read_v2
import numpy as np
import os

#Get the coefficient files for the DIFI
baseDir = os.getcwd() + "/"
filename_DIFI = baseDir+'DIFI/coefs/difi-coefs.txt'
filename_f107 = baseDir+'DIFI/coefs/f107.DBL'
s = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)
difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(filename_f107)

def main():
    print("Thanks for running the DIFI program")
    print("If you open the python file you can see some quick examples")
    B = getSQfield(21.3166, -157.9996, 2015, 3, 1)
    print("This is the SQ field for lat = 5, lon = 5, on May 1st, 2019")
    print(B)

    year = 2019
    month = np.array([1,2,3,4,5,6,7,8,9,10,11,12])
    day = 15
    B = getSQfield(20, -50, year, month, day)
    print "This is a time series of the field for the fifteenth day of each month in 2019" 
    print '\n'+str(B)+'\n' 
    print "Now attempting with an out of bounds altitude" 


    print("If you have any questions about the software please contact")
    print("Li-Yin Young at: liyin.young@noaa.gov")

def getSQfield(lat, lon, year, month, day, hour=0, minutes=0, h = 0):
    """
    Input:
        Latitude, lat (in WGS-84 coordinates)
        Longtitude, lon
        An array of year, year (Only good between 2014.0 and 2020.0
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
    end_time = 8765.5
    sq_t_arr = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)

    sq_t = sq_t_arr[0]

    if sq_t < start_time:
        raise Exception("Request before 2014.0 reached difi calculation improperly")
    elif sq_t > end_time:
        raise Exception("Request after 2024.0 reached difi calculation improperly")
    while sq_t < 5114.0:
        sq_t +=365

    j = 0
    while difi_t_f107[j] < sq_t:
        j += 1

    frac = sq_t - math.floor(sq_t)
    f107_1 = np.multiply(difi_f107[j] * frac + difi_f107[j - 1] * (1 - frac), np.ones(1))

    B_XYZ = {}
    if h < 20: #LIMIT ALTITUDE REQUESTS TO 20 KM.  Model invalid above!
        [r_gc, theta_gc] = geod2geoc.geod2geoc(np.radians(lat), h)
        r_gc = a + h
        theta_gc = np.degrees(theta_gc)

        [B_1, B_2] = forward_Sq_d_Re.forward_Sq_d_Re(a, theta_gc, lon, sq_t, f107_1, s)

        B_C = B_1+B_2
        B_XYZ['Z'] = -1*B_C[0]
        B_XYZ['Y'] = B_C[2]
        B_XYZ['X'] = -1*B_C[1]
    else:
        raise ValueError("Requested altitude {h} km is higher than 20 km".format(h=repr(h)))
        B_XYZ['Z'] = 0
        B_XYZ['Y'] = 0
        B_XYZ['X'] = 0

    return B_XYZ

def RotateMagneticVector(B_geoc, lat_geoc, latitude):
# Rotate the Magnetic Vectors to Geodetic Coordinates
# Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
# Equation 16, WMM Technical report
#INPUT : lat_geoc: geocentric latitude 
#                latitude:  geodetic latitude
#
#                B_geoc: Dictionary with the following elements
#                        X      North
#                        Y      East
#                        Z      Down
#
#OUTPUT: MagneticResultsGeo Dictionary with the following elements
#                        X      North
#                        Y      East
#                        Z      Down
#
#    /* Difference between the spherical and Geodetic latitudes */
    Psi = np.radians(lat_geoc - latitude)
    B={}

#    /* Rotate spherical field components to the Geodetic system */
    B['Z'] = B_geoc['X'] * np.sin(Psi) + B_geoc['Z'] * np.cos(Psi)
    B['X'] = B_geoc['X'] * np.cos(Psi) - B_geoc['Z'] * np.sin(Psi)
    B['Y'] = B_geoc['Y']
    return B
#MAG_RotateMagneticVector

if __name__ == "__main__":
    main()
