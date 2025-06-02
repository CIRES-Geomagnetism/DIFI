import os

# Get the coefficient files for the DIFI
baseDir = os.path.dirname(__file__) + "/"
import sys
sys.path.insert(1, baseDir + '/DIFI/')

from DIFI import forward_Sq_d_Re
# import geomaglib
# print(dir(geomaglib.util))
# from geomaglib.util import geod_to_geoc_lat as geod2geoc
from DIFI import jd2000_dt
from DIFI import SwarmL2_F107_Read
from DIFI import SwarmL2_MIO_SHA_Read_v2
import numpy as np
filename_DIFI = baseDir + 'DIFI/coefs/difi-coefs.txt'
filename_f107 = baseDir + 'DIFI/coefs/f107.DBL'
s = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)
difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(filename_f107)


def main():
    print("Thanks for running the DIFI program")
    print("If you open the python file you can see some quick examples")
    B = getSQfield(21.3166, -157.9996, 2023, 11, 1)
    print(
        "This is the SQ field for lat = 21.3166, lon = -157.9996, "
        "on November 1th, 2023"
    )
    print(B)

    year = 2023
    month = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    day = 15
    B, f107 = getSQfield(20, -50, year, month, day)
    print()
    print(
        "This is a time series of the field for the 5th day "
        "of each month in 2023"
    )
    print(
        "\n" + str(B) + "\n"
    )


    print("If you have any questions about the software please contact")
    print("Li-Yin Young at: liyin.young@noaa.gov")


# def getSQfield(lat, lon, year, month, day, hour=0, minutes=0, h=0):
#     """
#     Input:
#         Latitude, lat (in WGS-84 coordinates)
#         Longtitude, lon
#         An array of year, year (Only good between 2014.0 and 2025.0
#         An array of months, month
#         An array of days, day

#     Optional Input:
#         An array of hours, hour
#         An array of minutes, minutes
#         Height above WGS84 ellipsoid, h

#     Output:
#         B, the magnetic field due to the SQ in WGS-84 coordinates
#     """
#     a = 6371.2
#     start_time = 5114.0
#     end_time = 9131.5
#     sq_t = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)
#     frac_arr = sq_t - np.floor(sq_t)
#     f107_1 = np.array([])
#     for i in range(np.size(sq_t)):
#         if sq_t[i] < start_time:
#             raise Exception(
#                 "Request before 2014.0 reached difi calculation improper"
#             )
#         elif sq_t[i] > end_time:
#             raise Exception(
#                 "Request after 2025.0 reached difi calculation improperly"
#             )
#         while sq_t[i] < 5114.0:
#             sq_t[i] += 365
#         j = 0
#         while difi_t_f107[j] < sq_t[i]:
#             j += 1
#         f107_1 = np.append(
#             f107_1,
#             difi_f107[j] * frac_arr[i] + difi_f107[j - 1] * (1 - frac_arr[i])
#         )
#     # index_f107 = [
#     #     i for i in range(len(difi_t_f107))
#     #     if (
#     #         difi_t_f107[i] >= math.floor(sq_t)
#     #         and difi_t_f107[i] <= math.ceil(sq_t)
#     #     )
#     # ]
#     # frac = sq_t - math.floor(sq_t)
#     # f107_1 = np.multiply(
#     #     difi_f107[index_f107[0]]
#     #     * (1-frac) + difi_f107[index_f107[-1]]*(frac),
#     #     np.ones(1),
#     # )
#     B_XYZ = {}
#     if h < 20:  # LIMIT ALTITUDE REQUESTS TO 20 KM.  Model invalid above!
#         [r_gc, theta_gc] = geod2geoc.geod2geoc(np.radians(lat), h)
#         r_gc = a + h
#         theta_gc = np.degrees(theta_gc)
#         # print "Difi input", r_gc, theta_gc, RV['lon'], sq_t, f107_1
#         [B_1, B_2] = forward_Sq_d_Re.forward_Sq_d_Re(
#             r_gc,
#             theta_gc,
#             lon,
#             sq_t,
#             f107_1,
#             s,
#         )
#         # print "Difi output", B_1, B_2
#         B_C = B_1 + B_2
#         # B_C = B_1
#         B_XYZ['Z'] = -1 * B_C[0]
#         B_XYZ['Y'] = B_C[2]
#         B_XYZ['X'] = -1 * B_C[1]
#     else:
#         raise ValueError(
#             "Requested altitude {h} km is higher than 20 km".format(h=repr(h))
#         )
#         B_XYZ['Z'] = 0
#         B_XYZ['Y'] = 0
#         B_XYZ['X'] = 0
#     B = RotateMagneticVector(B_XYZ, 90 - theta_gc, lat)
#     # B = B_XYZ
#     return B, f107_1


# def RotateMagneticVector(B_geoc, lat_geoc, latitude):
#     # Rotate the Magnetic Vectors to Geodetic Coordinates
#     # Manoj Nair, June, 2009 Manoj.C.Nair@Noaa.Gov
#     # Equation 16, WMM Technical report
#     # INPUT : lat_geoc: geocentric latitude
#     #                latitude:  geodetic latitude
#     #
#     #                B_geoc: Dictionary with the following elements
#     #                        X      North
#     #                        Y      East
#     #                        Z      Down
#     #
#     # OUTPUT: MagneticResultsGeo Dictionary with the following elements
#     #                        X      North
#     #                        Y      East
#     #                        Z      Down
#     #
#     # Difference between the spherical and Geodetic latitudes
#     Psi = np.radians(lat_geoc - latitude)
#     B = {}

#     # Rotate spherical field components to the Geodetic system
#     B['Z'] = B_geoc['X'] * np.sin(Psi) + B_geoc['Z'] * np.cos(Psi)
#     B['X'] = B_geoc['X'] * np.cos(Psi) - B_geoc['Z'] * np.sin(Psi)
#     B['Y'] = B_geoc['Y']
#     return B
    # MAG_RotateMagneticVector


if __name__ == "__main__":
    main()
