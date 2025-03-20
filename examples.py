import os

# Get the coefficient files for the DIFI
baseDir = os.path.dirname(__file__) + "/"
import sys
sys.path.insert(1, baseDir + '/DIFI/')

from DIFI import forward_Sq_d_Re
from DIFI import jd2000_dt
from DIFI import SwarmL2_F107_Read
from DIFI import SwarmL2_MIO_SHA_Read_v2
from DIFI import getSQfield as new_sqfield

from output_results import vectorized_out_difi_new_write_file, vectorized_out_difi_new_return_dict

import numpy as np
filename_DIFI = baseDir + 'DIFI/coefs/difi-coefs.txt'
filename_f107 = baseDir + 'DIFI/coefs/f107.DBL'
s = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)
difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(filename_f107)


def main():
    print("Thanks for running the DIFI program")
    print("If you open the python file you can see some quick examples")
    B = new_sqfield(21.3166, -157.9996, 2023, 11, 1)
    print(
        "This is the SQ field for lat = 21.3166, lon = -157.9996, "
        "on November 1th, 2023"
    )
    print(B)

    year = 2023
    month = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    day = 15
    B = new_sqfield(20, -50, year, month, day)
    print()
    print(
        "This is a time series of the field for the 5th day "
        "of each month in 2023"
    )
    print(
        "\n" + str(B) + "\n"
    )

    Num_elements = 1000
    lats = np.linspace(-90.0, 90.0, Num_elements)
    lons = np.linspace(-180.0, 360.0, Num_elements)
    year = 2024
    month = 6
    day = 10
    hours = np.linspace(1,22,Num_elements)
    height = np.linspace(1,22,Num_elements)
    outfile = "example_out_difi_vectorized.csv"

    height = np.linspace(1,19,Num_elements)
    height = 11000
    vectorized_out_difi_new_write_file(outfile, height, lats, lons, hours, year, month, day, write_inputs= False)
    B = vectorized_out_difi_new_return_dict(height,lats, lons, hours, year, month, day)
    print()
    print(
        "\n" , type(B) ,np.shape(B['X']),np.shape(B['Y']),np.shape(B['Z']), "\n"
    )

    print("If you have any questions about the software please contact")
    print("Li-Yin Young at: liyin.young@noaa.gov")

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
    # MAG_RotateMagneticVector


if __name__ == "__main__":
    main()
