import math

from DIFI import SwarmL2_F107_Read
from DIFI import SwarmL2_MIO_SHA_Read_v2
from DIFI import getSQfield as new_sqfield


import numpy as np
import os

def vectorized_out_difi_new_return_dict(lats, lons, hours, year, month, day):
    B = new_sqfield(lats, lons, year, month, day, hour=hours)
    return B
def vectorized_out_difi_new_write_file(outfile, lats, lons, hours, year, month, day,write_inputs = False):
    with open(outfile, "w") as file:
        B = new_sqfield(lats, lons, year, month, day, hour=hours)
        for i in range(0, len(B['X'])):
            if(not write_inputs):
                file.write(str(B['X'][i]) + "," + str(B['Y'][i]) + "," + str(B['Z'][i]) + "\n")
            else:
                file.write(str(lats[i]) + "," + str(lons[i]) + "," + str(year[i]) + "," +str(month[i]) + "," +str(day[i]) + "," +str(hours[i]) + "," +str(B['X'][i]) + "," + str(B['Y'][i]) + "," + str(B['Z'][i]) + "\n")
  
def main():

    Num_elements = 1000
    lats = np.linspace(-90.0, 90.0, Num_elements)
    lons = np.linspace(-180.0, 360.0, Num_elements)
    year = 2024
    month = 6
    day = 10
    hours = np.linspace(1,22,Num_elements)
    outfile = "20250320_outputs_difi_vectorized.csv"

    vectorized_out_difi_new_write_file(outfile, lats, lons, hours, year, month, day, write_inputs= False)


if __name__=="__main__":
    main()



