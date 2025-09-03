import math

from DIFI import SwarmL2_F107_Read
from DIFI import SwarmL2_MIO_SHA_Read_v2
from DIFI import getSQfield
from DIFI import getSQfield as new_sqfield

import numpy as np
import os

#Get the coefficient files for the DIFI
#baseDir = os.getcwd() + "/"
#filename_DIFI = baseDir+'DIFI/coefs/difi-coefs.txt'
#filename_f107 = baseDir+'DIFI/coefs/f107.DBL'
#s = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)
#difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(filename_f107)

def out_difi_old(outfile, lats, lons, hours, years, month, day):


    count = 0

    with open(outfile, "w") as file:
        for lat in lats:
            for lon in lons:
                for hour in hours:
                    for year in years:
                        B = getSQfield(lat, lon, year, month, day, hour=hour, h = 110)



                        file.write(str(B['X'][0]) + "," + str(B['Y'][0]) + "," + str(B['Z'][0]) + "\n")
                        if count % 100 == 0:
                            print("Run " + str(count) + "sample")

                        count += 1

def out_difi_new(outfile, lats, lons, hours, year, month, day):

    count = 0
    with open(outfile, "w") as file:
    
        B = new_sqfield(lats, lons, year, month, day, hour=hours)
        for i in range(0, len(B['X'])):
            file.write(str(B['X'][i]) + "," + str(B['Y'][i]) + "," + str(B['Z'][i]) + "\n")
        # for lat in lats:
        #     for lon in lons:
        #         for hour in hours:
        #             B = new_sqfield(lat, lon, year, month, day, hour=hour)

        #             file.write(str(B['X'][0]) + "," + str(B['Y'][0]) + "," + str(B['Z'][0]) + "\n")
        #             if count % 100 == 0:
        #                 print("Run " + str(count) + "sample")

        #             count += 1

def main():

    N_datapoints = 1000
    lats = np.linspace(-90,90, N_datapoints)
    lons = np.linspace(-180, 180, N_datapoints)
    year = np.ones(N_datapoints)*2020
    month = np.ones(N_datapoints)*6
    days = [1,2,3,4]
    day_gatherer= []
    hour_gatherer = []
    for day in days:
        day_fill = np.ones(N_datapoints//len(days))*day
        day_gatherer.extend(day_fill)
        hour_gatherer.extend(np.linspace(0,23.99, len(day_fill)))

    hours = np.array(hour_gatherer)
    day = np.array(day_gatherer)

    outfile = "20250220_outputs_difi_new.csv"

    out_difi_new(outfile, lats, lons, hours, year, month, day)
    # #out_difi_old_call
    # outfile = "20250220_outputs_difi_old.csv"
    # lats = np.arange(-90.0, 90.0, 9.37)
    # lons = np.arange(-180.0, 360.0, 7.41)
    # year = np.linspace(2000,2024, 60)
    # month = 6
    # day = 10
    # hours = [3, 5, 7]
    # out_difi_old(outfile, lats, lons, hours, year, month, day)

if __name__=="__main__":
    main()



