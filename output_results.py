import math

from memory_profiler import profile
from DIFI import SwarmL2_F107_Read
from DIFI import SwarmL2_MIO_SHA_Read_v2
from DIFI import getSQfield as new_sqfield


import numpy as np
import os

#Get the coefficient files for the DIFI
#baseDir = os.getcwd() + "/"
#filename_DIFI = baseDir+'DIFI/coefs/difi-coefs.txt'
#filename_f107 = baseDir+'DIFI/coefs/f107.DBL'
#s = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)
#difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(filename_f107)

def out_difi_old(outfile, lats, lons, hours, year, month, day):


    count = 0

    with open(outfile, "w") as file:
        for lat in lats:
            for lon in lons:
                for hour in hours:
                    B, f107 = new_sqfield(lat, lon, year, month, day, hour=hour)



                    file.write(str(B['X'][0]) + "," + str(B['Y'][0]) + "," + str(B['Z'][0]) + "\n")
                    if count % 100 == 0:
                        print("Run " + str(count) + "sample")

                    count += 1

def out_difi_new(outfile, lats, lons, hours, year, month, day):

    count = 0

    with open(outfile, "w") as file:
        for lat in lats:
            for lon in lons:
                for hour in hours:
                    B = new_sqfield(lat, lon, year, month, day, hour=hour)

                    file.write(str(B['X'][0]) + "," + str(B['Y'][0]) + "," + str(B['Z'][0]) + "\n")
                    if count % 100 == 0:
                        print("Run " + str(count) + "sample")

                    count += 1

def out_difi_new_test(outfile, lats, lons, hours, year, month, day):

    count = 0

    with open(outfile, "w") as file:
        B = np.zeros((len(lats), 3))
        for i in range(0, len(lats)):
            tempB = new_sqfield(lats[i], lons[i], year, month, day, hour=hours[i])
            B[i,0], B[i,1], B[i,2] = tempB['X'][0], tempB['Y'][0] ,tempB['Z'][0]
        return B
            # file.write(str(B['X'][0]) + "," + str(B['Y'][0]) + "," + str(B['Z'][0]) + "\n")
            # if count % 100 == 0:
            #     print("Run " + str(count) + "sample")

            # count += 1
def vectorized_out_difi_new(outfile, lats, lons, hours, year, month, day):
    with open(outfile, "w") as file:
        B = new_sqfield(lats, lons, year, month, day, hour=hours)
        return B
        print(np.shape(B))
        for i in range(0, len(B['X'])):
            file.write(str(B['X'][i]) + "," + str(B['Y'][i]) + "," + str(B['Z'][i]) + "\n")
        
        # if count % 100 == 0:
        #     print("Run " + str(count) + "sample")
    
def compare_file_outputs():
    f=open("20250310_outputs_difi_new.csv")
    f1=open("20250310_outputs_difi_vectorized.csv")
    vec_lines = f1.readlines()
    lines = f.readlines()
    for i in range(0, len(lines)):
        ans = lines[i].split(",")
        pencil = vec_lines[i].split(",")
        for n in range(0,len(ans)):
            ans[n] = float(ans[n])
            pencil[n] = float(pencil[n])
        if(not np.all(np.isclose(ans,pencil))):
            print("idk man...")
            print(ans, pencil)

def main():
    for Num_elements in [int(100000)]:
        # print("beginning datapoints = ", 10**powers, " test")
        # Num_elements = 10**powers
        print("number_of_datapoints is ",Num_elements)
        import time
        lats = np.linspace(-90.0, 90.0, Num_elements)
        lons = np.linspace(-180.0, 360.0, Num_elements)
        year = 2024
        month = 6
        day = 10
        hours = np.linspace(1,22,Num_elements)
        time1 = time.time()
        outfile = "20250310_outputs_difi_vectorized.csv"

        vectorized_out_difi_new(outfile, lats, lons, hours, year, month, day)
        print("time for vectorized input", time.time() - time1)
        outfile = "20250310_outputs_difi_new.csv"

        # time1 = time.time()
        # out_difi_new_test(outfile, lats, lons, hours, year, month, day)
        # print("time for scalar input", time.time() - time1)

        # compare_file_outputs()
        # print("ALL TEST CASES PASSED")

if __name__=="__main__":
    main()



