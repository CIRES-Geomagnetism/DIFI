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
def vectorized_out_difi_new_return_arr(lats, lons, hours, year, month, day):
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
def verify_outputs(B):
    f=open("20250310_outputs_difi_new.csv")
    lines = f.readlines()
    all_match = True
    for i in range(0, len(lines)):
        ans = lines[i].split(",")
        for n in range(0,len(ans)):
            ans[n] = float(ans[n])

        if(not np.all(np.isclose(ans,B))):
            print("mismatch at index",i,ans, B)
            all_match = False
    if(all_match):
        print("All matched! :)")


def difi_vector_test_cases(which_case):
    #This test case attests that my additions of vector
    #Methods are consistent with the output from WMM previous to my changes
    #The cases are that alt, lat, lon (rtp) have
    f = open("20250319_outputs_inputs_difi_pypidifi7.csv")
    lats, lons, year, month, day, hours = [],[],[],[],[],[]
    for line in f.readlines():
        # print(line)
        split_line = line.split(',')
        # print(split_line)
        lats.append(float(split_line[0]))
        lons.append(float(split_line[1]))
        year.append(float(split_line[2]))
        month.append(float(split_line[3]))
        day.append(float(split_line[4]))
        hours.append(float(split_line[5]))
    lat = np.array(lats)
    lon = np.array(lons)
    year = np.array(year)
    hours = np.array(hours)
    month = np.array(month)
    day = np.array(day)
    f.close()
    #base case full vector call
    f=open("tests/20250226_outputs_difi_pypidifi7.csv")
    lines = f.readlines()
    ans = np.zeros((len(lines), 3))
    for i in range(0, len(lines)):
        splitted = lines[i].split(",")
        for j in range(0,len(splitted)):
            ans[i,j] = float(splitted[j])
    f.close()
    test_values = ans
    print("beginining test case", which_case)
    # print("tsest value shape", np.shape(test_values))
    if(which_case == 0):#1a) All scalar
        for i in range(0,len(lat)):
            vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i], year[i], month[i], day[i])
        
            assert(np.isclose(vecB['X'], test_values[i,0]))
            assert(np.isclose(vecB['Y'], test_values[i,1]))
            assert(np.isclose(vecB['Z'], test_values[i,2]))
    elif(which_case == 1):#1b) All vector
        vecB = vectorized_out_difi_new_return_arr(lat, lon, hours, year, month, day)
        assert(np.all(np.isclose(vecB['X'], test_values[:,0])))
        assert(np.all(np.isclose(vecB['Y'], test_values[:,1])))
        assert(np.all(np.isclose(vecB['Z'], test_values[:,2])))
        # assert(np.isclose(vecB, test_values))
    elif(which_case == 2):#1c) 1 vector 2 scalar pos
        # print("should produce 4 warnings due to broadcasted locations being in blackout zones")
        print("test case ", which_case)
        for i in range(0,len(lat)):
            if(not i%100):
                print("beginning datapoint ", i)
            vecB = vectorized_out_difi_new_return_arr(lat[i], lon, hours[i], year[i], month[i], day[i])
        
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
            vecB = vectorized_out_difi_new_return_arr(lat, lon[i], hours[i], year[i], month[i], day[i])
       
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
        
    elif(which_case == 3):#1d) vectors have different length
        # 1di) vectors have the same length
        # print("vestigial case ", which_case)
        vecB = vectorized_out_difi_new_return_arr(lat[i:1], lon[i], hours[i], year[i], month[i], day[i])
        vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i:1], hours[i], year[i], month[i], day[i])
        vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i:1], year[i], month[i], day[i])
        vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i], year[i:1], month[i], day[i])
        vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i], year[i], month[i:1], day[i])
        vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i], year[i], month[i], day[i:1])
    elif(which_case == 4):#2a) All scalar time
        print("test case ", which_case)

        
        
    elif(which_case == 5):#2b) All vector time
        print("vestigial ", which_case)

        # model.setup_env(lat, lon, alt)
        # year,month,day  = get_ymd(time, lat)
        # model.setup_time(year, month, day)
        # vec_ans = model.get_all()
        # for i in range(0,len(lat)):
        
        #     assert np.isclose(test_X[i], vec_ans['x'][i], rtol=0, atol=0.1)
        #     assert np.isclose(test_ddec[i], vec_ans['ddec'][i], rtol=0, atol=0.01)
        #     assert np.isclose(test_dec[i] ,vec_ans['dec'][i], rtol=0, atol=0.01)
        #     assert np.isclose(test_H[i] , vec_ans['h'][i], rtol=0, atol=0.1)
        #     assert np.isclose(test_ydot[i] , vec_ans['dy'][i], rtol=0, atol=0.1)
    elif(which_case == 6):#2c) 1 vector 2 scalar time
        print("vestigial ", which_case)
        for i in range(0,len(lat)):
            if(not i%100):
                print("beginning datapoint ", i)
            #Begin hour
            vecB = vectorized_out_difi_new_return_arr(lat, lon, hours[i:1], year[i], month[i], day[i])
        
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
            vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i:1], year[i], month[i], day[i])
       
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
            #Begin year
            vecB = vectorized_out_difi_new_return_arr(lat, lon, hours[i], year[i:1], month[i], day[i])
        
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
            vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i], year[i:1], month[i], day[i])
       
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
            #Begin month
            vecB = vectorized_out_difi_new_return_arr(lat, lon, hours[i], year[i], month[i:1], day[i])
        
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
            vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i], year[i], month[i:1], day[i])
       
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
            #Begin day
            vecB = vectorized_out_difi_new_return_arr(lat, lon, hours[i], year[i], month[i], day[i:1])
        
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))
            vecB = vectorized_out_difi_new_return_arr(lat[i], lon[i], hours[i], year[i], month[i], day[i:1])
       
            assert(np.isclose(vecB['X'][i], test_values[i,0]))
            assert(np.isclose(vecB['Y'][i], test_values[i,1]))
            assert(np.isclose(vecB['Z'][i], test_values[i,2]))

    elif(which_case == 7):#2d) vectors have different length
        print("test case ", which_case)

        # 1di) vectors have the same length
        # year, month, day = get_ymd(time, lat)
        # for i in range(0,len(lat)):

        
        #     model.setup_env(lat, lon, alt)
        #     model.setup_time(year, month, day[i])
        #     vec_ans = model.get_all()
        #     assert np.isclose(test_ddec[i], vec_ans['ddec'][i], rtol=0, atol=0.01)
        #     assert np.isclose(test_X[i], vec_ans['x'][i], rtol=0, atol=.1)
        #     assert np.isclose(test_dec[i] ,vec_ans['dec'][i], rtol=0, atol=0.01)
        #     assert np.isclose(test_H[i] , vec_ans['h'][i], rtol=0, atol=0.1)
        #     assert np.isclose(test_ydot[i] , vec_ans['dy'][i], rtol=0, atol=0.1)
        # for i in range(0,len(lat)):

        
        #     model.setup_env(lat, lon, alt)
        #     model.setup_time(year, month[i], day)
        #     vec_ans = model.get_all()
        #     assert np.isclose(test_ddec[i], vec_ans['ddec'][i], rtol=0, atol=0.01)
        #     assert np.isclose(test_X[i], vec_ans['x'][i], rtol=0, atol=.1)
        #     assert np.isclose(test_dec[i] ,vec_ans['dec'][i], rtol=0, atol=0.01)
        #     assert np.isclose(test_H[i] , vec_ans['h'][i], rtol=0, atol=0.1)
        #     assert np.isclose(test_ydot[i] , vec_ans['dy'][i], rtol=0, atol=0.1)
        # for i in range(0,len(lat)):

        
        #     model.setup_env(lat, lon, alt)
        #     model.setup_time(year[i], month, day)
        #     vec_ans = model.get_all()

        #     assert np.isclose(test_ddec[i], vec_ans['ddec'][i], rtol=0, atol=0.01)
        #     assert np.isclose(test_X[i], vec_ans['x'][i], rtol=0, atol=.1)
        #     assert np.isclose(test_dec[i] ,vec_ans['dec'][i], rtol=0, atol=0.01)
        #     assert np.isclose(test_H[i] , vec_ans['h'][i], rtol=0, atol=0.1)
        #     assert np.isclose(test_ydot[i] , vec_ans['dy'][i], rtol=0, atol=0.1)
        # 
    
        
    elif(which_case == 8):#3a) vector dyear scalar pos
        print("test case ", which_case)
       
    elif(which_case == 9):#3b) scalar dyear vector pos
        print("test case ", which_case)
        
    elif(which_case == 10):#3c) vector dyear vector position
        print("test case ", which_case)
        
    else:
         print(f'you havent created testcase {which_case} yet')
    print(f'passed test case {which_case}')
    return

def main():
    for tests in range(3,10):
        difi_vector_test_cases(tests)

        # lats = np.linspace(-90.0, 90.0, Num_elements)
        # lons = np.linspace(-180.0, 360.0, Num_elements)
        # year = 2024
        # month = 6
        # day = 10
        # hours = np.linspace(1,22,Num_elements)
        # time1 = time.time()
        # outfile = "20250310_outputs_difi_vectorized.csv"

        # vectorized_out_difi_new(outfile, lats, lons, hours, year, month, day)
        # print("time for vectorized input", time.time() - time1)
        # outfile = "20250319_outputs_difi_new.csv"
        
        # time1 = time.time()
        # out_difi_new_test(outfile, lats, lons, hours, year, month, day)
        # print("time for scalar input", time.time() - time1)

        # compare_file_outputs()
        # print("ALL TEST CASES PASSED")

if __name__=="__main__":
    main()



