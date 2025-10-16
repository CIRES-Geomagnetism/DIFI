import numpy as np
from DIFI.getSQfield import getSQfield
import pytest


def parse_geomagnetic_data(filepath):
    # Define column labels
    keys = [
        "year", "month", "day", "hour", "min", "sec",
        "r", "theta", "phi", "h", "lat", "lon", "F10.7",
        "B_r", "B_theta", "B_phi", "X", "Y", "Z"
    ]

    data = []
    with open(filepath, 'r') as f:
        for line in f:
            # Split the line into values and convert to floats
            values = list(map(float, line.strip().split()))
            # Create a dictionary for each row
            entry = dict(zip(keys, values))
            data.append(entry)
    
    return data
def test_above_and_below_SQ():
    models = ['xDIFI2', 'DIFI8']
    for model_name in models:
        data = parse_geomagnetic_data(f"tests/test_values_{model_name}_v1_20250928.txt")
        lat, lon, year, month, day, hour, minute, h, f107 = [],[],[],[],[],[],[],[],[]
        for i in range(0,len(data)):
            lat.append(data[i]['lat'])
            lon.append(data[i]['lon'])
            year.append(data[i]['year'])
            month.append(data[i]['month'])
            day.append(data[i]['day'])
            hour.append(data[i]['hour'])
            minute.append(data[i]['min'])
            h.append(data[i]['h'])
            f107.append(data[i]['F10.7'])
        #================================================================================================
        #Create test values
        #================================================================================================
        #Create values below SQ field
        h = np.linspace(80, 105, len(h))

        B = getSQfield(lat, lon, year, 
                        month, day, hour=hour,
                            minutes = minute, h = h,f107_1 = f107, 
                            model_name = model_name)
        B = np.array([B['X'], B['Y'], B['Z']])
        np.save(f"tests/{model_name}_below_SQ_test_values", B)
        h = np.linspace(115, 140, len(h))
        B = getSQfield(lat, lon, year, 
                        month, day, hour=hour,
                            minutes = minute, h = h,f107_1 = f107, 
                            model_name = model_name)
        B = np.array([B['X'], B['Y'], B['Z']])
        np.save(f"tests/{model_name}_after_SQ_test_values", B)
        
        max_diff = 0.0
        
        temp_lat = []
        temp_lon = []
        temp_year = []
        temp_month = []
        temp_day = []
        temp_hour = []
        temp_minute = []
        temp_h = []
        temp_f107 = []
        temp_model_name = []

        temp_lat.extend(lat)
        temp_lon.extend(lon)
        temp_year.extend(year)
        temp_month.extend(month)
        temp_day.extend(day)
        temp_hour.extend(hour)
        temp_minute.extend(minute)
        temp_h.extend(np.linspace(80, 105, len(h)))
        temp_f107.extend(f107)

        temp_lat.extend(lat)
        temp_lon.extend(lon)
        temp_year.extend(year)
        temp_month.extend(month)
        temp_day.extend(day)
        temp_hour.extend(hour)
        temp_minute.extend(minute)
        temp_h.extend(np.linspace(115, 140, len(h)))
        temp_f107.extend(f107)
        lat = temp_lat
        lon = temp_lon
        year = temp_year
        month = temp_month
        day = temp_day
        hour = temp_hour
        minute = temp_minute
        h = temp_h
        f107 = temp_f107


        B = getSQfield(lat, lon, year, 
                        month, day, hour=hour,
                            minutes = minute, h = h,f107_1 = f107, 
                            model_name = model_name)
        below_ans = np.load(f"tests/{model_name}_below_SQ_test_values.npy", allow_pickle= False)
        above_ans = np.load(f"tests/{model_name}_after_SQ_test_values.npy", allow_pickle= False)
        ansX = []
        ansY = []
        ansZ = []
        ansX.extend(below_ans[0])
        ansY.extend(below_ans[1])
        ansZ.extend(below_ans[2])
        
        ansX.extend(above_ans[0])
        ansY.extend(above_ans[1])
        ansZ.extend(above_ans[2])
        for i in range(0, len(ansX)):
            max_diff = max(max_diff, np.abs(ansX[i] - B["X"][i]))
            max_diff = max(max_diff, np.abs(ansY[i] - B["Y"][i]))
            max_diff = max(max_diff, np.abs(ansZ[i] - B["Z"][i]))
        
            # if not i // 10:
            #     print(np.abs(ansX[i] - B["X"][i]), ansX[i] ,  B["X"][i])
            #     print(np.abs(ansY[i] - B["Y"][i]), ansY[i] ,  B["Y"][i])
            #     print(np.abs(ansZ[i] - B["Z"][i]), ansZ[i] ,  B["Z"][i])
        B = getSQfield(lat[0], lon[0], year[0], 
                        month[0], day[0], hour=hour[0],
                            minutes = minute[0], h = h[0],f107_1 = f107[0], 
                            model_name = model_name)
        below_ans = np.load(f"tests/{model_name}_below_SQ_test_values.npy", allow_pickle= False)
        above_ans = np.load(f"tests/{model_name}_after_SQ_test_values.npy", allow_pickle= False)
        ansX = []
        ansY = []
        ansZ = []
        ansX.extend(below_ans[0])
        ansY.extend(below_ans[1])
        ansZ.extend(below_ans[2])
        
        ansX.extend(above_ans[0])
        ansY.extend(above_ans[1])
        ansZ.extend(above_ans[2])
        #Create inputs which have span above and below the SQ field
        
        # for i in range(0, len(data)):
        #     # B = getSQfield(data[i]['lat'], data[i]['lon'], data[i]['year'], data[i]['month'], data[i]['day'], hour=data[i]['hour'], minutes = data[i]['min'], h= data[i]['h'],f107_1 = data[i]['F10.7'], model_name='xdifi2')

        #     if B['X'][i] != data[i]['X']:
        #         max_diff = max(max_diff, np.abs(B['X'][i] - data[i]['X']))
        #     if B['Y'][i] != data[i]['Y']:
        #         max_diff = max(max_diff, np.abs(B['Y'][i] - data[i]['Y']))
        #     #     print("Y mismatch", B['Y'][i] , data[i]['Y'], B['Y'][i] - data[i]['Y'])
        #     if B['Z'][i] != data[i]['Z']:
        #         max_diff = max(max_diff, np.abs(B['Z'][i] - data[i]['Z']))

        #     if np.abs(B['X'][i] - data[i]['X']) > 1e-3:
        #         print("there is a mismatch in B['X'][i] - data[i]['X'] of size", B['X'][i] - data[i]['X'])
        #         raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the X component")
        #     if np.abs(B['Y'][i] - data[i]['Y']) > 1e-3:
        #         print("there is a mismatch in B['Y'][i] - data[i]['Y'] of size", B['Y'][i] - data[i]['Y'])
        #         raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
        #     if np.abs(B['Z'][i] - data[i]['Z']) > 1e-3:
        #         print("there is a mismatch in B['Z'][i] - data[i]['Z'][i] of size", B['Z'][i] - data[i]['Z'])
        #         raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
        print(f"Max difference of {model_name} when calculating above and below SQ field simultaneously vs individually input: {max_diff}. Should be zero")
        assert(max_diff < 1e-12)
