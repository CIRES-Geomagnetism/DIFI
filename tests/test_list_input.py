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
def test_list_input_accuracy():
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


        B = getSQfield(lat, lon, year, 
                        month, day, hour=hour,
                            minutes = minute, h = h,f107_1 = f107, 
                            model_name = model_name)
        max_diff = 0.0
        for i in range(0, len(data)):
            # B = getSQfield(data[i]['lat'], data[i]['lon'], data[i]['year'], data[i]['month'], data[i]['day'], hour=data[i]['hour'], minutes = data[i]['min'], h= data[i]['h'],f107_1 = data[i]['F10.7'], model_name='xdifi2')

            if B['X'][i] != data[i]['X']:
                max_diff = max(max_diff, np.abs(B['X'][i] - data[i]['X']))
            if B['Y'][i] != data[i]['Y']:
                max_diff = max(max_diff, np.abs(B['Y'][i] - data[i]['Y']))
            #     print("Y mismatch", B['Y'][i] , data[i]['Y'], B['Y'][i] - data[i]['Y'])
            if B['Z'][i] != data[i]['Z']:
                max_diff = max(max_diff, np.abs(B['Z'][i] - data[i]['Z']))

            if np.abs(B['X'][i] - data[i]['X']) > 1e-3:
                print("there is a mismatch in B['X'][i] - data[i]['X'] of size", B['X'][i] - data[i]['X'])
                raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the X component")
            if np.abs(B['Y'][i] - data[i]['Y']) > 1e-3:
                print("there is a mismatch in B['Y'][i] - data[i]['Y'] of size", B['Y'][i] - data[i]['Y'])
                raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
            if np.abs(B['Z'][i] - data[i]['Z']) > 1e-3:
                print("there is a mismatch in B['Z'][i] - data[i]['Z'][i] of size", B['Z'][i] - data[i]['Z'])
                raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
        print(f"Max difference of {model_name} list input: {max_diff}")
        assert(max_diff < 1e-3)
