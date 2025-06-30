import numpy as np
# from DIFI import get_SQ_field_TESTING_ONLY as getSQfield
from DIFI.getSQfield import getSQfield

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

data = parse_geomagnetic_data("tests/test_values_xDIFI2_v1_20250528.txt")
print(data[0])
for i in range (0, len(data)):
    B = getSQfield(data[i]['lat'], data[i]['lon'], data[i]['year'], data[i]['month'], data[i]['day'], hour=data[i]['hour'], minutes = data[i]['min'], h = data[i]['r']-6371.2,f107_1 = data[i]['F10.7'], model_name='xdifi2')
    # if B['X'] != data[i]['X']:
    #     print("x mismatch", B['X'] , data[i]['X'], B['X'] - data[i]['X'])
    # if B['Y'] != data[i]['Y']:
    #     print("Y mismatch", B['Y'] , data[i]['Y'], B['Y'] - data[i]['Y'])
    # if B['Z'] != data[i]['Z']:
    #     print("Z mismatch", B['Z'] , data[i]['Z'], B['Z'] - data[i]['Z'])

    # if np.abs(B['X'] - data[i]['X']) > 1e-3:
    #     print("there is a mismatch in B['X'] - data[i]['X'] of size", B['X'] - data[i]['X'], "answer = ", B['X'] ,"calculated = ",  data[i]['X'],"lat",  data[i]['lat'], "lon", data[i]['lon'], 'r',data[i]['r']-6371.2, 'f107', data[i]['F10.7'])
    # if np.abs(B['Y'] - data[i]['Y']) > 1e-3:
    #     print("there is a mismatch in B['Y'] - data[i]['Y'] of size", B['Y'] - data[i]['Y'], "answer = ", B['Y'] ,"calculated = ",  data[i]['Y'],"lat",  data[i]['lat'], "lon", data[i]['lon'], 'r',data[i]['r']-6371.2, 'f107', data[i]['F10.7'])
    # if np.abs(B['Z'] - data[i]['Z']) > 1e-3:
    #     print("there is a mismatch in B['Z'] - data[i]['Z'] of size", B['Z'] - data[i]['Z'], "answer = ", B['Z'] ,"calculated = ",  data[i]['Z'],"lat",  data[i]['lat'], "lon", data[i]['lon'], 'r',data[i]['r']-6371.2, 'f107', data[i]['F10.7'])

    if np.abs(B['X'] - data[i]['X']) > 1e-3:
        print("there is a mismatch in B['X'] - data[i]['X'] of size", B['X'] - data[i]['X'])
        raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the X component")
    if np.abs(B['Y'] - data[i]['Y']) > 1e-3:
        print("there is a mismatch in B['Y'] - data[i]['Y'] of size", B['Y'] - data[i]['Y'])
        raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
    if np.abs(B['Z'] - data[i]['Z']) > 1e-3:
        print("there is a mismatch in B['Z'] - data[i]['Z'] of size", B['Z'] - data[i]['Z'])
        raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
print("Values match to 3 decimal places")
    # print("All values match to 3 decimal places... which isn't perfect considering that's many of the values magnitude")