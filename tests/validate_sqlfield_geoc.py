import numpy as np
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

        for i in range(len(data)):
            data[i]["theta"] = 90 - data[i]["theta"]  # Convert theta to cotangent

    return data


data = parse_geomagnetic_data("tests/test_values_DIFI8_v1_20250528.txt")
print(data[0])
max_diff = 0.0

for i in range(0, len(data)):

    B = getSQfield(data[i]['theta'],  data[i]['lon'], data[i]['year'], data[i]['month'], data[i]['day'], hour=data[i]['hour'], minutes = data[i]['min'], r = data[i]['r'],f107_1 = data[i]['F10.7'], geoc = True, model_name='difi8')
    if B['X'] != data[i]['X']:
         max_diff = max(max_diff, np.abs(B['X'] - data[i]['X']))
    if B['Y'] != data[i]['Y']:
        max_diff = max(max_diff, np.abs(B['Y'] - data[i]['Y']))
    #     print("Y mismatch", B['Y'] , data[i]['Y'], B['Y'] - data[i]['Y'])
    if B['Z'] != data[i]['Z']:
        max_diff = max(max_diff, np.abs(B['Z'] - data[i]['Z']))
    #     print("Z mismatch", B['Z'] , data[i]['Z'], B['Z'] - data[i]['Z'])

    # print(f"lon: {data[i]['lon']}, phi: {data[i]['phi']} lat: {data[i]['lat']}, theta: {data[i]['theta']} lat_np_degree: {np.rad2deg(90 - data[i]['theta'])}")
    if np.abs(B['X'] - data[i]['X']) > 1e-3:
        print("there is a mismatch in B['X'] - data[i]['X'] of size", B['X'] - data[i]['X'])
    if np.abs(B['Y'] - data[i]['Y']) > 1e-3:
        print("there is a mismatch in B['Y'] - data[i]['Y'] of size", B['Y'] - data[i]['Y'])
    if np.abs(B['Z'] - data[i]['Z']) > 1e-3:
        print("there is a mismatch in B['Z'] - data[i]['Z'] of size", B['Z'] - data[i]['Z'])

print(f"Max difference: {max_diff}")