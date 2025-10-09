import numpy as np
from DIFI.getSQfield import getSQfield
import DIFI.getSQfield_par as par_SQ
import time

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
if __name__ == '__main__':
    #Create dataset or smth
    N_data = int(1e6)
    lat = np.linspace(1,89,N_data)
    lon = np.linspace(1,89,N_data)
    r = np.linspace(1,1,N_data)
    f107 = np.linspace(100,110,N_data)
    years = np.linspace(2000, 2025, N_data)
    months = np.linspace(1,12, N_data)
    days = np.linspace(1,31, N_data)
    hours = np.linspace(0,23, N_data)
    minutes = np.linspace(0,59, N_data)

    #Create test values with non_par implementation
    time1 = time.time()
    # truth = getSQfield(lat, lon, years, months, days, hour=hours, minutes = minutes, h = r,f107_1 = f107, model_name = 'DIFI8')
    truth = 0
    non_partime = time.time() - time1

    time1 = time.time()
    B = par_SQ.getSQfield(lat, lon, years, months, days, hour=hours, minutes = minutes, h = r,f107_1 = f107, model_name = 'DIFI8')
    partime = time.time() - time1

    
    # print("finished running, now onto comparing")
    max_diff = 0.0

    print("printing one pair of values to avoid any potential null result:\n", B['X'][0] , truth['X'][0])
    if np.max(np.abs(B['X'] - truth['X'])) > 1e-3:
        print("there is a mismatch in B['X'] - truth['X'] of size", B['X'] - truth['X'])
        raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the X component")
    if np.max(np.abs(B['Y'] - truth['Y'])) > 1e-3:
        print("there is a mismatch in B['Y'] - truth['Y'] of size", B['Y'] - truth['Y'])
        raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
    if np.max(np.abs(B['Z'] - truth['Z'])) > 1e-3:
        print("there is a mismatch in B['Z'] - truth['Z'] of size", B['Z'] - truth['Z'])
        raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
    print(f"Max difference: {max_diff}")
    print(f"----\nRuntimes:\nnon_par: {non_partime}\npar: {partime}")