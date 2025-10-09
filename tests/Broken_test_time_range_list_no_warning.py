"""
This test is supposed to affirm that warnings are not printed when they shouldn't. 
It seems like there's a problem with asserting that "record" is empty.
This test case gets left behind due to its unimportance and holding up development of CM6 and HDGMXT.
"""

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

def test_warning_noprint():
    with pytest.warns(UserWarning) as record:

        #all input dates are within time range, no warnings should print
        model_name = 'difi8'
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
        N = 40
        #Create time from 2000.0 to 2025.9
        lat, lon, year, month, day, hour, minute, h, f107 = lat[:N], lon[:N], year[:N], month[:N], day[:N], hour[:N], minute[:N], h[:N], f107[:N]
        years = np.linspace(2001, 2023, N)
        months = np.linspace(1,12, N)
        days = np.linspace(1,31, N)
        hours = np.linspace(0,23, N)
        minutes = np.linspace(0,59, N)

        B = getSQfield(lat, lon,years, 
                        months, days, hour= hours, minutes=minutes ,
                        h = 0,f107_1 = 100, 
                        model_name = 'xdifi2')
        
        years = np.linspace(2014, 2023, N)
        B = getSQfield(lat, lon,years, 
                        months, days, hour= hours, minutes=minutes ,
                        h = 0,f107_1 = 100, 
                        model_name = 'difi8')
    print(len(record))
    assert len(record) == 0

with pytest.warns(UserWarning) as record:
    #all input dates are within time range, no warnings should print
    model_name = 'difi8'
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
    N = 40
    #Create time from 2000.0 to 2025.9
    lat, lon, year, month, day, hour, minute, h, f107 = lat[:N], lon[:N], year[:N], month[:N], day[:N], hour[:N], minute[:N], h[:N], f107[:N]
    years = np.linspace(2001, 2023, N)
    months = np.linspace(1,12, N)
    days = np.linspace(1,31, N)
    hours = np.linspace(0,23, N)
    minutes = np.linspace(0,59, N)

    B = getSQfield(lat, lon,years, 
                    months, days, hour= hours, minutes=minutes ,
                    h = 0,f107_1 = 100, 
                    model_name = 'xdifi2')
    
    years = np.linspace(2014, 2023, N)
    B = getSQfield(lat, lon,years, 
                    months, days, hour= hours, minutes=minutes ,
                    h = 0,f107_1 = 100, 
                    model_name = 'difi8')
print("???")
print(len(record))
# assert len(record) == 0