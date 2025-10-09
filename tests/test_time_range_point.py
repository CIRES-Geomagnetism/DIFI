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
def test_single_time_point_warnings():
    with pytest.warns(UserWarning) as record:
        #early warning: input date is just before the minimum recomended time
        B = getSQfield(0, 0, 2000, 
                        1, 1, h = 0,f107_1 = 100, 
                            model_name = 'xdifi2')

        #late warning: input date is just after the maximum recomended time
        B = getSQfield(0, 0, 2024, 
                        1, 2, h = 0,f107_1 = 100, 
                            model_name = 'xdifi2')
        #early warning: input date is just before the minimum recomended time
        B = getSQfield(0, 0, 2000, 
                        1, 1, h = 0,f107_1 = 100, 
                            model_name = 'difi8')
        #late warning: input date is just after the maximum recomended time
        B = getSQfield(0, 0, 2024, 
                        1, 2, h = 0,f107_1 = 100, 
                            model_name = 'difi8')
    assert len(record) == 4

    # Extract the actual warning messages
    messages = [str(w.message) for w in record]

    expected_messages = [
        "Dataset contains date before 2001.0, outside xDIFI2's recommended range",
        "Dataset contains date after 2024.0, outside xDIFI2's recommended range",
        "Dataset contains date before 2014.0, outside DIFI8's recommended range",
        "Dataset contains date after 2024.0, outside DIFI8's recommended range",
    ]

    assert messages == expected_messages