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
    
    return data
def test_single_time_points():
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

if __name__ == '__main__':
    """Should print 4 warnings, for out of time range"""
    test_single_time_points()    