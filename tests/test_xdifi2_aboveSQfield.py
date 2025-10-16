import numpy as np
# from DIFI import get_SQ_field_TESTING_ONLY as getSQfield
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
def test_xdifi2_accuracy():
    with pytest.warns(UserWarning) as record:
        data = parse_geomagnetic_data("tests/test_values_xDIFI2_v2_aboveSQfield_20251015.txt")
        max_diff = 0.0

        for i in range(0, len(data)):
            B = getSQfield(data[i]['lat'], data[i]['lon'], data[i]['year'], data[i]['month'], data[i]['day'], hour=data[i]['hour'], minutes = data[i]['min'], h= data[i]['h'],f107_1 = data[i]['F10.7'], model_name='xdifi2')

            if B['X'] != data[i]['X']:
                max_diff = max(max_diff, np.abs(B['X'] - data[i]['X']))
            if B['Y'] != data[i]['Y']:
                max_diff = max(max_diff, np.abs(B['Y'] - data[i]['Y']))
            #     print("Y mismatch", B['Y'] , data[i]['Y'], B['Y'] - data[i]['Y'])
            if B['Z'] != data[i]['Z']:
                max_diff = max(max_diff, np.abs(B['Z'] - data[i]['Z']))

            if np.abs(B['X'] - data[i]['X']) > 1e-3:
                print("there is a mismatch in B['X'] - data[i]['X'] of size", B['X'] - data[i]['X'], B['X'] , data[i]['X'])
                raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the X component")
            if np.abs(B['Y'] - data[i]['Y']) > 1e-3:
                print("there is a mismatch in B['Y'] - data[i]['Y'] of size", B['Y'] - data[i]['Y'], B['Y'] , data[i]['Y'])
                raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
            if np.abs(B['Z'] - data[i]['Z']) > 1e-3:
                print("there is a mismatch in B['Z'] - data[i]['Z'] of size", B['Z'] - data[i]['Z'], B['Z'] , data[i]['Z'])
                raise ValueError("There was a mismatch of magnitude greater than 1e-3 in the Y component")
        # print(f"Max difference: {max_diff}")
        # messages = [str(w.message) for w in record]
        # print("The thing inside messages:",messages)
        # expected_messages = [
        #     "Warning: The altitude is within 30 km of 110 km, where the model switches between internal and external ionospheric current sources. The model is not expected to accurately represent the field within this transition zone.', 'Warning: The altitude is within 30 km of 110 km, where the model switches between internal and external ionospheric current sources. The model is not expected to accurately represent the field within this transition zone.', 'Warning: The altitude is within 30 km of 110 km, where the model switches between internal and external ionospheric current sources. The model is not expected to accurately represent the field within this transition zone."
        # ]

        # assert messages == expected_messages[0]
        assert(max_diff < 1e-3)

# test_xdifi2_accuracy()