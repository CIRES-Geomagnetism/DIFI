import numpy as np
from DIFI.getSQfield import getSQfield
def parse_sq_input_output(filepath):
    import re
    import numpy as np

    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            key = parts[0]

            # Handle keys like B_1_1(1)
            match = re.match(r'([A-Za-z0-9_]+)\((\d+)\)', key)
            if match:
                base_key, idx = match.groups()
                idx = int(idx) - 1  # convert to 0-based index
                if base_key not in data:
                    data[base_key] = []
                values = list(map(float, parts[1:]))
                # Ensure we can assign at the right index
                while len(data[base_key]) <= idx:
                    data[base_key].append([])
                data[base_key][idx] = values
            else:
                # Regular scalar or vector key
                values = list(map(float, parts[1:]))
                data[key] = values[0] if len(values) == 1 else values

    # Convert matrix-like lists to numpy arrays
    for k in ['B_1_1', 'B_2_1']:
        if k in data:
            data[k] = np.array(data[k])

    return data

inputs = parse_sq_input_output('tests/test_val_above_SQ.txt')
# print("inptuts", inputs)
inputs['theta_1'] = 90 - np.array(inputs['theta_1'])
B = getSQfield(inputs['theta_1'], inputs['phi_1'], year = 2014, month = 3, day = 21, hour = 12, h = inputs['r_1']-6371.2,f107_1 = inputs['f107_1'], model_name = 'difi8')#, geoc = True )
Bx = -B['X']
By = B['Y']
Bz = -B['Z']
True_Bz = inputs['B_1_1'][0] + inputs['B_2_1'][0]
True_Bx = inputs['B_1_1'][1] + inputs['B_2_1'][1]
True_By = inputs['B_1_1'][2] + inputs['B_2_1'][2]
error = []
for i in range(0, len(Bx)):
    # print("Bx diff, mine/mat/diff", Bx[i], True_Bx[i], Bx[i] - True_Bx[i])
    error.append(Bx[i] - True_Bx[i])
    # print("Bz diff, mine/mat/diff", Bz[i], True_Bz[i], Bz[i] - True_Bz[i])
    error.append(Bz[i] - True_Bz[i])
    # print("By diff, mine/mat/diff", By[i], True_By[i], By[i] - True_By[i])
    error.append(By[i] - True_By[i])
error = np.array(error)
print("average error in outputs between matlab test values and python is", np.average(error), "max is ", np.max(error))