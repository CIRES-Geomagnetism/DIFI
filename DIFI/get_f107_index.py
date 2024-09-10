import numpy as np
from DIFI import difi_t_f107, difi_f107

def get_f107_index(sq_t: list, start_time: float, end_time: float) -> list:

    frac_arr = sq_t - np.floor(sq_t)
    f107_1 = np.array([])
    for i in range(np.size(sq_t)):
        if sq_t[i] < start_time:
            raise Exception(
                "Request before 2014.0 reached difi calculation improper"
            )
        elif sq_t[i] > end_time:
            raise Exception(
                "Request after 2025.0 reached difi calculation improperly"
            )
        while sq_t[i] < 5114.0:
            sq_t[i] += 365
        j = 0
        while difi_t_f107[j] < sq_t[i]:
            j += 1
        f107_1 = np.append(
            f107_1,
            difi_f107[j] * frac_arr[i] + difi_f107[j - 1] * (1 - frac_arr[i])
        )

    return f107_1