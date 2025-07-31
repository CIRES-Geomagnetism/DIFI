import os
import numpy as np
from functools import lru_cache
from DIFI import getSQfield, SwarmL2_F107_Read, SwarmL2_MIO_SHA_Read_v2
from typing import Optional, Union


@lru_cache(maxsize=1)
def load_coefs() -> tuple[list, list]:
    baseDir = os.path.dirname(__file__)
    filename_f107 = os.path.join(baseDir, "coefs", "f107.DBL")
    difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(filename_f107)

    return difi_t_f107, difi_f107
1
@lru_cache(maxsize=1)
def load_swarm_DIFI7() -> dict:
    """
    This function is called in this file to import 
    DIFI7 coefficients
    """
    baseDir = os.path.dirname(__file__)
    filename_DIFI = os.path.join(baseDir, "coefs", "difi-coefs.txt")
    swarm_data = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)

    return swarm_data

@lru_cache(maxsize=1)
def load_swarm_DIFI8() -> dict:
    baseDir = os.path.dirname(__file__)
    filename_DIFI = os.path.join(baseDir, "coefs", "DIFI8.DBL")
    swarm_data = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)

    return swarm_data

@lru_cache(maxsize=1)
def load_swarm_xDIFI() -> dict:
    baseDir = os.path.dirname(__file__)
    filename_DIFI = os.path.join(baseDir, "coefs", "xDIFI2.DBL")
    swarm_data = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)

    return swarm_data


# difi_t_f107, difi_f107 = load_coefs()
# swarm_data = load_swarm_DIFI7()

def get_f107_index(sq_t: list, start_time: float, end_time: float, difi_f107: Union[float,list], difi_t_f107: Union[float,list]) ->np.ndarray:

    frac_arr = (sq_t + .5) - np.floor(sq_t + .5)
    # frac_arr = (sq_t + 0) - np.floor(sq_t + 0)
    f107_1 = np.array([])

    for i in range(np.size(sq_t)):
        if sq_t[i] < start_time:
            raise Exception(
                "This package does not contain f10.7 data before 01/01/2000. Input time data contains a date corresponding to f10.7 data not contained in this package."
            )
        elif sq_t[i] > end_time:
            raise Exception(
                "This package does not contain f10.7 data after noon 12/31/2025. Input time data contains a date corresponding to f10.7 data not contained in this package."
            )
        while sq_t[i] < 0:
            sq_t[i] += 365
        j = 0
        while difi_t_f107[j] <= sq_t[i]:
            j += 1
        if sq_t[i] < .5:
            f107_1 = np.append(
            f107_1,
            difi_f107[j]
            )
        else:
            f107_1 = np.append(
                f107_1,
                difi_f107[j] * frac_arr[i] + difi_f107[j - 1] * (1 - frac_arr[i])
            )
            # f107_1 = np.append(
            #     f107_1,
            #     difi_f107[j] * (1-frac_arr[i]) + difi_f107[j - 1] * (frac_arr[i])
            # )

    return f107_1

