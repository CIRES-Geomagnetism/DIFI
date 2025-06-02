import datetime
import numpy as np
from typing import Union


def jd2000_dt(year: Union[float, list], month: Union[float, list], day: Union[float, list], UT: Union[float, list]=0, minutes: Union[float, list]=0) -> np.ndarray:
    # set offset to 2000
    offset = datetime.datetime(2000, 1, 1)
    td_array = np.array([])
    size = (
        np.size(year),
        np.size(month),
        np.size(day),
        np.size(UT),
        np.size(minutes)
    )
    max_size = max(size)
    if np.size(year) == 1:
        year = np.ones(max_size, dtype=int) * year
    if np.size(month) == 1:
        month = np.ones(max_size, dtype=int) * month
    if np.size(day) == 1:
        day = np.ones(max_size, dtype=int) * day
    if np.size(UT) == 1:
        UT = np.ones(max_size, dtype=int) * UT
    if np.size(minutes) == 1:
        minutes = np.ones(max_size, dtype=int) * minutes
    for i in range(np.size(year)):
        td = datetime.datetime(
            int(year[i]),
            int(month[i]),
            int(day[i]),
            int(UT[i]),
            int(minutes[i])
        ) - offset
        td_f = td.days+td.seconds/(3600.0*24)
        td_array = np.append(td_array, td_f)

    return td_array
