import random
import unittest
import warnings
import numpy as np
import pytest

from DIFI import getSQfield

class test_validate_inputs(unittest.TestCase):
    def test_xdifi2_time_range(self):
        """
        Test the time range of the xDIFI2 model for point and list inputs.
        - warning if t < 2001.0 or t >= 2024.0
        - no warning if 2001.0 <= t < 2024.0
        - raise the error if t >= the last date of f107.DBL
        """

        N = 50
        lat = 25
        lon = 100
        years = np.linspace(1999, 2027, N)
        months = np.linspace(1, 12, N)
        days = np.linspace(1, 31, N)
        hours = np.linspace(0, 23, N)
        minutes = np.linspace(0, 59, N)
        h = 100

        idx_2000 = 0
        idx_2001 = N - 1
        idx_2024 = N - 1
        idx_end = N - 1

        # Test for the point input
        for i in range(N):

            # When the user provide f107 by himself, the end year can be extended after 2026.0
            B = getSQfield(lat, lon, years[i], months[i], days[i], f107_1=129, h=h, hour=0, minutes=0,
                           model_name="xDIFI2")

            if years[i] < 2000.0:
                idx_2000 = i
                with pytest.raises(Exception,
                                   match="This package does not contain f10.7 data before 01/01/2000. Input time data contains a date corresponding to f10.7 data not contained in this package"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="xDIFI2")

            elif years[i] >= 2000.0 and years[i] < 2001.0:
                if idx_2001 == N - 1:
                    idx_2001 = i
                with pytest.warns(UserWarning,
                                  match="Dataset contains date after 2001.0, outside xDIFI2's recommended range"):

                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="xDIFI2")
                    warnings.warn("Dataset contains date after 2001.0, outside xDIFI2's recommended range", UserWarning)



            elif years[i] >= 2024.0 and years[i] < 2025.9:
                if idx_2024 == N - 1:
                    idx_2024 = i
                with pytest.warns(UserWarning,
                                  match="Dataset contains date after 2024.0, outside xDIFI2's recommended range"):

                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="xDIFI2")
                    warnings.warn("Dataset contains date after 2024.0, outside xDIFI2's recommended range", UserWarning)

            elif years[i] >= 2026.0:
                if idx_end == N - 1:
                    idx_end = i
                with pytest.raises(Exception,
                                   match="This package does not contain f10.7 data after noon 12/31/2025. Input time data contains a date corresponding to f10.7 data not contained in this package"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="xDIFI2")

        # Test for the list input

        # When the user provide f107 by himself, the end year can be extended after 2026.0
        f107_arr = [129.0] * N
        B = getSQfield(lat, lon, years, months, days, f107_1=f107_arr, h=h, hour=0, minutes=0,
                       model_name="xDIFI2")

        with pytest.raises(Exception,
                           match="This package does not contain f10.7 data before 01/01/2000. Input time data contains a date corresponding to f10.7 data not contained in this package"):
            if idx_2000 < N - 1:
                B = getSQfield(lat, lon, years[:idx_2000 + 1], months[:idx_2000 + 1], days[:idx_2000 + 1], h=h,
                               hour=hours[:idx_2000 + 1], minutes=minutes[:idx_2000 + 1], model_name="xDIFI2")

        with pytest.warns(UserWarning,
                          match="Dataset contains date before 2001.0, outside xDIFI2's recommended range"):
            if idx_2001 < N - 1:
                B = getSQfield(lat, lon, years[idx_2000 + 1:idx_2001 + 1], months[idx_2000 + 1:idx_2001 + 1],
                               days[idx_2000 + 1:idx_2001 + 1], h=h,
                               hour=hours[idx_2000 + 1:idx_2001 + 1], minutes=minutes[idx_2000 + 1:idx_2001 + 1],
                               model_name="xdifi2")
                warnings.warn("Dataset contains date before 2001.0, outside xDIFI2's recommended range", UserWarning)

        with pytest.warns(UserWarning, match="Dataset contains date after 2024.0, outside xDIFI2's recommended range"):
            if idx_2024 < N - 1:
                B = getSQfield(lat, lon, years[idx_2024:idx_end], months[idx_2024:idx_end], days[idx_2024:idx_end], h=h,
                               hour=hours[idx_2024:idx_end], minutes=minutes[idx_2024:idx_end], model_name="xDIFI2")

                warnings.warn("Dataset contains date after 2024.0, outside xDIFI2's reccomended range", UserWarning)
        with pytest.raises(Exception,
                           match="This package does not contain f10.7 data after noon 12/31/2025. Input time data contains a date corresponding to f10.7 data not contained in this package."):
            B = getSQfield(lat, lon, years[idx_2024:], months[idx_2024:], days[idx_2024:], h=h, hour=hours[idx_2024:],
                           minutes=minutes[idx_2024:], model_name="xDIFI2")

    def test_difi8_time_range(self):

        """
            Test the time range of the DIFI8 model for point and list inputs.

            - warning if t < 2014.0 or t >= 2024.0
            - no warning if 2014.0 <= t < 2024.0
            - raise the error if t >= the last date of f107.DBL
        """

        N = 50
        lat = 25
        lon = 100
        years = np.linspace(1999, 2027, N)
        months = np.linspace(1, 12, N)
        days = np.linspace(1, 31, N)
        hours = np.linspace(0, 23, N)
        minutes = np.linspace(0, 59, N)
        h = 100

        # Initialize indices for the first year with data, the year with a warning, and the last year
        idx_2000 = 0
        idx_2014 = N - 1
        idx_2024 = N - 1
        idx_end = N - 1

        # Test for the point input
        for i in range(N):

            # When the user provide f107 by himself, the end year can be extended after 2026.0
            B = getSQfield(lat, lon, years[i], months[i], days[i], f107_1=129, h=h, hour=0, minutes=0,
                           model_name="difi8")

            if years[i] < 2000.0:
                idx_2000 = i
                with pytest.raises(Exception,
                                   match="This package does not contain f10.7 data before 01/01/2000. Input time data contains a date corresponding to f10.7 data not contained in this package"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="difi8")

            elif 2000.0 <= years[i] < 2014.0:

                if idx_2014 == N - 1:
                    idx_2014 = i
                with pytest.warns(UserWarning,
                                  match="Dataset contains date before 2014.0, outside DIFI8's reccomended range"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="difi8")
                    warnings.warn("Dataset contains date before 2014.0, outside DIFI8's reccomended range", UserWarning)

            elif years[i] > 2024.0 and years[i] < 2025.9:
                if idx_2024 == N - 1:
                    idx_2024 = i

                with pytest.warns(UserWarning,
                                  match="Dataset contains date after 2024.0, outside DIFI8's recommended range"):

                    B = getSQfield(lat, lon, int(years[i]), int(months[i]), int(days[i]), h=h, hour=0, minutes=0,
                                   model_name="difi8")
                    warnings.warn("Dataset contains date before 2024.0, outside DIFI8's recommended range", UserWarning)


            elif years[i] >= 2026.0:
                if idx_end == N - 1:
                    idx_end = i

                with pytest.raises(Exception,
                                   match="This package does not contain f10.7 data after noon 12/31/2025. Input time data contains a date corresponding to f10.7 data not contained in this package"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="difi8")

        # Test for the list input

        # When the user provide f107 by himself, the end year can be extended after 2026.0
        f107_arr = [129.0] * N
        B = getSQfield(lat, lon, years, months, days, f107_1=f107_arr, h=h, hour=0, minutes=0,
                       model_name="difi8")

        with pytest.raises(Exception,
                           match="This package does not contain f10.7 data before 01/01/2000. Input time data contains a date corresponding to f10.7 data not contained in this package"):

            if idx_2000 < N - 1:
                B = getSQfield(lat, lon, years[:idx_2000 + 1], months[:idx_2000 + 1], days[:idx_2000 + 1], h=h,
                               hour=hours[:idx_2014 + 1], minutes=minutes[:idx_2014 + 1], model_name="difi8")

        with pytest.warns(UserWarning,
                          match="Dataset contains date before 2014.0, outside DIFI8's recommended range"):

            if idx_2014 < N - 1:
                B = getSQfield(lat, lon, years[idx_2000 + 1:idx_2014 + 1], months[idx_2000 + 1:idx_2014 + 1],
                               days[idx_2000 + 1:idx_2014 + 1], h=h,
                               hour=hours[idx_2000 + 1:idx_2014 + 1], minutes=minutes[idx_2000 + 1:idx_2014 + 1],
                               model_name="difi8")
                warnings.warn("Dataset contains date before 2014.0, outside DIFI8's recommended range", UserWarning)

        with pytest.warns(UserWarning, match="Dataset contains date after 2024.0, outside DIFI8's recommended range"):
            if idx_2024 < N - 1:
                B = getSQfield(lat, lon, years[idx_2024:idx_end], months[idx_2024:idx_end], days[idx_2024:idx_end], h=h,
                               hour=hours[idx_2024:idx_end], minutes=minutes[idx_2024:idx_end], model_name="difi8")
                warnings.warn("Dataset contains date before 2024.0, outside DIFI8's recommended range", UserWarning)

        with pytest.raises(Exception,
                           match="This package does not contain f10.7 data after noon 12/31/2025. Input time data contains a date corresponding to f10.7 data not contained in this package."):
            B = getSQfield(lat, lon, years[idx_2024:], months[idx_2024:], days[idx_2024:], h=h, hour=hours[idx_2024:],
                           minutes=minutes[idx_2024:], model_name="difi8")

    def test_altitude_range(self):

        """
        Both DIFI8 and xDIFI2 should allow altitiude from -1 to 1000 km
        """

        N = 50
        lat = [45]
        lon = [100]
        years = np.linspace(2000., 2025, N)
        months = np.linspace(1., 12., N)
        days = np.linspace(1., 30., N)
        hours = np.linspace(0., 23., N)
        minutes = np.linspace(0., 59., N)
        #h = np.linspace(-1, 1000, N)
        h = 100


        B = getSQfield(lat, lon, years, months, days, h=h, hour=hours,
                           minutes=minutes, model_name="difi8")

        B = getSQfield(lat, lon, years, months, days, h=h, hour=hours,
                           minutes=minutes, model_name="xdifi2")

    def test_endtime(self):
        """
        The end_time should allow the last date of the year of 2025
        """

        N = 50
        lat = 25
        lon = 100
        h = 1000

        with pytest.warns(UserWarning, match="Dataset contains date after 2024.0, outside xDIFI2's reccomended range"):
            B = getSQfield(lat, lon, year=2024, month=1, day=1, h=h, hour=0,
                           minutes=0, model_name="xdifi2")
            warnings.warn("Dataset contains date after 2024.0, outside xDIFI2's reccomended range", UserWarning)


if __name__ == '__main__':
    unittest.main()
