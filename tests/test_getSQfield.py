import unittest
from pathlib import Path
import csv  # from Python stdlib

import numpy as np

from DIFI.getSQfield import getSQfield
from DIFI import jd2000_dt
from DIFI import get_f107_index_all
import warnings
import pytest

# The data files will be in the same directory as this test script
inputs_file_path = Path(__file__).parent / "test_values_xDIFI2_v1_20250528.txt"



class test_getSQfield(unittest.TestCase):

    def setUp(self):
        #Read in the test points as a list of dictionaries (one dict per each row)
        self.csvfields = ['lat','lon','year','month','day','hour','x','y','z']
        inputs_file = inputs_file_path.resolve()
        
        self.points = self.parse_geomagnetic_data(inputs_file)

    def parse_csv_files(self, inputs_file):

        points = []
        with open(inputs_file, "r") as fin:
            for row in csv.DictReader(fin, fieldnames=self.csvfields):
                # Will read in as str, convert to float
                point = {key: float(val) for key, val in row.items()}
                # Add in inputs that aren't in CSV
                point['minutes'] = 0.
                point['h'] = 0.
                points.append(point)

        # Pick a (repeatable) random subset so it doesn't take forever
        rng = np.random.default_rng(42)
        subset_inds = rng.integers(low=0, high=len(points), size=60).tolist()
        points = [points[ind] for ind in subset_inds]

        return points
  
    def get_f107(self,point):
        """Helper function to get f107 by date (test data doesn't specify)"""
        difi_t_f107, difi_f107 = get_f107_index_all.load_coefs()
        start_time = 5114.0
        end_time = float(difi_t_f107[-1])
        sq_t = jd2000_dt.jd2000_dt(point['year'], 
                                   point['month'], 
                                   point['day'], 
                                   point['hour'],
                                    point['minutes'])
        f107 = get_f107_index_all.get_f107_index(sq_t, start_time, end_time, difi_f107, difi_t_f107)
        return f107.tolist()[0] #returns an array and we just want float

    def parse_geomagnetic_data(self, filepath):
        # Define column labels
        keys = [
            "year", "month", "day", "hour", "minutes", "sec",
            "r", "theta", "phi", "h", "lat", "lon", "F10.7",
            "B_r", "B_theta", "B_phi", "x", "y", "z"
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

    def test_all_scalar(self):
        for point in self.points:
            Bxyz_test = getSQfield(lat=point['lat'],
                                lon=point['lon'],
                                year=point['year'],
                                month=point['month'],
                                day=point['day'],
                                hour=point['hour'],
                                minutes=point['minutes'],
                                h=point['h'],
                                f107_1=point["F10.7"])

            Bxyz_expected = {'X':point['x'],
                            'Y':point['y'],
                            'Z':point['z']}

            for cmpnt in ['X','Y','Z']:
                with self.subTest(f'B{cmpnt} {point}'):
                    B_test = Bxyz_test[cmpnt][0]
                    B_expect = Bxyz_expected[cmpnt]
                    self.assertAlmostEqual(B_test,B_expect,places=5)

    def test_all_lists(self):
        """Make sure we can handle all inputs as length 2 list"""
        point1,point2 = self.points[0],self.points[1]
        Bxyz_test = getSQfield(lat=[point1['lat'],point2['lat']],
                            lon=[point1['lon'],point2['lon']],
                            year=[point1['year'],point2['year']],
                            month=[point1['month'],point2['month']],
                            day=[point1['day'],point2['day']],
                            hour=[point1['hour'],point2['hour']],
                            minutes=[point1['minutes'],point2['minutes']],
                            h=[point1['h'],point2['h']],
                            f107_1=[point1['F10.7'],point2['F10.7']])

        Bxyz_expected = {'X':np.array([point1['x'],point2['x']]),
                        'Y':np.array([point1['y'],point2['y']]),
                        'Z':np.array([point1['z'],point2['z']])}

        for cmpnt in ['X','Y','Z']:
            with self.subTest(f'B{cmpnt} {point1} {point2}'):
                B_test = Bxyz_test[cmpnt].tolist()
                B_expect = Bxyz_expected[cmpnt].tolist()
                #point 1
                self.assertAlmostEqual(B_test[0],B_expect[0],places=6)
                #point 2
                self.assertAlmostEqual(B_test[1],B_expect[1],places=6)

    def test_inputs_scalar_except_one_list(self):
        """Check all possible permutations of scalars and one list
        (make sure all inputs can be a list when everything else is scalar)"""
        for i in range(len(self.points) // 3):
            point = self.points[i].copy() #make a copy so we don't overwrite original

            point['f107_1']=self.points[i]['F10.7']
            input_names = ['lat','lon','year','month','day','hour','minutes','h','f107_1']
            scalar_getSQfield_inputs = {key:point[key] for key in input_names}
            for input_name in input_names:
                with self.subTest(f'Test {input_name} is a list'):
                    getSQfield_inputs = scalar_getSQfield_inputs.copy()
                    
                    #Make one input a list
                    getSQfield_inputs[input_name] = [point[input_name],
                                                    point[input_name]]

                    Bxyz_test = getSQfield(**getSQfield_inputs)

                    #We expect to get 2 identical points back
                    Bxyz_expected = {'X':[point['x'],point['x']],
                                    'Y':[point['y'],point['y']],
                                    'Z':[point['z'],point['z']]}

                    for cmpnt in ['X','Y','Z']:
                        with self.subTest(f'B{cmpnt}'):
                            B_test = Bxyz_test[cmpnt]
                            B_expect = Bxyz_expected[cmpnt]
                            self.assertAlmostEqual(B_test[0],B_expect[0],places=5)
                            self.assertAlmostEqual(B_test[1],B_expect[1],places=5)

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

        idx_2001 = 0
        idx_2024 = N - 1
        idx_end = N - 1

        # Test for the point input
        for i in range(N):
            if years[i] < 2001.0:
                idx_2001 = i
                with pytest.warns(UserWarning, match="Dataset contains date before 2001.0, outside xDIFI2's reccomended range"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="xDIFI2")


            elif years[i] > 2024.0 and years[i] < 2026.0:
                if idx_2024 == N - 1:
                    idx_2024 = i
                with pytest.warns(UserWarning, match="Dataset contains date after 2024.0, outside xDIFI2's reccomended range"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="xDIFI2")

            elif years[i] >= 2026.0:
                if idx_end == N -1:
                    idx_end = i


                with pytest.raises(Exception, match="Request after year 2026, out of model's valid date range"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="xDIFI2")

        # Test for the list input
        with pytest.warns(UserWarning, match="Dataset contains date before 2001.0, outside xDIFI2's reccomended range"):
            if idx_2001 < N -1:
                B = getSQfield(lat, lon, years[:idx_2001+1], months[:idx_2001+1], days[:idx_2001+1], h=h, hour=hours[:idx_2001+1], minutes=minutes[:idx_2001+1], model_name="xDIFI2")

        with pytest.warns(UserWarning, match="Dataset contains date after 2024.0, outside xDIFI2's reccomended range"):
            if idx_2024 < N -1:
                B = getSQfield(lat, lon, years[idx_2024:idx_end], months[idx_2024:idx_end], days[idx_2024:idx_end], h=h, hour=hours[idx_2024:idx_end], minutes=minutes[idx_2024:idx_end], model_name="xDIFI2")

        with pytest.raises(Exception, match="Request after year 2026, out of model's valid date range"):
            B = getSQfield(lat, lon, years[idx_2024:], months[idx_2024:], days[idx_2024:], h=h, hour=hours[idx_2024:], minutes=minutes[idx_2024:], model_name="xDIFI2")

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

        idx_2014 = 0
        idx_2024 = N - 1
        idx_end = N - 1

        # Test for the point input
        for i in range(N):
            if years[i] < 2014.0:
                idx_2014 = i
                with pytest.warns(UserWarning,
                                  match="Dataset contains date before 2014.0, outside DIFI8's reccomended range"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="difi8")


            elif years[i] > 2024.0 and years[i] < 2026.0:
                if idx_2024 == N - 1:
                    idx_2024 = i
                with pytest.warns(UserWarning,
                                  match="Dataset contains date after 2024.0, outside DIFI8's reccomended range"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="difi8")

            elif years[i] >= 2026.0:
                if idx_end == N - 1:
                    idx_end = i

                with pytest.raises(Exception, match="Request after year 2026, out of model's valid date range"):
                    B = getSQfield(lat, lon, years[i], months[i], days[i], h=h, hour=0, minutes=0, model_name="difi8")

        # Test for the list input
        with pytest.warns(UserWarning, match="Dataset contains date before 2014.0, outside DIFI8's reccomended range"):
            if idx_2014 < N - 1:
                B = getSQfield(lat, lon, years[:idx_2014 + 1], months[:idx_2014 + 1], days[:idx_2014 + 1], h=h,
                               hour=hours[:idx_2014 + 1], minutes=minutes[:idx_2014 + 1], model_name="difi8")

        with pytest.warns(UserWarning, match="Dataset contains date after 2024.0, outside DIFI8's reccomended range"):
            if idx_2024 < N - 1:
                B = getSQfield(lat, lon, years[idx_2024:idx_end], months[idx_2024:idx_end], days[idx_2024:idx_end], h=h,
                               hour=hours[idx_2024:idx_end], minutes=minutes[idx_2024:idx_end], model_name="difi8")

        with pytest.raises(Exception, match="Request after year 2026, out of model's valid date range"):
            B = getSQfield(lat, lon, years[idx_2024:], months[idx_2024:], days[idx_2024:], h=h, hour=hours[idx_2024:],
                           minutes=minutes[idx_2024:], model_name="difi8")

    def test_altitude_range(self):

        """
        Both DIFI8 and xDIFI2 should allow altitiude from -1 to 1000 km
        """

        N = 50
        lat = 25
        lon = 100
        years = np.linspace(2000., 2025, N)
        months = np.linspace(1., 12., N)
        days = np.linspace(1., 31., N)
        hours = np.linspace(0., 23., N)
        minutes = np.linspace(0., 59., N)
        h = np.linspace(-1, 1000, N)

        for i in range(N):
            B = getSQfield(lat, lon, years, months, days, h=h[i], hour=hours,
                       minutes=minutes, model_name="difi8")

            B = getSQfield(lat, lon, years, months, days, h=h[i], hour=hours,
                       minutes=minutes, model_name="xdifi2")

    def test_endtime(self):
        """
        The end_time should allow the last date of the year of 2025
        """

        N = 50
        lat = 25
        lon = 100
        h = 1000


        B = getSQfield(lat, lon, year=2025, month=12, day=31, h=h, hour=23,
                       minutes=59, model_name="xdifi2")







if __name__ == '__main__':
    unittest.main()




                            
        
            
