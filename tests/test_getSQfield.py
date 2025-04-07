import unittest
from pathlib import Path
import csv  # from Python stdlib

import numpy as np

from DIFI.getSQfield import getSQfield

# The data files will be in the same directory as this test script
inputs_file_path = Path(__file__).parent / "20250319_outputs_inputs_difi_pypidifi7.csv"

class test_getSQfield(unittest.TestCase):

    def setUp(self):
        #Read in the test points as a list of dictionaries (one dict per each row)
        csvfields = ['lat','lon','year','month','day','hour','x','y','z']
        inputs_file = inputs_file_path.resolve()
        
        points = []
        with open(inputs_file, "r") as fin:
            for row in csv.DictReader(fin,fieldnames=csvfields):
                #Will read in as str, convert to float
                point = {key:float(val) for key,val in row.items()}
                #Add in inputs that aren't in CSV
                point['minutes']=0.
                point['h']=0.
                points.append(point)

        #Pick a (repeatable) random subset so it doesn't take forever
        rng = np.random.default_rng(42)
        subset_inds = rng.integers(low=0,high=len(points),size=100).tolist()
        self.points = [points[ind] for ind in subset_inds]
        

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
                                f107_1=None)

            Bxyz_expected = {'X':point['x'],
                            'Y':point['y'],
                            'Z':point['z']}

            for cmpnt in ['X','Y','Z']:
                with self.subTest(f'B{cmpnt} {point}'):
                    B_test = Bxyz_test[cmpnt][0]
                    B_expect = Bxyz_expected[cmpnt]
                    self.assertAlmostEqual(B_test,B_expect,places=6)

    def test_inputs_scalar_except_one_list(self):
        point = self.points[0]
        input_names = ['lat','lon','year','month','day','hour','minutes','h']
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
                        self.assertAlmostEqual(B_test[0],B_expect[0],places=6)
                        self.assertAlmostEqual(B_test[1],B_expect[1],places=6)
                        
        
            
