import unittest
from pathlib import Path
import csv  # from Python stdlib

import numpy as np

from DIFI.getSQfield import getSQfield
from DIFI import jd2000_dt
from DIFI import get_f107_index

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
        subset_inds = rng.integers(low=0,high=len(points),size=60).tolist()
        self.points = [points[ind] for ind in subset_inds]
  
    def get_f107(self,point):
        """Helper function to get f107 by date (test data doesn't specify)"""
        start_time = 5114.0
        end_time = float(get_f107_index.difi_t_f107[-1])
        sq_t = jd2000_dt.jd2000_dt(point['year'], 
                                   point['month'], 
                                   point['day'], 
                                   point['hour'],
                                    point['minutes'])
        f107 = get_f107_index.get_f107_index(sq_t, start_time, end_time)
        return f107.tolist()[0] #returns an array and we just want float 
    
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
                            f107_1=None)

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
        for i in range(len(self.points)):
            point = self.points[i].copy() #make a copy so we don't overwrite original
            point['f107_1']=self.get_f107(point)
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
                            self.assertAlmostEqual(B_test[0],B_expect[0],places=6)
                            self.assertAlmostEqual(B_test[1],B_expect[1],places=6)
                            
        
            
