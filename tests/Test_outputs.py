import unittest
import os
import numpy as np

import output_results



class Test_branch_results(unittest.TestCase):

    def setUp(self):

        cur_dir = os.path.dirname(__file__)
        main_filename = "20250220_outputs_difi_main.csv"
        branch_filename = "20250226_outputs_difi_pypidifi7.csv"

        self.difi_main_file = os.path.join(cur_dir, main_filename)
        self.difi_py3_file = os.path.join(cur_dir, branch_filename)


    def read_file(self, filename):

        Bx, By, Bz = [], [], []

        with open(filename, "r") as file:

            for line in file:
                vals = line.split(",")
                Bx.append(vals[0])
                By.append(vals[1])
                Bz.append(vals[2])

        Bx = np.array(Bx, dtype=float)
        By = np.array(By, dtype=float)
        Bz = np.array(Bz, dtype=float)

        return Bx, By, Bz


    def test_outputs(self):

        #out_filename = "difimax_outputs.csv"
        #output_results.output_results(out_filename)

        py2_x, py2_y, py2_z = self.read_file(self.difi_main_file)
        py3_x, py3_y, py3_z = self.read_file(self.difi_py3_file)

        N = len(py2_x)

        # Define tolerance levels
        rtol = 1e-20
        atol = 1e-20
        max_e = 0

        for i in range(N):

            self.assertAlmostEqual(py2_x[i], py3_x[i], places=10)
            self.assertAlmostEqual(py2_y[i], py3_y[i], places=10)
            self.assertAlmostEqual(py2_z[i], py3_z[i], places=10)








