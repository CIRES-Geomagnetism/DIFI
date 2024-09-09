import unittest
import os
import numpy as np

from geomaglib import util

from DIFI import SwarmL2_F107_Read, geod2geoc, legendre

class Test_sqfield(unittest.TestCase):

    def setUp(self):

        self.curr_dir = os.path.dirname(os.path.abspath(__file__))
        self.top_dir = os.path.dirname(self.curr_dir)
        self.f107_file = os.path.join(self.top_dir, "DIFI", "coefs", "f107.DBL")

    def test_get_end_time(self):
        # get end time from f107.DBL (the end of first column)
        difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(self.f107_file)

        self.assertEqual(difi_t_f107[-1], 9.1315000e+03)  # add assertion here

    def test_geod2geoc(self):

        # test Nils' vs Aranud's
        phi = 85 * np.pi / 180
        h = 0

        n_r, n_sph_phi = geod2geoc.geod2geoc(phi, h)
        n_theta = 90-n_sph_phi*180.0/np.pi

        a_r, a_theta = util.geod_to_geoc_lat(85, h)

        self.assertAlmostEqual(n_r, a_r, places=8)
        self.assertAlmostEqual(n_theta, a_theta, places=8)

    def test_legendre_for_sha_e(self):
        # test McupLow vs Collin's

        lat = 50.0

        self.assertAlmostEqual(True, False , places=6)

    def test_legendre_for_sha_i(self):
        # test McupLow vs Collin's

        lat = 50.0

        self.assertAlmostEqual(True, False , places=6)

    def test_getSQfield(self):

        # from examples.py vs new getSQfield
        self.assertAlmostEqual(True, False, places=6)






if __name__ == '__main__':
    unittest.main()
