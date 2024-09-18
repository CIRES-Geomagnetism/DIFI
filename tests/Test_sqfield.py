import unittest
import os
import numpy as np

from geomaglib import util

from DIFI import SwarmL2_F107_Read, geod2geoc, legendre, jd2000_dt, get_f107_index, forward_Sq_d_Re, gg2gm_2010, design_SHA_Sq_e_Re_v2, design_SHA_Sq_i_Re_v2, getSQfield


class Test_sqfield(unittest.TestCase):

    def setUp(self):

        self.curr_dir = os.path.dirname(os.path.abspath(__file__))
        self.top_dir = os.path.dirname(self.curr_dir)
        self.f107_file = os.path.join(self.top_dir, "DIFI", "coefs", "f107.DBL")

        self.difi_t_f107, self.difi_f107 = get_f107_index.load_coefs()
        self.swarm_data = get_f107_index.load_swarm()

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

    def test_gg2gm(self):

        # Test the results for latitude and colatitude
        lat, lon, h = 44.67, 55.5, 0
        lat_rad = lat * np.pi / 180
        r_gc, theta_gc = geod2geoc.geod2geoc(lat_rad, h)
        theta_gc = theta_gc * 180.0/np.pi

        a_r, a_theta = util.geod_to_geoc_lat(lat, h)
        a_cotheta = 90.0 - a_theta

        theta_d, phi_d, rotmat_d = gg2gm_2010.gg2gm_2010(np.array([theta_gc]), np.array([lon]), get_R=True)
        print(np.array([theta_gc]))
        theta_a, phi_a, rotmat_a = gg2gm_2010.gg2gm_2010(np.array([a_cotheta]), np.array([lon]), get_R=True)

        self.assertAlmostEqual(theta_d[0], theta_a[0], places=8)
        self.assertAlmostEqual(phi_d[0], phi_a[0], places=8)

    def test_legendre_for_sha_e(self):
        # test McupLow vs Collin's

        lat, lon, h = 44.67, 55.5, 0
        lat_rad = lat * np.pi / 180

        a = 6371.2
        r_gc, theta_gc = geod2geoc.geod2geoc(lat_rad, h)
        rho = r_gc / a
        theta_d, phi_d, rotmat = gg2gm_2010.gg2gm_2010(np.array([theta_gc]), np.array([lon]), get_R=True)

        print(f"nmax: {self.swarm_data['nmax']} mmax:{self.swarm_data['mmax']}")

        arr_internal = np.array(design_SHA_Sq_e_Re_v2.design_SHA_Sq_e_Re_v2(
            rho,
            theta_d,
            phi_d,
            self.swarm_data['nmax'],
            self.swarm_data['mmax'],
        ))

        print(arr_internal)

        #self.assertAlmostEqual(True, False , places=6)

    def test_legendre_for_sha_i(self):
        # test McupLow vs Collin's

        lat = 50.0

        self.assertAlmostEqual(True, False , places=6)


    def test_geod_to_geoc_lat_for_forward_sq(self):

        # replacing geod2geoc at forward_Sq_d_re with geomaglib.util.geod_to_geoc_lat()
        a = 6371.2
        start_time = 5114.0
        end_time = float(self.difi_t_f107[-1])
        year, month, day, hour, minutes = 2024, 9, 5, 0, 0
        sq_t = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)
        f107_1 = get_f107_index.get_f107_index(sq_t, start_time, end_time)
        lat, lon, h = 44.67, 55.5, 0
        r_gc, theta_gc = geod2geoc.geod2geoc(np.radians(lat), h)
        r_gc = a + h
        theta_gc = np.degrees(theta_gc)
            # print "Difi input", r_gc, theta_gc, RV['lon'], sq_t, f107_1

        # theta from geod2geoc
        B_1, B_2 = forward_Sq_d_Re.forward_Sq_d_Re(
            r_gc,
            theta_gc,
            lon,
            sq_t,
            f107_1,
            self.swarm_data,
        )

        # theta from geod_to_geoc_lat



        self.assertAlmostEqual()


    def test_getSQfield(self):

        # from examples.py vs new getSQfield
        self.assertAlmostEqual(True, False, places=6)






if __name__ == '__main__':
    unittest.main()
