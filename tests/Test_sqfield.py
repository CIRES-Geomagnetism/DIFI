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

    def test_jd2000_dt(self):
        inputs = {"year": 2025,
                  "month": 12,
                  "day": 31,
                  "UT": 12,
                  "minutes": 30
                  }
        output = jd2000_dt.jd2000_dt(**inputs)

        print(output)
        self.assertAlmostEqual(output[0], 9.4965000e+03 , places=1)
    def test_get_end_time(self):
        # get end time from f107.DBL (the end of first column)
        difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(self.f107_file)

        self.assertEqual(difi_t_f107[-1], 9.1315000e+03)  # add assertion here

    def test_geod2geoc(self):

        # test Nir's vs Aranud's
        phi = 85 * np.pi / 180
        h = 0

        n_r, n_sph_phi = geod2geoc.geod2geoc(phi, h)
        n_theta = 90-n_sph_phi*180.0/np.pi

        a_r, a_theta = util.geod_to_geoc_lat(85, h)
        a_to_n_theta = (90 - a_theta)*np.pi/180.0

        self.assertAlmostEqual(n_r, a_r, places=8)
        self.assertAlmostEqual(n_theta, a_theta, places=8)
        self.assertAlmostEqual(n_theta, a_theta, places=8)
        self.assertAlmostEqual(n_sph_phi, a_to_n_theta, places=8)

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

        arr_internal1 = np.array(design_SHA_Sq_e_Re_v2.design_SHA_Sq_e_Re_v2(
            rho,
            theta_d,
            phi_d,
            self.swarm_data['nmax'],
            self.swarm_data['mmax'],
        ))

        arr_internal2 = np.array(design_SHA_Sq_e_Re_v2.design_SHA_Sq_e_Re_v1(
            rho,
            theta_d,
            phi_d,
            self.swarm_data['nmax'],
            self.swarm_data['mmax'],
        ))

        N = len(arr_internal1)
        M = len(arr_internal1[0])

        for i in range(N):
            for j in range(M):

                self.assertAlmostEqual(arr_internal1[i][j][0], arr_internal2[i][j][0], places=6)

    def test_legendre_for_sha_i(self):
        # test McupLow vs Collin's

        lat, lon, h = 44.67, 55.5, 0
        lat_rad = lat * np.pi / 180

        a = 6371.2
        r_gc, theta_gc = geod2geoc.geod2geoc(lat_rad, h)
        rho = r_gc / a
        theta_d, phi_d, rotmat = gg2gm_2010.gg2gm_2010(np.array([theta_gc]), np.array([lon]), get_R=True)

        print(f"nmax: {self.swarm_data['nmax']} mmax:{self.swarm_data['mmax']}")

        arr_internal1 = np.array(design_SHA_Sq_i_Re_v2.design_SHA_Sq_i_Re_v2(
            rho,
            theta_d,
            phi_d,
            self.swarm_data['nmax'],
            self.swarm_data['mmax'],
        ))

        arr_internal2 = np.array(design_SHA_Sq_i_Re_v2.design_SHA_Sq_i_Re_v1(
            rho,
            theta_d,
            phi_d,
            self.swarm_data['nmax'],
            self.swarm_data['mmax'],
        ))

        print(np.shape(arr_internal1))
        N = len(arr_internal1)
        M = len(arr_internal1[0])

        for i in range(N):
            for j in range(M):

                self.assertAlmostEqual(arr_internal1[i][j][0], arr_internal2[i][j][0], places=6)

    def test_geod_to_geoc_lat_for_forward_sq(self):

        # replacing geod2geoc at forward_Sq_d_re with geomaglib.util.geod_to_geoc_lat()
        a = 6371.2
        start_time = 5114.0
        end_time = float(self.difi_t_f107[-1])
        year, month, day, hour, minutes = 2024, 9, 5, 0, 0
        sq_t = jd2000_dt.jd2000_dt(year, month, day, hour, minutes)
        f107_1 = get_f107_index.get_f107_index(sq_t, start_time, end_time)




        lats = np.linspace(57.75,77.75,10)
        print(lats)
        lon = 100.34
        h = 33.55

        for i in range(len(lats)):
            lat = lats[i]
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
            r_gc, theta_gc2 = util.geod_to_geoc_lat(lat, h)
            r_gc = a + h
            cotheta_gc = 90 - theta_gc2


            B_12, B_22 = forward_Sq_d_Re.forward_Sq_d_Re(
                r_gc,
                cotheta_gc,
                lon,
                sq_t,
                f107_1,
                self.swarm_data,
            )

            for i in range(3):
                self.assertAlmostEqual(B_1[i][0], B_12[i][0])
                self.assertAlmostEqual(B_2[i][0], B_22[i][0])






    def test_getSQfield(self):

        year = 2023
        month = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        day = 15

        B = getSQfield(21.3166, -157.9996, 2023, 11, 1)

        print(B)

        B = getSQfield(20, -50, year, month, day)

        print(B)

        Z = [-1.11753492, -2.17806267, -2.36170017, -1.55958987, -0.4729657 ,
        0.02314406, -0.62739855, -1.84115162, -2.41599221, -1.93664788,
       -0.91328872, -0.47807066]

        X = [-1.63857938, -0.99863218, -0.70018449, -1.01992701, -1.70769538,
       -2.239885  , -2.03735071, -1.2146875 , -0.53504296, -0.62733946,
       -1.33669145, -1.80754398]

        Y = [-0.53981766,  0.02365691,  1.0244612 ,  2.68627321,  4.66764832,
        6.15517171,  5.83115627,  3.98613964,  2.00863353,  0.73887103,
       -0.06879104, -0.51973947]

        for i in range(12):
            self.assertAlmostEqual(B["X"][i], X[i], places=6)
            self.assertAlmostEqual(B["Y"][i], Y[i], places=6)
            self.assertAlmostEqual(B["Z"][i], Z[i], places=6)











if __name__ == '__main__':
    unittest.main()
