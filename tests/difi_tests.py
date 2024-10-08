import unittest
import numpy as np
import os

from DIFI import jd2000_dt

from DIFI import legendre as l
from DIFI import sun_md2000 as s
from DIFI import getmut as gm
from DIFI import geod2geoc as gg
from DIFI import SwarmL2_MIO_SHA_Read_v2 as sr

from DIFI import forward_Sq_d_Re as fs
from DIFI import design_SHA_Sq_i_Re_v2 as dsi
from DIFI import design_SHA_Sq_e_Re_v2 as dse
from DIFI import gg2gm_2010 as gg2gm

def raiseException(e):
    raise e

class TestDifi(unittest.TestCase):
    def test_jd2000_dt(self):
        inputs = {"year": 2024,
                  "month": 12,
                  "day": 31,
                  "UT": 12,
                  "minutes": 30
        }
        output = jd2000_dt.jd2000_dt(**inputs)

        print(output)
        self.assertAlmostEqual(output[0], 9.1315000e+03, places=1)
    def test_jd2000_dt_at0(self):        
        inputs = {"year": 2000, 
                  "month": 1, 
                  "day": 1
        }
        output = jd2000_dt.jd2000_dt(**inputs)
        self.assertEqual(output[0], 0)
    def test_jd2000_dt_string(self):
        inputs = {"year": 2000, 
                  "month": "1", 
                  "day": 1
        }
        output = jd2000_dt.jd2000_dt(**inputs)
        self.assertEqual(output[0], 0)
    def test_jd2000_dt_string_error(self):
        inputs = {"year": 2000, 
                  "month": "january", 
                  "day": 1
        }
        self.assertRaises(ValueError, jd2000_dt.jd2000_dt, **inputs) 
    def test_jd2000_dt_array(self):
        inputs = {"year": [2000, 2001, 2002, 2003, 2004, 2005], 
                  "month": 1, 
                  "day": 1}
        output = jd2000_dt.jd2000_dt(**inputs)
        self.assertTrue((output == np.array([0,365+1,2*365+1, 3*365+1, 4*365+1, 5*365+2])).all(), 
                         msg="")
        
    def test_legendre(self):
        [P, dP] = l.legendre(50, 12)
        self.assertAlmostEqual(P[0], 1, places=10)
        self.assertAlmostEqual(0.000000000000000, dP[0], places=12)
        self.assertAlmostEqual(-0.319004346471378, P[10], places=12)
        self.assertAlmostEqual(1.363673909157531, dP[10], places=12)
        self.assertAlmostEqual(0.076984407859974, P[20], places=12)
        self.assertAlmostEqual(-0.458732223204308, dP[20], places=12)
        self.assertAlmostEqual(-0.019666817104087, P[30], places=12)
        self.assertAlmostEqual(3.681114093113790, dP[30], places=12)
        self.assertAlmostEqual(0.570611823069099, P[40], places=12)
        self.assertAlmostEqual(0.757367737114840, dP[40], places=12)
        self.assertAlmostEqual(0.550902041938949, P[50], places=12)
        self.assertAlmostEqual(-0.294630595643946, dP[50], places=12)
        self.assertAlmostEqual(0.509169749102922, P[60], places=12)
        self.assertAlmostEqual(1.354234355166279, dP[60], places=12)
    
    def test_r2r(self):
        x = 5.5
        rads = s.r2r(x)
        self.assertAlmostEqual(rads,np.pi, places=12)
        x = -0.83
        rads = s.r2r(x)
        self.assertAlmostEqual(rads, np.pi*0.34, places=12)
    
    def test_sun_md2000(self):
        day = 1945.1
        rasc, decl = s.sun_md2000(day)
        self.assertAlmostEqual(0.636963379271627, rasc, places=15)
        self.assertAlmostEqual(0.252378643254093, decl, places=15)
        day = np.array([300.25522, 3671.167])
        rasc, decl = s.sun_md2000(day)
        self.assertAlmostEqual(-2.585172002406036, rasc[0], places=15)
        self.assertAlmostEqual(-1.027903722933005, rasc[1], places=15)
        self.assertAlmostEqual(-0.225089538635647, decl[0], places=15)
        self.assertAlmostEqual(-0.355440330381075, decl[1], places=15)
        
    def test_getmut(self):
        t = 1945.1
        theta = 120
        phi = 70
        t_mut = gm.getmut(np.array([t]), theta, phi)
        self.assertAlmostEqual(t_mut[0], 4.627250370601096, places=15)
        t_mut = gm.getmut(np.array([t]), theta, -phi)
        self.assertAlmostEqual(t_mut[0], 16.741613977999435, places=15)

    def test_geod2geoc(self):
        phi = 85*np.pi/180
        h = 0
        r, sph_phi = gg.geod2geoc(phi, h)
        self.assertAlmostEqual(r, 6356.9161142608, places=10)
        self.assertAlmostEqual(90-sph_phi*180.0/np.pi, 84.966475056615, places=12)
        r, sph_phi = gg.geod2geoc(-phi, h)
        self.assertAlmostEqual(90-sph_phi*180.0/np.pi, -84.966475056615, places=12)
        self.assertAlmostEqual(r, 6356.9161142608, places=10)
    
    def test_MIO_SHA_Read(self):
        baseDir = os.getcwd()+ "/"
        filename_DIFI = baseDir+'coefs/difi-coefs.txt'
        s = sr.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)
        self.assertEqual(s["nmax"], 60)
        self.assertEqual(s["mmax"], 12)
        self.assertAlmostEqual(s["theta_NGP"], 9.6900, places=5)
        self.assertAlmostEqual(s["phi_NGP"], 287.3800, places=4)
        self.assertAlmostEqual(s["N"], 0.01485, places=5)
        self.assertEqual(s['p_vec'], [0,1,2,3,4])
        self.assertEqual(s['s_vec'], [-2,-1,0,1,2])
        self.assertAlmostEqual(s['m_e_d_Re'][19][11], 0.0210697187, places=10)
        self.assertAlmostEqual(s['m_e_d_Re'][27][16], -0.0232533088, places=10)
        self.assertEqual(s['h'], 110)

    def test_forward_Sq_d_Re(self):

        a = 6371.2
        theta = np.pi/3
        lon = 140.5
        t = 6325.3
        f107 = 100
        baseDir = os.getcwd()+ "/"
        filename_DIFI = baseDir+'coefs/difi-coefs.txt'
        s = sr.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)
        s['h'] = 3
        [B1, B2] = fs.forward_Sq_d_Re(a, theta, lon, t, f107, s)
        self.assertAlmostEqual(B1[0][0], 7.627417313149250, places=12)
        self.assertAlmostEqual(B1[1][0], 5.066420402691515, places=12)
        self.assertAlmostEqual(B1[2][0], -28.54155907663003, places=12)
        self.assertAlmostEqual(B2[0][0], -0.522463251909774, places=12)
        self.assertAlmostEqual(B2[1][0], -1.246802486890230, places=12)
        self.assertAlmostEqual(B2[2][0], -5.746270416836315, places=12)

    def test_design_SHA_Sq_i_Re_v2(self):
        rho = np.array([1.2])
        theta = np.array([20])
        phi = np.array([50])
        t = 2020.5
        t_ut = 14
        p_vec = [0, 1, 2, 3, 4]
        s_vec = [-2, -1, 0, 1, 2]
        nmax = 60
        mmax = 12
        A_r, A_theta, A_phi = dsi.design_SHA_Sq_i_Re_v2(rho, theta, phi, t, t_ut, nmax, mmax, p_vec, s_vec)
        self.assertEqual(np.size(A_r), 68400)
        self.assertEqual(np.size(A_theta), 68400)
        self.assertEqual(np.size(A_phi), 68400)
        self.assertAlmostEqual(A_r[250][0], -0.331262203079535)
        self.assertAlmostEqual(A_r[0][0], 1.08760719998369)
        self.assertAlmostEqual(A_phi[0][0], 0)
        self.assertAlmostEqual(A_phi[41732][0],-0.00647992486071359)
        self.assertAlmostEqual(A_phi[41731][0], 0.0112235590879846)
        self.assertAlmostEqual(A_phi[55828][0], 0.000321462827923509)
        self.assertAlmostEqual(A_theta[0][0], 0.197928323683836)
        self.assertAlmostEqual(A_theta[55827][0], -0.000141047508325803)
        self.assertAlmostEqual(A_theta[55828][0], -0.000148227572514236)
    def test_design_SHA_Sq_e_Re_v2(self):
        rho = np.array([1.2])
        theta = np.array([20])
        phi = np.array([50])
        t = 2020.5
        t_ut = 14
        p_vec = [0, 1, 2, 3, 4]
        s_vec = [-2, -1, 0, 1, 2]
        nmax = 60
        mmax = 12
        A_r, A_theta, A_phi = dse.design_SHA_Sq_e_Re_v2(rho, theta, phi, t, t_ut, nmax, mmax, p_vec, s_vec)
        self.assertEqual(np.size(A_r), 68400)
        self.assertEqual(np.size(A_theta), 68400)
        self.assertEqual(np.size(A_phi), 68400)
        self.assertAlmostEqual(A_r[250][0], 127.886312992629)
        self.assertAlmostEqual(A_r[0][0], -0.939692620785908)
        self.assertAlmostEqual(A_phi[0][0], 0)
        self.assertAlmostEqual(A_phi[41732][0],-1308.31841762352)
        self.assertAlmostEqual(A_theta[0][0], 0.342020143325669)
        self.assertAlmostEqual(A_theta[1][0], -0.604022773555054)
        self.assertAlmostEqual(A_theta[55827][0], -14017.5227411039) #disagrees with matlab by factor -1

    def test_gg2gm(self):
        theta_gm, phi_gm, R = gg2gm.gg2gm_2010(theta_gg = np.array([20.0]), phi_gg = np.array([40.0]), get_R = True)
        self.assertAlmostEqual(theta_gm,25.396461447872113)
        self.assertAlmostEqual(phi_gm,132.4172298377015)
        self.assertAlmostEqual(R[0][0][0], 0.9283, places=4)
        self.assertAlmostEqual(R[0][1][0], -0.3719, places=4)
        self.assertAlmostEqual(R[0][0][1], 0.3719, places=4)
        self.assertAlmostEqual(R[0][1][1], 0.9283, places=4)

def main():
    unittest.main()

if __name__ == "__main__":
    main()
