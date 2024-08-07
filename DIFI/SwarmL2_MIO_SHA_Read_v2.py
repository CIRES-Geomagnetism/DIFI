# Partially autogenerated with SMOP 0.29
import numpy as np

def SwarmL2_MIO_SHA_Read_v2(filename):

    # s = SwarmL2_MIO_SHA_Read(filename)

    # Purpose:
# Read Swarm L2 MIO_SHA model from "*.DBL" file.

    # Input:
#    filename           name of file to read

    # Output: s is a dictionary with the following elements
#    m_e_d_Re(:,:)    primary Sq model (REAL coefficients)
#    m_i_d_Re(:,:)    secondary Sq model (REAL coefficients)
#    nmax, mmax       max degree and order in dipole coord.
#    p_vec(:)         diurnal wavenumbers
#    s_vec(:)         seasonal wavenumbers
#    theta_NGP        co-latitude of North Geomagnetic dipole [deg]
#    phi_NGP          longitude of North Geomagnetic dipole [deg]
#    h                altitude of Sq current layer [km]
#    N                Wold ratio of the F10.7 dependence

    # A. Chulliat, 2016-09-22
# (from an earlier version dated 2011-03-27, with inputs from N. Olsen)
# Python version rewritten by A. Woods, 2016-11-29
    s ={}
    #if exist(filename,r"file") == 0:

    try:
        fid=open(filename, 'r')
    except IOError as err:
        print("Coefficient file ",filename,r" does not exist!")
        raise err

# SwarmL2_MIO_SHA_Read_v2.m:29
    # skip comment lines

    tline=fid.readline()

# SwarmL2_MIO_SHA_Read_v2.m:33
    while tline[0]==r"#":
        tline=fid.readline()
# SwarmL2_MIO_SHA_Read_v2.m:35


    # read header part (1 line)
    # "%i %i %i %i %i %i %f %f %f %f"
    f=tline.split()
    #try:
# SwarmL2_MIO_SHA_Read_v2.m:40
    s["nmax"] = int(f[0])
    # SwarmL2_MIO_SHA_Read_v2.m:41
    s["mmax"] = int(f[1])
    # SwarmL2_MIO_SHA_Read_v2.m:42
    pmin= int(f[2])
    # SwarmL2_MIO_SHA_Read_v2.m:43
    pmax= int(f[3])
    # SwarmL2_MIO_SHA_Read_v2.m:44
    smin= int(f[4])
    # SwarmL2_MIO_SHA_Read_v2.m:45
    smax= int(f[5])
    # SwarmL2_MIO_SHA_Read_v2.m:46
    s["theta_NGP"] = float(f[6])
    # SwarmL2_MIO_SHA_Read_v2.m:47
    s["phi_NGP"] = float(f[7])
    # SwarmL2_MIO_SHA_Read_v2.m:48
    s["h"] = float(f[8])
    # SwarmL2_MIO_SHA_Read_v2.m:49
    s["N"] = float(f[9])
    # SwarmL2_MIO_SHA_Read_v2.m:50
    s["p_vec"] = range(pmin,pmax+1)
    # SwarmL2_MIO_SHA_Read_v2.m:52
    s["s_vec"] = range(smin,smax+1)
    # SwarmL2_MIO_SHA_Read_v2.m:53
    N_nm=s["mmax"]*(s["mmax"] + 2) + (s["nmax"] - s["mmax"])*(2*s["mmax"] + 1)
    # SwarmL2_MIO_SHA_Read_v2.m:55
    N_sp=len(s["p_vec"])*len(s["s_vec"])
    # SwarmL2_MIO_SHA_Read_v2.m:56
        # read primary model coefficients
    #except:
        #print tline
        #quit()
    s["m_e_d_Re"] = np.zeros([N_nm,2*N_sp])
# SwarmL2_MIO_SHA_Read_v2.m:60
    for i in range(0,N_nm):
        tline=fid.readline()
# SwarmL2_MIO_SHA_Read_v2.m:63
        if tline.split() is None:
            error(r"Wrong number of model coefficients!")
        f=[float(j) for j in tline.split()]

# SwarmL2_MIO_SHA_Read_v2.m:67
        s["m_e_d_Re"][i,:]=f[2:2*N_sp + 2]
# SwarmL2_MIO_SHA_Read_v2.m:68

    # read secondary model coefficients

    s["m_i_d_Re"] = np.zeros([N_nm,2*N_sp])
# SwarmL2_MIO_SHA_Read_v2.m:73
    for i in range(0,N_nm):
        tline=fid.readline()
# SwarmL2_MIO_SHA_Read_v2.m:76
        if tline.split() is []:
            error(r"Wrong number of model coefficients!")
        f=[float(j) for j in tline.split()]
# SwarmL2_MIO_SHA_Read_v2.m:80
        s["m_i_d_Re"][i,:]=f[2:2*N_sp + 2]
# SwarmL2_MIO_SHA_Read_v2.m:81

    fid.close()
    return s
