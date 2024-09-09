import os
from . import SwarmL2_F107_Read, SwarmL2_MIO_SHA_Read_v2

baseDir = os.path.dirname(__file__)

filename_DIFI = os.path.join(baseDir, "coefs", "difi-coefs.txt")
filename_f107 = os.path.join(baseDir, "coefs", "f107.DBL")
swarm_data = SwarmL2_MIO_SHA_Read_v2.SwarmL2_MIO_SHA_Read_v2(filename_DIFI)
difi_t_f107, difi_f107 = SwarmL2_F107_Read.SwarmL2_F107_Read(filename_f107)
