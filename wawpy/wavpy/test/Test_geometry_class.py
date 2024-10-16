import numpy as np
import sys
import os
HOME = os.getenv("HOME")
wavpyDIR = HOME + "/gnssr_analysis/trunk/src/waveform_pylib/wavpy/"
sys.path.insert(1,wavpyDIR)
import wavpy

data_path = "../../src/data_share/"
receiver_file = data_path + "R_TXYZ_KM"
satellite_file = data_path + "igs16915.sp3"
inertials_file = data_path + "R_TINS";
mygeom = wavpy.Specular_geometry()
mygeom.set_Undulation(100.0)
LonLatHeight_R = np.array([114.558930, 22.481011, 0.121809])
mygeom.set_LongLatHeight_Receiver( LonLatHeight_R )
mygeom.dump_parameters()
pos_R = np.array([-2450.730138, 5363.012234, 2423.760669])
mygeom.set_ECEFpos_Receiver(pos_R)
LonLatHeight_T = np.array([114.558758, -10.481069, 20000.0])
mygeom.set_LongLatHeight_Transmitter(LonLatHeight_T)
mygeom.compute_specular_point(1)
mygeom.dump_parameters()
mygeom.read_ECEFpos_Receiver(receiver_file, 1691, 453295)
mygeom.read_ECEFpos_GNSS_Transmitter(satellite_file, 1691, 453295, 10, 'G')
mygeom.compute_specular_point(1)
mygeom.read_Inertials_Receiver(inertials_file, 1691, 453295)
mygeom.dump_parameters()
vector_r_a_BF = np.array([1.0, 0.0, 0.0], dtype=np.float32)
vector_r_t_BF = np.array([0.0, 1.0, 0.0], dtype=np.float32)
rvv = np.array([0.81447333, 0.05226387], dtype=np.float64)
rhh = np.array([-0.8186183, -0.05124843], dtype=np.float64)
windup_phase_R_L = mygeom.compute_Beyerle_windup_reflected(vector_r_a_BF, vector_r_t_BF, rvv, rhh, 1691, 453295)
vector_ant_BF = np.array([0.0, 0.0, 2.7], dtype=np.float32)
inertdel = mygeom.compute_inertial_delay(vector_ant_BF)
print "WINDUP reflected RHCP=%f LHCP=%f, INERTIAL DELAY %f" % (windup_phase_R_L[0]/(2*3.141592653589793238), windup_phase_R_L[1]/(2*3.141592653589793238), inertdel)
mygeom.set_Undulation(100.0)
mygeom.compute_specular_point(0)
mygeom.dump_parameters()
