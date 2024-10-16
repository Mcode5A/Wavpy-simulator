#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include "specular_geometry.h"
#include "waveform_pylib_config.h"

using namespace std;

//============================ PROGRAM =====================================
int main (int argc, char* argv[]) {
	double pos_R[3], LonLatHeight_R[3], LonLatHeight_T[3];
	double vector_r_a_BF[3], vector_r_t_BF[3], rvv[2], rhh[2], windup_phase_R_L[2], vector_ant_BF[3];
	double inertdel;
	string data_path = LOCAL_DATA_PATH;
    cout << "LOCAL_DATA_PATH: " << LOCAL_DATA_PATH << endl;
	string receiver_file = data_path + "R_TXYZ_KM";
	string satellite_file = data_path + "igs16915.sp3";
	string inertials_file = data_path + "R_TINS";
	Specular_geometry mygeom;
	mygeom.set_Undulation(100.0);
	LonLatHeight_R[0] = 114.558930;
	LonLatHeight_R[1] = 22.481011;
	LonLatHeight_R[2] = 0.121809;
	mygeom.set_LongLatHeight_Receiver(LonLatHeight_R);
	mygeom.dump_parameters();
	pos_R[0] = -2450.730138;
	pos_R[1] = 5363.012234;
	pos_R[2] = 2423.760669;
	mygeom.set_ECEFpos_Receiver(pos_R);
	LonLatHeight_T[0] = 114.558758;
	LonLatHeight_T[1] = -10.481069;
	LonLatHeight_T[2] = 20000.0;
	mygeom.set_LongLatHeight_Transmitter(LonLatHeight_T);
	mygeom.compute_specular_point(1);
	mygeom.dump_parameters();
	mygeom.read_ECEFpos_Receiver(receiver_file.c_str(), 1691, 453295);
	mygeom.read_ECEFpos_GNSS_Transmitter(satellite_file.c_str(), 1691, 453295, 10, 'G');
	mygeom.compute_specular_point(1);
	mygeom.read_Inertials_Receiver(inertials_file.c_str(), 1691, 453295);
	mygeom.dump_parameters();
	vector_r_a_BF[0] = 1.0;
	vector_r_a_BF[1] = 0.0;
	vector_r_a_BF[2] = 0.0;
	vector_r_t_BF[0] = 0.0;
	vector_r_t_BF[1] = 1.0;
	vector_r_t_BF[2] = 0.0;
	rvv[0] = 0.81447333;
	rvv[1] = 0.05226387;
	rhh[0] = -0.8186183;
	rhh[1] = -0.05124843;
	mygeom.compute_Beyerle_windup_reflected(vector_r_a_BF, vector_r_t_BF, rvv, rhh, 1691, 453295, windup_phase_R_L);
	vector_ant_BF[0] = 0.0;
	vector_ant_BF[1] = 0.0;
	vector_ant_BF[2] = 2.7;
	inertdel = mygeom.compute_inertial_delay(vector_ant_BF);
	printf("WINDUP reflected RHCP=%f LHCP=%f, INERTIAL DELAY %f\n", windup_phase_R_L[0]/(2*3.141592653589793238), windup_phase_R_L[1]/(2*3.141592653589793238), inertdel);
	mygeom.set_Undulation(100.0);
	mygeom.compute_specular_point(0);
	mygeom.dump_parameters();
	return 0;
}
