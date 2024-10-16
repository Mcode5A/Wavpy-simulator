#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include "ancillary_functions.h"
#include "specular_geometry.h"

using namespace std;

int main (int argc, char* argv[]) {
	double posR[3];
	double posT_ECEF[3];
	double **xx, **yy, **zz, **vxx, **vyy, **vzz;
	double *tt;
	double **week_sow_XYZ;
	double *undMat, *undLats, *undLongs;
	double resultsT[6];  //[x, y, z, vx, vy, vz]
	int i, prn, nsats, nlinesR, nlinesT, sow_diff, ref_GPSweek;
	bool correctData, sat_vel;
	Specular_geometry specGeom;
	float elev_min = 0.0;
	float elev_max = 90.0;
	float azim_min = 0.0;
	float azim_max = 360.0;
	string satellite_file;
	string receiver_file;
	double sow;
	string char_GNSS;
	//Analysis of command line
	if((argc==4)||(argc==8)||(argc==9)){
		receiver_file.assign(argv[1]);
		satellite_file.assign(argv[2]);
		char_GNSS.assign(argv[3]);
		if (argc>4){
			elev_min = atof(argv[4]);
			elev_max = atof(argv[5]);
			azim_min = atof(argv[6]);
			azim_max = atof(argv[7]);
		}
	}else{
		cout << "USAGE : [1] " << argv[0] << "  receiver_ascii_file  sp3_file  char_GNSS  or  [2] " << argv[0] << "    receiver_ascii_file  sp3_file  char_GNSS  elev_min  elev_max  azim_min  azim_max  or  [3] " << argv[0] << "    receiver_ascii_file  sp3_file  char_GNSS  elev_min  elev_max  azim_min  azim_max  1 (for LonLatH-rec-file)" << endl;
		return -1;
	}
	if((char(char_GNSS[0]) != 'G')&&(char(char_GNSS[0]) != 'E')&&(char(char_GNSS[0]) != 'C')&&(char(char_GNSS[0]) != 'J')){
		cout << "Read  char_GNSS: " << char_GNSS << endl;
		cout << "Valid char_GNSS:  G (GPS),  E (Galileo),  C (BeiDou),  J (QZSS)" << endl;
		return -1;
	}
	if(satellite_file.size()<5){
		cout << "Incorrect satellite file" << endl;
		return -1;
	}
	if(satellite_file.substr(satellite_file.size()-4,4) != ".sp3"){
		cout << "Incorrect satellite file" << endl;
		return -1;
	}
	//Read sp3 file
	CheckSP3File(satellite_file.c_str(), nsats, nlinesT, sow_diff);
	if((nsats == 0)||(nlinesT == 0)){
		cout << "ERROR! Empty or wrong GPS.sp3 file: " << satellite_file << endl;
		return -1;
	}
	undMat = (double*) malloc (681*1440*sizeof(double));
	undLats = (double*) malloc (681*sizeof(double));
	undLongs = (double*) malloc (1440*sizeof(double));
	ReadUndulationFile(undMat, undLats, undLongs);
	tt = (double *) malloc(nlinesT*sizeof(double));
	xx = (double **) malloc(nlinesT*sizeof(double *));
	yy = (double **) malloc(nlinesT*sizeof(double *));
	zz = (double **) malloc(nlinesT*sizeof(double *));
	vxx = (double **) malloc(nlinesT*sizeof(double *));
	vyy = (double **) malloc(nlinesT*sizeof(double *));
	vzz = (double **) malloc(nlinesT*sizeof(double *));
	for(i=0; i<nlinesT; i++){
		xx[i] = (double *) malloc(32*sizeof(double));
		yy[i] = (double *) malloc(32*sizeof(double));
		zz[i] = (double *) malloc(32*sizeof(double));
		vxx[i] = (double *) malloc(32*sizeof(double));
		vyy[i] = (double *) malloc(32*sizeof(double));
		vzz[i] = (double *) malloc(32*sizeof(double));
	}
	sat_vel = ReadSP3File(satellite_file.c_str(), tt, xx, yy, zz, vxx, vyy, vzz, nlinesT, sow_diff, ref_GPSweek, char(char_GNSS[0]));
	//Read receiver's file
	nlinesR = Check_week_sow_3values_file(receiver_file.c_str());
	if(nlinesR <= 0){
		cout << "ERROR! Empty or wrong week_sow_PosR file: " << receiver_file << endl;
		return -1;
	}
	week_sow_XYZ = (double **) malloc(5*sizeof(double *));
	for(i=0; i<5; i++){
		week_sow_XYZ[i] = (double *) malloc(nlinesR*sizeof(double));
	}
	Read_week_sow_3values_file(receiver_file.c_str(), nlinesR, week_sow_XYZ);
	//Undulation
	specGeom.set_Undulation(0.0);
	//specGeom.set_Undulation(3.2250);
	//
	for(i=0; i<nlinesR; i++){
		sow = week_sow_XYZ[1][i];
		posR[0] = week_sow_XYZ[2][i];
		posR[1] = week_sow_XYZ[3][i];
		posR[2] = week_sow_XYZ[4][i];
		if(argc != 9){
			specGeom.set_ECEFpos_Receiver(posR);
		}else{
			specGeom.set_LongLatHeight_Receiver(posR);
		}
		if((week_sow_XYZ[0][i] == ref_GPSweek)&&(sow >= tt[0])&&(sow <= tt[nlinesT - 1])){
			for(prn=1; prn<33; prn++){
				if(xx[0][prn - 1] != 999999.){
					correctData = Interpol_sat_sp3(sow, prn, nlinesT, resultsT, tt, xx, yy, zz, sat_vel, vxx, vyy, vzz);
					if(correctData){
						posT_ECEF[0] = resultsT[0];
						posT_ECEF[1] = resultsT[1];
						posT_ECEF[2] = resultsT[2];
						specGeom.set_ECEFpos_Transmitter(posT_ECEF);
						specGeom.compute_specular_point(0);
						//undulation = Interpol_und(specGeom.longitudeS, specGeom.latitudeS, undMat, undLats, undLongs);
						//specGeom.set_Undulation(undulation);
						//specGeom.compute_specular_point(0);
						if((specGeom.elevation >= elev_min)&&(specGeom.elevation <= elev_max)&&(specGeom.azimuthT >= azim_min)&&(specGeom.azimuthT <= azim_max)){
							printf("%f %d %f %f %f %f\n", fmod(sow, 86400), prn, specGeom.elevation, specGeom.azimuthT, specGeom.latitudeS, specGeom.longitudeS);
						}
					}
				}
			}
		}
	}
	free(undMat);
	free(undLats);
	free(undLongs);
	free(tt);
	free(xx);
	free(yy);
	free(zz);
	free(vxx);
	free(vyy);
	free(vzz);
	free(week_sow_XYZ);
	return 0;
}
