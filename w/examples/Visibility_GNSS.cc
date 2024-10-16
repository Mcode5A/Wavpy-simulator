#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include "recipes_cpp/nr3.h"
#include "ancillary_functions.h"
#include "specular_geometry.h"

using namespace std;

bool ReadSP3File( const char* namefile, VecDoub &tt, MatDoub &xx, MatDoub &yy, MatDoub &zz, MatDoub &vxx, MatDoub &vyy, MatDoub &vzz, int &nsats, int &weekGPS_ref, char gnss_identifier );
bool Interpol_sat_sp3( double sow, int prn, VecDoub &results, VecDoub &tt, MatDoub &xx, MatDoub &yy, MatDoub &zz, bool satvel, MatDoub &vxx, MatDoub &vyy, MatDoub &vzz );
void ReadUndulationFile( MatDoub &undMat, VecDoub &undLats, VecDoub &undLongs );
double Interpol_und( double latitude, double longitude, MatDoub &undMat, VecDoub &undLats, VecDoub &undLongs );

int main (int argc, char* argv[]) {
	double posR_ECEF[3];
	double LonLatHeight_R[3];
	double posT_ECEF[3];
	double lotlatHeight_R[3];
	//TIGRIS setup
	//posR_ECEF[0] = -2450.730138;
	//posR_ECEF[1] = 5363.012234;
	//posR_ECEF[2] = 2423.760669;
	//Svalbard
	//LonLatHeight_R[0] = 20.0;
	//LonLatHeight_R[1] = 80.0;
	//LonLatHeight_R[2] = 2.0;
	//Helsinki bay
	//LonLatHeight_R[0] = 26.0;
	//LonLatHeight_R[1] = 60.0;
	//LonLatHeight_R[2] = 3.0;
	//American tower, Dome-C
	posR_ECEF[0] = -903.842351442126;
	posR_ECEF[1] = 1375.99962379035;
	posR_ECEF[2] = -6144.73169128670;
	//Far de Sant Sebastià
	//LonLatHeight_R[0] = 3.0 + 12.0/60.0 + 10.15863/3600.0;
	//LonLatHeight_R[1] = 41.0 + 53.0/60.0 + 49.61786/3600.0;
	//LonLatHeight_R[2] = 178.354/1000.0;
	Specular_geometry specGeom;
	specGeom.set_ECEFpos_Receiver(posR_ECEF);
	//specGeom.set_LongLatHeight_Receiver(LonLatHeight_R);
	float elev_min = 10.0;
	float elev_max = 90.0;
	float azim_min = 270.0;
	float azim_max = 360.0;
	int interpol_time = 60;
	string satellite_file;
	int i, prn, iterations, number, nsats, ref_GPSweek, orbit_period, num_extra_sidays;
	double current_time;
	MatDoub xx, yy, zz, vxx, vyy, vzz;
	VecDoub tt;
	VecDoub resultsT(6);  //[x, y, z, vx, vy, vz]
	bool correctData, sat_vel;
	MatDoub undMat;
	VecDoub undLats, undLongs;
	double undulation;
	string char_GNSS;
	num_extra_sidays = 0;
	//Analysis of command line
	if(argc>=3){
		satellite_file.assign(argv[1]);
		char_GNSS.assign(argv[2]);
		if (argc==4){
			number = atoi(argv[3]);
			interpol_time = double(number);
		}
		if (argc==5){
			number = atoi(argv[3]);
			interpol_time = double(number);
			num_extra_sidays = atoi(argv[4]);
		}
	}else {
		cout << "USAGE : [1] " << argv[0] << "  sp3_file  char_GNSS  or  [2] " << argv[0] << "  sp3_file  char_GNSS  interpol_time  or  [3] " << argv[0] << "  sp3_file  G  interpol_time  num_extrapolation_sidays" << endl;
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
	if((argc==5)&&(char(char_GNSS[0]) != 'G')){
		cout << "Only GPS can be extrapolated" << endl;
		return -1;
	}
	orbit_period = 86164; //GPS
	//Satellite position
	sat_vel = ReadSP3File(satellite_file.c_str(), tt, xx, yy, zz, vxx, vyy, vzz, nsats, ref_GPSweek, char(char_GNSS[0]));  
	if(nsats == 0){
		cout << "ERROR! Empty or wrong GPS.sp3 file: " << satellite_file << endl;
		return -1;
	}
	//Undulation
	//ReadUndulationFile(undMat, undLats, undLongs);
	specGeom.set_Undulation(3.2250);
	//Number of iterations
	iterations = int((tt[tt.size()-1]-tt[0])/interpol_time);
	//And here... we... go!
	for(i=0; i<(iterations+1); i++){
		current_time = tt[0] + double(i*interpol_time) + interpol_time/2;
		if(current_time <= tt[tt.size()-1] + 900){
			for(prn=1; prn<33; prn++){
				if(xx[0][prn - 1] != 999999.){
					correctData = Interpol_sat_sp3(current_time, prn, resultsT, tt, xx, yy, zz, sat_vel, vxx, vyy, vzz);
					if(correctData){
						//printf("ECEF %f %d %f %f %f\n", current_time, prn, resultsT[0], resultsT[1], resultsT[2]);
						posT_ECEF[0] = resultsT[0];
						posT_ECEF[1] = resultsT[1];
						posT_ECEF[2] = resultsT[2];
						//specGeom.set_Undulation(0.0);
						specGeom.set_ECEFpos_Transmitter(posT_ECEF);
						specGeom.compute_specular_point(0);
						if((specGeom.elevation > elev_min)&&(specGeom.elevation < elev_max)&&(specGeom.azimuthT > azim_min)&&(specGeom.azimuthT < azim_max)){
							//undulation = Interpol_und(specGeom.latitudeS, specGeom.longitudeS, undMat, undLats, undLongs);
							//specGeom.set_Undulation(undulation);
							specGeom.compute_specular_point(0);
							specGeom.get_LongLatHeight_Receiver(lotlatHeight_R);
							//printf("%f %d %f %f %f %f\n", current_time, prn, specGeom.elevation, (1000.0*specGeom.geometric_delay), specGeom.longitudeS, specGeom.latitudeS);
							//printf("%f %d %f %f %f %f %f %f\n", fmod((current_time + double(num_extra_sidays*orbit_period)), 604800), prn, specGeom.elevation, specGeom.azimuthT, (1000.0*specGeom.geometric_delay), specGeom.longitudeS, specGeom.latitudeS, (1000.0*lotlatHeight_R[2]));
							printf("%f %d %f %f\n", fmod((current_time + double(num_extra_sidays*orbit_period)), 604800), prn, specGeom.elevation, specGeom.azimuthT);
						}
					}
				}
			}
		}
	}
	return 0;
}
