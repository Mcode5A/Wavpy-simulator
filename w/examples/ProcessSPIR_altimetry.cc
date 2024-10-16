#include "../recipes_cpp/nr3.h"
#include "../ancillary_functions.h"
#include "../wavpy_global_variables.h"
#include "../waveform_power.h"
#include "../waveform_complex.h"
#include "../specular_geometry.h"
#include "../twoDim_planar_array.h"


using namespace std;

void ReadUndulationFile( MatDoub &undMat, VecDoub &undLats, VecDoub &undLongs );
double Interpol_und( double latitude, double longitude, MatDoub &undMat, VecDoub &undLats, VecDoub &undLongs );

//============================ PROGRAM =====================================
int main (int argc, char* argv[]){

	int i, prn;
	string spir_file_list, sp3file, posRfile, inertfile, outfile;
	FILE *fpo_data;
	MatDoub undMat;
	VecDoub undLats, undLongs;
	double **xx, **yy, **zz, **vxx, **vyy, **vzz;
	double *tt;
	double resultsT[6];  //[x, y, z, vx, vy, vz]
	int nsats, nlines_sat, sow_diff, ref_GPSweek;
	bool correctData, sat_vel;
	double sow_diffWeek;
	double **week_sow_XYZ;
	double resultsPos[6];
	int nlines_posR;
	double **week_sow_RollPitchYaw;
	double resultsInert[6];
	int nlines_inert;
	double sampling_rate = 80000000.0;
	Waveform_complex_cluster mywavCluster;
	TwoDim_planar_array UParray;
	TwoDim_planar_array DWarray;
	Specular_geometry mygeom;
	Waveform_power wav_data;
	Waveform_power wav_data_accum_1msec;
	Waveform_power wav_data_accum_10msec;
	wav_data.set_sampling_rate(sampling_rate);
	wav_data.set_min_resolution_fft_interp(0.05);
	wav_data.set_fit_length(5.0);
	wav_data_accum_1msec.set_sampling_rate(sampling_rate);
	wav_data_accum_1msec.set_min_resolution_fft_interp(0.05);
	wav_data_accum_1msec.set_fit_length(5.0);
	wav_data_accum_10msec.set_sampling_rate(sampling_rate);
	wav_data_accum_10msec.set_min_resolution_fft_interp(0.05);
	wav_data_accum_10msec.set_fit_length(5.0);
	double phases_offset[16];
	phases_offset[0] = 0.00;
	phases_offset[1] = 146.92;
	phases_offset[2] = -146.96;
	phases_offset[3] = -147.47;
	phases_offset[4] = -176.24;
	phases_offset[5] = -155.03;
	phases_offset[6] = -152.11;
	phases_offset[7] = 171.89;
	phases_offset[8] = 0.00;
	phases_offset[9] = 84.24;
	phases_offset[10] = 76.38;
	phases_offset[11] = 118.13;
	phases_offset[12] = 117.05;
	phases_offset[13] = 100.82;
	phases_offset[14] = -23.70;
	phases_offset[15] = 76.34;
	double pos_elems[16];
	pos_elems[0] = -0.1651;
	pos_elems[1] = -0.03808;
	pos_elems[2] = -0.1651;
	pos_elems[3] = 0.07622;
	pos_elems[4] = -0.0508;
	pos_elems[5] = -0.03808;
	pos_elems[6] = -0.0508;
	pos_elems[7] = 0.07622;
	pos_elems[8] = 0.0635;
	pos_elems[9] = -0.03808;
	pos_elems[10] = 0.0635;
	pos_elems[11] = 0.07622;
	pos_elems[12] = 0.1778;
	pos_elems[13] = -0.03808;
	pos_elems[14] = 0.1778;
	pos_elems[15] = 0.07622;
	UParray.set_antenna_elements_pos_AF(pos_elems, 8, 2, 0);
	DWarray.set_antenna_elements_pos_AF(pos_elems, 8, 2, 0);
	double uparray_BF_E_vec[3], uparray_BF_H_vec[3];
	uparray_BF_E_vec[0] = 1.0;
	uparray_BF_E_vec[1] = 0.0;
	uparray_BF_E_vec[2] = 0.0;
	uparray_BF_H_vec[0] = 0.0;
	uparray_BF_H_vec[1] = -1.0;
	uparray_BF_H_vec[2] = 0.0;
	UParray.set_antenna_orientation_BF_EH(uparray_BF_E_vec, uparray_BF_H_vec);
	double posUParrayBF[3], posDWarrayBF[3];
	posUParrayBF[0] = -1.21; 
	posUParrayBF[1] = 0.122; 
	posUParrayBF[2] = -0.01;
	posDWarrayBF[0] = -1.23; 
	posDWarrayBF[1] = 0.093; 
	posDWarrayBF[2] = 2.093;
	double retracking_1msec[999];
	double retracking_10msec[100];
	stringstream lineRead;
	string line;
	int sod;
	string spirfile;
	double posUParrayLocal[3], posDWarrayLocal[3];
	double posT_ECEF[3], posR_ECEF[3], posR_LLH[3], inert_RPY[3]; 
	double phasesUP[8], phasesDW[8], phase_delaysUP[8], phase_delaysDW[8];
	double undu;
	int week, msec;
	double ref_inert, inertial_delay, ref_atmos, atmospheric_delay, ref_range;
	double start_window_delay, start_window_delay_int10, sow;
	double int_wav_1msec[512];
	double int_wav_10msec[512];
	double ref_range_int10, ref_inert_int10, ref_atmos_int10, elev_int10, undu_int10;
	bool first_int10 = true;
	double apriori_scattdel;
	int filter_number;
	string char_GNSS;
	////////////////////////////////////////////////////////////////////////////////////////
	if(argc==8){
		spir_file_list.assign(argv[1]);
		sp3file.assign(argv[2]);
		posRfile.assign(argv[3]);
		inertfile.assign(argv[4]);
		char_GNSS.assign(argv[5]);
		prn = atoi(argv[6]);
		outfile.assign(argv[7]);
	}else{
		cout << "USAGE:  " << argv[0] << "  spir_file_list  sp3file  posRfile  inertfile  char_GNSS  PRN  outfile" << endl;
		return 1;
	}
	//Undulation file
	ReadUndulationFile(undMat, undLats, undLongs);
	//Satellite file
	CheckSP3File(sp3file.c_str(), nsats, nlines_sat, sow_diff);
	tt = (double *) malloc(nlines_sat*sizeof(double));
	xx = (double **) malloc(nlines_sat*sizeof(double *));
	yy = (double **) malloc(nlines_sat*sizeof(double *));
	zz = (double **) malloc(nlines_sat*sizeof(double *));
	vxx = (double **) malloc(nlines_sat*sizeof(double *));
	vyy = (double **) malloc(nlines_sat*sizeof(double *));
	vzz = (double **) malloc(nlines_sat*sizeof(double *));
	for(i=0; i<nlines_sat; i++){
		xx[i] = (double *) malloc(32*sizeof(double));
		yy[i] = (double *) malloc(32*sizeof(double));
		zz[i] = (double *) malloc(32*sizeof(double));
		vxx[i] = (double *) malloc(32*sizeof(double));
		vyy[i] = (double *) malloc(32*sizeof(double));
		vzz[i] = (double *) malloc(32*sizeof(double));
	}
	sat_vel = ReadSP3File(sp3file.c_str(), tt, xx, yy, zz, vxx, vyy, vzz, nlines_sat, sow_diff, ref_GPSweek, char(char_GNSS[0]));
	if(xx[0][prn-1]==999999.){
		cout << "ERROR! No data with PRN" << prn << " in .sp3 file" << endl;
		return 1;
	}
	//Receiver's position file
	nlines_posR = Check_week_sow_3values_file(posRfile.c_str());
	if(nlines_posR <= 0){
		cout << "ERROR! Empty or wrong week_sow_ECEFxyz file: " << posRfile << endl;
		return 1;
	}
	week_sow_XYZ = (double **) malloc(5*sizeof(double *));
	for(i=0; i<5; i++){
		week_sow_XYZ[i] = (double *) malloc(nlines_posR*sizeof(double));
	}
	Read_week_sow_3values_file(posRfile.c_str(), nlines_posR, week_sow_XYZ);
	//Receiver's inertials file
	nlines_inert = Check_week_sow_3values_file(inertfile.c_str());
	if(nlines_inert <= 0){
		cout << "ERROR! Empty or wrong week_sow_RollPitchYaw file: " << inertfile << endl;
		return 1;
	}
	week_sow_RollPitchYaw = (double **) malloc(5*sizeof(double *));
	for(i=0; i<5; i++){
		week_sow_RollPitchYaw[i] = (double *) malloc(nlines_inert*sizeof(double));
	}
	Read_week_sow_3values_file(inertfile.c_str(), nlines_inert, week_sow_RollPitchYaw);
	//
	if(char(char_GNSS[0]) == 'E'){
		apriori_scattdel = 25.0;
		filter_number = 1;
		wav_data.set_normtail_length(100.0);
	}else{
		apriori_scattdel = 10.0;
		filter_number = 0;
		wav_data.set_normtail_length(50.0);
	}
	////////////////////////////////////////////////////////////////////////////////////////
	//START PROCESS
	fpo_data = fopen(outfile.c_str(), "w");
	ifstream filelist(spir_file_list.c_str());
	if(!filelist){
		cout << "ERROR! Unable to read SPIR-files list file: " << spir_file_list << endl;
		return false;
	}
	if(filelist.is_open()){
		getline(filelist,line);
		while(!filelist.eof()){
			//lineRead << line;
			lineRead.clear();
			lineRead.str(line);
			lineRead >> sod >> spirfile;
			sow = double(sod) + 4.0*86400.0;
			cout << "At SoD:" << sod << " File:" << spirfile << endl;
			week = 1873;
			for(msec=0; msec<999; msec++){
				sow_diffWeek = sow + double(msec)/1000.0 + (week - ref_GPSweek)*604800.0;
				correctData = Interpol_sat_sp3(sow_diffWeek, prn, nlines_sat, resultsT, tt, xx, yy, zz, sat_vel, vxx, vyy, vzz);
				correctData = Interpol_ECEFpos_file(week, sow + double(msec)/1000.0, nlines_posR, resultsPos, week_sow_XYZ);
				correctData = Interpol_ECEFpos_file(week, sow + double(msec)/1000.0, nlines_inert, resultsInert, week_sow_RollPitchYaw);
				for(i=0; i<3; i++){
					posT_ECEF[i] = resultsT[i];
					posR_LLH[i] = resultsPos[i];
					inert_RPY[i] = resultsInert[i];
				}
				mygeom.set_ECEFpos_Transmitter(posT_ECEF);
				mygeom.set_LongLatHeight_Receiver(posR_LLH);
				mygeom.set_inertials(inert_RPY[0], inert_RPY[1], inert_RPY[2]);
				mygeom.compute_specular_point(0);
				mygeom.rotate_vector_BF_to_local(posUParrayBF, posUParrayLocal);
				mygeom.rotate_vector_BF_to_local(posDWarrayBF, posDWarrayLocal);
				inertial_delay = posUParrayLocal[1]*(-cos(mygeom.elevation*PI_NUM/180.0)) + posUParrayLocal[2]*(-sin(mygeom.elevation*PI_NUM/180.0)) + posDWarrayLocal[1]*cos(mygeom.elevation*PI_NUM/180.0) +  posDWarrayLocal[2]*(-sin(mygeom.elevation*PI_NUM/180.0));
				atmospheric_delay = (4.6/sin(mygeom.elevation*PI_NUM/180.0))*(1.0 - exp(-posR_LLH[2]/8.621));
				if(msec == 0){
					ref_inert = inertial_delay;
					ref_atmos = atmospheric_delay;
					ref_range = mygeom.geometric_delay*1000.0 - inertial_delay + atmospheric_delay;
					mygeom.get_ECEFpos_Receiver(posR_ECEF);
					UParray.compute_phase_delays_pos_ECEF_RT(inert_RPY, posR_ECEF, posT_ECEF);
					UParray.get_phase_delays(phase_delaysUP, 8);
					DWarray.compute_phase_delays_pos_ECEF_RT(inert_RPY, posR_ECEF, posT_ECEF);
					DWarray.get_phase_delays(phase_delaysDW, 8);
					for(i=0; i<8; i++){
						phasesUP[i] = (phase_delaysUP[i]*180.0/PI_NUM) - phases_offset[i];
						phasesDW[i] = (phase_delaysDW[i]*180.0/PI_NUM) - phases_offset[i + 8];
					}
					undu = Interpol_und(mygeom.latitudeS, mygeom.longitudeS, undMat, undLats, undLongs);
				}
				retracking_1msec[msec] = ref_range - (mygeom.geometric_delay*1000.0 - inertial_delay + atmospheric_delay);
				if(msec%10 == 0){
					retracking_10msec[msec/10] = ref_range - (mygeom.geometric_delay*1000.0 - inertial_delay + atmospheric_delay);
				}
			}
			start_window_delay = mywavCluster.load_ITF_waveforms_SPIR(spirfile.c_str(), ref_range, phasesUP, phasesDW, filter_number);
			mywavCluster.integrate_waveforms_retracking(1, sampling_rate, retracking_1msec, 999, int_wav_1msec, 512);
			//mywavCluster.integrate_waveforms(1, int_wav_1msec, 512);
			wav_data.set_init_range(start_window_delay);
			wav_data.set_waveform(int_wav_1msec, 512);
			//wav_data.compute_delays();
			wav_data.compute_delays_wlimits((ref_range - 16.5*2.0*sin(mygeom.elevation*PI_NUM/180.0)), 10.0, apriori_scattdel);
			fprintf(fpo_data, "%d 1 1 %f %f %f %f %f %f %f %f %f %f %f %f\n", sod, wav_data.positionDer, ref_range, ref_inert, ref_atmos, (ref_range - wav_data.positionDer)/(2.0*sin(mygeom.elevation*PI_NUM/180.0)), undu*1000.0, mygeom.elevation, wav_data.positionMax, wav_data.powerMax, wav_data.power_posDer, wav_data.floorNoise, wav_data.slope_normTail);
			mywavCluster.integrate_waveforms_retracking(10, sampling_rate, retracking_10msec, 99, int_wav_10msec, 512);
			//mywavCluster.integrate_waveforms(10, int_wav_10msec, 512);
			wav_data.set_init_range(start_window_delay);
			wav_data.set_waveform(int_wav_10msec, 512);
			//wav_data.compute_delays();
			wav_data.compute_delays_wlimits((ref_range - 16.5*2.0*sin(mygeom.elevation*PI_NUM/180.0)), 10.0, apriori_scattdel);
			fprintf(fpo_data, "%d 10 1 %f %f %f %f %f %f %f %f %f %f %f %f\n", sod, wav_data.positionDer, ref_range, ref_inert, ref_atmos, (ref_range - wav_data.positionDer)/(2.0*sin(mygeom.elevation*PI_NUM/180.0)), undu*1000.0, mygeom.elevation, wav_data.positionMax, wav_data.powerMax, wav_data.power_posDer, wav_data.floorNoise, wav_data.slope_normTail);
			if(sod%10 == 0){
				if(first_int10){
					first_int10 = false;
				}else{
					//wav_data_accum_1msec.compute_delays();
					wav_data_accum_1msec.compute_delays_wlimits((ref_range_int10 - 16.5*2.0*sin(elev_int10*PI_NUM/180.0)), 10.0, apriori_scattdel);
					fprintf(fpo_data, "%d 1 10 %f %f %f %f %f %f %f %f %f %f %f %f\n", sod - 10, wav_data_accum_1msec.positionDer, ref_range_int10, ref_inert_int10, ref_atmos_int10, (ref_range_int10 - wav_data_accum_1msec.positionDer)/(2.0*sin(elev_int10*PI_NUM/180.0)), undu_int10*1000.0, elev_int10, wav_data_accum_1msec.positionMax, wav_data_accum_1msec.powerMax, wav_data_accum_1msec.power_posDer, wav_data_accum_1msec.floorNoise, wav_data_accum_1msec.slope_normTail);
					//wav_data_accum_10msec.compute_delays();
					wav_data_accum_10msec.compute_delays_wlimits((ref_range_int10 - 16.5*2.0*sin(elev_int10*PI_NUM/180.0)), 10.0, apriori_scattdel);
					fprintf(fpo_data, "%d 10 10 %f %f %f %f %f %f %f %f %f %f %f %f\n", sod - 10, wav_data_accum_10msec.positionDer, ref_range_int10, ref_inert_int10, ref_atmos_int10, (ref_range_int10 - wav_data_accum_10msec.positionDer)/(2.0*sin(elev_int10*PI_NUM/180.0)), undu_int10*1000.0, elev_int10, wav_data_accum_10msec.positionMax, wav_data_accum_10msec.powerMax, wav_data_accum_10msec.power_posDer, wav_data_accum_10msec.floorNoise, wav_data_accum_10msec.slope_normTail);
				}
				ref_range_int10 = ref_range;
				ref_inert_int10 = ref_inert;
				ref_atmos_int10 = ref_atmos;
				elev_int10 = mygeom.elevation;
				undu_int10 = undu;
				start_window_delay_int10 = start_window_delay;
				wav_data_accum_1msec.set_init_range(start_window_delay);
				wav_data_accum_1msec.set_waveform(int_wav_1msec, 512);
				wav_data_accum_10msec.set_init_range(start_window_delay);
				wav_data_accum_10msec.set_waveform(int_wav_10msec, 512);
			}
			if((!first_int10)&&(sod%10 != 0)){
				wav_data_accum_1msec.add_waveform_retracking(int_wav_1msec, 512, ref_range_int10 - ref_range + start_window_delay - start_window_delay_int10, 0.1);
				//wav_data_accum_1msec.add_waveform_retracking(int_wav_1msec, 512, start_window_delay - start_window_delay_int10, 0.1);
				wav_data_accum_10msec.add_waveform_retracking(int_wav_10msec, 512, ref_range_int10 - ref_range + start_window_delay - start_window_delay_int10, 0.1);
				//wav_data_accum_10msec.add_waveform_retracking(int_wav_10msec, 512, start_window_delay - start_window_delay_int10, 0.1);
			}
			getline(filelist,line);
		}
		filelist.close();
	}
	////////////////////////////////////////////////////////////////////////////////////////
	fclose(fpo_data);
	free(tt);
	for(i=0; i<nlines_sat; i++){
		free(xx[i]);
		free(yy[i]);
		free(zz[i]);
		free(vxx[i]);
		free(vyy[i]);
		free(vzz[i]);
	}
	free(xx);
	free(yy);
	free(zz);
	free(vxx);
	free(vyy);
	free(vzz);
	for(i=0; i<5; i++){
		free(week_sow_XYZ[i]);
		free(week_sow_RollPitchYaw[i]);
	}
	free(week_sow_XYZ);
	free(week_sow_RollPitchYaw);
	return 0;
}


