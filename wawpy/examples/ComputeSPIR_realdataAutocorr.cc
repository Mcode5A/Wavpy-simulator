#include "ancillary_functions.h"
#include "wavpy_global_variables.h"
#include "waveform_power.h"
#include "waveform_complex.h"
#include "specular_geometry.h"
#include "rf_front_end.h"

using namespace std;

//============================ PROGRAM =====================================
int main (int argc, char* argv[]){
	int size_wav = 512;
	int i, prn;
	string spir_file_list, sp3file, posRfile, inertfile, outfile;
	FILE *fpo_data;
	double **xx, **yy, **zz, **vxx, **vyy, **vzz;
	double *tt;
	double resultsT[6];  //[x, y, z, vx, vy, vz]
	int nsats, nlines_sat, sow_diff, ref_GPSweek;
	bool sat_vel;
	double sow_diffWeek;
	double **week_sow_XYZ;
	double resultsPos[6];
	int nlines_posR;
	double **week_sow_RollPitchYaw;
	double resultsInert[6];
	int nlines_inert;
	double sampling_rate = 80000000.0;
	Waveform_complex_cluster mywavCluster;
	Waveform_power wavint_data;
	wavint_data.set_sampling_rate(sampling_rate);
	RF_FrontEnd myreceiver;
	Specular_geometry mygeom;
	double phases_offset[16];
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
	myreceiver.set_antenna_elements_pos_AF(pos_elems, 8, 2, 0);
	double uparray_BF_E_vec[3], uparray_BF_H_vec[3];
	uparray_BF_E_vec[0] = 1.0;
	uparray_BF_E_vec[1] = 0.0;
	uparray_BF_E_vec[2] = 0.0;
	uparray_BF_H_vec[0] = 0.0;
	uparray_BF_H_vec[1] = -1.0;
	uparray_BF_H_vec[2] = 0.0;
	myreceiver.set_antenna_orientation_BF_EH(uparray_BF_E_vec, uparray_BF_H_vec);
	stringstream lineRead;
	string line;
	int sod;
	string spirfile;
	double posT_ECEF[3], posR_ECEF[3], posR_LLH[3], inert_RPY[3]; 
	double phases_BF[16], phase_delaysUP[8];
	int week, msec;
	double start_window_delay, sow;
	int filter_number;
	string char_GNSS;
	//NEW VARIABLES
	int coh_msecs;
	string freq_band;
	double int_wav[size_wav], range_wav[size_wav];
	signed char elements_1[16], elements_2[16];
	for(i=0; i<16; i++){
		elements_1[i] = 0;
		elements_2[i] = 0;
	}
	int combination_type;
	////////////////////////////////////////////////////////////////////////////////////////
	if(argc==11){
		coh_msecs = atoi(argv[1]);
		combination_type = atoi(argv[2]);
		spir_file_list.assign(argv[3]);
		sp3file.assign(argv[4]);
		posRfile.assign(argv[5]);
		inertfile.assign(argv[6]);
		char_GNSS.assign(argv[7]);
		freq_band.assign(argv[8]);
		prn = atoi(argv[9]);
		outfile.assign(argv[10]);
	}else{
		cout << "USAGE:  " << argv[0] << " coh_msecs  combination_type  spir_file_list  sp3file  posRfile  inertfile  char_GNSS  freq_band  PRN  outfile" << endl;
		return 1;
	}
	//Check input vales
	if((coh_msecs <= 0)||(coh_msecs > 999)){
		cout << "ERROR! Wrong value for coh_msecs: " << coh_msecs << endl;
		return 1;
	}
	switch(combination_type){
		case 1: {
			elements_1[0] = 1;
			elements_1[1] = 1;
			elements_1[2] = 1;
			elements_1[3] = 1;
			elements_2[4] = 1;
			elements_2[5] = 1;
			elements_2[6] = 1;
			elements_2[7] = 1;
			break;
		}
		case 2: {
			elements_1[0] = 1;
			elements_1[1] = 1;
			elements_1[4] = 1;
			elements_1[5] = 1;
			elements_2[2] = 1;
			elements_2[3] = 1;
			elements_2[6] = 1;
			elements_2[7] = 1;
			break;
		}
		case 3: {
			elements_1[0] = 1;
			elements_1[1] = 1;
			elements_1[6] = 1;
			elements_1[7] = 1;
			elements_2[2] = 1;
			elements_2[3] = 1;
			elements_2[4] = 1;
			elements_2[5] = 1;
			break;
		}
		case 4: {
			elements_1[0] = 1;
			elements_1[2] = 1;
			elements_1[4] = 1;
			elements_1[6] = 1;
			elements_2[1] = 1;
			elements_2[3] = 1;
			elements_2[5] = 1;
			elements_2[7] = 1;
			break;
		}
		default: {
			cout << "ERROR! Wrong value for combination_type: " << coh_msecs << endl;
			return 1;
			break;
		}
	}
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
	if((freq_band[0] == 'L') && (freq_band[1] == '5')){
		filter_number = 2;
		myreceiver.set_frequency(FREQ_GPS_L5);
		phases_offset[0] = 0.00;
		phases_offset[1] = 111.88;
		phases_offset[2] = 145.78;
		phases_offset[3] = 146.22;
		phases_offset[4] = 140.00;
		phases_offset[5] = 139.11;
		phases_offset[6] = 171.94;
		phases_offset[7] = 129.75;
		phases_offset[8] = 0.00;
		phases_offset[9] = 154.58;
		phases_offset[10] = 130.89;
		phases_offset[11] = 125.09;
		phases_offset[12] = 103.49;
		phases_offset[13] = 111.18;
		phases_offset[14] = -0.85;
		phases_offset[15] = 82.12;
	}else{
		if(char(char_GNSS[0]) == 'E'){
			filter_number = 1;
		}else{
			filter_number = 0;
		}
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
	}
	////////////////////////////////////////////////////////////////////////////////////////
	//START PROCESS
	if((fpo_data = fopen(outfile.c_str(), "wb")) == NULL){
		printf("ERROR! outfile can't be opened.\n");
		return 1;
	}
	ifstream filelist(spir_file_list.c_str());
	if(!filelist){
		cout << "ERROR! Unable to read SPIR-files list file: " << spir_file_list << endl;
		return 1;
	}
	if(filelist.is_open()){
		getline(filelist,line);
		while(!filelist.eof()){
			lineRead.clear();
			lineRead.str(line);
			lineRead >> sod >> spirfile;
			sow = double(sod) + 4.0*86400.0;
			cout << "At SoD:" << sod << " File:" << spirfile << endl;
			week = 1873;
			msec = 500;
			sow_diffWeek = sow + double(msec)/1000.0 + (week - ref_GPSweek)*604800.0;
			Interpol_sat_sp3(sow_diffWeek, prn, nlines_sat, resultsT, tt, xx, yy, zz, sat_vel, vxx, vyy, vzz);
			Interpol_ECEFpos_file(week, sow + double(msec)/1000.0, nlines_posR, resultsPos, week_sow_XYZ);
			Interpol_ECEFpos_file(week, sow + double(msec)/1000.0, nlines_inert, resultsInert, week_sow_RollPitchYaw);
			for(i=0; i<3; i++){
				posT_ECEF[i] = resultsT[i];
				posR_LLH[i] = resultsPos[i];
				inert_RPY[i] = resultsInert[i];
			}
			mygeom.set_LongLatHeight_Receiver(posR_LLH);
			mygeom.get_ECEFpos_Receiver(posR_ECEF);
			myreceiver.compute_phase_delays_pos_ECEF_RT(inert_RPY, posR_ECEF, posT_ECEF);
			myreceiver.get_phase_delays(phase_delaysUP, 8);
			for(i=0; i<16; i++){
				if(i<8){
					phases_BF[i] = (phase_delaysUP[i]*180.0/PI_NUM) - phases_offset[i];
				}else{
					phases_BF[i] = 0.0;
				}
			}
			start_window_delay = mywavCluster.load_ITF_waveforms_SPIR_selected_signals(spirfile.c_str(), 0.0, phases_BF, elements_1, elements_2, filter_number);
			mywavCluster.integrate_waveforms(coh_msecs, int_wav, size_wav);
			wavint_data.set_init_range(start_window_delay);
			wavint_data.set_waveform(int_wav, size_wav);
			wavint_data.get_range_waveform(range_wav, size_wav);
			//Write results
			fwrite(&sod, sizeof(int), 1, fpo_data);
			fwrite(&prn, sizeof(int), 1, fpo_data);
			for(i=0; i<size_wav; i++){
				fwrite((char*)&range_wav[i], sizeof(double), 1, fpo_data);
				fwrite((char*)&int_wav[i], sizeof(double), 1, fpo_data);
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


