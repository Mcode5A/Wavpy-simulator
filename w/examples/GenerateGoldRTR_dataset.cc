#include "ancillary_functions.h"
#include "wavpy_global_variables.h"
#include "waveform_power.h"
#include "waveform_complex.h"
#include "specular_geometry.h"

struct wav_header {
  //General info
  int prn;
  int gps_week;
  int gps_sow;
  int num_wavs_L1;
  int num_wavs_L2;
  int num_wavs_L3;
  //WAV info
  int range_samples;
  char dummy[4];
  double range_start;
  double range_resolution;
  //Data measurements at zero doppler (waveform)
  double max_wav_delay_L1;
  double pnr_L1_dB;
  double max_der_delay_L2;
  double max_wav_delay_L2;
  double pnr_L2_dB;
  double max_der_delay_L3;
  double max_wav_delay_L3;
  double pnr_L3_dB;
  //Delay models
  double geometric_delay;
  double eccentricity_delay_L2;
  double eccentricity_delay_L3;
  double atmospheric_delay;
  //Geometry from ancillary data
  double height_rcv;
  double azimuth_txr;
  double elevation_txr;
  double longitude_spec;
  double latitude_spec;
  double undulation_egm96_spec;
  double rcv_position[3];
  double txr_position[3];
  double spec_position[3];
  double rcv_velocity[3];
  double txr_velocity[3];
  double attitude[3];
};

// waveform complete structure definition
struct waveform_in {
  int weeksow;                  // 0
  short millisecond;            // 4
  char status_numcorr;          // 6
  char link_updw;               // 7
  char prn;                     // 8
  char max_pos;                 // 9
  char amplitude;               // 10
  char phase;                   // 11
  int range_model;              // 12
  int doppler_avion;            // 16
  int sampling_freq_int;        // 20
  short sampling_freq_frac;     // 24
  short d_freq;                 // 26
  short d_tao;                  // 28
  short sin_elevation;          // 30
  char data[128];               // 32
};

void georef( int sow, int msec, int week, int ref_week, int prn, bool load_info, wav_header *wavinfo, Specular_geometry geom, int nlines_sat, bool sat_vel, double* tt, double** xx, double** yy, double** zz, double** vxx, double** vyy, double** vzz, int nlines_posR, double** week_sow_XYZ, int nlines_inert, double** week_sow_RollPitchYaw, double* undMat, double* undLats, double* undLongs, double &range_mod_L2, double &range_mod_L3 );

using namespace std;

//============================ PROGRAM =====================================
int main (int argc, char* argv[]){
	int size_wav = 64;
	int i, j, prn, coh_msecs, incoh_secs, week, sow, msec, msec_ref_sow, isec_ref_sow, i_rtrack;
	string goldrtr_file_list, sp3file, posRfile, inertfile, outfile, file_input_name, line;
	stringstream lineRead;
	FILE *fpi, *fpo_data;
	signed char *wav_i, *wav_q;
	int num_pow_wavs[3], num_wavs_msec[3], num_wavs_isec[3];
	double *undMat, *undLats, *undLongs, *retracking_reord, *retracking_L1, *retracking_L2, *retracking_L3, *wav;
	double **xx, **yy, **zz, **vxx, **vyy, **vzz, *tt, **week_sow_XYZ, **week_sow_RollPitchYaw;
	int size_read, status, numcorr, link;
	int nsats, nlines_sat, sow_diff, ref_GPSweek, nlines_posR, nlines_inert;
	bool sat_vel, msec_ref_req, isec_ref_req, isec_swd_ref_req;
	bool *valid_rtrack_L1, *valid_rtrack_L2, *valid_rtrack_L3; 
	double sowmilli, prev_sowmilli, swd; 
	double range_mod_L2, range_mod_L3, ref_msec_range_mod_L2, ref_msec_range_mod_L3, ref_isec_range_mod_L2, ref_isec_range_mod_L3;
	double sampling_rate = 20000000.0;
	Waveform_complex_cluster** mywavCluster = new Waveform_complex_cluster*[3];
	for(i=0; i<3; i++){
		mywavCluster[i] = new Waveform_complex_cluster;
		mywavCluster[i]->initialize(1000, size_wav);
	}
	wav = (double *) calloc(size_wav, sizeof(double));
	Waveform_power** wavs_pow = new Waveform_power*[3];
	for(i=0; i<3; i++){
		wavs_pow[i] = new Waveform_power;
		wavs_pow[i]->set_sampling_rate(sampling_rate);
		wavs_pow[i]->set_min_resolution_fft_interp(0.05);
		wavs_pow[i]->set_waveform(wav, size_wav);
		//wavs_pow[i]->set_fit_length(5.0);
		num_pow_wavs[i] = 0;
		num_wavs_isec[i] = 0;
		num_wavs_msec[i] = 0;
	}
	wavs_pow[0]->set_init_range(-32.0*C_LIGHT/sampling_rate);
	Waveform_power wav_pow;
	wav_pow.set_sampling_rate(sampling_rate);
	wav_pow.set_init_range(-32.0*C_LIGHT/sampling_rate);
	Specular_geometry mygeom;
	wav_header wav_dump_info;
	waveform_in wav_raw;
	////////////////////////////////////////////////////////////////////////////////////////
	if(argc==9){
		coh_msecs = atoi(argv[1]);
		incoh_secs = atoi(argv[2]);
		goldrtr_file_list.assign(argv[3]);
		sp3file.assign(argv[4]);
		posRfile.assign(argv[5]);
		inertfile.assign(argv[6]);
		prn = atoi(argv[7]);
		outfile.assign(argv[8]);
	}else{
		cout << "USAGE:  " << argv[0] << " coh_msecs  incoh_secs  goldrtr_file_list  sp3file  posRfile  inertfile  PRN  outfile" << endl;
		return 1;
	}
	//Check input vales
	if(coh_msecs != 1){
		cout << "ERROR! Wrong value for coh_msecs: " << coh_msecs << endl;
		return 1;
	}
	if(incoh_secs <= 0){
		cout << "ERROR! Wrong value for incoh_secs: " << incoh_secs << endl;
		return 1;
	}
	retracking_L1 = (double *) calloc(1000, sizeof(double));
	retracking_L2 = (double *) calloc(1000, sizeof(double));
	retracking_L3 = (double *) calloc(1000, sizeof(double));
	valid_rtrack_L1 = (bool *) calloc(1000, sizeof(bool));
	valid_rtrack_L2 = (bool *) calloc(1000, sizeof(bool));
	valid_rtrack_L3 = (bool *) calloc(1000, sizeof(bool));
	prev_sowmilli = -1.0;
	msec_ref_req = true;
	isec_ref_req = true;
	isec_swd_ref_req = true;
	msec_ref_sow = 8*86400;
	isec_ref_sow = 8*86400;
	wav_i = (signed char *) malloc(size_wav*sizeof(signed char));
	wav_q = (signed char *) malloc(size_wav*sizeof(signed char));
	//Undulation file
	undMat = (double*) malloc (681*1440*sizeof(double));
	undLats = (double*) malloc (681*sizeof(double));
	undLongs = (double*) malloc (1440*sizeof(double));
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
	sat_vel = ReadSP3File(sp3file.c_str(), tt, xx, yy, zz, vxx, vyy, vzz, nlines_sat, sow_diff, ref_GPSweek, 'G');
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
	wav_dump_info.prn = prn;
	wav_dump_info.range_samples = size_wav;
	wav_dump_info.range_resolution = C_LIGHT/sampling_rate;
	////////////////////////////////////////////////////////////////////////////////////////
	//START PROCESS
	if((fpo_data = fopen(outfile.c_str(), "wb")) == NULL){
		printf("ERROR! outfile can't be opened.\n");
		return 1;
	}
	ifstream filelist(goldrtr_file_list.c_str());
	if(!filelist){
		cout << "ERROR! Unable to read GOLD-RTR-files list file: " << goldrtr_file_list << endl;
		return 1;
	}
	if(filelist.is_open()){
		getline(filelist,line);
		while(!filelist.eof()){
			//lineRead << line;
			lineRead.clear();
			lineRead.str(line);
			lineRead >> file_input_name;
			cout << "=======================================" << endl;
			cout << "Analyzing " << file_input_name << endl;
			fpi = fopen(file_input_name.c_str(),"r");
			if((fpi = fopen(file_input_name.c_str(),"r"))== NULL){
				printf("Cannot open input file \n");
				size_read = 0;
			}else{
				size_read = 1;
			}
			while(size_read == 1){
				size_read = fread(&wav_raw, sizeof(wav_raw), 1, fpi);
				if(size_read == 1) {
					status = wav_raw.status_numcorr&0xF0;
					numcorr = wav_raw.status_numcorr&0x0F;
					if((status == 0x00)&&(wav_raw.prn == prn)&&(numcorr > 0)&&(numcorr < 10)){
						week = 1024 + ((wav_raw.weeksow >> 20)&0x00000FFF);
						sow  = wav_raw.weeksow&0x000FFFFF;
						msec = int(wav_raw.millisecond);
						sowmilli = double(sow) + double(msec)/1000.0;
//==============================================================================================================================================
						if(sow > msec_ref_sow){ //1sec integration of clusters
							//Link 1
							if(num_wavs_msec[0] > 0){
								retracking_reord = (double *) calloc(1000, sizeof(double));
								i_rtrack = 0;
								for(i=0; i<1000; i++){
									if(valid_rtrack_L1[i]){
										retracking_reord[i_rtrack] = retracking_L1[i];
										i_rtrack ++;
									}
								}
								mywavCluster[0]->integrate_waveforms_retracking(coh_msecs, sampling_rate, retracking_reord, 1000, wav, size_wav);
								wav_pow.set_waveform(wav, size_wav);
								wav_pow.compute_delays();
								free(retracking_reord);
								num_pow_wavs[0] ++;
								wavs_pow[0]->add_waveform_retracking(wav, size_wav, -wav_pow.positionMax, 1.0/double(num_pow_wavs[0]), false);
								mywavCluster[0]->initialize(1000, size_wav);
							}
							//Link 2
							if((num_wavs_msec[1] > 0)&&(num_wavs_msec[0] > 0)){
								retracking_reord = (double *) calloc(1000, sizeof(double));
								i_rtrack = 0;
								for(i=0; i<1000; i++){
									if(valid_rtrack_L2[i]){
										if(valid_rtrack_L1[i]){
											retracking_reord[i_rtrack] = retracking_L1[i] + retracking_L2[i];
										}else{
											retracking_reord[i_rtrack] = retracking_L2[i];
										}
										i_rtrack ++;
									}
								}
								mywavCluster[1]->integrate_waveforms_retracking(coh_msecs, sampling_rate, retracking_reord, 1000, wav, size_wav);
								free(retracking_reord);
								num_pow_wavs[1] ++;
								wavs_pow[1]->add_waveform_retracking(wav, size_wav, ref_isec_range_mod_L2 - range_mod_L2 + swd - wav_dump_info.range_start - wav_pow.positionMax, 1.0/double(num_pow_wavs[1]), false);
								mywavCluster[1]->initialize(1000, size_wav);
							}
							//Link 3
							if((num_wavs_msec[2] > 0)&&(num_wavs_msec[0] > 0)){
								retracking_reord = (double *) calloc(1000, sizeof(double));
								i_rtrack = 0;
								for(i=0; i<1000; i++){
									if(valid_rtrack_L3[i]){
										if(valid_rtrack_L1[i]){
											retracking_reord[i_rtrack] = retracking_L1[i] + retracking_L3[i];
										}else{
											retracking_reord[i_rtrack] = retracking_L3[i];
										}
										i_rtrack ++;
									}
								}
								mywavCluster[2]->integrate_waveforms_retracking(coh_msecs, sampling_rate, retracking_reord, 1000, wav, size_wav);
								free(retracking_reord);
								num_pow_wavs[2] ++;
								wavs_pow[2]->add_waveform_retracking(wav, size_wav, ref_isec_range_mod_L3 - range_mod_L3 + swd - wav_dump_info.range_start - wav_pow.positionMax, 1.0/double(num_pow_wavs[2]), false);
								mywavCluster[2]->initialize(1000, size_wav);
							}
							for(i=0; i<1000; i++){
								valid_rtrack_L1[i] = false;
								valid_rtrack_L2[i] = false;
								valid_rtrack_L3[i] = false;
							}
							for(i=0; i<3; i++){
								num_wavs_isec[i] = num_wavs_isec[i] + num_wavs_msec[i];
								num_wavs_msec[i] = 0;
							}
							msec_ref_req = true;
						}
//==============================================================================================================================================
						if(sow >= (isec_ref_sow + incoh_secs)){ //dump of waveforms
							if((num_pow_wavs[0] > 0)&&(num_pow_wavs[1] > 0)&&(num_pow_wavs[2] > 0)){
								//Link 1
								wav_dump_info.num_wavs_L1 = num_wavs_isec[0];
								wavs_pow[0]->compute_delays();
								wav_dump_info.max_wav_delay_L1 = wavs_pow[0]->positionMax;
								wav_dump_info.pnr_L1_dB = 0.0;
								if(wavs_pow[0]->floorNoise > 0.0){
									wav_dump_info.pnr_L1_dB = 10.0*log10(wavs_pow[0]->powerMax/wavs_pow[0]->floorNoise);
								}
								//Link 2
								wav_dump_info.num_wavs_L2 = num_wavs_isec[1];
								wavs_pow[1]->compute_delays();
								wav_dump_info.max_der_delay_L2 = wavs_pow[1]->positionDer;
								wav_dump_info.max_wav_delay_L2 = wavs_pow[1]->positionMax;
								wav_dump_info.pnr_L2_dB = 0.0;
								if(wavs_pow[1]->floorNoise > 0.0){
									wav_dump_info.pnr_L2_dB = 10.0*log10(wavs_pow[1]->power_posDer/wavs_pow[1]->floorNoise);
								}
								//Link 3
								wav_dump_info.num_wavs_L3 = num_wavs_isec[2];
								wavs_pow[2]->compute_delays();
								wav_dump_info.max_der_delay_L3 = wavs_pow[2]->positionDer;
								wav_dump_info.max_wav_delay_L3 = wavs_pow[2]->positionMax;
								wav_dump_info.pnr_L3_dB = 0.0;
								if(wavs_pow[2]->floorNoise > 0.0){
									wav_dump_info.pnr_L3_dB = 10.0*log10(wavs_pow[2]->power_posDer/wavs_pow[2]->floorNoise);
								}
								//WRITE wav_dump_info AND wavs L2/L3 IN BINARY FORMAT
								fwrite(&wav_dump_info, sizeof(wav_header), 1, fpo_data);
								for(i=1; i<3; i++){
									wavs_pow[i]->get_waveform(wav, size_wav);
									for(j=0; j<size_wav; j++){
										fwrite((char*)&wav[j], sizeof(double), 1, fpo_data);
									}
								}
							}
							free(wav);
							wav = (double *) calloc(size_wav, sizeof(double));
							for(i=0; i<3; i++){
								wavs_pow[i]->set_waveform(wav, size_wav);
								num_pow_wavs[i] = 0;
								num_wavs_isec[i] = 0;
							}
							isec_ref_req = true;
							isec_swd_ref_req = true;
						}
//==============================================================================================================================================
						if(sowmilli > prev_sowmilli){ //Georef
							if(isec_ref_req){
								georef(sow - sow%incoh_secs, 0, week, ref_GPSweek, prn, true, &wav_dump_info, mygeom, nlines_sat, sat_vel, tt, xx, yy, zz, vxx, vyy, vzz, nlines_posR, week_sow_XYZ, nlines_inert, week_sow_RollPitchYaw, undMat, undLats, undLongs, range_mod_L2, range_mod_L3);
								ref_isec_range_mod_L2 = range_mod_L2;
								ref_isec_range_mod_L3 = range_mod_L3;
								isec_ref_sow = sow - sow%incoh_secs;
								isec_ref_req = false;
							}
							if(msec_ref_req){
								if(sow%incoh_secs != 0){
									georef(sow, 0, week, ref_GPSweek, prn, false, &wav_dump_info, mygeom, nlines_sat, sat_vel, tt, xx, yy, zz, vxx, vyy, vzz, nlines_posR, week_sow_XYZ, nlines_inert, week_sow_RollPitchYaw, undMat, undLats, undLongs, range_mod_L2, range_mod_L3);
								}
								ref_msec_range_mod_L2 = range_mod_L2;
								ref_msec_range_mod_L3 = range_mod_L3;
								msec_ref_sow = sow;
								msec_ref_req = false;
							}
							if(msec != 0){
								georef(sow, msec, week, ref_GPSweek, prn, false, &wav_dump_info, mygeom, nlines_sat, sat_vel, tt, xx, yy, zz, vxx, vyy, vzz, nlines_posR, week_sow_XYZ, nlines_inert, week_sow_RollPitchYaw, undMat, undLats, undLongs, range_mod_L2, range_mod_L3);
							}
						}
//==============================================================================================================================================
						//Read complex waveform
						link = int((wav_raw.link_updw >> 5)&0x03);
						for(i=0; i<size_wav; i++){
							wav_i[i] = wav_raw.data[i*2 + 1];
							wav_q[i] = wav_raw.data[i*2];
						}
						mywavCluster[link - 1]->add_waveform_GOLD(wav_i, size_wav, wav_q, size_wav, msec);
						num_wavs_msec[link - 1] ++;
						//Retracking values (msec-level) and swd
						if(link == 1){
							for(i=0; i<size_wav; i++){
								wav[i] = double(wav_i[i])*double(wav_i[i]) + double(wav_q[i])*double(wav_q[i]);
							}
							wav_pow.set_waveform(wav, size_wav);
							wav_pow.compute_delays();
							retracking_L1[msec] = wav_pow.positionMax;
							valid_rtrack_L1[msec] = true;
						}else{
							swd = double(wav_raw.range_model) - 32.0*C_LIGHT/sampling_rate;
							if(isec_swd_ref_req){
								wav_dump_info.range_start = swd;
								wavs_pow[1]->set_init_range(wav_dump_info.range_start);
								wavs_pow[2]->set_init_range(wav_dump_info.range_start);
								isec_swd_ref_req = false;
							}
							if(link == 2){
								retracking_L2[msec] = ref_msec_range_mod_L2 - range_mod_L2;
								valid_rtrack_L2[msec] = true;
							}else{
								retracking_L3[msec] = ref_msec_range_mod_L3 - range_mod_L3;
								valid_rtrack_L3[msec] = true;
							}
						}
						prev_sowmilli = sowmilli;
//printf("END %d %f %d %d %d %d %d %d %d %d %d %d\n", isec_ref_sow, sowmilli, link, num_wavs_msec[0], num_wavs_msec[1], num_wavs_msec[2], num_wavs_isec[0], num_wavs_isec[1], num_wavs_isec[2], num_pow_wavs[0], num_pow_wavs[1], num_pow_wavs[2]);
					}
				}
			}
			cout << "=======================================" << endl;
			fclose(fpi);
			getline(filelist, line);
		}
		filelist.close();
	}
	fclose(fpo_data);
	free(retracking_L1);
	free(retracking_L2);
	free(retracking_L3);
	free(valid_rtrack_L1);
	free(valid_rtrack_L2);
	free(valid_rtrack_L3);
	free(wav);
	free(wav_i);
	free(wav_q);
	free(undMat);
	free(undLats);
	free(undLongs);
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
  
void georef( int sow, int msec, int week, int ref_week, int prn, bool load_info, wav_header *wavinfo, Specular_geometry geom, int nlines_sat, bool sat_vel, double* tt, double** xx, double** yy, double** zz, double** vxx, double** vyy, double** vzz, int nlines_posR, double** week_sow_XYZ, int nlines_inert, double** week_sow_RollPitchYaw, double* undMat, double* undLats, double* undLongs, double &range_mod_L2, double &range_mod_L3 ){
	int i;
	double sow_diffWeek, inertial_delay_L2, inertial_delay_L3, atmospheric_delay;
	double resultsT[6], resultsPos[6], resultsInert[6];
	double posT_ECEF[3], posR_ECEF[3], posS_ECEF[3], posR_LLH[3], inert_RPY[3], velT_ECEF[3], velR_ECEF[3], temp_posR_LLH[3];
	double ant_L1_BF[3], ant_L2_BF[3], ant_L3_BF[3], ant_L1_Local[3], ant_L2_Local[3], ant_L3_Local[3];
	ant_L1_BF[0] = -1.21;
	ant_L1_BF[1] = 0.122;
	ant_L1_BF[2] = -0.01;
	ant_L2_BF[0] = -1.23;
	ant_L2_BF[1] = 0.093;
	ant_L2_BF[2] = 2.093;
	ant_L3_BF[0] = 0.19;
	ant_L3_BF[1] = 1.098;
	ant_L3_BF[2] = 2.351;
	sow_diffWeek = double(sow) + double(msec)/1000.0 + double(week - ref_week)*604800.0;
	Interpol_sat_sp3(sow_diffWeek, prn, nlines_sat, resultsT, tt, xx, yy, zz, sat_vel, vxx, vyy, vzz);
	Interpol_ECEFpos_file(week, double(sow) + double(msec)/1000.0, nlines_posR, resultsPos, week_sow_XYZ);
	Interpol_ECEFpos_file(week, double(sow) + double(msec)/1000.0, nlines_inert, resultsInert, week_sow_RollPitchYaw);
	for(i=0; i<3; i++){
		posT_ECEF[i] = resultsT[i];
		velT_ECEF[i] = resultsT[i + 3];
		posR_LLH[i] = resultsPos[i];
		inert_RPY[i] = resultsInert[i];
	}
	geom.set_ECEFpos_Transmitter(posT_ECEF);
	geom.set_ECEFvel_Transmitter(velT_ECEF);
	geom.set_LongLatHeight_Receiver(posR_LLH);
	geom.set_inertials(inert_RPY[0], inert_RPY[1], inert_RPY[2]);
	geom.compute_specular_point(0);
	geom.get_ECEFpos_Specular(posS_ECEF);
	geom.rotate_vector_BF_to_local(ant_L1_BF, ant_L1_Local);
	geom.rotate_vector_BF_to_local(ant_L2_BF, ant_L2_Local);
	geom.rotate_vector_BF_to_local(ant_L3_BF, ant_L3_Local);
	inertial_delay_L2 = ant_L1_Local[1]*(-cos(geom.elevation*PI_NUM/180.0)) + ant_L1_Local[2]*(-sin(geom.elevation*PI_NUM/180.0)) + ant_L2_Local[1]*cos(geom.elevation*PI_NUM/180.0) + ant_L2_Local[2]*(-sin(geom.elevation*PI_NUM/180.0));
	inertial_delay_L3 = ant_L1_Local[1]*(-cos(geom.elevation*PI_NUM/180.0)) + ant_L1_Local[2]*(-sin(geom.elevation*PI_NUM/180.0)) + ant_L3_Local[1]*cos(geom.elevation*PI_NUM/180.0) + ant_L3_Local[2]*(-sin(geom.elevation*PI_NUM/180.0));
	atmospheric_delay = (4.6/sin(geom.elevation*PI_NUM/180.0))*(1.0 - exp(-posR_LLH[2]/8.621));
	range_mod_L2 = geom.geometric_delay*1000.0 - inertial_delay_L2 + atmospheric_delay;
	range_mod_L3 = geom.geometric_delay*1000.0 - inertial_delay_L3 + atmospheric_delay; 
	if(load_info){
		//Compute receiver's velocity
		Interpol_ECEFpos_file(week, double(sow) + 0.5, nlines_posR, resultsPos, week_sow_XYZ);
		for(i=0; i<3; i++){
			temp_posR_LLH[i] = resultsPos[i];
		}
		geom.set_LongLatHeight_Receiver(temp_posR_LLH);
		geom.get_ECEFpos_Receiver(posR_ECEF);
		for(i=0; i<3; i++){
			velR_ECEF[i] = posR_ECEF[i];
		}
		Interpol_ECEFpos_file(week, double(sow) - 0.5, nlines_posR, resultsPos, week_sow_XYZ);
		for(i=0; i<3; i++){
			temp_posR_LLH[i] = resultsPos[i];
		}
		geom.set_LongLatHeight_Receiver(temp_posR_LLH);
		geom.get_ECEFpos_Receiver(posR_ECEF);
		for(i=0; i<3; i++){
			velR_ECEF[i] = velR_ECEF[i] - posR_ECEF[i];
		}
		geom.set_ECEFvel_Receiver(velR_ECEF);
		//============================
		geom.set_LongLatHeight_Receiver(posR_LLH);
		wavinfo->gps_week = week;
		wavinfo->gps_sow = sow;
		wavinfo->geometric_delay = geom.geometric_delay*1000.0;
		wavinfo->eccentricity_delay_L2 = inertial_delay_L2;
		wavinfo->eccentricity_delay_L3 = inertial_delay_L3;
		wavinfo->atmospheric_delay = atmospheric_delay;
		wavinfo->height_rcv = posR_LLH[2]*1000.0;
		wavinfo->azimuth_txr = geom.azimuthT;
		wavinfo->elevation_txr = geom.elevation;
		wavinfo->longitude_spec = geom.longitudeS;
		wavinfo->latitude_spec = geom.latitudeS;
		wavinfo->undulation_egm96_spec = 1000.0*Interpol_und(geom.latitudeS, geom.longitudeS, undMat, undLats, undLongs);
		for(i=0; i<3; i++){
			wavinfo->attitude[i] = inert_RPY[i];
		}
		geom.get_ECEFpos_Receiver(wavinfo->rcv_position);
		geom.get_ECEFvel_Receiver(wavinfo->rcv_velocity);
		geom.get_ECEFpos_Transmitter(wavinfo->txr_position);
		geom.get_ECEFvel_Transmitter(wavinfo->txr_velocity);
		geom.get_ECEFpos_Specular(wavinfo->spec_position);
	}
	return;
}
  
