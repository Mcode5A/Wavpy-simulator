#include "ancillary_functions.h"
#include "wavpy_global_variables.h"
#include "waveform_power.h"
#include "waveform_complex.h"
#include "specular_geometry.h"
#include "rf_front_end.h"

struct wav_header {
  //General info
  char gnss_info[4];
  int prn;
  int gps_week;
  int gps_sow;
  int msec_coh;
  int sec_incoh;
  //WAV info
  int range_samples;
  char dummy[4];
  double range_start;
  double range_resolution;
  //Data measurements at zero doppler (waveform)
  double max_der_delay;
  double sigma_max_der_delay;
  double max_der_pow;
  double max_der_pow_der;
  double max_wav_delay;
  double sigma_max_wav_delay;
  double max_wav_pow;
  double noise_floor;
  //Delay models
  double nominal_delay;
  double geometric_delay;
  double eccentricity_delay;
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

void get_theta_phi_beam_optimization( double elev, double azim, double iRPY[3], RF_FrontEnd recUP, RF_FrontEnd recDW, double tp_UP[2], double tp_DW[2] );

using namespace std;

//============================ PROGRAM =====================================
int main (int argc, char* argv[]){
	int size_wav = 512;
	int i, j, prn;
	string spir_file_list, sp3file, posRfile, inertfile, outfile;
	FILE *fpo_data;
	double *undMat, *undLats, *undLongs;
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
	RF_FrontEnd receiver_UP;
	RF_FrontEnd receiver_DW;
	receiver_UP.set_receiver_params(3.0, 20.0, 3.0, 12000000.0, 0);
	receiver_DW.set_receiver_params(3.0, 300.0, 3.0, 12000000.0, 0);
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
	receiver_UP.set_antenna_elements_pos_AF(pos_elems, 8, 2, 0);
	receiver_DW.set_antenna_elements_pos_AF(pos_elems, 8, 2, 0);
	double uparray_BF_E_vec[3], uparray_BF_H_vec[3];
	uparray_BF_E_vec[0] = 1.0;
	uparray_BF_E_vec[1] = 0.0;
	uparray_BF_E_vec[2] = 0.0;
	uparray_BF_H_vec[0] = 0.0;
	uparray_BF_H_vec[1] = -1.0;
	uparray_BF_H_vec[2] = 0.0;
	receiver_UP.set_antenna_orientation_BF_EH(uparray_BF_E_vec, uparray_BF_H_vec);
	double posUParrayBF[3], posDWarrayBF[3];
	posUParrayBF[0] = -1.21; 
	posUParrayBF[1] = 0.122; 
	posUParrayBF[2] = -0.01;
	posDWarrayBF[0] = -1.23; 
	posDWarrayBF[1] = 0.093; 
	posDWarrayBF[2] = 2.093;
	double *retracking_msec, *retracking_direct_msec;
	stringstream lineRead;
	string line;
	int sod;
	string spirfile;
	double posUParrayLocal[3], posDWarrayLocal[3];
	double posT_ECEF[3], posR_ECEF[3], posS_ECEF[3], posR_LLH[3], inert_RPY[3], velT_ECEF[3], velR_ECEF[3], temp_posR_LLH[3];
	double phasesBF[16], phase_delaysUP[8], phase_delaysDW[8];
	double uRT[3], uRS[3], uST[3], vR_min_vT[3], min_vT[3];
	double undu;
	int week, msec;
	double inertial_delay, atmospheric_delay, ref_range;
	double start_window_delay, sow, norm_uRT, norm_uRS, norm_uST, doppler_up, doppler_dw;
	bool first_int = true;
	//NEW VARIABLES
	int coh_msecs, incoh_secs, retrack_samples;
	double *wav, *wav_i, *wav_q;
	wav_header wav_dump_info;
	int wavs_updws;
	//Beaming optimization
	bool apply_beam_optimization = true;
	double theta_phi_UP[2], theta_phi_DW[2];
	////////////////////////////////////////////////////////////////////////////////////////
	if(argc==9){
		coh_msecs = atoi(argv[1]);
		incoh_secs = atoi(argv[2]);
		spir_file_list.assign(argv[3]);
		sp3file.assign(argv[4]);
		posRfile.assign(argv[5]);
		inertfile.assign(argv[6]);
		prn = atoi(argv[7]);
		outfile.assign(argv[8]);
	}else{
		cout << "USAGE:  " << argv[0] << " coh_msecs  incoh_secs  spir_file_list  sp3file  posRfile  inertfile  PRN  outfile" << endl;
		return 1;
	}
	//Check input vales
	if((coh_msecs <= 0)||(coh_msecs > 999)){
		cout << "ERROR! Wrong value for coh_msecs: " << coh_msecs << endl;
		return 1;
	}
	if(incoh_secs <= 0){
		cout << "ERROR! Wrong value for incoh_secs: " << incoh_secs << endl;
		return 1;
	}
	if(coh_msecs == 1){
		retrack_samples = 999;
	}else{
		retrack_samples = 1000/coh_msecs;
	}
	retracking_msec = (double *) malloc(retrack_samples*sizeof(double));
	retracking_direct_msec = (double *) malloc(retrack_samples*sizeof(double));
	wav = (double *) malloc(size_wav*sizeof(double));
	wav_i = (double *) malloc(size_wav*sizeof(double));
	wav_q = (double *) malloc(size_wav*sizeof(double));
	Waveform_power wav_pow;
	wav_pow.set_sampling_rate(sampling_rate);
	Waveform_power** wavs_updw = new Waveform_power*[2];
	for(i=0; i<2; i++){
		wavs_updw[i] = new Waveform_power;
		wavs_updw[i]->set_sampling_rate(sampling_rate);
		wavs_updw[i]->set_min_resolution_fft_interp(0.05);
		wavs_updw[i]->set_fit_length(5.0);
	}
	//Read Antenna gain files
	double antGain[32400];
	double total_gain[181][360];
	string antfileUP = "/home/multivac/fabra/CAMPAIGNS/SPIR/FLIGHT_DEC15/ANALYSIS/TEST/mean_gain_UP.txt";
	ifstream fpi_antUP (antfileUP.c_str());
	if (!fpi_antUP){
		cout << "ERROR! Unable to read UP antenna mean pattern file." << endl ;
		return 1;
	}
	if(fpi_antUP.is_open()){
		for(i=0;i<32400;i++){
			getline(fpi_antUP,line);
			lineRead.clear();
			lineRead.str(line);
			lineRead >> antGain[i];
		}
		fpi_antUP.close();
	}
	for(i=0; i<90; i++){
		for(j=0; j<360; j++){
			total_gain[i][j] = antGain[i*360 + j];
		}
	}
	receiver_UP.set_antenna_whole_pattern(total_gain);
	string antfileDW = "/home/multivac/fabra/CAMPAIGNS/SPIR/FLIGHT_DEC15/ANALYSIS/TEST/mean_gain_DW.txt";
	ifstream fpi_antDW (antfileDW.c_str());
	if (!fpi_antDW){
		cout << "ERROR! Unable to read DW antenna mean pattern file." << endl ;
		return 1;
	}
	if(fpi_antDW.is_open()){
		for(i=0;i<32400;i++){
			getline(fpi_antDW,line);
			lineRead.clear();
			lineRead.str(line);
			lineRead >> antGain[i];
		}
		fpi_antDW.close();
	}
	for(i=0; i<90; i++){
		for(j=0; j<360; j++){
			total_gain[i][j] = antGain[i*360 + j];
		}
	}
	receiver_DW.set_antenna_whole_pattern(total_gain);
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
	wav_dump_info.gnss_info[0] = 'G';
	wav_dump_info.gnss_info[1] = 'P';
	wav_dump_info.gnss_info[2] = 'S';
	wav_dump_info.gnss_info[3] = '1';
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
	wav_dump_info.prn = prn;
	wav_dump_info.msec_coh = coh_msecs;
	wav_dump_info.sec_incoh = incoh_secs;
	wav_dump_info.range_samples = size_wav;
	wav_dump_info.range_resolution = C_LIGHT/sampling_rate;
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
			//lineRead << line;
			lineRead.clear();
			lineRead.str(line);
			lineRead >> sod >> spirfile;
			sow = double(sod) + 4.0*86400.0;
			cout << "At SoD:" << sod << " File:" << spirfile << endl;
			week = 1873;
			if(sod%incoh_secs == 0){
				if(first_int){
					first_int = false;
				}else{
					wavs_updw[1]->compute_delays();
					wav_dump_info.max_der_delay = wavs_updw[1]->positionDer;
					wav_dump_info.sigma_max_der_delay = wavs_updw[1]->sigma_posDer;
					wav_dump_info.max_der_pow = wavs_updw[1]->power_posDer;
					wav_dump_info.max_der_pow_der = wavs_updw[1]->powerDer_posDer;
					wav_dump_info.max_wav_delay = wavs_updw[1]->positionMax;
					wav_dump_info.sigma_max_wav_delay = wavs_updw[1]->sigma_posMax;
					wav_dump_info.max_wav_pow = wavs_updw[1]->powerMax;
					wav_dump_info.noise_floor = wavs_updw[1]->floorNoise;
					//WRITE wav_dump_info AND ddm IN BINARY FORMAT
					fwrite(&wav_dump_info, sizeof(wav_header), 1, fpo_data);
					for(i=0; i<2; i++){
						wavs_updw[i]->get_waveform(wav, size_wav);
						for(j=0; j<size_wav; j++){
							fwrite((char*)&wav[j], sizeof(double), 1, fpo_data);
						}
					}
				}
			}
			if(!first_int){
				for(msec=0; msec<999; msec++){
					sow_diffWeek = sow + double(msec)/1000.0 + (week - ref_GPSweek)*604800.0;
					Interpol_sat_sp3(sow_diffWeek, prn, nlines_sat, resultsT, tt, xx, yy, zz, sat_vel, vxx, vyy, vzz);
					Interpol_ECEFpos_file(week, sow + double(msec)/1000.0, nlines_posR, resultsPos, week_sow_XYZ);
					Interpol_ECEFpos_file(week, sow + double(msec)/1000.0, nlines_inert, resultsInert, week_sow_RollPitchYaw);
					for(i=0; i<3; i++){
						posT_ECEF[i] = resultsT[i];
						velT_ECEF[i] = resultsT[i + 3];
						posR_LLH[i] = resultsPos[i];
						inert_RPY[i] = resultsInert[i];
					}
					mygeom.set_ECEFpos_Transmitter(posT_ECEF);
					mygeom.set_ECEFvel_Transmitter(velT_ECEF);
					mygeom.set_LongLatHeight_Receiver(posR_LLH);
					mygeom.set_inertials(inert_RPY[0], inert_RPY[1], inert_RPY[2]);
					mygeom.compute_specular_point(0);
					mygeom.get_ECEFpos_Specular(posS_ECEF);
					mygeom.rotate_vector_BF_to_local(posUParrayBF, posUParrayLocal);
					mygeom.rotate_vector_BF_to_local(posDWarrayBF, posDWarrayLocal);
					inertial_delay = posUParrayLocal[1]*(-cos(mygeom.elevation*PI_NUM/180.0)) + posUParrayLocal[2]*(-sin(mygeom.elevation*PI_NUM/180.0)) + posDWarrayLocal[1]*cos(mygeom.elevation*PI_NUM/180.0) +  posDWarrayLocal[2]*(-sin(mygeom.elevation*PI_NUM/180.0));
					atmospheric_delay = (4.6/sin(mygeom.elevation*PI_NUM/180.0))*(1.0 - exp(-posR_LLH[2]/8.621));
					if(msec == 0){
						//Compute receiver's velocity
						Interpol_ECEFpos_file(week, sow + 0.5, nlines_posR, resultsPos, week_sow_XYZ);
						for(i=0; i<3; i++){
							temp_posR_LLH[i] = resultsPos[i];
						}
						mygeom.set_LongLatHeight_Receiver(temp_posR_LLH);
						mygeom.get_ECEFpos_Receiver(posR_ECEF);
						for(i=0; i<3; i++){
							velR_ECEF[i] = posR_ECEF[i];
						}
						Interpol_ECEFpos_file(week, sow - 0.5, nlines_posR, resultsPos, week_sow_XYZ);
						for(i=0; i<3; i++){
							temp_posR_LLH[i] = resultsPos[i];
						}
						mygeom.set_LongLatHeight_Receiver(temp_posR_LLH);
						mygeom.get_ECEFpos_Receiver(posR_ECEF);
						for(i=0; i<3; i++){
							velR_ECEF[i] = velR_ECEF[i] - posR_ECEF[i];
						}
						mygeom.set_ECEFvel_Receiver(velR_ECEF);
						//============================
						mygeom.set_LongLatHeight_Receiver(posR_LLH);
						ref_range = mygeom.geometric_delay*1000.0 - inertial_delay + atmospheric_delay;
						mygeom.get_ECEFpos_Receiver(posR_ECEF);
						if(apply_beam_optimization){
							get_theta_phi_beam_optimization(mygeom.elevation, mygeom.azimuthT, inert_RPY, receiver_UP, receiver_DW, theta_phi_UP, theta_phi_DW);
							receiver_UP.compute_phase_delays_UPA(theta_phi_UP[0], theta_phi_UP[1]);
							receiver_DW.compute_phase_delays_UPA(theta_phi_DW[0], theta_phi_DW[1]);
						}else{
							receiver_UP.compute_phase_delays_pos_ECEF_RT(inert_RPY, posR_ECEF, posT_ECEF);
							receiver_DW.compute_phase_delays_pos_ECEF_RT(inert_RPY, posR_ECEF, posS_ECEF);
						}
						receiver_UP.get_phase_delays(phase_delaysUP, 8);
						receiver_DW.get_phase_delays(phase_delaysDW, 8);
						for(i=0; i<8; i++){
							phasesBF[i] = (phase_delaysUP[i]*180.0/PI_NUM) - phases_offset[i];
							phasesBF[i + 8] = (phase_delaysDW[i]*180.0/PI_NUM) - phases_offset[i + 8];
						}
						undu = Interpol_und(mygeom.latitudeS, mygeom.longitudeS, undMat, undLats, undLongs);
						if(sod%incoh_secs == 0){
							wav_dump_info.gps_week = week;
							wav_dump_info.gps_sow = sow;
							wav_dump_info.nominal_delay = ref_range;
							wav_dump_info.geometric_delay = mygeom.geometric_delay*1000.0;
							wav_dump_info.eccentricity_delay = inertial_delay;
							wav_dump_info.atmospheric_delay = atmospheric_delay;
							wav_dump_info.height_rcv = posR_LLH[2]*1000.0;
							wav_dump_info.azimuth_txr = mygeom.azimuthT;
							wav_dump_info.elevation_txr = mygeom.elevation;
							wav_dump_info.longitude_spec = mygeom.longitudeS;
							wav_dump_info.latitude_spec = mygeom.latitudeS;
							wav_dump_info.undulation_egm96_spec = undu*1000.0;
							for(i=0; i<3; i++){
								wav_dump_info.attitude[i] = inert_RPY[i];
							}
							mygeom.get_ECEFpos_Receiver(wav_dump_info.rcv_position);
							mygeom.get_ECEFvel_Receiver(wav_dump_info.rcv_velocity);
							mygeom.get_ECEFpos_Transmitter(wav_dump_info.txr_position);
							mygeom.get_ECEFvel_Transmitter(wav_dump_info.txr_velocity);
							mygeom.get_ECEFpos_Specular(wav_dump_info.spec_position);
						}
						//COMPUTE-DOPPLER
						for(i=0; i<3; i++){
							uRT[i] = (posT_ECEF[i] - posR_ECEF[i]);
							uRS[i] = (posS_ECEF[i] - posR_ECEF[i]);
							uST[i] = (posT_ECEF[i] - posS_ECEF[i]);
						}
						norm_uRT = norm3vec(uRT);
						norm_uRS = norm3vec(uRS);
						norm_uST = norm3vec(uST);
						for(i=0; i<3; i++){
							uRT[i] = uRT[i]/norm_uRT;
							uRS[i] = uRS[i]/norm_uRS;
							uST[i] = uST[i]/norm_uST;
						}
						for(i=0; i<3; i++){
							vR_min_vT[i] = velR_ECEF[i] - velT_ECEF[i];
							min_vT[i] = -velT_ECEF[i];
						}
						doppler_up = 1000.0*scalar3Prod(vR_min_vT, uRT)/(C_LIGHT/FREQ_GPS_L1);
						doppler_dw = 1000.0*(scalar3Prod(velR_ECEF, uRS) + scalar3Prod(min_vT, uST))/(C_LIGHT/FREQ_GPS_L1);
						start_window_delay = mywavCluster.load_CR_waveforms_SPIR(spirfile.c_str(), 0.0, doppler_up, 0.0, phasesBF, true, prn);
						wav_pow.set_init_range(start_window_delay);
					}
					if(msec%coh_msecs == 0){
						retracking_msec[msec/coh_msecs] = ref_range - (mygeom.geometric_delay*1000.0 - inertial_delay + atmospheric_delay);
						mywavCluster.get_waveform(wav_i, size_wav, wav_q, size_wav, msec);
						for(i=0; i<size_wav; i++){
							wav[i] = wav_i[i]*wav_i[i] + wav_q[i]*wav_q[i];
						}
						wav_pow.set_waveform(wav, size_wav);
						wav_pow.compute_delays();
						retracking_direct_msec[msec/coh_msecs] = wav_pow.positionMax;
					}
				}
				for(i=0; i<retrack_samples; i++){
					retracking_direct_msec[i] = retracking_direct_msec[500] - retracking_direct_msec[i]; //VALID FOR 1MSEC COH TIME!!!!!
					retracking_msec[i] = retracking_msec[i] + retracking_direct_msec[i];
				}
				//Direct signal
				mywavCluster.integrate_waveforms_retracking(coh_msecs, sampling_rate, retracking_direct_msec, retrack_samples, wav, size_wav);
				wav_pow.set_waveform(wav, size_wav);
				wav_pow.set_init_range(start_window_delay);
				wav_pow.compute_delays();
				if(sod%incoh_secs == 0){
					wavs_updws = 1;
					wavs_updw[0]->set_init_range(start_window_delay);
					wavs_updw[0]->set_waveform(wav, size_wav); //This will not have any effect (just initialization)
					wavs_updw[0]->add_waveform_retracking(wav, size_wav, -wav_pow.positionMax, 1.0, false);
				}else{
					wavs_updws ++;
					wavs_updw[0]->add_waveform_retracking(wav, size_wav, -wav_pow.positionMax, 1.0/double(wavs_updws), false);
				}
				//Reflected signal
				start_window_delay = mywavCluster.load_CR_waveforms_SPIR(spirfile.c_str(), ref_range, doppler_up, doppler_dw - doppler_up, phasesBF, false, prn);
				mywavCluster.integrate_waveforms_retracking(coh_msecs, sampling_rate, retracking_msec, retrack_samples, wav, size_wav);
				if(sod%incoh_secs == 0){
					wav_dump_info.range_start = start_window_delay;
					wavs_updw[1]->set_init_range(start_window_delay);
					wavs_updw[1]->set_waveform(wav, size_wav); //This will not have any effect (just initialization)
					wavs_updw[1]->add_waveform_retracking(wav, size_wav, -wav_pow.positionMax, 1.0, false);
				}else{
					wavs_updw[1]->add_waveform_retracking(wav, size_wav, wav_dump_info.nominal_delay - ref_range + start_window_delay - wav_dump_info.range_start - wav_pow.positionMax, 1.0/double(wavs_updws), false);
				}
			}
			getline(filelist,line);
		}
		filelist.close();
	}
	////////////////////////////////////////////////////////////////////////////////////////
	fclose(fpo_data);
	free(retracking_msec);
	free(retracking_direct_msec);
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


void get_theta_phi_beam_optimization( double elev, double azim, double iRPY[3], RF_FrontEnd recUP, RF_FrontEnd recDW, double tp_UP[2], double tp_DW[2] ){
	int i, elev_iter;
	double vector[3];
	double tp_ref_UP[2], tp_ref_DW[2];
	double delta_elev, theta, phi, gain_diff, gain_diff_ref_UP, gain_diff_ref_DW;
	delta_elev = 1.0;
	elev_iter = min(int(elev/delta_elev), 30);
	gain_diff_ref_UP = 1000.0;
	gain_diff_ref_DW = 1000.0;
	for(i=0; i<elev_iter; i++){
		//UP
		vector[0] = cos(azim*PI_NUM/180.0)*cos((elev - double(i)*delta_elev)*PI_NUM/180.0);
		vector[1] = sin(azim*PI_NUM/180.0)*cos((elev - double(i)*delta_elev)*PI_NUM/180.0);
		vector[2] = -sin((elev - double(i)*delta_elev)*PI_NUM/180.0);
		//NED2BF_antiRotation
		Rot3DaxisZ(vector, -iRPY[2]);
		Rot3DaxisY(vector, -iRPY[1]);
		Rot3DaxisX(vector, -iRPY[0]);
		theta = atan2(sqrt(vector[0]*vector[0] + vector[1]*vector[1]), -vector[2])*180.0/PI_NUM;
		phi = atan2(-vector[1], vector[0])*180.0/PI_NUM;
		recUP.compute_phase_delays_UPA(theta, phi);
		recUP.compute_array_factor();
		if(phi < 0.0){
			phi = phi + 360.0;
		}
		if(i == 0){
			tp_ref_UP[0] = theta;
			tp_ref_UP[1] = phi;
		}
		gain_diff = (recUP.get_PhiTheta_gain_dB(tp_ref_UP[1], tp_ref_UP[0]) - recUP.get_PhiTheta_gain_dB(tp_ref_UP[1], tp_ref_UP[0] + 1.0))*(recUP.get_PhiTheta_gain_dB(tp_ref_UP[1], tp_ref_UP[0]) - recUP.get_PhiTheta_gain_dB(tp_ref_UP[1], tp_ref_UP[0] + 1.0)) + (recUP.get_PhiTheta_gain_dB(tp_ref_UP[1], tp_ref_UP[0]) - recUP.get_PhiTheta_gain_dB(tp_ref_UP[1], tp_ref_UP[0] - 1.0))*(recUP.get_PhiTheta_gain_dB(tp_ref_UP[1], tp_ref_UP[0]) - recUP.get_PhiTheta_gain_dB(tp_ref_UP[1], tp_ref_UP[0] - 1.0));
		if(gain_diff < gain_diff_ref_UP){
			tp_UP[0] = theta;
			tp_UP[1] = phi;
			gain_diff_ref_UP = gain_diff;
		}
		//DW
		vector[0] = cos(azim*PI_NUM/180.0)*cos((elev - double(i)*delta_elev)*PI_NUM/180.0);
		vector[1] = sin(azim*PI_NUM/180.0)*cos((elev - double(i)*delta_elev)*PI_NUM/180.0);
		vector[2] = sin((elev - double(i)*delta_elev)*PI_NUM/180.0);
		//NED2BF_antiRotation
		Rot3DaxisZ(vector, -iRPY[2]);
		Rot3DaxisY(vector, -iRPY[1]);
		Rot3DaxisX(vector, -iRPY[0]);
		theta = atan2(sqrt(vector[0]*vector[0] + vector[1]*vector[1]), vector[2])*180.0/PI_NUM;
		phi = atan2(vector[1], vector[0])*180.0/PI_NUM;
		recDW.compute_phase_delays_UPA(theta, phi);
		recDW.compute_array_factor();
		if(phi < 0.0){
			phi = phi + 360.0;
		}
		if(i == 0){
			tp_ref_DW[0] = theta;
			tp_ref_DW[1] = phi;
		}
		gain_diff = (recDW.get_PhiTheta_gain_dB(tp_ref_DW[1], tp_ref_DW[0]) - recDW.get_PhiTheta_gain_dB(tp_ref_DW[1], tp_ref_DW[0] + 1.0))*(recDW.get_PhiTheta_gain_dB(tp_ref_DW[1], tp_ref_DW[0]) - recDW.get_PhiTheta_gain_dB(tp_ref_DW[1], tp_ref_DW[0] + 1.0)) + (recDW.get_PhiTheta_gain_dB(tp_ref_DW[1], tp_ref_DW[0]) - recDW.get_PhiTheta_gain_dB(tp_ref_DW[1], tp_ref_DW[0] - 1.0))*(recDW.get_PhiTheta_gain_dB(tp_ref_DW[1], tp_ref_DW[0]) - recDW.get_PhiTheta_gain_dB(tp_ref_DW[1], tp_ref_DW[0] - 1.0));
		if(gain_diff < gain_diff_ref_DW){
			tp_DW[0] = theta;
			tp_DW[1] = phi;
			gain_diff_ref_DW = gain_diff;
		}
	}
	return;
}
