#include "ancillary_functions.h"
#include "wavpy_global_variables.h"
#include "waveform_power.h"
#include "waveform_complex.h"
#include "specular_geometry.h"
#include "rf_front_end.h"

struct ddm_header {
  //General info
  char gnss_info[4];
  int prn;
  int gps_week;
  int gps_sow;
  int msec_coh;
  int sec_incoh;
  //DDM info
  int range_samples;
  int doppler_samples;
  double range_start;
  double range_resolution;
  double doppler_start;
  double doppler_resolution;
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
	double posreceiver_UPBF[3], posreceiver_DWBF[3];
	posreceiver_UPBF[0] = -1.21; 
	posreceiver_UPBF[1] = 0.122; 
	posreceiver_UPBF[2] = -0.01;
	posreceiver_DWBF[0] = -1.23; 
	posreceiver_DWBF[1] = 0.093; 
	posreceiver_DWBF[2] = 2.093;
	double *retracking_msec;
	stringstream lineRead;
	string line;
	int sod;
	string spirfile;
	double posreceiver_UPLocal[3], posreceiver_DWLocal[3];
	double posT_ECEF[3], posR_ECEF[3], posS_ECEF[3], posR_LLH[3], inert_RPY[3], velT_ECEF[3], velR_ECEF[3], temp_posR_LLH[3]; 
	double phasesUP[8], phasesDW[8], phase_delaysUP[8], phase_delaysDW[8];
	double undu;
	int week, msec;
	double inertial_delay, atmospheric_delay, ref_range;
	double start_window_delay, sow;
	bool first_int = true;
	double apriori_scattdel;
	int filter_number;
	string char_GNSS;
	//NEW VARIABLES
	int coh_msecs, incoh_secs, doppler_samples, retrack_samples;
	double delta_doppler;
	string freq_band;
	double **ddm;
	double ddm_freq_slice[size_wav];
	ddm_header ddm_dump_info;
	int accum_ddms;
	//Beaming optimization
	bool apply_beam_optimization = false;
	double theta_phi_UP[2], theta_phi_DW[2];
	////////////////////////////////////////////////////////////////////////////////////////
	if(argc==13){
		coh_msecs = atoi(argv[1]);
		incoh_secs = atoi(argv[2]);
		doppler_samples = atoi(argv[3]);
		delta_doppler = atof(argv[4]);
		spir_file_list.assign(argv[5]);
		sp3file.assign(argv[6]);
		posRfile.assign(argv[7]);
		inertfile.assign(argv[8]);
		char_GNSS.assign(argv[9]);
		freq_band.assign(argv[10]);
		prn = atoi(argv[11]);
		outfile.assign(argv[12]);
	}else{
		cout << "USAGE:  " << argv[0] << " coh_msecs  incoh_secs  doppler_samples  delta_doppler  spir_file_list  sp3file  posRfile  inertfile  char_GNSS  freq_band  PRN  outfile" << endl;
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
	if(doppler_samples <= 1){
		cout << "ERROR! Wrong value for doppler_samples: " << doppler_samples << endl;
		return 1;
	}
	if(delta_doppler < 0.0){
		cout << "ERROR! Wrong value for delta_doppler: " << delta_doppler << endl;
		return 1;
	}
	if(coh_msecs == 1){
		retrack_samples = 999;
	}else{
		retrack_samples = 1000/coh_msecs;
	}
	retracking_msec = (double *) malloc(retrack_samples*sizeof(double));
	ddm = (double **) malloc(size_wav*sizeof(double *));
	for(i=0; i<size_wav; i++){
		ddm[i] = (double *) malloc(doppler_samples*sizeof(double));
	}
	Waveform_power** accum_ddm = new Waveform_power*[doppler_samples];
	for(i=0; i<doppler_samples; i++){
		accum_ddm[i] = new Waveform_power;
		accum_ddm[i]->set_sampling_rate(sampling_rate);
		accum_ddm[i]->set_min_resolution_fft_interp(0.05);
		accum_ddm[i]->set_fit_length(5.0);
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
		if(char(char_GNSS[0]) == 'E'){
			ddm_dump_info.gnss_info[0] = 'G';
			ddm_dump_info.gnss_info[1] = 'A';
			ddm_dump_info.gnss_info[2] = 'L';
			ddm_dump_info.gnss_info[3] = '5';
		}else{
			ddm_dump_info.gnss_info[0] = 'G';
			ddm_dump_info.gnss_info[1] = 'P';
			ddm_dump_info.gnss_info[2] = 'S';
			ddm_dump_info.gnss_info[3] = '5';
		}
		filter_number = 2;
		receiver_UP.set_frequency(FREQ_GPS_L5);
		receiver_DW.set_frequency(FREQ_GPS_L5);
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
			ddm_dump_info.gnss_info[0] = 'G';
			ddm_dump_info.gnss_info[1] = 'A';
			ddm_dump_info.gnss_info[2] = 'L';
			ddm_dump_info.gnss_info[3] = '1';
			apriori_scattdel = 25.0;
			filter_number = 1;
			accum_ddm[doppler_samples/2]->set_normtail_length(100.0);
		}else{
			ddm_dump_info.gnss_info[0] = 'G';
			ddm_dump_info.gnss_info[1] = 'P';
			ddm_dump_info.gnss_info[2] = 'S';
			ddm_dump_info.gnss_info[3] = '1';
			apriori_scattdel = 10.0;
			filter_number = 0;
			accum_ddm[doppler_samples/2]->set_normtail_length(50.0);
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
	ddm_dump_info.prn = prn;
	ddm_dump_info.msec_coh = coh_msecs;
	ddm_dump_info.sec_incoh = incoh_secs;
	ddm_dump_info.range_samples = size_wav;
	ddm_dump_info.doppler_samples = doppler_samples;
	ddm_dump_info.range_resolution = C_LIGHT/sampling_rate;
	ddm_dump_info.doppler_start = double(-doppler_samples/2)*delta_doppler;
	ddm_dump_info.doppler_resolution = delta_doppler;
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
					//wav_data_accum_10msec.compute_delays();
					accum_ddm[doppler_samples/2]->compute_delays_wlimits((ddm_dump_info.nominal_delay - 16.5*2.0*sin(ddm_dump_info.elevation_txr*PI_NUM/180.0)), 10.0, apriori_scattdel);
					ddm_dump_info.max_der_delay = accum_ddm[doppler_samples/2]->positionDer;
					ddm_dump_info.sigma_max_der_delay = accum_ddm[doppler_samples/2]->sigma_posDer;
					ddm_dump_info.max_der_pow = accum_ddm[doppler_samples/2]->power_posDer;
					ddm_dump_info.max_der_pow_der = accum_ddm[doppler_samples/2]->powerDer_posDer;
					ddm_dump_info.max_wav_delay = accum_ddm[doppler_samples/2]->positionMax;
					ddm_dump_info.sigma_max_wav_delay = accum_ddm[doppler_samples/2]->sigma_posMax;
					ddm_dump_info.max_wav_pow = accum_ddm[doppler_samples/2]->powerMax;
					ddm_dump_info.noise_floor = accum_ddm[doppler_samples/2]->floorNoise;
					//WRITE ddm_dump_info AND ddm IN BINARY FORMAT
					fwrite(&ddm_dump_info, sizeof(ddm_header), 1, fpo_data);
					for(i=0; i<doppler_samples; i++){
						accum_ddm[i]->get_waveform(ddm_freq_slice, size_wav);
						for(j=0; j<size_wav; j++){
							fwrite((char*)&ddm_freq_slice[j], sizeof(double), 1, fpo_data);
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
					mygeom.rotate_vector_BF_to_local(posreceiver_UPBF, posreceiver_UPLocal);
					mygeom.rotate_vector_BF_to_local(posreceiver_DWBF, posreceiver_DWLocal);
					//inertial_delay = posreceiver_UPLocal[1]*(-cos(mygeom.elevation*PI_NUM/180.0)) + posreceiver_UPLocal[2]*(-sin(mygeom.elevation*PI_NUM/180.0)) + posreceiver_DWLocal[1]*cos(mygeom.elevation*PI_NUM/180.0) + posreceiver_DWLocal[2]*(-sin(mygeom.elevation*PI_NUM/180.0));
                    inertial_delay = (posreceiver_DWLocal[1] - posreceiver_UPLocal[1])*(cos(mygeom.elevation*PI_NUM/180.0)) + (posreceiver_DWLocal[2] - posreceiver_UPLocal[2])*(-sin(mygeom.elevation*PI_NUM/180.0)); //Projection over reflected signal path
                    //inertial_delay = (posreceiver_DWLocal[1] - posreceiver_UPLocal[1])*(-cos(mygeom.elevation*PI_NUM/180.0)) + (posreceiver_DWLocal[2] - posreceiver_UPLocal[2])*(-sin(mygeom.elevation*PI_NUM/180.0)); //Projection over direct signal path
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
							phasesUP[i] = (phase_delaysUP[i]*180.0/PI_NUM) - phases_offset[i];
							phasesDW[i] = (phase_delaysDW[i]*180.0/PI_NUM) - phases_offset[i + 8];
						}
						undu = Interpol_und(mygeom.latitudeS, mygeom.longitudeS, undMat, undLats, undLongs);
						if(sod%incoh_secs == 0){
							ddm_dump_info.gps_week = week;
							ddm_dump_info.gps_sow = sow;
							ddm_dump_info.nominal_delay = ref_range;
							ddm_dump_info.geometric_delay = mygeom.geometric_delay*1000.0;
							ddm_dump_info.eccentricity_delay = inertial_delay;
							ddm_dump_info.atmospheric_delay = atmospheric_delay;
							ddm_dump_info.height_rcv = posR_LLH[2]*1000.0;
							ddm_dump_info.azimuth_txr = mygeom.azimuthT;
							ddm_dump_info.elevation_txr = mygeom.elevation;
							ddm_dump_info.longitude_spec = mygeom.longitudeS;
							ddm_dump_info.latitude_spec = mygeom.latitudeS;
							ddm_dump_info.undulation_egm96_spec = undu*1000.0;
							for(i=0; i<3; i++){
								ddm_dump_info.attitude[i] = inert_RPY[i];
							}
							mygeom.get_ECEFpos_Receiver(ddm_dump_info.rcv_position);
							mygeom.get_ECEFvel_Receiver(ddm_dump_info.rcv_velocity);
							mygeom.get_ECEFpos_Transmitter(ddm_dump_info.txr_position);
							mygeom.get_ECEFvel_Transmitter(ddm_dump_info.txr_velocity);
							mygeom.get_ECEFpos_Specular(ddm_dump_info.spec_position);
						}
					}
					if(msec%coh_msecs == 0){
						retracking_msec[msec/coh_msecs] = ref_range - (mygeom.geometric_delay*1000.0 - inertial_delay + atmospheric_delay);
					}
				}
				start_window_delay = mywavCluster.load_ITF_waveforms_SPIR(spirfile.c_str(), ref_range, phasesUP, phasesDW, filter_number);
				mywavCluster.compute_whole_DDM_retracking(coh_msecs, delta_doppler, ddm, size_wav, doppler_samples, sampling_rate, retracking_msec, retrack_samples);
				if(sod%incoh_secs == 0){
					ddm_dump_info.range_start = start_window_delay;
					accum_ddms = 1;
					for(i=0; i<doppler_samples; i++){
						for(j=0; j<size_wav; j++){
							ddm_freq_slice[j] = ddm[j][i];
						}
						accum_ddm[i]->set_init_range(start_window_delay);
						accum_ddm[i]->set_waveform(ddm_freq_slice, size_wav);
					}
				}else{
					accum_ddms ++;
					for(i=0; i<doppler_samples; i++){
						for(j=0; j<size_wav; j++){
							ddm_freq_slice[j] = ddm[j][i];
						}
						accum_ddm[i]->add_waveform_retracking(ddm_freq_slice, size_wav, ddm_dump_info.nominal_delay - ref_range + start_window_delay - ddm_dump_info.range_start, 1.0/double(accum_ddms), false);
					}
				}
			}
			getline(filelist,line);
		}
		filelist.close();
	}
	////////////////////////////////////////////////////////////////////////////////////////
	fclose(fpo_data);
	free(retracking_msec);
	for(i=0; i<size_wav; i++){
		free(ddm[i]);
	}
	free(ddm);
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
