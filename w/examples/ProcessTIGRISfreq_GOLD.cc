#include <dirent.h>
#include <wav_netcdf.h>
#include "../recipes_cpp/nr3.h"
#include "../ancillary_functions.h"
#include "../waveform_power.h"
#include "../waveform_complex.h"
#include "../specular_geometry.h"
#include "../reflecting_surface.h"
#include "../rf_front_end.h"
#include "../gnss_composite.h"
#include "../ZavorotnyVoronovichModel.h"


using namespace std;

bool readAntennaPatterns( double angles[187], double gain[187], double correctedPattern[360][360] );
void get_ddm_parameters( float** normalized_ddm, double delta_freq, int freq_samples, int lag_spec, float ddm_parameters[13] );
bool ReadSP3File( const char* namefile, VecDoub &tt, MatDoub &xx, MatDoub &yy, MatDoub &zz, MatDoub &vxx, MatDoub &vyy, MatDoub &vzz, int &nsats, int &weekGPS_ref, char gnss_identifier );
bool Interpol_sat_sp3( double sow, int prn, VecDoub &results, VecDoub &tt, MatDoub &xx, MatDoub &yy, MatDoub &zz, bool satvel, MatDoub &vxx, MatDoub &vyy, MatDoub &vzz );
void ReadUndulationFile( MatDoub &undMat, VecDoub &undLats, VecDoub &undLongs );
double Interpol_und( double latitude, double longitude, MatDoub &undMat, VecDoub &undLats, VecDoub &undLongs );


//============================ PROGRAM =====================================
int main (int argc, char* argv[]){
	int msec_coh, sec_uncoh, total_millisecs, cluster_pos, n_wav, i, ncid, lag_pos_max, total_wavs, disc_wavs;
	int size, size_ext, n_files, wav_sow, index_group, index_link, wav_corr, wav_link, wav_prn;
	struct dirent **namelist;
	int n = 0;
	string dir, strcoh, struncoh, file_input, file_output, sat_file;
	FILE *fpo_data;
	string ext = ".WAV.nc";
	size_ext = ext.size();
	struct wav_raw_nc s;
	struct wav_raw_nc *w=&s;
	short wav_GPSweek, wav_msec;
	signed char I_component[64], Q_component[64];
	float pow_wav[64];
	double absolute_time_sec, mid_sow_REF, max_val;
	double delMax_link1, pow_link1, noise_link1, cohtime_link1, antgain_linkUP, antgain_linkDW;
	double delMax_link2, pow_link2, noise_link2, cohtime_link2, delDer_link2, sigMax_link2, sigDer_link2, nom_delay_link2;
	double delMax_link3, pow_link3, noise_link3, cohtime_link3, delDer_link3, sigMax_link3, sigDer_link3, nom_delay_link3;
	double absolute_time_sec_REF[3], accum_nom_delay[3][2];
	int prn_group[3];
	bool first_wav[3];
	double delay_zero_GOLD = (299792458.0/20000000)*33.0; //GOLD-RTR delay 0 during TIGRIS = Lag 33
	double pi_num = 3.141592653589793238462643383279502884197;
	MatDoub xx, yy, zz, vxx, vyy, vzz;
	VecDoub tt;
	VecDoub resultsT(6);  //[x, y, z, vx, vy, vz]
	int nsats, ref_GPSweek;
	bool correctData, sat_vel;
	double sow_diffWeek;
	MatDoub undMat;
	VecDoub undLats, undLongs;
	ZaVoModel_GNSSR wavmodel;
	double elev_ant = 30.0;
	double azim_ant = 110.0;
	Waveform_power waveform;
	double coast_pattern[360][360];
	double angles[187], pattern[187];
	double vectorE[3], vectorH[3], posR_ECEF[3], posT_ECEF[3], velT_ECEF[3];
	double reflecRHCP[2], reflecLHCP[2];
	int ddm_coh_int = 100;
	int freq_samples = 400;
	int pos_max_ddm[2];
	double delta_freq = 0.1;
	double ddm_max_link2, ddm_max_link3, lagHolo_ratio_link2, lagHolo_ratio_link3, central_LagHolo, sides_LagHolo;
	float **normalized_ddm;
	normalized_ddm = (float **) malloc(64*sizeof(float *));
	for(i=0;i<64;i++)
		normalized_ddm[i] = (float *) malloc(freq_samples*sizeof(float));
	int fft_samples = 512;
	float powerLagHolo[fft_samples];
	float ddm_parameters_link2[13], ddm_parameters_link3[13];
	double epsilon_air[2];
	epsilon_air[0] = 1.0;
	epsilon_air[1] = 0.0;
	////////////////////////////////////////////////////////////////////////////////////////
	for(i=0; i<3; i++){
		accum_nom_delay[i][0] = 0.0;
		accum_nom_delay[i][1] = 0.0;
		first_wav[i] = true;
	}
	dir = ".";
	if(argc==4){
		msec_coh = atoi(argv[1]);
		sec_uncoh = atoi(argv[2]);
		strcoh.assign(argv[1]);
		struncoh.assign(argv[2]);
		sat_file.assign(argv[3]);
	}else{
		cout << "USAGE:  " << argv[0] << "  msec_coherent  sec_uncoherent  sp3_sat_file" << endl;
		return 1;
	}
	total_millisecs = sec_uncoh*1000;
	Waveform_complex_cluster** wavCluster = new Waveform_complex_cluster*[9];
	for(i = 0; i < 9; i++){
		wavCluster[i] = new Waveform_complex_cluster;
		wavCluster[i]->initialize(total_millisecs, 64);
	}
	signed char valid_phasor[total_millisecs];
	float phases[total_millisecs], navBit[total_millisecs];
	file_output = "GOLDfreqdel_" + strcoh + "msec_" + struncoh + "sec.txt";
	sat_vel = ReadSP3File(sat_file.c_str(), tt, xx, yy, zz, vxx, vyy, vzz, nsats, ref_GPSweek, 'G'); 
	if(nsats == 0){
		cout << "ERROR! Empty or wrong GPS.sp3 file: " << sat_file << endl;
		return 1;
	}
	ReadUndulationFile(undMat, undLats, undLongs);
	if(!readAntennaPatterns(angles, pattern, coast_pattern)){
		return 1;
	}
	//Waveform model
	waveform.set_sampling_rate(20000000.0);
	wavmodel.receiver_Up.set_receiver_params(0.0, 10.0, 3.0, 12000000.0, 0);
	wavmodel.receiver_Up.set_antenna_pattern_interp(angles, 187, pattern, 187, -11.0);
	vectorE[0] = -sin(elev_ant*pi_num/180.0);
	vectorE[1] = 0.0;
	vectorE[2] = -cos(elev_ant*pi_num/180.0);
	vectorH[0] = 0.0;
	vectorH[1] = 1.0;
	vectorH[2] = 0.0;
	wavmodel.receiver_Up.set_antenna_orientation_EH(vectorE, vectorH);
	wavmodel.receiver_Down.set_receiver_params(0.0, 200.0, 3.0, 12000000.0, 0);
	wavmodel.receiver_Down.set_antenna_whole_pattern(coast_pattern);
	vectorE[0] = cos((90.0 - elev_ant)*pi_num/180.0);
	vectorE[2] = -sin((90.0 - elev_ant)*pi_num/180.0);
	wavmodel.receiver_Down.set_antenna_orientation_EH(vectorE, vectorH);
	posR_ECEF[0] = -2450.730138;
	posR_ECEF[1] = 5363.012234;
	posR_ECEF[2] = 2423.760669;
	wavmodel.geometry.set_ECEFpos_Receiver(posR_ECEF);
	wavmodel.geometry.set_inertials(0.0, 0.0, azim_ant);
	wavmodel.surface.mss_1 = 0.01;
	wavmodel.surface.mss_2 = 0.01;
	wavmodel.gnss_signal.weight_CA = 1.0;
	wavmodel.gnss_signal.weight_PY = 0.0;
	wavmodel.gnss_signal.weight_M = 0.0;
	wavmodel.sampling_rate = 20000000.0;
	wavmodel.wav_length = 64;
	wavmodel.coherent_integration = double(msec_coh)/1000.0;
	////////////////////////////////////////////////////////////////////////////////////////
	n_files = scandir(dir.c_str(), &namelist, 0, alphasort);
	total_wavs = 0;
	disc_wavs = 0;
	if(n_files < 0) perror("scandir");
	else{
		while(n < n_files){
			file_input = namelist[n]->d_name;
			size = file_input.size();
			if(size>size_ext && int(file_input.find(ext,size-size_ext))==(size-size_ext)){
				cout << "=======================================" << endl;
				cout << "Analyzing " << file_input << endl;
				if(wav_open(file_input.c_str(), NULL, WAV_NC_TYPE, "r", &ncid) != 1) return 1;
				if(first_wav[0] && first_wav[1] && first_wav[2]){
					fpo_data = fopen(file_output.c_str(), "w");
				}else{
					fpo_data = fopen(file_output.c_str(), "a");
				}
				n_wav = 0;
				while(wav_read(ncid, WAV_NC_TYPE, n_wav++, (void*)w, NULL)){
					wav_GPSweek = w->week;
					wav_msec = w->millisecond;   
					wav_sow = w->sow;
					wav_corr = int(w->channel);
					wav_link = int(w->link);
					wav_prn = int(w->prn);
					index_group = int((wav_corr - 1)/3);
					if((wav_link > 0)&&(wav_link < 4)&&(wav_corr > 0)&&(wav_corr < 10)&&(wav_GPSweek > 1750)&&(wav_GPSweek < 1760)&&(wav_msec < 1000)&&(wav_msec > -1)&&(wav_sow < 604800)&&(wav_sow > -1)){
						total_wavs ++;
						for(i=0;i<64;i++){
							I_component[i] = w->I[i];
							Q_component[i] = w->Q[i];
						}
						if(first_wav[index_group]){
							absolute_time_sec_REF[index_group] = double(wav_GPSweek)*604800.0 + double(wav_sow) - double(wav_sow%sec_uncoh);
							prn_group[index_group] = wav_prn;
							first_wav[index_group] = false;
						}
						absolute_time_sec = double(wav_GPSweek)*604800.0 + double(wav_sow) + double(wav_msec)/1000.0;
						cluster_pos = int(round((absolute_time_sec - absolute_time_sec_REF[index_group])*1000.0));
						if(cluster_pos < 0){
							cout << "WARNING! Bad timing in wav " << n_wav << endl;
							disc_wavs ++;
							//return 1;
						}
						if((cluster_pos >= 0)&&(cluster_pos < total_millisecs)&&(wav_prn == prn_group[index_group])){
							wavCluster[index_group*3 + wav_link - 1]->add_waveform_GOLD(I_component, 64, Q_component, 64, cluster_pos);
							if(wav_link != 1){
								accum_nom_delay[index_group][wav_link - 2] = accum_nom_delay[index_group][wav_link - 2] + double(w->nominal_delay);
							}
						}
						if(((cluster_pos == (total_millisecs - 1))&&(wav_link == 3))||(cluster_pos >= total_millisecs)||(wav_prn != prn_group[index_group])){
							if((wavCluster[index_group*3]->num_valid_wavs > 0)&&(wavCluster[index_group*3 + 1]->num_valid_wavs > 0)&&(wavCluster[index_group*3 + 2]->num_valid_wavs > 0)){
								mid_sow_REF = absolute_time_sec_REF[index_group] - floor(absolute_time_sec_REF[index_group]/604800.0)*604800.0 + double(sec_uncoh)/2.0;
								//Geometry Model
								sow_diffWeek = absolute_time_sec_REF[index_group] - ref_GPSweek*604800.0 + double(sec_uncoh)/2.0;
								correctData = Interpol_sat_sp3(sow_diffWeek, wav_prn, resultsT, tt, xx, yy, zz, sat_vel, vxx, vyy, vzz);
								for(i=0; i<3; i++){
									posT_ECEF[i] = resultsT[i];
									velT_ECEF[i] = resultsT[i + 3];
								}
								wavmodel.geometry.set_ECEFpos_Transmitter(posT_ECEF);
								wavmodel.geometry.set_ECEFvel_Transmitter(velT_ECEF);
								wavmodel.geometry.set_Undulation(0.0);
								wavmodel.geometry.compute_specular_point(0);
								wavmodel.geometry.set_Undulation(Interpol_und(wavmodel.geometry.latitudeS, wavmodel.geometry.longitudeS, undMat, undLats, undLongs));
								wavmodel.geometry.compute_specular_point(0);
								wavmodel.surface.compute_Rfresnel_circular((90.0-wavmodel.geometry.elevation)*pi_num/180.0, epsilon_air, reflecRHCP, reflecLHCP);
								antgain_linkUP = wavmodel.receiver_Up.get_angular_gain_dB( (wavmodel.geometry.elevation - elev_ant), (wavmodel.geometry.azimuthT - azim_ant));
								antgain_linkDW = wavmodel.receiver_Down.get_angular_gain_dB( (elev_ant - wavmodel.geometry.elevation), (azim_ant - wavmodel.geometry.azimuthT));
								//Link 1 RHCP-D
								wavCluster[index_group*3]->correct_navigation_bit(33, 1);
								cohtime_link1 = wavCluster[index_group*3]->compute_coherence_time(33, 0);
								max_val = wavCluster[index_group*3]->integrate_waveforms(msec_coh, pow_wav, 64);
								waveform.set_norm_waveform(pow_wav, 64, max_val);
								waveform.compute_delays();
								pow_link1 = waveform.powerMax;
								noise_link1 = waveform.floorNoise;
								delMax_link1 = waveform.positionMax;
								wavCluster[index_group*3]->get_phasor(navBit, total_millisecs, phases, total_millisecs, valid_phasor, total_millisecs);
								for(i=0; i<total_millisecs; i++){
									phases[i] = atan2(0.0, navBit[i]);
									valid_phasor[i] = 1;
								}
								//Link 2 LHCP-R
								nom_delay_link2 = accum_nom_delay[index_group][0]/double(wavCluster[index_group*3 + 1]->num_valid_wavs);
								wavCluster[index_group*3 + 1]->counterrot_waveforms(phases, total_millisecs, valid_phasor, total_millisecs);
								cohtime_link2 = wavCluster[index_group*3 + 1]->compute_coherence_time(33, 0);
								max_val = wavCluster[index_group*3 + 1]->integrate_waveforms_remdir(msec_coh, total_millisecs, pow_wav, 64);
								waveform.set_norm_waveform(pow_wav, 64, max_val);
								waveform.compute_delays();
								delMax_link2 = waveform.positionMax;
								sigMax_link2 = waveform.sigma_posMax;
								pow_link2 = waveform.powerMax;
								delDer_link2 = waveform.positionDer;
								sigDer_link2 = waveform.sigma_posDer;
								noise_link2 = waveform.floorNoise;
								//
								lag_pos_max = int(delMax_link2/(C_LIGHT/wavmodel.sampling_rate));
								wavCluster[index_group*3 + 1]->compute_LagHologram( lag_pos_max, powerLagHolo, fft_samples);
								sides_LagHolo = 0.0;
								central_LagHolo = 0.0;
								for(i=0;i<fft_samples/4;i++){
									sides_LagHolo = sides_LagHolo + powerLagHolo[i] + powerLagHolo[i + (fft_samples*3/4)];
									central_LagHolo = central_LagHolo + powerLagHolo[i + (fft_samples/4)] + powerLagHolo[i + (fft_samples/2)]; 
								}
								if(sides_LagHolo > 0.0){
									lagHolo_ratio_link2 = central_LagHolo/sides_LagHolo;
								}else{
									lagHolo_ratio_link2 = 0.0;
								}
								//
								ddm_max_link2 = wavCluster[index_group*3 + 1]->compute_whole_DDM_remdir(ddm_coh_int, total_millisecs, delta_freq, normalized_ddm, 64, freq_samples, pos_max_ddm);
								get_ddm_parameters(normalized_ddm, delta_freq, freq_samples, lag_pos_max, ddm_parameters_link2);
								//Link 3 RHCP-R
								nom_delay_link3 = accum_nom_delay[index_group][1]/double(wavCluster[index_group*3 + 2]->num_valid_wavs);
								wavCluster[index_group*3 + 2]->counterrot_waveforms(phases, total_millisecs, valid_phasor, total_millisecs);
								cohtime_link3 = wavCluster[index_group*3 + 2]->compute_coherence_time(33, 0);
								max_val = wavCluster[index_group*3 + 2]->integrate_waveforms_remdir(msec_coh, total_millisecs, pow_wav, 64);
								waveform.set_norm_waveform(pow_wav, 64, max_val);
								waveform.compute_delays();
								delMax_link3 = waveform.positionMax;
								sigMax_link3 = waveform.sigma_posMax;
								pow_link3 = waveform.powerMax;
								delDer_link3 = waveform.positionDer;
								sigDer_link3 = waveform.sigma_posDer;
								noise_link3 = waveform.floorNoise;
								//
								lag_pos_max = int(delMax_link3/(C_LIGHT/wavmodel.sampling_rate));
								wavCluster[index_group*3 + 2]->compute_LagHologram( lag_pos_max, powerLagHolo, fft_samples);
								sides_LagHolo = 0.0;
								central_LagHolo = 0.0;
								for(i=0;i<fft_samples/4;i++){
									sides_LagHolo = sides_LagHolo + powerLagHolo[i] + powerLagHolo[i + (fft_samples*3/4)];
									central_LagHolo = central_LagHolo + powerLagHolo[i + (fft_samples/4)] + powerLagHolo[i + (fft_samples/2)]; 
								}
								if(sides_LagHolo > 0.0){
									lagHolo_ratio_link3 = central_LagHolo/sides_LagHolo;
								}else{
									lagHolo_ratio_link3 = 0.0;
								}
								//
								ddm_max_link3 = wavCluster[index_group*3 + 2]->compute_whole_DDM_remdir(ddm_coh_int, total_millisecs, delta_freq, normalized_ddm, 64, freq_samples, pos_max_ddm);
								get_ddm_parameters(normalized_ddm, delta_freq, freq_samples, lag_pos_max, ddm_parameters_link3);
								//Print Results
								fprintf(fpo_data, "%f %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %f %f\n",
								mid_sow_REF, //1
								prn_group[index_group], //2
								(delMax_link1 - delay_zero_GOLD), //3
								10.0*log10(pow_link1), //4
								10.0*log10(noise_link1), //5
								cohtime_link1, //6
								antgain_linkUP, //7
								(delMax_link2 - delay_zero_GOLD + nom_delay_link2), //8
								10.0*log10(pow_link2), //9
								10.0*log10(noise_link2), //10
								cohtime_link2, //11
								(delDer_link2 - delay_zero_GOLD + nom_delay_link2), //12
								sigMax_link2, //13
								sigDer_link2, //14
								double(ddm_parameters_link2[0])*ddm_max_link2, //15
								ddm_parameters_link2[1], //16
								ddm_parameters_link2[2], //17
								double(ddm_parameters_link2[3])*ddm_max_link2, //18
								ddm_parameters_link2[4], //19
								ddm_parameters_link2[5], //20
								double(ddm_parameters_link2[6])*ddm_max_link2, //21
								ddm_parameters_link2[7], //22
								ddm_parameters_link2[8], //23
								double(ddm_parameters_link2[9])*ddm_max_link2, //24
								ddm_parameters_link2[10], //25
								ddm_parameters_link2[11], //26
								ddm_parameters_link2[12], //27
								(delMax_link3 - delay_zero_GOLD + nom_delay_link3), //28
								10.0*log10(pow_link3), //29
								10.0*log10(noise_link3), //30
								cohtime_link3, //31
								(delDer_link3 - delay_zero_GOLD + nom_delay_link3), //32
								sigMax_link3, //33
								sigDer_link3, //34
								double(ddm_parameters_link3[0])*ddm_max_link3, //35
								ddm_parameters_link3[1], //36
								ddm_parameters_link3[2], //37
								double(ddm_parameters_link3[3])*ddm_max_link3, //38
								ddm_parameters_link3[4], //39
								ddm_parameters_link3[5], //40
								double(ddm_parameters_link3[6])*ddm_max_link3, //41
								ddm_parameters_link3[7], //42
								ddm_parameters_link3[8], //43
								double(ddm_parameters_link3[9])*ddm_max_link3, //44
								ddm_parameters_link3[10], //45
								ddm_parameters_link3[11], //46
								ddm_parameters_link3[12], //47
								antgain_linkDW, //48
								10.0*log10(reflecLHCP[0]*reflecLHCP[0] + reflecLHCP[1]*reflecLHCP[1]), //49
								10.0*log10(reflecRHCP[0]*reflecRHCP[0] + reflecRHCP[1]*reflecRHCP[1]), //50
								(wavmodel.geometry.geometric_delay*1000.0), //51
								wavmodel.geometry.elevation, //52
								wavmodel.geometry.azimuthT, //53
								wavmodel.geometry.longitudeS, //54
								wavmodel.geometry.latitudeS, //55
								wavCluster[index_group*3]->num_valid_wavs, //56
								wavCluster[index_group*3 + 1]->num_valid_wavs, //57
								wavCluster[index_group*3 + 2]->num_valid_wavs, //58
								lagHolo_ratio_link2, //59
								lagHolo_ratio_link3 //60
								);
							}
							wavCluster[index_group*3]->reset();
							wavCluster[index_group*3 + 1]->reset();
							wavCluster[index_group*3 + 2]->reset();
							absolute_time_sec_REF[index_group] = absolute_time_sec_REF[index_group] + double(sec_uncoh);
							accum_nom_delay[index_group][0] = 0.0;
							accum_nom_delay[index_group][1] = 0.0;
							if((cluster_pos > (total_millisecs - 1))||(wav_prn != prn_group[index_group])){
								absolute_time_sec_REF[index_group] = double(wav_GPSweek)*604800.0 + double(wav_sow) - double(wav_sow%sec_uncoh);
								cluster_pos = int((absolute_time_sec - absolute_time_sec_REF[index_group])*1000.0);
								prn_group[index_group] = wav_prn;
								wavCluster[index_group*3 + index_link]->add_waveform_GOLD(I_component, 64, Q_component, 64, cluster_pos);
								if(wav_link != 1){
									accum_nom_delay[index_group][wav_link - 2] = double(w->nominal_delay);
								}
							}
						}
					}
				}
				fclose(fpo_data);
				if(wav_close(ncid) != 1) return 1;
			}
			free(namelist[n]);
			n ++;
		}
		free(namelist);
	}
	if(first_wav[0] && first_wav[1] && first_wav[2]){
		cout << "ERROR! Any output file generated!" << file_input << endl;
	}else{
		cout << "FINISH! Total valid wavs: " << total_wavs << ". Discarded wavs: " << disc_wavs << " (" << (100.0*float(disc_wavs)/float(total_wavs)) << "%)" << endl; 
	}
	return 0;
}


bool readAntennaPatterns( double angles[187], double gain[187], double correctedPattern[360][360] )
{
	stringstream lineRead;
	string line;
	int row, col;
	string namefile1 = "../ANCILLARY/angles_gain_AntPattern.txt";
	string namefile2 = "../ANCILLARY/coast_antenna_pattern.txt";
	ifstream fpi1(namefile1.c_str());
	ifstream fpi2(namefile2.c_str());
	if(!fpi1){
		cout << "ERROR! Unable to read file: " << namefile1 << endl ;
		return false;
	}
	if(!fpi2){
		cout << "ERROR! Unable to read file: " << namefile2 << endl ;
		return false;
	}
	if(fpi1.is_open()){
		row = 0;
		while(!fpi1.eof()){
			getline(fpi1,line);
			lineRead << line;
			lineRead.clear();
			lineRead.str(line);
			lineRead >> angles[row] >> gain[row];
			row ++;
		}
		fpi1.close();
	}
	if(fpi2.is_open()){
		row = 0;
		while(!fpi2.eof()){
			getline(fpi2,line);
			lineRead << line;
			lineRead.clear();
			lineRead.str(line);
			col = 0;
			while(lineRead >> correctedPattern[row][col]) col++;
			row ++;
		}
		fpi2.close();
	}
	return true;
}


void get_ddm_parameters( float** normalized_ddm, double delta_freq, int freq_samples, int lag_spec, float ddm_parameters[13] )
{
	int i, k, i_left, i_right, ref_lag, i_freq_max;
	int lag_offset[4];
	double sum_left, sum_right;
	float ref_val, coef_left, coef_right, freq_left, freq_right;
	memset(ddm_parameters, 0, sizeof(float)*13);
	lag_offset[0] = 0;
	lag_offset[1] = 10;
	lag_offset[2] = 20;
	lag_offset[3] = -10;
	for(k=0; k<4; k++){
		ref_lag = lag_spec + lag_offset[k];
		if((ref_lag >= 0)&&(ref_lag < 64)){
			for(i=0; i<freq_samples; i++){
				if(normalized_ddm[ref_lag][i] > ddm_parameters[k*3]){
					ddm_parameters[k*3] = normalized_ddm[ref_lag][i];
					ddm_parameters[k*3 + 1] = float(i - freq_samples/2)*float(delta_freq);
					i_freq_max = i;
				}
			}
			ref_val = ddm_parameters[k*3]/2.0;
			if((normalized_ddm[ref_lag][0] < ref_val)&&(normalized_ddm[ref_lag][63] < ref_val)){
				i = 0;
				while((normalized_ddm[ref_lag][i] < ref_val)&&(i < i_freq_max)) i++;
				coef_left = (normalized_ddm[ref_lag][i] - ref_val)/(normalized_ddm[ref_lag][i] - normalized_ddm[ref_lag][i-1]);
				coef_right = (ref_val - normalized_ddm[ref_lag][i-1])/(normalized_ddm[ref_lag][i] - normalized_ddm[ref_lag][i-1]);
				freq_left = (float(i -1 - freq_samples/2)*coef_left + float(i - freq_samples/2)*coef_right)*float(delta_freq);
				i = freq_samples;
				while((normalized_ddm[ref_lag][i] < ref_val)&&(i > i_freq_max)) i--;
				coef_left = (ref_val - normalized_ddm[ref_lag][i+1])/(normalized_ddm[ref_lag][i] - normalized_ddm[ref_lag][i+1]);
				coef_right = (normalized_ddm[ref_lag][i] - ref_val)/(normalized_ddm[ref_lag][i] - normalized_ddm[ref_lag][i+1]);
				freq_right = (float(i - freq_samples/2)*coef_left + float(i +1 - freq_samples/2)*coef_right)*float(delta_freq);
				ddm_parameters[k*3 + 2] = freq_right - freq_left;
			}
		}
	}
	i_left = lag_spec - 1;
	i_right = lag_spec + 1;
	sum_left = 0.0;
	sum_right = 0.0;
	while((i_left >= 0)&&(i_right <= 63)){
		for(i=0; i<freq_samples; i++){
			sum_left = sum_left + double(normalized_ddm[i_left][i]);
			sum_right = sum_right + double(normalized_ddm[i_right][i]);
		}
		i_left --;
		i_right ++;
	}
	if(sum_left > 0.0){
		ddm_parameters[12] = float(sum_right/sum_left);
	}
	return;
}

