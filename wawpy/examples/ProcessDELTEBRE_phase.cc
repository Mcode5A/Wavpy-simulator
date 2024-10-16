#include <dirent.h>
#include <wav_netcdf.h>
#include "../recipes_cpp/nr3.h"
#include "../ancillary_functions.h"
#include "../waveform_power.h"
#include "../waveform_complex.h"

using namespace std;

//============================ PROGRAM =====================================
int main (int argc, char* argv[]){
	int msec_coh, sec_uncoh, total_millisecs, cluster_pos, n_wav, i, ncid, lag_pos_max, total_wavs, disc_wavs;
	int size, size_ext, n_files, wav_sow, index_group, wav_corr, wav_link, wav_prn;
	double pow_link1, noise_link1, delMax_link1, pow_link2, noise_link2, delMax_link2, sigMax_link2, delDer_link2, sigDer_link2, nom_delay_link2;
	struct dirent **namelist;
	int n = 0;
	string dir, file_input, file_output, file_output_phase, prn_str;
	stringstream number;
	FILE *fpo_data;
	FILE *fpo_data_phase;
	string ext = ".nc";
	size_ext = ext.size();
	struct wav_raw_nc s;
	struct wav_raw_nc *w=&s;
	short wav_GPSweek, wav_msec;
	signed char I_component[64], Q_component[64];
	float pow_wav[64];
	double absolute_time_sec, mid_sow_REF, max_val;
	double absolute_time_sec_REF[5], accum_nom_delay[5];
	int prn_group[5];
	bool first_wav[5];
	double pi_num = 3.141592653589793238462643383279502884197;
	Waveform_power waveform;
	////////////////////////////////////////////////////////////////////////////////////////
	for(i=0; i<5; i++){
		accum_nom_delay[i] = 0.0;
		first_wav[i] = true;
	}
	dir = ".";
	msec_coh = 1;
	sec_uncoh = 1;
	total_millisecs = sec_uncoh*1000;
	Waveform_complex_cluster** wavCluster = new Waveform_complex_cluster*[10];
	for(i = 0; i < 10; i++){
		wavCluster[i] = new Waveform_complex_cluster;
		wavCluster[i]->initialize(total_millisecs, 64);
	}
	signed char valid_phasor[total_millisecs];
	float phases[total_millisecs], phasor_I[total_millisecs], phasor_Q[total_millisecs];
	file_output = "sow_prn_pow1_noise1_dmax1_pow2_noise2_dmax2_sigmax2_deriv2_sigderiv2_nomdelay.txt";
	system("rm Phasor_stopL1_PRN*.txt");
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
				if(first_wav[0] && first_wav[1] && first_wav[2] && first_wav[3] && first_wav[4]){
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
					index_group = int((wav_corr - 1)/2);
					if((wav_link > 0)&&(wav_link < 3)&&(wav_corr > 0)&&(wav_corr < 10)&&(wav_msec < 1000)&&(wav_msec > -1)&&(wav_sow < 604800)&&(wav_sow > -1)){
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
							wavCluster[index_group*2 + wav_link - 1]->add_waveform_GOLD(I_component, 64, Q_component, 64, cluster_pos);
							if(wav_link != 1){
								accum_nom_delay[index_group] = accum_nom_delay[index_group] + double(w->nominal_delay);
							}
						}
						if(((cluster_pos == (total_millisecs - 1))&&(wav_link == 2))||(cluster_pos >= total_millisecs)||(wav_prn != prn_group[index_group])){
							if((wavCluster[index_group*2]->num_valid_wavs > 0)&&(wavCluster[index_group*2 + 1]->num_valid_wavs > 0)){
								mid_sow_REF = absolute_time_sec_REF[index_group] - floor(absolute_time_sec_REF[index_group]/604800.0)*604800.0 + double(sec_uncoh)/2.0;
								//Link 1 RHCP-D
								max_val = wavCluster[index_group*2]->integrate_waveforms(1, pow_wav, 64);
								waveform.set_norm_waveform(pow_wav, 64, max_val);
								waveform.compute_delays();
								pow_link1 = waveform.powerMax;
								noise_link1 = waveform.floorNoise;
								delMax_link1 = waveform.positionMax;
								wavCluster[index_group*2]->store_phasor_wavs(20);
								wavCluster[index_group*2]->get_phasor(phasor_I, total_millisecs, phasor_Q, total_millisecs, valid_phasor, total_millisecs);
								for(i=0; i<total_millisecs; i++){
									phases[i] = atan2(phasor_Q[i], phasor_I[i]);
									valid_phasor[i] = 1;
								}
								//Link 2 LHCP-R
								nom_delay_link2 = accum_nom_delay[index_group]/double(wavCluster[index_group*2 + 1]->num_valid_wavs);
								wavCluster[index_group*2 + 1]->counterrot_waveforms(phases, total_millisecs, valid_phasor, total_millisecs);
								max_val = wavCluster[index_group*2 + 1]->integrate_waveforms(msec_coh, pow_wav, 64);
								waveform.set_norm_waveform(pow_wav, 64, max_val);
								waveform.compute_delays();
								delMax_link2 = waveform.positionMax;
								sigMax_link2 = waveform.sigma_posMax;
								pow_link2 = waveform.powerMax;
								delDer_link2 = waveform.positionDer;
								sigDer_link2 = waveform.sigma_posDer;
								noise_link2 = waveform.floorNoise;
								lag_pos_max = int(delMax_link2/15.0);
								if((lag_pos_max<0)||(lag_pos_max>63)){
									lag_pos_max = 32;
								}
								wavCluster[index_group*2 + 1]->store_phasor_wavs(lag_pos_max);
								wavCluster[index_group*2 + 1]->get_phasor(phasor_I, total_millisecs, phasor_Q, total_millisecs, valid_phasor, total_millisecs);
								number << prn_group[index_group];
								number >> prn_str;
								number.clear();
								file_output_phase = "Phasor_stopL1_PRN" + prn_str + ".txt";
								fpo_data_phase = fopen(file_output_phase.c_str(), "a");
								for(i=0; i<total_millisecs; i++){
									if(valid_phasor[i] == 1){
										fprintf(fpo_data_phase, "%f %f %f\n", (mid_sow_REF - (double(sec_uncoh)/2.0) + (double(i)/1000.0)), phasor_I[i], phasor_Q[i]);
									}
								}
								fclose(fpo_data_phase);
								//Print Results
								fprintf(fpo_data, "%f %d %f %f %f %f %f %f %f %f %f %f\n",
								mid_sow_REF, //1
								prn_group[index_group], //2
								pow_link1, //3
								noise_link1, //4
								delMax_link1, //5
								pow_link2, //6
								noise_link2, //7
								delMax_link2, //8
								sigMax_link2, //9
								delDer_link2, //10
								sigDer_link2, //11
								nom_delay_link2 //12
								);
							}
							wavCluster[index_group*2]->reset();
							wavCluster[index_group*2 + 1]->reset();
							absolute_time_sec_REF[index_group] = absolute_time_sec_REF[index_group] + double(sec_uncoh);
							accum_nom_delay[index_group] = 0.0;
							if((cluster_pos > (total_millisecs - 1))||(wav_prn != prn_group[index_group])){
								absolute_time_sec_REF[index_group] = double(wav_GPSweek)*604800.0 + double(wav_sow) - double(wav_sow%sec_uncoh);
								cluster_pos = int((absolute_time_sec - absolute_time_sec_REF[index_group])*1000.0);
								prn_group[index_group] = wav_prn;
								wavCluster[index_group*2 + wav_link - 1]->add_waveform_GOLD(I_component, 64, Q_component, 64, cluster_pos);
								if(wav_link != 1){
									accum_nom_delay[index_group] = double(w->nominal_delay);
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
	if(first_wav[0] && first_wav[1] && first_wav[2] && first_wav[3] && first_wav[4]){
		cout << "ERROR! Any output file generated!" << file_input << endl;
	}else{
		cout << "FINISH! Total valid wavs: " << total_wavs << ". Discarded wavs: " << disc_wavs << " (" << (100.0*float(disc_wavs)/float(total_wavs)) << "%)" << endl; 
	}
	return 0;
}
