#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <wav_netcdf.h>
#include "../waveform_power.h"

using namespace std;

bool check_wav_to_retrack(struct wav_int_nc *wav_ref, struct wav_int_nc *wav_new, int sec_interv);
void assign_retracked_waveform(struct wav_int_nc *wav_out, Waveform_power retracked_wav);

//============================ PROGRAM =====================================
int main (int argc, char* argv[]){
	int i, ncidr, ncidw, i_chn, n_wav_in, n_wav_out, secs_interval, min_wavs, used_wavs;
	double diff_delay;
	double max_diff_delay = 90.0; //6 GOLD-RTR lags
	double single_wavPow[64];
	bool first_wav[10];
	int nwavs[10];
	string file_input, file_output;
	struct wav_int_nc wav_ref[10];
	struct wav_int_nc s;
	struct wav_int_nc *w=&s;
	//===========================================================================
	for(i=0; i<10; i++){
		first_wav[i] = true;
	}
	Waveform_power** wavsRetracked = new Waveform_power*[10];
	for(i = 0; i < 10; i++){
		wavsRetracked[i] = new Waveform_power;
		wavsRetracked[i]->set_sampling_rate(20000000.0); //GOLD-RTR
	}
	if(argc==4){
		file_input.assign(argv[1]);
		file_output.assign(argv[2]);
		secs_interval = atoi(argv[3]);
	}else{
		cout << "USAGE:  " << argv[0] << "  input_file  output_file  secs_interval" << endl;
		return 1;
	}
	cout << "|===========> ProcessLongIntegration " << endl;
	min_wavs = secs_interval/2;
	n_wav_in = 0;
	n_wav_out = 0;
	used_wavs = 0;
	if(wav_open(file_input.c_str(), NULL, WAV_INT_NC_TYPE, "r", &ncidr) != 1) return 1;
	if(wav_open(file_output.c_str(), NULL, WAV_INT_NC_TYPE, "w", &ncidw) != 1) return 1;
	while(wav_read(ncidr, WAV_INT_NC_TYPE, n_wav_in++, (void*)w, NULL)){
		i_chn = int(w->channel) - 1;
		if((i_chn > -1)&&(i_chn < 10)&&(w->sow < 604800)&&(w->sow > -1)){
			if(first_wav[i_chn]){
				wav_ref[i_chn] = *w;
				nwavs[i_chn] = 0;
				first_wav[i_chn] = false;
			}
			if(check_wav_to_retrack(&wav_ref[i_chn], w, secs_interval)){
				//Differential delay to be retracked
				//1st correction: GeomCorr_REF - GeomCorr_NEW
				//2nd correction: SWD_NEW - SWD_REF
				diff_delay = (wav_ref[i_chn].geometric_delay - wav_ref[i_chn].eccentricity_delay - wav_ref[i_chn].atmospheric_delay) - (w->geometric_delay - w->eccentricity_delay - w->atmospheric_delay) + w->start_window_delay - wav_ref[i_chn].start_window_delay;
				//Comparison with max_diff_delay due to FFT-shift periodic effect
				if(fabs(diff_delay) < max_diff_delay){
					nwavs[i_chn] ++;
					used_wavs ++;
					for(i=0; i<64; i++){
						single_wavPow[i] = double(w->waveform[i]);
					}
					wavsRetracked[i_chn]->add_waveform_retracking(single_wavPow, 64, diff_delay, 1.0/nwavs[i_chn]);
				}
			}else{
				if(nwavs[i_chn] > min_wavs){
					//Assign retracked waveform and computed delays
					assign_retracked_waveform(&wav_ref[i_chn], *wavsRetracked[i_chn]);
					//Write integrated waveform
					wav_write(ncidw, WAV_INT_NC_TYPE, n_wav_out++, (void*)&wav_ref[i_chn], NULL);
				}else{
					used_wavs = used_wavs - nwavs[i_chn];
				}
				wav_ref[i_chn] = *w;
				nwavs[i_chn] = 1;
				for(i=0; i<64; i++){
					single_wavPow[i] = double(w->waveform[i]);
				}
				wavsRetracked[i_chn]->add_waveform_retracking(single_wavPow, 64, 0.0, 1.0);
			}
		}
	}
	for(i=0; i<10; i++){
		if(nwavs[i_chn] > min_wavs){
			//Assign retracked waveform and computed delays
			assign_retracked_waveform(&wav_ref[i_chn], *wavsRetracked[i_chn]);
			//Write integrated waveform
			wav_write(ncidw, WAV_INT_NC_TYPE, n_wav_out++, (void*)&wav_ref[i_chn], NULL);
		}
	}
	if(wav_close(ncidr) != 1) return 1;
	if(wav_close(ncidw) != 1) return 1;
	cout << "|=> END! Waveforms discarded: " << (n_wav_in - used_wavs) << " of " << n_wav_in << " (" << (100.0*(1.0 - float(used_wavs)/float(n_wav_in))) << "%)" << endl;
	return 0;
}

bool check_wav_to_retrack(wav_int_nc *wav_ref, wav_int_nc *wav_new, int sec_interv){
	int sow_diff;
	sow_diff = wav_new->sow - wav_ref->sow + int(wav_new->week - wav_ref->week)*604800;
	if((sow_diff >= sec_interv)||(sow_diff < 0)){
		return false;
	}
	if((wav_new->d_freq)!=(wav_new->d_freq)){
		return false;
	}
	if((wav_new->prn)!=(wav_new->prn)){
		return false;
	}
	if((wav_new->link)!=(wav_new->link)){
		return false;
	}
	return true;
}

void assign_retracked_waveform(wav_int_nc *wav_out, Waveform_power retracked_wav){
	int i;
	double wavPow[64];
	retracked_wav.get_waveform(wavPow, 64);
	for(i=0; i<64; i++){
		wav_out->waveform[i] = float(wavPow[i]);
	}
	retracked_wav.compute_delays();
	wav_out->max_wav = float(retracked_wav.powerMax);
	wav_out->max_wav_delay = retracked_wav.positionMax + wav_out->start_window_delay;
	wav_out->sig_max_wav_delay = float(retracked_wav.sigma_posMax);
	wav_out->max_der_delay = retracked_wav.positionDer + wav_out->start_window_delay;
	wav_out->sig_max_der_delay = float(retracked_wav.sigma_posDer);
	return;
}

