#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <dirent.h>
#include <wav_netcdf.h>
#include "../waveform_power.h"
#include "../waveform_complex.h"

using namespace std;


//============================ PROGRAM =====================================
int main (int argc, char* argv[]){
	int msec_coh, sec_uncoh, total_millisecs, cluster_pos, n_wav, i, ncid;
	int size, size_ext, n_files, wav_sow, prev_sow, frozen_sow;
	struct dirent **namelist;
	int n = 0;
	string dir, strcoh, struncoh, file_input, file_output;
	FILE *fpo_data;
	string ext = ".WAV-PIR.nc";
	size_ext = ext.size();
	struct wav_pir s;
	struct wav_pir *w=&s;
	short wav_GPSweek, wav_msec, prev_msec;
	short wavXiYi[320], wavXqYq[320], wavXiYq[320], wavXqYi[320];
	float pow_wav[320];
	float coh_time;
	double absolute_time_sec, absolute_time_sec_REF, mid_sow_REF;
	bool first_wav = true;
	double posMax, sigmaMax, powMax, posDer, sigmaDer, powDer, maxwav;
	double delay_zero_PIR = (299792458.0/80000000)*80.0; //PIR delay 0 during TIGRIS = Lag 80
	dir = ".";
	if(argc==3){
		msec_coh = atoi(argv[1]);
		sec_uncoh = atoi(argv[2]);
		strcoh.assign(argv[1]);
		struncoh.assign(argv[2]);
	}else{
		cout << "USAGE:  " << argv[0] << "  msec_coherent  sec_uncoherent" << endl;
		return 1;
	}
	total_millisecs = sec_uncoh*1000;
	file_output = "PIRdelays_" + strcoh + "msec_" + struncoh + "sec.txt";
	Waveform_complex_cluster wavCluster;
	wavCluster.initialize(total_millisecs, 320);
	Waveform_power waveform;
	waveform.set_sampling_rate(80000000.0);
	prev_sow = 0;
	prev_msec = 0;
	frozen_sow = -1;
	n_files = scandir(dir.c_str(), &namelist, 0, alphasort);
	if(n_files < 0) perror("scandir");
	else{
		while(n < n_files){
			file_input = namelist[n]->d_name;
			size = file_input.size();
			if(size>size_ext && int(file_input.find(ext,size-size_ext))==(size-size_ext)){
				cout << "=======================================" << endl;
				cout << "Analyzing " << file_input << endl;
				if(wav_open(file_input.c_str(), NULL, WAV_PIR_TYPE, "r", &ncid) != 1) return 1;
				if(first_wav){
					fpo_data = fopen (file_output.c_str(),"w");
				}else{
					fpo_data = fopen (file_output.c_str(),"a");
				}
				n_wav = 0;
				while(wav_read(ncid, WAV_PIR_TYPE, n_wav++, (void*)w, NULL)){
					wav_GPSweek = w->week;
					wav_msec = w->millisecond;   
					wav_sow = w->sow;
					if((wav_GPSweek > 1750)&&(wav_GPSweek < 1760)&&(wav_msec < 1000)&&(wav_msec > -1)&&(wav_sow < 604800)&&(wav_sow > -1)){
						for(i=0;i<320;i++){
							wavXiYi[i] = w->XiYi[i];
							wavXqYq[i] = w->XqYq[i];
							wavXiYq[i] = w->XiYq[i];
							wavXqYi[i] = w->XqYi[i];
						}
						if(first_wav){
							absolute_time_sec_REF = double(wav_GPSweek)*604800.0 + double(wav_sow) - double(wav_sow%sec_uncoh);
							first_wav = false;
						}
						if((wav_sow == prev_sow)&&(wav_msec < prev_msec)){ //Due to PIR error during TIGRIS
							wav_sow = prev_sow + 1;
							frozen_sow = prev_sow;
							if(wav_sow >= 604800){
								wav_sow = wav_sow - 604800;
								wav_GPSweek++;
							}
						}
						if(wav_sow == frozen_sow){ //Due to PIR error during TIGRIS
							if(wav_msec < prev_msec){
								wav_sow = prev_sow + 1;
								if(wav_sow >= 604800){
									wav_sow = wav_sow - 604800;
									wav_GPSweek++;
								}
							}else{
								wav_sow = prev_sow;
							}
						}
						absolute_time_sec = double(wav_GPSweek)*604800.0 + double(wav_sow) + double(wav_msec)/1000.0;
						cluster_pos = int(round((absolute_time_sec - absolute_time_sec_REF)*1000.0));
						if(cluster_pos < 0){
							cout << "ERROR! Bad timing in wav " << n_wav << endl;
							return 1;
						}
						if(cluster_pos < total_millisecs){
							if(((wav_sow%60)<30) || ((wav_sow%60)>35)){ //During seconds 30 to 35 there is noise autocorr
								wavCluster.add_waveform_PIR(wavXiYi, 320, wavXqYq, 320, wavXiYq, 320, wavXqYi, 320, cluster_pos);
							}
						}
						if(cluster_pos >= (total_millisecs - 1)){
							if(wavCluster.num_valid_wavs > 0){
								maxwav = wavCluster.integrate_waveforms_remdir(msec_coh,  total_millisecs, pow_wav, 320);
								waveform.set_norm_waveform(pow_wav, 320, maxwav);
								waveform.compute_delays();
								coh_time = wavCluster.compute_coherence_time(int(waveform.positionMax/3.75), 0);
								mid_sow_REF = absolute_time_sec_REF - floor(absolute_time_sec_REF/604800.0)*604800.0 + double(sec_uncoh)/2.0;
								printf("Printing sow=%f\n", mid_sow_REF);
								fprintf(fpo_data, "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", mid_sow_REF, (waveform.positionMax - delay_zero_PIR), waveform.sigma_posMax, 10*log10(waveform.powerMax), 10*log10(pow_wav[0]*maxwav), coh_time, (waveform.positionDer - delay_zero_PIR), waveform.sigma_posDer);
							}
							wavCluster.reset();
							if(cluster_pos == (total_millisecs - 1)){
								absolute_time_sec_REF = absolute_time_sec_REF + double(sec_uncoh);
							}else{
								absolute_time_sec_REF = double(wav_GPSweek)*604800.0 + double(wav_sow) - double(wav_sow%sec_uncoh);
								cluster_pos = int((absolute_time_sec - absolute_time_sec_REF)*1000.0);
								if(((wav_sow%60)<30) || ((wav_sow%60)>35)){ //During seconds 30 to 35 there is noise autocorr
									wavCluster.add_waveform_PIR(wavXiYi, 320, wavXqYq, 320, wavXiYq, 320, wavXqYi, 320, cluster_pos);
								}
							}
						}
						prev_sow = wav_sow;
						prev_msec = wav_msec;
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
	if(first_wav){
		cout << "ERROR! Any output file generated!" << file_input << endl;
	}
	return 0;
}