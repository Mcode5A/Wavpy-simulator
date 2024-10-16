#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include "waveform_structures.h"
#include "waveform_power.h"
#include "waveform_pylib_config.h"

using namespace std;

//============================ PROGRAM =====================================
int main (int argc, char* argv[]){
	int size_read, n_wav, i;
	FILE *fpi;
	struct wav_int waveform_int;
	Waveform_power mywaveform;
	float pow_wav[64];
	float pow_wav_ext[320];
	string data_path = LOCAL_DATA_PATH;
	string data_file = data_path + "GEOHALO_sample.WAVINTGEO";
	if((fpi=fopen(data_file.c_str(),"r"))==NULL){
		printf("Can't open input file \n");
		return 0;
	}
	size_read = 1;
	n_wav     = 0;
	while(size_read == 1){
		size_read=fread(&waveform_int,sizeof(waveform_int),1,fpi);
		if(size_read == 1){
			n_wav = n_wav + 1;
			for(i=0;i<64;i++){
				pow_wav[i] = waveform_int.data[i]*waveform_int.data[i];
				pow_wav_ext[i] = waveform_int.data[i]*waveform_int.data[i];
			}
			if(n_wav==2){
				mywaveform.set_float_waveform(pow_wav, 64);
				mywaveform.set_sampling_rate(20000000.0);
				mywaveform.dump_norm_waveform();
				mywaveform.compute_delays();
				mywaveform.dump_delays();
				printf("%f %f %f %f %f %f\n", mywaveform.positionMax, mywaveform.sigma_posMax, mywaveform.powerMax, mywaveform.positionDer, mywaveform.sigma_posDer, mywaveform.power_posDer);
				for(i=64;i<320;i++){
					pow_wav_ext[i] = 0.0;
				}
				mywaveform.set_float_waveform(pow_wav_ext, 320);
				mywaveform.set_sampling_rate(20000000.0);
				mywaveform.compute_delays();
				mywaveform.dump_delays();
				return 0;
			}
		}
	}
	fclose(fpi);
	return 0;
}

