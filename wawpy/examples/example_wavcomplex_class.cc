#include "waveform_structures.h"
#include "wavpy_global_variables.h"
#include "waveform_power.h"
#include "waveform_complex.h"
#include "waveform_pylib_config.h"

using namespace std;

//============================ PROGRAM =====================================
int main (int argc, char* argv[]) {
	int size_read, i;
	FILE *fpi;
	struct wav_raw waveform_raw;
	Waveform_complex_cluster mywavClusterL1;
	mywavClusterL1.initialize(1000, 64);
	Waveform_complex_cluster mywavClusterL2;
	mywavClusterL2.initialize(1000, 64);
	double phasorI[1000];
	double phasorQ[1000];
	double phasesL1[1000];
	signed char valid_phasor[1000];
	signed char wavcomplexI[64];
	signed char wavcomplexQ[64];
	int weeksow;
	short millisecond;
	char link_updw;
	int week;
	int sow;
	int link;
	int numcorr;
	int sow_ref;
	bool first_wav = true;
	string data_path = LOCAL_DATA_PATH;
	string data_file = data_path + "Greenland_sample.WAV";
	if((fpi=fopen(data_file.c_str(),"r"))==NULL){
		printf("Can't open input file \n");
		return 0;
	}
	size_read = 1;
	while (size_read == 1){
		size_read = fread(&waveform_raw,sizeof(waveform_raw),1,fpi);
		if(size_read == 1){
			weeksow        = waveform_raw.weeksow;
			millisecond    = waveform_raw.millisecond;
			link_updw      = waveform_raw.link_updw;
			numcorr = (waveform_raw.status_numcorr)&0x0F;
			//status  = (waveform_raw.status_numcorr)&0xF0;
			// Decode weeksow
			week = weeksow;
			// added 1024 in the following sentence
			week = 1024+((week >> 20)&0xFFFF);
			sow  = weeksow&(0x0000FFFFF);
			// Decode link_updw
			link = (link_updw >> 5)&0x03;
			//updw = (link_updw >> 4)&0x01;
			if((numcorr == 7)&&first_wav){
				sow_ref = sow;
				first_wav = false;
			}
			if(link == 1){
				if(sow>sow_ref){
					mywavClusterL1.store_phasor_wavs(33);
					mywavClusterL1.get_phasor( phasorI, 1000, phasorQ, 1000, valid_phasor, 1000 );
					for(i=0;i<1000;i++){
						if(valid_phasor[i] == 1){
							phasesL1[i] = atan2(phasorQ[i], phasorI[i]);
						}else{
							phasesL1[i] = 0.0;
						}
					}
					mywavClusterL1.initialize(1000, 64);
				}
				for(i=0;i<64;i++){
					wavcomplexI[i] = waveform_raw.data[i*2+1];
					wavcomplexQ[i] = waveform_raw.data[i*2];
				}
				mywavClusterL1.add_waveform_GOLD( wavcomplexI, 64, wavcomplexQ, 64, int(millisecond) );
			}
			if(link == 2){
				if(sow>sow_ref){
					mywavClusterL2.store_phasor_wavs(33);
					mywavClusterL2.counterrot_phasor(phasesL1, 1000, valid_phasor, 1000);
					mywavClusterL2.get_phasor( phasorI, 1000, phasorQ, 1000, valid_phasor, 1000 );
					for(i=0;i<1000;i++){
						if(valid_phasor[i] == 1){
							printf("PHASE %d %d %f\n",sow_ref,i,(atan2(phasorQ[i], phasorI[i])));
						}
					}
					sow_ref = sow;
					mywavClusterL2.initialize(1000, 64);
				}
				for(i=0;i<64;i++){
					wavcomplexI[i] = waveform_raw.data[i*2+1];
					wavcomplexQ[i] = waveform_raw.data[i*2];
				}
				mywavClusterL2.add_waveform_GOLD( wavcomplexI, 64, wavcomplexQ, 64, int(millisecond) );
			}
		}
	}
	fclose(fpi);
	return 0;
}
