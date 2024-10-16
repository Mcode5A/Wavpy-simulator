/**
    \file waveform_complex.cc
    \brief Implementation of the Waveform_complex_cluster class.
*/

#include "wavpy_global_variables.h"
#include "waveform_power.h"
#include "waveform_complex.h"
#include "ancillary_functions.h"
#include <fstream>
//#include "gsl/gsl_fft_complex.h"
//#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <fftw3.h>
//#include <time.h>
#include "waveform_pylib_config.h"

void Waveform_complex_cluster::prepare_arrays( int cluster_in_L, int wav_in_L ){
	int i, j;
	//Simple reset
	if((cluster_length == cluster_in_L)&&(wav_length == wav_in_L)){
		for(i=0; i<cluster_length; i++){
			for(j=0; j<wav_length; j++){
				Icomponents[i][j] = 0;
				Qcomponents[i][j] = 0;
			}
			valid_wavs[i] = false;
			phasorI[i] = 0.0;
			phasorQ[i] = 0.0;
			valid_phasor[i] = false;
		}
		num_valid_wavs = 0;
		num_phasor_iter = 0;
		return;
	}
	//Free memory
	if(cluster_length > 0){
		if((wav_length != wav_in_L)||(cluster_length != cluster_in_L)){
			for(i=0; i<cluster_length; i++){
				free(Icomponents[i]);
				free(Qcomponents[i]);
			}
		}
		if(cluster_length != cluster_in_L){
			free(Icomponents);
			free(Qcomponents);
			free(valid_wavs);
			free(phasorI);
			free(phasorQ);
			free(valid_phasor);
		}
	}
	//Memory allocation and zero-initialization
	if((cluster_length == 0)||(cluster_length != cluster_in_L)){
		Icomponents = (double **) calloc(cluster_in_L, sizeof(double *));
		Qcomponents = (double **) calloc(cluster_in_L, sizeof(double *));
		for(i=0; i<cluster_in_L; i++){
			Icomponents[i] = (double *) calloc(wav_in_L, sizeof(double));
			Qcomponents[i] = (double *) calloc(wav_in_L, sizeof(double));
		}
		valid_wavs = (bool*) calloc(cluster_in_L, sizeof(bool));
		phasorI = (double*) calloc(cluster_in_L, sizeof(double));
		phasorQ = (double*) calloc(cluster_in_L, sizeof(double));
		valid_phasor = (bool*) calloc(cluster_in_L, sizeof(bool));
	}
	if((cluster_length == cluster_in_L)&&(wav_length != wav_in_L)){
		for(i=0; i<cluster_in_L; i++){
			Icomponents[i] = (double *) calloc(wav_in_L, sizeof(double));
			Qcomponents[i] = (double *) calloc(wav_in_L, sizeof(double));
		}
	}
	cluster_length = cluster_in_L;
	wav_length = wav_in_L;
	num_valid_wavs = 0;
	num_phasor_iter = 0;
	return;
}

double Waveform_complex_cluster::compute_ITF_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double phases_BF[16], bool elements_UP[16], bool elements_DW[16], int filter_num ){
	int msec, spir_samples_msec, start_window_sample;
	double sampling_rate, lag_res, start_window_delay, interm_freq;
	int *spirdata;
	float *filter;
	const int spir_bytes = 320012288;
	sampling_rate = 80000000.0;
	lag_res = C_LIGHT/sampling_rate;
	spir_samples_msec = int(sampling_rate)/1000;
	if(sizeof(int) != 4){
		printf("ERROR! Size of int must be equal to 4 bytes.\n");
		return 0.0;
	}
	std::ifstream spirfile(namefile);
	spirfile.seekg(0, std::ios_base::end);
	std::size_t size = spirfile.tellg();
	if(size != spir_bytes){
		std::cout << "ERROR! File <" << namefile << "> is not a complete SPIR file." << std::endl ;
		spirfile.close();
		return 0.0;
	}
	spirfile.seekg(0, std::ios_base::beg);
	spirdata = (int *) malloc(spir_bytes);
	spirfile.read((char*) &spirdata[0], size);
	spirfile.close();
	if(!check_SPIR_blocks(spirdata)){
		free(spirdata);
		return 0.0;
	}
	//Start Process
	start_window_sample = int(peak_delay_estimate/lag_res) - wav_length/2;
	start_window_delay = double(start_window_sample)*lag_res;
	filter = (float *) malloc(spir_samples_msec*sizeof(float));
	switch(filter_num){
		case 1:
			//Galileo at E1
			compute_RF_filter_GalileoE1(filter);
			interm_freq = 19420000.0;
			break;
		case 2:
			//GPS or Galileo at L5
			compute_RF_filter_GPS_L5(filter);
			interm_freq = 19307100.0;
			break;
		default:
			//GPS at L1 without C/A code
			compute_RF_filter(filter);
			interm_freq = 19420000.0;
			break;
	}
	//clock_t t;
	//t = clock();
	for(msec=0; msec<cluster_length; msec++){
		compute_ITFwav_SPIR(spirdata, interm_freq, filter, &Icomponents[msec][0], &Qcomponents[msec][0], start_window_sample, wav_length, phases_BF, elements_UP, elements_DW, msec);
		valid_wavs[msec] = true;
		num_valid_wavs ++;
	}
	//t = clock() - t;
	//printf ("Buffer time: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	free(filter);
	free(spirdata);
	return start_window_delay;
}

double Waveform_complex_cluster::compute_CR_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double doppler_estimate, double delta_doppler, double phases_BF[16], bool elements_CORR[16], int code_ref, bool selected_signals, bool searching, bool assist_doppler ){
	int msec, spir_samples_msec, start_window_sample, i, channel, prn;
	double sampling_rate, lag_res, interm_freq, final_doppler;
	int *spirdata;
	float *filter;
    float *code_CR;
	int code_length = 1023;
	const int spir_bytes = 320012288;
	sampling_rate = 80000000.0;
	interm_freq = 19420000.0; //GPS L1
	lag_res = C_LIGHT/sampling_rate;
	spir_samples_msec = int(sampling_rate)/1000;
	if(sizeof(int) != 4){
		printf("ERROR! Size of int must be equal to 4 bytes.\n");
		return 0.0;
	}
	if(!searching){
        if(((code_ref - 1)<0)||((code_ref - 1)>32)){
            printf("ERROR! Reference code not valid.\n");
            return 0.0;
        }
        code_CR = (float *) malloc(code_length*sizeof(float));
        read_PRN_code(code_ref, code_CR, code_length);
    }
	std::ifstream spirfile(namefile);
	spirfile.seekg(0, std::ios_base::end);
	std::size_t size = spirfile.tellg();
	if(size != spir_bytes){
		std::cout << "ERROR! File <" << namefile << "> is not a complete SPIR file." << std::endl ;
		spirfile.close();
		return 0.0;
	}
	spirfile.seekg(0, std::ios_base::beg);
	spirdata = (int *) malloc(spir_bytes);
	spirfile.read((char*) &spirdata[0], size);
	spirfile.close();
	if(!check_SPIR_blocks(spirdata)){
		free(spirdata);
		return 0.0;
	}
	filter = (float *) malloc(spir_samples_msec*sizeof(float));
///////////////////////////
	double freq;
	double df = sampling_rate/double(spir_samples_msec);
	for(i=0; i<spir_samples_msec; i++){
		if(i > spir_samples_msec/2){
			freq = df*double(i - spir_samples_msec);
		}else{
			freq = df*double(i);
		}
		filter[i] = float(1.0/(1.0 + pow(freq/2000000.0, 40.0)));
	}
///////////////////////////
	//Start Process
	if(searching){
        code_CR = (float *) malloc(code_length*sizeof(float));
        for(channel=0; channel<16; channel++){
            if(elements_CORR[channel]){
                printf("=========> Channel %d \n", channel + 1);
                for(prn=1; prn<33; prn++){
                    printf("PRN %d ", prn);
                    read_PRN_code(prn, code_CR, code_length);
                    get_doppler_initsamp_CR_SPIR(spirdata, interm_freq, code_CR, doppler_estimate, code_length, filter, 512, phases_BF, final_doppler, start_window_sample, selected_signals, elements_CORR, true);
                }
            }
        }
        return 0.0;
    }
    if(assist_doppler){
        get_doppler_initsamp_CR_SPIR(spirdata, interm_freq, code_CR, doppler_estimate, code_length, filter, wav_length, phases_BF, final_doppler, start_window_sample, selected_signals, elements_CORR, false);
    }else{
        final_doppler = doppler_estimate;
    }
	start_window_sample = start_window_sample + int(peak_delay_estimate/lag_res);
	//Start Process
	//clock_t t;
	//t = clock();
	for(msec=0; msec<cluster_length; msec++){
		compute_CRwav_SPIR(spirdata, interm_freq, code_CR, final_doppler + delta_doppler, code_length, filter, &Icomponents[msec][0], &Qcomponents[msec][0], wav_length, phases_BF, elements_CORR, msec, start_window_sample, false);
		valid_wavs[msec] = true;
		num_valid_wavs ++;
	}
	//t = clock() - t;
	//printf ("Buffer time: %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	free(filter);
	free(code_CR);
	free(spirdata);
	return lag_res*double(int(peak_delay_estimate/lag_res) - wav_length/2); //start window delay
}

Waveform_complex_cluster::Waveform_complex_cluster( void ){
	cluster_length = 0;
	wav_length = 0;
	return;
}

void Waveform_complex_cluster::initialize( int in_cluster_length, int wav_in_length ){
	if(in_cluster_length <= 0){
		printf("ERROR! cluster_length must be greater than zero.");
		return;
	}
	if(wav_in_length <= 0){
		printf("ERROR! wav_length must be greater than zero.");
		return;
	}
	prepare_arrays(in_cluster_length, wav_in_length);
	return;
}

void Waveform_complex_cluster::add_waveform( double* Icomp_in, int len_I, double* Qcomp_in, int len_Q, int cluster_pos ){
	int j;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(len_I != len_Q){
		printf("ERROR-1! I and Q components must have equal size (size-I %d, size-Q %d)\n", len_I, len_Q);
		return;
	}
	if(len_I != wav_length){
		printf("ERROR-2! Waveform size not valid: %d (it has to be %d)\n", len_I, wav_length);
		return;
	}
	if((cluster_pos >=  cluster_length)||(cluster_pos < 0)){
		printf("ERROR-3! Position out of range (max = %d, min = 0)\n", (cluster_length - 1));
		return;
	}
	for(j=0; j<len_I; j++){
		Icomponents[cluster_pos][j] = Icomp_in[j];
		Qcomponents[cluster_pos][j] = Qcomp_in[j];
	}
	valid_wavs[cluster_pos] = true;
	num_valid_wavs ++;
	return;
}

void Waveform_complex_cluster::add_waveform_scale( double* Icomp_in, int len_I, double* Qcomp_in, int len_Q, int cluster_pos, double scale_factor ){
	int j;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(len_I != len_Q){
		printf("ERROR-1! I and Q components must have equal size (size-I %d, size-Q %d)\n", len_I, len_Q);
		return;
	}
	if(len_I != wav_length){
		printf("ERROR-2! Waveform size not valid: %d (it has to be %d)\n", len_I, wav_length);
		return;
	}
	if((cluster_pos >=  cluster_length)||(cluster_pos < 0)){
		printf("ERROR-3! Position out of range (max = %d, min = 0)\n", (cluster_length - 1));
		return;
	}
	for(j=0; j<len_I; j++){
		Icomponents[cluster_pos][j] = Icomp_in[j]*scale_factor;
		Qcomponents[cluster_pos][j] = Qcomp_in[j]*scale_factor;
	}
	valid_wavs[cluster_pos] = true;
	num_valid_wavs ++;
	return;
}

void Waveform_complex_cluster::add_waveform_GOLD( signed char* Icomponents_in, int wav_in_length_I, signed char* Qcomponents_in, int wav_in_length_Q, int cluster_pos ){
	int j;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(wav_in_length_I != wav_in_length_Q){
		printf("ERROR-1! I and Q components must have equal size (size-I %d, size-Q %d)\n", wav_in_length_I, wav_in_length_Q);
		return;
	}
	if(wav_in_length_I != wav_length){
		printf("ERROR-2! Waveform size not valid: %d (it has to be %d)\n", wav_in_length_I, wav_length);
		return;
	}
	if((cluster_pos >=  cluster_length)||(cluster_pos < 0)){
		printf("ERROR-3! Position out of range (max = %d, min = 0)\n", (cluster_length - 1));
		return;
	}
	for(j=0; j<wav_in_length_I; j++){
		Icomponents[cluster_pos][j] = double(Icomponents_in[j]);
		Qcomponents[cluster_pos][j] = double(Qcomponents_in[j]);
	}
	valid_wavs[cluster_pos] = true;
	num_valid_wavs ++;
	return;
}

void Waveform_complex_cluster::add_waveform_PIR( short* XiYi, int wav_in_length_II, short* XqYq, int wav_in_length_QQ, short* XiYq, int wav_in_length_IQ, short* XqYi, int wav_in_length_QI, int cluster_pos ){
	int j;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((wav_in_length_II != wav_in_length_QQ)||(wav_in_length_II != wav_in_length_IQ)||(wav_in_length_II != wav_in_length_QI)){
		printf("ERROR-1! II, QQ, IQ and QI components must have equal size (size-II %d, size-QQ %d, size-IQ %d, size-QI %d)\n", wav_in_length_II, wav_in_length_QQ, wav_in_length_IQ, wav_in_length_QI);
		return;
	}
	if(wav_in_length_II != wav_length){
		printf("ERROR-2! Waveform size not valid: %d (it has to be %d)\n", wav_in_length_II, wav_length);
		return;
	}
	if((cluster_pos >=  cluster_length)||(cluster_pos < 0)){
		printf("ERROR-3! Position out of range (max = %d, min = 0)\n", (cluster_length - 1));
		return;
	}
	for(j=0; j<wav_length; j++){
		Icomponents[cluster_pos][j] = double(XiYi[j]) + double(XqYq[j]);
		Qcomponents[cluster_pos][j] = double(XqYi[j]) - double(XiYq[j]);
	}
	valid_wavs[cluster_pos] = true;
	num_valid_wavs ++;
	return;
}

double Waveform_complex_cluster::load_ITF_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double BF_phases_UP[8], double BF_phases_DW[8], int filter_num ){
	int i, wav_spir_length, msecs_spir_file;
	double phasesBF[16];
	bool elementsUP[16], elementsDW[16];
	msecs_spir_file = 999;
	wav_spir_length = 512;
	//Prepare array variables
	prepare_arrays(msecs_spir_file, wav_spir_length);
	for(i=0; i<16; i++){
		if(i<8){
			phasesBF[i] = BF_phases_UP[i];
			elementsUP[i] = true;
			elementsDW[i] = false;
		}else{
			phasesBF[i] = BF_phases_DW[i - 8];
			elementsUP[i] = false;
			elementsDW[i] = true;
		}
	}
	return compute_ITF_waveforms_SPIR(namefile, peak_delay_estimate, phasesBF, elementsUP, elementsDW, filter_num);
}

double Waveform_complex_cluster::load_ITF_waveforms_SPIR_selected_signals( const char* namefile, double peak_delay_estimate, double BF_phases[16], signed char signal_elements_1[16], signed char signal_elements_2[16], int filter_num ){
	int i, wav_spir_length, msecs_spir_file;
	bool elementsUP[16], elementsDW[16];
	msecs_spir_file = 999;
	wav_spir_length = 512;
	//Prepare array variables
	prepare_arrays(msecs_spir_file, wav_spir_length);
	for(i=0; i<16; i++){
		if(signal_elements_1[i] == 1){
			elementsUP[i] = true;
		}else{
			elementsUP[i] = false;
		}
		if(signal_elements_2[i] == 1){
			elementsDW[i] = true;
		}else{
			elementsDW[i] = false;
		}
	}
	return compute_ITF_waveforms_SPIR(namefile, peak_delay_estimate, BF_phases, elementsUP, elementsDW, filter_num);
}

double Waveform_complex_cluster::load_CR_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double doppler_estimate, double delta_doppler, double BF_phases[16], bool uplooking_channel, int code_ref ){
	int i, wav_spir_length, msecs_spir_file;
	bool elements_CORR[16];
	msecs_spir_file = 999;
	wav_spir_length = 512;
	//Prepare array variables
	prepare_arrays(msecs_spir_file, wav_spir_length);
	for(i=0; i<16; i++){
		if(i<8){
			if(uplooking_channel){
				elements_CORR[i] = true;
			}else{
				elements_CORR[i] = false;
			}
		}
		else{
			if(uplooking_channel){
				elements_CORR[i] = false;
			}else{
				elements_CORR[i] = true;
			}
		}
	}
	return compute_CR_waveforms_SPIR(namefile, peak_delay_estimate, doppler_estimate, delta_doppler, BF_phases, elements_CORR, code_ref, false, false, true);
}

double Waveform_complex_cluster::load_CR_waveforms_SPIR_selected_signals( const char* namefile, double peak_delay_estimate, double doppler_estimate, bool assist_doppler, double BF_phases[16], signed char signal_elements[16], int code_ref ){
	int i, wav_spir_length, msecs_spir_file;
	bool elements_CORR[16];
	msecs_spir_file = 999;
	wav_spir_length = 512;
	//Prepare array variables
	prepare_arrays(msecs_spir_file, wav_spir_length);
	for(i=0; i<16; i++){
        if(signal_elements[i] == 1){
            elements_CORR[i] = true;
        }else{
            elements_CORR[i] = false;
        }
	}
	return compute_CR_waveforms_SPIR(namefile, peak_delay_estimate, doppler_estimate, 0.0, BF_phases, elements_CORR, code_ref, true, false, assist_doppler);
}

void Waveform_complex_cluster::searching_CR_waveforms_SPIR( const char* namefile, signed char signal_elements[16] ){
    int i;
    bool elements_CORR[16];
    double *phasesBF;
    phasesBF = (double*) calloc (16, sizeof(double));
    for(i=0; i<16; i++){
        if(signal_elements[i] == 1){
            elements_CORR[i] = true;
        }else{
            elements_CORR[i] = false;
        }
	}
	compute_CR_waveforms_SPIR(namefile, 0.0, 0.0, 0.0, phasesBF, elements_CORR, 0, true, true, true);
}

void Waveform_complex_cluster::get_waveform( double* Icomp_out, int len_I, double* Qcomp_out, int len_Q, int cluster_pos ){
	int j;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(len_I != len_Q){
		printf("ERROR-1! I and Q components must have equal size (size-I %d, size-Q %d)\n", len_I, len_Q);
		return;
	}
	if(len_I != wav_length){
		printf("ERROR-2! Waveform size not valid: %d (it has to be %d)\n", len_I, wav_length);
		return;
	}
	if((cluster_pos >=  cluster_length)||(cluster_pos < 0)){
		printf("ERROR-3! Position out of range (max = %d, min = 0)\n", (cluster_length - 1));
		return;
	}
	if(!valid_wavs[cluster_pos]){
		printf("ERROR-4! Not valid waveform at position %d\n", cluster_pos);
		return;
	}
	for(j=0; j<len_I; j++){
		Icomp_out[j] = Icomponents[cluster_pos][j];
		Qcomp_out[j] = Qcomponents[cluster_pos][j];
	}
	return;
}

void Waveform_complex_cluster::integrate_waveforms( int coherent_int, double* wav_out, int wav_out_length ){
	double *no_retrack=NULL; // NULL, because the integrate_wavs function called with retracking=false
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(wav_out_length != wav_length){
		printf("ERROR! Waveform size not valid: %d (it has to be %d)\n", wav_out_length, wav_length);
		return;
	}
	if(coherent_int > cluster_length){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if(coherent_int == 0){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	integrate_wavs(coherent_int, (cluster_length + 1), cluster_length, wav_length, 0.0, valid_wavs, Icomponents, Qcomponents, wav_out, false, no_retrack);
	return;
}

void Waveform_complex_cluster::integrate_waveforms_remdir( int coherent_int, int coherent_int_dir, double* wav_out, int wav_out_length ){
	double *no_retrack=NULL; // NULL, because the integrate_wavs function called with retracking=false
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(wav_out_length != wav_length){
		printf("ERROR! Waveform size not valid: %d (it has to be %d)\n", wav_out_length, wav_length);
		return;
	}
	if((coherent_int > cluster_length)||(coherent_int_dir > cluster_length)){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if((coherent_int == 0)||(coherent_int_dir == 0)){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	integrate_wavs(coherent_int, coherent_int_dir, cluster_length, wav_length, 0.0, valid_wavs, Icomponents, Qcomponents, wav_out, false, no_retrack);
	return;
}

void Waveform_complex_cluster::integrate_waveforms_retracking( int coherent_int, double sampling_rate, double* retracking_meters, int retracking_meters_length, double* wav_out, int wav_out_length ){
	int i;
	double *retrack_samples;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(wav_out_length != wav_length){
		printf("ERROR! Waveform size not valid: %d (it has to be %d)\n", wav_out_length, wav_length);
		return;
	}
	if(coherent_int > cluster_length){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if(coherent_int == 0){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	if((cluster_length/coherent_int) != retracking_meters_length){
		printf("ERROR! The number of retracking values (%d) must be equal to cluster_size/coherent_in (%d)\n", retracking_meters_length, (cluster_length/coherent_int));
		return;
	}
	retrack_samples = (double *) malloc(retracking_meters_length*sizeof(double));
	for(i=0; i<retracking_meters_length; i++){
		retrack_samples[i] = retracking_meters[i]*sampling_rate/C_LIGHT;
	}
	integrate_wavs(coherent_int, (cluster_length + 1), cluster_length, wav_length, 0.0, valid_wavs, Icomponents, Qcomponents, wav_out, true, retrack_samples);
	free(retrack_samples);
	return;
}

void Waveform_complex_cluster::dump_phase( int lag_pos ){
	int i;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((lag_pos >=  wav_length)||(lag_pos < 0)){
		printf("ERROR! Lag position out of range (max = %d, min = 0)\n", (wav_length - 1));
		return;
	}
	for(i=0; i<cluster_length; i++){
		if(valid_wavs[i]){
			printf("PHASE-L%d %d %f\n",lag_pos,i,atan2(Qcomponents[i][lag_pos],Icomponents[i][lag_pos]));
		}
	}
	return;
}

void Waveform_complex_cluster::dump_phase_peak( void ){
	int i, j, max_pos; 
	double max_val, power;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	for(i=0; i<cluster_length; i++){
		if(valid_wavs[i]){
			max_pos = 0;
			max_val = 0;
			for(j=0; j<wav_length; j++){
				power = Icomponents[i][j]*Icomponents[i][j] + Qcomponents[i][j]*Qcomponents[i][j];
				if(power > max_val){
					max_pos = j;
					max_val = power;
				}
			}
			printf("PHASE-PEAK %d %f\n",i,atan2(Qcomponents[i][max_pos],Icomponents[i][max_pos]));
		}
	}
	return;
}

void Waveform_complex_cluster::store_phasor_wavs( int lag_pos ){
	int i;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((lag_pos >=  wav_length)||(lag_pos < 0)){
		printf("ERROR! Lag position out of range (max = %d, min = 0)\n", (wav_length - 1));
		return;
	}
	for(i=0; i<cluster_length; i++){
		valid_phasor[i] = false;
		if(valid_wavs[i]){
			phasorI[i] = Icomponents[i][lag_pos];
			phasorQ[i] = Qcomponents[i][lag_pos];
			valid_phasor[i] = true;
		}
	}
	num_phasor_iter = 1;
	return;
}

void Waveform_complex_cluster::get_phasor( double* phasorI_out, int phasor_length_I, double* phasorQ_out, int phasor_length_Q, signed char* valid_phasor_out, int phasor_length_bool ){
	int i;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((phasor_length_I != cluster_length)||(phasor_length_Q != cluster_length)||(phasor_length_bool != cluster_length)){
		printf("ERROR! I, Q and BOOL components must have size = %d (size-I %d, size-Q %d, size-BOOL %d)\n", cluster_length, phasor_length_I, phasor_length_Q, phasor_length_bool);
		return;
	}
	for(i=0; i<cluster_length; i++){
		phasorI_out[i] = phasorI[i];
		phasorQ_out[i] = phasorQ[i];
		if(valid_phasor[i]){
			valid_phasor_out[i] = 1;
		}else{
			valid_phasor_out[i] = 0;
		}
	}
	return;
}

double Waveform_complex_cluster::get_sigma_phase_phasor( int min_valid_samples ){
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return -1.0;
	}
	if(min_valid_samples <= 0){
		printf("ERROR! min_valid_samples has to be > 0 in get_sigma_phase_phasor\n");
		return -1.0;
	}
	if(min_valid_samples > cluster_length){
		printf("ERROR! min_valid_samples > cluster_length\n");
		return -1.0;
	}
	return compute_sigma_phase_phasor(0, cluster_length, min_valid_samples, phasorI, phasorQ, valid_phasor);
}

double Waveform_complex_cluster::get_sigma_phase_phasor_interv( int init_sample, int interv_samples, int min_valid_samples ){
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return -1.0;
	}
	if(min_valid_samples <= 0){
		printf("ERROR! min_valid_samples has to be > 0 in get_sigma_phase_phasor\n");
		return -1.0;
	}
	if(min_valid_samples > cluster_length){
		printf("ERROR! min_valid_samples > cluster_length\n");
		return -1.0;
	}
	if(init_sample < 0){
		printf("ERROR! init_sample < 0\n");
		return -1.0;
	}
	if(interv_samples < 2){
		printf("ERROR! interv_samples < 2\n");
		return -1.0;
	}
	if((init_sample + interv_samples) > cluster_length){
		printf("ERROR! (init_sample + interv_samples) > cluster_length\n");
		return -1.0;
	}
	return compute_sigma_phase_phasor(init_sample, interv_samples, min_valid_samples, phasorI, phasorQ, valid_phasor);
}

void Waveform_complex_cluster::counterrot_phasor( double* phases_rad, int phases_length, signed char* valid_phases, int phases_length_bool ){
	int i;
	double phasorI_temp, phasorQ_temp;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((phases_length != cluster_length)||(phases_length_bool != cluster_length)){
		printf("ERROR! Phases and BOOL must have size = %d (size-Phases %d, size-BOOL %d)\n", cluster_length, phases_length, phases_length_bool);
		return;
	}
	for(i=0; i<cluster_length; i++){
		if(valid_phasor[i] && (valid_phases[i]==1)){
			phasorI_temp = phasorI[i]*cos(phases_rad[i]) + phasorQ[i]*sin(phases_rad[i]);
			phasorQ_temp = phasorQ[i]*cos(phases_rad[i]) - phasorI[i]*sin(phases_rad[i]);
			phasorI[i] = phasorI_temp;
			phasorQ[i] = phasorQ_temp;
			valid_phasor[i] = true;
		}else{
			valid_phasor[i] = false;
		}
	}
	num_phasor_iter ++;
	return;
}

void Waveform_complex_cluster::counterrot_waveforms( double* phases_rad, int phases_length, signed char* valid_phases, int phases_length_bool ){
	int i, j;
	double phasorI_temp, phasorQ_temp, phaseI, phaseQ;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((phases_length != cluster_length)||(phases_length_bool != cluster_length)){
		printf("ERROR! Phases and BOOL must have size = %d (size-Phases %d, size-BOOL %d)\n", cluster_length, phases_length, phases_length_bool);
		return;
	}
	num_valid_wavs = 0;
	for(i=0; i<cluster_length; i++){
		if(valid_wavs[i] && (valid_phases[i]==1)){
			phaseI = cos(phases_rad[i]);
			phaseQ = sin(phases_rad[i]);
			for(j=0; j<wav_length; j++){
				phasorI_temp = Icomponents[i][j]*phaseI + Qcomponents[i][j]*phaseQ;
				phasorQ_temp = Qcomponents[i][j]*phaseI - Icomponents[i][j]*phaseQ;
				Icomponents[i][j] = phasorI_temp;
				Qcomponents[i][j] = phasorQ_temp;
			}
			valid_wavs[i] = true;
			num_valid_wavs ++;
		}else{
			valid_wavs[i] = false;
		}
	}
	num_phasor_iter ++;
	return;
}

void Waveform_complex_cluster::correct_navigation_bit( int lag_pos, int store_navbit_phasorI ){
	int i, j;
	int half_nav_period = 10;
	int samples_left, samples_right, bit_sign;
	bool no_near_peaks;
	double accum_phase_left, accum_phase_right, phase_diff;
	double accum_phase_diff, prev1_accum_phase_diff, prev2_accum_phase_diff;
	double phase_bound = double(PI_NUM/2.0);
	double deriv_phase_v1[cluster_length];
	double deriv_phase_v2[cluster_length];
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((lag_pos >=  wav_length)||(lag_pos < 0)){
		printf("ERROR! Lag position out of range (max = %d, min = 0)\n", (wav_length - 1));
		return;
	}
	memset(deriv_phase_v1, 0, sizeof(double)*cluster_length);
	memset(deriv_phase_v2, 0, sizeof(double)*cluster_length);
	accum_phase_diff = 0.0;
	prev1_accum_phase_diff = 0.0;
	prev2_accum_phase_diff = 0.0;
	for(i=1; i<(cluster_length - 1); i++){
		accum_phase_diff = 0.0;
		if(valid_wavs[i]){
			samples_left = 0;
			accum_phase_left = 0.0;
			samples_right = 0;
			accum_phase_right = 0.0;
			for(j=1; j<(half_nav_period + 1); j++){
				if((i - j)>=0){
					if(valid_wavs[i - j]){
						phase_diff = fabs(atan2((Qcomponents[i][lag_pos]*Icomponents[i - j][lag_pos] - Icomponents[i][lag_pos]*Qcomponents[i - j][lag_pos]), (Icomponents[i][lag_pos]*Icomponents[i - j][lag_pos] + Qcomponents[i][lag_pos]*Qcomponents[i - j][lag_pos])));
						accum_phase_left = accum_phase_left + phase_diff;
						samples_left ++;
					}
				}
				if((i + j)<cluster_length){
					if(valid_wavs[i + j]){
						phase_diff = fabs(atan2((Qcomponents[i][lag_pos]*Icomponents[i + j][lag_pos] - Icomponents[i][lag_pos]*Qcomponents[i + j][lag_pos]), (Icomponents[i][lag_pos]*Icomponents[i + j][lag_pos] + Qcomponents[i][lag_pos]*Qcomponents[i + j][lag_pos])));
						accum_phase_right = accum_phase_right + phase_diff;
						samples_right ++;
					}
				}
			}
			if((samples_left > 0)&&(samples_right > 0)){
				accum_phase_diff = (accum_phase_left/double(samples_left)) - (accum_phase_right/double(samples_right));
				deriv_phase_v1[i] = accum_phase_diff - prev1_accum_phase_diff;
				if(i > 1){
					deriv_phase_v2[i - 1] = accum_phase_diff - prev2_accum_phase_diff;
				}
			}
		}
		prev2_accum_phase_diff = prev1_accum_phase_diff;
		prev1_accum_phase_diff = accum_phase_diff;
	}
	//First sample
	bit_sign = 1;
	if(store_navbit_phasorI == 1){
		phasorI[0] = 1.0;
	}
	for(i=1; i<(cluster_length - 1); i++){
		if((deriv_phase_v1[i] > phase_bound)&&(deriv_phase_v2[i] > phase_bound)){
			no_near_peaks =  true;
			for(j=1; j<(half_nav_period + 1); j++){
				if((i - j)>=0){
					if((deriv_phase_v1[i - j] > deriv_phase_v1[i])&&(deriv_phase_v2[i - j] > phase_bound)){
						no_near_peaks = false;
					}
				}
				if((i + j)<cluster_length){
					if((deriv_phase_v1[i + j] > deriv_phase_v1[i])&&(deriv_phase_v2[i + j] > phase_bound)){
						no_near_peaks = false;
					}
				}
			}
			if(no_near_peaks){ //CHANGE OF NAVIGATION BIT
				bit_sign = -1*bit_sign;
			}
		}
		if(store_navbit_phasorI == 1){
			phasorI[i] = double(bit_sign);
		}
		if(bit_sign == -1){
			for(j=0; j<wav_length; j++){
				Icomponents[i][j] = -Icomponents[i][j];
				Qcomponents[i][j] = -Qcomponents[i][j];
			}
		}
	}
	//Last sample
	if(store_navbit_phasorI == 1){
		phasorI[cluster_length - 1] = double(bit_sign);
	}
	if(bit_sign == -1){
		for(j=0; j<wav_length; j++){
			Icomponents[cluster_length - 1][j] = -Icomponents[cluster_length - 1][j];
			Qcomponents[cluster_length - 1][j] = -Qcomponents[cluster_length - 1][j];
		}
	}
	return;
}

double Waveform_complex_cluster::compute_coherence_time( int lag_peak, int store_acf_phasorQ ){
	int i, j, prev_pos, num_samples;
	double acf[cluster_length];
	bool valid_acf[cluster_length];
	double accum_acf, phase_diff, prev_acf, coh_sample, coh_bound;
	bool val_not_found;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return -1.0;
	}
	memset(acf, 0, sizeof(double)*cluster_length);
	acf[0] = 0.0;
	valid_acf[0] = true;
	if(store_acf_phasorQ == 1){
		phasorQ[0] = 0.0;
		valid_phasor[0] = true;
	}
	for(i=1; i<cluster_length; i++){
		accum_acf = 0.0;
		num_samples = 0;
		for(j=0; j<(cluster_length - i); j++){
			if(valid_wavs[i]&&valid_wavs[j]){
				phase_diff = fabs(atan2((Qcomponents[j + i][lag_peak]*Icomponents[j][lag_peak] - Icomponents[j + i][lag_peak]*Qcomponents[j][lag_peak]), (Icomponents[j + i][lag_peak]*Icomponents[j][lag_peak] + Qcomponents[j + i][lag_peak]*Qcomponents[j][lag_peak])));
				accum_acf = accum_acf + phase_diff;
				num_samples ++;
			}
		}
		if(num_samples > 0){
			acf[i] = accum_acf/double(num_samples);
			valid_acf[i] = true;
			if(store_acf_phasorQ == 1){
				phasorQ[i] = acf[i];
				valid_phasor[i] = true;
			}
		}else{
			valid_acf[i] = false;
			if(store_acf_phasorQ == 1){
				valid_phasor[i] = false;
			}
		}
	}
	coh_bound = 1; //In radians
	prev_pos = 0;
	val_not_found = true;
	for(i=1; i<cluster_length; i++){
		if((valid_acf[i])&&val_not_found){
			if(acf[i] >= coh_bound){
				coh_sample = prev_pos + ((i - prev_pos)*(coh_bound - prev_acf)/(acf[i] - prev_acf)); //Linear interpolation
				val_not_found = false;
			}else{
				prev_acf = acf[i];
				prev_pos = i;
			}
		}
	}
	if(val_not_found){
		return double(cluster_length);
	}else{
		return coh_sample;
	}
}

void Waveform_complex_cluster::compute_singlefreq_DDM( int coherent_int, double doppler_freq, double* freq_ddm, int delay_samples ){
	double *no_retrack=NULL; // NULL, because the integrate_wavs function called with retracking=false
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(delay_samples != wav_length){
		printf("ERROR! DDM delay-size not valid: %d (it has to be %d)\n", delay_samples, wav_length);
		return;
	}
	if(coherent_int > cluster_length){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if(coherent_int == 0){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	integrate_wavs(coherent_int, (cluster_length + 1), cluster_length, wav_length, doppler_freq, valid_wavs, Icomponents, Qcomponents, freq_ddm, false, no_retrack);
	return;
}

void Waveform_complex_cluster::compute_singlefreq_DDM_remdir( int coherent_int, int coherent_int_dir, double doppler_freq, double* freq_ddm, int delay_samples ){
	double *no_retrack=NULL; // NULL, because the integrate_wavs function called with retracking=false
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(delay_samples != wav_length){
		printf("ERROR! DDM delay-size not valid: %d (it has to be %d)\n", delay_samples, wav_length);
		return;
	}
	if((coherent_int > cluster_length)||(coherent_int_dir > cluster_length)){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if((coherent_int == 0)||(coherent_int_dir == 0)){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	integrate_wavs(coherent_int, coherent_int_dir, cluster_length, wav_length, doppler_freq, valid_wavs, Icomponents, Qcomponents, freq_ddm, false, no_retrack);
	return;
}

void Waveform_complex_cluster::compute_singlelag_DDM( int coherent_int, int lag_pos, double delta_freq, double* lag_ddm, int freq_samples ){
	int i;
	double doppler_freq;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((lag_pos < 0)||(lag_pos >= wav_length)){
		printf("ERROR! DDM lag position not valid: %d (it has to be < %d)\n", lag_pos, wav_length);
		return;
	}
	if(coherent_int > cluster_length){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if(coherent_int == 0){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	for(i=0; i<freq_samples; i++){
		doppler_freq = double(i - freq_samples/2)*delta_freq;
		lag_ddm[i] = integrate_wav_lag(coherent_int, (cluster_length + 1), cluster_length, lag_pos, doppler_freq, valid_wavs, Icomponents, Qcomponents);
	}
	return;
}

void Waveform_complex_cluster::compute_singlelag_DDM_remdir( int coherent_int, int coherent_int_dir, int lag_pos, double delta_freq, double* lag_ddm, int freq_samples ){
	int i;
	double doppler_freq;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if((lag_pos < 0)||(lag_pos >= wav_length)){
		printf("ERROR! DDM lag position not valid: %d (it has to be < %d)\n", lag_pos, wav_length);
		return;
	}
	if((coherent_int > cluster_length)||(coherent_int_dir > cluster_length)){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if((coherent_int == 0)||(coherent_int_dir == 0)){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	for(i=0; i<freq_samples; i++){
		doppler_freq = double(i - freq_samples/2)*delta_freq;
		lag_ddm[i] = integrate_wav_lag(coherent_int, coherent_int_dir, cluster_length, lag_pos, doppler_freq, valid_wavs, Icomponents, Qcomponents);
	}
	return;
}

double Waveform_complex_cluster::compute_DopplerMap_BW( int coherent_int, int lag_pos, int freq_samples, double delta_freq, double pos_pow_Max[2] ){
	int i;
	double doppler_freqs[freq_samples];
	double doppler_map[freq_samples];
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return 0.0;
	}
	if((lag_pos < 0)||(lag_pos >= wav_length)){
		printf("ERROR! DDM lag position not valid: %d (it has to be < %d)\n", lag_pos, wav_length);
		return 0.0;
	}
	if(coherent_int > cluster_length){
		printf("ERROR! Coherent integration time > cluster length\n");
		return 0.0;
	}
	if(coherent_int == 0){
		printf("ERROR! Coherent integration time = 0\n");
		return 0.0;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return 0.0;
	}
	for(i=0; i<freq_samples; i++){
		doppler_freqs[i] = double(i - freq_samples/2)*delta_freq;
		doppler_map[i] = integrate_wav_lag(coherent_int, (cluster_length + 1), cluster_length, lag_pos, doppler_freqs[i], valid_wavs, Icomponents, Qcomponents);
	}
	return Compute_Doppler_bandwidth(freq_samples, doppler_freqs, doppler_map, pos_pow_Max);
}

double Waveform_complex_cluster::compute_DopplerMap_BW_remdir( int coherent_int, int coherent_int_dir, int lag_pos, int freq_samples, double delta_freq, double pos_pow_Max[2] ){
	int i;
	double doppler_freqs[freq_samples];
	double doppler_map[freq_samples];
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return 0.0;
	}
	if((lag_pos < 0)||(lag_pos >= wav_length)){
		printf("ERROR! DDM lag position not valid: %d (it has to be < %d)\n", lag_pos, wav_length);
		return 0.0;
	}
	if((coherent_int > cluster_length)||(coherent_int_dir > cluster_length)){
		printf("ERROR! Coherent integration time > cluster length\n");
		return 0.0;
	}
	if((coherent_int == 0)||(coherent_int_dir == 0)){
		printf("ERROR! Coherent integration time = 0\n");
		return 0.0;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return 0.0;
	}
	for(i=0; i<freq_samples; i++){
		doppler_freqs[i] = double(i - freq_samples/2)*delta_freq;
		doppler_map[i] = integrate_wav_lag(coherent_int, coherent_int_dir, cluster_length, lag_pos, doppler_freqs[i], valid_wavs, Icomponents, Qcomponents);
	}
	return Compute_Doppler_bandwidth(freq_samples, doppler_freqs, doppler_map, pos_pow_Max);
}

void Waveform_complex_cluster::compute_whole_DDM( int coherent_int, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples ){
	// THIS FUNCTION DOES NOT WORK IN PYTHON DUE TO 2D-ARRAY LIMITATIONS IN NUMPY-SWIG
	int j, i_freq;
	double doppler_freq;
	double pow_wav[wav_length];
	double *no_retrack=NULL; // NULL, because the integrate_wavs function called with retracking=false
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(ddm_delay_samples != wav_length){
		printf("ERROR! DDM delay-size not valid: %d (it has to be %d)\n", ddm_delay_samples, wav_length);
		return;
	}
	if(coherent_int > cluster_length){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if(coherent_int == 0){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	for(i_freq=0; i_freq<ddm_freq_samples; i_freq++){
		doppler_freq = double(i_freq - ddm_freq_samples/2)*ddm_delta_freq;
		integrate_wavs(coherent_int, (cluster_length + 1), cluster_length, wav_length, doppler_freq, valid_wavs, Icomponents, Qcomponents, pow_wav, false, no_retrack);
		for(j=0; j<wav_length; j++){
			ddm[j][i_freq] = pow_wav[j];
		}
	}
	return;
}

void Waveform_complex_cluster::compute_whole_DDM_retracking( int coherent_int, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples, double sampling_rate, double* retracking_meters, int retracking_meters_length ){
	// THIS FUNCTION DOES NOT WORK IN PYTHON DUE TO 2D-ARRAY LIMITATIONS IN NUMPY-SWIG
	int j, i_freq;
	double doppler_freq;
	double pow_wav[wav_length];
	double *retrack_samples;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(ddm_delay_samples != wav_length){
		printf("ERROR! DDM delay-size not valid: %d (it has to be %d)\n", ddm_delay_samples, wav_length);
		return;
	}
	if(coherent_int > cluster_length){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if(coherent_int == 0){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	retrack_samples = (double *) malloc(retracking_meters_length*sizeof(double));
	for(j=0; j<retracking_meters_length; j++){
		retrack_samples[j] = retracking_meters[j]*sampling_rate/C_LIGHT;
	}
	for(i_freq=0; i_freq<ddm_freq_samples; i_freq++){
		doppler_freq = double(i_freq - ddm_freq_samples/2)*ddm_delta_freq;
		integrate_wavs(coherent_int, (cluster_length + 1), cluster_length, wav_length, doppler_freq, valid_wavs, Icomponents, Qcomponents, pow_wav, true, retrack_samples);
		for(j=0; j<wav_length; j++){
			ddm[j][i_freq] = pow_wav[j];
		}
	}
	free(retrack_samples);
	return;
}

void Waveform_complex_cluster::compute_whole_DDM_remdir( int coherent_int, int coherent_int_dir, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples ){
	// THIS FUNCTION DOES NOT WORK IN PYTHON DUE TO 2D-ARRAY LIMITATIONS IN NUMPY-SWIG
	int j, i_freq;
	double doppler_freq;
	double pow_wav[wav_length];
	double *no_retrack=NULL; // NULL, because the integrate_wavs function called with retracking=false
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(ddm_delay_samples != wav_length){
		printf("ERROR! DDM delay-size not valid: %d (it has to be %d)\n", ddm_delay_samples, wav_length);
		return;
	}
	if((coherent_int > cluster_length)||(coherent_int_dir > cluster_length)){
		printf("ERROR! Coherent integration time > cluster length\n");
		return;
	}
	if((coherent_int == 0)||(coherent_int_dir == 0)){
		printf("ERROR! Coherent integration time = 0\n");
		return;
	}
	if(num_valid_wavs == 0){
		printf("ERROR! No valid waveforms in cluster\n");
		return;
	}
	for(i_freq=0; i_freq<ddm_freq_samples; i_freq++){
		doppler_freq = double(i_freq - ddm_freq_samples/2)*ddm_delta_freq;
		integrate_wavs(coherent_int, coherent_int_dir, cluster_length, wav_length, doppler_freq, valid_wavs, Icomponents, Qcomponents, pow_wav, false, no_retrack);
		for(j=0; j<wav_length; j++){
			ddm[j][i_freq] = pow_wav[j];
		}
	}
	return;
}

void Waveform_complex_cluster::compute_whole_LagHologram( double** powerLagHolo, int lagHolo_lags, int lagHolo_freqs ){
	int power2, lag, i;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(lagHolo_lags != wav_length){
		printf("ERROR! Lag-Hologram lags != waveform's length\n");
		return;
	}
	if(lagHolo_freqs < 2){
		printf("ERROR! Lag-Hologram freqs < 2\n");
		return;
	}
	power2 = 2;
	while(lagHolo_freqs > power2) power2 = power2*2;
	if(lagHolo_freqs != power2){
		printf("ERROR! Lag-Hologram freqs must be a power of 2\n");
		return;
	}
	double real_part[lagHolo_freqs], imag_part[lagHolo_freqs], power_spectrum[lagHolo_freqs];
	if(lagHolo_freqs > num_valid_wavs){
		printf("WARNING! Zero-padding applied (%d/%d)\n", (lagHolo_freqs - num_valid_wavs), lagHolo_freqs);
	}
	for(lag=0; lag<wav_length; lag++){
		for(i=0; i<lagHolo_freqs; i++){
			real_part[i] = 0.0;
			imag_part[i] = 0.0;
			if(i<cluster_length){
				if(valid_wavs[i]){
					real_part[i] = Icomponents[i][lag];
					imag_part[i] = Qcomponents[i][lag];
				}
			}
		}
		Compute_power_spectrum(lagHolo_freqs, real_part, imag_part, power_spectrum);
		for(i=0; i<lagHolo_freqs; i++){
			powerLagHolo[lag][i] = power_spectrum[i];
		}
	}
	return;
}

void Waveform_complex_cluster::compute_LagHologram( int lag, double* powerLagHolo_singleLag, int fft_samples ){
	int power2, i;
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return;
	}
	if(fft_samples < 2){
		printf("ERROR! Lag-Hologram freqs < 2\n");
		return;
	}
	power2 = 2;
	while(fft_samples > power2) power2 = power2*2;
	if(fft_samples != power2){
		printf("ERROR! Lag-Hologram freqs must be a power of 2\n");
		return;
	}
	double real_part[fft_samples], imag_part[fft_samples];
	if(fft_samples > num_valid_wavs){
		printf("WARNING! Zero-padding applied (%d/%d)\n", (fft_samples - num_valid_wavs), fft_samples);
	}
	for(i=0; i<fft_samples; i++){
		real_part[i] = 0.0;
		imag_part[i] = 0.0;
		if(i<cluster_length){
			if(valid_wavs[i]){
				real_part[i] = Icomponents[i][lag];
				imag_part[i] = Qcomponents[i][lag];
			}
		}
	}
	Compute_power_spectrum(fft_samples, real_part, imag_part, powerLagHolo_singleLag);
	return;
}

int Waveform_complex_cluster::get_wav_length( void ){
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return 0;
	}
	return wav_length;
}

int Waveform_complex_cluster::get_cluster_length( void ){
	if(cluster_length == 0){
		printf("ERROR! Cluster not initialized.\n");
		return 0;
	}
	return cluster_length;
}


void integrate_wavs( int coh_int, int coh_int_dir, int cl_length, int w_length, double freq, bool* valid, double** Icomp, double** Qcomp, double* wav, bool retracking, double* retrack_samples ){
	int i, j, num_coh_wavs, num_incoh_wavs, num_coh_dir_wavs;
	double scale1, arg;
	double *accum_I, *accum_Q, *pow_wav, *accum_I_dir, *accum_Q_dir, *pow_wav_dir, *wav_retrack;
	accum_I = (double*) calloc(w_length, sizeof(double));
	accum_Q = (double*) calloc(w_length, sizeof(double));
	pow_wav = (double*) calloc(w_length, sizeof(double));
	accum_I_dir = (double*) calloc(w_length, sizeof(double));
	accum_Q_dir = (double*) calloc(w_length, sizeof(double));
	pow_wav_dir = (double*) calloc(w_length, sizeof(double));
	if(retracking){
		wav_retrack = (double *) malloc(w_length*sizeof(double));
	}
	memset(wav, 0, sizeof(double)*w_length);
	num_coh_wavs = 0;
	num_coh_dir_wavs = 0;
	num_incoh_wavs = 0;
	for(i=0; i<cl_length; i++){
		if(valid[i]){
			arg = 2.0*PI_NUM*freq*double(i)/double(cl_length);
			for(j=0; j<w_length; j++){
				accum_I[j] = accum_I[j] + double(Icomp[i][j])*cos(arg) + double(Qcomp[i][j])*sin(arg);
				accum_Q[j] = accum_Q[j] + double(Qcomp[i][j])*cos(arg) - double(Icomp[i][j])*sin(arg);
				if(coh_int_dir <= cl_length){
					accum_I_dir[j] = accum_I_dir[j] + double(Icomp[i][j])*cos(arg) + double(Qcomp[i][j])*sin(arg);
					accum_Q_dir[j] = accum_Q_dir[j] + double(Qcomp[i][j])*cos(arg) - double(Icomp[i][j])*sin(arg);
				}
			}
			num_coh_wavs ++;
			if(coh_int_dir <= cl_length){
				num_coh_dir_wavs ++;
			}
		}
		if((i%coh_int) == (coh_int - 1)){
			if(num_coh_wavs > 0){
				//scale1 = pow(1./(8.e+4), 2)/(num_coh_wavs*num_coh_wavs);
				scale1 = 1.0/double(num_coh_wavs);
				if(retracking){
					for(j=0; j<w_length; j++){
						wav_retrack[j] = scale1*((accum_I[j]*accum_I[j]) + (accum_Q[j]*accum_Q[j]));
					}
					retrack_real_waveform(wav_retrack, w_length, retrack_samples[num_incoh_wavs]);
				}
				for(j=0; j<w_length; j++){
					if(retracking){
						pow_wav[j] = pow_wav[j] + wav_retrack[j];
					}else{
						pow_wav[j] = pow_wav[j] + scale1*((accum_I[j]*accum_I[j]) + (accum_Q[j]*accum_Q[j]));
					}
				}
				num_incoh_wavs ++;
			}
			memset(accum_I, 0, sizeof(double)*w_length);
			memset(accum_Q, 0, sizeof(double)*w_length);
			num_coh_wavs = 0;
		}
		if((i%coh_int_dir) == (coh_int_dir - 1)){
			if(num_coh_dir_wavs > 0){
				//scale1 = pow(1./(8.e+4), 2)/(num_coh_wavs*num_coh_wavs);
				scale1 = 1.0/double(num_coh_dir_wavs);
				for(j=0; j<w_length; j++){
					pow_wav_dir[j] = pow_wav_dir[j] + scale1*((accum_I_dir[j]*accum_I_dir[j]) + (accum_Q_dir[j]*accum_Q_dir[j]));
				}
			}
			memset(accum_I_dir, 0, sizeof(double)*w_length);
			memset(accum_Q_dir, 0, sizeof(double)*w_length);
			num_coh_dir_wavs = 0;
		}
	}
	if(num_incoh_wavs > 0){
		for(j=0; j<w_length; j++){
			wav[j] = fabs(pow_wav[j] - pow_wav_dir[j]);
		}
	}
	free(accum_I);
	free(accum_Q);
	free(pow_wav);
	free(accum_I_dir);
	free(accum_Q_dir);
	free(pow_wav_dir);
	if(retracking){
		free(wav_retrack);
	}
	return;
}

double integrate_wav_lag( int coh_int, int coh_int_dir, int cl_length, int lag, double freq, bool* valid, double** Icomp, double** Qcomp ){
	int i, num_coh_wavs, num_incoh_wavs, num_coh_dir_wavs;
	double scale1, arg;
	double accum_I = 0.0;
	double accum_Q = 0.0;
	double pow_wav = 0.0;
	double accum_I_dir = 0.0;
	double accum_Q_dir = 0.0;
	double pow_wav_dir = 0.0;
	num_coh_wavs = 0;
	num_coh_dir_wavs = 0;
	num_incoh_wavs = 0;
	for(i=0; i<cl_length; i++){
		if(valid[i]){
			arg = 2.0*PI_NUM*freq*double(i)/double(cl_length);
			accum_I = accum_I + double(Icomp[i][lag])*cos(arg) + double(Qcomp[i][lag])*sin(arg);
			accum_Q = accum_Q + double(Qcomp[i][lag])*cos(arg) - double(Icomp[i][lag])*sin(arg);
			if(coh_int_dir <= cl_length){
				accum_I_dir = accum_I_dir + double(Icomp[i][lag])*cos(arg) + double(Qcomp[i][lag])*sin(arg);
				accum_Q_dir = accum_Q_dir + double(Qcomp[i][lag])*cos(arg) - double(Icomp[i][lag])*sin(arg);
			}
			num_coh_wavs ++;
			if(coh_int_dir <= cl_length){
				num_coh_dir_wavs ++;
			}
		}
		if((i%coh_int) == (coh_int - 1)){
			if(num_coh_wavs > 0){
				//scale1 = pow(1./(8.e+4), 2)/(num_coh_wavs*num_coh_wavs);
				scale1 = 1.0/double(num_coh_wavs);
				pow_wav = pow_wav + scale1*((accum_I*accum_I) + (accum_Q*accum_Q));
				num_incoh_wavs ++;
			}
			accum_I = 0.0;
			accum_Q = 0.0;
			num_coh_wavs = 0;
		}
		if((i%coh_int_dir) == (coh_int_dir - 1)){
			if(num_coh_dir_wavs > 0){
				//scale1 = pow(1./(8.e+4), 2)/(num_coh_wavs*num_coh_wavs);
				scale1 = 1.0/double(num_coh_dir_wavs);
				pow_wav_dir = pow_wav_dir + scale1*((accum_I_dir*accum_I_dir) + (accum_Q_dir*accum_Q_dir));
			}
			accum_I_dir = 0.0;
			accum_Q_dir = 0.0;
			num_coh_dir_wavs = 0;
		}
	}
	if(num_incoh_wavs == 0){
		return 0.0;
	}else{
		return fabs(pow_wav - pow_wav_dir)/double(num_incoh_wavs);
	}
}

double compute_sigma_phase_phasor( int init, int interval, int min_samples, double* phasI, double* phasQ, bool* valid ){
	int i, num_valid_samples;
	double accum_I, accum_Q, mean_I, mean_Q, accum_sigma2_I, accum_sigma2_Q;
	num_valid_samples = 0;
	accum_I = 0.0;
	accum_Q = 0.0;
	for(i=init; i<(init + interval); i++){
		if(valid[i]){
			num_valid_samples ++;
			accum_I = accum_I + phasI[i];
			accum_Q = accum_Q + phasQ[i];
		}
	}
	if(num_valid_samples < min_samples){
		return -1.0;
	}
	if((accum_I == 0.0)&&(accum_Q == 0.0)){
		return -1.0;
	}
	mean_I = accum_I/num_valid_samples;
	mean_Q = accum_Q/num_valid_samples;
	accum_sigma2_I = 0.0;
	accum_sigma2_Q = 0.0;
	for(i=init; i<(init + interval); i++){
		if(valid[i]){
			accum_sigma2_I = accum_sigma2_I + ((phasI[i] - mean_I)*(phasI[i] - mean_I));
			accum_sigma2_Q = accum_sigma2_Q + ((phasQ[i] - mean_Q)*(phasQ[i] - mean_Q));
		}
	}
	return atan2(sqrt((accum_sigma2_I + accum_sigma2_Q)/num_valid_samples), sqrt(mean_I*mean_I + mean_Q*mean_Q));
}

bool check_SPIR_blocks( int* spirdata ){
	int i, word, blockID, anciID;
	const int spir_data_blocks = 4882;
	const int spir_block_ints = 16384;
	for(i=0; i<spir_data_blocks; i++){
		word = spirdata[i*spir_block_ints];
		blockID = word&0x0000FFFF;
		anciID = (word>>30)&0x00000003;
		if((anciID!=0)||(blockID!=i)){
			printf("ERROR! Unexpected value at block label %d\n", i);
			return false;
		}
	}
	return true;
}

void compute_RF_filter( float* filter ){
	int i, j, pos, pos_min, pos_max;
	const int n_interf = 11;
	const int n_0 = n_interf/2;
	double sampling_rate = 80000000.0;
	const int num_samples = int(sampling_rate)/1000; //1 msec
	double freqList[n_interf];
	double main_freq = FREQ_GPS_L1;
	double df = sampling_rate/double(num_samples);
	double freq_comb = 4000000.0;
	double deltaF = 100000.0;
	double freq_cut = 12000000.0;
	double freq_low = 1000000.0;
	double BW_deg = 20.0;
	double freq, freq_offset;
	for(i=0; i<num_samples; i++){
		filter[i] = 1.0;
	}
	freq_offset = int(main_freq/freq_comb)*freq_comb - main_freq;
	for(i=0; i<n_interf; i++){
			freqList[i] = (i - n_0)*freq_comb + freq_offset;
	}
	for(i=0; i<n_interf; i++){
		pos = int(freqList[i]/df);
		pos_min = pos - int(deltaF/df);
		pos_max = pos + int(deltaF/df);
		if((pos_min<0)&&(pos_max<0)){
			for(j=(pos_min + num_samples); j<=(pos_max + num_samples); j++){
				filter[j] = 0.0;
			}
		}
		if((pos_min<0)&&(pos_max>=0)){
			for(j=(pos_min + num_samples); j<num_samples; j++){
				filter[j] = 0.0;
			}
			for(j=0; j<=pos_max; j++){
				filter[j] = 0.0;
			}
		}
		if((pos_min>=0)&&(pos_max>0)){
			for(j=pos_min; j<=pos_max; j++){
				filter[j] = 0.0;
			}
		}
	}
	for(i=0; i<num_samples; i++){
		if(i > num_samples/2){
			freq = df*double(i - num_samples);
		}else{
			freq = df*double(i);
		}
		filter[i] = filter[i]*float(1.0/(1.0 + pow(freq/freq_cut, 2.0*BW_deg)));
		filter[i] = filter[i]*float(1.0 - (1.0/(1.0 + pow(freq/freq_low, 2.0*BW_deg))));
	}
	return;
}

void compute_RF_filter_GalileoE1( float* filter ){
	int i, j, pos, pos_min, pos_max;
	const int n_interf = 81;
	double sampling_rate = 80000000.0;
	const int num_samples = int(sampling_rate)/1000; //1 msec
	double freqList[n_interf];
	double freqWidths[n_interf];
	freqList[0] = 750;
	freqWidths[0] = 750;
	freqList[1] = 319249;
	freqWidths[1] = 750;
	freqList[2] = 2321;
	freqWidths[2] = 16;
	freqList[3] = 18321;
	freqWidths[3] = 16;
	freqList[4] = 34321;
	freqWidths[4] = 16;
	freqList[5] = 50321;
	freqWidths[5] = 16;
	freqList[6] = 61944;
	freqWidths[6] = 16;
	freqList[7] = 66321;
	freqWidths[7] = 16;
	freqList[8] = 98323;
	freqWidths[8] = 16;
	freqList[9] = 29609;
	freqWidths[9] = 21109;
	freqList[10] = 290388;
	freqWidths[10] = 21109;
	freqList[11] = 160000;
	freqWidths[11] = 88000;
	freqList[12] = 194323;
	freqWidths[12] = 16;
	freqList[13] = 210323;
	freqWidths[13] = 16;
	freqList[14] = 218289;
	freqWidths[14] = 16;
	freqList[15] = 222693;
	freqWidths[15] = 16;
	freqList[16] = 226323;
	freqWidths[16] = 16;
	freqList[17] = 258321;
	freqWidths[17] = 16;
	freqList[18] = 274321;
	freqWidths[18] = 16;
	freqList[19] = 290321;
	freqWidths[19] = 16;
	freqList[20] = 306321;
	freqWidths[20] = 16;
	freqList[21] = 308674;
	freqWidths[21] = 16;
	freqList[22] = 318216;
	freqWidths[22] = 16;
	freqList[23] = 887;
	freqWidths[23] = 16;
	freqList[24] = 923;
	freqWidths[24] = 16;
	freqList[25] = 847;
	freqWidths[25] = 16;
	freqList[26] = 22267;
	freqWidths[26] = 16;
	freqList[27] = 22281;
	freqWidths[27] = 16;
	freqList[28] = 22310;
	freqWidths[28] = 16;
	freqList[29] = 22345;
	freqWidths[29] = 16;
	freqList[30] = 22409;
	freqWidths[30] = 16;
	freqList[31] = 22663;
	freqWidths[31] = 16;
	freqList[32] = 23045;
	freqWidths[32] = 16;
	freqList[33] = 23299;
	freqWidths[33] = 16;
	freqList[34] = 23427;
	freqWidths[34] = 16;
	freqList[35] = 23533;
	freqWidths[35] = 16;
	freqList[36] = 23681;
	freqWidths[36] = 16;
	freqList[37] = 23809;
	freqWidths[37] = 16;
	freqList[38] = 318260;
	freqWidths[38] = 16;
	freqList[39] = 318381;
	freqWidths[39] = 16;
	freqList[40] = 318431;
	freqWidths[40] = 16;
	freqList[41] = 318743;
	freqWidths[41] = 16;
	freqList[42] = 319219;
	freqWidths[42] = 16;
	freqList[43] = 319083;
	freqWidths[43] = 16;
	freqList[44] = 319159;
	freqWidths[44] = 16;
	freqList[45] = 319583;
	freqWidths[45] = 16;
	freqList[46] = 319503;
	freqWidths[46] = 16;
	freqList[47] = 319667;
	freqWidths[47] = 16;
	freqList[48] = 319671;
	freqWidths[48] = 16;
	freqList[49] = 319710;
	freqWidths[49] = 16;
	freqList[50] = 319759;
	freqWidths[50] = 16;
	freqList[51] = 319798;
	freqWidths[51] = 16;
	freqList[52] = 319939;
	freqWidths[52] = 16;
	freqList[53] = 319879;
	freqWidths[53] = 16;
	freqList[54] = 319978;
	freqWidths[54] = 16;
	freqList[55] = 20000;
	freqWidths[55] = 20;
	freqList[56] = 299999;
	freqWidths[56] = 20;
	freqList[57] = 268850;
	freqWidths[57] = 16;
	freqList[58] = 267584;
	freqWidths[58] = 16;
	freqList[59] = 266318;
	freqWidths[59] = 16;
	freqList[60] = 265052;
	freqWidths[60] = 16;
	freqList[61] = 263786;
	freqWidths[61] = 16;
	freqList[62] = 262520;
	freqWidths[62] = 16;
	freqList[63] = 261254;
	freqWidths[63] = 16;
	freqList[64] = 259988;
	freqWidths[64] = 16;
	freqList[65] = 258722;
	freqWidths[65] = 16;
	freqList[66] = 257456;
	freqWidths[66] = 16;
	freqList[67] = 256190;
	freqWidths[67] = 16;
	freqList[68] = 254924;
	freqWidths[68] = 16;
	freqList[69] = 253658;
	freqWidths[69] = 16;
	freqList[70] = 252392;
	freqWidths[70] = 16;
	freqList[71] = 251126;
	freqWidths[71] = 16;
	freqList[72] = 251126;
	freqWidths[72] = 16;
	freqList[73] = 249860;
	freqWidths[73] = 16;
	freqList[74] = 248594;
	freqWidths[74] = 16;
	freqList[75] = 318137;
	freqWidths[75] = 16;
	freqList[76] = 1574;
	freqWidths[76] = 16;
	freqList[77] = 58215;
	freqWidths[77] = 16;
	freqList[78] = 58272;
	freqWidths[78] = 16;
	freqList[79] = 58200;
	freqWidths[79] = 16;
	freqList[80] = 317726;
	freqWidths[80] = 32;
	for(i=0; i<num_samples; i++){
		filter[i] = 1.0;
	}
	for(i=0; i<n_interf; i++){
		pos = freqList[i]/4;
		pos_min = pos - freqWidths[i]/4;
		pos_max = pos + freqWidths[i]/4;
		for(j=pos_min; j<=pos_max; j++){
			filter[j] = 0.0;
		}
	}
	return;
}

void compute_RF_filter_GPS_L5( float* filter ){
	int i, j, pos, pos_min, pos_max;
	const int n_interf = 22;
	double sampling_rate = 80000000.0;
	const int num_samples = int(sampling_rate)/1000; //1 msec
	double freqList[n_interf];
	double freqWidths[n_interf];
	freqList[0] = 5;
	freqWidths[0] = 5; 
	freqList[1] = 3800; 
	freqWidths[1] = 600; 
	freqList[2] = 6542;
	freqWidths[2] = 200;
	freqList[3] = 9264;
	freqWidths[3] = 100;
	freqList[4] = 10231;
	freqWidths[4] = 15;
	freqList[5] = 11582;
	freqWidths[5] = 15;
	freqList[6] = 12516;
	freqWidths[6] = 300;
	freqList[7] = 17836;
	freqWidths[7] = 10;
	freqList[8] = 40000;
	freqWidths[8] = 30000;
	freqList[9] = 62603;
	freqWidths[9] = 300;
	freqList[10] = 65576;
	freqWidths[10] = 500;
	freqList[11] = 68642;
	freqWidths[11] = 100;
	freqList[12] = 69771;
	freqWidths[12] = 15;
	freqList[13] = 70335;
	freqWidths[13] = 15;
	freqList[14] = 70442;
	freqWidths[14] = 15;
	freqList[15] = 70659;
	freqWidths[15] = 30;
	freqList[16] = 70610;
	freqWidths[16] = 15;
	freqList[17] = 70506;
	freqWidths[17] = 200;
	freqList[18] = 74580;
	freqWidths[18] = 300;
	freqList[19] = 74979;
	freqWidths[19] = 15;
	freqList[20] = 79550;
	freqWidths[20] = 400;
	freqList[21] = 79994;
	freqWidths[21] = 5;
	for(i=0; i<num_samples; i++){
		filter[i] = 1.0;
	}
	for(i=0; i<n_interf; i++){
		pos = freqList[i];
		pos_min = pos - freqWidths[i];
		pos_max = pos + freqWidths[i];
		for(j=pos_min; j<=pos_max; j++){
			filter[j] = 0.0;
		}
	}
	return;
}

// Based on FFT from GSL version. It takes ~35 seconds for 1 second of data at multivac (FFTW-based version takes ~30 second)
/*
void compute_ITFwav_SPIR( int* spirdata, float* filter, double* wav_real, double* wav_imag, int init_sample, int wav_size, double* phases_UP, double* phases_DW, int msec ){
	int i, j, ind;
	double sampling_rate = 80000000.0;
	const int num_samples = int(sampling_rate)/1000; //1 msec
	double interm_freq = 19420000.0; //Must be >= 0.0 and < sampling_rate/2
	const int interm_freq_sample = int(interm_freq)/(int(sampling_rate)/num_samples);
	const int spir_block_ints = 16384;
	double phasors_UP[8][2];
	double phasors_DW[8][2];
	double *wavUP, *wavDW, *wavCROSS;
	gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc(num_samples);
	workspace = gsl_fft_complex_workspace_alloc(num_samples);
	wavUP = (double*) calloc((2*num_samples), sizeof(double));
	wavDW = (double*) calloc((2*num_samples), sizeof(double));
	wavCROSS = (double *) malloc(2*num_samples*sizeof(double));
	//Get signals
	for(i=0; i<8; i++){
		phasors_UP[i][0] = cos(phases_UP[i]*PI_NUM/180.0);
		phasors_UP[i][1] = sin(phases_UP[i]*PI_NUM/180.0);
		phasors_DW[i][0] = cos(phases_DW[i]*PI_NUM/180.0);
		phasors_DW[i][1] = sin(phases_DW[i]*PI_NUM/180.0);
	}
	ind = msec*num_samples + msec*num_samples/spir_block_ints;
	for(i=0; i<num_samples; i++){
		if(ind%spir_block_ints == 0){
			ind ++;
		}
		for(j=0; j<8; j++){
			wavUP[i*2] = wavUP[i*2] + getSPIRdatabit(spirdata[ind], 2*j)*phasors_UP[j][0] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_UP[j][1];
			wavUP[i*2 + 1] = wavUP[i*2 + 1] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_UP[j][0] - getSPIRdatabit(spirdata[ind], 2*j)*phasors_UP[j][1];
			wavDW[i*2] = wavDW[i*2] + getSPIRdatabit(spirdata[ind], (2*j + 16))*phasors_DW[j][0] + getSPIRdatabit(spirdata[ind], (2*j + 17))*phasors_DW[j][1];
			wavDW[i*2 + 1] = wavDW[i*2 + 1] + getSPIRdatabit(spirdata[ind], (2*j + 17))*phasors_DW[j][0] - getSPIRdatabit(spirdata[ind], (2*j + 16))*phasors_DW[j][1];
		}
		ind ++;
	}
	gsl_fft_complex_forward(wavUP, 1, num_samples, wavetable, workspace);
	gsl_fft_complex_forward(wavDW, 1, num_samples, wavetable, workspace);
	//Frequency Domain
	for(i=(num_samples/2); i<(num_samples - interm_freq_sample); i++){
		j = i + interm_freq_sample;
		wavCROSS[2*i] = (wavDW[2*j]*wavUP[2*j] + wavDW[2*j + 1]*wavUP[2*j + 1])*double(filter[i]);
		wavCROSS[2*i + 1] = (wavDW[2*j + 1]*wavUP[2*j] - wavDW[2*j]*wavUP[2*j + 1])*double(filter[i]);
	}
	for(i=(num_samples - interm_freq_sample); i<num_samples; i++){
		j = i - num_samples + interm_freq_sample;
		wavCROSS[2*i] = (wavDW[2*j]*wavUP[2*j] + wavDW[2*j + 1]*wavUP[2*j + 1])*double(filter[i]);
		wavCROSS[2*i + 1] = (wavDW[2*j + 1]*wavUP[2*j] - wavDW[2*j]*wavUP[2*j + 1])*double(filter[i]);
	}
	for(i=0; i<(num_samples/2 - interm_freq_sample); i++){
		j = i + interm_freq_sample;
		wavCROSS[2*i] = (wavDW[2*j]*wavUP[2*j] + wavDW[2*j + 1]*wavUP[2*j + 1])*double(filter[i]);
		wavCROSS[2*i + 1] = (wavDW[2*j + 1]*wavUP[2*j] - wavDW[2*j]*wavUP[2*j + 1])*double(filter[i]);
	}
	for(i=(num_samples/2 - interm_freq_sample); i<(num_samples/2); i++){
		wavCROSS[2*i] = 0.0;
		wavCROSS[2*i + 1] = 0.0;
	}
	free(wavUP);
	free(wavDW);
	gsl_fft_complex_inverse(wavCROSS, 1, num_samples, wavetable, workspace);
	//Time Domain
	for(i=0; i<wav_size; i++){
		wav_real[i] = wavCROSS[2*i + 2*init_sample];
		wav_imag[i] = wavCROSS[2*i + 1 + 2*init_sample];
	}
	free(wavCROSS);
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);
	return;
}
*/
// Based on FFT from FFTW version. It takes ~30 seconds for 1 second of data at multivac.
/*
void compute_ITFwav_SPIR( int* spirdata, double interm_freq, float* filter, double* wav_real, double* wav_imag, int init_sample, int wav_size, double* phases_UP, double* phases_DW, int msec ){
	//interm_freq Must be >= 0.0 and < sampling_rate/2
	int i, j, ind;
	double sampling_rate = 80000000.0;
	const int num_samples = int(sampling_rate)/1000; //1 msec
	const int interm_freq_sample = int(interm_freq)/(int(sampling_rate)/num_samples);
	const int spir_block_ints = 16384;
	double phasors_UP[8][2];
	double phasors_DW[8][2];
	fftw_complex *wavUP_time, *wavUP_freq, *wavDW_time, *wavDW_freq, *wavCROSS_time, *wavCROSS_freq;
	fftw_plan up_fw, dw_fw, cross_bw;
	wavUP_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavUP_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	up_fw = fftw_plan_dft_1d(num_samples, wavUP_time, wavUP_freq, FFTW_FORWARD, FFTW_ESTIMATE);
	wavDW_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavDW_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	dw_fw = fftw_plan_dft_1d(num_samples, wavDW_time, wavDW_freq, FFTW_FORWARD, FFTW_ESTIMATE);
	wavCROSS_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavCROSS_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	cross_bw = fftw_plan_dft_1d(num_samples, wavCROSS_freq, wavCROSS_time, FFTW_BACKWARD, FFTW_ESTIMATE);
	//Get signals
	for(i=0; i<8; i++){
		phasors_UP[i][0] = cos(phases_UP[i]*PI_NUM/180.0);
		phasors_UP[i][1] = sin(phases_UP[i]*PI_NUM/180.0);
		phasors_DW[i][0] = cos(phases_DW[i]*PI_NUM/180.0);
		phasors_DW[i][1] = sin(phases_DW[i]*PI_NUM/180.0);
	}
	ind = msec*num_samples + msec*num_samples/spir_block_ints;
	for(i=0; i<num_samples; i++){
		if(ind%spir_block_ints == 0){ //Beginning of block = ID block
			ind ++;
		}
		wavUP_time[i][0] = 0.0;
		wavUP_time[i][1] = 0.0;
		wavDW_time[i][0] = 0.0;
		wavDW_time[i][1] = 0.0;
		for(j=0; j<8; j++){
			wavUP_time[i][0] = wavUP_time[i][0] + getSPIRdatabit(spirdata[ind], 2*j)*phasors_UP[j][0] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_UP[j][1];
			wavUP_time[i][1] = wavUP_time[i][1] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_UP[j][0] - getSPIRdatabit(spirdata[ind], 2*j)*phasors_UP[j][1];
			wavDW_time[i][0] = wavDW_time[i][0] + getSPIRdatabit(spirdata[ind], (2*j + 16))*phasors_DW[j][0] + getSPIRdatabit(spirdata[ind], (2*j + 17))*phasors_DW[j][1];
			wavDW_time[i][1] = wavDW_time[i][1] + getSPIRdatabit(spirdata[ind], (2*j + 17))*phasors_DW[j][0] - getSPIRdatabit(spirdata[ind], (2*j + 16))*phasors_DW[j][1];
		}
		ind ++;
	}
	fftw_execute(up_fw);
	fftw_execute(dw_fw);
	//Frequency Domain
	for(i=(num_samples/2); i<(num_samples - interm_freq_sample); i++){
		j = i + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavDW_freq[j][0]*wavUP_freq[j][0] + wavDW_freq[j][1]*wavUP_freq[j][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavDW_freq[j][1]*wavUP_freq[j][0] - wavDW_freq[j][0]*wavUP_freq[j][1])*double(filter[i]);
	}
	for(i=(num_samples - interm_freq_sample); i<num_samples; i++){
		j = i - num_samples + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavDW_freq[j][0]*wavUP_freq[j][0] + wavDW_freq[j][1]*wavUP_freq[j][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavDW_freq[j][1]*wavUP_freq[j][0] - wavDW_freq[j][0]*wavUP_freq[j][1])*double(filter[i]);
	}
	for(i=0; i<(num_samples/2 - interm_freq_sample); i++){
		j = i + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavDW_freq[j][0]*wavUP_freq[j][0] + wavDW_freq[j][1]*wavUP_freq[j][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavDW_freq[j][1]*wavUP_freq[j][0] - wavDW_freq[j][0]*wavUP_freq[j][1])*double(filter[i]);
	}
	for(i=(num_samples/2 - interm_freq_sample); i<(num_samples/2); i++){
		wavCROSS_freq[i][0] = 0.0;
		wavCROSS_freq[i][1] = 0.0;
	}
	fftw_destroy_plan(up_fw);
	fftw_destroy_plan(dw_fw);
	fftw_free(wavUP_time);
	fftw_free(wavUP_freq);
	fftw_free(wavDW_time);
	fftw_free(wavDW_freq);
	fftw_execute(cross_bw);
	//Time Domain
	for(i=0; i<wav_size; i++){
		wav_real[i] = wavCROSS_time[i + init_sample][0];
		wav_imag[i] = wavCROSS_time[i + init_sample][1];
	}
	fftw_destroy_plan(cross_bw);
	fftw_free(wavCROSS_time);
	fftw_free(wavCROSS_freq);
	return;
}
*/

void compute_ITFwav_SPIR( int* spirdata, double interm_freq, float* filter, double* wav_real, double* wav_imag, int init_sample, int wav_size, double* phases_BF, bool* elementsUP, bool* elementsDW, int msec ){
	//interm_freq Must be >= 0.0 and < sampling_rate/2
	int i, j, ind;
	double sampling_rate = 80000000.0;
	const int num_samples = int(sampling_rate)/1000; //1 msec
	const int interm_freq_sample = int(interm_freq)/(int(sampling_rate)/num_samples);
	const int spir_block_ints = 16384;
	double phasors_BF[16][2];
	fftw_complex *wavUP_time, *wavUP_freq, *wavDW_time, *wavDW_freq, *wavCROSS_time, *wavCROSS_freq;
	fftw_plan up_fw, dw_fw, cross_bw;
	wavUP_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavUP_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	up_fw = fftw_plan_dft_1d(num_samples, wavUP_time, wavUP_freq, FFTW_FORWARD, FFTW_ESTIMATE);
	wavDW_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavDW_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	dw_fw = fftw_plan_dft_1d(num_samples, wavDW_time, wavDW_freq, FFTW_FORWARD, FFTW_ESTIMATE);
	wavCROSS_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavCROSS_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	cross_bw = fftw_plan_dft_1d(num_samples, wavCROSS_freq, wavCROSS_time, FFTW_BACKWARD, FFTW_ESTIMATE);
	//Get signals
	for(i=0; i<16; i++){
		phasors_BF[i][0] = cos(phases_BF[i]*PI_NUM/180.0);
		phasors_BF[i][1] = sin(phases_BF[i]*PI_NUM/180.0);
	}
	ind = msec*num_samples + msec*num_samples/spir_block_ints;
	for(i=0; i<num_samples; i++){
		if(ind%spir_block_ints == 0){ //Beginning of block = ID block
			ind ++;
		}
		wavUP_time[i][0] = 0.0;
		wavUP_time[i][1] = 0.0;
		wavDW_time[i][0] = 0.0;
		wavDW_time[i][1] = 0.0;
		for(j=0; j<16; j++){
			if(elementsUP[j]){
				wavUP_time[i][0] = wavUP_time[i][0] + getSPIRdatabit(spirdata[ind], 2*j)*phasors_BF[j][0] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_BF[j][1];
				wavUP_time[i][1] = wavUP_time[i][1] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_BF[j][0] - getSPIRdatabit(spirdata[ind], 2*j)*phasors_BF[j][1];
			}
			if(elementsDW[j]){
				wavDW_time[i][0] = wavDW_time[i][0] + getSPIRdatabit(spirdata[ind], 2*j)*phasors_BF[j][0] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_BF[j][1];
				wavDW_time[i][1] = wavDW_time[i][1] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_BF[j][0] - getSPIRdatabit(spirdata[ind], 2*j)*phasors_BF[j][1];
			}
		}
		ind ++;
	}
	fftw_execute(up_fw);
	fftw_execute(dw_fw);
	//Frequency Domain
	for(i=(num_samples/2); i<(num_samples - interm_freq_sample); i++){
		j = i + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavDW_freq[j][0]*wavUP_freq[j][0] + wavDW_freq[j][1]*wavUP_freq[j][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavDW_freq[j][1]*wavUP_freq[j][0] - wavDW_freq[j][0]*wavUP_freq[j][1])*double(filter[i]);
	}
	for(i=(num_samples - interm_freq_sample); i<num_samples; i++){
		j = i - num_samples + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavDW_freq[j][0]*wavUP_freq[j][0] + wavDW_freq[j][1]*wavUP_freq[j][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavDW_freq[j][1]*wavUP_freq[j][0] - wavDW_freq[j][0]*wavUP_freq[j][1])*double(filter[i]);
	}
	for(i=0; i<(num_samples/2 - interm_freq_sample); i++){
		j = i + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavDW_freq[j][0]*wavUP_freq[j][0] + wavDW_freq[j][1]*wavUP_freq[j][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavDW_freq[j][1]*wavUP_freq[j][0] - wavDW_freq[j][0]*wavUP_freq[j][1])*double(filter[i]);
	}
	for(i=(num_samples/2 - interm_freq_sample); i<(num_samples/2); i++){
		wavCROSS_freq[i][0] = 0.0;
		wavCROSS_freq[i][1] = 0.0;
	}
	fftw_destroy_plan(up_fw);
	fftw_destroy_plan(dw_fw);
	fftw_free(wavUP_time);
	fftw_free(wavUP_freq);
	fftw_free(wavDW_time);
	fftw_free(wavDW_freq);
	fftw_execute(cross_bw);
	//Time Domain
	for(i=0; i<wav_size; i++){
		//wav_real[i] = wavCROSS_time[i + init_sample][0];
		//wav_imag[i] = wavCROSS_time[i + init_sample][1];
		wav_real[i] = wavCROSS_time[((i + init_sample)%num_samples + num_samples)%num_samples][0];
		wav_imag[i] = wavCROSS_time[((i + init_sample)%num_samples + num_samples)%num_samples][1];
	}
	fftw_destroy_plan(cross_bw);
	fftw_free(wavCROSS_time);
	fftw_free(wavCROSS_freq);
	return;
}

bool compute_CRwav_SPIR( int* spirdata, double interm_freq, float* model_CR, double doppler_CR, int samples_CR, float* filter, double* wav_real, double* wav_imag, int wav_size, double* phases_BF, bool* elements_CORR, int msec, int &init_sample, bool compute_init_sample ){
	//interm_freq Must be >= 0.0 and < sampling_rate/2
	int i, j, ind, ref_zero_sample;
	bool valid_init_sample;
	double max_val;
	double sampling_rate = 80000000.0;
	const int num_samples = int(sampling_rate)/1000; //1 msec
	const int interm_freq_sample = int(interm_freq)/(int(sampling_rate)/num_samples);
	const int spir_block_ints = 16384;
	double phasors_BF[16][2];
	fftw_complex *wavModel_time, *wavModel_freq, *wavData_time, *wavData_freq, *wavCROSS_time, *wavCROSS_freq;
	fftw_plan model_fw, data_fw, cross_bw;
	wavModel_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavModel_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	model_fw = fftw_plan_dft_1d(num_samples, wavModel_time, wavModel_freq, FFTW_FORWARD, FFTW_ESTIMATE);
	wavData_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavData_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	data_fw = fftw_plan_dft_1d(num_samples, wavData_time, wavData_freq, FFTW_FORWARD, FFTW_ESTIMATE);
	wavCROSS_time = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	wavCROSS_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_samples);
	cross_bw = fftw_plan_dft_1d(num_samples, wavCROSS_freq, wavCROSS_time, FFTW_BACKWARD, FFTW_ESTIMATE);
	//Get signals
	for(i=0; i<16; i++){
		phasors_BF[i][0] = cos(phases_BF[i]*PI_NUM/180.0);
		phasors_BF[i][1] = sin(phases_BF[i]*PI_NUM/180.0);
	}
	ind = msec*num_samples + msec*num_samples/spir_block_ints;
	for(i=0; i<num_samples; i++){
		if(ind%spir_block_ints == 0){ //Beginning of block = ID block
			ind ++;
		}
		wavModel_time[i][0] = model_CR[(i*samples_CR/num_samples)%samples_CR]*cos(doppler_CR*(double(i)/sampling_rate)*2.0*PI_NUM);
		wavModel_time[i][1] = model_CR[(i*samples_CR/num_samples)%samples_CR]*sin(doppler_CR*(double(i)/sampling_rate)*2.0*PI_NUM);
		//
		wavData_time[i][0] = 0.0;
		wavData_time[i][1] = 0.0;
		for(j=0; j<16; j++){
			if(elements_CORR[j]){
				wavData_time[i][0] = wavData_time[i][0] + getSPIRdatabit(spirdata[ind], 2*j)*phasors_BF[j][0] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_BF[j][1];
				wavData_time[i][1] = wavData_time[i][1] + getSPIRdatabit(spirdata[ind], (2*j + 1))*phasors_BF[j][0] - getSPIRdatabit(spirdata[ind], 2*j)*phasors_BF[j][1];
			}
		}
		ind ++;
	}
	fftw_execute(model_fw);
	fftw_execute(data_fw);
	//Frequency Domain
	for(i=(num_samples/2); i<(num_samples - interm_freq_sample); i++){
		j = i + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavData_freq[j][0]*wavModel_freq[i][0] + wavData_freq[j][1]*wavModel_freq[i][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavData_freq[j][1]*wavModel_freq[i][0] - wavData_freq[j][0]*wavModel_freq[i][1])*double(filter[i]);
	}
	for(i=(num_samples - interm_freq_sample); i<num_samples; i++){
		j = i - num_samples + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavData_freq[j][0]*wavModel_freq[i][0] + wavData_freq[j][1]*wavModel_freq[i][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavData_freq[j][1]*wavModel_freq[i][0] - wavData_freq[j][0]*wavModel_freq[i][1])*double(filter[i]);
	}
	for(i=0; i<(num_samples/2 - interm_freq_sample); i++){
		j = i + interm_freq_sample;
		wavCROSS_freq[i][0] = (wavData_freq[j][0]*wavModel_freq[i][0] + wavData_freq[j][1]*wavModel_freq[i][1])*double(filter[i]);
		wavCROSS_freq[i][1] = (wavData_freq[j][1]*wavModel_freq[i][0] - wavData_freq[j][0]*wavModel_freq[i][1])*double(filter[i]);
	}
	for(i=(num_samples/2 - interm_freq_sample); i<(num_samples/2); i++){
		wavCROSS_freq[i][0] = 0.0;
		wavCROSS_freq[i][1] = 0.0;
	}
	fftw_destroy_plan(model_fw);
	fftw_destroy_plan(data_fw);
	fftw_free(wavModel_time);
	fftw_free(wavModel_freq);
	fftw_free(wavData_time);
	fftw_free(wavData_freq);
	fftw_execute(cross_bw);
	//Time Domain
	valid_init_sample = true;
	if(compute_init_sample){
		max_val = wavCROSS_time[0][0]*wavCROSS_time[0][0] + wavCROSS_time[0][1]*wavCROSS_time[0][1];
		ref_zero_sample = 0;
		for(i=1; i<num_samples; i++){
			if((wavCROSS_time[i][0]*wavCROSS_time[i][0] + wavCROSS_time[i][1]*wavCROSS_time[i][1]) > max_val){
				ref_zero_sample = i;
				max_val = wavCROSS_time[i][0]*wavCROSS_time[i][0] + wavCROSS_time[i][1]*wavCROSS_time[i][1];
			}
		}
		if(max_val/(wavCROSS_time[(ref_zero_sample + 80)%num_samples][0]*wavCROSS_time[(ref_zero_sample + 80)%num_samples][0] + wavCROSS_time[(ref_zero_sample + 80)%num_samples][1]*wavCROSS_time[(ref_zero_sample + 80)%num_samples][1]) < 2.0){
			valid_init_sample = false;
			init_sample = 0;
		}else{
			init_sample = ref_zero_sample - wav_size/2;
		}
	}
	for(i=0; i<wav_size; i++){
		wav_real[i] = wavCROSS_time[((i + init_sample)%num_samples + num_samples)%num_samples][0];
		wav_imag[i] = wavCROSS_time[((i + init_sample)%num_samples + num_samples)%num_samples][1];
	}
	fftw_destroy_plan(cross_bw);
	fftw_free(wavCROSS_time);
	fftw_free(wavCROSS_freq);
	return valid_init_sample;
}

double getSPIRdatabit( int word, int pos ){
	int val;
	val = (word>>pos)&0x00000001;
	if(val == 1){
		return 1.0;
	}else{
		return -1.0;
	}
}

void get_doppler_initsamp_CR_SPIR( int* spirdata, double interm_freq, float* model_CR, double doppler_CR, int samples_CR, float* filter, int wav_size, double* phases_BF, double &doppler_peak, int &init_sample, bool selected_signals, bool elements_CORR[16], bool searching ){
	int i, doppler_samples, msec_ref, iterations;
	double doppler_res, max_pow;
	double *wav_i, *wav_q, *wav_pow, *dop_fit, *pow_dop_fit;
	bool elements_CORR_UP[16];
	Waveform_power wavpow;
	wav_i = (double*) malloc (wav_size*sizeof(double));
	wav_q = (double*) malloc (wav_size*sizeof(double));
	wav_pow = (double*) malloc (wav_size*sizeof(double));
    if(selected_signals){
        for(i=0; i<16; i++){
            elements_CORR_UP[i] = elements_CORR[i];
        }
    }else{
        for(i=0; i<16; i++){
            if(i<8){
                elements_CORR_UP[i] = true;
            }else{
                elements_CORR_UP[i] = false;
            }
        }
    }
	//Coarse acquisition
	iterations = 0;
	msec_ref = 500;
	while((!compute_CRwav_SPIR(spirdata, interm_freq, model_CR, doppler_CR, samples_CR, filter, wav_i, wav_q, wav_size, phases_BF, elements_CORR_UP, msec_ref, init_sample, true))&&(iterations < 200)){
		msec_ref = (msec_ref + 5)%998;
		iterations ++;
	}
	wavpow.set_amp_waveform(wav_i, wav_size, wav_q, wav_size);
	wavpow.compute_delays();
	max_pow = wavpow.powerMax;
	doppler_peak = doppler_CR;
	//
    if(searching){
        doppler_samples = 100;
        doppler_res = 100.0;
    }else{
        doppler_samples = 20;
        doppler_res = 200.0;
    }
	for(i=0; i<doppler_samples; i++){
		if(i != doppler_samples/2){ //Already done
			compute_CRwav_SPIR(spirdata, interm_freq, model_CR, doppler_CR + double(i - doppler_samples/2)*doppler_res, samples_CR, filter, wav_i, wav_q, wav_size, phases_BF, elements_CORR_UP, msec_ref, init_sample, true);
			wavpow.set_amp_waveform(wav_i, wav_size, wav_q, wav_size);
			wavpow.compute_delays();
			if(wavpow.powerMax > max_pow){
				max_pow = wavpow.powerMax;
				doppler_peak = doppler_CR + double(i - doppler_samples/2)*doppler_res;
			}
		}
	}
	//Fine acquisition
	doppler_samples = 20;
	doppler_res = 25.0;
	dop_fit = (double*) malloc (doppler_samples*sizeof(double));
	pow_dop_fit = (double*) malloc (doppler_samples*sizeof(double));
	dop_fit[doppler_samples/2] = doppler_peak;
	pow_dop_fit[doppler_samples/2] = max_pow;
	for(i=0; i<doppler_samples; i++){
		if(i != doppler_samples/2){ //Already done
			compute_CRwav_SPIR(spirdata, interm_freq, model_CR, doppler_peak + double(i - doppler_samples/2)*doppler_res, samples_CR, filter, wav_i, wav_q, wav_size, phases_BF, elements_CORR_UP, msec_ref, init_sample, true);
			wavpow.set_amp_waveform(wav_i, wav_size, wav_q, wav_size);
			wavpow.compute_delays();
			dop_fit[i] = doppler_peak + double(i - doppler_samples/2)*doppler_res;
			pow_dop_fit[i] = wavpow.powerMax;
		}
	}
	//Polynomial fit to obtain Doppler max
	double chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *c;
	X = gsl_matrix_alloc (doppler_samples, 3);
	y = gsl_vector_alloc (doppler_samples);
	c = gsl_vector_alloc (3);
	cov = gsl_matrix_alloc (3, 3);
	for(i=0; i<doppler_samples; i++){
		gsl_matrix_set (X, i, 0, 1.0);
		gsl_matrix_set (X, i, 1, dop_fit[i]);
		gsl_matrix_set (X, i, 2, dop_fit[i]*dop_fit[i]);
		gsl_vector_set (y, i, pow_dop_fit[i]);
	}
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(doppler_samples, 3);
	gsl_multifit_linear(X, y, c, cov, &chisq, work);
	gsl_multifit_linear_free(work);
	doppler_peak = -gsl_vector_get(c,(1))/(2.0*gsl_vector_get(c,(2)));
	compute_CRwav_SPIR(spirdata, interm_freq, model_CR, doppler_peak, samples_CR, filter, wav_i, wav_q, wav_size, phases_BF, elements_CORR_UP, msec_ref, init_sample, true);
    if(searching){
        wavpow.set_amp_waveform(wav_i, wav_size, wav_q, wav_size);
        wavpow.compute_delays();
        printf(" Doppler %5.1f Hz, 1-msec PNR %3.1f dB\n", doppler_peak, 10.0*log10(wavpow.powerMax/wavpow.floorNoise));
    }
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
	free(wav_i);
	free(wav_q);
	free(wav_pow);
	free(dop_fit);
	free(pow_dop_fit);
	return;
}

void read_PRN_code( int code_ref, float* code_CR, int code_length ){
    std::string data_path = LOCAL_DATA_PATH;
	std::string codefilepath, line, code_ref_str;
	std::stringstream number;
	number << code_ref;
        number >> code_ref_str;
	if(code_ref_str.length()==1){
		codefilepath = data_path + "GPS_L1_PRN0" + code_ref_str + ".txt";
	}else{
		codefilepath = data_path + "GPS_L1_PRN" + code_ref_str + ".txt";
	}
	std::ifstream codefile(codefilepath.c_str());
	if(!codefile){
		std::cout << "ERROR! Unable to read code file: " << codefilepath << std::endl;
		return;
	}
	getline(codefile,line);
	int i;
	for(i=0; i<code_length; i++){
		if(line[i] == '1'){
			code_CR[i] = 1.0;
		}else{
			code_CR[i] = -1.0;
		}
	}
	codefile.close();
    return;
}


