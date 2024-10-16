/**
    \file compute_delays_wav.h
    \brief TODO
*/
void compute_delays_wav( double* waveform_in, int wav_length, double sampling_rate, double relFactor, double &positionMax, double &posSampleMax, double &sigma_posMax, double &power_posMax, double &positionDer, double &posSampleDer, double &powSampleDer, double &sigma_posDer, double &power_posDer, double &powerDer_posDer, double &positionRel, double &posSampleRel, double &sigma_posRel, double &noiseFloor, double &slope_NormTail, double &sigma_slope_NormTail, double normTail_length, double min_resolution_fft_interp, double fit_length, bool apply_speckle_weights, int num_incoh, bool apply_limits, double limits_center, double limits_width, double apriori_scattdel, double bandwidth_LPF, int tail_factor, int noise_lags, int tail_lags, bool get_deriv_tracks, double* tau_deriv, double* wav_track_deriv, double* wav_first_deriv, double* wav_second_deriv );
