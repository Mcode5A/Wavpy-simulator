/**
  \file waveform_power.h
  \brief Header of the Waveform_power class.
*/

/**
  \brief This class provides means to analyze and characterize power (real values) GNSS+R waveforms.
 
  TODO: Add a definition/formula of what is a power waveform.
*/

class Waveform_power
{
  private:
	double *waveform;
	int wav_length, tail_factor, noise_lags, tail_lags;
	double sampling_rate, scale_factor, rel_factor, init_range, min_resolution_fft_interp, fit_length, normTail_length;
	void scale_delays( void );
  public:
    double positionMax;          /**< \brief Range position of waveform's peak (computed by means of a linear fit at the first waveform's derivative) in meters. Default value: 0.0 */
    double posSampleMax;         /**< \brief Range position of waveform's peak (computed by taking the maximum sample value of the interpolated waveform) in meters. Default value: 0.0 */
    double sigma_posMax;         /**< \brief Formal standard deviation of the position of the waveform's peak (as a result of the linear fir applied) in meters. Default value: -1.0 */
    double powerMax;             /**< \brief Waveform's peak power (computed by taking the maximum sample value of the interpolated waveform) in the stored waveform's units. Default value: 0.0 */
    double positionDer;          /**< \brief Range position of waveform's peak derivative (computed by means of a linear fit at the second waveform's derivative) in meters. Default value: 0.0 */
    double posSampleDer;         /**< \brief Range position of waveform's peak derivative (computed by taking the maximum sample value of the interpolated waveform's first derivative) in meters. Default value: 0.0 */
    double powSampleDer;         /**< \brief Formal standard deviation of the position of the waveform's peak derivative (as a result of the linear fir applied) in meters. Default value: -1.0 */
    double sigma_posDer;         /**< \brief Waveform's power at position of maximum's first derivative (computed by taking the corresponding sample value of the interpolated waveform) in the stored waveform's units. Default value: 0.0 */
    double power_posDer;         /**< \brief Waveform's first derivative power at position of maximum's first derivative (computed by taking the corresponding sample value of the interpolated waveform) in the stored waveform's units divided by meters. Default value: 0.0 */
    double powerDer_posDer;      /**< \brief Waveform's noise level (computed by averaging the first samples) in the stored waveform's units. Default value: 0.0 */
    double floorNoise;           /**< \brief Range position of waveform's peak multiplied by __rel_factor__ at the leading edge (computed by means of a linear fit at the waveform) in meters. Default value: 0.0 */
    double positionRel;          /**< \brief Range position of waveform's peak multiplied by __rel_factor__ at the leading edge (computed by taking the corresponding sample value of the interpolated waveform) in meters. Default value: 0.0 */
    double posSampleRel;         /**< \brief Formal standard deviation of the position of waveform's peak multiplied by __rel_factor__ at the leading edge (as a result of the linear fir applied) in meters. Default value: -1.0 */
    double sigma_posRel;         /**< \brief Slope of the trailing edge of the normalized waveform (computed by means of a linear fit at the waveform) in meters$^{-1}$. Default value: 0.0 */
    double slope_normTail;       /**< \brief Formal standard deviation of the slope of the trailing edge of the normalized waveform (as a result of the linear fir applied) in meters$^{-1}$. Default value: -1.0 */
    double sigma_slope_normTail; /**< \brief Formal standard deviation of the slope of the trailing edge of the normalized waveform (as a result of the linear fir applied) in meters$^{-1}$. Default value: -1.0 */
	// Methods
    /** \brief Constructor
    */
    Waveform_power( void );

	/** \brief Set a power waveform

    	@param waveform_in        Power waveform in arbitrary units
	    @param wav_in_length      Number of samples of the input power waveform
	*/
	void set_waveform( double* waveform_in, int wav_in_length );

	/** \brief Set a power waveform of type float

    	@param float_waveform_in    Power waveform in arbitrary units
    	@param wav_in_length        Number of samples of the input power waveform
	*/
	void set_float_waveform( float* float_waveform_in, int wav_in_length );

	/** \brief Set a normalized power waveform of type float

    	@param norm_waveform_in     Normalized power waveform
    	@param norm_wav_in_length   Number of samples of the input power waveform
	    @param scale_factor_in      Maximum value of the corresponding de-normalized waveform in arbitrary units
	*/
	void set_norm_waveform( float* norm_waveform_in, int norm_wav_in_length, double scale_factor_in );

	/** \brief 	Set a power waveform from in-phase and quadrature amplitude components

        @param wav_i_in         Amplitude waveform (in-phase component) in arbitrary units
        @param wav_i_in_length  Number of samples of the input amplitude waveform
        @param wav_q_in         Amplitude waveform (quadrature component) in arbitrary units
        @param wav_q_in_length  Number of samples of the input amplitude waveform
	*/
	void set_amp_waveform( double* wav_i_in, int wav_i_in_length, double* wav_q_in, int wav_q_in_length );

	/** \brief Get the stored power waveform

        @param waveform_out     Power waveform in arbitrary units
        @param wav_out_length   Number of samples of the ouput power waveform
	*/
	void get_waveform( double* waveform_out, int wav_out_length );

	/** \brief Update the stored waveform by weighty averaging its values with an input waveform after being re-tracked a given delay (by means of FFT). The resultant stored waveform will be given by:

    	\f[
             new_stored_waveform[range] = (1.0 - wav_weight) * waveform[range] + wav_weight * waveform_in[range - retrack_delay]
    	\f]
	    @param waveform_in             Normalized power waveform
    	@param wav_in_length           Number of samples of the input power waveform
	    @param retrack_delay           Retracking delay applied to waveform_in before being averaged with the stored waveform in meters.
    	@param wav_weight              Weight applied to waveform_in for its averaging with the stored waveform
	    @param apply_safety_margin     Boolean that controls the application of a safety margin during the re-tracking
	*/
	void add_waveform_retracking( double* waveform_in, int wav_in_length, double retrack_delay, float wav_weight, bool apply_safety_margin );

	/** \brief Set the sampling rate of the stored waveform

    	@param samplingRate     Sampling rate of the waveform in samples/sec
	*/
	void set_sampling_rate( double samplingRate );

	/** \brief 	Set the scaling factor used for the computation of the relative delay located at the leading edge with a power equal to the product of such factor with the peak of the waveform

    	@param relFactor   Scaling factor used for the computation of the relative delay located at the leading edge with a power equal to the product of such factor with the peak of the waveform
	*/
	void set_rel_factor( double relFactor );

	/** \brief 	Get the scaling factor used for the computation of the relative delay located at the leading edge with a power equal to the product of such factor with the peak of the waveform

    	@return  Scaling factor used for the computation of the relative delay located at the leading edge with a power equal to the product of such factor with the peak of the waveform
	*/
	double get_rel_factor( void );

	/** \brief Set the minimum range resolution of the FFT interpolation for the computation of the delays

    	@param resolution_in     Minimum range resolution of the FFT interpolation for the computation of the delays in meters
	*/
	void set_min_resolution_fft_interp( double resolution_in );

	/** \brief Set the range length of the segment employed for the linear fitting when computing the delays

    	@param fit_length_in     Range length of the segment employed for the linear fitting when computing the delays in meters
	*/
	void set_fit_length( double fit_length_in );

	/** \brief 	Set range length of the segment employed for the linear fitting when computing the slope of the trailing edge in meters

    	@param normtail_length_in     Range length of the segment employed for the linear fitting when computing the slope of the trailing edge in meters
	*/
	void set_normtail_length( double normtail_length_in );

	/** \brief Set factor to increase waveform's length with the purpose of connecting both ends during computation of delays

    	@param tail_factor_in     Integer factor (to be multiplied by wav_length) used to compute the extended length of the waveform for connecting both ends before FFT-based interpolation.
	*/
	void set_tail_factor( int tail_factor_in );
    
    /** \brief Set the number of lags employed to compute the noise floor

        @param noise_lags_in     Number of lags, starting from zero, employed to compute the noise floor
	*/
	void set_noise_lags( int noise_lags_in );
    
    /** \brief 	Set the number of lags employed to extrapolate the tail of the waveform until the noise floor

        @param tail_lags_in     Number of lags employed to extrapolate the tail of the waveform until the noise floor
	*/
	void set_tail_lags( int tail_lags_in );
    
    /** \brief Compute relevant delay positions and power levels along the waveform
	*/
	void compute_delays( void );
    
    /** \brief Compute relevant delay positions and power lever along the waveform by assuming speckle noise distribution for the estimation of the standard deviations

        @param num_incoh     Number of independent waveforms employed during the incoherent integration
	*/
	void compute_delays_wspeckle( int num_incoh );
    
    /** \brief Compute relevant delay positions and power levels around a limited interval of the waveform

        @param limits_center      Center location of the range interval to constraint the initial estimation of the specular delay
        @param limits_width       Range interval around limits_center to constraint the initial estimation of the specular delay
        @param apriori_scattdel   A-priori estimation of the scatterometric delay (delay distance from specular to peak of the waveform)
	*/
	void compute_delays_wlimits( double limits_center, double limits_width, double apriori_scattdel );
    
    /** \brief Compute relevant delay positions and power levels around a limited interval of the waveform after applying a box low-pass filter

        @param limits_center        Center location of the range interval to constraint the initial estimation of the specular delay
	    @param limits_width         Range interval around limits_center to constraint the initial estimation of the specular delay
        @param apriori_scattdel     A-priori estimation of the scatterometric delay (delay distance from specular to peak of the waveform)
        @param bandwidth_LPF_Hz     Bandwidth (in Hz units) of the low-pass filter applied
	*/
	void compute_delays_wlimits_LPF( double limits_center, double limits_width, double apriori_scattdel, double bandwidth_LPF_Hz );
    
    /** \brief 	Set initial range value

    	@param init_range_in     Initial range value
	*/
	void set_init_range( double init_range_in );
    
    /** \brief Get range values for each lag of the waveform

        @param range        Range values for each lag of the waveform
        @param size_range   Size of the waveform
	*/
	void get_range_waveform( double* range, int size_range );
    
    /** \brief Get size (number of samples) of the intervals where the fits are applied on first and second waveform's derivatives 
	*/
	int get_size_deriv_waveform_tracks( void );
    
    /** \brief 	Get the intervals where the fits are applied (around specular delay) from waveform and its first and second derivatives after computing relevant delay positions and power levels

        @param tau_der         Range of the waveform and its first and second derivatives for the interval around the specular delay
        @param size_tau_der    Number of samples of range interval
        @param wav_der         Waveform interval around the specular delay where the fits are applied
        @param size_wav_der    Number of samples of waveform interval
        @param wav_der1        Waveform's first derivative interval around the specular delay where the fits are applied
        @param size_wav_der1   Number of samples of waveform's first derivative interval
        @param wav_der2        Waveform's second derivative interval around the specular delay where the fits are applied
        @param size_wav_der2   Number of samples of waveform's second derivative interval
	*/
	void get_deriv_waveform_tracks( double* tau_der, int size_tau_der, double* wav_der, int size_wav_der, double* wav_der1, int size_wav_der1, double* wav_der2, int size_wav_der2 );
    
    /** \brief Get the intervals where the fits are applied (around specular delay) from waveform and its first and second derivatives after computing relevant delay positions and power levels after applying limits

        @param limits_center      Center location of the range interval to constraint the initial estimation of the specular delay
        @param limits_width       Range interval around limits_center to constraint the initial estimation of the specular delay
        @param apriori_scattdel   A-priori estimation of the scatterometric delay (delay distance from specular to peak of the waveform)
        @param tau_der            Range of the waveform and its first and second derivatives for the interval around the specular delay
        @param size_tau_der       Number of samples of range interval
        @param wav_der            Waveform interval around the specular delay where the fits are applied
        @param size_wav_der       Number of samples of waveform interval
        @param wav_der1           Waveform's first derivative interval around the specular delay where the fits are applied
        @param size_wav_der1      Number of samples of waveform's first derivative interval
        @param wav_der2           Waveform's second derivative interval around the specular delay where the fits are applied
        @param size_wav_der2      Number of samples of waveform's second derivative interval
	*/
	void get_deriv_waveform_tracks_wlimits( double limits_center, double limits_width, double apriori_scattdel, double* tau_der, int size_tau_der, double* wav_der, int size_wav_der, double* wav_der1, int size_wav_der1, double* wav_der2, int size_wav_der2 );
    
    /** \brief 	Dump normalized waveform
	*/
	void dump_norm_waveform( void );
    
    /** \brief 	Dump relevant delay positions and power levels from the waveform (they have to be computed first)
	*/
	void dump_delays( void );
    
    /** \brief 	Get size (number of samples) of the waveform
	*/
	int get_wav_length( void );

  private:
    double* allocate_der();
};
