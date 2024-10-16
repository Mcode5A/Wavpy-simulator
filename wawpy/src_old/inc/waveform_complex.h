/**
    \file waveform_complex.h
    \brief Header file of the Waveform_complex_cluster class.
*/

/**
   \brief This class provides means to analyze and characterize time series of complex GNSS+R waveforms.

*/ 
class Waveform_complex_cluster
{
	double **Icomponents;   /**< \brief In-phase components of waveforms stored in the cluster. Array of `cluster_length x wav_length`. Default value: void  */ 
	double **Qcomponents;   /**< \brief Quadrature components of waveforms stored in the cluster. Array of `cluster_length x wav_length`. Default value: void  */ 
	bool *valid_wavs;       /**< \brief Boolean array indicating the validity of the stored complex waveforms in the cluster. Array of length `cluster_length` elements. Default value: void  */ 
	double *phasorI;        /**< \brief In-phase component of a complex phasor stored in the cluster. Array of length `cluster_length` elements. Default value: void  */ 
	double *phasorQ;        /**< \brief Quadrature component of a complex phasor stored in the cluster. Array of length `cluster_length` elements. Default value: void  */
	bool *valid_phasor;     /**< \brief Boolean array indicating the validity of the stored complex phasor in the cluster. Array of length `cluster_length` elements. Default value: void  */ 
	int wav_length;         /**< \brief Size of the waveforms stored in the cluster. Default value: 0  */ 
	int cluster_length;     /**< \brief Length of the cluster. Default value: 0  */ 

    /** \brief TODO

        \param cluster_in_L  TODO
        \param wav_in_L      TODO
    */
	void prepare_arrays( int cluster_in_L, int wav_in_L );

    /** \brief TODO

        \param namefile                TODO
        \param peak_delay_estimate     TODO
        \param phases_BF               TODO
        \param elements_UP             TODO
        \param elements_DW             TODO
        \param filter_num              TODO
    */
	double compute_ITF_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double phases_BF[16], bool elements_UP[16], bool elements_DW[16], int filter_num );

	/** \brief TODO

        \param namefile              TODO
        \param peak_delay_estimate   TODO
        \param doppler_estimate      TODO
        \param delta_doppler         TODO
        \param phases_BF             TODO
        \param elements_CORR         TODO
        \param code_ref              TODO
        \param selected_signals      TODO
        \param searching             TODO
        \param assist_doppler        TODO
    */
    double compute_CR_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double doppler_estimate, double delta_doppler, double phases_BF[16], bool elements_CORR[16], int code_ref, bool selected_signals, bool searching, bool assist_doppler );
  public:
	int num_valid_wavs;    /**< \brief Number of valid waveforms stored in the cluster. Default value: 0  */ 
	int num_phasor_iter;   /**< \brief Number of operations made with the stored phasor. Default value: 0  */ 
    // METHODS
    /** \brief Constructor */
	Waveform_complex_cluster( void );

    /** \brief Initialize the size of the cluster to allocate the required memory.

        \param in_cluster_length   Length of the cluster.                          
        \param wav_in_length       Length of the waveforms to be stored in the cluster.
    */
	void initialize( int in_cluster_length, int wav_in_length );

    /** \brief Add a complex waveform to the cluster at a given position.

        \param Icomp_in     In-phase components of complex waveform in arbitrary units.
        \param len_I        Number of elements of the `Icomp_in` array
        \param Qcomp_in     Quadrature components of complex waveform in arbitrary units. 
        \param len_Q        Number of elements of the `Qcomp_in` array
        \param cluster_pos  Position to store the complex waveform inside the cluster (from 0 to `cluster_length - 1`).
    */
	void add_waveform( double* Icomp_in, int len_I, double* Qcomp_in, int len_Q, int cluster_pos );

    /** \brief Add a complex waveform to the cluster at a given position and multiplied by a given scaling factor. 

        \param Icomp_in       In-phase components of complex waveform in arbitrary units.    
        \param len_I          Number of elements of the `Icomp_in` array
        \param Qcomp_in       Quadrature components of complex waveform in arbitrary units. 
        \param len_Q          Number of elements of the `Qcomp_in` array
        \param cluster_pos    Position to store the complex waveform inside the cluster (from 0 to `cluster_length - 1`).
        \param scale_factor   Scaling factor applied to the input complex waveform before being stored in the cluster.
    */
	void add_waveform_scale( double* Icomp_in, int len_I, double* Qcomp_in, int len_Q, int cluster_pos, double scale_factor );

    /** \brief Add a complex waveform with the format of IEEC's GOLD-RTR dataset to the cluster at a given position.

        \param Icomponents_in     In-phase components of complex waveform in arbitrary units.
        \param wav_in_length_I    Number of elements in `Icomponents_in` array. 
        \param Qcomponents_in     Quadrature components of complex waveform in arbitrary units.
        \param wav_in_length_Q    Number of elements in `Qcomponents_in` array. Must be the same as `wav_in_length_I`
        \param cluster_pos        Position to store the complex waveform inside the cluster (from 0 to `cluster_length - 1`).
    */
	void add_waveform_GOLD( signed char* Icomponents_in, int wav_in_length_I, signed char* Qcomponents_in, int wav_in_length_Q, int cluster_pos );

    /** \brief Add a complex waveform with the format of IEEC's PIR dataset to the cluster at a given position.

        \param XiYi               Correlation between in-phase components of direct and reflected signals in arbitrary units.                              
        \param wav_in_length_II   Number of elements of the `XiYi` array.
        \param XqYq               Correlation between quadrature components of direct and reflected signals in arbitrary units.
        \param wav_in_length_QQ   Number of elements of the `XqYq` array. 
        \param XiYq               Correlation between in-phase component of direct signal and quadrature component of reflected signal in arbitrary units.
        \param wav_in_length_IQ   Number of elements of the `XiYq` array.
        \param XqYi               Correlation between quadrature component of direct signal and in-phase component of reflected signal in arbitrary units.
        \param wav_in_length_QI   Number of elements of the `XqYi` array.
        \param cluster_pos        Position to store the complex waveform inside the cluster (from 0 to `cluster_length - 1`).
    */
	void add_waveform_PIR( short* XiYi, int wav_in_length_II, short* XqYq, int wav_in_length_QQ, short* XiYq, int wav_in_length_IQ, short* XqYi, int wav_in_length_QI, int cluster_pos );

    /** \brief Load a cluster of 999 complex interferometric waveforms from a 1-second IEEC's SPIR binary file.

        \param namefile             Filename of SPIR 1-second binary file. 
        \param peak_delay_estimate  Raw estimation of the range delay of the waveform's peak in meters.
        \param BF_phases_UP         Beamformer phases (in degrees) to be applied at each of the 8 up-looking antenna elements of the SPIR setup to get maximum antenna gain towards a desired direction.
        \param BF_phases_DW         Beamformer phases (in degrees) to be applied at each of the 8 down-looking antenna elements of the SPIR setup to get maximum antenna gain towards a desired direction.
        \param filter_num           Code to select the frequency-domain filter to be applied during the processing.
                                    Value | Meaning
                                    :----:|:-------
                                      1   | Galileo at E1
                                      2   | GPS/Galileo at L5
                                    other | GPS at L1 removing the C/A code band. 
    */
	double load_ITF_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double BF_phases_UP[8], double BF_phases_DW[8], int filter_num );

    /** \brief Load a cluster of 999 complex interferometric waveforms from a 1-second IEEC's SPIR binary file by selecting a set of input signals (16 available).

        \param namefile              Filename of SPIR 1-second binary file.
        \param peak_delay_estimate   Raw estimation of the range delay of the waveform's peak in meters.
        \param BF_phases             Beamformer phases (in degrees) to be applied at each of the 16 antenna elements of the SPIR setup to get maximum antenna gain towards a desired direction. Note that only those elements indicated on input variables `signal_elements_1` and `signal_elements_2` will be employed.
        \param signal_elements_1     Array indicator of which antenna elements (those indexes with a value equal to '1') will be combined to generate signal 1 (to be cross-correlated against signal 2).
        \param signal_elements_2     Array indicator of which antenna elements (those indexes with a value equal to '1') will be combined to generate signal 2 (to be cross-correlated against signal 1).
        \param filter_num            Code to select the frequency-domain filter to be applied during the processing. 
                                     Value | Meaning
                                    :-----:|:--------
                                       1   | Galileo at E1
                                       2   | GPS/Galileo at L5
                                     other | GPS at L1 removing the C/A code band. 
    */
	double load_ITF_waveforms_SPIR_selected_signals( const char* namefile, double peak_delay_estimate, double BF_phases[16], signed char signal_elements_1[16], signed char signal_elements_2[16], int filter_num );

    /** \brief Load a cluster of 999 complex GPS L1-CA clean-replica waveforms from a 1-second IEEC's SPIR binary file.

        \param namefile              Filename of SPIR 1-second binary file.
        \param peak_delay_estimate   Raw estimation of the range delay of the waveform's peak in meters.
        \param doppler_estimate      Estimation of the Doppler frequency of the direct signal in Hz (also required when computing waveforms from reflected signals).
        \param delta_doppler         Estimation of the Doppler frequency difference of the reflected signal with respect to the direct one in Hz (set to zero in case of computing waveforms from direct signals).
        \param BF_phases             Beamformer phases (in degrees) to be applied at each of the 16 antenna elements of the SPIR setup to get maximum antenna gain towards a desired direction.
        \param uplooking_channel     Select which antennas have to be cross-correlated
                                     Value | Meaning
                                     :----:|:-------
                                      True | Cross-correlation against direct signals using up-looking antenna elements. 
                                      False| Cross-correlation against direct or reflected signals using down-looking antenna elements.
        \param code_ref              GPS PRN number. 
        \return                      Range delay of the first waveform sample in meters.        
    */
	double load_CR_waveforms_SPIR( const char* namefile, double peak_delay_estimate, double doppler_estimate, double delta_doppler, double BF_phases[16], bool uplooking_channel, int code_ref );

    /** \brief Load a cluster of 999 complex GPS L1-CA clean-replica waveforms from a 1-second IEEC's SPIR binary file by selecting a set of input signals (16 available).

        \param namefile              Filename of SPIR 1-second binary file.
        \param peak_delay_estimate   Raw estimation of the range delay of the waveform's peak in meters.
        \param doppler_estimate      Estimation of the Doppler frequency of the direct signal in Hz (also required when computing waveforms from reflected signals).
        \param assist_doppler        Enable search for Doppler
                                     Value | Meaning
                                    :-----:|:--------
                                     True  | Enable search of Doppler frequency (having `doppler_estimate` as initial reference). 
                                     False | Disable search of Doppler frequency (rely on `doppler_estimate`).
        \param BF_phases             Beamformer phases (in degrees) to be applied at each of the 16 antenna elements of the SPIR setup to get maximum antenna gain towards a desired direction. Note that only those elements indicated on input variable `signal_elements` will be employed.
        \param signal_elements       Array indicator of which antenna elements (those indexes with a value equal to '1') will be combined to generate the signal.
        \param code_ref              GPS PRN number. 
        \return                      Range delay of the first waveform sample in meters.
    */
	double load_CR_waveforms_SPIR_selected_signals( const char* namefile, double peak_delay_estimate, double doppler_estimate, bool assist_doppler, double BF_phases[16], signed char signal_elements[16], int code_ref );

    /** \brief Search of GPS L1-CA clean-replica signals from a 1-second IEEC's SPIR binary file by selecting a set of input signals (16 available). Basically, it scans the Doppler domain cross-correlating against a clean-replica model of each PRN in order to get the maximum response. Then, it prints on the screen the corresponding Doppler frequency and peak-to-noise ratio values obtained.
       
        \param namefile          Filename of SPIR 1-second binary file.
        \param signal_elements   Array indicator of which antenna elements (those indexes with a value equal to '1') will be evaluated.
    */
    void searching_CR_waveforms_SPIR( const char* namefile, signed char signal_elements[16] );

    /** \brief Get one of the stored complex waveforms.

        \param Icomp_out     Output - In-phase component of the complex waveform in arbitrary units.
        \param len_I         Number of samples of the in-phase component of the stored waveform (it has to be equal to `wav_length`). Recommendation: use `my_wavcluster.get_wav_length()` as `len_I`.
        \param Qcomp_out     Output - Quadrature component of the complex waveform in arbitrary units. 
        \param len_Q         Number of samples of the quadrature component of the stored waveform (it has to be equal to `wav_length`). Recommendation: use `my_wavcluster.get_wav_length()` as `len_Q`.
        \param cluster_pos   Position of the complex waveform to be got.
    */
	void get_waveform( double* Icomp_out, int len_I, double* Qcomp_out, int len_Q, int cluster_pos );

    /** \brief Integrate the cluster of complex waveforms, both coherent and incoherent, to get a power waveform.

        \param coherent_int     Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged.
        \param wav_out          Integrated power waveform in arbitrary units. 
        \param wav_out_length   Number of samples of the stored waveforms (it has to be equal to `wav_length`). Recommendation: use `my_wavcluster.get_wav_length()` as `wav_len`.
    */
	void integrate_waveforms( int coherent_int, double* wav_out, int wav_out_length );

    /** \brief Integrate the cluster of complex waveforms, both coherent and incoherent, to get a power waveform after removing a longer coherent component to mitigate direct signal interference. 

        \param coherent_int       Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged.
        \param coherent_int_dir   Number of samples of coherent integration for mitigation of direct signal interference (from 1 to `cluster_length`).
        \param wav_out            Output - Integrated power waveform in arbitrary units.
        \param wav_out_length     Number of samples of the stored waveforms (it has to be equal to `wav_length`). Recommendation: use `my_wavcluster.get_wav_length()` as `wav_len`.
    */
	void integrate_waveforms_remdir( int coherent_int, int coherent_int_dir, double* wav_out, int wav_out_length );

    /** \brief Integrate the cluster of complex waveforms, both coherent and incoherent, to get a power waveform after applying a given complex retracking along the cluster. In this context, "to apply a retracking" means to rotate and move in delay each complex waveform given a range value with the purpose of better align them before integration.

        \param coherent_int               Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged.
        \param sampling_rate              Sampling rate of the waveforms from the cluster in samples/sec.
        \param retracking_meters          Retracking (in meters) to be applied to each waveform of the cluster before the integration.
        \param retracking_meters_length   Length of the `retracking_meters` array.
        \param wav_out                    Output - Integrated power waveform in arbitrary units.
        \param wav_out_length             Number of samples of the stored waveforms (it has to be equal tob `wav_length`. Recommendation: use `my_wavcluster.get_wav_length()` as `wav_len`.
    */
	void integrate_waveforms_retracking( int coherent_int, double sampling_rate, double* retracking_meters, int retracking_meters_length, double* wav_out, int wav_out_length );

    /** \brief  Print a time series of the phase from the stored waveform's cluster at a given lag position (only where there are valid values). 

        \param lag_pos   Lag position at the waveform's cluster selected to print out a time series of its phase.
    */
	void dump_phase( int lag_pos );

    /** \brief Print a time series of the phase from the stored waveform's cluster at a the peak position (only where there are valid values).
    */ 
	void dump_phase_peak( void );

    /** \brief Store the contents of the waveform's cluster at a given lag in a parallel phasor for analysis purposes. 

        \param lag_pos   Lag position at the waveform's cluster selected to store its contents in a parallel phasor.
    */
	void store_phasor_wavs( int lag_pos );

    /** \brief Get the stored phasor.

        \param phasorI_out          In-phase component of the stored phasor.
        \param phasor_length_I      Number of samples of the in-phase component of the stored phasor (it has to be equal to `cluster_length`). Recommendation: use `my_wavcluster.get_cluster_length()` as `len_Iphasor`.
        \param phasorQ_out          Quadrature component of the stored phasor.
        \param phasor_length_Q      Number of samples of the quadrature component of the stored phasor (it has to be equal to `cluster_length`). Recommendation: use `my_wavcluster.get_cluster_length()` as `len_Qphasor`.              
        \param valid_phasor_out     Array that indicates the validity of the stored phasor 
                                    Value | Meaning
                                    :----:|:-------
                                      1   |  valid
                                      0   | non-valid
        \param phasor_length_bool   Number of samples of the array that indicates the validity of the stored phasor (it has to be equal to `cluster_length`). Recommendation: use `my_wavcluster.get_cluster_length()` as `phasor_length_bool`.
    */
	void get_phasor( double* phasorI_out, int phasor_length_I, double* phasorQ_out, int phasor_length_Q, signed char* valid_phasor_out, int phasor_length_bool );
    
    /** \brief Compute the standard deviation of the phase from the stored phasor with the requirement of a given minimum number of valid samples.

        \param min_valid_samples   Minimum number of valid samples in stored phasor to compute the standard deviation of the phase.
        \return                    -1 if not enough valid samples.
    */
	double get_sigma_phase_phasor( int min_valid_samples );

    /** \brief Compute the standard deviation of the phase from the stored phasor within a given interval and with the requirement of a given minimum number of valid samples.

        \param init_sample         Initial sample of the selected interval in the stored phasor to compute the standard deviation of the phase.
        \param interv_samples      Number of samples of the selected interval in the stored phasor to compute the standard deviation of the phase.
        \param min_valid_samples   Minimum number of valid samples within the selected interval of the stored phasor to compute the standard deviation of the phase.
        \return                    -1 if not enough valid samples.
    */
	double get_sigma_phase_phasor_interv( int init_sample, int interv_samples, int min_valid_samples );

    /** \brief Counter-rotate the stored phasor with a given array of phases.

        \param phases_rad          Time series of input phases (in radians) employed to counter-rotate the stored phasor.    
        \param phases_length       Length of array `ohases_rad`
        \param valid_phases        Array that indicates the validity of `phases_rad`
                                     Value | Meaning
                                    :-----:|:-------
                                      1    | valid
                                      0    | non-valid
        \param phases_length_bool  Length of array `valid_phases`
    */
	void counterrot_phasor( double* phases_rad, int phases_length, signed char* valid_phases, int phases_length_bool );

    /** \brief Counter-rotate the stored complex waveforms with a given array of phases. 

        \param phases_rad           Time series of input phases (in radians) employed to counter-rotate the stored complex waveforms.
        \param phases_length        Length of array `phases_rad` 
        \param valid_phases         Array that indicates the validity of `phases_rad`
                                    Value | Meaning
                                    :----:|:-------
                                      1   | valid
                                      0   | non-valid
        \param phases_length_bool   Length of array `valid_phases`
    */
	void counterrot_waveforms( double* phases_rad, int phases_length, signed char* valid_phases, int phases_length_bool );

    /** \brief Apply an algorithm to correct the navigation bit (for GPS L1 C/A code) in all the waveforms from the stored cluster. Such algorithm analyzes the phase variations at a given lag and decides where there is a bit transition. Since coherence and proper SNR are required, this function is recommended for waveform clusters containing direct GPS L1 signals; however, the information obtained could be later applied to the clusters containing their corresponding reflections (if available).

        \param lag_pos                Lag position where the algorithm is applied (peak is recommended).                                                                                           
        \param store_navbit_phasorI   Tells if the navigation bit shall be stored.
                                      Value | Meaning
                                      :----:|:-------
                                        1   | store the time series of the navigation bit estimated at the in-phase component of the stored phasor
                                        0   | do not store the navigation bit
    */
	void correct_navigation_bit( int lag_pos, int store_navbit_phasorI );

    /** \brief Apply an algorithm to estimate the coherence time from the stored cluster. Such algorithm checks the accumulated phase variation at a given lag and sets the coherence time at the moment when such variation reaches 1~radian.

        \param lag_peak            Lag position where the algorithm is applied (peak is recommended).
        \param store_acf_phasorQ   Store or not the time series of the accumulated phase difference at the quadrature component of the stored phasor.
                                    Value | Meaning
                                    :----:|:-------
                                      1   | Store the time series of the accumulated phase difference at the quadrature component of the stored phasor.
                                      0   | Do not store the navigation bit.
    */
	double compute_coherence_time( int lag_peak, int store_acf_phasorQ );

    /** \brief Compute a DDM's Doppler slice from the waveform cluster. The procedure is the same as doing an integration, but after modulating the waveforms using a phasor at the given frequency. 

        \param coherent_int  Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged. 
        \param doppler_freq   Doppler frequency of the output DDM slice in inverse cluster sample units.
        \param freq_ddm       DDM's Doppler slice in arbitrary units. 
        \param delay_samples  Number of delay samples of the DDM (it has to be equal to `wav_length`). Recommendation: use `my_wavcluster.get_wav_length()` as `ddm_lag_len`.
    */
	void compute_singlefreq_DDM( int coherent_int, double doppler_freq, double* freq_ddm, int delay_samples );

    /** \brief Compute a DDM's Doppler slice from the waveform cluster after removing a longer coherent component to mitigate direct signal interference. The procedure is the same as doing an integration, but after modulating the waveforms using a phasor at the given frequency.

        \param coherent_int       Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged. 
        \param coherent_int_dir   Number of samples of coherent integration for mitigation of direct signal interference (from 1 to `cluster_length).
        \param doppler_freq       Doppler frequency of the output DDM slice in inverse cluster sample units.
        \param freq_ddm           DDM's Doppler slice in arbitrary units.
        \param delay_samples      Number of delay samples of the DDM (it has to be equal to `wav_length`). Recommendation: use `my_wavcluster.get_wav_length() as `ddm_lag_len`.
    */
	void compute_singlefreq_DDM_remdir( int coherent_int, int coherent_int_dir, double doppler_freq, double* freq_ddm, int delay_samples );

    /** \brief  Compute a DDM's delay slice from the waveform cluster. The procedure is the same as doing an integration, but after modulating the waveforms using a phasor at a given set of frequencies.

        \param coherent_int  Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged. 
        \param lag_pos       Delay lag of the output DDM slice.
        \param delta_freq    Doppler frequency resolution of the DDM in inverse cluster sample units.
        \param lag_ddm       
        \param freq_samples  Number of Doppler bins in the DDM, starting from `(-ddm_freq_len/2)*delta_freq`.
    */
	void compute_singlelag_DDM( int coherent_int, int lag_pos, double delta_freq, double* lag_ddm, int freq_samples );

    /** \brief Compute a DDM's delay slice from the waveform cluster after removing a longer coherent component to mitigate direct signal interference. The procedure is the same as doing an integration, but after modulating the waveforms using a phasor at a given set of frequencies.

        \param coherent_int       Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged.
        \param coherent_int_dir   Number of samples of coherent integration for mitigation of direct signal interference (from 1 to `cluster_length`).
        \param lag_pos            Delay lag of the output DDM slice.
        \param delta_freq         Doppler frequency resolution of the DDM in inverse cluster sample units.
        \param lag_ddm            DDM's delay slice in arbitrary units.
        \param freq_samples       Number of Doppler bins in the DDM, starting from `(-ddm_freq_len/2)*delta_freq`.
    */
	void compute_singlelag_DDM_remdir( int coherent_int, int coherent_int_dir, int lag_pos, double delta_freq, double* lag_ddm, int freq_samples );

    /** \brief Compute the 3~dB bandwidth from the DDM at the Doppler domain for a given lag.

        \param coherent_int   Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged.
        \param lag_pos        Delay lag of the DDM slice.
        \param freq_samples   Number of Doppler bins in the DDM, starting from `(-ddm_freq_len/2)*delta_freq`.
        \param delta_freq     Doppler frequency resolution of the DDM in inverse cluster sample units.
        \param pos_pow_Max    Output - Peak location of the DDM's delay slice in inverse sample units (first element). Peak value of the DDM's delay slice in arbitrary units (second element).
        \return               TODO
    */
	double compute_DopplerMap_BW( int coherent_int, int lag_pos, int freq_samples, double delta_freq, double pos_pow_Max[2] );

    /** \brief Compute the 3 dB bandwidth from the DDM at the Doppler domain for a given lag and after removing a longer coherent component to mitigate direct signal interference. 

        \param coherent_int       Number of samples of coherent integration (from 1 to `cluster_length`). The complex waveforms are coherently integrated in blocks of `coherent_int` samples and the results are incoherently averaged.
        \param coherent_int_dir   Number of samples of coherent integration for mitigation of direct signal interference (from 1 to `cluster_length`).
        \param lag_pos            Delay lag of the DDM slice.
        \param freq_samples       Number of Doppler bins in the DDM, starting from `-(ddm_freq_len/2)*delta_freq`.
        \param delta_freq         Doppler frequency resolution of the DDM in inverse cluster sample units.
        \param pos_pow_Max        Output - Peak location of the DDM's delay slice in inverse cluster sample units (first element). Peak value of the DDM's delay slice in arbitrary units (second element).
        \return                   TODO
    */
	double compute_DopplerMap_BW_remdir( int coherent_int, int coherent_int_dir, int lag_pos, int freq_samples, double delta_freq, double pos_pow_Max[2] );

    /** \brief TODO 

        \param coherent_int        TODO
        \param ddm_delta_freq      TODO
        \param ddm                 TODO
        \param ddm_delay_samples   TODO
        \param ddm_freq_samples    TODO
    */
	void compute_whole_DDM( int coherent_int, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples );

    /** \brief TODO

        \param coherent_int        TODO
        \param coherent_int_dir    TODO
        \param ddm_delta_freq      TODO
        \param ddm                 TODO
        \param ddm_delay_samples   TODO
        \param ddm_freq_samples    TODO
    */
	void compute_whole_DDM_remdir( int coherent_int, int coherent_int_dir, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples );

    /** \brief TODO

        \param coherent_int               TODO
        \param ddm_delta_freq             TODO
        \param ddm                        TODO
        \param ddm_delay_samples          TODO
        \param ddm_freq_samples           TODO
        \param sampling_rate              TODO
        \param retracking_meters          TODO
        \param retracking_meters_length   TODO
    */
	void compute_whole_DDM_retracking( int coherent_int, double ddm_delta_freq, double** ddm, int ddm_delay_samples, int ddm_freq_samples, double sampling_rate, double* retracking_meters, int retracking_meters_length );

    /** \brief TODO

        \param powerLagHolo
        \param lagHolo_lags
        \param lagHolo_freqs
    */
	void compute_whole_LagHologram( double** powerLagHolo, int lagHolo_lags, int lagHolo_freqs );

    /** \brief Compute a lag-hologram delay slice from the waveform cluster. The procedure consist in computing the FFT algorithm from the time series at the given lag position from the cluster of waveforms.

        \param lag                      Delay lag of the output lag-hologram slice.
        \param powerLagHolo_singleLag   Lag-hologram delay slice in arbitrary units. The corresponding frequencies range from `-(int(fft_len/2) - 1)` to `int(fft_len/2)` in inverse cluster sample units.  
        \param fft_samples              Number of FFT samples (from 2 to `cluster_length`, being a power of 2). 
    */
	void compute_LagHologram( int lag, double* powerLagHolo_singleLag, int fft_samples );

    /** \brief Get the length of the waveforms stored. 

        \return  Length of the waveforms stored.
    */
	int get_wav_length( void );

    /** \brief  Get the number of waveforms stored in the cluster.

        \return   Number of waveforms stored in the cluster.
    */
	int get_cluster_length( void );
};

void integrate_wavs( int coh_int, int coh_int_dir, int cl_length, int w_length, double freq, bool* valid, double** Icomp, double** Qcomp, double* wav, bool retracking, double* retrack_samples );
double integrate_wav_lag( int coh_int, int coh_int_dir, int cl_length, int lag, double freq, bool* valid, double** Icomp, double** Qcomp );
double compute_sigma_phase_phasor( int init, int interval, int min_samples, double* phasI, double* phasQ, bool* valid );
bool check_SPIR_blocks( int* spirdata );
void compute_RF_filter( float* filter );
void compute_RF_filter_GalileoE1( float* filter );
void compute_RF_filter_GPS_L5( float* filter );
void compute_ITFwav_SPIR( int* spirdata, double interm_freq, float* filter, double* wav_real, double* wav_imag, int init_sample, int wav_size, double* phases_BF, bool* elementsUP, bool* elementsDW, int msec );
bool compute_CRwav_SPIR( int* spirdata, double interm_freq, float* model_CR, double doppler_CR, int samples_CR, float* filter, double* wav_real, double* wav_imag, int wav_size, double* phases_BF, bool* elements_CORR, int msec, int &init_sample, bool compute_init_sample );
double getSPIRdatabit( int word, int pos );
void get_doppler_initsamp_CR_SPIR( int* spirdata, double interm_freq, float* model_CR, double doppler_CR, int samples_CR, float* filter, int wav_size, double* phases_BF, double &doppler_peak, int &init_sample, bool selected_signals, bool elements_CORR[16], bool searching );
void read_PRN_code( int code_ref, float* code_CR, int code_length );

