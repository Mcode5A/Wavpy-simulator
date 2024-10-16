/**
  \file ZavorotnyVoronovichModel.h
  \brief Public header of the ZaVoModel_GNSSR class.
*/

/**
  \brief The ZaVoModel_GNSSR class allows to simulate GNSS-R waveforms and DDMs.

  This class models GNSS-R waveforms and DDM's based on \cite zavorotny2000sgs and their corresponding covariance matrix based on \cite Li2017Revisiting, which allows the simulation of realistic noise (thermal and speckle) realizations with a proper statistical characterization in both range and Doppler domains.
  

*/

class ZaVoModel_GNSSR
{
  private: // variables
	bool dump_isolines_data;             /**< \brief Indicates if isolines (delay and Doppler) shall be written to a binary file.
                                              Value | Description
                                              :----:|:-----------
                                              True  | Store a binary file containing information of range, Doppler, sigma_0 and gain of the antenna over the surface. 
                                              False (default) | Do not store surface data. 
                                         */
    bool stare_processing_flag;          /**< \brief Indicate if the simulation for the stare processing shall be done or not.
                                              Value | Description
                                              :----:|:-----------
                                              True  | Make simulation for a stare processing case, which will produce a DDM with a PDF of slopes set to 1.
                                              False (default) | Make the simulation using the standard computation of the PDF of slopes.
                                         */
    bool use_transmitter_eirp;           /**< \brief Indicate which value for the EIRP of the transmitter has to be used.
                                              Value | Description
                                              :----:|:-----------
                                              True  | Use txr_eirp for the computation of transmitted power. 
                                              False (default) |  Standard computation of transmitted power by applying the code weights in gnss\_signal to the nominal values of the different code-components. 
                                         */
    double txr_eirp;                     /**< \brief Equivalent Isotropic Radiated Power (EIRP) in dB units employed for the computation of transmitted power when use_transmitter_eirp is True. Default value 0.0 */
	std::string isolines_data_namefile;  /**< \brief Name of the binary output file containing the information of range, Doppler, sigma_0 and gain of the antenna over the surface when making the simulation.
                                              Default value: void 
                                         */
    std::string stareProc_data_namefile; /**< \brief Name of the text output file containing the information related to stare processing in an eight-columns format:
                                              range | Doppler | mean X-component of scattering vector | mean Y-component of scattering vector | mean Z-component of scattering vector | mean longitude | mean latitude | ambiguity status
                                              :----:|:-------:|:-------------------------------------:|:-------------------------------------:|:-------------------------------------:|:--------------:|:-------------:|:----------------
                                               ...  | ...     | ...                                   | ...                                   | ...                                   | ...            | ...           | ...     
                                             ("1" if there is no ambiguity, "0" else). Default value: void
                                         */ 
	int size_ddm_stored[2];              /**< \brief (integer, array of 2 elements): Number of samples of the stored DDM without taking into account the central Doppler. Default value: [0, 0] */
	double **ddm;                        /**< \brief 2D array of size_ddm_stored[0] x size_ddm_stored[1] elements. Simulated power DDM. Default value: void */
	int len_cov_stored;                  /**< \brief Number of samples (single dimension) of the stored covariance. Default value: 0 */
	double **cov;                        /**< \brief Simulated waveform or DDM covariance. Default value: void*/
	double **chol;                       /**< \brief TBC */
  public: 
    // VARIABLES
	bool interferometric_flag;           /**< \brief Indicates if computation of noise level is based on interferometric/clean-replica approach
                                              Value | Description
                                              :----:|:-------------
                                              True  | Computation of noise level based on the interferometric approach. 
                                              False (default) | Computation of noise level based on the clean-replica approach.
                                         */
	bool curvature_approx_flag;          /**< \brief Apply or not Earth curvature approximation when integrating.
                                              Value | Description
                                              :----:|:-----------
                                              True  | Apply Earth curvature approximation when integrating the reflected power over the surface (recommended when the receiver is at high altitude).
                                              False (defaule) | Apply planar reflected surface. Default value: False
                                         */
	bool covariance_wav_flag;            /**< \brief Indicates if the mean waveform has been computed using the covariance approach or not.
                                              Value | Description
                                              :----:|:-----------
                                               True | Computation of mean waveform based on the covariance approach (the corresponding covariance matrix is also stored).
                                               False (default)| Computation of mean waveform based on the regular approach (if covariance\_ddm\_flag is also at False).
                                         */
	bool covariance_ddm_flag;            /**< \brief Indicates if the mean waveform and DDM has been computed based on the covariance approach or not.
                                              Value | Description
                                              :----:|:-----------
                                               True | Computation of mean waveform and DDM based on the covariance approach (the corresponding covariance matrix is also stored).
                                               False (default) | Computation of mean DDM based on the regular approach.
                                         */
	bool recompute_Lambda_flag;          /**< \brief Indicates if the the autocorrelation is computed before the simulation or not.
                                              Value | Description
                                              :----:|:-----------
                                              True (default) | Computation of autocorrelation function from gnss\_signal before simulation.
                                              False | Use of already stored autocorrelation function in gnss\_signal during simulation.
                                         */
    bool sigma0_to_one_flag;             /**< \brief Indicates if sigma_0 has been set to 1 for computation fo the DDM or not.
                                              Value | Description
                                              :----:|:-----------
                                              True  | Computation of mean waveform or DDM by setting sigma\_0 to 1. 
                                              False (default) | Computation of mean waveform or DDM based on the regular approach.
                                         */
	char polarization;                   /**< \brief Polarization of the reflected signal. 
                                              Valid values | Description 
                                              :-----------:|:------------
                                                'R'        |  RHCP 
                                                'L' (default) | LHCP. 
                                         */ 
	int num_angles;                      /**< \brief Number of integration points over the surface ellipse for each range sample.
                                              It determines the angular resolution of the surface integral during the simulation. Default value is 120.
                                         */
	int wav_length;                      /**< \brief Length of modelled waveform. Default value is 256. */
	int ddm_half_dopplers;               /**< \brief Number of additional Doppler slices. The total number of Doppler lines in the DDM will be: 2x ddm\_half\_dopplers + 1. Default value: 0*/
	double sampling_rate;                /**< \brief Sampling rate of the modelled waveform in samples/sec. Default value is 80e6 */ 
	double delta_doppler;                /**< \brief Doppler increment between consecutive slices of the DDM in Hz. Default value is 0.0 */
	double delta_freq;                   /**< \brief Doppler offset of the DDM in Hz. Default value is 0.0 */
	double coherent_integration;         /**< \brief Coherent integration time of the modelled waveform in secs. Default value is 1e-3 */
	double weight_cohPower;              /**< \brief Weight applied to the coherent component of the reflected power. Default value: 0.0 */
    Specular_geometry geometry;          /**< \brief Defines the geometry used during the simulation. */
	Reflecting_surface surface;          /**< \brief Defines the properties of the reflecting surface during simulation. */
	RF_FrontEnd receiver_Up;             /**< \brief Defines the properties of the up-looking receiver front end, that collects the direct GNSS signals, used during the simulation. */
	RF_FrontEnd receiver_Down;           /**< \brief Defines the properties of the down-looking reciever front end, that collects the reflected GNSS signals, used during the simulation. */
	GNSS_composite gnss_signal;          /**< \brief Defines the properties of the GNSS signal used during the simulation. */
	Waveform_power waveform_POW;         /**< \brief Stores the simulated waveform. */  
    // METHODS
    // Constructors 
    /** \brief Constructor for defining a new object.
      @return      the simulated GNSSR object
    */
    ZaVoModel_GNSSR( void );

    // Other methods
    /** \brief Enable the storage of a binary file containing information of range, Doppler, sigma_0 and gain of the antenna over the surface during the simulation. 

      @param namefile   The name of the file to which the data shall be written to.

      The binary file contains blocks of 8 float numbers:
        Field   | Description
      :--------:|:-----------
         x      | x coordinate of the point location in the local frame expressed in m
         y      | y coordinate of the point location in the local frame expressed in m
         lon    | longitude of the point location expressed in degrees
         lat    | latitude of the point location expressed in degrees
         tau    | the range corresponding to the point expressed in m
       sigma_0  | the sigma_0 at the point expressed in dB
       gain     | the receiver antenna gain at the point expressed in dB
       doppler  | the Doppler frequency at the point expressed in Hz

      The next lines show a python example on how to read the contents of file and load them into a set of numpy arrays:

      \verbatim
        import numpy as np
        
        data_surf = np.fromfile(binfile, dtype=np.float32)
        x = np.zeros(len(data_surf)/8, dtype=np.float32)
        y = np.zeros(len(data_surf)/8, dtype=np.float32)
        lon = np.zeros(len(data_surf)/8, dtype=np.float32)
        lat = np.zeros(len(data_surf)/8, dtype=np.float32)
        tau = np.zeros(len(data_surf)/8, dtype=np.float32)
        sigma_0 = np.zeros(len(data_surf)/8, dtype=np.float32)
        gain = np.zeros(len(data_surf)/8, dtype=np.float32)
        doppler = np.zeros(len(data_surf)/8, dtype=np.float32)
        
        for index in range(len(data_surf)/8):
           x[index] = data_surf[index*8]
           y[index] = data_surf[index*8 + 1]
           lon[index] = data_surf[index*8 + 2]
           lat[index] = data_surf[index*8 + 3]
           tau[index] = data_surf[index*8 + 4]
           sigma_0[index] = data_surf[index*8 + 5]
           gain[index] = data_surf[index*8 + 6]
           doppler[index] = data_surf[index*8 + 7]
      \endverbatim
    */
	void enable_isolines_data_dump( const char* namefile );

    /** \brief Disable the storage of a binary file containing information of range, Doppler, sigma_0 and gain of the antenna over the surface during the simulation. */
	void disable_isolines_data_dump( void );
  
    /** \brief Set simulation into stare processing mode. Under this condition, the waveform/DDM is simulated by setting the PDF of the slopes to 1.
               In addition, a text file is stored with further information of the DDM process, such as mean scattering vector and ambiguity status of each delay-Doppler cell.

      @param namefile     Name of the text output file containing the information related to stare processing.

      The text file has an eight-columns format:
                Column  | Parameter 
                :------:|:---------
                  1     | range
                  2     | Doppler
                  3     | mean X-component of scattering vector
                  4     | mean Y-component of scattering vector
                  5     | mean Z-component of scattering vector
                  6     | mean longitude
                  7     | mean latitude 
                  8     | ambiguity status ("1" if there is no ambiguity, "0" else)
    */
    void set_stare_processing_mode( const char* namefile );

    /** \brief Unset the stare processing mode from the simulation. */
	void unset_stare_processing_mode( void );

    /** \brief Use of a given value for the computation of transmitted power. 

        @param txr_eirp_in    Equivalent Isotropic Radiated Power (EIRP) in dB units employed for the computation of transmitted power.
    */
    void set_transmitter_EIRP( double txr_eirp_in );

    /** \brief Unset the use of a previously loaded value for the computation of transmitted power. Then, such computation is made by 
        applying the code weights in gnss_signal to the nominal values of the different code-components.
    */
    void unset_transmitter_EIRP( void );

    /** \brief Simulate the waveform and the (DDM if ddm_half_dopplers and delta_doppler are greater than zero) for the given characterization. */
	void compute_waveform( void );

    /** \brief Get a Doppler slice from the stored DDM (if it has been computed previously).

      @param doppler_index   Index of the Doppler slice (0 for the central frequency)
      @param dm_out          Output Doppler slice from the stored DDM
      @param dm_out_length  The number of elements of the doppler slice array 
    */
	void get_DDM_doppler_slice( int doppler_index, double* dm_out, int dm_out_length );
	
    /** \brief Get a slice of the stored covariance matrix (if it has been computed).
 
      @param cov_index       Index of the covariance slice
      @param cov_out         Output covariance slice from the stored covariance.
      @param cov_out_length  Number of range samples of the stored covariance (it has to be equal to len_cov_stored). For a waveform covariance, this number shuld be equal to the number of range samples of the corresponding waveform. In the case of a DDM, the previous number should be additionally multiplied by the number of Doppler samples.
    */
    void get_cov_slice( int cov_index, double* cov_out, int cov_out_length );

    /** \brief Get a noisy waveform computed from the stored covariance and the mean waveform.

      @param wav_out        Noisy output power waveform in arbitrary units.
      @param wav_out_length Number of samples of the stored mean waveform (it has to be equal to wav_length). Recommendation: use  my_model.waveform_POW.get_wav_length() as wav_len.
      @param seed_in        Seed used for the internal random Gaussian noise generator. 
                            TODO: Explain what random generator is used internally.
    */ 
	void get_noisy_waveform( double* wav_out, int wav_out_length, unsigned long int seed_in );

    /** \brief Get a noisy DDM computed from the stored covariance and the mean DDM.
      
      @param ddm_row_out          Noisy output DDM in arbitrary units.
                                  Recommendation: use ddm_out.reshape((my_model.ddm_half_dopplers*2 + 1), my_model.wav_length)
                                  to convert it into a 2-dimensional array.

      @param ddm_row_out_length   Number of samples of the stored mean DDM (it has to be equal to size_ddm_stored[0] multiplied by size_ddm_stored[1] + 1). 
                                  Recommendation: use my_model.waveform_POW.get_wav_length()*(my_model.ddm_half_dopplers*2 + 1) as ddm_len.
      @param seed_in              Seed used for the internal random Gaussian generator.
                                  TODO: Explain what random generator is used internally.
    */
	void get_noisy_DDM( double* ddm_row_out, int ddm_row_out_length, unsigned long int seed_in );
};

/** \brief Computes the transmitted GNSS power

    \param  signal    The GNSS signal under consideration
    \param  elev      The elevation of the signal
    \param  atm_loss  The atmospheric losses expressed in [TODO]
    \param  posTnorm  TODO
    \return           The transmitted power
*/
static double Compute_power_trans( GNSS_composite signal, double elev, double atm_loss, double posTnorm );  //TODO this function should be static. TODO Confirm if it should also be private.

/** \brief Computes the Doppler frequency of the reflections path. ECEF frame is used in the calculation TODO. Should this function be moved to some geometry class TODO?

    \param n_scat      Unit (TODO) vector indicating the propagation direction of the scattering direction (TODO: confirm direction, from Rx to S, or from S to Rx?)
    \param velR        Velocity vector of the transmitter (TODO: which units?)
    \param n_inc       Unit (TODO) vector indicating the propagation direction of the incidence direction (TODO: confirm direction, from Tx to S, or from S to Tx?)
    \param velT        Velocity vector of the transmitter (TODO: which units?)
    \param frequency   Carrier frequency of the signal expressed in Hz
    \return            The computed Doppler frequency expressed in Hz
*/
double Doppler_func( double n_scat[3], double velR[3], double n_inc[3], double velT[3], double frequency );


/** \brief Computes the sigma_0 for a given surface and geometry. Should this function be moved to the surface class TODO?

    \param surf            The reflecting surface
    \param n_scat          Unit (TODO) vector indicating the propagation direction of the scattering (TODO: confirm direction, from Rx to S or from S to Rx?)
    \param n_inc           Unit (TODO) vector indication the propagation direction of the incidence (TODO: confirm direction, from Tx to S or from S to Tx?)
    \param pol             Polarization of the wave
    \param azimuthT        TODO
    \param set_PDF_to_one  Indicate if the slopes pdf has to set to 1 (TODO: Does it make sense to have it to 1? Can we have normalization problems in the sigma_0?)
    \return                The computed sigma_0 in [m^2/m^2]
*/
double Sigma0_func( Reflecting_surface surf, double n_scat[3], double n_inc[3], char pol, double azimuthT, bool set_PDF_to_one );

/** \brief                 Computes the PD (TODO) of the slopes distribution for a given surface (TODO: confirm)

    \param surf            The reflecting surface
    \param q               TODO
    \param azimuthT        TODO
    \return                TODO
*/
double PDFunction( Reflecting_surface surf, double q[3], double azimuthT );

/** \brief TODO

    \param ant_k             TODO
    \param ant_E             TODO
    \param ant_H             TODO
    \param n_target_to_ant   TODO
    \param out_angles        TODO  is the return data
*/
void AntGain_EH_angles( double ant_k[3], double ant_E[3], double ant_H[3], double n_target_to_ant[3], double out_angles[2] );


/** \brief TODO

    \param ant_k             TODO
    \param ant_E             TODO
    \param ant_H             TODO
    \param n_target_to_ant   TODO
    \param out_angles        TODO  is the return data
*/
void AntGain_PhiTheta_angles( double ant_k[3], double ant_E[3], double ant_H[3], double n_target_to_ant[3], double out_angles[2] );

/** \brief  Computes the DMM using the covariance method TODO

    \param power_direct_surf      TODO
    \param power_direct_rec       TODO
    \param sampling_rate          TODO
    \param coh_time               TODO
    \param BW                     TODO
    \param pow_nd                 TODO
    \param pow_nr                 TODO
    \param delta_doppler          TODO
    \param power_surf             TODO
    \param doppler_len            TODO
    \param delay_len              TODO
    \param lambda_func            TODO
    \param lambda_len             TODO
    \param wav_out                TODO 
    \param ddm_out                TODO
    \param cov_out                TODO 
    \param chol_out               TODO
    \param wav_len                TODO
    \param interferometric        TODO
    \param ddm_freq_factor        TODO
    \param compute_ddm            TODO
*/
void Compute_Covariance_DDM( double power_direct_surf, double power_direct_rec, double sampling_rate, double coh_time, double BW, double pow_nd, double pow_nr, double delta_doppler, double** power_surf, int doppler_len, int delay_len, double* lambda_func, int lambda_len, double* wav_out, double** ddm_out, double** cov_out, double** chol_out, int wav_len, int interferometric, int ddm_freq_factor, bool compute_ddm );

/** \brief TODO 

    TODO say what random generator is used.
    \param out_string            Why are almost all input parameters called xxx_string? Output value TODO
    \param mean_string           TODO
    \param len_string            TODO
    \param chol_mat              TODO
    \param seed_rn               Some random seed. TODO
*/
void Get_noise_string( double* out_string, double* mean_string, int len_string, double** chol_mat, int seed_rn );

/** \brief TODO Some kind of auxiliary function. By its behaviour, it should probably defined as a static and private method. 

    \param current_cell          TODO
    \param prev_cell             TODO
    \param status                
*/
void Check_status_Doppler_cell( bool current_cell, bool prev_cell, int &status );
