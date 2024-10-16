/**
  \file rf_front_end.h
  \brief Header of the RF_FrontEnd class.
*/

/**
  \brief The RF_FrontEnd class characterizes the main aspects of a receiver front end.

  The following definitions have been used in the implementation of this class.

  __Receiver frame:__ The receiver body frame is a Cartesian coordinates system of the structure containing the receiver, typically a satellite or an aircraft, where the X-axis points towards the front, the Y-axis points towards the right-side (XY define the horizontal plane) and the Z-axis towards Nadir. In absence of inertial rotation of the body frame, we will assume that the X-axis points towards the Earth's North and the Z-axis points towards the Earth's center.

  __Antenna frame:__ the antenna frame is a Cartesian coordinates system with its origin at the antenna's physical center, X- and Y-axis defining the physical plane of the antenna and the Z-axis pointing perpendicularly towards the propagation direction by complying the right-hand rule. Typically, this system is also represented with spherical coordinates, using radial distance, polar angle (\f$\theta\f$) and azimuth angle (\f$\phi\f$) as commonly used in physics (ISO convention). The antenna gain pattern is given using this type of representation. We define the E-plane and the H-plane as XZ-plane and YZ-plane respectively.

  __Relationship between variables:__ Whenever a user changes the value of any variable by means of the available functions, the rest of them are updated to keep consistency with the following equations:


    antenna_Aeff = (C_LIGHT/frequency)^2 * 10.0**(antenna_Gain_dB/10.0)/(4.0 * pi)

    noise_T = antenna_T + 290.0 * (10.0**(noise_F_dB}/10.0) - 1.0)

    noise_pow_dBW = 10.0 * log10(K_BOLTZ * noise_T * filter_BB_BW)

where C\_LIGHT is the speed of light in vacuum and K_BOLTZ is the Boltzmann's constant.
*/
class RF_FrontEnd
{
	double antenna_pattern_dB[181][360];    /**< \brief Antenna pattern (for a single element) as a function of \f$\theta\f$ and \f$\phi\f$ in the antenna frame in dB. Default value: 0.0 for all elements */
	double antenna_vector_BF_E[3];          /**< \brief Antenna frame's X-axis in the receiver body frame (defines the orientation of the antenna in the receiver frame). Default value: [1.0, 0.0, 0.0] */
	double antenna_vector_BF_H[3];          /**< \brief Antenna frame's Y-axis in the receiver body frame (defines the orientation of the antenna in the receiver frame). Default value: [0.0, 1.0, 0.0] */  
	double antenna_vector_BF_k[3];          /**< \brief Antenna frame's Z-axis in the receiver body frame (defines the orientation of the antenna in the receiver frame). Default value: [0.0, 0.0, 1.0]   */  
	bool isotropic;                         /**< \brief Indicates if an isotropic antenna is used
                                                        Value | Meaning
                                                        :----:|:-------
                                                         True | An isotropic antenna is employed (antenna pattern neglected)
                                               False (Default)| A non-isotropic antenna is employed.
                                             */ 
	double frequency;                       /**< \brief Frequency of the receiver in Hz. Default value: 1575420000.0 (GPS L1)  */  
	double antenna_Gain_dB;                 /**< \brief Antenna gain in dB (this value is added to the antenna pattern). Default value: 3.0  */  
	double antenna_Aeff;                    /**< \brief Effective area of the antenna in \f$ m^2\f$. Default value: 0.00575  */  
	double antenna_T;                       /**< \brief Antenna temperature in K. Default value: 200.0  */  
	double noise_T;                         /**< \brief Noise temperature in K. Default value: 488.626  */  
	double noise_pow_dBW;                   /**< \brief Noise power in dBW. Default value: -134.719  */  
	double noise_F_dB;                      /**< \brief Noise figure in dB. Default value: 3.0  */   
	double filter_BB_BW;                    /**< \brief Base-band bandwidth of the receiver in Hz. Default value: 5000000.0  */ 
	//2D PLANAR ARRAY
	int array_num_elements;                 /**< \brief If \f$>1\f$, the antenna is a 2D planar array of such number of elements. Default value: 1 (single antenna)  */ 
	double **element_pos_AF;                /**< \brief Position of the array elements as 2D coordinates in plane XY of the antenna frame in meters. Legth of `array_num_elements`x 2. Default value: void  */ 
	double *phase_delay;                    /**< \brief Phase applied to each element to obtain a desired array factor in radians. Length `array_num_elements`. Default value: void  */ 
	double array_factor_dB[181][360];       /**< \brief Array factor as a function of \f$\theta\f$ and \f$\phi\f$ in the antenna frame in dB. Default value: void  */ 
	bool array_factor_ready;                /**< \brief Indicates if the array factor has been comupted.
                                                        Value | Meaning
                                                        :----:|:-------
                                                        True  | Array factor is computed
                                              False (Default) | Array factor is not computed.
                                             */ 
  public:
    /** \brief Constructor
    */
	RF_FrontEnd( void );

    /** \brief Print relevant information about the object's content (public and private variables).
    */
	void dump_parameters( void );

    /** \brief Set the orientation of the antenna in the receiver body frame by editing private variables `antenna_vector_BF_E` and `antenna_vector_BF_H` (`antenna_vector_BF_k` is constructed using the right hand rule).

        \param antenna_vector_BF_E_in    Antenna frame's X-axis in the receiver body frame.
        \param antenna_vector_BF_H_in    Antenna frame's Y-axis in the receiver body frame.
    */
	void set_antenna_orientation_BF_EH( double antenna_vector_BF_E_in[3], double antenna_vector_BF_H_in[3] );

   /** \brief Set the pointing direction of the antenna in the receiver body frame by editing the private variable `antenna_vector_BF_k`. In this case, `antenna_vector_BF_E` and `antenna_vector_BF_H` are arbitrary constructed using the right hand rule (this method is valid when there is symmetry between both axis).

        \param antenna_vector_BF_k_in  Antenna frame's Z-axis in the receiver body frame.
    */
	void set_antenna_orientation_BF_k( double antenna_vector_BF_k_in[3] );

    /** \brief Get the antenna orientation.

        \param antenna_vector_BF_E_out   Output - Antenna frame's X-axis in the receiver body frame.
        \param antenna_vector_BF_H_out   Output - Antenna frame's Y-axis in the receiver body frame.
        \param antenna_vector_BF_k_out   Output - Antenna frame's Z-axis in the receiver body frame.
    */
	void get_antenna_orientation_BF( double antenna_vector_BF_E_out[3], double antenna_vector_BF_H_out[3], double antenna_vector_BF_k_out[3] );

    /** \brief Set the antenna pattern.

        \param antenna_whole_pattern_dB_in    Antenna pattern as a function of \f$\theta\f$ and \f$\phi\f$ in the antenna frame in dB. // TODO define better theta and phi
    */
	void set_antenna_whole_pattern( double antenna_whole_pattern_dB_in[181][360] );

    /** \brief Set an individual value of the antenna pattern.

        \param phi_AntFrame        Sample at coordinate \f$\phi\f$ of the antenna pattern (coincides with the integer value of angle \f$\phi\f$ in the antenna frame)
        \param theta_AntFrame      Sample at coordinate \f$\theta\f$ of the antenna pattern (coincides with the integer value of angle \f$\theta\f$ in the antenna frame).
        \param pattern_dB_value    Value of the antenna pattern in dB.
    */
	void set_val_antenna_pattern( int phi_AntFrame, int theta_AntFrame, double pattern_dB_value );

    /** \brief Set antenna pattern from a single full (for all \f$\theta\f$ angles) full (for \f$\phi\f$=0 and \f$\phi\f$=180) cut at the E-plane. It is assumed that both E- and H-planes are equal and interpolation is applied for the remaining points.

        \param antenna_full_pattern_dB_in    Antenna pattern in the E-plane as a function of \f$\theta\f$ angle in dB.
    */
	void set_antenna_pattern_FF( double antenna_full_pattern_dB_in[360] );

    /** \brief Set antenna pattern from a single full (for all \f$\theta\f$ angles) half (only for \f$\phi\f$=0) cut at the E-plane. It is assumed that both E- and H-planes are equal and interpolation is applied for the remaining points.

        \param antenna_half_pattern_dB_in    Antenna pattern in the E-plane as a function of \f$\theta\f$ angle in dB.
    */
	void set_antenna_pattern_FH( double antenna_half_pattern_dB_in[181] );

    /** \brief Set the antenna pattern by providing a set of points in the E-plane and a minimum level. The whole shape of the E-plane pattern is then built by means of a piece-wise spline interpolation strategy (spline interpolation for each segment between two minimum-level values). If \f$\theta\f$ values range from 0 to 180 degrees, symmetry is applied. Finally, it is assumed that both E- and H-planes are equal and interpolation is applied for the remaining points.

        \param antenna_angles_deg       Theta angle values in the E-plane in degrees
        \param angles_length            Length of the `antenna_angles_deg` array // TODO confirm
        \param antenna_pattern_dB_in    Antenna pattern values at `antenna_angles_deg` in dB.
        \param pattern_length           Length of the `antenna_pattern_dB_in` array
        \param min_level_dB             Minimum level to split the spline interpolation between contiguous segments.
    */
	void set_antenna_pattern_interp( double* antenna_angles_deg, int angles_length, double* antenna_pattern_dB_in, int pattern_length, double min_level_dB );

    /** \brief Set antenna pattern from single full (for all \f$\theta\f$ angles) full (for \f$\phi\f$=0/90 and \f$\phi\f$=180/270) cuts at both E- and H-plane. Interpolation is applied for the remaining points.

        \param antenna_full_pattern_E_dB_in   Antenna pattern in the E-plane as a function of \f$\theta\f$ angle in dB.
        \param antenna_full_pattern_H_dB_in   Antenna pattern in the H-plane as a function of \f$\theta\f$ angle in dB.
    */
	void set_antenna_patterns_FF( double antenna_full_pattern_E_dB_in[360], double antenna_full_pattern_H_dB_in[360] );

    /** \brief  Set antenna pattern from single full (for all \f$\theta\f$ angles) half (only for \f$\phi\f$=0/90) cuts at both E- and H-plane. Interpolation is applied for the remaining points.

        \param antenna_half_pattern_E_dB_in   Antenna pattern in the E-plane as a function of \f$\theta\f$ angle in dB.
        \param antenna_half_pattern_H_dB_in   Antenna pattern in the H-plane as a function of \f$\theta\f$ angle in dB.
    */
	void set_antenna_patterns_FH( double antenna_half_pattern_E_dB_in[181], double antenna_half_pattern_H_dB_in[181] );

    /** \brief Set the antenna pattern by providing a set of points in both E- and H-planes and a minimum level. The whole shape of both E- and H-plane patterns is then built by means of a piece-wise spline interpolation strategy (spline interpolation for each segment between two minimum-level values). If \f$\theta\f$ values range from 0 to 180 degrees, symmetry is applied. Finally, interpolation is applied for the remaining points.

        \param antenna_angles_E_deg      Theta angle values in the E-plane in degrees
        \param angles_E_length           Number of elements in array `antenna_angles_E_deg`
        \param antenna_pattern_E_dB_in   Antenna pattern values at `antenna_angles_E_deg` in dB
        \param pattern_E_length          Number of elements in array `antenna_pattern_E_dB_in`
        \param antenna_angles_H_deg      Theta angle values in the H-plane in degrees
        \param angles_H_length           Number of elements in array `antenna_angles_H_deg`
        \param antenna_pattern_H_dB_in   Antenna pattern values at `antenna_angles_H_deg` in dB
        \param pattern_H_length          Number of elements in array `antenna_pattern_H_dB_in`
        \param min_level_dB              Minimum level (in dB) to split the spline interpolation between contiguous segments.
    */ 
	void set_antenna_patterns_interp( double* antenna_angles_E_deg, int angles_E_length, double* antenna_pattern_E_dB_in, int pattern_E_length, double* antenna_angles_H_deg, int angles_H_length, double* antenna_pattern_H_dB_in, int pattern_H_length, double min_level_dB );

    /** \brief Get the antenna pattern.

        \param antenna_pattern_dB_out   Output - Antenna pattern as a function of \f$\theta\f$ and \f$\phi\f$ in the antenna frame in dB.
    */
	void get_antenna_whole_pattern( double antenna_pattern_dB_out[181][360] );

    /** \brief Get the E- and H-plane cuts of the antenna pattern.

        \param antenna_pattern_E_dB_out   E-plane cut of the antenna pattern in dB.
        \param antenna_pattern_H_dB_out   H-plane cut of the antenna pattern in dB.
    */
	void get_antenna_patterns( double antenna_pattern_E_dB_out[360], double antenna_pattern_H_dB_out[360] );

    /** \brief Set the main parameters of the receiver (by editing the corresponding private variables).

        \param antenna_Gain_dB_in  Antenna gain in dB.                          
        \param antenna_T_in        Antenna temperature in K.
        \param noise_F_dB_in       Noise figure in dB.
        \param filter_BB_BW_in     Base-band bandwidth of the receiver in Hz.
        \param isotropic_antenna   Indicates if the antenna is isotropic. 
                                   Value | Meaning
                                   :----:|:-------
                                     1   | Antenna is isotropic
                                     0   | Antenna is non-isotropic

                                   When the antenna is isotropic the internally stored antenna pattern is ignored.
    */
	void set_receiver_params( double antenna_Gain_dB_in, double antenna_T_in, double noise_F_dB_in, double filter_BB_BW_in, signed char isotropic_antenna );

    /** \brief Set the antenna effective area.

        \param antenna_Aeff_in    Effective area of the antenna in \f$\mathrm{m}^2\f$
    */
    
	void set_antenna_eff_area( double antenna_Aeff_in );

    /** \brief Set the noise temperature. //TODO check code to see if it is receiver noise temperature or what

        \param noise_T_in Noise temperature in K.
    */
	void set_noise_T( double noise_T_in );

    /** \brief Set the noise power.  // TODO check code to see if it is the receiver noise power at the input and how it is related to the receiver noise temperature

        \param noise_pow_dBW_in    Noise power in dBW
    */
	void set_noise_pow_dBW( double noise_pow_dBW_in );

    /** \brief Set the frequency of the receiver.

        \param frequency_in  Frequency in Hz.
    */ 
	void set_frequency( double frequency_in );

    /** \brief This function is not documented and in the code there is a comment warning that the functions does not work properly. TODO: Review function.
    */
	double get_anglesEH_gain_dB( double angle_E_plane, double angle_H_plane );

    /** \brief Get the antenna gain as a function of \f$\phi\f$ and \f$\theta\f$ angles in the antenna frame.

        \param phi_AntFrame     Angle \f$\phi\f$ in the antenna frame in degrees.
        \param theta_AntFrame   Angle \f$\theta\f$ in the antenna frame in degrees.
        \return                 The antenna gain in this direction, expressed in dB.
    */
	double get_PhiTheta_gain_dB( double phi_AntFrame, double theta_AntFrame );

    /** \brief Get the antenna gain as a function of the incidence or transmitting vector in the receiver body frame. 

        \param incvector    Incidence or transmitting vector in the receiver body frame.
        \return             Antenna gain in this direction, expressed in dB.
    */
	double get_incvector_gain_dB( double incvector[3] );

    /** \brief Get the frequency of the receiver

        \return      The frequency of the receiver in Hz.
    */
	double get_frequency( void );

    /** \brief Get the antenna gain (for isotropic antennas). // TODO This does not make any sense to me, as isotropic antennas have always 0 dB directivity.

        \return    Antenna gain in dB
    */
	double get_antenna_Gain_dB( void );

    /** \brief Get the effective area of the antenna.

        \return Effective area of the antenna in \f$ \mathrm{m}^2\f$.
    */
	double get_antenna_Aeff( void );

    /** \brief  Get the antenna temperature.

        \return   Antenna temperature in K
    */
	double get_antenna_T( void );

    /** \brief Get the noise temperature. // TODO: check if it is the noise temperature of the receiver at the input/output

        /return   Noise temperature in K.
    */
	double get_noise_T( void );

    /** \brief  Get the noise power.

        \return  Noise power in dBW
    */
	double get_noise_pow_dBW( void );

    /** \brief Get the noise figure.

        \return Noise figure in dB  // TODO check if it is only for the receiver or it include the noise figure of the antenna.
    */ 
	double get_noise_F_dB( void );

    /** \brief Get the base-band bandwidth of the receiver.

        \return     Base-band bandwidth of the receiver in Hz.
    */
	double get_filter_BB_BW( void );
	//2D PLANAR ARRAY

    /** \brief Set a distribution of antenna elements over the antenna frame to have a 2D planar array.

        \param element_pos_in    Positions of the array elements as 2D coordinates in the XY plane of the antenna frame in meters or in wavelengths (units of lambda `lambda_units` = 1).
        \param num_elem_in       Length of the array `element_pos_in`.
        \param plane_dim         ? TODO  Clarify what is the meaning of the different input parameters.
        \param lambda_units      Choose if distances are expressed in meters or in wavelengths
                                 Value | Meaning
                                 :----:|:-------
                                   1   | units expressed in wavelengths (lambda)
                                   0   | units expressed in m 
    */
	void set_antenna_elements_pos_AF( double* element_pos_in, int num_elem_in, int plane_dim, char lambda_units );

    /** \brief Set the phases applied to each element to obtain a desired array factor.

        \param phase_delay_in   Phases applied to each element to obtain a desired array factor in radians.
        \param num_elems_in     Number of elements in array `phase_delay_in`
    */
	void set_phase_delays( double* phase_delay_in, int num_elems_in );

    /** \brief Get the phases applied to each element to obtain a desired array factor.

        \param phase_delay_out   Phases applied to each element to obtain a desired array factor in radians.
        \param num_elems_out     Number of elements in array `phase_delay_out`
    */
	void get_phase_delays( double* phase_delay_out, int num_elems_out );

    /** \brief Compute and store the array factor based on the phases and the distribution of the antenna elements.
    */
	void compute_array_factor( void );

    /** \brief Get the array factor stored.

        \param array_factor_dB_out   Array factor as a function of \f$\theta\f$ and \f$\phi\f$ in the antenna frame in dB.  
    */
	void get_array_factor( double array_factor_dB_out[181][360] );

    /** \brief Compute the phases to be applied to each element to get an array factor with a desired pointing direction, based on uniformly-distributed planar array (UPA) theory.

        \param theta_max   Angle \f$\theta\f$ with maximum array factor gain in the antenna frame in degrees.
        \param phi_max     Angle \f$\phi\f$ with maximum array factor gain in the antenna frame in degrees.
    */
	void compute_phase_delays_UPA( double theta_max, double phi_max );

    /** \brief Compute the phases to be applied to each element to get an array factor with a desired pointing direction, based on the ECEF position of a receiver (where the array is placed), a transmitter (where the array is pointing at) and the inertial information of the receiver.

        \param inertials   Inertial rotation of the X-axis (positive clockwise) of the receiver body frame in degrees (first element). Inertial rotation of the Y-axis (positive clockwise) of the receiver body frame in degrees(second element). Inertial rotation of the Z-axis (positive clockwise) of the receiver body frame in degrees (third element).
        \param posR_km     Position of the receiver in ECEF coordinates in km.
        \param posT_km     Position of the transmitter in ECEF coordinates in km.                     
    */ 
	void compute_phase_delays_pos_ECEF_RT( double inertials[3], double posR_km[3], double posT_km[3] );
};

// TODO these functions shall be moved to statuc method functions. Probably private functions
// TODO document these functions
void Set_whole_pattern_from_EH_planes( double ant_pattern_E_dB[360], double ant_pattern_H_dB[360], double ant_pattern_dB[181][360] );
void Set_whole_pattern_from_EH_planes_RevLinInterp( double ant_pattern_E_dB[360], double ant_pattern_H_dB[360], double ant_pattern_dB[181][360] );
void Compute_PhiTheta_from_anglesEH( double angleE, double angleH, double* phi, double* theta );
double Get_gain_pattern( double theta, double phi, double ant_pattern_dB[181][360] );
bool Check_if_UPA_distribution( int &dimX, int &dimY, double &distX, double &distY, int num, double** positions );
