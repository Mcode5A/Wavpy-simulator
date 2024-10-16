/**
    \file MultiRaySingleRefl_model.h 
    \brief Header of the MRSR_Model class
*/

/**
    \brief This class models GNSS-R complex waveforms in a scenario with several reflecting layers. 

    This model, referred as Multiple Ray - Single Reflection (MRSR), constructs the reflected signal as a complex sum of single coherent reflections coming from a set of layer interfaces (one reflection from each layer). Such approach was successfully tested in an Antarctic campaign \cite Cardellach12.
    
    Additional information can be found below.

    __Surface frame:__
    The surface frame is a Cartesian coordinates system with its origin at the vertical projection of the receiver's location towards the surface level, X- and Y-axis defining the horizontal plane parallel to the surface, with the X-axis pointing towards North, and the Z-axis pointing to Zenith by complying the right-hand rule.
*/
class MRSR_Model
{
	int num_layers;                       /**< \brief Number of layers. Default value: 0  */ 
	//XYZ system centered at the surface level (plane XY) under the receiver and Z positive towards it
	double height_z;                      /**< \brief Height of the receiver with respect of the surface level (or its position in the Z-axis of the surface frame) in meters. Default value: 0   */ 
	double *depth_layer;                  /**< \brief Initial depth level of the different layers at the Z-axis of the surface frame in meters. Note that positive values refer to negative values in the Z-axis' positions (e.g. 10 meters depth means -10 meters in the Z-axis of the surface frame). Default value: void  */ 
	double *alpha_x;                      /**< \brief Rotation angle in the X-axis of the surface frame of the different layers in degrees (following right-hand rule with the surface frame) (pitch). Default value: void  */ 
	double *alpha_y;                      /**< \brief Rotation angle in the Y-axis of the surface frame of the different layers in degrees (following right-hand rule with the surface frame) (roll). Default value: void  */  
	double *epsilon_r;                    /**< \brief Real part of relative permittivity of the different layers without units. Default value: void  */ 
	double *epsilon_i;                    /**< \brief Imaginary part of relative permittivity of the different layers without units. Default value: void  */ 
  public:
	RF_FrontEnd receiver;                 /**< \brief  Object characterizing the receiver front-end (collecting reflected GNSS signals) used during the simulation. */ 
	GNSS_composite gnss_signal;           /**< \brief  Object characterizing the GNSS signal used during the simulation.  */ 
	Waveform_complex_cluster waveforms;   /**< \brief  Object that stores the simulated cluster of complex waveforms with functions for data analysis. */ 

	/** \brief Constructor
    */
    MRSR_Model( void );
    /** \brief Characterize a multiple-layer scenario. 

        \param height_in       Height of the receiver with respect of the surface level (or its position in the Z-axis of the surface frame) in meters.
        \param depths_in       Initial depth level of the different layers at the Z-axis of the surface frame in meters.
        \param num_depths      Length of array `depths_in`
        \param alpha_x_in      Rotation angle in the X-axis of the surface frame of the different layers in degrees (following right-hand rule with the surface frame).
        \param num_alpha_x     Length of array `alpha_x_in`
        \param alpha_y_in      Rotation angle in the Y-axis of the surface frame of the different layers in degrees (following right-hand rule with the surface frame).
        \param num_alpha_y     Length of array `alpha_y_in`
        \param epsilon_r_in    Real part of relative permittivity of the different layers without units.
        \param num_epsilon_r   Length of array `epsilon_r_in`
        \param epsilon_i_in    Imaginary part of relative permittivity of the different layers without units.
        \param num_epsilon_i   Length of array `epsilon_i_in`
    */
	void set_general_scenario( double height_in, double* depths_in, int num_depths, double* alpha_x_in, int num_alpha_x, double* alpha_y_in, int num_alpha_y, double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i );

    /** \brief Characterize a multiple-layer scenario where all layers are parallel to the XY-plane of the surface frame.

        \param height_in      Height of the receiver with respect of the surface level (or its position in the Z-axis of the surface frame) in meters. 
        \param depths_in      Initial depth level of the different layers at the Z-axis of the surface frame in meters.
        \param num_depths     Length of array `depths_in`
        \param epsilon_r_in   Real part of relative permittivity of the different layers without units.
        \param num_epsilon_r  Length of array `epsilon_r_in`
        \param epsilon_i_in   Imaginary part of relative permittivity of the different layers without units.
        \param num_epsilon_i  Length of array `epsilon_i_in`
    */
	void set_planar_layers_scenario( double height_in, double* depths_in, int num_depths, double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i );

    /** \brief Characterize a dry snow multiple-layer scenario where all layers are parallel to the XY-plane of the surface frame.

        \param height_in       Height of the receiver with respect of the surface level (or its position in the Z-axis of the surface frame) in meters.
        \param depths_in       Initial depth level of the different layers at the Z-axis of the surface frame in meters.
        \param num_depths      Length of array `depths_in` 
        \param snow_dens_in    Snow density of the different layers in gr/cm$^3$ units.
        \param num_snow_dens   Length of array `num_snow_dens`
    */
	void set_dry_snow_planar_layers_scenario( double height_in, double* depths_in, int num_depths, double* snow_dens_in, int num_snow_dens );

    /** \brief Modify the height of the receiver and the depths of the layers in the already characterized scenario.

        \param height_in    Height of the receiver with respect of the surface level (or its position in the Z-axis of the surface frame) in meters.       
        \param depths_in    Initial depth level of the different layers at the Z-axis of the surface frame in meters.
        \param num_depths   Length of array `depths_in`
    */
	void mod_height_depths( double height_in, double* depths_in, int num_depths );

    /** \brief Modify the rotation angles of the layers in the already characterized scenario.

        \param alpha_x_in    Rotation angle in the X-axis of the surface frame of the different layers in degrees (following right-hand rule with the surface frame).
        \param num_alpha_x   Length of array `alpha_x_in` 
        \param alpha_y_in    Rotation angle in the Y-axis of the surface frame of the different layers in degrees (following right-hand rule with the surface frame).
        \param num_alpha_y   Length of array `alpha_y_in`
    */
	void mod_alphas( double* alpha_x_in, int num_alpha_x, double* alpha_y_in, int num_alpha_y );

    /** \brief Modify the permittivity of the layers in the already characterized scenario.

        \param epsilon_r_in     Real part of relative permittivity of the different layers without units.
        \param num_epsilon_r    Length of array `epsilon_r_in` 
        \param epsilon_i_in     Imaginary part of relative permittivity of the different layers without units.
        \param num_epsilon_i    Length of array `epsilon_i_in`
    */
	void mod_epsilon( double* epsilon_r_in, int num_epsilon_r, double* epsilon_i_in, int num_epsilon_i );

    /** \brief Simulate a cluster of complex waveforms for the characterized scenario.

        \param wav_lags         Size of the waveforms stored in the simulated cluster.
        \param lag_direct_pos   Lag position corresponding to the direct signal's peak in the simulated waveform's range.
        \param sampling_rate    Sampling rate of the waveforms from the simulated cluster in samples/sec.
        \param elevations       Elevation of the incident GNSS signal at the surface level for each sample of the simulated cluster in degrees.
        \param size_elevs       Azimuth of the incident GNSS signal at the surface level for each sample of the simulated cluster in degrees.
        \param yaws             ? TODO
        \param size_yaws        ? TODO
    */
	void compute_GNSS_wavcluster( int wav_lags, int lag_direct_pos, double sampling_rate, double* elevations, int size_elevs, double* yaws, int size_yaws );

    /** \brief Compute the relationship between interferometric frequency and depth in the resultant lag-hologram (obtained from a simulated waveforms cluster) for the characterized scenario. 

        \param elev_range           Elevation range (initial and final values) of the incident GNSS signal at the surface level for the simulated cluster in degrees.
        \param azim_range           Azimuth range (initial and final values) of the incident GNSS signal at the surface level for the simulated cluster in degrees.
        \param time_range           Time range (initial and final values) of the simulated cluster in seconds.
        \param freq_LH              Interferometric frequencies (in cycles/degree of elevation units) of the resultant lag-hologram (obtained from a simulated waveforms cluster) for the characterized scenario.
        \param samples_freq_LH      Number of frequency samples of the resultant lag-hologram.                                                
        \param depth_LH             Depths (in meters) of the resultant lag-hologram (obtained from a simulated waveforms cluster) for the characterized scenario.
        \param samples_depth_LH     Number of depths samples of the resultant lag-hologram (it has to be equal than `samples_freq_LH`). 
    */
	void compute_LH_freqs_and_depths( double elev_range[2], double azim_range[2], double time_range[2], double* freq_LH, int samples_freq_LH, double* depth_LH, int samples_depth_LH );

    /** \brief Compute the relative received power (not waveforms) assuming an incoming signal (with power~=~1 at the surface level) at a given frequency and linear polarizations for the characterized scenario.

        \param elevation  Elevation of the incident signal at the surface level in degrees.
        \param yaw        Azimuth of the incident signal at the surface level in degrees.
        \param freq       Frequency of the incident signal in Hz.
        \param pow_HV     Output - Relative received power at horizontal polarization (first element). Relative received power at vertical polarization (second elements).
    */
	void compute_pow_linearPol( double elevation, double yaw, double freq, double pow_HV[2] );
};

int compute_elevs_out( double elev_in, double *elevs_out, bool *valid_elevs_out, double *depth_layer, double *alpha_plane, double *epsilon_r, double *epsilon_i, int num_layers );
double theta_layer2_Snell( double theta_layer1, double eps_r_layer1, double eps_i_layer1, double eps_r_layer2, double eps_i_layer2 );
int get_layered_powrange( double elev_in, double height_rec, double lambda, double *depth_layer, double *alpha_plane, double *epsilon_r, double *epsilon_i, double *elevs_out, bool *valid_elevs_out, int num_layers, double *range_out, double *amp_pol1, double *amp_pol2, bool lin_pol );
bool traverse_single_layer( double theta_normal, double orig_thickness, double alpha_layer, double alpha_beneath_layer, double epsilon_r, double epsilon_i, double lambda, bool up2down, double &horizontal_pos, double &delta_range, double &amp_atten );
