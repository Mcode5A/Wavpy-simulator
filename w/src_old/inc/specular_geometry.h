/** \file specular_geometry.h 
    \brief Header of the specular_geometry class.
*/


/**
  \brief The specular_gemoetry class to define the geometry of a reflectometry scenario.

  This class provides functions to characterize a reflectometry scenario composed by a transmitter, a receiver and their corresponding specular point over a model of the Earth surface (ellipsoid WGS84 + a given undulation).

  The following definitions apply.

  __Receiver frame:__
  The receiver body frame is a Cartesian coordinates system of the structure containing the receiver, typically a satellite or an aircraft, where the X-axis points towards the front, the Y-axis points towards the right-side (XY define the horizontal plane) and the Z-axis towards Nadir. In abscense of inertial rotation of the body frame, we will assume that the X-axis points towards the Earth's North and the Z-axis points towards the Earth's center.

  __Local frame:__  The local frame is a Cartesian coordinates system with its origin at the specular point, X- and Y-axis defining the horizontal plane parallel to the surface, with the Y-axis pointing towards the transmitter, and the Z-axis pointing to Zenith by complying the right-hand rule.

*/


class Specular_geometry
{
  //IMPORTANT: Distances in km, angles in deg
  private: // Variables
	double posR_ECEF[3];                             /**< \brief Position of the receiver in ECEF coordinates in km. Default value: [6381.137, 0.0, 0.0] */ 
	double posT_ECEF[3];                             /**< \brief Position of the transmitter in ECEF coordinates in km. Default value: [26378.137, 0.0, 0.0] */ 
	double velR_ECEF[3];                             /**< \brief Velocity vector of the receiver in ECEF coordinates in km/s. Default value: [0.0, 0.0, 0.0] */ 
	double velT_ECEF[3];                             /**< \brief Velocity vector of the transmitter in ECEF coordinates in km/s. Default value: [0.0, 0.0, 0.0] */ 
	double posS_ECEF[3];                             /**< \brief Position of the specular point in ECEF coordinates in km. Default value: [6378.137, 0.0, 0.0] */
	double longitudeR;                               /**< \brief Longitude coordinate of the receiver in degrees. Default value: 0.0 */ 
    double latitudeR;                                /**< \brief Latitude coordinate of the receiver in degrees. Default value: 0.0 */ 
	double longitudeT;                               /**< \brief Latitude coordinate of the receiver in degrees. Default value: 0.0 */
    double latitudeT;                                /**< \brief Latitude coordinate of the transmitter in degrees. Default value: 0.0 */
	double heightR;                                  /**< \brief Height of the receiver with respect to ellipsoid WGS84 in km. Default value: 3.0 */ 
    double heightT;                                  /**< \brief Height of the transmitter with respect to ellipsoid WGS84 in km. Default value: 20000.0 */ 
	double local_heightR;                            /**< \brief Height of the receiver in the local frame in km. Default value: 3.0  */ 
    double local_heightT;                            /**< \brief Height of the transmitter in the local frame in km. Default value: 3.0 */ 
	double undulation;                               /**< \brief Vertical height offset of the specular point with respect to ellipsoid WGS84 in km. Default value: 0.0 */ 
	double roll;                                     /**< \brief Inertial rotation of the X-axis (positive clockwise) of the receiver body frame in degrees. Default value: 0.  */ 
    double pitch;                                    /**< \brief Inertial rotation of the Y-axis (positive clockwise) of the receiver body frame in degrees. Default value: 0.0  */ 
    double heading;                                  /**< \brief Inertial rotation of the Z-axis (positive clockwise) of the receiver body frame with respect to North in degrees. Default value: 0.0  */ 
    int nlines_sp3;
    double **xx_sp3, **yy_sp3, **zz_sp3, **vxx_sp3, **vyy_sp3, **vzz_sp3;
	double *tt_sp3;
    int sow_diff_sp3;
    int ref_GPSweek_sp3;
    char gnss_sp3;
    bool sat_vel_sp3;
  public:
	double longitudeS;                                /**< \brief Longitude coordinate of the specular point in degrees. Default value: 0.0 */
    double latitudeS;                                 /**< \brief Latitude coordinate of the specular point in degrees. Default value: 0.0 */
	double elevation;                                 /**< \brief Elevation angle of the transmitter at the specular point in degrees. Default value: 90.0  */ 
	double azimuthR;                                  /**< \brief Azimuth angle of the receiver at the specular point in degrees. Default value: 0.0  */ 
    double azimuthT;                                  /**< \brief Azimuth angle of the transmitter at the specular point in degrees. Default value: 0.0  */ 
	double geometric_delay;                           /**< \brief Delay path difference between transmitter-specular-receiver and transmitter-receiver in km. Default value: 6.0  */ 
	Specular_geometry( void );                        /**< \brief Constructor  */ 
	void dump_parameters( void );                     /**< \brief Print relevant information about the object's content (public and private variables). TODO: Why should we dump the contents of private variables?  */ 

    /** \brief Set the position of the receiver in ECEF coordinates. 
               The orientation of the Z-axis of the receiver body frame with respect to North (__heading__) is updated with the receiver's position and flight direction.  

        \param posR_in   Position of the receiver in ECEF coordinates in km.
    */  
	void set_ECEFpos_Receiver( double posR_in[3] );   

    /** \brief Get the position of the receiver in ECEF coordinates.

        \param posR_out   Output - Position of the receiver in ECEF coordinates in km
   */ 
	void get_ECEFpos_Receiver( double posR_out[3] );  

    /** \brief Set the velocity vector of the receiver in ECEF coordinates. The orientation of the Z-axis of the receiver body frame with respect to North (_heading_) is updated with the receiver's position and flight direction.

        \param velR_in    Velocity vector of the receiver in ECEF coordinates in km/s.
    */
	void set_ECEFvel_Receiver( double velR_in[3] );

    /** \brief Get the position of the receiver in ECEF coordinates.

        \param velR_out   Output -  Velocity vector of the receiver in ECEF coordinates in km/s.
    */
	void get_ECEFvel_Receiver( double velR_out[3] );

    /** \brief Set the position of the transmitter in ECEF coordinates.
    
        \param posT_in    Position of the transmitter in ECEF coordinates in km.
    */
	void set_ECEFpos_Transmitter( double posT_in[3] );

    /** \brief Get the position of the transmitter in ECEF coordinates

        \param posT_out   Position of the transmitter in ECEF coordinates in km.
    */
	void get_ECEFpos_Transmitter( double posT_out[3] );

     /** \brief Set the velocity vector of the transmitter in ECEF coordinates.

         \param velT_in   Velocity vector of the transmitter in ECEF coordinates in km/s.
    */ 
	void set_ECEFvel_Transmitter( double velT_in[3] );

    /** \brief  Get the position of the transmitter in ECEF coordinates.

        \param velT_out   Velocity vector of the transmitter in ECEF coordinates in km/s.
    */ 
	void get_ECEFvel_Transmitter( double velT_out[3] );

    /** \brief Get the position of the specular point in ECEF coordinates.

        \param posS_out   Position of the specular point in ECEF coordinates in km.
    */ 
	void get_ECEFpos_Specular( double posS_out[3] );

    /** \brief Set the position of the receiver with Longitude-Latitude-Height coordinates. The orientation of the Z-axis of the receiver body frame with respect to North (_heading_) is updated with the receiver's position and flight direction.

        \param LonLatHeight_R_in   Longitude coordinate of the receiver in degrees (first element). Latitude coordinate of the receiver in degrees (second element). Height of the receiver with respect to ellipsoid WGS84 in km (third element). 
    */
	void set_LongLatHeight_Receiver( double LonLatHeight_R_in[3] );

    /** \brief Get the position of the receiver with Longitude-Latitude-Height coordinates.

        \param LonLatHeight_R_out  Longitude coordinate of the receiver in degrees (first element). Latitude coordinate of the receiver in degrees (second element). Height of the receiver with respect to ellipsoid WGS84 in km (third element).
    */
	void get_LongLatHeight_Receiver( double LonLatHeight_R_out[3] );

    /** \brief Set the position of the transmitter with Longitude-Latitude-Height coordinates.

        \param LonLatHeight_T_in   Longitude coordinate of the transmitter in degrees (first element).  Latitude coordinate of the transmitter in degrees (second element).  Height of the transmitter with respect to ellipsoid WGS84 in km (third element).
    */
	void set_LongLatHeight_Transmitter( double LonLatHeight_T_in[3] );

    /** \brief Get the position of the transmitter with Longitude-Latitude-Height coordinates.

        \param LonLatHeight_T_out   Longitude coordinate of the transmitter in degrees (first parameter).  Latitude coordinate of the transmitter in degrees (second parameter).  Height of the transmitter with respect to ellipsoid WGS84 in km (third parameter).
    */
	void get_LongLatHeight_Transmitter( double LonLatHeight_T_out[3] );

    /** \brief Set the geometry of the different elements from their heights, the elevation and azimuth angles, and the location of the specular point in Longitude-Latitude coordinates. The orientation of the Z-axis of the receiver body frame with respect to North (_heading_) is updated with the receiver's position and flight direction.

        \param elev_in       Elevation angle of the transmitter at the specular point in degrees.
        \param heightR_in    Height of the receiver with respect to ellipsoid WGS84 in km.
        \param heightT_in    Height of the transmitter with respect to ellipsoid WGS84 in km.
        \param lonS_in       Longitude coordinate of the specular point in degrees.
        \param latS_in       Latitude coordinate of the specular point in degrees.
        \param azimT_in      Azimuth angle of the transmitter at the specular point in degrees.
        \param heightS_in    Vertical height offset of the specular point with respect to ellipsoid WGS84 in km. 
        \param computeUndu   TODO
    */
	void set_geometry_from_ElevHeightsSpec( double elev_in, double heightR_in, double heightT_in, double lonS_in, double latS_in, double azimT_in, double heightS_in, char computeUndu );

    /** \brief  Set a velocity vector for the receiver tangential to the Earth surface. The orientation of the Z-axis of the receiver body frame with respect to North (__heading__) is updated with the receiver's position and flight direction.

        \param velocity       Speed of the receiver in km/s.
        \param specAzim_deg   Clockwise azimuth angle with respect to the pointing direction towards the specular point in degrees.
    */
	void set_tangEarthVel_Receiver( double velocity, double specAzim_deg );

    /** \brief  Set a velocity vector for the transmitter tangential to the Earth surface.

        \param velocity       Speed of the transmitter in km/s.
        \param specAzim_deg   Clockwise azimuth angle with respect to the pointing direction towards the specular point in degrees.
    */ 
	void set_tangEarthVel_Transmitter( double velocity, double specAzim_deg );

    /** \brief Set a vertical height offset of the specular point with respect to ellipsoid WGS84.

        \param undu    Vertical height offset of the specular point with respect to ellipsoid WGS84 in km.
    */
	void set_Undulation( double undu );

    /** \brief Get the vertical height offset of the specular point with respect to ellipsoid WGS84.

        \return    Vertical height offset of the specular point with respect to ellipsoid WGS84 in km.
    */
	double get_Undulation( void );

    /** \brief Set the ECEF position and velocity of the receiver from an ASCII file of five columns containing the following variables: GPS_week - Second_of_week - Position_X_km - Position_Y_km - Position_Z_km. Interpolation is applied when required. The orientation of the Z-axis of the receiver body frame with respect to North (__heading__) is updated with the receiver's position and flight direction.

        \param namefile   Filename of the ASCII file containing the time series of the receiver's ECEF position.
        \param week       GPS week.
        \param sow        GPS second of week.
    */
	void read_ECEFpos_Receiver( const char* namefile, int week, double sow );

    /** \brief  Set the ECEF position and velocity of the transmitter from an ASCII file of five columns containing the following variables: GPS\_week - Second\_of\_week - Position\_X\_km - Position\_Y\_km - Position\_Z\_km. Interpolation is applied when required.

        \param namefile   Filename of the ASCII file containing the time series of the transmitter's ECEF position.
        \param week       GPS week.
        \param sow        GPS second of the week.
   */ 
	void read_ECEFpos_Transmitter( const char* namefile, int week, double sow );

    /** \brief Set the ECEF position and velocity of a GNSS transmitter from a SP3 file. If several reads have to be done at the same SP3 file, it is better to use _load_sp3File_, _read_ECEFpos_GNSS_Transmitter_sp3Loaded_ and _free_sp3File_.

        \param namefile           Filename of the SP3 file containing ECEF positions of several GNSS satellites.
        \param week               GPS week.
        \param sow                GPS second of the week.
        \param prn                PRN of the desired GNSS satellite.
        \param gnss_identifier    GNSS identifier 
                                  Identifier | Constellation
                                  :---------:|:-------------
                                     'G'     |  GPS
                                     'E'     |  Galileo
                                     'C'     |  Beidou
                                     'R'     |  GLONASS
                                     'J'     |  QZSS          
    */
	void read_ECEFpos_GNSS_Transmitter( const char* namefile, int week, double sow, int prn, char gnss_identifier );

    /** \brief Load the orbits of a single GNSS from a SP3 file.

        \param namefile           Filename of the SP3 file containing ECEF positions of several GNSS satellites.
        \param gnss_identifier    GNSS identifier 
                                  Identifier | Constellation
                                  :---------:|:-------------
                                     'G'     |  GPS
                                     'E'     |  Galileo
                                     'C'     |  Beidou
                                     'R'     |  GLONASS
                                     'J'     |  QZSS          
    */
    void load_sp3File( const char* namefile, char gnss_identifier );

    /** \brief Free from memory previous loaded orbits of a single GNSS from a SP3 file.
    */
    void free_sp3File( void );

    /** \brief Set the ECEF position and velocity of a GNSS transmitter from a previously loaded SP3 file.

        \param week               GPS week.
        \param sow                GPS second of the week.
        \param prn                PRN of the desired GNSS satellite.
        \param gnss_identifier    GNSS identifier 
                                  Identifier | Constellation
                                  :---------:|:-------------
                                     'G'     |  GPS
                                     'E'     |  Galileo
                                     'C'     |  Beidou
                                     'R'     |  GLONASS
                                     'J'     |  QZSS          
    */
    void read_ECEFpos_GNSS_Transmitter_sp3Loaded( int week, double sow, int prn, char gnss_identifier );

    /** \brief Compute the specular point from the positions of transmitter and receiver over the ellipsoid WGS84 plus a given undulation by increasing both semi-axis with such value.

        \param computeUndu    Indicate if undulation of EGM96 has to be computed
                              Value | Meaning
                              :----:|:-------
                               1    | interpolate undulation from EGM96 stored grid
                               0    | do not compute undulation.
    */
	void compute_specular_point( char computeUndu );

    /** \brief Compute the specular point from the positions of transmitter and receiver over the ellipsoid WGS84 plus a given undulation by increasing the radius of curvature at initial specular position over WGS84.
    */
	void compute_specular_point_Undu_Spherical_Earth( void );

    /** \brief  Compute elevation and azimuth angles of the transmitter as seen from the receiver's point of view.

        \param elevAzimT_out    Output - Elevation angle of the transmitter from the receiver's point of view in degrees (first element). Azimuth angle of the transmitter from the receiver's point of view in degree (second element).
    */ 
	void compute_ElevAzimT_from_receiver( double elevAzimT_out[2] );

    /** \brief Set the inertial rotation of the receiver, including its orientation with respect to North.

        \param roll_in      Inertial rotation of the X-axis (positive clockwise) of the receiver body frame in degrees.
        \param pitch_in     Inertial rotation of the Y-axis (positive clockwise) of the receiver body frame in degrees.
        \param heading_in   Inertial rotation of the Z-axis (positive clockwise) of the receiver body frame with respect to North in degrees.
    */
	void set_inertials( double roll_in, double pitch_in, double heading_in );

    /** \brief Get the stored inertial rotation of the receiver.

        \param vector_RPY_out   Output - Inertial rotation of the X-axis (positive clockwise) of the receiver body frame in degrees (first element). Inertial rotation of the Y-axis (positive clockwise) of the receiver body frame in degrees (second element). Inertial rotation of the Z-axis (positive clockwise) of the receiver body frame with respect to North in degrees (third element).
    */
	void get_inertials( double vector_RPY_out[3] );

    /** \brief Rotate an input vector in the receiver body frame to the local frame's orientation (it is not a change of coordinates).

        \param  vector_BF_in       Input vector in the receiver body frame.
        \param  vector_local_out   Output vector with the local frame's orientation.
    */
	void rotate_vector_BF_to_local( double vector_BF_in[3], double vector_local_out[3] );

    /** \brief Rotate an input vector in the receiver body frame to ECEF orientation (it is not a change of coordinates).

        \param vector_BF_in      Input vector in the receiver body frame
        \param vector_ECEF_out    Output vector with ECEF orientation.
    */
	void rotate_vector_BF_to_ECEF( double vector_BF_in[3], double vector_ECEF_out[3] );

    /** \brief Compute the projection of an input vector in the receiver body frame into the reflected signal's delay path.

        \param vector_BF_in   Input vector in the receiver body frame. As output: Projection of the input vector into the reflected signal's delay path in the same units as vector_BF_in.
    */
	double compute_inertial_delay( double vector_BF_in[3] );

    /** \brief Set the inertial rotation of the receiver from an ASCII file of five columns containing the following variables: GPS_week - Second_of_week - roll_deg - pitch_deg - yaw_deg. Interpolation is applied when required.

        \param namefile   Filename of the ASCII file containing the time series of the receiver's inertial rotation.
        \param week       GPS week.
        \param sow        GPS second of the week.
    */
	void read_Inertials_Receiver( const char* namefile, int week, double sow );

    /** \brief Compute the carrier phase wind-up of the up-looking antenna (collecting direct GNSS signals) based on \cite Beyerle2009.

        \param vector_r_a_BF      Antenna vector r_a (as in \cite Beyerle2009) in the receiver body frame.
        \param vector_r_t_BF      Antenna vector r_t (as in \cite Beyerle2009) in the receiver body frame.
        \param week               GPS week.
        \param sow                GPS second of the week.
        \param windup_phase_R_L   Output - Carrier phase wind-up for RHCP polarization in radians (first element). Carrier phase wind-up for LHCP polarization in radians (second element).
    */
	void compute_Beyerle_windup_direct( double vector_r_a_BF[3], double vector_r_t_BF[3], int week, double sow, double windup_phase_R_L[2] );

    /** \brief Compute the carrier phase wind-up of the down-looking antenna (collecting direct GNSS signals) based on \cite Beyerle2009.

        \param vector_r_a_BF      Antenna vector r_a (as in \cite Beyerle2009) in the receiver body frame.
        \param vector_r_t_BF      Antenna vector r_t (as in \cite Beyerle2009) in the receiver body frame.
        \param rvv                Complex reflection coefficient for vertical polarization (real and imaginary parts).
        \param rhh                Complex reflection coefficient for horizontal polarization (real and imaginary parts).
        \param week               GPS week.
        \param sow                GPS second of the week.
        \param windup_phase_R_L   Output - Carrier phase wind-up for RHCP polarization in radians (first element). Carrier phase wind-up for LHCP polarization in radians (second element).
    */
	void compute_Beyerle_windup_reflected( double vector_r_a_BF[3], double vector_r_t_BF[3], double rvv[2], double rhh[2], int week, double sow, double windup_phase_R_L[2] );
};
