/**
  \file reflecting_surface.h
  \brief Header of the Reflecting_surface class.
*/

using namespace std;

/**
  \brief The Reflecting_surface class allows to define the parameters of a reflecting surface

  This class provides functions to characterize a reflecting surface in a radar scenario, in particular regarding its dielectric properties \cite UlabyI and its roughness.

*/
class Reflecting_surface
{
  private: // variables
	double freq_GHz;                  /**< \brief Frequency of the system in GHz. Default value: 1.57542 (GPS L1). */ 
	double **surface_spectrum;        /**< \brief Spectrum of the surface. 2D array of nx_spec x ny_spec elements. Default value: void. */ 
	double *kx_spec;                  /**< \brief Wavenumber's range of stored spectrum at upwind direction over the reflecting surface in \f$\mathrm{m}^{-1}\f$.array of nx_spec elements.  Default value: void. */
	double *ky_spec;                  /**< \brief Wavenumber's range of stored spectrum at cross-wind direction over the reflecting surface in \f$\mathrm{m}^{-1}\f$.array of ny_spec elements.  Default value: void. */ 
	int nx_spec;                      /**< \brief Number of samples of stored spectrum at upwind direction over the reflecting surface. Default value: 0. */
	int ny_spec;                      /**< \brief Number of samples of stored spectrum at cross-wind direction over the reflecting surface. Default value: 0. */
	double k_threshold;               /**< \brief Wavenumber's limit for the computation of MSS from the stored spectrum. Default value: \f$2\pi/3\lambda_\mathrm{L1}\f$ */
	double **wind_U10_speed_grid;     /**< \brief Grid of wind speeds at 10 meters above the surface in m/sec. 2D array of size_lon_wgrid x size_lat_wgrid elements. Default value: void. */
	double **wind_U10_azimuth_grid;   /**< \brief Grid of wind azimuths at 10 meters above the surface in m/sec. 2D array of size_lon_wgrid x size_lat_wgrid elements. Default value: void. */ 
	double *lon_wgrid;                /**< \brief Longitude values of wind grid in degrees. Array of size_lon_wgrid elements. Default value: void. */
	double *lat_wgrid;                /**< \brief Latitude values of wind grid in degrees. Array of size_lat_wgrid elements. Default value: void. */ 
	int size_lon_wgrid;               /**< \brief Number of samples of longitudes in stored wind grid. Default value: 0. */
	int size_lat_wgrid;               /**< \brief Number of samples of latitudes in stored wind grid. Default value: 0. */
	bool use_wind_grid;               /**< \brief Indicates if a wind grid is being used
                                                  Value | Meaning
                                                  :----:|:--------
                                                   True | A wind grid is stored
                                                  False | There is no wind grid stored
                                      */
  public:
	//Relative permittivity
	double epsilon_real;     /**< \brief Real part of relative permittivity without units. Default value: 73.423        // TODO: Why this value? */
	double epsilon_imag;     /**< \brief Imaginary part of relative permittivity without units. Default value: 56.067   // TODO: Why this value? */
	//Roughness
	double mss_x;            /**< \brief Mean square slope at upwind direction over the reflecting surface without units. Default value: 0.0075  //TODO define what is upwind. Why this value? */
	double mss_y;            /**< \brief Mean square slope at crosswind direction over the reflecting surface without units. Default value: 0.0075 //TODO defube what is crosswind. Why this value? */
	double sigma_z;          /**< \brief Standard deviation of surface height in meters. Default value: 0.069. TODO: Why this value? */
	double c21_coeff;        /**< \brief Gram-Charlier C21 coefficient for 5 m/s wind. Default value: -0.033. TODO: Why this value?  This might only apply to sea-surface*/ 
	double c03_coeff;        /**< \brief Gram-Charlier C03 coefficient for 5 m/s wind. Default value: -0.125. TODO: Why this value? This might only apply to sea-surface*/
	double wind_U10_speed;   /**< \brief Wind speed magnitude at 10 meters above the surface in m/s. Default value: 5.0. TODO: Why this value? This might only apply to sea-surface */
	double wind_U10_azimuth; /**< \brief  Wind azimuth (clockwise starting from North) at 10 m above the surface in degrees. The direction follows the meteorological convention (0 deg means wind coming from North). Default value: 0.0 */
	string medium;           /**< \brief Brief description of the medium. Default value: "Sea water with T=15C and sal=35psu" */

    // Methods
    /** \brief Constructor. Sets default parameters
    */
	Reflecting_surface( void );

    /** \brief Print relevant information about the object's content (public and private variables).
    */
	void dump_parameters( void );

    /** \brief Set the working frequency, expressed in GHz.
    */
	void set_frequency( double frequency_GHz );

    /** \brief Set wavenumber's limit for the computation of MSS from the stored spectrum (\f$2\pi/3\lambda_\mathrm{L1}\f$ by default).
    */
	void set_k_threshold( double k_threshold_in );

    /** \brief Set wavenumber's limit for the computation of MSS from the stored spectrum using \cite Brown78. 
    
        \param incidence_deg  Incidence angle of the incoming signal in degrees.
    */
	void set_k_threshold_Brown( double incidence_deg );

    /** \brief Get wavenumber's limit internally stored

        \return wavenumber's limit internally stored in \f$\mathrm{m}^{-1}\f$
    */ 
	double get_k_threshold( void );

    /** \brief Set relative permittivity for sea water as a function of salinity and temperature. TODO: Say which model is being used (Reference and equations)

        \param salinity      Salinity of sea water in psu.
        \param temperature   Temperature of sea water in Celsius.
    */
	void epsilon_sea_water( double salinity, double temperature );

    /** \brief Set relative permittivity of sea ice as a function of brine.

        \param brine_volume    Brine volume of sea ice in 1/1000 units.
    */
	void epsilon_sea_ice( double brine_volume );

    /** \brief Set relative permittivity of dry snow as a function of snow density.

        \param snow_density  Snow density of dry snow in \f$\mathrm{g}/\mathrm{cm^3}\f$
    */
	void epsilon_dry_snow( double snow_density );

    /** \brief Set relative permittivity of wet snow as a function of snow density and water volume.

        \param snow_density   Snow density of wet snow in \f$\mathrm{g}/\mathrm{cm^3}\f$
        \param water_volume   Water volume of wet snow in \%
    */
	void epsilon_wet_snow( double snow_density, double water_volume );

    /** \brief Compute the reflection Fresnel coefficients for linear polarizations.

        \param incidence_deg     Incidence angle of incoming signal in degrees.
        \param epsilon_upLayer   Complex relative permittivity (real and imaginary parts) of layer above the reflecting surface. For air, simply set epsilon_upLayer = [1.0, 0.0].
        \param rvv               Output. Complex reflection coefficient for vertical polarization (real and imaginary parts).
        \param rhh               Output. Complex reflection coefficient for horizontal polarization (real and imaginary parts).
    */
	void compute_Rfresnel_linear( double incidence_deg, double epsilon_upLayer[2], double rvv[2], double rhh[2] );

    /** \brief Compute the reflection Fresnel coefficients for circular polarizations.

        \param incidence_deg    Incidence angle of incoming signal in degrees.
        \param epsilon_upLayer  Complex relative permittivity (real and imaginary parts) of layer above the reflecting surface. For air, simply set epsilon_upLayer} = [1.0, 0.0].
        \param rco              Output - Complex reflection coefficient for co-polar [RHCP for GPS] (real and imaginary parts).
        \param rcross           Output - Complex reflection coefficient for cross-polar [LHCP for GPS] (real and imaginary parts).
    */
	void compute_Rfresnel_circular( double incidence_deg, double epsilon_upLayer[2], double rco[2], double rcross[2] );

    /** \brief  Compute the transmission Fresnel coefficients for linear polarizations

        \param incidence_deg    Incidence angle of incoming signal in degrees.
        \param epsilon_upLayer  Complex relative permittivity (real and imaginary parts) of layer above the reflecting surface. For air, simply set epsilon_up_layer = [1.0, 0.0].
        \param tvv              Output - Complex transmission coefficient for vertical polarization (real and imaginary parts).
        \param thh              Output - Complex transmission coefficient for horizontal polarization (real and imaginary parts).
    */
	void compute_Tfresnel_linear( double incidence_deg, double epsilon_upLayer[2], double tvv[2], double thh[2] );

    /** \brief Compute the transmission Fresnel coefficients for circular polarizations.

        \param incidence_deg    Incidence angle of incoming signal in degrees
        \param epsilon_upLayer  Complex relative permittivity (real and imaginary parts) of layer above the reflecting surface. For air, simply set epsilon_up_layer = [1.0, 0.0].
        \param tco              Output - Complex transmission coefficient for co-polar [RHCP for GPS] (real and imaginary parts)
        \param tcross           Output - Complex transmission coefficient for cross-polar [LHCP for GPS] (real and imaginary parts).
    */
	void compute_Tfresnel_circular( double incidence_deg, double epsilon_upLayer[2], double tco[2], double tcross[2] );

    /** \brief  Compute the sea spectrum based on \cite Elfouhaily97.

        \param num_samples  Number of samples of the generated spectrum in a single dimension.
        \param delta_k_in   Wavenumber's resolution of the generated spectrum in \f$ \mathrm{m}^{-1}\f$
        \param theta_deg    Wind waves angle in degrees
        \param omega        Wave age
    */
	void compute_sea_spectrum( int num_samples, double delta_k_in, double theta_deg, double omega );

    /** \brief Set sea surface spectrum with a constant wavenumber's resolution.

        \param spectrum_in    Sea surface spectrum. (2D array of MxN elements)
        \param dim_x          Wavenumber's range of spectrum at upwind direction over the reflecting surface in meter$^{-1}$ (a non-repetitive ascending sequence is required for a proper computation of MSS).
        \param dim_y          Wavenumber's range of spectrum at crosswind direction over the reflecting surface in meter$^{-1}$ (a non-repetitive ascending sequence is required for a proper computation of MSS)
        \param kx_spec_in     ? TODO
        \param dim_kx         ? TODO
        \param ky_spec_in     ? TODO
        \param dim_ky         ? TODO
    */
	void set_surf_spectrum( double* spectrum_in, int dim_x, int dim_y, double* kx_spec_in, int dim_kx, double* ky_spec_in, int dim_ky );

    /** \brief Set omnidirectional sea surface with a constant wavenumber's resolution.

        \param spectrum_in_omnidir  Sea surface spectrum. (array of N elements)
        \param dim_x_omnidir        Wavenumber's range of omnidirectional spectrum in \f$ \mathrm{m}^{-1}\f$ (a non-repetitive ascending sequence is required for a proper computation of MSS).
        \param k_spec_in_omnidir    ? TODO
        \param dim_k_omnidir        ? TODO
    */
	void set_surf_spectrum_omnidir( double* spectrum_in_omnidir, int dim_x_omnidir, double* k_spec_in_omnidir, int dim_k_omnidir );

    /** \brief Get value from stored sea surface spectrum as a function of grid's position.

        \param x               Sample index at upwind direction over the reflecting surface.
        \param y               Sample index at crosswind direction over the reflecting surface
        \param kx_ky_spectrum  Output - Wavenumber at sample x, y in \f$ m^{-1}\f$ (frist two elements of the array). The spectrum value at this point (last element of the vector). TODO: Confirm
    */
	void get_surf_spectrum( int x, int y, double kx_ky_spectrum[3] );

    /** \brief Get value from stored sea surface spectrum as a function of array's position.

        \param x             Sample index in the omnidirectional array.
        \param k_spectrum    Output - Spectrum value at sample x in \f$ \mathrm{m}^{-1}\f$ (first element of the array). Spectrum value at that sample (second element of the array).
    */ 
	void get_surf_spectrum_omnidir( int x, double k_spectrum[2] );

    /** \brief Compute MSS from the stored spectrum.
    */
	void compute_mss_from_spectrum( void );

    /** \brief Compute MSS from the wind speed parameters based on \cite Katzberg2006.

    */
	void compute_mss_from_wind( void );

    /** \brief Store an inhomogeneous wind grid as a function of latitude and longitude.

        \param wspeed_grid_in   Grid of wind speeds at 10 meters above the surface in m/sec. (2D array of MxN elements)
        \param dim_wsx          Number of elements if first dimension of wind speed grid (M) (confirm TODO)
        \param dim_wsy          Number of elements in second dimension of wind speed grid (N) (confirm TODO)
        \param wazim_grid_in    Grid of wind azimuths at 10 meters above the surface in degrees. (2D array of MxN elements)
        \param dim_wax          Number of elements in first dimension of wind azimuths grid (M) (confirm TODO)
        \param dim_way          Number of elements in second dimension of wind azimuths grid (N) (confirm TODO)
        \param lon_in           Longitude values of wind grid in degrees.
        \param dim_lon          Number of elements in latitude array (confirm TODO)
        \param lat_in           Latitude values of wind grid in degrees.
        \param dim_lat          Number of elements in latitude array (confirm TODO)
    */
	void set_wind_grid( double* wspeed_grid_in, int dim_wsx, int dim_wsy, double* wazim_grid_in, int dim_wax, int dim_way, double* lon_in, int dim_lon, double* lat_in, int dim_lat );

    /** \brief Interpolate wind speed and wind azimuth from stored wind grid and store the results obtained in public variables wind_U10_speed and wind_U10_azimuth respectively.

        \param lon_in   Longitude coordinate for interpolation of the wind grid in degrees.
        \param lat_in   Latitude coordinate for interpolation of the wind grid in degrees.
    */
	void interp_wind_grid( double lon_in, double lat_in );

    /** \brief Remove stored wind grid
    */
    void disable_wind_grid();

    /** \brief Return current status of stored wind grid.

        \return   Value | Meaning
                  :----:|:-------
                   True | The is a wind grid stored
                  False | There is no wind grid stored
    */
	bool get_wind_grid_status();
};
