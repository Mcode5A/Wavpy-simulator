/** 
    \file gnss_composite.h
    \brief Header file of the GNSS_composite class
*/

/**
    \brief This class generates the autocorrelation function of a set of GNSS signals.
*/
class GNSS_composite
{
    double *lambda_func;   /**< \brief Normalized autocorrelation function stored. Default value: GPS-L1 autocorrelation function with all code components. */   
	double sampling_rate;  /**< \brief Sampling rate of the receiver in samples/sec. Default value: 80000000.0 */  
	double filter_BB_BW;   /**< \brief Base-band bandwidth of the filter applied in Hz. Default value: 12000000.0  */  
	bool lambda_updated;   /**< \brief TODO  */ 
  public:
	int lambda_size;      /**< \brief Number of samples of the autocorrelation function. Default value: 157  */ 
	float weight_CA;      /**< \brief Weight for the C/A code (GPS L1). Default value: 1.0  */ 
	float weight_PY;      /**< \brief Weight for the PY code (GPS L1). Default value: 1.0  */ 
	float weight_M;       /**< \brief Weight for the M code (GPS L1). Default value: 1.0  */ 
	float weight_IM;      /**< \brief Weight for the Inter-Modulation component (GPS L1). Default value: 1.0  */ 
	float weight_E1A;     /**< \brief Weight for the A code (Galileo E1). Default value: 0.0  */ 
	float weight_E1B;     /**< \brief Weight for the B code (Galileo E1). Default value: 0.0  */ 
	float weight_E1C;     /**< \brief Weight for the C code (Galileo E1). Default value: 0.0  */ 
	float weight_B1I;     /**< \brief Weight for the C code (BeiDou B1). Default value: 0.0  */ 
	float weight_L1C;     /**< \brief Weight for the L1C code (QZSS). Default value: 0.0  */ 
    double frequency;     /**< \brief Frequency of the GNSS signal in Hz. Default value: 1575420000.0 (GPS L1)  */ 
    // METHODS
    /** \brief Constructor 
    */
	GNSS_composite( void );

    /** \brief Print relevant information about the object's content (public and private variables).
    */
	void dump_parameters( void );

    /** \brief Set the instrumental parameters that have to be taken into account for the computation of the autocorrelation function.

        \param input_sampling_rate   Input sampling rate in samples/sec.
        \param input_filter_BW       Input base-band bandwidth of the filter applied in Hz.
        \param computeLambda         Choose to compute the autocorrelation function or only set the input parameters
                                     Value | Meaning
                                     :----:|:-------
                                      1    | Compute the autocorrelation function with this function call
                                      0    | Only set the input parameters without computing the autocorrelation function. 
    */
	void set_instrumental_params( double input_sampling_rate, double input_filter_BW, char computeLambda );

    /** \brief Compute and store the autocorrelation function according to current configuration.
    */
	void compute_lambda_func( void );

    /** \brief Provide the range and values of the autocorrelation function internally stored.

        \param range_vec      Output - Range of the autocorrelation function in meters. 
        \param range_length   Number of samples of the range of the autocorrelation function (it has to be equal to `lambda_size`).   
        \param lambda_out     Output - Normalized autocorrelation function. 
        \param lambda_length  Number of samples of the autocorrelation function (it has to be equal to  `lambda_size`).
    */
	void get_lambda_func( double* range_vec, int range_length, double* lambda_out, int lambda_length );

    /** \brief Set a given autocorrelation function.

        \param lambda_in       Normalized autocorrelation function.
        \param lambda_length   The number of elements in array `lambda_in`
    */
	void set_lambda_func( double* lambda_in, int lambda_length );
};



