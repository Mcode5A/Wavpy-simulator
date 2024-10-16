#ifndef COMPLEXSIGNAL_HPP_
#define COMPLEXSIGNAL_HPP_

/**
    \file   ComplexSignal.hpp
    \brief  Header file for the ComplexSignal class.
*/

/**
    \brief  The ComplexSignal class allows to easily work with complex signals.

    The class allows to define commonly used complex signals, such as phasors, delta function, pulses, 
    GNSS PRN signals, etc.,
    and manipulated them in several ways like amplifying, adding, multiplying or computing the cross-correlation
    between different signals.
    It is also possible to print or plot the signals.
*/

#include <iostream>
#include <string>
#include <gsl/gsl_complex.h>

#include "Observable.hpp"
#include "ComplexWaveform.hpp"

class ComplexSignal : public Observable
{
  // Member variables
#ifndef TEST_COMPLEXSIGNAL
  private:
#else  
  public:
#endif
    unsigned int  mLength;        // number of samples
    double*       mReal=NULL;     // real part of the signal
    double*       mImag=NULL;     // imaginary part of the signal
    double        mFs;            // sampling frequency
    long int      mInitialSample; // initial sample respect to origin

  // Public member functions
  public:
    // Constructors and destructor
    ComplexSignal(unsigned int length=0, double Fs=80e6, long int initial_sample=0); //
    ComplexSignal(unsigned int length, double* real, double* imag, double Fs=80e6, long int initial_sample=0); //
    ComplexSignal(const ComplexSignal& cs); ///
    ~ComplexSignal(); ///

    // Set things
    void SetFs(double Fs); ///
    void SetInitialSample(long int initial_sample); ///
    void SetComplex(unsigned int length, double* real, double* imag); ///
    void SetReal(unsigned int length, double* real); ///
    void SetImag(unsigned int length, double* imag); ///
    void SetSample(unsigned int pos, double vreal, double vimag); ///

    // Get things
    unsigned int GetLength()        const; ///
    double       GetFs()            const; ///
    int          GetInitialSample() const; ///

    // Operators
    bool operator== (const ComplexSignal& s) const;  ///
    bool operator!= (const ComplexSignal& s) const;  ///
    ComplexSignal& operator= (const ComplexSignal& s); ///
    ComplexSignal operator+ (const ComplexSignal& s) const; ///
    ComplexSignal operator+ (gsl_complex v) const; ///
    ComplexSignal operator- (const ComplexSignal& s) const; ///
    ComplexSignal operator- (gsl_complex v) const; ///
    ComplexSignal operator* (const ComplexSignal& s) const; ///
    ComplexSignal operator* (double a) const; ///
    ComplexSignal operator/ (const ComplexSignal& s) const; ///
    ComplexSignal operator/ (double a) const; ///
    ComplexSignal Conj(void) const; ///
    friend std::ostream& operator<< (std::ostream &out, const ComplexSignal& s); ///

    // Manipulate Signals
    ComplexSignal& Crop(unsigned int n0, unsigned int n1) const;

    // Estimators
    gsl_complex Mean() const;
    double      Var() const;
    double      Std() const;

    // Special Signals

    /**
      Generates a phasor

      A phasor of frequency f  is described with 
      \f[
        s(t) = e^{j 2 \pi f t + \phi_0}
      \f]
      This function generates a signal that is the sampled version of a phasor:
      \f[
        s[n] = s(n T_s) = e^{ j 2 \pi f n / F_s + \phi_0}
      \f]
      where 
      \f[
        T_s = 1/F_s
      \f]
      @param Freq     the frequency of the phasor in [Hz]
      @param Fs       the sampling frequency in [Hz]
      @param length   the length of the signal in number of samples
      @param Phi0     the initial phase of the phasor in radians
      @bool parallel  Set to 'true' to use paralle code, to 'false' for sequential code execution
      @return         the pointer to the phasor
    */
    static ComplexSignal Phasor(double Freq, double Fs, unsigned int length, double Phi0=0, bool parallel=false); ///

    /**
      Generates Additive White Gaussian Noise (AWGN)
      @param sigma   the noise's standard deviation for each (in-phase and quadrature) component
      @param length  the length of the signal in number of samples
      @param Fs      the sampling frequency in [Hz]
      @param seed    the seed to be used to generate the AWGN. If no argument is passed, the current
                     seed of the random number generator is used.
      @return        The pointer to the noise signal
    */
    static ComplexSignal AWGN(double sigma, unsigned int length, double Fs, int seed=-1); ///

    /**
      Generates a GPS CA code at baseband
      @param PRN    The PRN number
      @param Fs     The sampling frequency
      @param length The length (number of samples) of the signal to be generated
    */
    ComplexSignal GPS_CA(unsigned char PRN, double Fs, unsigned int length);

    //-----------------------------------------------------
    // Functions related to CA code generation

    /**
      Generates a sequence of GPS CA code chipping sequences
      @param ca    Is a return value of an array of integers (+1, -1) representing the chipping sequence
      @param PRN   The PRN number
      @return      1 on success, 0 otherwise (NULL ca pointer, or wrong PRN number)
    */
    static int CACode(int *ca_code, unsigned int PRN);


  //-------------------------------------------------------
  // Functions related to correlations
#ifndef TEST_COMPLEXSIGNAL
  private:
#else
  public:
#endif
    /**
      Checks that the cross-correlation limits are correct so that arrays are not accessed outside their range.

      @param length_x  Length of the first signal
      @param length_y  Length of the last signal
      @param N         The number of signal samples to be integrated during the cross-correlation.
      @param n0        The initial saple from which the cross-correlation has to be compted.
      @param m_min     The minimum lag at which the cross-correlation shall be computed.
      @param m_max     The maximum lag at which the cross-correlatio shall be computed.
      @return          The real correlation value.
    */
    static bool CorrelationLimitsCorrect(unsigned int length_x, unsigned int length_y, unsigned int N, unsigned int n0, int min_m, int max_m); //

    /**
      Computes the cross-correlation of two real signals at a single lag.
      @param x        The first input array.
      @param y        The second input array.
      @param N        The number of signal samples to be integrated during the cross-correlation.
      @param n0       The intial sample from which the cross-correlation has to be computed.
      @param m        The lag position at which the cross-correlation is computed.
      @param parallel Set to 1 if parallel code execution is desired. 0 otherwise.
      @return         The complex value of the cross-correlation.
    */
    static double CorrelationRealSingleLag(const double* x, const double* y, unsigned int N, unsigned int n0, int m, bool parallel=0); //

    /**
      Computes the cross-correlation of two complex signals at a single lag.
     
      @param x        The first input signal.
      @param y        The second input signal.
      @param N        The number of signal samples to be integrated during the cross-correlation.
      @param n0       The intial sample from which the cross-correlation has to be computed.
      @param m        The lag position at which the cross-correlation is computed.
      @param parallel Set to 1 if parallel code execution is desired. 0 otherwise.
      @return         The complex value of the cross-correlation.
    */
    static gsl_complex CorrelationComplexSingleLag(const ComplexSignal& x, const ComplexSignal& y, unsigned int N, unsigned int n0, int m, bool parallel=0); //

  public:
    /**
      Computes the cross-correlation of two signals. The computation is carried out in the time domain.
 
      @param x       the first input signal
      @param y       the second input signal
      @param N       the number of samples to be integrated
      @param n0      the initial sample at which the correlation integral is computed
      @param min_m   the minimum lag at which the correlation is computed
      @param step_m  the separation between lags (in samples)
      @param max_m   the maximum lag at which the correlation is computed
      @parallel      to chose for parallel computation
      @return        a complex signal that represent the cross-correlation function of the input signals
    */
    static ComplexWaveform Correlation(const ComplexSignal& x, const ComplexSignal& y, unsigned int N, unsigned int n0, int min_m, unsigned int step_m, int max_m, bool parallel=false);

   
    // ----------------------------------------------------
    // Load BIBASPIR files
#ifndef TEST_COMPLEXSIGNAL
  private:
#else
  public:
#endif


    static bool IsBIBASPIRFlag(const int* bibaspirdata);
    /**
      Checks the IDs of all blocks within a BIBA-SPIR file.

      @param bibaspirdata   Memory buffer holding a complete BIBA-SPIR file
      @return               Return true if all blocks are correct, false otherwise.
      @todo                 Check if this is true for the input parameter.
    */
    static bool CheckBIBASPIRBlocks(int* bibaspirdata);
  
    /**
      Extract one real signal sample (UP/DW, L1/L5, I/Q) from a four-bytes word that corresponds to a sampling instant in a BIBA-SPIR file

      @param word     the 32 bit word that corresponds to a sampling interval
      @param UP_nDW   Select *up* or *down* channels
                          UP_nDW  | selection
                          :------:|:--------:
                             0    |   down
                             1    |   up
      @param L1_nL5   Select *L1* or *L5* frequency band
                          L1_nL5  | selection
                          :------:|:--------:
                             0    |   L5
                             1    |   L1
      @param I_nQ     Select the *in-phase* or *quadrature* component
                            I_nQ  | selection
                          :------:|:------------:
                             0    |   quadrature
                             1    |   in-phase
      @return         the selected real signal sample
    */
    static double ExtractBIBASPIRSample(int word, bool UP_nDW, bool L1_nL5, bool I_nQ);

    /**
      Extract one in-phase signal sample (UP/DW, L1/L5, I/Q) from a four-bytes word that corresponds to a sampling instant in a BIBA-SPIR file

      @param word     the 32 bit word that corresponds to a sampling interval
      @param UP_nDW   Select *up* or *down* channels
                          UP_nDW  | selection
                          :------:|:--------:
                             0    |   down
                             1    |   up
      @param L1_nL5   Select *L1* or *L5* frequency band
                          L1_nL5  | selection
                          :------:|:--------:
                             0    |   L5
                             1    |   L1
      @return         the selected real signal sample
    */
    static double ExtractBIBASPIRISample(int word, bool UP_nDW, bool L1_nL5);
    /**
      Extract one in-quadrature signal sample (UP/DW, L1/L5, I/Q) from a four-bytes word that corresponds to a sampling instant in a BIBA-SPIR file

      @param word     the 32 bit word that corresponds to a sampling interval
      @param UP_nDW   Select *up* or *down* channels
                          UP_nDW  | selection
                          :------:|:--------:
                             0    |   down
                             1    |   up
      @param L1_nL5   Select *L1* or *L5* frequency band
                          L1_nL5  | selection
                          :------:|:--------:
                             0    |   L5
                             1    |   L1
      @return         the selected real signal sample
    */
    static double ExtractBIBASPIRQSample(int word, bool UP_nDW, bool L1_nL5);

  public:

    /** Read a complex signal from a BIBA-SPIR file.

      @param filename     the name of the file to be read.
      @param UP_nDP       Select an UP or DW signal
                          UP_nDW  | selection
                          :------:|:--------:
                             0    |   down
                             1    |   up
      @param L1_nL5       Select a frequency band
                          L1_nL5  | selection
                          :------:|:--------:
                             0    |   L5
                             1    |   L1
      @param parallel     Set to 'true' to use parallel code execution, to 'false' for sequential code
      @return             Returns a pointer to the complex signal read.
    */
    static ComplexSignal LoadFromBIBASPIRFile(std::string filename, bool UP_nDW, bool L1_nL5, bool parallel=true);
 

    // ----------------------------------------------------
    // Virtual Functions of the parent class 
    ComplexSignal& Load(std::string fname);
    void Save(std::string fname);
    void Print();
    void Plot(const std::string fname="");
    void PlotMagPhase(std::string fname="", unsigned int n0=0, unsigned int n1=80e6) const;
#ifndef TEST_COMPLEXSIGNAL
  private:
#else
  public:
#endif
    // Set things
    void SetLength(unsigned int length);
    // Memory Management
    void AllocateMemory(unsigned int length);
    void FreeMemory();
};


#endif  //COMPLEXSIGNAL_HPP_
