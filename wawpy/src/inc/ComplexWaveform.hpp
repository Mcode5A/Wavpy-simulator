#ifndef COMPLEX_WAVEFORM_HPP_
#define COMPLEX_WAVEFORM_HPP_


/**
    \file   ComplexWaveform.hpp
    \brief  Public header file for the ComplexWaveform class.
*/

#include "Observable.hpp"

#include <gsl/gsl_complex_math.h>

class ComplexWaveform : public Observable
{
  // Member variables
#ifndef TEST_COMPLEXWAVEFORM
  private:
#else
  public:
#endif
    unsigned int mLength;       /**< number of lags of the complex waveform */ 
    double*      mReal=NULL;    /**< pointer to the real part of the waveform */
    double*      mImag=NULL;    /**< pointer to the imaginary part of the waveform */
    double       mTs;           /**< sampling period in [s] */
    int          mFirstLag;     /**< position of the initial lag of the waveform respect to the origin [sample number] */

  // Public member functions
  public:
    // Constructors & destructor
    ComplexWaveform(unsigned int length=0, double Ts=12.5e-9, int first_lag=0); ///
    ComplexWaveform(unsigned int length, const double* real, const double* imag, double Ts=12.5e-9, int initial_lag=0); //
    ComplexWaveform(const ComplexWaveform& cw); ///
    ~ComplexWaveform(); ///

  public:
    // Set things
    void SetTs(double Ts); ///
    void SetFirstLag(double first_lag); ///
    void SetWaveform(unsigned int length, const double* real, const double* imag); ///
    void SetValue(unsigned int pos, double vreal, double vimag); ///
    void SetValue(unsigned int pos, gsl_complex vcmplx); ///
    
    // Get things
    unsigned int GetLength() const; //
    double       GetTs() const; ///
    int          GetFirstLag() const; //
    gsl_complex  GetValue(unsigned int index) const; //
    double       GetReal(unsigned int index) const;  // 
    double       GetImag(unsigned int index) const;  //
    double       GetMag(unsigned int index) const;   //
    double       GetPhase(unsigned int index) const; //
 
    // Print things
    //void Print(unsigned int length=5) const;
    //void PrintLength() const;
    //void PrintTs() const;
    //void PrintInitialLag() const;
    
    // Operators
    bool operator== (const ComplexWaveform& cw) const; ///
    bool operator!= (const ComplexWaveform& cw) const; ///
    ComplexWaveform& operator= (const ComplexWaveform& w); ///
    ComplexWaveform operator+ (const ComplexWaveform& w) const; ///
    ComplexWaveform operator+ (gsl_complex v) const; ///
    ComplexWaveform operator- (const ComplexWaveform& w) const; ///
    ComplexWaveform operator- (gsl_complex v) const; ///
    ComplexWaveform operator* (const ComplexWaveform& w) const; ///
    ComplexWaveform operator* (double a) const; ///
    ComplexWaveform operator/ (const ComplexWaveform& w) const; ///
    ComplexWaveform operator/ (double a) const; ///
    ComplexWaveform Conj(void) const; ///
    friend std::ostream& operator<< (std::ostream &out, const ComplexWaveform& s); ///
    // Virtual Functions of the parent class
    ComplexWaveform& Load(std::string fname); /// TODO
    void Save(std::string fname);             /// TODO
    void Print();                             ///
    void Plot(const std::string fname="");    ///
    void PlotMagPhase(const std::string fname="");
    void PlotMag(const std::string fname="");
    void PlotPhase(const std::string fname="");
    void PlotPow(const std::string fname="");
    void PlotReIm(const std::string fname="");

#ifndef TEST_COMPLEXWAVEFORM
  private:
#else
  public:
#endif
    // Private (auxiliary) member functions
    // manage memory
    void AllocateMemory(unsigned int length); ///
    void FreeMemory(); ///
    
    // Set things
    void SetLength(unsigned int length); ///
};



#endif
