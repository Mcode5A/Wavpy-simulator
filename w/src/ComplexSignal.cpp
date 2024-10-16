#ifndef COMPLEXSIGNAL_CPP_
#define COMPLEXSIGNAL_CPP_
/**
    \file   ComplexSignal.cpp
    \brief  Implementation of the ComplexSignal class
*/

#include "ComplexSignal.hpp"
#include "BIBASPIRfile.hpp"
#include "Plots.hpp"


#include <iostream>
//#include <iomanip>
#include <fstream>
#include <string>
#include <string.h>
//#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex_math.h>


namespace plt = matplotlibcpp;

// PUBLIC MEMBER FUNCTIONS
// Constructors & destructor
//---------------------------------------------------------
ComplexSignal::ComplexSignal(unsigned int length, double Fs, long int initial_sample)
{
  AllocateMemory(length);
  SetLength(length);
  SetFs(Fs);
  SetInitialSample(initial_sample);
}

//---------------------------------------------------------
ComplexSignal::ComplexSignal(unsigned int length, double* real, double* imag, double Fs, long int initial_sample)
{
  AllocateMemory(length);
  SetLength(length);
  SetComplex(length, real, imag);
  SetFs(Fs);
  SetInitialSample(initial_sample);
}

//---------------------------------------------------------
ComplexSignal::ComplexSignal(const ComplexSignal& cs)
{
  mLength = cs.GetLength();
  AllocateMemory(mLength);
  memcpy(mReal, cs.mReal, sizeof(*mReal)*mLength);
  memcpy(mImag, cs.mImag, sizeof(*mImag)*mLength);
  mFs = cs.GetFs();
  mInitialSample = cs.GetInitialSample();
}

//---------------------------------------------------------
ComplexSignal::~ComplexSignal()
{
  FreeMemory();
  SetLength(0);
}


// SET THINGS
//---------------------------------------------------------
void ComplexSignal::SetFs(double Fs)
{
  mFs = Fs;
  return;
}

//---------------------------------------------------------
void ComplexSignal::SetInitialSample(long int initial_sample)
{
  mInitialSample = initial_sample;
  return;
}

//---------------------------------------------------------
void ComplexSignal::SetComplex(unsigned int length, double* real, double* imag)
{
  if (mLength != length)
  {
    std::string s;
    s = "ERROR: assigning complex data of length" + std::to_string(length) + " to a ComplexSignal of length " + std::to_string(mLength) + ".";
    throw(s);
  }
  SetReal(length, real);
  SetImag(length, imag);
  return;
}

//---------------------------------------------------------
void ComplexSignal::SetReal(unsigned int length, double* real)
{
  if (mLength != length)
  {
    std::string s;
    s = "ERROR: assigning real data of length" + std::to_string(length) + " to a ComplexSignal of length " + std::to_string(mLength) + ".";
    throw(s);
  }
  if (NULL == real)
  {
    // throw("ERROR: NULL pointer when assigning real data to ComplexSignal.");
    // putting all zeros when passing a NULL data pointer
    for (unsigned int k=0; k<mLength; k++)
    {
      mReal[k] = 0;
    }
  }
  else
  {
    memcpy(mReal, real, length*sizeof(double));
  }
  return;
}

//---------------------------------------------------------
void ComplexSignal::SetImag(unsigned int length, double* imag)
{
  if (mLength != length)
  {
    std::string s;
    s = "ERROR: assigning imaginary data of length" + std::to_string(length) + " to a ComplexSignal of length " + std::to_string(mLength) + ".";
    throw(s);
  }
  if (NULL == imag)
  {
    //throw("ERROR: NULL pointer when assigning imaginary data to ComplexSignal.");
    // putting all zeros when passing a NULL data pointer
    for (unsigned int k=0; k<mLength; k++)
    {
      mImag[k] = 0;
    }
  }
  else
  {
    memcpy(mImag, imag, length*sizeof(double));
  }
  return;
}

//---------------------------------------------------------
void ComplexSignal::SetSample(unsigned int pos, double vreal, double vimag)
{
  if (pos >= mLength)
  {
    std::string s;
    s = "ERROR: Trying to assign values at position " + std::to_string(pos) + " in a ComplexSignal of length " + std::to_string(mLength) + ".";
    throw(s);
  }
  mReal[pos] = vreal;
  mImag[pos] = vimag;
  return;
}

// GET THINGS
//---------------------------------------------------------
unsigned int ComplexSignal::GetLength() const 
{
  return mLength;
}

//---------------------------------------------------------
double ComplexSignal::GetFs() const 
{
  return mFs;
}

//---------------------------------------------------------
int ComplexSignal::GetInitialSample() const
{
  return mInitialSample;
}

// OPERATORS
//---------------------------------------------------------
bool ComplexSignal::operator== (const ComplexSignal& s) const
{
  if (this == &s)                         return true;  // Self-comparision
  if (mLength != s.mLength)               return false;
  if (mFs != s.mFs)                       return false;
  if (mInitialSample != s.mInitialSample) return false;
  for (unsigned int k=0; k < mLength; k++)
  {
    //std::cout << "k= " << k << " re= " << mReal[k] << " " << s.mReal[k] << " im= " << mImag[k] << " " << s.mImag[k] << std::endl;
    if (mReal[k] != s.mReal[k]) return false;
    if (mImag[k] != s.mImag[k]) return false;
  }
  return true;
}

//---------------------------------------------------------
bool ComplexSignal::operator!= (const ComplexSignal& s) const
{
  return (!(*this==s));
}

//---------------------------------------------------------
ComplexSignal& ComplexSignal::operator= (const ComplexSignal &s)
{
  // Check for self-assignment
  if (this != &s)
  {
    FreeMemory();
    mLength        = s.mLength;
    mFs            = s.mFs;
    mInitialSample = s.mInitialSample;
    AllocateMemory(mLength);
    memcpy(mReal, s.mReal, mLength*sizeof(*mReal));
    memcpy(mImag, s.mImag, mLength*sizeof(*mImag));
  }
  return *this;
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::operator+ (const ComplexSignal& s) const
{
  if (mLength!= s.mLength)
  {
    throw("Cannot add ComplexSignals of different length.");
  }
  if (mFs != s.mFs)
  {
    std::cout << "WARNING: Adding ComplexSignals with different sampling frequencies." << std::endl;
  }
  if (mInitialSample != s.mInitialSample)
  {
    std::cout << "WARNING: Adding ComplexSignals with different initial samples." << std::endl;
  }
  
  ComplexSignal cs(mLength, mFs, mInitialSample);
  for (unsigned int k=0; k<mLength; k++)
  {
    cs.mReal[k] = mReal[k] + s.mReal[k];
    cs.mImag[k] = mImag[k] + s.mImag[k];
  }
  return cs;
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::operator+ (gsl_complex v) const
{
  ComplexSignal cs(mLength, mFs, mInitialSample);
  double re = GSL_REAL(v);
  double im = GSL_IMAG(v);

  for (unsigned int k=0; k<mLength; k++)
  {
    cs.mReal[k] = mReal[k] + re;
    cs.mImag[k] = mImag[k] + im;
  }
  return cs;
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::operator- (const ComplexSignal& s) const
{
  if (mLength!= s.mLength)
  {
    throw("Cannot subtract ComplexSignals of different length.");
  }
  if (mFs != s.mFs)
  {
    std::cout << "WARNING: Subtracting ComplexSignals with different sampling frequencies." << std::endl;
  }
  if (mInitialSample != s.mInitialSample)
  {
    std::cout << "WARNING: Subtracting ComplexSignals with different initial samples." << std::endl;
  }
  
  ComplexSignal cs(mLength, mFs, mInitialSample);
  for (unsigned int k=0; k<mLength; k++)
  {
    cs.mReal[k] = mReal[k] - s.mReal[k];
    cs.mImag[k] = mImag[k] - s.mImag[k];
  }
  return cs;
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::operator- (gsl_complex v) const
{
  double re = GSL_REAL(v);
  double im = GSL_IMAG(v);
  gsl_complex  w = gsl_complex_rect(-re, -im);
  return (this->operator+(w)); 
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::operator* (const ComplexSignal& s) const
{
  if (mLength!= s.mLength)
  {
    throw("Cannot multiply ComplexSignals of different length.");
  }
  if (mFs != s.mFs)
  {
    std::cout << "WARNING: Multiplying ComplexSignals with different sampling frequencies." << std::endl;
  }
  if (mInitialSample != s.mInitialSample)
  {
    std::cout << "WARNING: Multiplying ComplexSignals with different initial samples." << std::endl;
  }
  
  ComplexSignal cs(mLength, mFs, mInitialSample);
  for (unsigned int k=0; k<mLength; k++)
  {
    cs.mReal[k] = mReal[k]*s.mReal[k] - mImag[k]*s.mImag[k];
    cs.mImag[k] = mImag[k]*s.mReal[k] + mReal[k]*s.mImag[k];
  }

  // Manage the name, timetag and origin
  cs.SetName(this->GetName() + " * " + s.GetName());
  cs.SetTimeTag(this->GetTimeTag() + " & " + s.GetTimeTag());
  cs.SetOrigin(this->GetOrigin() + " &\n" + s.GetOrigin());

  return cs;
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::operator* (double a) const
{
  ComplexSignal cs(mLength, mFs, mInitialSample);
  for (unsigned int k=0; k<cs.mLength; k++)
  {
    cs.mReal[k] = mReal[k]*a;
    cs.mImag[k] = mImag[k]*a;
  }
  return cs;
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::operator/ (const ComplexSignal& s) const 
{
  if (mLength!= s.mLength)
  {
    throw("Cannot divide ComplexSignals of different length.");
  }
  if (mFs != s.mFs)
  {
    std::cout << "WARNING: Dividing ComplexSignals with different sampling frequencies." << std::endl;
  }
  if (mInitialSample != s.mInitialSample)
  {
    std::cout << "WARNING: Dividing ComplexSignals with different initial samples." << std::endl;
  }
  ComplexSignal sconj = s.Conj();
  ComplexSignal cs = this->operator*(sconj);
  for (unsigned int k=0; k<mLength; k++)
  {
    double div = s.mReal[k]*s.mReal[k] + s.mImag[k]*s.mImag[k];
    if (div!=0)
    {
      cs.mReal[k] = cs.mReal[k] / div;
      cs.mImag[k] = cs.mImag[k] / div;
    }
    else
    {
      if (cs.mReal[k] > 0)  cs.mReal[k] =  INFINITY;
      if (cs.mReal[k] == 0) cs.mReal[k] =  SNAN;
      if (cs.mReal[k] < 0)  cs.mReal[k] = -INFINITY;
      if (cs.mImag[k] > 0)  cs.mImag[k] =  INFINITY;
      if (cs.mImag[k] == 0) cs.mImag[k] =  SNAN;
      if (cs.mImag[k] < 0)  cs.mImag[k] = -INFINITY;
      std::cout << "WARNING: Division by zero." <<  std::endl;
    }
  }
  
  //ComplexSignal* cs_ret = new ComplexSignal(0);
  //*cs_ret = cs; 
  return cs;
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::operator/ (double a) const
{
  double b = 1/a;
  return (this->operator*(b));
}


//---------------------------------------------------------
ComplexSignal ComplexSignal::Conj(void) const
{
  ComplexSignal cs=*this;
  //memcpy(cs->mReal, mReal, mLength*sizeof(*mReal));
  for (unsigned int k=0; k<mLength; k++)
  {
    cs.mImag[k] = -mImag[k];
  }
  return cs;
}

//---------------------------------------------------------
ComplexSignal& ComplexSignal::Crop(unsigned int n0, unsigned int n1) const
{
  // Check input paramters
  if (n1 < n0)
  {
    throw("n1 cannot be smaller than n0 when cropping a ComplexSignal.");
  }
  if (n0 >= GetLength())
  {
    throw("Cannot crop signal for n0.");
  }
  if (n1 >= GetLength())
  {
    throw("Cannot crop signal with current value of n1.");
  }

  // Length of the new signal
  unsigned int L = n1 - n0 + 1;
  double *re, *im;
  re = mReal+n0;
  im = mImag+n0;
  ComplexSignal* cs_ret = new ComplexSignal(L, re, im, GetFs(), GetInitialSample()+n0);

  return *cs_ret;
}

//---------------------------------------------------------
gsl_complex ComplexSignal::Mean(void) const
{
  gsl_complex m;
  double re=0;
  double im=0;
  for (unsigned int k=0; k<mLength; k++)
  {
    re = re + mReal[k];
    im = im + mImag[k];
  }
  re = re/mLength;
  im = im/mLength;
  m = gsl_complex_rect(re, im);
  return m;
}

//---------------------------------------------------------
double ComplexSignal::Var(void) const
{
  gsl_complex m = Mean();
  double m_re = GSL_REAL(m);
  double m_im = GSL_IMAG(m);

  double v = 0;
  for (unsigned int k=0; k<mLength; k++)
  {
    double c, d;
    c = mReal[k] - m_re;
    d = mImag[k] - m_im;
    v = v + c*c + d*d;
  }
  v = v / (mLength-1);

  return v;
}

//---------------------------------------------------------
double ComplexSignal::Std(void) const
{
  return sqrt(Var());
}

//---------------------------------------------------------
std::ostream& operator<< (std::ostream &out, const ComplexSignal& s)
{
  out << std::endl;
  out << "ComplexSignal:" << std::endl;
  out << "  Length:     " << s.mLength        << std::endl;
  out << "  Fs:         " << s.mFs            << std::endl;
  out << "  InitSample: " << s.mInitialSample << std::endl;
  unsigned int kmax = 10;
  std::string endstr = "...";
  if (s.mLength < kmax)
  {
    kmax = s.mLength;
    endstr = "";
  }
  for (unsigned int k=0; k<kmax; k++)
  {
    out << "cs[" << k << "]=" << s.mReal[k];
    if (s.mImag[k]>0)
    {
      out << "+j" << s.mImag[k] << "  ";
    }
    else
    {
      out << "-j" << -s.mImag[k] << "  ";
    }
  }
  out << endstr << std::endl;
  return out;
}


// SPECIAL SIGNALS
//---------------------------------------------------------
ComplexSignal ComplexSignal::Phasor(double Freq, double Fs, unsigned int length, double phi_0, bool parallel)
{
  // the initial phase phi_0 is in radians !
  long int initial_sample = 0;
  ComplexSignal P(length, Fs, initial_sample);
  double Fdig;

  Fdig = Freq/Fs;

  if (true==parallel)
  {
    throw("Parallelization not implemented in ComplexSignal::Phasor.");
  }
  else
  {
    for (unsigned int k=0; k<length; k++)
    {
      double arg = 2*M_PI*Fdig*k + phi_0;
      P.mReal[k] = cos(arg);
      P.mImag[k] = sin(arg);
    }
  }
  return P;
}

//----------------------------------------------------------
ComplexSignal ComplexSignal::AWGN(double sigma, unsigned int length, double Fs, int seed)
{
  static gsl_rng* r = NULL;

  if (NULL==r)
  { 
    //std::cout << "Setting up rng." << std::endl;
    r = gsl_rng_alloc(gsl_rng_ranlxs0);  
    if (NULL==r)
    {
      throw("ERROR! Could not obtain memory when creating random number generator when creating AWGN.");
    }
  }
  if (seed!=-1)
  {
    //std::cout << "Setting seed: " << seed << std::endl;
    gsl_rng_set(r, (unsigned long int) seed);
  }
  
  double* IQ;
  IQ = (double*)malloc(sizeof(double)*length*2);
  if (NULL==IQ)
  {
    throw("ERROR! Could not obtain memory when generating AWGN.");
  }

  sigma = abs(sigma); // Make sigma positive for the case it is negative
  for (unsigned int k=0;k<length;k++)
  {
    IQ[k]        = gsl_ran_gaussian(r, sigma);
    IQ[k+length] = gsl_ran_gaussian(r, sigma);
  }
  ComplexSignal s(length, &(IQ[0]), &(IQ[length]), Fs);
  free(IQ);
  gsl_rng_free(r);

  return s;
}


//---------------------------------------------------------
ComplexSignal ComplexSignal::GPS_CA(unsigned char PRN, double Fs, unsigned int length)
{
  /** C/A PRN signal as defined in *INTERFACE SPECIFICATION IS-GPS-200, REVISION E, Navstar GPS Space Segment/Navigation User Interfaces, 8 June, 2010*.
  */
  int* ca_code;  
  double Fchip = 1.023e6;
  double Isample, Qsample;
  double chip_step;
  int chip;

  ca_code = (int*)calloc(1024, sizeof(int));
  if (NULL==ca_code)
  {
    std::cerr << "Could not allocate memory for CA code generation." << std::endl; // TODO
    return (EXIT_FAILURE);
  }
  CACode(ca_code, PRN);
  ComplexSignal S = ComplexSignal(length, Fs, 0);

  chip_step = Fchip / Fs;

  Isample = 0;
  for (unsigned int k=0;k<length;k++)
  {
    chip = ((int)floor((chip_step*k))) %1023;
    //std::cout << "k: " << k << "  chip #: " << chip << "  ca_code[" << chip << "]: " << ca_code[chip] << std::endl;
    if (ca_code[chip] == 1)
    {
      Qsample = 1;
    }
    else
    {
      Qsample = -1;
    }
    S.SetSample(k, Isample, Qsample);
  }
  
  free(ca_code);
  return S;
}


//==========================================================
// SPECIAL OPERATIONS
//==========================================================
//-----------------------------------------------------
// Functions related to CA code generation

//---------------------------------------------------------
int ComplexSignal::CACode(int *ca_code, unsigned int PRN)
{
  unsigned int G1, G2;
  unsigned int G1bit, G2bit;  // bits to update registers
  unsigned int G1out, G2out;  // output bits
 
  const unsigned int G1MASK = 0x00000408;
  const unsigned int G2MASK = 0x0000074C;

  const unsigned int G1init = 0x000007FE;
  const unsigned int G2init = 0x000007FE;
  
  const unsigned int MASK11BITS = 0x7FF;
  const unsigned int MASK10BITS = 0x7FE;

  const unsigned int G2PRNMASK[38] = { 0x000,  // PRN  0
                                       0x044,  // PRN  1
                                       0x088,  // PRN  2
                                       0x110,  // PRN  3
                                       0x220,  // PRN  4
                                       0x202,  // PRN  5
                                       0x404,  // PRN  6
                                       0x102,  // PRN  7
                                       0x204,  // PRN  8
                                       0x408,  // PRN  9
                                       0x00C,  // PRN 10
                                       0x018,  // PRN 11
                                       0x060,  // PRN 12
                                       0x0C0,  // PRN 13
                                       0x180,  // PRN 14
                                       0x300,  // PRN 15
                                       0x600,  // PRN 16
                                       0x012,  // PRN 17
                                       0x024,  // PRN 18
                                       0x048,  // PRN 19
                                       0x090,  // PRN 20
                                       0x120,  // PRN 21
                                       0x240,  // PRN 22
                                       0x00A,  // PRN 23
                                       0x050,  // PRN 24
                                       0x0A0,  // PRN 25
                                       0x140,  // PRN 26
                                       0x280,  // PRN 27
                                       0x500,  // PRN 28
                                       0x042,  // PRN 29
                                       0x084,  // PRN 30
                                       0x108,  // PRN 31
                                       0x210,  // PRN 32
                                       0x420,  // PRN 33
                                       0x410,  // PRN 34 
                                       0x082,  // PRN 35
                                       0x104,  // PRN 36
                                       0x410}; // PRN 37
  

  if (PRN<1 || PRN>37)  // TODO: handle error properly
  {
    std::cerr << "ERROR: PRN must be in the range [1,37]." << std::endl;
    return 0;
    exit(EXIT_FAILURE);
  }
  
  if (NULL==ca_code)  // TODO: handle error properly
  {
    std::cerr << "ERROR: The pointer to store the CA code must not be NULL." << std::endl;
    return 0;
    exit(EXIT_FAILURE);
  }

  // Init registers
  G1 = G1init;
  G2 = G2init;

  for (unsigned int k=0;k<1023;k++)
  {
    //std::cout << "  G1[" << std::dec << k << "]: " << std::hex << (G1 & MASK10BITS);
    //std::cout << "  G2[" << std::dec << k << "]: " << std::hex << (G2 & MASK10BITS);
    //std::cout << std::dec;
    
    // output bits
    G1out = (G1 >> 10) & 0x01;
    G2out = (__builtin_popcount(G2 & G2PRNMASK[PRN]))%2;

    //std::cout << "  G1out[" << std::dec << k << "]: " << G1out;
    //std::cout << "  G2out[" << std::dec << k << "]: " << G2out;
    //std::cout << "  OUT["   << std::dec << k << "]: " << (G1out ^ G2out);

    // Compute output (chipping code)
    if ((G1out ^ G2out)==1)
    {
      ca_code[k] = 1;
    }
    else
    {
      ca_code[k] = 0;
    }
    // Feedback bits of the polynomial 
    G1bit  = (__builtin_popcount(G1 & G1MASK))%2;
    G2bit  = (__builtin_popcount(G2 & G2MASK))%2;
    //std::cout << "  G1bit[" << std::dec << k << "]: " << G1bit;
    //std::cout << "  G2bit[" << std::dec << k << "]: " << G2bit;

    // Compute contents of register for next clock cycle
    G1 = (((G1 & MASK10BITS) | G1bit) << 1) & MASK11BITS;
    G2 = (((G2 & MASK10BITS) | G2bit) << 1) & MASK11BITS;
    //std::cout << std::endl;

  } 
  return 1; 
}




// Functions related to correlations
//---------------------------------------------------------
bool ComplexSignal::CorrelationLimitsCorrect(unsigned int length_x, unsigned int length_y, unsigned int N, unsigned int n0, int min_m, int max_m)
{
  int n_start, n_stop;

  // Check the number of samples
  if (N==0)
  {
    std::cerr << "WARNING: At least one sample needs to be correlated, but N=0." << std::endl;
    return false;
  }

  // Check limits for x
  n_start = (int)n0;    
  n_stop  = (int)n0+(int)N-1;
  // start limit
  if (n_start<0 || n_start>=(int)length_x)
  {
    std::cerr << "WARNING: Correlation interval out of range for signal x. n_start = " << n_start << " [0," << length_x << ")." << std::endl;
    return false;
  }
  // stop limit
  if (n_stop<0 || n_stop>=(int)length_x)
  {
    std::cerr << "WARNING: Correlation interval out of range for signal x. n_stop = " << n_stop << " [0," << length_x << ")." << std::endl;
    return false;
  }

  // Check limits for y when m=m_min 
  n_start = (int)n0-min_m;    
  n_stop  = (int)n0+N-1-min_m;
  // start limit
  if (n_start<0 || n_start>=(int)length_y)
  {
    std::cerr << "WARNING: Correlation interval out of range for signal y when m=" << min_m << ". n_start = " << (int) n_start << " [0," << length_y << ")." << std::endl;
    return false;
  }
  // stop limit
  if (n_stop<0 || n_stop>=(int)length_y)
  {
    std::cerr << "WARNING: Correlation interval out of range for signal y when m=" << min_m << ". n_stop = " << (int) n_stop << " [0," << length_y << ")." << std::endl;
    return false;
  }

  // Check limits for y when m=m_max 
  n_start = (int)n0-max_m;    
  n_stop  = (int)n0+N-1-max_m;
  // start limit
  if (n_start<0 || n_start>=(int)length_y)
  {
    std::cerr << "WARNING: Correlation interval out of range for signal y when m=" << max_m << ". n_start = " << (int) n_start << " [0," << length_y << ")." << std::endl;
    return false;
  }
  // stop limit
  if (n_stop<0 || n_stop>=(int)length_y)
  {
    std::cerr << "WARNING: Correlation interval out of range for signal y when m=" << max_m << ". n_stop = " << (int) n_stop << " [0," << length_y << ")." << std::endl;
    return false;
  }
  return true;
}




//---------------------------------------------------------
double ComplexSignal::CorrelationRealSingleLag(const double* x, const double* y, unsigned int N, unsigned int n0, int m, bool parallel)
{
  /**
    The real cross-correlation is computed according to:
    \f[
      R[m] = \frac{1}{N} \sum_{n=n_0}^{n_0 + N -1} x[n] y[n-m]
    \f]
  */
  double Rxy = 0;

  if (parallel==1) // parallel code execution
  { 
    /*
    std::cout << "    double* CorrelationRealSingleLagTimeDomain(...) is parallel" << std::endl;
    //double start_time, end_time;
    //start_time = omp_get_wtime();

    omp_set_num_threads(12);
    #pragma omp parallel
    {
      double Rxy_partial = 0;
      int id, nthrds;
      id = omp_get_thread_num();
      nthrds = omp_get_num_threads();
      
      for (unsigned int n=n0+id; n<(n0+N); n=n+nthrds)
      {
        Rxy_partial = Rxy_partial + x[n]*y[n-m];
      } 
      #pragma omp critical
      {
        Rxy = Rxy+Rxy_partial;
      }
    }
    //end_time = omp_get_wtime();
    //std::cout << "Single lag time: " << end_time - start_time << std::endl;
    */
    std::cout << "WARNING: Parallel correlation not implemented. Using sequential implementation." << std::endl;
  }
  { // Sequential implementation of the correlation.
    //std::cout << "    double* CorrelationRealSingleLagTimeDomain(...) is sequential" << std::endl;
    for (unsigned int n=n0; n<(n0+N); n++)
    {
      //std::cout << "  n=" << n << std::endl;
      Rxy = Rxy + x[n]*y[n-m];
    } 
  
    //std::cout << "Rxy = " << Rxy << std::endl;
  }
  Rxy = Rxy / N;
  return Rxy;
}

//----------------------------------------------------------
gsl_complex ComplexSignal::CorrelationComplexSingleLag(const ComplexSignal& x, const ComplexSignal& y, unsigned int N, unsigned int n0, int m, bool parallel)
{
  /**
    The complex cross-correlation is computed according to:
    \f[
      R[m] = \frac{1}{N} \sum_{n=n_0}^{n_0 + N -1} x[n] y^*[n-m]
    \f]
  */
  gsl_complex Rxy;
  double II=0;
  double QI=0;
  double QQ=0;
  double IQ=0;

  //std::cout << "N=" << N << "   n0=" << n0 << "   m=" << m << std::endl;    
  II = CorrelationRealSingleLag(x.mReal, y.mReal, N, n0, m, parallel);
  QI = CorrelationRealSingleLag(x.mImag, y.mReal, N, n0, m, parallel);
  QQ = CorrelationRealSingleLag(x.mImag, y.mImag, N, n0, m, parallel);
  IQ = CorrelationRealSingleLag(x.mReal, y.mImag, N, n0, m, parallel);

  Rxy = gsl_complex_rect(II+QQ, QI-IQ);

  return Rxy;  
}

//---------------------------------------------------------
ComplexWaveform ComplexSignal::Correlation(const ComplexSignal& x, const ComplexSignal& y, unsigned int N, unsigned int n0, int min_m, unsigned int step_m, int max_m, bool parallel)
{
  /** The cross-correlation of the signals is calculated according to
    \f[
      R[m] = \frac{1}{N} \sum_{n=n_0}^{n_0 + N -1} x[n] y^*[n-m]
    \f]
  */
  unsigned int length_R; 
  int initial_sample_R;
  double Fs_R;
  if (!CorrelationLimitsCorrect(x.GetLength(), y.GetLength(), N, n0, min_m, max_m))
  {
    throw("ERROR: Incorrect correlation limits.");
  }

  // Check that the sampling frequencies of both signals are equal
  if (x.GetFs()!=y.GetFs())
  {
    std::cerr << "WARNING: Sampling frequencies of the signals are different. Fsx = " << x.GetFs() << " Hz. Fsy = " << y.GetFs() << " Hz." << std::endl;
  }
  
  length_R = floor((max_m-min_m)/step_m) + 1;
  initial_sample_R = (min_m + x.mInitialSample - y.mInitialSample)/(int)step_m; 
  Fs_R = x.GetFs()/step_m;
  ComplexWaveform Rxy(length_R, 1/Fs_R, initial_sample_R);

  if (parallel==true)  // parallel code
  {
    std::cout << "    Parallel computation of cross-correlation not implemented yet. Using sequential computation." << std::endl;
    /*
    std::cout << "    ComplexWaveform* CorrelationInTimeDomain(...) is parallel" << std::endl;
    // I parallelize for each lag
    double start_time, end_time;
    start_time = omp_get_wtime();
    omp_set_num_threads(12);
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      int nthrds = omp_get_num_threads();
  
      for (int m=min_m+(id*step_m), k=id;m<=max_m; m=m+(step_m*nthrds), k=k+nthrds)
      {
        gsl_complex r = CorrelationSingleLagTimeDomain(x, y, N, n0, m);
        Rxy->SetLag(k, GSL_REAL(r), GSL_IMAG(r));
      }
    }
    end_time = omp_get_wtime();
    std::cout << "Time elapsed to compute one waveform: " << end_time - start_time << std::endl;
    */
  }
  //else
  {
    std::cout << "    ComplexWaveform* CorrelationInTimeDomain(...) is sequential" << std::endl;
    for (int m=min_m, k=0;m<=max_m; m=m+step_m, k++)
    {
      //std::cout << "m=" << m << std::endl;
      gsl_complex r = CorrelationComplexSingleLag(x, y, N, n0, m);
      Rxy.SetValue(k, GSL_REAL(r), GSL_IMAG(r));
    }
  }

  // Manage the name, timetag and origin
  Rxy.SetName(x.GetName() + " (XCOR) " + y.GetName());
  Rxy.SetTimeTag(x.GetTimeTag() + " & " + y.GetTimeTag());
  Rxy.SetOrigin(x.GetOrigin() + " &\n " + y.GetOrigin());
  return Rxy;
}


//---------------------------------------------------------
// Load BIBASPIR files
//---------------------------------------------------------
bool ComplexSignal::IsBIBASPIRFlag(const int* bibaspirdata)
{
  unsigned char* buf;
  buf = (unsigned char*) bibaspirdata;
  unsigned char flag;
  flag = buf[BIBASPIR_FLAG_POS];

  if (BIBASPIR_FLAG==flag)
  {
    return true;
  }
  else
  {
    return false;
  }
}

//---------------------------------------------------------
bool ComplexSignal::CheckBIBASPIRBlocks(int* bibaspirdata)
{
  int word, blockID, anciID;

  for (int i=0; i<BIBASPIR_NBLOCKS; i++)
  {
    word = bibaspirdata[i*BIBASPIR_BLOCK_SIZE_INTS];
    blockID = word & BIBASPIR_BLOCK_ID_MASK;
    anciID  = (word>>BIBASPIR_ANCI_ID_BIT) & BIBASPIR_ANCI_ID_MASK;
    if ((anciID!=0) || (blockID!=i))
    {
      std::cout << "WARNING: Wrong block ID at block number: " << i << std::endl;
      return false;
    }
  }
  return true;
}

//---------------------------------------------------------
double ComplexSignal::ExtractBIBASPIRSample(int word, bool UP_nDW, bool L1_nL5, bool I_nQ)
{
  int  maskUPDW;
  int  maskL1L5;
  int  maskIQ;
  int  mask;
  int  shift;
  char aux, sample;

  switch ((int)UP_nDW)
  {
    case true:
      maskUPDW = BIBASPIR_UP_MASK; 
      shift    = 0;
      break;
    case false:
      maskUPDW = ~BIBASPIR_UP_MASK;
      shift = 16;
      break;
  }

  switch ((int)L1_nL5)
  {
    case true:
      maskL1L5 = BIBASPIR_L1_MASK;
      shift    = shift + 0;
      break;
    case false:
      maskL1L5 = ~BIBASPIR_L1_MASK;
      shift    = shift + 8;
      break; 
  }
 
  switch ((int)I_nQ)
  {
    case false:  // The I/Q (Q/I) definition in MAX2112 may be swapped
                 // For this reason I switch the case statements (true/false)
                 // to swap the I and Q components so that the phasors 
                 // rotate in the correct direction
      maskIQ = BIBASPIR_I_MASK;
      shift  = shift + 4;
      break;
    case true:
      maskIQ = ~BIBASPIR_I_MASK;
      shift  = shift + 0;
      break;
  }
  
  mask   = maskUPDW & maskL1L5 & maskIQ;
  aux   = (char)((word & mask)>>shift);  // nibble shifted maximum to the right
  sample = aux<<4;                       // shift four bits to the left to get correct sign 
                      
  return double(sample); 
}

//---------------------------------------------------------
double ComplexSignal::ExtractBIBASPIRISample(int word, bool UP_nDW, bool L1_nL5)
{
  return ExtractBIBASPIRSample(word, UP_nDW, L1_nL5, 1);
}

//---------------------------------------------------------
double ComplexSignal::ExtractBIBASPIRQSample(int word, bool UP_nDW, bool L1_nL5)
{
  return ExtractBIBASPIRSample(word, UP_nDW, L1_nL5, 0);
}

//---------------------------------------------------------
ComplexSignal ComplexSignal::LoadFromBIBASPIRFile(std::string filename, bool UP_nDW, bool L1_nL5, bool parallel)
{
  std::string s;
  int* bibaspirdata=NULL;
  //// ComplexSignal* cs=NULL;

  bibaspirdata = (int *) malloc(BIBASPIR_FILESIZE);  // check allocation went well
  if (NULL==bibaspirdata)
  {
    throw("Error: could not allocate memory to load BIBASPIR file.");
  }

  // Open file and read data
  // Check architecture
  if(sizeof(int) != 4)
  {
    free(bibaspirdata);
    bibaspirdata=NULL;
    throw("ERROR! Size of int must be equal to 4 bytes. Computer architecture not supported.");
  }

  // Check file length
  std::ifstream bibaspirfile;
  bibaspirfile.open(filename.c_str(), std::ios_base::in);
  if (!bibaspirfile)
  {
    free(bibaspirdata);
    bibaspirdata=NULL;
    throw("ERROR! Could not open file.");
  }
  ///
  bibaspirfile.seekg(0, std::ios_base::end);
  std::size_t size = bibaspirfile.tellg();
  if(size != BIBASPIR_FILESIZE)
  {
    ///
  	bibaspirfile.close();
    free(bibaspirdata);
    bibaspirdata=NULL;
    s = "ERROR! File <" + filename + "> is not a complete BIBASPIR file.";
    throw(s);
  }
  
  // read data from file
  bibaspirfile.seekg(0, std::ios_base::beg);
  bibaspirfile.read((char*) bibaspirdata, size);  /* @todo check reading went well */
  bibaspirfile.close();
  if (!bibaspirfile)
  {
    free(bibaspirdata);
    bibaspirdata=NULL;
    s = s + "Error: could only read " + std::to_string(bibaspirfile.gcount()) + " bytes from BIBASPIRFILE.";
    throw(s);
  }
  
  // Check that there is the BIBASPIR flag on the file.
  if (!IsBIBASPIRFlag(bibaspirdata))
  {
    free(bibaspirdata);
    bibaspirdata=NULL;
    throw("Error: This is not a BIBASPIR file, it might be a SPIR file.");
  }
 
  // Check if the blocks are consistently numbered
  if (!CheckBIBASPIRBlocks(bibaspirdata))
  {
    free(bibaspirdata);
    bibaspirdata=NULL;
    throw("Error: Data read from BIBASPIR file: " + filename + " is not consistent.");
  }
 
  // Now extract the samples and fill up the ComplexSignal object
  ComplexSignal cs(BIBASPIR_NSAMPLES);

  if (parallel==true) // Use parallel code
  {
    std::cout << "Parallel BIBASPIR file reading not implemented yet. Using sequential implementation." << std::endl;
  }
    /*
    // I parallelize the outer loop
    double start_time, end_time;
    start_time = omp_get_wtime();
  
    omp_set_num_threads(12);
    #pragma omp parallel 
    {
      int id, nthrds;
  
      id = omp_get_thread_num();
      nthrds = omp_get_num_threads();
  
      for (unsigned int kblock=id;kblock<BIBASPIR_NBLOCKS;kblock=kblock+nthrds)
      {
        unsigned int firstinblock = ((kblock*BIBASPIR_BLOCK_SIZE)/4+1);   // +1 to skip the ID at the begining of the block
        unsigned int lastinblock  = ((kblock+1)*BIBASPIR_BLOCK_SIZE)/4-1; 
        unsigned int l = kblock*(BIBASPIR_BLOCK_SIZE/4-1); 
         
        for (unsigned int k=firstinblock;k<=lastinblock;k++)
        {
          int word = bibaspirdata[k];
          double Isample = ExtractBIBASPIRISample(word,UP_nDW,L1_nL5);
          double Qsample = ExtractBIBASPIRQSample(word,UP_nDW,L1_nL5);
          cs->SetSample(l,Isample, Qsample);
          l++;
        }
      }
    }
    
    end_time = omp_get_wtime();
    std::cout << "Time elapsed when reading block of data: " << end_time - start_time << std::endl;
  }
  else // Use sequential code
*/
  {
    std::cout << "Using sequential code to read BIBASPIR file." << std::endl;
    unsigned int l=0;
    for (unsigned int kblock=0;kblock<BIBASPIR_NBLOCKS;kblock++)
    {
      unsigned int firstinblock = ((kblock*BIBASPIR_BLOCK_SIZE)/4+1);   // +1 to skip the ID at the begining of the block
      unsigned int lastinblock  = ((kblock+1)*BIBASPIR_BLOCK_SIZE)/4-1; 
     
      for (unsigned int k=firstinblock;k<=lastinblock;k++)
      {
        int word = bibaspirdata[k];
        double Isample = ExtractBIBASPIRISample(word,UP_nDW,L1_nL5);
        double Qsample = ExtractBIBASPIRQSample(word,UP_nDW,L1_nL5);
        cs.SetSample(l,Isample, Qsample); // maybe this could be modified to make it faster (SetComplex, after the loop)
        l++;
      }
    }
  }

  // Add labels to the complex signal
  std::string channel;
  std::string band;
  if (UP_nDW)
  {
    channel = "UP";
  }
  else
  {
    channel = "DW";
  }
  if (L1_nL5)
  {
    band = "L1";
  }
  else
  {
    band = "L5";
  }

  cs.SetName("BIBASPIR-" + band + "-" + channel);
  cs.SetOrigin(filename);
  free(bibaspirdata);
  return cs;
}


// VIRTUAL FUNCTIONS OF THE PARENT CLASS
//---------------------------------------------------------
ComplexSignal& ComplexSignal::Load(std::string fname)
{
  std::cout << "ComplexSignal::Load() TBD." << std::endl;
  ComplexSignal* cs = new ComplexSignal();
  return *cs;
}

//---------------------------------------------------------
void ComplexSignal::Save(std::string fname)
{
  std::cout << "ComplexSignal::Save() TBD." << std::endl;
  return;
}

//---------------------------------------------------------
void ComplexSignal::Print()
{
  std::cout << *this;
  return;
}


//---------------------------------------------------------
void ComplexSignal::Plot(std::string fname)
{
  unsigned int length = mLength;

  if (length>PLT_MAX_N)
  {
    std::cout << "Very long signal. Plotting first " << PLT_MAX_N << " samples only." << std::endl;
    length = PLT_MAX_N;
  }

  // Horizontal axis data
  std::vector<double> t(length);
  double Ts = 1/mFs;
  for (unsigned int k = 0; k<length; k++)
  {
    t.at(k) = (double) (k+mInitialSample) * Ts;
  }
  
  std::vector<std::vector<double>> tt(2);
  tt.at(0) = t;
  tt.at(1) = t;
  
  // Vertical axis data
  std::vector<double> re(mReal, mReal + length);
  std::vector<double> im(mImag, mImag + length);

  std::vector<std::vector<double>> datay(2);
  datay.at(0) = re;
  datay.at(1) = im;

  // Labels
  std::vector<std::string> xlabels(2);
  std::vector<std::string> ylabels(2);
  std::vector<std::string> titles(2);

  xlabels.at(0) = "";
  xlabels.at(1) = "Time [s]";
  
  ylabels.at(0) = "Real";
  ylabels.at(1) = "Imag";

  titles.at(0) = GetName() + " - " + GetTimeTag() + "\n" + GetOrigin();
  titles.at(1) = "";

  Plots::PlotMultiple (fname, tt, datay, xlabels, ylabels, titles);

  return;
}

//---------------------------------------------------------
void ComplexSignal::PlotMagPhase(std::string fname, unsigned int n0, unsigned int n1) const
{
  if (n1>=mLength)
  {
    n1 = mLength-1;
  }
  unsigned int length = n1-n0 + 1;
  
  if (length<2)
  {
    throw("ERROR: Cannot print less that two signal samples.");
  }
  if (n0>n1)
  {
    throw("ERROR: Cannot n0 must be smaller or equal n1.");
  }

  // Horizontal axis data
  double Ts = 1/mFs;
  std::vector<double> tau(length);
  for (unsigned int k=0;k<length; k++)
  {
    tau.at(k) = (double) ((int)k + mInitialSample) * Ts;
  }
  std::vector<std::vector<double>> tautau(2);
  tautau.at(0) = tau;
  tautau.at(1) = tau;

  // Vertical axis data
  std::vector<double> mag(length);
  std::vector<double> phase(length);
  for (unsigned int k=0; k<length; k++)
  {
    mag.at(k) = sqrt(mReal[k]*mReal[k] + mImag[k]*mImag[k]);
    phase.at(k) = atan2(mImag[k], mReal[k]) * 180/M_PI;
  }
  std::vector<std::vector<double>>data(2);
  data.at(0) = mag;
  data.at(1) = phase;

  // Labels
  std::vector<std::string> xlabels(2);
  std::vector<std::string> ylabels(2);
  std::vector<std::string> titles(2);

  xlabels.at(0) = "Tau [s]";
  xlabels.at(1) = "Tau [s]";

  ylabels.at(0) = "Magnitude"; 
  ylabels.at(1) = "Phase [deg]";

  titles.at(0) = GetName() + " - " + GetTimeTag() + "\n" + GetOrigin();
  titles.at(1) = "";
  
  Plots::PlotMultiple(fname, tautau, data, xlabels, ylabels, titles);

  return;
}  

// PRIVATE FUNCTIONS
// Set Things
//---------------------------------------------------------
void ComplexSignal::SetLength(unsigned int length)
{
  mLength = length;
}

// MEMORY MANAGEMENT
//---------------------------------------------------------
void ComplexSignal::AllocateMemory(unsigned int length)
{
  if (NULL!=mReal)
  {
    throw("ERROR: Trying to allocate new memory to already allocated real data in ComplexSignal.");
  }
  if (NULL!=mImag)
  {
    throw("ERROR: Trying to allocate new memory to already allocated imaginary data in ComplexSignal.");
  }
  mReal = (double*) calloc(length, sizeof(double));
  if (NULL==mReal)
  {
    throw("ERROR: Could not allocate memory for real data in ComplexSignal.");
  }

  mImag = (double*) calloc(length, sizeof(double));
  if (NULL==mImag)
  {
    throw("ERROR: Could not allocate memory for imaginary data in ComplexSignal.");
  }
  return;
}

//---------------------------------------------------------
void ComplexSignal::FreeMemory()
{
  free(mReal);
  free(mImag);
  mReal = NULL;
  mImag = NULL;
}



#endif  //COMPLEXSIGNAL_CPP_
