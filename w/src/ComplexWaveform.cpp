#ifndef COMPLEX_WAVEFORM_CPP_
#define COMPLEX_WAVEFORM_CPP_

/**
    \file   ComplexWaveform.cpp
    \brief  Contains the definition for the class ComplexWaveform
*/

#include "ComplexWaveform.hpp"
#include "Plots.hpp"

#include <string.h>
#include <gsl/gsl_math.h>
#include <iostream>

//---------------------------------------------------------
ComplexWaveform::ComplexWaveform(unsigned int length, double Ts, int first_lag)
{
  AllocateMemory(length);
  SetLength(length);
  SetTs(Ts);
  SetFirstLag(first_lag);
}

//---------------------------------------------------------
ComplexWaveform::ComplexWaveform(unsigned int length, const double* real, const double* imag, double Ts, int first_lag)
{
  AllocateMemory(length);
  SetLength(length);
  SetWaveform(length, real, imag);
  SetTs(Ts);
  SetFirstLag(first_lag);
}

//---------------------------------------------------------
ComplexWaveform::ComplexWaveform(const ComplexWaveform& cw)
{
  AllocateMemory(cw.mLength);
  SetLength(cw.mLength);
  SetWaveform(cw.mLength, cw.mReal, cw.mImag);
  SetTs(cw.mTs);
  SetFirstLag(cw.mFirstLag);
}

//---------------------------------------------------------
ComplexWaveform::~ComplexWaveform()
{
  FreeMemory();
  SetLength(0);
}

//---------------------------------------------------------
// SET THINGS
//---------------------------------------------------------
void ComplexWaveform::SetWaveform(unsigned int length, const double* real, const double* imag)
{
  if (mLength != length)
  {
    throw("ERROR: assigning data of length " + std::to_string(length) + " to a complex waveform of length " + std::to_string(mLength));
  }
  memcpy(mReal, real, length*sizeof(double));
  memcpy(mImag, imag, length*sizeof(double));
  return;
}

//---------------------------------------------------------
void ComplexWaveform::SetValue(unsigned int pos, double vreal, double vimag)
{
  if (pos >= mLength)
  {
    throw("ERROR: Trying to assign values to position " + std::to_string(pos) + " in a complex waveform of length " + std::to_string(mLength) +". ");
    exit(EXIT_FAILURE);
  }
  mReal[pos] = vreal;
  mImag[pos] = vimag;

  return;
}

//---------------------------------------------------------
void ComplexWaveform::SetValue(unsigned int pos, gsl_complex vcmplx)
{
  double re, im;
  re = GSL_REAL(vcmplx);
  im = GSL_IMAG(vcmplx);
  SetValue(pos, re, im);
  return;
}


//---------------------------------------------------------
void ComplexWaveform::SetTs(double Ts)
{
  mTs = Ts;
  return;
}

//---------------------------------------------------------
void ComplexWaveform::SetFirstLag(double first_lag)
{
  mFirstLag = first_lag;
  return;
}

//---------------------------------------------------------
// Get things
//---------------------------------------------------------
unsigned int ComplexWaveform::GetLength() const
{
  return mLength;
}

//---------------------------------------------------------
double ComplexWaveform::GetTs() const
{
  return mTs;
}

//---------------------------------------------------------
int ComplexWaveform::GetFirstLag() const
{
  return mFirstLag;
}

//---------------------------------------------------------
gsl_complex ComplexWaveform::GetValue(unsigned int index) const
{
  if (index >= mLength)
  {
    throw("ERROR: Trying to read values at index " + std::to_string(index) + " in a complex waveform of length " + std::to_string(mLength) + ".");
  }
  gsl_complex cplx;
  cplx = gsl_complex_rect(mReal[index], mImag[index]);
  return (cplx);
}

//---------------------------------------------------------
double ComplexWaveform::GetReal(unsigned int index) const
{
  if (index >= mLength)
  {
    throw("ERROR: Trying to read values at index " + std::to_string(index) + " in a complex waveform of length " + std::to_string(mLength) + ".");
  }
  double v;
  v = mReal[index];
  return v;
}

//---------------------------------------------------------
double ComplexWaveform::GetImag(unsigned int index) const
{
  if (index >= mLength)
  {
    throw("ERROR: Trying to read values at index " + std::to_string(index) + " in a complex waveform of length " + std::to_string(mLength) + ".");
  }
  double v;
  v = mImag[index];
  return v;
}

//---------------------------------------------------------
double ComplexWaveform::GetMag(unsigned int index) const
{
  if (index >= mLength)
  {
    throw("ERROR: Trying to read values at index " + std::to_string(index) + " in a complex waveform of length " + std::to_string(mLength) + ".");
  }
  double re, im, mag;
  re = mReal[index];
  im = mImag[index];
  mag = sqrt(re*re + im*im);
  return mag;
}

//---------------------------------------------------------
double ComplexWaveform::GetPhase(unsigned int index) const
{
  if (index >= mLength)
  {
    throw("ERROR: Trying to read values at index " + std::to_string(index) + " in a complex waveform of length " + std::to_string(mLength) + ".");
  }
  double re, im, phase;
  re = mReal[index];
  im = mImag[index];
  phase = atan2(im, re) * 180/M_PI;
  return phase;
}


//---------------------------------------------------------
// OPERATORS
//---------------------------------------------------------

bool ComplexWaveform::operator== (const ComplexWaveform& cw) const
{
  if (this == &cw)                        return true; // self-comparison
  if (mLength != cw.mLength)              return false;
  if (mTs != cw.mTs)                      return false;
  if (mFirstLag != cw.mFirstLag)          return false;
  for (unsigned int k=0; k<mLength; k++)
  {
    if (mReal[k] != cw.mReal[k]) return false;
    if (mImag[k] != cw.mImag[k]) return false;
  }
  return true;
}

//---------------------------------------------------------
bool ComplexWaveform::operator!= (const ComplexWaveform& cw) const
{
  return (!(*this==cw));
}


//---------------------------------------------------------
ComplexWaveform& ComplexWaveform::operator= (const ComplexWaveform &cw)
{
  // Check for self-assignment
  if (this != &cw)
  {
    FreeMemory();
    mLength      = cw.mLength;
    mTs          = cw.mTs;
    mFirstLag    = cw.mFirstLag;
    AllocateMemory(mLength);
    memcpy(mReal, cw.mReal, mLength*sizeof(*mReal));
    memcpy(mImag, cw.mImag, mLength*sizeof(*mImag));
  }
  return *this;
}

//---------------------------------------------------------
ComplexWaveform ComplexWaveform::operator+ (const ComplexWaveform& cw) const
{
  if (mLength!= cw.mLength)
  {
    throw("Cannot add ComplexWaveforms of different lengths.");
  }
  if (mTs != cw.mTs)
  {
    std::cout << "WARNING: Adding ComplexWaveforms with different sampling periods." << std::endl;
  }
  if (mFirstLag != cw.mFirstLag)
  {
    std::cout << "WARNING: Adding ComplexWaveforms with different initial lag." << std::endl;
  }
  
  ComplexWaveform cwo(mLength, mTs, mFirstLag);
  for (unsigned int k=0; k<mLength; k++)
  {
    cwo.mReal[k] = mReal[k] + cw.mReal[k];
    cwo.mImag[k] = mImag[k] + cw.mImag[k];
  }
  return cwo;
}

//---------------------------------------------------------
ComplexWaveform ComplexWaveform::operator+ (gsl_complex v) const
{
  ComplexWaveform cw(mLength, mTs, mFirstLag);
  double re = GSL_REAL(v);
  double im = GSL_IMAG(v);

  for (unsigned int k=0; k<mLength; k++)
  {
    cw.mReal[k] = mReal[k] + re;
    cw.mImag[k] = mImag[k] + im;
  }
  return cw;
}

//---------------------------------------------------------
ComplexWaveform ComplexWaveform::operator- (const ComplexWaveform& w) const
{
  if (mLength!= w.mLength)
  {
    throw("Cannot subtract ComplexWaveforms of different length.");
  }
  if (mTs != w.mTs)
  {
    std::cout << "WARNING: Subtracting ComplexWaveforms with different sampling periods." << std::endl;
  }
  if (mFirstLag != w.mFirstLag)
  {
    std::cout << "WARNING: Subtracting ComplexWaveforms with different initial lags." << std::endl;
  }
  
  ComplexWaveform cw(mLength, mTs, mFirstLag);
  for (unsigned int k=0; k<mLength; k++)
  {
    cw.mReal[k] = mReal[k] - w.mReal[k];
    cw.mImag[k] = mImag[k] - w.mImag[k];
  }
  return cw;
}

//---------------------------------------------------------
ComplexWaveform ComplexWaveform::operator- (gsl_complex v) const
{
  double re = GSL_REAL(v);
  double im = GSL_IMAG(v);
  gsl_complex  w = gsl_complex_rect(-re, -im);
  return (this->operator+(w)); 
}

//---------------------------------------------------------
ComplexWaveform ComplexWaveform::operator* (const ComplexWaveform& w) const
{
  if (mLength!= w.mLength)
  {
    throw("Cannot multiply ComplexWaveforms of different length.");
  }
  if (mTs != w.mTs)
  {
    std::cout << "WARNING: Multiplying ComplexWaveforms with different sampling periods." << std::endl;
  }
  if (mFirstLag != w.mFirstLag)
  {
    std::cout << "WARNING: Multiplying ComplexWaveforms with different first lags." << std::endl;
  }
  
  ComplexWaveform cw(mLength, mTs, mFirstLag);
  for (unsigned int k=0; k<mLength; k++)
  {
    cw.mReal[k] = mReal[k]*w.mReal[k] - mImag[k]*w.mImag[k];
    cw.mImag[k] = mImag[k]*w.mReal[k] + mReal[k]*w.mImag[k];
  }
  return cw;
}

//---------------------------------------------------------
ComplexWaveform ComplexWaveform::operator* (double a) const
{
  ComplexWaveform cw(mLength, mTs, mFirstLag);
  for (unsigned int k=0; k<cw.mLength; k++)
  {
    cw.mReal[k] = mReal[k]*a;
    cw.mImag[k] = mImag[k]*a;
  }
  return cw;
}

//---------------------------------------------------------
ComplexWaveform ComplexWaveform::operator/ (const ComplexWaveform& w) const 
{
  if (mLength!= w.mLength)
  {
    throw("Cannot divide ComplexWaveforms of different length.");
  }
  if (mTs != w.mTs)
  {
    std::cout << "WARNING: Dividing ComplexWaveforms with different sampling periods." << std::endl;
  }
  if (mFirstLag != w.mFirstLag)
  {
    std::cout << "WARNING: Dividing ComplexWaveforms with different first lags." << std::endl;
  }
  ComplexWaveform wc = w.Conj();
  ComplexWaveform cw = this->operator*(wc);
  for (unsigned int k=0; k<mLength; k++)
  {
    double div = w.mReal[k]*w.mReal[k] + w.mImag[k]*w.mImag[k];
    if (div!=0)
    {
      cw.mReal[k] = cw.mReal[k] / div;
      cw.mImag[k] = cw.mImag[k] / div;
    }
    else
    {
      if (cw.mReal[k] > 0)  cw.mReal[k] =  INFINITY;
      if (cw.mReal[k] == 0) cw.mReal[k] =  SNAN;
      if (cw.mReal[k] < 0)  cw.mReal[k] = -INFINITY;
      if (cw.mImag[k] > 0)  cw.mImag[k] =  INFINITY;
      if (cw.mImag[k] == 0) cw.mImag[k] =  SNAN;
      if (cw.mImag[k] < 0)  cw.mImag[k] = -INFINITY;
      std::cout << "WARNING: Division by zero." <<  std::endl;
    }
  }
  
  return cw;
}

//---------------------------------------------------------
ComplexWaveform ComplexWaveform::operator/ (double a) const
{
  double b = 1/a;
  return (this->operator*(b));
}


//---------------------------------------------------------
ComplexWaveform ComplexWaveform::Conj(void) const
{
  ComplexWaveform cw(mLength, mTs, mFirstLag);
  memcpy(cw.mReal, mReal, mLength*sizeof(*mReal));
  for (unsigned int k=0; k<mLength; k++)
  {
    cw.mImag[k] = -mImag[k];
  }
  return cw;
}



//---------------------------------------------------------
std::ostream& operator<< (std::ostream &out, const ComplexWaveform& w)
{
  out << std::endl;
  out << "ComplexWaveform:" << std::endl;
  out << "  Length:       " << w.mLength   << std::endl;
  out << "  Ts:           " << w.mTs       << std::endl;
  out << "  FirstLag:     " << w.mFirstLag << std::endl;
  unsigned int kmax = 10;
  std::string endstr = "...";
  if (w.mLength < kmax)
  {
    kmax = w.mLength;
    endstr = "";
  }
  for (unsigned int k=0; k<kmax; k++)
  {
    out << "cw[" << k << "]=" << w.mReal[k];
    if (w.mImag[k]>0)
    {
      out << "+j" << w.mImag[k] << "  ";
    }
    else
    {
      out << "-j" << -w.mImag[k] << "  ";
    }
  }
  out << endstr << std::endl;
  return out;
}

// VIRTUAL FUNCTIONS OF THE PARENT CLASS
//---------------------------------------------------------
ComplexWaveform& ComplexWaveform::Load(std::string fname)
{
  std::cout << "ComplexWaveform::Load() TBD." << std::endl;
  ComplexWaveform* cw = new ComplexWaveform();
  return *cw;
}

//---------------------------------------------------------
void ComplexWaveform::Save(std::string fname)
{
  std::cout << "ComplexWaveform::Save() TBD." << std::endl;
  return;
}

//---------------------------------------------------------
void ComplexWaveform::Print()
{
  std::cout << *this;
  return;
}

//---------------------------------------------------------
void ComplexWaveform::Plot(std::string fname)
{
  ComplexWaveform::PlotMagPhase(fname);
}

//---------------------------------------------------------
void ComplexWaveform::PlotMagPhase(std::string fname)
{
  unsigned int length = mLength;
  //std::cout << mLength << std::endl;
  if (length>PLT_MAX_N)
  {
    std::cout << "Very long waveform. Plotting first " << PLT_MAX_N << " lags only." << std::endl;
    length = PLT_MAX_N;
  }
  // Horizontal axis data
  std::vector<double> tau(length);
  for (unsigned int k=0; k<length; k++)
  {
    tau.at(k) = (double) ((int)k + mFirstLag) * mTs;
    //std::cout << "tau(" << k << ")= " << tau.at(k) << std::endl ;
  }
  std::vector<std::vector<double>> tautau(2);
  tautau.at(0) = tau;
  tautau.at(1) = tau;

  // Vertical axis data
  std::vector<double> mag(length);
  std::vector<double> phase(length);
  for (unsigned int k=0; k<length; k++)
  {
    mag.at(k)   = sqrt(mReal[k]*mReal[k] + mImag[k]*mImag[k]);
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

  Plots::PlotMultiple (fname, tautau, data, xlabels, ylabels, titles);

  return;
}

//---------------------------------------------------------
void ComplexWaveform::PlotPhase(std::string fname)
{
  unsigned int length = mLength;
  std::cout << mLength << std::endl;
  if (length>PLT_MAX_N)
  {
    std::cout << "Very long waveform. Plotting first " << PLT_MAX_N << " lags only." << std::endl;
    length = PLT_MAX_N;
  }
  // Horizontal axis data
  std::vector<double> tau(length);
  for (unsigned int k=0; k<length; k++)
  {
    tau.at(k) = (double) ((int)k + mFirstLag) * mTs;
    //std::cout << "tau(" << k << ")= " << tau.at(k) << std::endl ;
  }
  std::vector<std::vector<double>> tautau(1);
  tautau.at(0) = tau;

  // Vertical axis data
  std::vector<double> phase(length);
  for (unsigned int k=0; k<length; k++)
  {
    phase.at(k) = atan2(mImag[k], mReal[k]) * 180/M_PI;
  }
  std::vector<std::vector<double>>data(1);
  data.at(0) = phase;

  // Labels
  std::vector<std::string> xlabels(1);
  std::vector<std::string> ylabels(1);
  std::vector<std::string> titles(1);

  xlabels.at(0) = "Tau [s]";

  ylabels.at(0) = "Phase [deg]";

  titles.at(0) = GetName() + " - " + GetTimeTag() + "\n" + GetOrigin();

  Plots::PlotMultiple (fname, tautau, data, xlabels, ylabels, titles);

  return;
}

//---------------------------------------------------------
void ComplexWaveform::PlotMag(std::string fname)
{
  unsigned int length = mLength;
  if (length>PLT_MAX_N)
  {
    std::cout << "Very long waveform. Plotting first " << PLT_MAX_N << " lags only." << std::endl;
    length = PLT_MAX_N;
  }
  // Horizontal axis data
  std::vector<double> tau(length);
  for (unsigned int k=0; k<length; k++)
  {
    tau.at(k) = (double) ((int)k + mFirstLag) * mTs;
    //std::cout << "tau(" << k << ")= " << tau.at(k) << std::endl ;
  }
  std::vector<std::vector<double>> tautau(1);
  tautau.at(0) = tau;

  // Vertical axis data
  std::vector<double> mag(length);
  for (unsigned int k=0; k<length; k++)
  {
    mag.at(k)   = sqrt(mReal[k]*mReal[k] + mImag[k]*mImag[k]);
  }
  std::vector<std::vector<double>>data(1);
  data.at(0) = mag;

  // Labels
  std::vector<std::string> xlabels(1);
  std::vector<std::string> ylabels(1);
  std::vector<std::string> titles(1);

  xlabels.at(0) = "Tau [s]";

  ylabels.at(0) = "Magnitude";

  titles.at(0) = GetName() + " - " + GetTimeTag() + "\n" + GetOrigin();

  Plots::PlotMultiple (fname, tautau, data, xlabels, ylabels, titles);

  return;
}

//---------------------------------------------------------
void ComplexWaveform::PlotPow(std::string fname)
{
  unsigned int length = mLength;
  std::cout << mLength << std::endl;
  if (length>PLT_MAX_N)
  {
    std::cout << "Very long waveform. Plotting first " << PLT_MAX_N << " lags only." << std::endl;
    length = PLT_MAX_N;
  }
  // Horizontal axis data
  std::vector<double> tau(length);
  for (unsigned int k=0; k<length; k++)
  {
    tau.at(k) = (double) ((int)k + mFirstLag) * mTs;
    //std::cout << "tau(" << k << ")= " << tau.at(k) << std::endl ;
  }
  std::vector<std::vector<double>> tautau(1);
  tautau.at(0) = tau;

  // Vertical axis data
  std::vector<double> pow(length);
  for (unsigned int k=0; k<length; k++)
  {
    pow.at(k)   = mReal[k]*mReal[k] + mImag[k]*mImag[k];
  }
  std::vector<std::vector<double>>data(1);
  data.at(0) = pow;

  // Labels
  std::vector<std::string> xlabels(1);
  std::vector<std::string> ylabels(1);
  std::vector<std::string> titles(1);

  xlabels.at(0) = "Tau [s]";

  ylabels.at(0) = "Power";

  titles.at(0) = GetName() + " - " + GetTimeTag() + "\n" + GetOrigin();

  Plots::PlotMultiple (fname, tautau, data, xlabels, ylabels, titles);

  return;
}

//---------------------------------------------------------
void ComplexWaveform::PlotReIm(std::string fname)
{
  unsigned int length = mLength;
  std::cout << mLength << std::endl;
  if (length>PLT_MAX_N)
  {
    std::cout << "Very long waveform. Plotting first " << PLT_MAX_N << " lags only." << std::endl;
    length = PLT_MAX_N;
  }
  // Horizontal axis data
  std::vector<double> tau(length);
  for (unsigned int k=0; k<length; k++)
  {
    tau.at(k) = (double) ((int)k + mFirstLag) * mTs;
    //std::cout << "tau(" << k << ")= " << tau.at(k) << std::endl ;
  }
  std::vector<std::vector<double>> tautau(2);
  tautau.at(0) = tau;
  tautau.at(1) = tau;

  // Vertical axis data
  std::vector<double> re(mReal, mReal + length);
  std::vector<double> im(mImag, mImag + length);

  std::vector<std::vector<double>>data(2);
  data.at(0) = re;
  data.at(1) = im;

  // Labels
  std::vector<std::string> xlabels(2);
  std::vector<std::string> ylabels(2);
  std::vector<std::string> titles(2);

  xlabels.at(0) = "Tau [s]";
  xlabels.at(1) = "Tau [s]";

  ylabels.at(0) = "Real";
  ylabels.at(1) = "Imag";

  titles.at(0) = GetName() + " - " + GetTimeTag() + "\n" + GetOrigin();
  titles.at(1) = "";

  Plots::PlotMultiple (fname, tautau, data, xlabels, ylabels, titles);

  return;
}

//---------------------------------------------------------
// PRIVATE FUNCTIONS
//---------------------------------------------------------
void ComplexWaveform::AllocateMemory(unsigned int length)
{
  if (NULL!=mReal)
  {
    throw("ERROR: Trying to allocate new memory to already allocated m_real in complex waveform.");
  }

  if (NULL!=mImag)
  {
    throw("ERROR: Trying to allocate new memory to already allocated m_imag in complex waveform.");
  }

  mReal = (double*) calloc(length, sizeof(double));
  if (NULL==mReal)
  {
    throw("ERROR: Could not allocate memory for m_real in complex waveform.");
  }

  mImag = (double*) calloc(length, sizeof(double));
  if (NULL==mImag)
  {
    free(mReal);
    throw("ERROR: Could not allocate memory for m_imag in complex waveform.");
  }
  return;
}

//---------------------------------------------------------
void ComplexWaveform::FreeMemory()
{
  free(mReal);
  mReal = NULL;
  free(mImag);
  mImag = NULL;
  return;
}


//---------------------------------------------------------
void ComplexWaveform::SetLength(unsigned int length)
{
  mLength = length;
}



#endif


