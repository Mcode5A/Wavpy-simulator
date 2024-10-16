#ifndef DDM_CPP_
#define DDM_CPP_

/**
    \file  DDM.cpp
    \brief Contains the definition for the class DDM
*/

#include "DDM.hpp"
#include <math.h>
#include "Plots.hpp"
//---------------------------------------------------------
DDM::DDM()
{
  mNcuts = 0;
  mfcuts = NULL;
  mcuts = NULL;
}

//---------------------------------------------------------
DDM::DDM(ComplexSignal& a, ComplexSignal& b, unsigned int Nsamples, unsigned int n0, int min_m, unsigned int step_m, int max_m, double fmin, double fstep, double fmax)
{
  struct timeval tval_before, tval_after, tval_result;
  gettimeofday(&tval_before, NULL);

  // Check that input signals are compatible
  if (a.GetFs() != b.GetFs())
  {
    throw("ERROR: signals with different sampling frequency cannot be used to generate DDM.");
  }

  // take the chunk of the signals that we will need for computing the DDM
  ComplexSignal a_crop = a.Crop(n0-max_m, n0+Nsamples-min_m+1);
  ComplexSignal b_crop = b.Crop(n0-max_m, n0+Nsamples-min_m+1);

  // obtain Doppler cuts
  mNcuts = floor((fmax - fmin)/fstep) + 1;
  std::cout << "mNcuts: " << mNcuts << std::endl;
  mfcuts = (double*)calloc(sizeof(double), mNcuts);
  if (NULL == mfcuts)
  {
    throw("ERROR: could not allocate memory for doppler cut frequencies.");
  } 
  for (unsigned int k=0;k<mNcuts;k++)
  {
    mfcuts[k] = fmin + fstep*k;
    std::cout << "mfcuts[" << k << "] = " << mfcuts[k] << std::endl;
  }

  mcuts = (ComplexWaveform**)calloc(sizeof(ComplexWaveform*), mNcuts);
  if (NULL== mcuts)
  {
    throw("ERROR: could not allocate memory for Doppler cuts array.");
  }

  // Compute waveform for each Doppler cut
  std::cout << "Compute waveform for each Doppler cut" << std::endl;
  for (unsigned int k=0; k<mNcuts; k++)
  {
    std::cout << std::endl << "Cut " << k << "/" << mNcuts << std::endl;
    std::cout << "Obtain phasor..." << std::endl;
    ComplexSignal phasor = ComplexSignal::Phasor(mfcuts[k], a_crop.GetFs(), a_crop.GetLength());
    ///ComplexSignal phasor = ComplexSignal::Phasor(mfcuts[k], a.GetFs(), a.GetLength());
    std::cout << "Multiply direct signal with phasor..." << std::endl;
    ComplexSignal aux = a_crop*phasor;
    ///ComplexSignal aux = a*phasor;
    ComplexWaveform cwf = ComplexSignal::Correlation(aux, b_crop, Nsamples, n0-n0+max_m, min_m, step_m, max_m, false);
    ///ComplexWaveform cwf = ComplexSignal::Correlation(aux, b, Nsamples, n0, min_m, step_m, max_m, false);
    ComplexWaveform* cwf2 = new ComplexWaveform(cwf);;
    mcuts[k] = cwf2;
  }
  gettimeofday(&tval_after, NULL);
  timersub(&tval_after, &tval_before, &tval_result);
  std::cout << "Time to compute DDM: " << (long int)tval_result.tv_sec << "." << (long int)tval_result.tv_usec << std::endl;
}

//---------------------------------------------------------
DDM::~DDM()
{
  free(mfcuts);
  for (unsigned int k=0; k<mNcuts; k++)
  {
    delete(mcuts[k]);
  }
  free(mcuts);
}



//---------------------------------------------------------
DDM& DDM::Load(std::string fname)
{
  std::cout << "DDM::Load() TBD." << std::endl;
  DDM* ddm = new DDM();
  return *ddm;
}

//---------------------------------------------------------
void DDM::Save(std::string fname)
{
  std::cout << "DDM::Save() TBD." << std::endl;
  return;
}

//---------------------------------------------------------
void DDM::Print()
{
  std::cout << "DDM::Print() TBD." << std::endl;
  return;
}

//---------------------------------------------------------
void DDM::Plot(const std::string fname)
{
  std::cout << "DDM::Plot() TBD." << std::endl;
  std::vector<std::vector<double>> x, y, z;
//  for (double i=-5; i <= 5; i += 0.25) 
//  {
//    std::vector<double> x_row, y_row, z_row;
//    for (double j=-5; j<=5; j+= 0.25)
//    {
//      x_row.push_back(i);
//      y_row.push_back(j);
//      z_row.push_back(::std::sin(::std::hypot(i, j)));
//    } 
//    x.push_back(x_row);
//    y.push_back(y_row);
//    z.push_back(
//  }

  unsigned int nLags = mcuts[0]->GetLength();

  for (unsigned int k=0; k<mNcuts; k++)
  {
    std::vector<double> x_row, y_row, z_row; 
    for (unsigned int n=0; n < nLags; n++)
    {
      x_row.push_back(mfcuts[k]);
      y_row.push_back(mcuts[k]->GetFirstLag()+n*mcuts[k]->GetTs()*3e8);
      //z_row.push_back(10*log10(mcuts[k]->GetMag(n)));
      z_row.push_back(mcuts[k]->GetMag(n));
    }
    x.push_back(x_row);
    y.push_back(y_row);
    z.push_back(z_row);
  }

  std::string xlabel = "delay";
  std::string ylabel = "Doppler";
  std::string zlabel = "";
  std::string title = "DDM";
 
  Plots::PlotImage(fname, x, y, z, xlabel, ylabel, zlabel, title);


  return;
}

#endif
