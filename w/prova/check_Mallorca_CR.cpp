#include "ComplexSignal.hpp"
#include "gsl/gsl_math.h"

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cerr << "ERROR: please indicate filename to be analyzed." << std::endl;
    return 0;
  }
  std::string  fname(argv[1]);
  std::cout << "=====================================" << std::endl;
  std::cout << "  " << argv[0]                         << std::endl;
  std::cout << "-------------------------------------" << std::endl;
  std::cout << "  Checking file: " << fname << std::endl;
 
  std::cout << "    ... loading data." << std::endl;
  ComplexSignal L1R = ComplexSignal::LoadFromBIBASPIRFile(fname, true, true);
  ComplexSignal L5R = ComplexSignal::LoadFromBIBASPIRFile(fname, true, false);
  ComplexSignal L1L = ComplexSignal::LoadFromBIBASPIRFile(fname, false, true);
  ComplexSignal L5L = ComplexSignal::LoadFromBIBASPIRFile(fname, false, false);

  // Calibrate DC offset
  std::cout << "Calibrating DC offsets..." << std::endl;
  L1R = L1R - L1R.Mean();
  L5R = L5R - L5R.Mean();
  L1L = L1L - L1L.Mean();
  L5L = L5L - L5L.Mean();

  std::cout << "    ... plotting signals." << std::endl;
  L1R.Plot("plots/L1R-signal.png");
  L1L.Plot("plots/L1L-signal.png");
  L5R.Plot("plots/L5R-signal.png");
  L5L.Plot("plots/L5L-signal.png");

  std::cout << "    ... printing signals." << std::endl;
  std::cout << L1R;
  std::cout << L1L;
  std::cout << L5R;
  std::cout << L5L;

  //std::cout << "    ... computing statistics from signals." << std::endl;

  
  //return 0;

  double Fs = 80e6;
  double f_loL1 = Fs / 4 * 78.75;
  double f_loL5 = Fs / 4 * 58.8125;

  unsigned int Nsamples = 80000*1;
  unsigned int n0 = 250;
  int min_m = -250;
  unsigned int step_m = 1;
  int max_m = 250;
  
  std::cout << "Computing cross-correlation." << std::endl;
  std::cout << "  Nsamples= " << Nsamples << std::endl;
  std::cout << "  n0=       " << n0       << std::endl;
  std::cout << "  min_m=    " << min_m    << std::endl;
  std::cout << "  step_m=   " << step_m   << std::endl;
  std::cout << "  max_m=    " << max_m    << std::endl;

   
  int K = 2;
  for (int k=0; k<K; k++)
  {
    std::cout << k << "/" << K << std::endl;
    // L1
    ComplexWaveform wfL1 = ComplexSignal::Correlation(L1R, L1L, Nsamples ,n0+k*Nsamples , min_m, step_m, max_m, false);
    std::string str;
    str = "plots/correlation-magphase-L1-" + std::to_string(k) + ".png";
    wfL1.PlotMagPhase(str);
    str = "plots/correlation-magnitude-L1-" + std::to_string(k) + ".png";
    wfL1.PlotMag (str);
    // L5
    ComplexWaveform wfL5 = ComplexSignal::Correlation(L5R, L5L, Nsamples ,n0+k*Nsamples , min_m, step_m, max_m, false);
    str = "plots/correlation-magphase-L5" + std::to_string(k) + ".png";
    wfL5.PlotMagPhase(str);
    str = "plots/correlation-magnitude-L5" + std::to_string(k) + ".png";
    wfL5.PlotMag (str);
  }
  
  return 0;
}
