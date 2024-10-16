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
  ComplexSignal A1 = ComplexSignal::LoadFromBIBASPIRFile(fname, true, true);
  ComplexSignal A2 = ComplexSignal::LoadFromBIBASPIRFile(fname, true, false);
  ComplexSignal B1 = ComplexSignal::LoadFromBIBASPIRFile(fname, false, true);
  ComplexSignal B2 = ComplexSignal::LoadFromBIBASPIRFile(fname, false, false);

  // Calibrate DC offset
  std::cout << "Calibrating DC offsets..." << std::endl;
  A1 = A1 - A1.Mean();
  A2 = A2 - A2.Mean();
  B1 = B1 - B1.Mean();
  B2 = B2 - B2.Mean();

  std::cout << "    ... plotting signals." << std::endl;
  A1.Plot("plots/A1-signal.png");
  A2.Plot("plots/A2-signal.png");
  B1.Plot("plots/B1-signal.png");
  B2.Plot("plots/B2-signal.png");

  std::cout << "    ... printing signals." << std::endl;
  std::cout << A1;
  std::cout << A2;
  std::cout << B1;
  std::cout << B2;

  //std::cout << "    ... computing statistics from signals." << std::endl;

  
  //return 0;

  unsigned int Nsamples = 80000;
  unsigned int n0 = 50;
  int min_m = -50;
  unsigned int step_m = 1;
  int max_m = 50;
  
  std::cout << "Computing cross-correlation." << std::endl;
  std::cout << "  Nsamples= " << Nsamples << std::endl;
  std::cout << "  n0=       " << n0       << std::endl;
  std::cout << "  min_m=    " << min_m    << std::endl;
  std::cout << "  step_m=   " << step_m   << std::endl;
  std::cout << "  max_m=    " << max_m    << std::endl;

   
  int K = 10;
  for (int k=0; k<K; k++)
  {
    std::cout << k << "/" << K << std::endl;
    // band 1
    ComplexWaveform wf1 = ComplexSignal::Correlation(A1, B1, Nsamples ,n0+k*Nsamples , min_m, step_m, max_m, false);
    std::string str;
    str = "plots/correlation-magphase-band1-" + std::to_string(k) + ".png";
    wf1.PlotMagPhase(str);
    str = "plots/correlation-magnitude-band1-" + std::to_string(k) + ".png";
    wf1.PlotMag (str);
    // band 2
    ComplexWaveform wf2 = ComplexSignal::Correlation(A2, B2, Nsamples ,n0+k*Nsamples , min_m, step_m, max_m, false);
    str = "plots/correlation-magphase-band2-" + std::to_string(k) + ".png";
    wf2.PlotMagPhase(str);
    str = "plots/correlation-magnitude-band2-" + std::to_string(k) + ".png";
    wf2.PlotMag (str);
  }
  
  return 0;
}
