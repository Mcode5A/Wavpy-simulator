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
  std::cout << "Processing file: " << fname << std::endl;
  // std::string fname("/home/ribo/projects/TX_postioning/data/satpos_20200930_171633/00Hour/00Min/36Sec");
  // std::string fname("/home/ribo/projects/TX_postioning/data/satpos_20201222_143730/00Hour/00Min/18Sec");
  // std::string fname("/home/ribo/projects/TX_postioning/data/satpos_20201223_114432/00Hour/00Min/12Sec");
  // std::string fname("/home/ribo/projects/TX_postioning/data/satpos_20201223_121915/00Hour/00Min/55Sec");
  // std::string fname("/home/ribo/projects/TX_postioning/data/satpos_20201223_122533/00Hour/01Min/28Sec");
  
  ComplexSignal up = ComplexSignal::LoadFromBIBASPIRFile(fname, true, true);
  ComplexSignal dw = ComplexSignal::LoadFromBIBASPIRFile(fname, false, true);
  up.Plot("plots/up-signal.png");
  dw.Plot("plots/dw-signal.png");

  unsigned int Nsamples = 80000;
  unsigned int n0 = 50;
  int min_m = -20;
  unsigned int step_m = 1;
  int max_m = 20;
  
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
    ComplexWaveform wf = ComplexSignal::Correlation(dw, up, Nsamples ,n0+k*Nsamples , min_m, step_m, max_m, false);
    std::string str;
    str = "plots/correlation-magphase-" + std::to_string(k) + ".png";
    wf.PlotMagPhase(str);
    str = "plots/correlation-magnitude-" + std::to_string(k) + ".png";
    wf.PlotMag (str);
  }
  
  return 0;
}
