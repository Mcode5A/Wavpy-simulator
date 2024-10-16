#include <iostream>

#include "ComplexSignal.hpp"

int main(int argc, char *argv[])
{
  std::cout << "prova_bibaspirfile.cpp" << std::endl;
  std::cout << argv[0] << std::endl;  

  std::string filename("/home/ribo/projects/COMSIRES/GelsaCampaign/data_samples/13Sec");
  
  std::cout << "Reading file " << filename << std::endl;

  ComplexSignal up = ComplexSignal::LoadFromBIBASPIRFile(filename, true, true, false);
  ComplexSignal dw = ComplexSignal::LoadFromBIBASPIRFile(filename, false, true, false);


  std::cout << "Computing cross-correlation." << std::endl;
  unsigned int N = 80000;
  unsigned int n0 = 8000;
  int min_m = -5;
  int step_m = 1;
  int max_m = 5;

  ComplexWaveform R = ComplexSignal::Correlation(up, dw, N, n0, min_m, step_m, max_m, false);
  R.Print();
  std::cout << "End of program." << std::endl;
  return 0;
}
