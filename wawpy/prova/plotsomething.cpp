#include "ComplexSignal.hpp"
#include "gsl/gsl_math.h"

int main(int argc, char *argv[])
{
  int N = 10;
  double t[N];
  double real[N];
  double imag[N];
  double Fs = 80e6;
  double Ts = 1/Fs;
  double Fc = 10e6;

  for (int n=0; n<10; n++)
  {
    t[n] = n*Ts*2*M_PI*Fc;
    real[n] = cos(t[n]);
    imag[n] = sin(t[n]);
  }
  //ComplexSignal s(N, real, imag, Fs, 123);
  ComplexSignal up = ComplexSignal::LoadFromBIBASPIRFile("/home/ribo/projects/TX_postioning/data/satpos_20200930_171633/00Hour/00Min/36Sec", true, true);
  //s.Print();
  ComplexSignal dw = ComplexSignal::LoadFromBIBASPIRFile("/home/ribo/projects/TX_postioning/data/satpos_20200930_171633/00Hour/00Min/36Sec", true, true);
  up.Plot("up-signal.png");
  dw.Plot("dw-signal.png");
  std::cout << "Computing cross-correlation." << std::endl;

  unsigned int Nsamples = 80000;
  unsigned int n0 = 50;
  int min_m = -20;
  unsigned int step_m = 1;
  int max_m = 20;

  ComplexWaveform wf = ComplexSignal::Correlation(dw, up, Nsamples ,n0, min_m, step_m, max_m, false);
  std::cout << "waveform length: " << wf.GetLength() << std::endl;
  wf.PlotMagPhase("correlation-magphase.png");
  wf.PlotReIm("correlation-reim.png");
  wf.PlotMag("correlation-mag.png");
  wf.PlotPow("correlation-pow.png");
  wf.PlotPhase("correlation-phase.png");
  return 0;
}
