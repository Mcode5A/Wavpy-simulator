#ifndef DDM_HPP_
#define DDM_HPP_

/**
    \file  DDM.hpp
    \brief Public header file for the DDM class.
*/

#include "Observable.hpp"
#include "ComplexWaveform.hpp"
#include "ComplexSignal.hpp"

class DDM : public Observable
{
  // Member variables
#ifndef TEST_DDM
  private:
#else
  public:
#endif
  unsigned int      mNcuts;   /**< number of Doppler cuts */
  double*           mfcuts;   /**< pointer holding the frequency of the Doppler cuts */
  ComplexWaveform** mcuts;    /**< pointer array holding the Doppler cuts */


  // Public member functions
  public:
    // Constructors & destructor
    DDM();
    DDM(ComplexSignal& a, ComplexSignal &b, unsigned int Nsamples, unsigned int n0, int min_m, unsigned int step_m, int max_m, double fmin, double fstep, double fmax);
    ~DDM();


    // Virtual functions of the parent class
    DDM& Load(std::string fname);
    void Save(std::string fname);
    void Print();
    void Plot(const std::string fname="");
    


};

#endif
