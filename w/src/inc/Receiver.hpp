#ifndef     RECEIVER_HPP_
#define     RECEIVER_HPP_

/**
    \file   Receiver.hpp
    \brief  Header file for the Receiver class.
*/

/**
    \brief  The Receiver class allows to the define a GNSS-R/SoOP-R
            receiver: its antenna, its noise figure, its filter 
            shape, its trajectory and its attitude.
*/

#include "ReferenceFrame.hpp"
#include "Antenna.hpp"
#include "DigitalSignalProcessor.hpp"
#include <gsl/gsl_vector.h>

class Receiver
{
    private:
        ReferenceFrame*         mBodyFrame;
        Antenna*                mAntenna;
        double                  mNF;
        gsl_vector_complex      mFilter;
        DigitalSignalProcessor* mProcessor;
    public:
        Receiver();
        ~Receiver();
        void Run(class ObservableProcessor* Proc);
        void SetNF(double NFi=3);
        void SetFilter();
        void SetAntenna();
        void SetBodyFrame();
        void SetDSP();
};

#endif  // TRANSMITTER_HPP_


