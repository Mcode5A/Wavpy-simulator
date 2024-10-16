#ifndef     TRANSMITTER_HPP_
#define     TRANSMITTER_HPP_

/**
    \file   Transmitter.hpp
    \brief  Header file for the Transmitter class.
*/

/**
    \brief  The Transmitter class allows to the define a GNSS-R/SoOP-R
            transmitter: its signal structure, its transmitted power,
            its transmitter antenna diagram, its trajectory and its
            attitude.
*/
#include "ReferenceFrame.hpp"
#include "Antenna.hpp"
#include "SourceType.hpp"

class Transmitter
{
    private:
        ReferenceFrame* mBodyFrame;
        Antenna*        mAntenna;
        SourceType*     mSrcType;
    public:
        Transmitter();
        ~Transmitter();
        void SetSourceType();
        void SetBodyFrame();
        void SetAntenna();
        void Run(class ObservableProcessor* Proc);
};

#endif  // TRANSMITTER_HPP_

