#ifndef     SCATTERER_HPP_
#define     SCATTERER_HPP_

/**
    \file   Scatterer.hpp
    \brief  Header file for the Scatterer abstract class.
*/

/**
    \brief  The Scatterer class defines an abstrct class in order to be
            able to define a set of scatterer child classes to handle
            point scatterers, surfce scatterers and volumn scatterers.
*/

#include "ReferenceFrame.hpp"

class Scatterer
{
    private:
        ReferenceFrame mScatFrame;
    public:
       Scatterer();
       ~Scatterer();
       void Run(class ObservableProcessor* Proc);
       void SetScattererFrame(ReferenceFrame ScatFrame);
       ReferenceFrame GetScattererFrame();
};



#endif  // SCATTERER_HPP_
