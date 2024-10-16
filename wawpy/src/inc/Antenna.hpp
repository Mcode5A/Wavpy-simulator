#ifndef     ANTENNA_HPP_
#define     ANTENNA_HPP_

/**
    \file   Antenna.hpp
    \brief  Header file for the Antenna class.
*/

/**
    \brief  The Antenna class allows to define an antenna (antenna pattern)
            its antenna frame, and the position/velocity/attitute of
            the antenna frame respect to a body frame.
*/

#include "ReferenceFrame.hpp"
#include <gsl/gsl_matrix.h>

class Antenna
{
    private:
        ReferenceFrame  AntennaFrame;
    public:
        Antenna();
        ~Antenna();
        void            SetPattern();
        gsl_matrix      GetPattern();
        void            SetAntennaFrame();
        ReferenceFrame  GetAntennaFrame();
};

#endif  // ANTENNA_HPP_

