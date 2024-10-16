#ifndef     ANTENNA_CPP_
#define     ANTENNA_CPP_

/**
    \file   Antenna.cpp
    \brief  Implementation of the Antenna class.
*/

#include "Antenna.hpp"

Antenna::Antenna()
{
    SetAntennaFrame();
    SetPattern();
}

Antenna::~Antenna()
{}

void Antenna::SetPattern()
{}

gsl_matrix Antenna::GetPattern()
{}

void Antenna::SetAntennaFrame()
{}

ReferenceFrame Antenna::GetAntennaFrame()
{}



#endif  //ANTENNA_CPP_
