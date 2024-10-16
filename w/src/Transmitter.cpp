#ifndef     TRANSMITTER_CPP_
#define     TRANSMITTER_CPP_

/**
    \file   Transmitter.cpp
    \brief  Implentation of the Transmitter class.
*/

#include "Transmitter.hpp"
#include "ObservableProcessor.hpp"
#include <iostream>


Transmitter::Transmitter()
{
    SetSourceType();
    SetAntenna();
    SetBodyFrame();
}

Transmitter::~Transmitter()
{
}

void Transmitter::SetSourceType()
{
}

void Transmitter::SetBodyFrame()
{
}

void Transmitter::SetAntenna()
{
}

void Transmitter::Run(ObservableProcessor* Proc)
{
    std::cout << " Run() Accepting processing request in Transmitter." << std::endl;
}


#endif  // TRANSMITTER_CPP_

