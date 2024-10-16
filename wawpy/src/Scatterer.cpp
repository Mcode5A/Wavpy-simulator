#ifndef     SCATTERER_CPP_
#define     SCATTERER_CPP_

/**
    \file   Scatterer.cpp
    \brief  Implemntation of the Scatterer abstract class.
*/

#include "Scatterer.hpp"
#include "ObservableProcessor.hpp"
#include <iostream>

Scatterer::Scatterer()
{
    // SetReferenceFrame()
}

Scatterer::~Scatterer()
{}

void Scatterer::Run(ObservableProcessor* Proc)
{
    std::cout << " Run() Accepting processing request in Scatterer." << std::endl;
    Proc->Process(this);
}

void Scatterer::SetScattererFrame(ReferenceFrame ScatFrame)
{}

ReferenceFrame Scatterer::GetScattererFrame()
{}


#endif  // SCATTERER_CPP_
