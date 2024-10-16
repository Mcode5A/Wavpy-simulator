#ifndef     RECEIVER_CPP_
#define     RECEIVER_CPP_

/**
    \file   Receiver.cpp
    \brief  Implentation of the Receiver class.
*/

#include "Receiver.hpp"
#include "ObservableProcessor.hpp"

Receiver::Receiver()
{
    SetNF();
    SetAntenna();
    SetBodyFrame();
    SetFilter();
    SetDSP();
}

Receiver::~Receiver()
{}

void Receiver::SetNF(double NF)
{
    // TODO check input value
    mNF = NF;
}

void Receiver::Run(ObservableProcessor* Proc)
{
    std::cout << "Run() Accepting processing request in Receiver." << std::endl;
    Proc->Process(this);
}

void Receiver::SetFilter()
{

}

void Receiver::SetAntenna()
{

}

void Receiver::SetBodyFrame()
{

}

void Receiver::SetDSP()
{}

#endif  // TRANSMITTER_CPP_


