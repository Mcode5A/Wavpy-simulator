#ifndef     OBSERVABLEPROCESSOR_CPP_
#define     OBSERVABLEPROCESSOR_CPP_

/**
    \file   ObservableProcessor.cpp
    \brief  Implementation of the virtual class Processor to process the data
*/

#include "ObservableProcessor.hpp"

#include <iostream>


/////////////////////////////////////// Processor

void ObservableProcessor::Show()
{
    std::cout << " --- " << mName << " --- " << std::endl;
    return;
}

/////////////////////////////////////// ProcessorA
ObservableProcessorA::ObservableProcessorA()
{
    mName = std::string("A");
}

void ObservableProcessorA::Process(Scenario* scen)
{
    Show();
    mScen = scen;
    std::cout << "The name of the Scenario is: " << mScen->GetName() << std::endl; 
    std::cout << "Changing Scenario name to 'Hola'." << std::endl;
    mScen->SetName("Hola");
}

void ObservableProcessorA::Process(Geometry* geom)
{
    mGeom = geom;
    std::cout << "Processing geometry." << std::endl;
}

void ObservableProcessorA::Process(Transmitter* trans)
{
    mTx = trans;
    std::cout << "Processing transmitter." << std::endl;
}

void ObservableProcessorA::Process(Receiver* rec)
{
    mRx = rec;
    std::cout << "Processing receiver." << std::endl;
}

void ObservableProcessorA::Process(Scatterer* scat)
{
    mScat = scat;
    std::cout << "Processing scatterer." << std::endl;
}

/////////////////////////////////////// ProcessorB
/*
ProcessorB::ProcessorB()
{
    mName = std::string("B");
}

void ProcessorB::Process(Scenario* scen)
{
    Show();
    mScen = scen;
    std::cout << "The name of the Scenario is: " << mScen->GetName() << std::endl;
}

*/
#endif  // PROCESSOR_CPP_
