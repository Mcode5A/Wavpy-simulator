#ifndef     GEOMETRY_CPP_
#define     GEOMETRY_CPP_

/**
    \file   Geometry.cpp
    \brief  Implementation of the Geometry class.
*/

#include "Geometry.hpp"
#include "ObservableProcessor.hpp"

Geometry::Geometry()
{
}

Geometry::~Geometry()
{
}

void Geometry::Compute()
{
    std::cout << "Computing Geometry..." << std::endl;
}


void Geometry::Run(ObservableProcessor* Proc)
{
    std::cout << " Run() Accepting processing request in Geometry." << std::endl;
    Proc->Process(this);
}

void Geometry::GetTx2RxDelay()
{
}

void Geometry::GetTx2Scatterer2RxDelay()
{
}

void Geometry::GetTx2RxDelayRate()     
{
}

void Geometry::GetTx2Scatterer2RxDelayRate()
{
}

void Geometry::GetScatterer2TxUnitVector()
{
}

void Geometry::GetRx2ScattererUnitVector()
{
}

void Geometry::GetTx2RxUnitVector()
{
}

void Geometry::GetTxAntennaPosition()
{
}

void Geometry::GetRxAntennaPosition()
{
}

void Geometry::GetTxAntennaVelocity()
{
}

void Geometry::GetRxAntennaVelocity()
{
}

void Geometry::GetTxAntennaAcceleration()
{
}

void Geometry::GetRxAntennaAcceleration()
{
}

void Geometry::GetTxDirectionInRxAntennaFrame()
{
}

void Geometry::GetScattererDirectionInRxAntennaFrame()
{
}

void Geometry::GetIsoDelayMapOverSurfaceGrid()
{
}

void Geometry::GetIsoDopplerMapOverSurfaceGrid()
{
}


#endif
