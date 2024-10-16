#ifndef     SCENARIO_CPP_
#define     SCENARIO_CPP_

/**
    \file   Scenario.cpp
    \brief  Implementation of the Scenario class.
*/

#include "Scenario.hpp"
#include <string.h>
#include <iostream>
#include <exception>
#include "ObservableProcessor.hpp"

Scenario::Scenario()
{
    SetTime();
    SetRx();
    SetTx();
    SetScatterer();
    SetGeometry();

    mModified=1;
    
}

Scenario::~Scenario()
{
    delete mRx;
    delete mTx;
    delete mScatterer;
    delete mGeometry;
}


void Scenario::SetName(std::string name)
{
    mName = name;
}
std::string Scenario::GetName()
{
    return mName;
}

void Scenario::SetTime()
{

}


void Scenario::SetRx()
{
    if (NULL!=mRx) // if we already have a Rx we delete it
    {
        delete mRx;
        mRx = NULL;
    }
    mRx = new Receiver;
}


void Scenario::SetTx()
{
    if (NULL!=mTx) // if we already have a Tx we delete it
    {
        delete mTx;
        mTx = NULL;
    }
    mTx = new Transmitter();
}

void Scenario::SetScatterer()
{
    if (NULL!=mScatterer) // if we already have a Scatterer we delete it
    {
        delete mScatterer;
        mScatterer = NULL;
    }
    mScatterer = new Scatterer();
}

void Scenario::SetGeometry()
{
}

void Scenario::ComputeGeometry()
{
    mGeometry->Compute();
}


void Scenario::AddObservable(Observable* Obs_in)
{
    Observable** aux;
    unsigned int n = mNumObservables+1;
    unsigned int b = n*sizeof(aux);
    try
    {
        aux = (Observable**)malloc(n*b);
        if (NULL==aux)
        {
            throw "Could not allocate memory when adding observable. Observable not added.";
        }
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << std::endl;
    }
    memcpy((void*)aux, (void*)mObservable, mNumObservables*b);
    free(mObservable);
    mObservable=aux;
    mObservable[mNumObservables] = Obs_in;
    mNumObservables++;
}



void Scenario::RemoveObservable()
{
}




void Scenario::AddRetrieval()
{
}

void Scenario::RemoveRetrieval(Retrieval*)
{
}

void Scenario::SaveScenario()
{
}

void  Scenario::LoadScenario()
{
}

void Scenario::Run(ObservableProcessor *Proc)
{
    std::cout << GetName() << ": Run() Accepting processing request in Scenario." << std::endl;
    mGeometry->Run(Proc);
    mTx->Run(Proc);
    mRx->Run(Proc);
    mScatterer->Run(Proc);
    Proc->Process(this);
}


#endif  // SCENARIO_CPP_
