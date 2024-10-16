#ifndef     OBSERVABLEPROCESSOR_HPP_
#define     OBSERVABLEPROCESSOR_HPP_


/**
    \file   ObservableProcessor.hpp
    \brief  Header file for the ObservableProcessor virtual class
*/


/**
    \brief  The ObservableProcessor virtual class is a _visitor_ class to implement
            the signal processing (simulation) of a scenario in order to generate
            observables.
*/

#include "Scenario.hpp"
#include "Geometry.hpp"
#include "Transmitter.hpp"
#include "Receiver.hpp"
#include "Scatterer.hpp"

#include <iostream>

class ObservableProcessor
{
    protected:
        std::string mName;
    public:
        virtual void Process(Scenario* scen) = 0;
        virtual void Process(Geometry* geom) = 0;
        virtual void Process(Transmitter* trans) = 0;
        virtual void Process(Receiver* rec) = 0;
        virtual void Process(Scatterer* scat) = 0;    
        void Show();
};

class ObservableProcessorA : public ObservableProcessor
{
    private:
        Scenario*       mScen;
        Geometry*       mGeom;
        Transmitter*    mTx;
        Receiver*       mRx;
        Scatterer*      mScat;
    public:
        ObservableProcessorA();
        ~ObservableProcessorA();
        void Process(Scenario* scen);
        void Process(Geometry* geom);
        void Process(Transmitter* trans);
        void Process(Receiver* rec);
        void Process(Scatterer* scat);
};

/*
class ProcessorB : public Processor
{
    private:
        Scenario *mScen;
        Geometry *mGeom;
    public:
        ProcessorB();
        ~ProcessorB();
        void Process(Scenario* scen);
};
*/

#endif  // PROCESSOR_HPP_
