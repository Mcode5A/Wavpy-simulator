#ifndef     SCENARIO_HPP_
#define     SCENARIO_HPP_

/**
    \file   Scenario.hpp
    \brief  Header file for the Scenario class.
*/

/**
    \brief  The Scenario abstract class implements a GNSS-R/SoOP-R scenario
            including one transmitter, one receiver and one scatterer.
            It allows to generate observables (simulated or real ones)
            and perform analyses.
*/

#include "Transmitter.hpp"
#include "Receiver.hpp"
#include "Scatterer.hpp"
#include "Geometry.hpp"
#include "Observable.hpp"
#include "Retrieval.hpp"
//#include "Processor.hpp"

#include <iostream>

class Scenario
{
    private: 
        std::string     mName = "NoName";
        bool            mModified=0;
        Transmitter*    mTx=NULL;
        Receiver*       mRx=NULL;
        Scatterer*      mScatterer=NULL;
        Geometry*       mGeometry=NULL;
        unsigned int    mNumObservables=0;
        Observable**    mObservable=NULL;
        unsigned int    mNumRetrievals=0;
        Retrieval*      mRetrievals=NULL;
        unsigned int    mNumTimeSamples=0; 
        double*         mTime=NULL;
        // TODO: Potser millor fer el vector de temps amb un gsl_array!
    public: 
        Scenario();
        ~Scenario();
        void SetName(std::string name);
        std::string GetName();
        void SetTime();
        void SetRx();
        void SetTx();
        void SetScatterer();
        void SetGeometry();
        void ComputeGeometry();
        void AddObservable(Observable* Obs_in);
        void RemoveObservable();
        void AddRetrieval();
        void RemoveRetrieval(Retrieval*);
        void LoadScenario();
        void SaveScenario();
        void Run(class ObservableProcessor* Proc);
};


#endif  // SCENARIO_HPP_
