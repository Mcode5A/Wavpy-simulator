#ifndef     RAWDATA_HPP_
#define     RAWDATA_HPP_

/**
    \file   RawData.hpp
    \brief  Header file for the RawData class.
*/

/**
    \brief  RawData class (child class of Observable) allow to define and use raw data (IF-sampled data).
*/
#include "Observable.hpp"
#include <gsl/gsl_vector.h>
#include <iostream>

class RawData : public Observable
{
    private:
        unsigned int        mLength = 0;
        gsl_vector_complex* mData;
        double              mFs;
    public:
        RawData();
        ~RawData();
        void SetFs(double Fs);
        double GetFs();
        void SetLength(unsigned int Length);
        unsigned int GetLength();
        void Show();
};


#endif  // RAWDATA_HPP_

