#ifndef     SOURCETYPE_HPP_
#define     SOURCETYPE_HPP_

/**
    \file   SourceType.hpp
    \brief  Header file for the SourceType class.
*/

/**
    \brief  The SourceType class allows to define the properties of a
            transmitted signal: its frequency, its autocorrelation
            function.
*/
#include <gsl/gsl_vector.h>

class SourceType
{
    private:
        double      mFreq;
        gsl_vector  mAutocorr;
    public:
        SourceType();
        ~SourceType();
        void SetFrequency();
        void SetAutocorrelation();
};

#endif  // SOURCETYPE_HPP_
