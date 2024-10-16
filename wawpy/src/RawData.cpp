#ifndef     RAWDATA_CPP_
#define     RAWDATA_CPP_

/**
    \file   RawData.cpp
    \brief  Implementation of the RawData class to define and operate on raw 
            data. It is a child class of Observable.
*/

#include "RawData.hpp"

void RawData::SetFs(double Fs)
{
    mFs = Fs;
    return;
} 

double RawData::GetFs(void)
{
    return mFs;
}

void RawData::SetLength(unsigned int Length)
{
    mLength = Length;
    return;
}

unsigned int RawData::GetLength()
{
    return mLength;
}

void RawData::Show()
{
    std::cout << "Fs="     << mFs     << std::endl;
    std::cout << "Length=" << mLength << std::endl;
}



#endif  // RAWDATA_CPP_
