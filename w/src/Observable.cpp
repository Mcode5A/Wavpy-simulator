#ifndef     OBSERVABLE_CPP_
#define     OBSERVABLE_CPP_

/**
    \file   Observable.cpp
    \brief  Implementation of the Observable virtual class.
*/

#include "Observable.hpp"

//Observable::Observable()
//{
//}

//----------------------------------------------------------
Observable::Observable(std::string name, std::string timetag, std::string origin)
{
  mName = name;
  mTimeTag = timetag;
  mOrigin = origin;
}

//----------------------------------------------------------
Observable::~Observable()
{
}

//----------------------------------------------------------
void Observable::SetName(std::string name)
{
  mName = name;
  return;
}

//----------------------------------------------------------
void Observable::SetTimeTag(std::string timetag)
{
  mTimeTag = timetag;
  return;
}

//----------------------------------------------------------
void Observable::SetOrigin(std::string origin)
{
  mOrigin = origin;
  return;
}

//----------------------------------------------------------
std::string Observable::GetName() const
{
  return mName;
}

//----------------------------------------------------------
std::string Observable::GetTimeTag() const
{
  return mTimeTag;
}

//----------------------------------------------------------
std::string Observable::GetOrigin() const
{
  return mOrigin;
}



#endif  //OBSERVABLE_CPP_
