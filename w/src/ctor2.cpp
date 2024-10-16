#include "ctor2.hpp"


#include <string>
#include <sstream>

Ctor2::Ctor2(std::string msg)
{
  mMsg = msg;
}

Ctor2::Ctor2(double x,double y)
{
  std::stringstream os;
  os << x << ":" << y << std::ends;
  set(os.str());
}

void Ctor2::set(std::string msg) 
{
  mMsg = msg;
}

std::string Ctor2::greet()
{
  return mMsg;
};


