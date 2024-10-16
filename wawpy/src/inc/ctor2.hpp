#include <string>

class Ctor2
{
  private:
    std::string mMsg;
  public:
    Ctor2(std::string msg);
    Ctor2(double x, double y);
    void set(std::string);
    std::string greet();
};
