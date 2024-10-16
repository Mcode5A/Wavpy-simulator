#ifndef     OBSERVABLE_HPP_
#define     OBSERVABLE_HPP_

/**
    \file   Observable.hpp
    \brief  Header file for the Observable virtual class.
*/

/**
    \brief  The abstract Observable class provides a general class to define 
            derived classes that fit to different types of observables.
*/

#include <string>
class Observable
{
  private:
    std::string mName;    /**<  Name of the Observable */
    std::string mTimeTag; /**<  Time tag of the Observable */
    std::string mOrigin;  /**<  Origin of the Observable (name of the campaign, simulation, etc*/
    
  protected:
    Observable(std::string name="NoName", std::string timetag="NoTimeTag", std::string origin="unknown");
  public:
    virtual ~Observable();

    void SetName(std::string name);
    void SetTimeTag(std::string timetag);
    void SetOrigin(std::string origin);

    std::string GetName() const;
    std::string GetTimeTag() const;
    std::string GetOrigin() const;


    virtual Observable& Load(std::string fname) = 0;
    virtual void Save(std::string fname) = 0;
    virtual void Print() = 0;
    virtual void Plot(const std::string) = 0;

    // TODO: Com ho faig per a  indicar que hi ha d'haver funcions virtuals per als operadors?
};


#endif  //OBSERVABLE_HPP_

