#ifndef     RETRIVAL_HPP_
#define     RETRIVAL_HPP_

/**
    \file   Retrieval.hpp
    \brief  Header file for the Retrieval abstract class.
*/

/**
    \brief  The Retrieval abstract class creates a frame in order to be
            able to define a set of child classes for the different types
            of retrievals that can be defined.
*/

class Retrieval
{
    public:
        Retrieval();
        ~Retrieval();
        void Show();

};



#endif  // RETRIEVAL_HPP_
