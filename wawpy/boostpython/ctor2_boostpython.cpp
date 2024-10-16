#include "ctor2.hpp"
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(ctor2)
{
    class_<Ctor2>("Ctor2", init<std::string>())
        .def(init<double, double>())
        .def("greet", &Ctor2::greet)
        .def("set", &Ctor2::set)
    ;

}


