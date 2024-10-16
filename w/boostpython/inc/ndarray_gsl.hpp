#ifndef NDARRAY_GSL_HPP
#define NDARRAY_GSL_HPP

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace py = boost::python;
namespace np = boost::python::numpy;

extern bool gNDARRAY_GSL_initialized;

class NdarrayGsl
{
  public:
    static void             initialize(void);
    static np::ndarray      gslvector2ndarray(gsl_vector* v);
    static np::ndarray      gslmatrix2ndarray(gsl_matrix* m);
    static gsl_vector_view  ndarray2gslvectorview(np::ndarray np);
    static gsl_matrix_view  ndarrat2gslmatrix(np::ndarray np);
};


#endif 
