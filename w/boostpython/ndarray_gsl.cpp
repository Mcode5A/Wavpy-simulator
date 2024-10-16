#ifndef NDARRAY_GSL_CPP
#define NDARRAY_GPS_CPP

#include "inc/ndarray_gsl.hpp"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <iostream>

namespace py = boost::python;
namespace np = boost::python::numpy;

bool gNDARRAY_GSL_initialized = false;

void NdarrayGsl::initialize(void)
{
  if (false == gNDARRAY_GSL_initialized)
  {
    //std::cout << "Initializing..." << std::endl;
    Py_Initialize();
    np::initialize();
    gNDARRAY_GSL_initialized=true;
  }
  else
  //{
  //  std::cout << "already initialized..." << std::endl;
  //}
  return; 
}

//////////// gslvector2ndarray //////////////
np::ndarray NdarrayGsl::gslvector2ndarray(gsl_vector* v)
{
  if (NULL == v)
  {
    throw("Cannot convert NULL gsl_vector to ndarray.");
  }
  // Initialize python and numpy 
  //NdarrayGsl::initialize();

  // Obtain size of vector
  unsigned int len;
  len = v->size;
  
  // Create the ndarray
  np::dtype   datatype = np::dtype::get_builtin<double>();
  py::tuple   shape    = py::make_tuple(len);
  py::tuple   stride   = py::make_tuple(sizeof(double));
  py::object  own;

  np::ndarray a        = np::from_data(v->data, datatype, shape, stride, own);
  return a;
}

//////////// gslmatrix2ndarray //////////////
np::ndarray NdarrayGsl::gslmatrix2ndarray(gsl_matrix* m)
{
  if (NULL == m)
  {
    throw("Cannot convert NULL gsl_matrix to ndarray.");
  }
  // Initialize python and numpy
  //NdarrayGsl::initialize();

  // Obtain the size of the matrix
  unsigned int rows, cols;
  rows = m->size1;
  cols = m->size2;

  // Create the ndarray
  np::dtype   datatype  = np::dtype::get_builtin<double>();
  py::tuple   shape     = py::make_tuple(rows,cols);
  py::tuple   stride    = py::make_tuple(sizeof(double), cols*sizeof(double));
  py::object  own;

  np::ndarray a         = np::from_data(m->data, datatype, shape, stride, own);
  return a;
}

//////////// ndarray2gslvectorview //////////////
gsl_vector_view NdarrayGsl::ndarray2gslvectorview(np::ndarray arr)
{
  if (arr.get_nd()!=1)
  {
    std::cerr << "The ndarray must have only one dimension." << std::endl;
    throw("The ndarray must have only one dimension.");
  }
  unsigned int length;
  length = arr.shape(0);

  if (length == 0)
  {
    std::cerr << "The length of the ndarray cannot be zero." << std::endl;
    throw("The length of the ndarray cannot be zero.");
  }
  // Check the data type
  std::string str1(py::extract<char const *>(py::str(arr.get_dtype())));
  std::string str2("float64");
  if (str1.compare(str2) != 0)
  {
    std::cerr << "The ndarray must be of type float64 (double), but it is "
              << py::extract<char const *>(py::str(arr.get_dtype())) << std::endl;
    throw("The ndarray must be of type float64 (double)");
  }
  //std::cout << "dtype:  " << py::extract<char const *>(py::str(arr.get_dtype())) << std::endl;
  //std::cout << "arr.get_dtype().get_itemsize()" << arr.get_dtype().get_itemsize() << std::endl;

  // Get the data
  double *data;
  data = reinterpret_cast<double *>(arr.get_data());
  
  gsl_vector_view vv;
  vv = gsl_vector_view_array(data, length);

  return vv;
}


#endif
