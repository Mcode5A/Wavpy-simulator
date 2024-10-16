#include "ReferenceFrame.hpp"

#include "inc/ndarray_gsl.hpp"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <iostream>

namespace py = boost::python;
namespace np = boost::python::numpy;

///////////////////////////////////////////////
// Define a child class of the main one, which
// is a wrapper class for specific methods
// conversion from gsl to ndarray, and viceversa,...
//
class ReferenceFrameWrap : public ReferenceFrame
{
  public:
    int testfun_wrap(int a, int b);
    np::ndarray GetMyReferenceFrame_wrap();
    //int GetMyReferenceFrame_wrap();
    np::ndarray GetPosition_wrap(unsigned int L=0) const;
    np::ndarray GetAttitude_wrap(unsigned int L=0) const;
    np::ndarray GetVelocity_wrap(unsigned int L=0) const;
    np::ndarray PosInParentFrame2PosInBF_wrap(np::ndarray PosInParent, unsigned int L=0) const;
    np::ndarray PosInBF2PosInParentFrame_wrap(np::ndarray PosInBF, unsigned int L=0) const;
    np::ndarray PosInBF2PosInECEF_wrap(np::ndarray PosInBF, unsigned int L=0) const;
    np::ndarray PosInECEF2PosInBF_wrap(np::ndarray PosInBF, unsigned int L=0) const;
};

int ReferenceFrameWrap::testfun_wrap(int a, int b)
{
  int retval;
  retval = this->testfun(a+b);
  return retval;
}

np::ndarray ReferenceFrameWrap::GetMyReferenceFrame_wrap()
{
  Py_Initialize();
  np::initialize();

  py::tuple   shape    = py::make_tuple(3,4);
  np::dtype   datatype = np::dtype::get_builtin<double>();
  np::ndarray a        = np::zeros(shape, datatype);
 
  return a;
}

/////////////////////////////// GetPosition_wrap ///////////////
np::ndarray ReferenceFrameWrap::GetPosition_wrap(unsigned int L) const
{
  gsl_vector *Pos;
  Pos = GetPosition(L);
  // Convert to np::ndarray
  np::ndarray Pos_ndarray = NdarrayGsl::gslvector2ndarray(Pos);
  return Pos_ndarray;
}

/////////////////////////////// GetAttitude_wrap ///////////////
np::ndarray ReferenceFrameWrap::GetAttitude_wrap(unsigned int L) const
{
  gsl_matrix *Att;
  Att = GetAttitude(L);
  // convert to np::ndarray
  np::ndarray Att_ndarray = NdarrayGsl::gslmatrix2ndarray(Att);
  return Att_ndarray;
}

/////////////////////////////// GetVelocity_wrap ///////////////
np::ndarray ReferenceFrameWrap::GetVelocity_wrap(unsigned int L) const
{
  gsl_vector *Vel;
  Vel = GetVelocity(L);
  // convert to np::ndarray
  np::ndarray Vel_ndarray = NdarrayGsl::gslvector2ndarray(Vel);
  return Vel_ndarray;
}

/////////////////////////////// PosInParentFrame2PosInBF_wrap //
np::ndarray ReferenceFrameWrap::PosInParentFrame2PosInBF_wrap(np::ndarray PosInParent, unsigned int L) const
{
  // Define gsl_vector pointers to use the actual library function
  gsl_vector_view PosInParent_gslview;
  gsl_vector_view PosInBF_gslview;
 
  // convert the input ndarray to gsl_vector
  PosInParent_gslview = NdarrayGsl::ndarray2gslvectorview(PosInParent);
  /*for (unsigned int k=0;k<(&PosInParent_gsl.vector)->size;k++)
    std::cout << "vector[" << k << "]= " << (&PosInParent_gsl.vector)->data[k] << std::endl;*/

  gsl_vector* v;
  v = &PosInBF_gslview.vector;
  v = PosInParentFrame2PosInBF(&PosInParent_gslview.vector, L);

  np::ndarray PosInBF = NdarrayGsl::gslvector2ndarray(v);

  return PosInBF;
}

/////////////////////////////// PosInBF2PosInParentFrame_wrap //
np::ndarray ReferenceFrameWrap::PosInBF2PosInParentFrame_wrap(np::ndarray PosInBF, unsigned int L) const
{
  //TODO
  // Define gsl_vector pointers to use the actual library function
  gsl_vector_view PosInParent_gslview;
  gsl_vector_view PosInBF_gslview;
 
  // convert the input ndarray to gsl_vector
  PosInBF_gslview = NdarrayGsl::ndarray2gslvectorview(PosInBF);
  /*for (unsigned int k=0;k<(&PosInParent_gsl.vector)->size;k++)
    std::cout << "vector[" << k << "]= " << (&PosInParent_gsl.vector)->data[k] << std::endl;*/

  gsl_vector* v;
  v = &PosInParent_gslview.vector;
  v = PosInBF2PosInParentFrame(&PosInBF_gslview.vector, L);

  np::ndarray PosInParent = NdarrayGsl::gslvector2ndarray(v);

  return PosInParent;
}


/////////////////////////////// PosInBF2PosInParentFrame_wrap //
np::ndarray ReferenceFrameWrap::PosInBF2PosInECEF_wrap(np::ndarray PosInBF, unsigned int L) const
{
  //TODO
  // Define gsl_vector pointers to use the actual library function
  gsl_vector_view PosInECEF_gslview;
  gsl_vector_view PosInBF_gslview;
 
  // convert the input ndarray to gsl_vector
  PosInBF_gslview = NdarrayGsl::ndarray2gslvectorview(PosInBF);
  /*for (unsigned int k=0;k<(&PosInParent_gsl.vector)->size;k++)
    std::cout << "vector[" << k << "]= " << (&PosInParent_gsl.vector)->data[k] << std::endl;*/

  gsl_vector* v;
  v = &PosInECEF_gslview.vector;
  v = PosInBF2PosInECEF(&PosInBF_gslview.vector, L);

  np::ndarray PosInECEF = NdarrayGsl::gslvector2ndarray(v);

  return PosInECEF;
}

/////////////////////////////// PosInECEF2PosInBF_wrap //
np::ndarray ReferenceFrameWrap::PosInECEF2PosInBF_wrap(np::ndarray PosInECEF, unsigned int L) const
{
  //TODO
  // Define gsl_vector pointers to use the actual library function
  gsl_vector_view PosInECEF_gslview;
  gsl_vector_view PosInBF_gslview;
 
  // convert the input ndarray to gsl_vector
  PosInECEF_gslview = NdarrayGsl::ndarray2gslvectorview(PosInECEF);
  /*for (unsigned int k=0;k<(&PosInParent_gsl.vector)->size;k++)
    std::cout << "vector[" << k << "]= " << (&PosInParent_gsl.vector)->data[k] << std::endl;*/

  gsl_vector* v;
  v = &PosInBF_gslview.vector;
  v = PosInECEF2PosInBF(&PosInECEF_gslview.vector, L);

  np::ndarray PosInBF = NdarrayGsl::gslvector2ndarray(v);

  return PosInBF;
}


// Automatic wrapper for overloaded functions
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ReferenceFrameWrap_GetPosition_wrap_overloads, GetPosition_wrap, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ReferenceFrameWrap_GetAttitude_wrap_overloads, GetAttitude_wrap, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ReferenceFrameWrap_GetVelocity_wrap_overloads, GetVelocity_wrap, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ReferenceFrameWrap_PosInParentFrame2PosInBF_wrap_overloads, PosInParentFrame2PosInBF_wrap, 1, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ReferenceFrameWrap_PosInBF2PosInParentFrame_wrap_overloads, PosInBF2PosInParentFrame_wrap, 1, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ReferenceFrameWrap_PosInBF2PosInECEF_wrap_overloads, PosInBF2PosInECEF_wrap, 1, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ReferenceFrameWrap_PosInECEF2PosInBF_wrap_overloads, PosInECEF2PosInBF_wrap, 1, 2);


BOOST_PYTHON_MODULE(referenceframe)
{
  //Initialize boost::python and boost::python::numpy
  NdarrayGsl::initialize();

  // expose the actual C++ class as if it would be a virtual class:
  //  add the boost::noncopyable 
  //  add the no_init
  //  put the underscore before the class (my naming convention)
  //  And define all access methods
  //  As the class is virtual as seen from python, it cannot be instatiated.
  //  Then,... (look below)
  py::class_<ReferenceFrame, boost::noncopyable>("_ReferenceFrame", py::no_init)
    .def(py::init<ReferenceFrame*, bool, gsl_matrix*, gsl_matrix*, gsl_matrix*, std::string>())
    .def("IsSolidary"              , &ReferenceFrame::IsSolidary)
    .def("IsMainReferenceFrame"    , &ReferenceFrame::IsMainReferenceFrame)
    //.def("GetMyReferenceFrame"     , &ReferenceFrame::GetMyReferenceFrame)
    .def("GetLength"               , &ReferenceFrame::GetLength)
    //.def("GetPosition"             , &ReferenceFrame::GetPosition) // in wrapper
    //.def("GetAttitude"             , &ReferenceFrame::GetAttitude) // in wrapper 
    //.def("GetVelocity"             , &ReferenceFrame::GetVelocity) // in wrapper
    .def("Print"                   , &ReferenceFrame::Print)
    //.def("PosInParentFrame2PosInBF", &ReferenceFrame::PosInParentFrame2PosInBF) // in wrapper
    //.def("PosInBF2PosInParentFrame", &ReferenceFrame::PosInBF2PosInParentFrame) // in wrapper
    //.def("PosInBF2PosInECEF"       , &ReferenceFrame::PosInBF2PosInECEF)        // in wrapper
    //.def("PosInECEF2PosInBF"       , &ReferenceFrame::PosInECEF2PosInBF)        // in wrapper
    .def("TestFun"                 , &ReferenceFrame::testfun)  // to be remove!
  ;
  
  // expose a wrapper (child) class, but using the name of the actual parent class
  // expose only the wrapped methods, but using the name of the original functions 
  py::class_<ReferenceFrameWrap, py::bases<ReferenceFrame> >("ReferenceFrame")
    .def("TestFun",                  &ReferenceFrameWrap::testfun_wrap)
    .def("GetMyReferenceFrame",      &ReferenceFrameWrap::GetMyReferenceFrame_wrap)
    .def("GetPosition",              &ReferenceFrameWrap::GetPosition_wrap,ReferenceFrameWrap_GetPosition_wrap_overloads())
    .def("GetAttitude",              &ReferenceFrameWrap::GetAttitude_wrap,ReferenceFrameWrap_GetAttitude_wrap_overloads())
    .def("GetVelocity",              &ReferenceFrameWrap::GetVelocity_wrap,ReferenceFrameWrap_GetVelocity_wrap_overloads())
    .def("PosInParentFrame2PosInBF", &ReferenceFrameWrap::PosInParentFrame2PosInBF_wrap, ReferenceFrameWrap_PosInParentFrame2PosInBF_wrap_overloads())
    .def("PosInBF2PosInParentFrame", &ReferenceFrameWrap::PosInBF2PosInParentFrame_wrap, ReferenceFrameWrap_PosInBF2PosInParentFrame_wrap_overloads())
    .def("PosInBF2PosInECEF",        &ReferenceFrameWrap::PosInBF2PosInECEF_wrap,        ReferenceFrameWrap_PosInBF2PosInECEF_wrap_overloads())
    .def("PosInECEF2PosInBF",        &ReferenceFrameWrap::PosInECEF2PosInBF_wrap,        ReferenceFrameWrap_PosInECEF2PosInBF_wrap_overloads())
  ;
};




