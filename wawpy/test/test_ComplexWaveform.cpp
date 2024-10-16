#define BOOST_TEST_MODULE test_ComplexWaveform
#include <boost/test/unit_test.hpp>

#define TEST_COMPLEXWAVEFORM

//#include <math.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <ComplexWaveform.hpp>


BOOST_AUTO_TEST_SUITE (test_ComplexWaveform)

BOOST_AUTO_TEST_CASE (test_ComplexWaveform_constructor)
{
  unsigned int len = 45;
  double       Ts = 40e6; 
  long int     init_sample = 345;

  // Constructor without data
  {
    {
      // Constructor without input data, one input parameter
      ComplexWaveform* cs = new ComplexWaveform(len);
      BOOST_CHECK_MESSAGE(cs->GetLength()   == len,     "Result 1.1");
      BOOST_CHECK_MESSAGE(cs->GetTs()       == 12.5e-9, "Result 1.2");
      BOOST_CHECK_MESSAGE(cs->GetFirstLag() == 0,       "Result 1.3");
      delete cs;
    }
    {
      // Constructor without input data, two input parameters
      ComplexWaveform* cs = new ComplexWaveform(len, Ts);
      BOOST_CHECK_MESSAGE(cs->GetLength()        == len, "Result 1.4");
      BOOST_CHECK_MESSAGE(cs->GetTs()            == Ts,  "Result 1.5");
      BOOST_CHECK_MESSAGE(cs->GetFirstLag()      == 0,   "Result 1.6");
      delete cs;
    }
    {
      // Constructor without input data, three input parameters
      ComplexWaveform* cs = new ComplexWaveform(len, Ts, init_sample);
      BOOST_CHECK_MESSAGE(cs->GetLength()        == len,         "Result 1.7");
      BOOST_CHECK_MESSAGE(cs->GetTs()            == Ts,          "Result 1.8");
      BOOST_CHECK_MESSAGE(cs->GetFirstLag() == init_sample, "Result 1.9");
      delete cs;
    }
  }

  // Constructor with data, all input parameters
  {
    double *re, *im;
    re = (double*) malloc(sizeof(double)*len);
    im = (double*) malloc(sizeof(double)*len);

    // Creating data
    for (unsigned int k=0; k<len; k++)
    {
      re[k] = k;
      im[k] = len - 2*k;
    }
    ComplexWaveform* cs = new ComplexWaveform(len, re, im, Ts, init_sample);
    BOOST_CHECK_MESSAGE(cs->GetLength()        == len,         "Result 1.10");
    BOOST_CHECK_MESSAGE(cs->GetTs()            == Ts,          "Result 1.11");
    BOOST_CHECK_MESSAGE(cs->GetFirstLag() == init_sample, "Result 1.12");

    for (unsigned int k=0; k<len; k++)
    {
      BOOST_CHECK_MESSAGE(cs->mReal[k] == k,         "Result 1.13");
      BOOST_CHECK_MESSAGE(cs->mImag[k] == len - 2*k, "Result 1.14");
    }
    delete cs;
    free(re);
    free(im);
  }

  // Copy constructor and comparisons 
  {
    double *re, *im;
    re = (double*) malloc(sizeof(double)*len);
    im = (double*) malloc(sizeof(double)*len);
    
    // Creating data
    for (unsigned int k=0; k<len; k++)
    {
      re[k] = 2*k;
      im[k] = len+k/3;
    }
    ComplexWaveform cs1(len, re, im, Ts, init_sample);
    // Check self-comaparison
    BOOST_CHECK_MESSAGE(cs1==cs1, "Result 1.15");
    // Check when equal
    ComplexWaveform cs2(cs1);
    BOOST_CHECK_MESSAGE(cs2==cs1, "Result 1.16");
    // Check when different Ts
    ComplexWaveform cs3(cs1);
    cs3.SetTs(Ts*2);
    BOOST_CHECK_MESSAGE(!(cs3==cs1), "Result 1.17");
    // Check when different InitSample
    ComplexWaveform cs4(cs1);
    cs4.SetFirstLag(init_sample*2);
    BOOST_CHECK_MESSAGE(!(cs4==cs1), "Result 1.18");
    // Check when different real data
    ComplexWaveform cs5(cs1);
    cs5.mReal[cs5.mLength-1] = cs5.mLength-1;
    BOOST_CHECK_MESSAGE(!(cs5==cs1), "Result 1.19");
    // Check when different imag data
    ComplexWaveform cs6(cs1);
    cs6.mImag[cs6.mLength-1] = 0;
    BOOST_CHECK_MESSAGE(!(cs6==cs1), "Result 1.20");
    // Check when different length
    ComplexWaveform cs7(cs1);
    cs7.SetLength(len-1);
    BOOST_CHECK_MESSAGE(!(cs7==cs1), "Result 1.21");
    // Check !=
    BOOST_CHECK_MESSAGE(cs3!=cs1, "Result 1.22");
    free(re);
    free(im);
  }
}


BOOST_AUTO_TEST_CASE (test_ComplexWaveform_assignment_operator)
{
  unsigned int Length = 12;
  double       Ts = 10e6;
  long int     InitialSample = 0;
  // Assignment operator
  ComplexWaveform s1(Length, Ts, InitialSample);
  ComplexWaveform s2(2*Length, 2*Ts, InitialSample+1);
  // Check that both signals are different
  BOOST_CHECK_MESSAGE(s2!=s1, "Result 2.1");
  // Assign 
  s2 = s1;
  // Check that both signals are equal
  BOOST_CHECK_MESSAGE(s2==s1, "Result 2.2");
  // Self assign
  s2 = s2;
  // Check that both signals are still equal
  BOOST_CHECK_MESSAGE(s2==s1, "Result 2.3");
}

BOOST_AUTO_TEST_CASE (test_ComplexWaveform_addition_operator)
{
  unsigned int Length = 128;
  double Ts = 10e6;
  long int InitialSample = 445;
  double re[Length], im[Length];
  for(int k=0;k<(int)Length;k++)
  {
    re[k] =  k;
    im[k] = -k;
  }
  ComplexWaveform s1(Length, re, im, Ts, InitialSample);
  ComplexWaveform s2(Length, Ts, InitialSample);
  // Addition
  ComplexWaveform s3 = s1 + s2;
  BOOST_CHECK_MESSAGE(s3==s1, "Result 3.1");
  // Force throw
  ComplexWaveform s4(2*Length, Ts, InitialSample);
  BOOST_CHECK_THROW(s4+s1, const char*);

  gsl_complex cmplx = gsl_complex_rect(2, 1);
  s3 = s1 + cmplx;
  for (int k=0; k< (int) Length; k++)
  {
    BOOST_CHECK_MESSAGE(s1.mReal[k] + GSL_REAL(cmplx) == s3.mReal[k], "Result 3.2");
    BOOST_CHECK_MESSAGE(s1.mImag[k] + GSL_IMAG(cmplx) == s3.mImag[k], "Result 3.2");
  }
}

BOOST_AUTO_TEST_CASE (test_ComplexWaveform_subtraction_operator)
{
  unsigned int Length = 128;
  double Ts = 10e6;
  long int InitialSample = 445;
  double a[Length], b[Length];
  for(int k=0;k< (int)Length;k++)
  {
    a[k] =  k;
    b[k] = -k;
  }
  ComplexWaveform s1(Length, a, b, Ts, InitialSample);
  ComplexWaveform s2(Length, a, b, Ts, InitialSample);
  ComplexWaveform s3(Length, Ts, InitialSample);
  // Subtraction
  ComplexWaveform s4 = s1-s2;
  BOOST_CHECK_MESSAGE(s4==s3, "Result 4.1");
  // Force throw
  ComplexWaveform s5(2*Length, Ts, InitialSample);
  BOOST_CHECK_THROW(s5+s1, const char*);

  gsl_complex cmplx = gsl_complex_rect(3, -1);
  s3 = s1 - cmplx;
  for (int k=0; k< (int) Length; k++)
  {
    BOOST_CHECK_MESSAGE(s1.mReal[k] - GSL_REAL(cmplx) == s3.mReal[k], "Result 4.2");
    BOOST_CHECK_MESSAGE(s1.mImag[k] - GSL_IMAG(cmplx) == s3.mImag[k], "Result 4.2");
  }
}

BOOST_AUTO_TEST_CASE (test_ComplexWaveform_multiplication_operator)
{
  unsigned int Length = 12;
  double Ts = 10e6;
  long int InitialSample = 445;
  double a[Length], b[Length], c[Length], d[Length];
  for(int k=0;k<(int)Length;k++)
  {
    a[k] =  k;
    b[k] =  2*k;
    c[k] =  k-5;
    d[k] = 2*k-5;
  }
  ComplexWaveform s1(Length, a, b, Ts, InitialSample);
  ComplexWaveform s2(Length, c, d, Ts, InitialSample);
  // Subtraction
  ComplexWaveform s3 = s1*s2;
  for (unsigned int k=0;k<Length;k++)
  {
    double re = a[k]*c[k] - b[k]*d[k];
    double im = b[k]*c[k] + a[k]*d[k];
    BOOST_CHECK_MESSAGE(s3.mReal[k]==re, "Result 5.1");
    BOOST_CHECK_MESSAGE(s3.mImag[k]==im, "Result 5.2");
  }
  // Force throw
  ComplexWaveform s4(2*Length, Ts, InitialSample);
  BOOST_CHECK_THROW(s4+s1, const char*);

  double dbl = 6.98;
  s3 = s1 * dbl;
  for (int k=0; k< (int) Length; k++)
  {
    BOOST_CHECK_MESSAGE(s1.mReal[k] * dbl == s3.mReal[k], "Result 5.2");
    BOOST_CHECK_MESSAGE(s1.mImag[k] * dbl == s3.mImag[k], "Result 5.2");
  }
}


BOOST_AUTO_TEST_CASE (test_ComplexWaveform_division_operator)
{
  unsigned int Length = 10;
  double Ts = 20;
  long int InitSample = 22;
  double a[Length], b[Length], c[Length], d[Length];
  for (int k=0; k< (int) Length; k++)
  {
    a[k] = k;
    b[k] = 2*k;
    c[k] = k-5;
    d[k] = 2*k-10;
  }
  ComplexWaveform s1(Length, a, b, Ts, InitSample);
  ComplexWaveform s2(Length, c, d, Ts, InitSample);
  // Division
  ComplexWaveform s3 = s1/s2;
  for (int k=0; k< (int) Length; k++)
  {
    double num_re = a[k]*c[k] + b[k]*d[k];
    double num_im = b[k]*c[k] - a[k]*d[k];
    double den = c[k]*c[k] + d[k]*d[k];
    if (den!=0)
    {
      double re = num_re/den;
      double im = num_im/den;
      BOOST_CHECK_MESSAGE(s3.mReal[k] == re, "Result 6.1");
      BOOST_CHECK_MESSAGE(s3.mImag[k] == im, "Result 6.2");
    }
    else
    {
      BOOST_CHECK_MESSAGE( !isnormal(s3.mReal[k]), "Result 6.3");
      BOOST_CHECK_MESSAGE( !isnormal(s3.mImag[k]), "Result 6.4");
    }
  }
  
  double dbl = 7.42;
  s3 = s1/dbl;
  for (int k=0; k<(int) Length; k++)
  {
    double inv = 1/dbl;
    double re = a[k] * inv;
    double im = b[k] * inv;
    BOOST_CHECK_MESSAGE(s3.mReal[k] == re, "Result 6.5");
    BOOST_CHECK_MESSAGE(s3.mImag[k] == im, "Result 6.6");
  }
}

BOOST_AUTO_TEST_CASE (test_ComplexWaveform_conjugate)
{
  unsigned int Length=16;
  double Ts = 4e6;
  long int InitSample = 1800;
  double a[Length], b[Length], c[Length];
  for (int k=0; k< (int)Length; k++)
  {
    a[k] = 2*k;
    b[k] =  k;
    c[k] = -k;
  }
  ComplexWaveform s1(Length, a, b, Ts, InitSample);
  ComplexWaveform s2(Length, a, c, Ts, InitSample);
  // Conjugate
  ComplexWaveform s3 = s1.Conj();
  BOOST_CHECK_MESSAGE(s3 == s2, "Result 7.1");
}

BOOST_AUTO_TEST_SUITE_END()






  // six ways to detect and report the same error
  /*
  BOOST_CHECK( add(2,2)==4);            // continues on error

  BOOST_REQUIRE( add(2,2)==4);          // throws on error

  if (add(2,2) !=4)
    BOOST_ERROR("Here is an error..."); // continues on error
  
  if (add(2,2) !=4)
    BOOST_FAIL("Here is an error...");  // throws on error

  BOOST_CHECK_MESSAGE( add(2,2) ==4,    // continues on error  
                       "add(...) result: " << add(2,2)); 

  BOOST_CHECK_EQUAL( add(2,2), 4);      // continues on error
    */


