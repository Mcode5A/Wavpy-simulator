#define BOOST_TEST_MODULE test_ComplexSignal
#include <boost/test/unit_test.hpp>

#define TEST_COMPLEXSIGNAL

//#include <math.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <ComplexSignal.hpp>


BOOST_AUTO_TEST_SUITE (test_ComplexSignal)

BOOST_AUTO_TEST_CASE (test_ComplexSignal_constructor)
{
  unsigned int len = 45;
  double       Fs = 40e6; 
  long int     init_sample = 345;

  // Constructor without data
  {
    {
      // Constructor without input data, one input parameter
      ComplexSignal* cs = new ComplexSignal(len);
      BOOST_CHECK_MESSAGE(cs->GetLength()        == len,  "Result 1.1");
      BOOST_CHECK_MESSAGE(cs->GetFs()            == 80e6, "Result 1.2");
      BOOST_CHECK_MESSAGE(cs->GetInitialSample() == 0,    "Result 1.3");
      delete cs;
    }
    {
      // Constructor without input data, two input parameters
      ComplexSignal* cs = new ComplexSignal(len, Fs);
      BOOST_CHECK_MESSAGE(cs->GetLength()        == len, "Result 1.4");
      BOOST_CHECK_MESSAGE(cs->GetFs()            == Fs,  "Result 1.5");
      BOOST_CHECK_MESSAGE(cs->GetInitialSample() == 0,   "Result 1.6");
      delete cs;
    }
    {
      // Constructor without input data, three input parameters
      ComplexSignal* cs = new ComplexSignal(len, Fs, init_sample);
      BOOST_CHECK_MESSAGE(cs->GetLength()        == len,         "Result 1.7");
      BOOST_CHECK_MESSAGE(cs->GetFs()            == Fs,          "Result 1.8");
      BOOST_CHECK_MESSAGE(cs->GetInitialSample() == init_sample, "Result 1.9");
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
    ComplexSignal* cs = new ComplexSignal(len, re, im, Fs, init_sample);
    BOOST_CHECK_MESSAGE(cs->GetLength()        == len,         "Result 1.10");
    BOOST_CHECK_MESSAGE(cs->GetFs()            == Fs,          "Result 1.11");
    BOOST_CHECK_MESSAGE(cs->GetInitialSample() == init_sample, "Result 1.12");

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
    ComplexSignal cs1(len, re, im, Fs, init_sample);
    // Check self-comaparison
    BOOST_CHECK_MESSAGE(cs1==cs1, "Result 1.15");
    // Check when equal
    ComplexSignal cs2(cs1);
    BOOST_CHECK_MESSAGE(cs2==cs1, "Result 1.16");
    // Check when different Fs
    ComplexSignal cs3(cs1);
    cs3.SetFs(Fs*2);
    BOOST_CHECK_MESSAGE(!(cs3==cs1), "Result 1.17");
    // Check when different InitSample
    ComplexSignal cs4(cs1);
    cs4.SetInitialSample(init_sample*2);
    BOOST_CHECK_MESSAGE(!(cs4==cs1), "Result 1.18");
    // Check when different real data
    ComplexSignal cs5(cs1);
    cs5.mReal[cs5.mLength-1] = cs5.mLength-1;
    BOOST_CHECK_MESSAGE(!(cs5==cs1), "Result 1.19");
    // Check when different imag data
    ComplexSignal cs6(cs1);
    cs6.mImag[cs6.mLength-1] = 0;
    BOOST_CHECK_MESSAGE(!(cs6==cs1), "Result 1.20");
    // Check when different length
    ComplexSignal cs7(cs1);
    cs7.SetLength(len-1);
    BOOST_CHECK_MESSAGE(!(cs7==cs1), "Result 1.21");
    // Check !=
    BOOST_CHECK_MESSAGE(cs3!=cs1, "Result 1.22");
    free(re);
    free(im);
  }
}


BOOST_AUTO_TEST_CASE (test_ComplexSignal_assignment_operator)
{
  unsigned int Length = 12;
  double       Fs = 10e6;
  long int     InitialSample = 0;
  // Assignment operator
  ComplexSignal s1(Length, Fs, InitialSample);
  ComplexSignal s2(2*Length, 2*Fs, InitialSample+1);
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

BOOST_AUTO_TEST_CASE (test_ComplexSignal_addition_operator)
{
  unsigned int Length = 128;
  double Fs = 10e6;
  long int InitialSample = 445;
  double re[Length], im[Length];
  for(int k=0;k<(int)Length;k++)
  {
    re[k] =  k;
    im[k] = -k;
  }
  ComplexSignal s1(Length, re, im, Fs, InitialSample);
  ComplexSignal s2(Length, Fs, InitialSample);
  // Addition
  ComplexSignal s3 = s1+s2;
  //ComplexSignal s3 = (s1).operator+(s2);
  BOOST_CHECK_MESSAGE(s3==s1, "Result 3.1");
  // Force throw
  ComplexSignal s4(2*Length, Fs, InitialSample);
  BOOST_CHECK_THROW(s4+s1, const char*);

  gsl_complex cmplx = gsl_complex_rect(2.0, 3.0);
  s3 = s1+cmplx;
  for (int k=0; k<(int)Length; k++)
  {
    double add_re = re[k] + GSL_REAL(cmplx);
    double add_im = im[k] + GSL_IMAG(cmplx);
    BOOST_CHECK_MESSAGE( add_re == s3.mReal[k], "Result 3.2");
    BOOST_CHECK_MESSAGE( add_im == s3.mImag[k], "Result 3.3");
  }
}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_subtraction_operator)
{
  unsigned int Length = 128;
  double Fs = 10e6;
  long int InitialSample = 445;
  double a[Length], b[Length];
  for(int k=0;k< (int)Length;k++)
  {
    a[k] =  k;
    b[k] = -k;
  }
  ComplexSignal s1(Length, a, b, Fs, InitialSample);
  ComplexSignal s2(Length, a, b, Fs, InitialSample);
  ComplexSignal s3(Length, Fs, InitialSample);
  // Subtraction
  ComplexSignal s4 = s1-s2;
  BOOST_CHECK_MESSAGE(s4==s3, "Result 4.1");
  // Force throw
  ComplexSignal s5(2*Length, Fs, InitialSample);
  BOOST_CHECK_THROW(s5+s1, const char*);

  gsl_complex cmplx = gsl_complex_rect(2.0, 3.0);
  s3 = s1-cmplx;
  for (int k=0; k<(int)Length; k++)
  {
    double sub_re = a[k] - GSL_REAL(cmplx);
    double sub_im = b[k] - GSL_IMAG(cmplx);
    BOOST_CHECK_MESSAGE( sub_re == s3.mReal[k], "Result 4.2");
    BOOST_CHECK_MESSAGE( sub_im == s3.mImag[k], "Result 4.2");

  }
}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_multiplication_operator)
{
  unsigned int Length = 12;
  double Fs = 10e6;
  long int InitialSample = 445;
  double a[Length], b[Length], c[Length], d[Length];
  for(int k=0;k<(int)Length;k++)
  {
    a[k] =  k;
    b[k] =  2*k;
    c[k] =  k-5;
    d[k] = 2*k-5;
  }
  ComplexSignal s1(Length, a, b, Fs, InitialSample);
  ComplexSignal s2(Length, c, d, Fs, InitialSample);
  // Subtraction
  ComplexSignal s3 = s1*s2;
  for (unsigned int k=0;k<Length;k++)
  {
    double re = a[k]*c[k] - b[k]*d[k];
    double im = b[k]*c[k] + a[k]*d[k];
    BOOST_CHECK_MESSAGE(s3.mReal[k]==re, "Result 5.1");
    BOOST_CHECK_MESSAGE(s3.mImag[k]==im, "Result 5.2");
  }
  // Force throw
  ComplexSignal s4(2*Length, Fs, InitialSample);
  BOOST_CHECK_THROW(s4+s1, const char*);

  double dbl = 2.0;
  s3 = s1*dbl;
  for (int k=0; k< (int) Length; k++)
  {
    double re = a[k] * dbl;
    double im = b[k] * dbl;
    BOOST_CHECK_MESSAGE(s3.mReal[k] == re, "Result 5.3");
    BOOST_CHECK_MESSAGE(s3.mImag[k] == im, "Result 5.4");
  }
}


BOOST_AUTO_TEST_CASE (test_ComplexSignal_division_operator)
{
  unsigned int Length = 10;
  double Fs = 20;
  long int InitSample = 22;
  double a[Length], b[Length], c[Length], d[Length];
  for (int k=0; k< (int) Length; k++)
  {
    a[k] = k;
    b[k] = 2*k;
    c[k] = k-5;
    d[k] = 2*k-10;
  }
  ComplexSignal s1(Length, a, b, Fs, InitSample);
  ComplexSignal s2(Length, c, d, Fs, InitSample);
  // Division signal/signal
  ComplexSignal s3 = s1/s2;
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
  double dbl = 2.0;
  s3 = s1/dbl;
  for (int k=0; k< (int) Length; k++)
  {
    double re = a[k]/dbl;
    double im = b[k]/dbl;
    BOOST_CHECK_MESSAGE(s3.mReal[k] == re, "Result 6.5");
    BOOST_CHECK_MESSAGE(s3.mImag[k] == im, "Result 6.6");
  }
}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_conjugate)
{
  unsigned int Length=16;
  double Fs = 4e6;
  long int InitSample = 1800;
  double a[Length], b[Length], c[Length];
  for (int k=0; k< (int)Length; k++)
  {
    a[k] = 2*k;
    b[k] =  k;
    c[k] = -k;
  }
  ComplexSignal s1(Length, a, b, Fs, InitSample);
  ComplexSignal s2(Length, a, c, Fs, InitSample);
  // Conjugate
  ComplexSignal s3 = s1.Conj();
  BOOST_CHECK_MESSAGE(s3 == s2, "Result 7.1");
}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_phasor)
{
  double Freq = 1234567.89;
  double Fs   = 80e6;
  double phi0 = M_SQRTPI;
  unsigned int length = 80e3;

  ComplexSignal P = ComplexSignal::Phasor(Freq, Fs, length, phi0);
  unsigned int k;
  for (k=0; k<length; k++)
  {
    BOOST_CHECK_MESSAGE(P.mReal[k] == cos(2*M_PI*Freq/Fs*k + phi0), "Result 8.1");
    BOOST_CHECK_MESSAGE(P.mImag[k] == sin(2*M_PI*Freq/Fs*k + phi0), "Result 8.2");
  }
}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_CorrelationLimits)
{
  unsigned int lenx, leny, N, n0;
  int minm, maxm;
  bool b;
  lenx = 80000;
  leny = 80000;
  N = 8000;
  n0 = 10;
  minm = -10;
  maxm =  10;

  b = ComplexSignal::CorrelationLimitsCorrect(lenx, leny, N, n0, minm, maxm);
  BOOST_CHECK_MESSAGE(b, "Result 9.1");
  
  b = ComplexSignal::CorrelationLimitsCorrect(lenx, leny, 0, n0, minm, maxm);
  BOOST_CHECK_MESSAGE(!b, "Result 9.2");
  
  b = ComplexSignal::CorrelationLimitsCorrect(lenx, leny, N, -1, minm, maxm);
  BOOST_CHECK_MESSAGE(!b, "Result 9.3");

  b = ComplexSignal::CorrelationLimitsCorrect(lenx, leny, N, lenx, minm, maxm);
  BOOST_CHECK_MESSAGE(!b, "Result 9.4");
  
  b = ComplexSignal::CorrelationLimitsCorrect(lenx, leny, N, n0, leny-1, leny);
  BOOST_CHECK_MESSAGE(!b, "Result 9.5");

  b = ComplexSignal::CorrelationLimitsCorrect(lenx, leny, N, lenx-N, minm, maxm);
  BOOST_CHECK_MESSAGE(!b, "Result 9.6");

  b = ComplexSignal::CorrelationLimitsCorrect(lenx, leny, N, n0, minm, 20);
  BOOST_CHECK_MESSAGE(!b, "Result 9.7");

  b = ComplexSignal::CorrelationLimitsCorrect(lenx, 8000+maxm-1, N, 10, 0, maxm);
  BOOST_CHECK_MESSAGE(!b, "Result 9.8");

}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_CorrelationRealSingleLag)
{
  unsigned int len = 10000;
  unsigned int N = 1000;
  unsigned int n0 = 10;

  double x[len];
  double y[len];
  for (unsigned int k=0; k<len; k++)
  {
    x[k] = pow(-1,k);
    y[k] = pow(-1,k);
    //std::cout << "x[" << k << "]=" << x[k] << "  ";
    //std::cout << "y[" << k << "]=" << y[k] << "  ";
  }
  
  // lag = 0
  {
    int m = 0;
    double Rxy = 1;
    double Rxy_serial = ComplexSignal::CorrelationRealSingleLag(x, y, N, n0, m, false);
    double Rxy_parallel = ComplexSignal::CorrelationRealSingleLag(x, y, N, n0, m, true);
    BOOST_CHECK_MESSAGE( Rxy_serial   == Rxy, "Result 10.1");
    BOOST_CHECK_MESSAGE( Rxy_parallel == Rxy, "Result 10.2");
  }
  // lag = -1
  {
    int m = -1;
    double Rxy = -1;
    double Rxy_serial = ComplexSignal::CorrelationRealSingleLag(x, y, N, n0, m, false);
    double Rxy_parallel = ComplexSignal::CorrelationRealSingleLag(x, y, N, n0, m, true);
    BOOST_CHECK_MESSAGE( Rxy_serial   == Rxy, "Result 10.3");
    BOOST_CHECK_MESSAGE( Rxy_parallel == Rxy, "Result 10.4");
  }
  // lag = 1
  {
    int m = 1;
    double Rxy = -1;
    double Rxy_serial = ComplexSignal::CorrelationRealSingleLag(x, y, N, n0, m, false);
    double Rxy_parallel = ComplexSignal::CorrelationRealSingleLag(x, y, N, n0, m, true);
    BOOST_CHECK_MESSAGE( Rxy_serial   == Rxy, "Result 10.5");
    BOOST_CHECK_MESSAGE( Rxy_parallel == Rxy, "Result 10.6");
  }
}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_CorrelationComplexSingleLag)
{
  unsigned int len = 4000;
  unsigned int N = 400;
  unsigned int n0 = 40;

  double re[len];
  double im[len];
  
  for (unsigned int k=0; k<len; k++)
  {
    if (k%4 == 0)
    {
      re[k] = 1;
      im[k] = 0;
    }
    else if (k%4 == 1)
    {
      re[k] = 0;
      im[k] = 1;
    }
    else if (k%4 == 2)
    {
      re[k] = -1;
      im[k] =  0;
    }
    else if (k%4 == 3)
    {
      re[k] = 0;
      im[k] = -1;
    }
  }
  { // lag 0;
    int m = 0;
    bool b;
    gsl_complex Rxy = gsl_complex_rect(1,0);
    ComplexSignal x(len, re, im);
    gsl_complex Rxy_serial = ComplexSignal::CorrelationComplexSingleLag(x, x, N, n0, m, false);
    gsl_complex Rxy_parallel = ComplexSignal::CorrelationComplexSingleLag(x, x, N, n0, m, true);
    b = ((GSL_REAL(Rxy) == GSL_REAL(Rxy_serial)) && (GSL_IMAG(Rxy) == GSL_IMAG(Rxy_serial)));
    BOOST_CHECK_MESSAGE(b, "Result 11.1");
    b = ((GSL_REAL(Rxy) == GSL_REAL(Rxy_parallel)) && (GSL_IMAG(Rxy) == GSL_IMAG(Rxy_parallel)));
    BOOST_CHECK_MESSAGE(b, "Result 11.2");
  }
  { // lag 1;
    int m = 1;
    bool b;
    gsl_complex Rxy = gsl_complex_rect(0,1);
    ComplexSignal x(len, re, im);
    gsl_complex Rxy_serial = ComplexSignal::CorrelationComplexSingleLag(x, x, N, n0, m, false);
    gsl_complex Rxy_parallel = ComplexSignal::CorrelationComplexSingleLag(x, x, N, n0, m, true);
    b = ((GSL_REAL(Rxy) == GSL_REAL(Rxy_serial)) && (GSL_IMAG(Rxy) == GSL_IMAG(Rxy_serial)));
    BOOST_CHECK_MESSAGE(b, "Result 11.3");
    b = ((GSL_REAL(Rxy) == GSL_REAL(Rxy_parallel)) && (GSL_IMAG(Rxy) == GSL_IMAG(Rxy_parallel)));
    BOOST_CHECK_MESSAGE(b, "Result 11.4");
  }
  { // lag 2;
    int m = 2;
    bool b;
    gsl_complex Rxy = gsl_complex_rect(-1,0);
    ComplexSignal x(len, re, im);
    gsl_complex Rxy_serial = ComplexSignal::CorrelationComplexSingleLag(x, x, N, n0, m, false);
    gsl_complex Rxy_parallel = ComplexSignal::CorrelationComplexSingleLag(x, x, N, n0, m, true);
    b = ((GSL_REAL(Rxy) == GSL_REAL(Rxy_serial)) && (GSL_IMAG(Rxy) == GSL_IMAG(Rxy_serial)));
    BOOST_CHECK_MESSAGE(b, "Result 11.5");
    b = ((GSL_REAL(Rxy) == GSL_REAL(Rxy_parallel)) && (GSL_IMAG(Rxy) == GSL_IMAG(Rxy_parallel)));
    BOOST_CHECK_MESSAGE(b, "Result 11.6");
  }
  { // lag 3;
    int m = 3;
    bool b;
    gsl_complex Rxy = gsl_complex_rect(0,-1);
    ComplexSignal x(len, re, im);
    gsl_complex Rxy_serial = ComplexSignal::CorrelationComplexSingleLag(x, x, N, n0, m, false);
    gsl_complex Rxy_parallel = ComplexSignal::CorrelationComplexSingleLag(x, x, N, n0, m, true);
    b = ((GSL_REAL(Rxy) == GSL_REAL(Rxy_serial)) && (GSL_IMAG(Rxy) == GSL_IMAG(Rxy_serial)));
    BOOST_CHECK_MESSAGE(b, "Result 11.7");
    b = ((GSL_REAL(Rxy) == GSL_REAL(Rxy_parallel)) && (GSL_IMAG(Rxy) == GSL_IMAG(Rxy_parallel)));
    BOOST_CHECK_MESSAGE(b, "Result 11.8");
  }
}



BOOST_AUTO_TEST_CASE (test_ComplexSignal_Correlation)
{
  unsigned int len = 4000;
  unsigned int N = 400;
  unsigned int n0 = 40;

  double re[len];
  double im[len];
  
  for (unsigned int k=0; k<len; k++)
  {
    if (k%4 == 0)
    {
      re[k] = 1;
      im[k] = 0;
    }
    else if (k%4 == 1)
    {
      re[k] = 0;
      im[k] = 1;
    }
    else if (k%4 == 2)
    {
      re[k] = -1;
      im[k] =  0;
    }
    else if (k%4 == 3)
    {
      re[k] = 0;
      im[k] = -1;
    }
  }

  { 
    int min_m = 0;
    int step_m = 1;
    int max_m = 3;
    ComplexSignal x(len, re, im);
    double real[max_m+1] = {1, 0, -1, 0};
    double imag[max_m+1] = {0, 1, 0, -1};
    ComplexWaveform Rxy(max_m+1, real, imag);
    ComplexWaveform CW_serial = ComplexSignal::Correlation(x, x, N, n0, min_m, step_m, max_m, false);
    ComplexWaveform CW_parallel = ComplexSignal::Correlation(x, x, N, n0, min_m, step_m, max_m, true);
    BOOST_CHECK_MESSAGE(Rxy == CW_serial, "Result 12.1");
    BOOST_CHECK_MESSAGE(Rxy == CW_parallel, "Result 12.2");
  }
}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_AWGN)
{
  double Fs   = 80e6;
  unsigned int length = 80e3;
  double  sigma = 2;
  int seed = 12345;

  ComplexSignal P = ComplexSignal::AWGN(sigma, length, Fs, seed);
  unsigned int k;
  double sigma2_est = 0;
  double aux;
  // Compute std. deviation
  for (k=0; k<length; k++)
  {
    aux = P.mReal[k] * P.mReal[k] + P.mImag[k] * P.mImag[k];
    sigma2_est = sigma2_est + aux;
    //BOOST_CHECK_MESSAGE(P.mReal[k] == cos(2*M_PI*Freq/Fs*k + phi0), "Result 8.1");
    //BOOST_CHECK_MESSAGE(P.mImag[k] == sin(2*M_PI*Freq/Fs*k + phi0), "Result 8.2");
  }
  sigma2_est = sigma2_est/(length-1);
  double sigma_est = sqrt(sigma2_est)/sqrt(2);
  double d = abs(sigma - sigma_est);
  BOOST_CHECK_MESSAGE( d < 0.01 , "Result 13.1"); // Check that the sigma is the expected one
}

BOOST_AUTO_TEST_CASE (test_ComplexSignal_Crop)
{
  double Fs   = 80e6;
  unsigned int length = 100;

  double real[length];
  double imag[length];

  for (unsigned int k=0;k<length;k++)
  {
    real[k] = k;
    imag[k] = 2*k;
  }

  ComplexSignal cs = ComplexSignal(length, real, imag, Fs, 0);

  unsigned int n0 = 30;
  unsigned int n1 = 40;

  // Test that input parameters checkings are OK
  BOOST_CHECK_THROW(ComplexSignal cr1 = cs.Crop(n1, n0), const char*);
  BOOST_CHECK_THROW(ComplexSignal cr2 = cs.Crop(length, n1), const char*);
  BOOST_CHECK_THROW(ComplexSignal cr3 = cs.Crop(n0, length), const char*); 
  
  // Check that the signal is properly cropped
  ComplexSignal cr = cs.Crop(n0, n1);
  BOOST_CHECK_MESSAGE( cr.GetLength() == n1-n0+1, "Result 14.1");
  BOOST_CHECK_MESSAGE( cr.GetFs() == cs.GetFs(), "Result 14.2");
  BOOST_CHECK_MESSAGE( cr.GetInitialSample() == cs.GetInitialSample() + (long int) n0, "Result 14.3");

  for (unsigned int k=0; k<n1-n0+1; k++)
  {
    std::cout << "k=" << k << "   cr.mReal[k]=" << cr.mReal[k] << std::endl;
    std::cout << "      "  <<    "cr.mImag[k]=" << cr.mImag[k] << std::endl;
    BOOST_CHECK_MESSAGE( cr.mReal[k] == n0 + k,      "Result 14.4");
    BOOST_CHECK_MESSAGE( cr.mImag[k] == 2*(n0 + k), "Result 14.5");
  } 
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


