#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <iomanip>
#include <string> 
#include "SurfaceGeneration.hpp"



using namespace std;

namespace plt = matplotlibcpp;

unsigned int seed;

inline int myRandom() {
    seed = (214013 * seed + 2531011);
    return  (seed >> 13);
}








int main(int argc, char *argv[]){


if ((argc>1) && strcmp(argv[1], "-h")==0) {
  std::cout<<std::endl;
  std::cout<<"Input values: filename Wind theta Omega Rand_mod SIZE_X SIZE_Y Nx Ny"<<std::endl;
  std::cout<<"Example: ./main name1 100.1 0.01 0.84 5 200 200 128 128"<<std::endl;
return 0;
}

char title[256];
strcpy(title, "SIZE:Nx x Ny");

if(argc>2){
  strcpy(title+strlen(title), " Wind:");
  strcat(title, argv[2]);
  }  
if(argc>4){
  strcpy(title+strlen(title), " Omega:");
  strcat(title, argv[4]);
  } 
if(argc>5){
  strcpy(title+strlen(title), " Rand:");
  strcat(title, argv[5]);
  }  
if(argc>6){
  strcpy(title+strlen(title), " WIDTH X/Y:");
  strcat(title, argv[6]);
  } 
if(argc>8){
  strcpy(title+strlen(title), " N X/Y:");
  strcat(title, argv[8]);
  } 

std::cout << title<<std::endl;


char name[15] = {"Output"};
if (argc>1) strcpy(name,argv[1]);
double wind  = (argc>2)?strtod(argv[2], NULL): 100.0;
double theta = (argc>3)?strtod(argv[3], NULL): 0.0;
double Omega = (argc>4)?strtod(argv[4], NULL): 0.85;
double Rand_mod= (argc>5)?strtod(argv[5], NULL): 1.0;
int Dx = (argc>6)?atoi(argv[6]): 200;
int Dy = (argc>7)?atoi(argv[7]): 200;
int Nx = (argc>8)?atoi(argv[8]): 128;
int Ny = (argc>9)?atoi(argv[9]): 128;

double deltakx = (double)2.0*M_PI/Dx;
double deltaky = (double)2.0*M_PI/Dy;
double deltaX = (double)Dx/Nx;
double deltaY = (double)Dy/Ny;

std::cout << std::setprecision(15) << std::fixed;
 
double *PSI = (double*)malloc(Nx*Ny*sizeof(double));  

double *PSI_COMPLEX  = (double*)malloc(2*Nx*Ny*sizeof(double));

elfouhaily_spectrum(wind, theta, Omega, PSI,  Nx, Ny, deltakx, deltaky);


// Set seed to the random number generator
srandom(Rand_mod); 


Generate_Random_fase(PSI_COMPLEX, PSI,  Nx,  Ny);

matrix_transform(PSI_COMPLEX, Nx, Ny*2);

matrix_2d_FFT(PSI_COMPLEX, Nx, Ny);


check_MSS_spectrum(PSI,  Nx,  Ny, deltakx, deltaky);


Complex_to_Real_Matrix(PSI_COMPLEX, PSI, Nx, Ny);


check_MSS_Derivative(PSI_COMPLEX, Nx, Ny, deltaX, deltaY);

check_MSS_surface(PSI_COMPLEX, Nx, Ny);




//std::cout << "Plotting to file..." << std::endl;
//plot_complex(Nx,Ny, PSI_COMPLEX, name, title, wind, Rand_mod, deltakx, deltaky );






 
free(PSI_COMPLEX);
free(PSI);


return  0;
}






