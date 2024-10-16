#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
//#define pi 3.141592  -> use the definition of pi in math.h -> M_PI
#include "../matplotlib-cpp-master/matplotlibcpp.h"
#include <cmath>
#include <vector>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <iomanip>
#include <string> 
using namespace std;

namespace plt = matplotlibcpp;

void check_MSS_Derivative(double *PSI_COMPLEX, int Nx, int Ny, double deltaX, double deltaY){


double sum_derv=0.0;

for(int i=0; i<Nx-1; i++)
  for(int j=0; j<Ny-1; j++){
  
    double ax = (double)i*deltaX-((double)i+1.0)*deltaX;
    double ay = (double)j*deltaY-((double)j)*deltaY;
    double az = PSI_COMPLEX[i*2*Ny+j*2]-PSI_COMPLEX[(i+1)*2*Ny+j*2];
    
    double bx = (double)i*deltaX-((double)i)*deltaX;
    double by = (double)j*deltaY-((double)j+1.0)*deltaY;
    double bz = PSI_COMPLEX[i*2*Ny+j*2]-PSI_COMPLEX[i*2*Ny+(j+1)*2];
    
    
    
    double cx = ay*bz-az*by;
    double cy = -(ax*bz-az*bx);
    double cz = ax*by-ay*bx;
    
    
    double c_mod = sqrt(cx*cx+cy*cy+cz*cz);
    
    cx/=c_mod;
    cy/=c_mod;
    cz/=c_mod;


    double beta = acos(cz);    // scalar product between C and 0,0,1
    
    double pendent = tan(beta);
    

    if(!(pendent != pendent))
      sum_derv+=(pendent*pendent);
     
    
  }

sum_derv/=((Nx-1)*(Ny-1));

std::cout << "derivada: "<<sum_derv<<std::endl;
}

void check_MSS_spectrum(double * PSI, int Nx, int Ny,double deltakx,double deltaky){

double MSSx=0.0;
double MSSy=0.0;
for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++){
  
    double kx_vec = (double)(-Nx/2+i)*deltakx;
    double ky_vec = (double)(-Ny/2+j)*deltaky;
  
    
    MSSx+=kx_vec*kx_vec*PSI[i*Ny+j]*deltakx*deltaky;
    MSSy+=ky_vec*ky_vec*PSI[i*Ny+j]*deltakx*deltaky;
    
  }
std::cout << "MSS: "<<MSSx+MSSy<<std::endl;
std::cout << "MSSx: "<<MSSx<<std::endl;
std::cout << "MSSy: "<<MSSy<<std::endl;
}

void check_MSS_surface(double *PSI_COMPLEX, int Nx, int Ny){
double checksum_complex=0.0;
 
for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++){
    checksum_complex+=PSI_COMPLEX[i*2*Ny+j*2]*PSI_COMPLEX[i*2*Ny+j*2];
  }  
  
  
std::cout << "Surface: "<<checksum_complex/(Nx*Ny)<<std::endl;
}







void plot(int Nx,int Ny, double *PSI){
 
  std::vector<std::vector<double>> x, y, z;
  for (int i = 0; i < Nx;  i++) {
      std::vector<double> x_row, y_row, z_row;
      for (int j = 0; j < Ny; j++) {
          x_row.push_back(i);
          y_row.push_back(j);
          z_row.push_back(PSI[i*Ny+j]);
        
      }
      
      x.push_back(x_row);
      y.push_back(y_row);
      z.push_back(z_row);
  }
  
  plt::plot_surface(x, y, z);
  //plt::show();
  plt::title("PSI ");
   const char* filename = "./k.png";
  std::cout << "Saving result to " << filename << std::endl;;
  plt::save(filename);



}
unsigned int seed;

inline int myRandom() {
    seed = (214013 * seed + 2531011);
    return  (seed >> 13);
}


int  elfouhaily_spectrum(double U10,double theta,double Omega,  double *PSI,  int Nx,int Ny, double deltakx, double deltaky)
{
    
// [S, PSI] = elfouhaily_spectrum(U10, theta, Omega)
//
//  Output:
//     S     : Omnidirectional spectrum
//     PSI   : Directional spectrum
//
//  Input:
//     U10   : windspeed at 10m above surface [m/s]
//     theta : angle between wind and waves [deg]
//     Omega : inverse dimensionless wave age > 0.83
//     kx    : wavenumber in x direction [rad/m]
//     ky    : wavenumber in y direction [rad/m]
//

// ----- Checking validity range of input and derived parameters

if(Omega < 0.84){
  std::cout<<"Omega must be larger than 0.84, but Omega="<< Omega<<std::endl;
  return 1;
}
double Omegac = Omega*cos(theta*M_PI/180);

if (Omegac<0.84 || Omegac>5){
  std::cout<<"Omega_c must be larger than 0.83 and smaller than 5, but Omega_c="<< Omegac<<std::endl;
  return 1;
}
double k_vec;
double phi_vec;
double c_vec; 
double Gmma_vec;  
double Jp_vec;  
double Lpm_vec;  
double Fp_vec;  
double Bl_vec;  
double Fm_vec;  
double Bh_vec;  
double Delta_vec;  
double g = 9.81; // [m/s^2] gravity acceleration
double cp = U10/Omega;    // Phase velocity at the spectral peak
double kp = g/(cp*cp);     // Wave number of the spectral peak (assumtion of deep waters. TODO ref. needed)

double alphap = 0.006*sqrt(Omega);
double gmma = ((0.84<=Omegac) && (Omegac<=1.0))? 1.7 : 1.7 + 6.0*log10(Omegac);

double sigma = 0.08*(1.0+4.0*pow(Omegac,(-3)));
double cm = 0.23; // [m/s] Elfouhaily 1997.
double ufric = sqrt(1e-3 * (0.8 + 0.065*U10))*U10;
double alpham = (ufric <= cm)? 1.0+log(ufric/cm) : 1.0 + 2.0*log(ufric/cm);
alpham*=0.01;
double km =  (2.0*g/(cm*cm)); // assuming deep waters => =370 rad/m
double a0 = 0.25*log(2.0); // Natural logarithm
double ap = 4.0;
double am = 0.13*ufric/cm;


double S_vec;
double PHI_vec;



double kx_vec;
double ky_vec;
 
double S=0.0;


for(int i=0; i<Nx; i++){
  for(int j=0; j<Ny; j++){

    kx_vec = (double)(-Nx/2+i)*deltakx;
    ky_vec = (double)(-Ny/2+j)*deltaky;

    k_vec = sqrt(kx_vec*kx_vec + ky_vec*ky_vec);
    
    phi_vec = atan2(ky_vec, kx_vec); 

    c_vec = sqrt(g/k_vec *(1.0+pow(k_vec/km, 2)) );    

    Gmma_vec = exp(-(pow(sqrt(k_vec/kp)-1,2))/(2*sigma*sigma));  // -> matrix

    Jp_vec = pow(gmma, Gmma_vec);    // -> matrix

    Lpm_vec = exp(-1.25*pow(kp/k_vec,2));// -> matrix
      
    Fp_vec = exp(-(Omega/sqrt(10.0))*(sqrt(k_vec/kp)-1));// -> matrix
    
    Bl_vec = 0.5*alphap*(cp/c_vec)*Fp_vec;// -> matrix

    Fm_vec = exp(-0.25*(pow((k_vec/km) -1 ,2)));  // -> matrix

 
    Bh_vec = 0.5 * alpham *(cm/c_vec)*Fm_vec; // -> matrix

    S_vec = pow(k_vec,(-3))*(Bl_vec + Bh_vec)*Lpm_vec*Jp_vec;  // -> matrix

    Delta_vec = tanh(a0 + ap*pow(c_vec/cp,2.5) + am*pow(cm/c_vec,2.5)); // -> matrix


  // ----- Computing the directional spectrum
    PHI_vec = 1.0/(2.0*M_PI) * (1.0 + Delta_vec * cos(2.0*phi_vec));  // -> matrix

    PSI[i*Ny+j] = S_vec*PHI_vec/k_vec;  //division by 0 on i=Nx/2 j =Ny/2


}}

  PSI[Nx*Ny/2+Ny/2]=0.0;
  
  return 0;


}


 
 
void plot_complex(int Nx, int Ny, double *PSI, char *name,char *title, double wind, double Rand_mod, double deltakx, double deltaky){

  
  std::vector<std::vector<double>> x, y, z;
  for (int i = 0; i < Nx;  i++) {
      std::vector<double> x_row, y_row, z_row;
      for (int j = 0; j < Ny; j++) {
          x_row.push_back(2.0*M_PI*i/(deltakx*Nx));
          y_row.push_back(2.0*M_PI*j/(deltaky*Ny));
          z_row.push_back(PSI[i*2*Ny+2*j]);
          
      }
      
      x.push_back(x_row);
      y.push_back(y_row);
      z.push_back(z_row);
  }
  
  plt::plot_surface(x, y, z);
  plt::title(title);
  char* filename = (char *) malloc(20*sizeof(char));
  strcat(filename, "./"); 
  strcpy(filename, name);
  strcat(filename, ".png"); 
  std::cout << "Saving result to " << filename << std::endl;;
  plt::save(filename);


}


void Generate_Random_fase(double *PSI_COMPLEX, double *PSI, int Nx, int Ny){

  // MAKING FIRST ROW IMAG PART =0// ALSO NOT MULTIPILY BY COS BC THERE IS NO IMAG PART// ALSO THERE IS NO SIMETRIC PART
  for(int i=0; i < Ny; i++){   //fila 0
     
      PSI_COMPLEX[2*i] = sqrt(PSI[i]);
      PSI_COMPLEX[2*i+1]=0.0;
      }
  
  // MAKING FIRST COLUMN IMAG PART =0// ALSO NOT MULTIPILY BY COS BC THERE IS NO IMAG PART// ALSO THERE IS NO SIMETRIC PART
  for(int i=0; i < Nx; i++){   //columna 0
     
      PSI_COMPLEX[2*i*Ny] = sqrt(PSI[i*Ny]); 
      PSI_COMPLEX[2*i*Ny+1]=0.0;
      }
  
  for(int i=1; i < Ny/2; i++){   //MAKING MIDDLE ROW SIMETRIC
  
      // The phase is a double. Otherwise it will be trucated to integers, and the phase will loose randomness 
      // The phase has to be a random number between 0 and 2*pi.
      // We use standard random number generator functions. If needed, we can optimize this later
  
      double randomphase = 2.0*M_PI*(random()/(double)RAND_MAX);
      PSI_COMPLEX[Nx*Ny+2*i] = sqrt(PSI[Ny*Nx/2 +i])*cos(randomphase);
      PSI_COMPLEX[Nx*Ny+2*i+1] = sqrt(PSI[Ny*Nx/2 +i])*sin(randomphase);
      
       //simetric values as above
      PSI_COMPLEX[Nx*Ny+2*Ny-2*i] = PSI_COMPLEX[Nx*Ny+2*i];
      PSI_COMPLEX[Nx*Ny+2*Ny-2*i+1] = -PSI_COMPLEX[Nx*Ny+2*i+1];
    
    }
  
  PSI_COMPLEX[Nx*Ny+Ny]=sqrt(PSI[Nx*Ny/2+Ny/2]); //MAKING MIDDLE POINT REAL PART = MODULE
  PSI_COMPLEX[Nx*Ny+Ny+1]=0.0;  //MAKING IMAG PART OF MIDDLE POINT =0
   
     
  for(int i = 1; i < Nx/2; i++){ //MAKING THE REST OF THE SURFACE SIMETRIC
    for (int j = 1; j < Ny; j++) {
           
  
            double randomphase = 2.0*M_PI*(random()/(double)RAND_MAX);
            PSI_COMPLEX[2*i*Ny+2*j]=sqrt(PSI[i*Ny+j])*cos(randomphase);
            PSI_COMPLEX[2*i*Ny+2*j+1]=sqrt(PSI[i*Ny+j])*sin(randomphase); 
            
            //simetric values as above
            PSI_COMPLEX[2*Nx*Ny-2*(i-1)*Ny-2*j]=PSI_COMPLEX[2*i*Ny+2*j];
            PSI_COMPLEX[2*Nx*Ny-2*(i-1)*Ny-2*j+1]=-PSI_COMPLEX[2*i*Ny+2*j+1];
  
  	}
  }
    
}

void matrix_transform(double * PSI, int Nx, int Ny){
  
  double  aux_table;
  
  for(int i=0; i < Nx/2; i++)
    for(int j=0; j < Ny/2; j++){
        aux_table = PSI[i*Ny+j];
        PSI[i*Ny+j] = PSI[(i+Nx/2)*Ny+j+Ny/2];
        PSI[(i+Nx/2)*Ny+j+Ny/2] = aux_table;
  }
  

  
  for(int i=0; i < Nx/2; i++)
    for(int j=0; j < Ny/2; j++){
        aux_table = PSI[i*Ny+j+Ny/2];
        PSI[i*Ny+j+Ny/2] = PSI[(i+Nx/2)*Ny+j];
        PSI[(i+Nx/2)*Ny+j] = aux_table;
        
  }
  
  
  
}


void matrix_2d_FFT(double *PSI_COMPLEX, int Nx, int Ny){

for (int i = 0; i < Ny; i++)
  gsl_fft_complex_radix2_backward(PSI_COMPLEX+2*i,Ny,Nx);
for (int i = 0; i < Nx; i++)
 gsl_fft_complex_radix2_backward(PSI_COMPLEX+2*i*Ny,1,Ny);
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



elfouhaily_spectrum(wind, theta, Omega, PSI,  Nx, Ny, deltakx, deltaky);


check_MSS_spectrum(PSI,  Nx,  Ny, deltakx, deltaky);



// Set seed to the random number generator
srandom(Rand_mod); 



double *PSI_COMPLEX  = (double*)malloc(2*Nx*Ny*sizeof(double));

Generate_Random_fase(PSI_COMPLEX, PSI,  Nx,  Ny);


free(PSI);


matrix_transform(PSI_COMPLEX, Nx, Ny*2);

matrix_2d_FFT(PSI_COMPLEX, Nx, Ny);


check_MSS_Derivative(PSI_COMPLEX, Nx, Ny, deltaX, deltaY);

check_MSS_surface(PSI_COMPLEX, Nx, Ny);


//std::cout << "Plotting to file..." << std::endl;
//plot_complex(Nx,Ny, PSI_COMPLEX, name, title, wind, Rand_mod, deltakx, deltaky );






 
free(PSI_COMPLEX);



return  0;
}






