#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#define pi 3.141592
#include "../../matplotlibcpp.h"
#include <cmath>
#include <vector>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <iomanip>

using namespace std;

namespace plt = matplotlibcpp;


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
    return (seed >> 13);
}
 
 
 int elfouhaily_spectrum(double U10,double theta,double Omega,double *kx,double *ky,double * S, double *PSI, double *PHI, int Nx,int Ny)
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


double Omegac = Omega*cos(theta*pi/180);

if (Omegac<0.84 || Omegac>5){
  std::cout<<"Omega_c must be larger than 0.83 and smaller than 5, but Omega_c="<< Omegac<<std::endl;
  return 1;
}



double *k=(double*)malloc(Nx*Ny*sizeof(double));
double *phi=(double*)malloc(Nx*Ny*sizeof(double));

for(int i=0; i<Nx; i++){
  for(int j=0; j<Ny; j++){

    k[i*Ny+j] = sqrt(kx[i]*kx[i] + ky[j]*ky[j]);
    phi[i*Ny+j] = atan2(ky[j], kx[i]); 
 }
}


// ----- constants
double g = 9.81; // [m/s^2] gravity acceleration




// ----- Computing long wave part of the spectrum

 
double alphap; 
double *c; // Phase velocity (deep waters assumption)
double cp = U10/Omega;    // Phase velocity at the spectral peak
double kp = g/(cp*cp);     // Wave number of the spectral peak (assumtion of deep waters. TODO ref. needed)

c  =(double*)malloc(Nx*Ny*sizeof(double));  

for(int i=0; i<Nx*Ny; i++){

    c[i] = sqrt(g/k[i]);   
}

alphap = 0.006*sqrt(Omegac);

double gmma;

if ((0.84<=Omegac) && (Omegac<=1))
  gmma = 1.7;
else //((1<Omegac) && (Omegac<=5))
  gmma = 1.7 + 6*log10(Omegac);



double *Gmma =(double*)malloc(Nx*Ny*sizeof(double));  
double *Jp =(double*)malloc(Nx*Ny*sizeof(double));  
double *Lpm =(double*)malloc(Nx*Ny*sizeof(double));  
double *Fp =(double*)malloc(Nx*Ny*sizeof(double));  
double *Bl =(double*)malloc(Nx*Ny*sizeof(double));  

double sigma = 0.08*(1.0+4.0*pow(Omegac,(-3)));


for(int i=0; i<Nx; i++){
  for(int j=0; j<Ny; j++){
      Gmma[i*Ny+j] = exp(-(pow(sqrt(kp/k[i*Ny+j])-1,2))/(2*sigma*sigma));  // -> matrix

      Jp[i*Ny+j] = pow(gmma, Gmma[i*Ny+j]);    // -> matrix

      Lpm[i*Ny+j] = exp(-5.0/4.0*pow(kp/k[i*Ny+j],2));// -> matrix
      

      Fp[i*Ny+j] = Lpm[i*Ny+j]*Jp[i*Ny+j]*exp(-(Omega/sqrt(10))*(sqrt(kp/k[i*Ny+j])-1));// -> matrix
      
 
      Bl[i*Ny+j] = 0.5*alphap*(cp/c[i*Ny+j])*Fp[i*Ny+j];// -> matrix

  }
}

free(Fp);
free(Jp);
free(Lpm);
free(Gmma);





// ----- Computing the short wave part of the spectrum
double cm = 0.23; // [m/s] Elfouhaily 1997.

double ufric = sqrt(1e-3 * (0.8 + 0.065*U10))*U10;
double alpham;
if (ufric <= cm)
  alpham = 1+log(ufric/cm);
else 
  alpham = 1 + 2*log(ufric/cm);


double km = sqrt(2*g/cm); // assuming deep waters => =370 rad/m

double *Fm =(double*)malloc(Nx*Ny*sizeof(double));  
double *Bh =(double*)malloc(Nx*Ny*sizeof(double));  

for(int i=0; i<Nx; i++){
  for(int j=0; j<Ny; j++){
    //FALTAN () PARA LA RESTA -1 ??
    Fm[i*Ny+j] = exp(-0.25*(pow((k[i*Ny+j]/km) -1 ,2)));  // -> matrix

    //FALTAN () PARA saber si * vantes de /
    Bh[i*Ny+j] = 0.5 * alpham *(cm/c[i*Ny+j])*Fm[i*Ny+j]; // -> matrix

    S[i*Ny+j] = pow(k[i*Ny+j],(-3))*(Bl[i*Ny+j] + Bh[i*Ny+j]);  // -> matrix
}}


free(Bl);
free(Bh);
free(Fm);


// ----- Computing the Spreading function

double a0 = 0.25*log(2); // Natural logarithm
double ap = 4;
double am = 0.13*ufric/cm;
double *Delta =(double*)malloc(Nx*Ny*sizeof(double));  



for(int i=0;i<Nx; i++){
  for(int j=0; j<Ny; j++){
    
    Delta[i*Ny+j] = tanh(a0 + ap*pow(c[i*Ny+j]/cp,2.5) + am*pow(cm/c[i*Ny+j],2.5)); // -> matrix


    // ----- Computing the directional spectrum
    PHI[i*Ny+j] = 1/(2*pi) * (1 + Delta[i*Ny+j] * cos(2*phi[i*Ny+j]));  // -> matrix

    PSI[i*Ny+j] = S[i*Ny+j]*PHI[i*Ny+j]/k[i*Ny+j];  //division by 0 on i=Nx/2 j =Ny/2

    
}}

free(k);
free(Delta);
free(phi);
free(c);

 
 

return 0;
}

int Powerof2(int n, int* m, int* twopm)
{
    if (n <= 1) {
        *m = 0;
        *twopm = 1;
        return(0);
    }

    *m = 1;
    *twopm = 2;
    do {
        (*m)++;
        (*twopm) *= 2;
    } while (2 * (*twopm) <= n);

    if (*twopm != n)
        return(0);
    else
        return(1);
}

int FFT(int dir, int m, double* x, double* y)
{
    long nn, i, i1, j, k, i2, l, l1, l2;
    double c1, c2, tx, ty, t1, t2, u1, u2, z;

    /* Calculate the number of points */
    nn = 1;
    for (i = 0; i < m; i++)
        nn *= 2;

    /* Do the bit reversal */
    i2 = nn >> 1;
    j = 0;
    for (i = 0; i < nn - 1; i++) {
        if (i < j) {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    /* Compute the FFT */
    c1 = -1.0;
    c2 = 0.0;
    l2 = 1;
    for (l = 0; l < m; l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j = 0; j < l1; j++) {
            for (i = j; i < nn; i += l2) {
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            z = u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }

    /* Scaling for forward transform */
    if (dir == 1) {
        for (i = 0; i < nn; i++) {
            x[i] /= (double)nn;
            y[i] /= (double)nn;
        }
    }

    return(1);
}



void plot_complex(int Nx,int Ny, double *PSI, char *name ){

  
  std::vector<std::vector<double>> x, y, z;
  for (int i = 0; i < Nx;  i++) {
      std::vector<double> x_row, y_row, z_row;
      for (int j = 0; j < Ny; j++) {
          x_row.push_back(i);
          y_row.push_back(j);
          z_row.push_back(PSI[i*2*Ny+2*j]);
          
      }
      
      x.push_back(x_row);
      y.push_back(y_row);
      z.push_back(z_row);
  }
  
  plt::plot_surface(x, y, z);
  //plt::show();
  plt::title("PSI time resolution 128, wind=10, fase Random:25");
  char* filename = (char *) malloc(20*sizeof(char));
  strcat(filename, "./"); 
  strcpy(filename, name);
  strcat(filename, ".png"); 
  std::cout << "Saving result to " << filename << std::endl;;
  plt::save(filename);


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


int main(){



int Nx,Ny;
Nx=128;
Ny=128;
 
double *kx =(double*)malloc(Nx*sizeof(double));  
double *ky =(double*)malloc(Ny*sizeof(double));  


double deltakx = (double)200.0/Nx;

double deltaky = (double)200.0/Ny;


for(int i=0; i<Nx; i++){
  kx[i]= (double)(-Nx/2+i)*deltakx;
}


 
for(int i=0; i<Ny; i++){
  ky[i]= (double)(-Ny/2+i)*deltaky;
}
 


double * S; 
double *PSI;
double *PHI;

S =(double*)malloc(Nx*Ny*sizeof(double));  
PHI =(double*)malloc(Nx*Ny*sizeof(double));   
PSI =(double*)malloc(Nx*Ny*sizeof(double));  


std::cout << std::setprecision(15) << std::fixed;
elfouhaily_spectrum(100.0,0.0,0.85,kx, ky, S, PSI, PHI, Nx, Ny);

 PSI[Nx*Ny/2+Ny/2]=0.0; //avoid nan

 
 
double checksum=0.0;
 
for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++){
    
    
    checksum+=PSI[i*Ny+j];
    
  
  }
  
std::cout << "checksum: "<<checksum<<std::endl;

free(PHI);
free(S);
free(kx);
free(ky);

int num= 1.0;
 
 

double *PSI_COMPLEX  = (double*)malloc(2*Nx*Ny*sizeof(double));

for(int i=0; i < Ny; i++){   //fila 0
   
    PSI_COMPLEX[2*i] = PSI[i];
    PSI_COMPLEX[2*i+1]=0;
    }
    
for(int i=0; i < Nx; i++){   //columna 0
   
    PSI_COMPLEX[2*i*Ny] = PSI[i*Ny]; 
    PSI_COMPLEX[2*i*Ny+1]=0;
    }

for(int i=1; i < Ny/2; i++){   //centro simetrico
    unsigned int random = num*myRandom();
    PSI_COMPLEX[Nx*Ny+2*i] = PSI[Ny*Nx/2 +i]*cos(random);
    PSI_COMPLEX[Nx*Ny+2*i+1] = PSI[Ny*Nx/2 +i]*sin(random);
    
    PSI_COMPLEX[Nx*Ny+2*Ny-2*i] = PSI_COMPLEX[Nx*Ny+2*i];
    PSI_COMPLEX[Nx*Ny+2*Ny-2*i+1] = -PSI_COMPLEX[Nx*Ny+2*i+1];
  
  }
   PSI_COMPLEX[Nx*Ny+Ny]=PSI[Nx*Ny/2+Ny/2]; //real
  PSI_COMPLEX[Nx*Ny+Ny+1]=0;  //centro 0 //img
 
   
for(int i = 1; i < Nx/2; i++)
  for (int j = 1; j < Ny; j++) {
          unsigned int random = num*myRandom();
          PSI_COMPLEX[2*i*Ny+2*j]=PSI[i*Ny+j]*cos(random);
          PSI_COMPLEX[2*i*Ny+2*j+1]=PSI[i*Ny+j]*sin(random); 
          
          //simetric values
          PSI_COMPLEX[2*Nx*Ny-2*(i-1)*Ny-2*j]=PSI_COMPLEX[2*i*Ny+2*j];
          PSI_COMPLEX[2*Nx*Ny-2*(i-1)*Ny-2*j+1]=-PSI_COMPLEX[2*i*Ny+2*j+1];
          //std::cout<<PSI_COMPLEX[2*Nx*Ny-2*j+1];
          //return 0;
}
  

double checksum_complex=0.0;
 
for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++){
    
    
    checksum_complex+=PSI_COMPLEX[i*2*Ny+j*2];
    
  
  }  
  
  
std::cout << "checksum_complex: "<<checksum_complex<<std::endl;

matrix_transform(PSI_COMPLEX, Nx, Ny*2);


  
free(PSI);
for (int i = 0; i < Ny; i++)
  gsl_fft_complex_radix2_backward(PSI_COMPLEX+2*i,Ny,Nx);
for (int i = 0; i < Nx; i++)
 gsl_fft_complex_radix2_backward(PSI_COMPLEX+2*i*Ny,1,Ny);


 char name3[]="Wave_FFT_V1";

//plot_complex(Nx,Ny, PSI_COMPLEX, name3);

 

free(PSI_COMPLEX);
///////END verion 3////////
/*for(int i = 0; i < Nx; i++){
  for (int j = 0; j < Ny; j++) 
  {
      std::cout<<PSI_COMPLEX[2*i*Ny+2*j+1]<<" ";
  
  }
  std::cout<<endl;
  }
  */
return  0;
}


















