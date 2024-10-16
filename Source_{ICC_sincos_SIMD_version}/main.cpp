#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <iomanip>
#include <string> 
#include <complex>
#include <cstring>
#include <limits>
#include <mpi.h>
#include "SurfaceGeneration.hpp"
#include "SurfaceReflection.hpp"
#include <unistd.h>
#include "plot_semisphere.hpp"



using namespace std;



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
std::cout << std::setprecision(15) << std::fixed;

char name[15] = {"Output"};
if (argc>1) strcpy(name,argv[1]);
double wind  = (argc>2)?strtod(argv[2], NULL): 3.0;
double theta = (argc>3)?strtod(argv[3], NULL): 0.0;
double Omega = (argc>4)?strtod(argv[4], NULL): 0.85;
double Rand_mod= (argc>5)?strtod(argv[5], NULL): 1.0;
int Dx = (argc>6)?atoi(argv[6]): 150;
int Dy = (argc>7)?atoi(argv[7]): 150;//ha de ser un REAL
int Nx = (argc>8)?atoi(argv[8]): 2048;
int Ny = (argc>9)?atoi(argv[9]): 2048;



//HAY QUE MIRAR Dx y Dy son int o REAL

struct coor sat_coor_transmitter;
struct coor sat_coor_receiver; 
struct coor polarization_inc;
struct coor polarization_ref;

double elev = 45.0;

REAL RX_RADIO = 1000.0;
//pendiente inizializar
init_Sat_ubi(&sat_coor_transmitter,20000000*cos(elev*M_PI/180.0), 0, 20000000*sin(elev*M_PI/180.0));    //done!
//init_Sat_ubi(&sat_coor_receiver,0.0, 0.0, 0.0);  //done!
//init_polarization(&polarization_inc, 0,-1,0);
//init_polarization(&polarization_ref, 1,1,1);



char pol_dir_TX = 'H';
char pol_dir_RX = 'H';



if(pol_dir_TX == 'V')
  init_polarization(&polarization_inc, -sin(elev*M_PI/180.0),0,cos(elev*M_PI/180.0));
  
if(pol_dir_TX == 'H')
  init_polarization(&polarization_inc, 0,-1,0);

//// INCIALIAZAR DESDE CONSOLA, O PARAMETRIZAR LA LLAMADA
REAL LANDA = 0.19; //pendiente inicializar
REAL Amp = 10.0; //pendiente inicializar
REAL Temp = 20.0;
REAL Salinity = 32.54; 

// EL 10 DEBE SER MODIFICABLE, DEPENDIENTE DE LA GRANULARIDAD DEL BARRIDO
int azimut = 360/10;
int inclination = 90/10;
REAL maxfield = Amp; //???????
REAL minfield = -Amp; //???????
REAL init_point_XX = -((REAL)Dx/2.0);
REAL init_point_YY = -((REAL)Dy/2.0);

double deltaX = (double)Dx/Nx;
double deltaY = (double)Dy/Ny;
  


int rank, size_procs;
MPI_Init( &argc, &argv );
MPI_Comm_rank( MPI_COMM_WORLD, &rank );
MPI_Comm_size( MPI_COMM_WORLD, &size_procs );
MPI_Status  status;
MPI_Request request;



double * surface  = (double*)malloc((2*Nx*Ny/size_procs + 2*Nx)*sizeof(double));  
double * original_surface_ptr = surface; 
surface += Nx;




generate_surface_merge_mpi(surface, Nx/size_procs,  Ny,  rank,  wind,  theta,  Omega,  Rand_mod,  Dx,  Dy,  deltaX,  deltaY, Rand_mod, size_procs);



// Serni: Té sentit utilitzar "struct sea_surface" i afegir-hi els paràmetres que la defineixen: Nx, Ny, Dx, Dy, wind theta, Omega, Rand_mod? Potser en una versió futura que ja sigui operativa es poden fer aquests canvis.



//removing imag part to recover memory
for(int i=0; i<Nx*Ny/size_procs; i++){  
    surface[i]=surface[i*2];
}
 
original_surface_ptr = (double *) realloc(original_surface_ptr, ((Nx*Ny/size_procs)+2*Nx)*sizeof(double));



//Sending last & firstrow of each for normal vec calc
int prev=  (rank-1)%size_procs;
int next =  (rank+1)%size_procs;
if(prev == -1) prev = size_procs-1;

MPI_Isend(surface+(Nx*Ny/size_procs)-Nx, Nx, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &request);
MPI_Recv(surface-Nx, Nx, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &status);

 
MPI_Isend(surface, Nx, MPI_DOUBLE, prev, 1, MPI_COMM_WORLD, &request);
MPI_Recv(surface+Nx*Ny/size_procs, Nx, MPI_DOUBLE, next, 1, MPI_COMM_WORLD, &status);




//assignar cordenada inicial despues de la division de la superficie //solo X porque cada uno tiene n filas completas
if(size_procs>1){
init_point_XX += ((REAL)Dy/size_procs)*rank;
int aux = Nx;
Nx=Ny; 
Ny=aux;
}


//BARREL OUTPUT DATA
double *barrel_result = (double*)malloc(azimut*inclination*sizeof(double)); 
double *barrel_result_Reduce = (double*)malloc(azimut*inclination*sizeof(double)); 


compute_surface_reflection_Receiver_barrel(Nx/size_procs,Ny, Dx, Dy, deltaX, deltaY, init_point_XX, init_point_YY, LANDA, Amp, Temp, Salinity,  &sat_coor_transmitter, &sat_coor_receiver, &polarization_inc, &polarization_ref, pol_dir_RX, surface, barrel_result, azimut, inclination, RX_RADIO);


MPI_Allreduce(barrel_result, barrel_result_Reduce, azimut*inclination, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

free(barrel_result);

MPI_Finalize();



maxfield = DBL_MIN;
minfield = DBL_MAX;
for (int i = 0; i < azimut*inclination; i++)
{
  
 if(maxfield<barrel_result_Reduce[i]) maxfield=barrel_result_Reduce[i];
 if(minfield>barrel_result_Reduce[i]) minfield=barrel_result_Reduce[i];

}


std::cout<<" max: "<<  maxfield<<std::endl;
std::cout<<" min: "<<  minfield<<std::endl;
plot_data_semisphere(azimut, inclination,  maxfield, minfield,  barrel_result_Reduce, &sat_coor_transmitter);


        

 

free(original_surface_ptr);
free(barrel_result_Reduce);
 // Serni




return  0;
}
