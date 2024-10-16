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
#include "SurfaceGeneration.hpp"
#include "SurfaceReflection.hpp"
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
int Dx = (argc>6)?atoi(argv[6]): 2000;
int Dy = (argc>7)?atoi(argv[7]): 2000;//ha de ser un REAL
int Nx = (argc>8)?atoi(argv[8]): 1024;
int Ny = (argc>9)?atoi(argv[9]): 1024;



struct coor sat_coor_transmitter;
struct coor sat_coor_receiver; 
struct coor polarization_inc;
struct coor polarization_ref;

double elev = 0.0;

REAL RX_RADIO = 1500.0;
//pendiente inizializar
init_Sat_ubi(&sat_coor_transmitter,20000000*cos(elev*M_PI/180.0), 0, 20000000*sin(elev*M_PI/180.0));    //done!
//init_Sat_ubi(&sat_coor_receiver,0.0, 0.0, 0.0);  //done!
init_polarization(&polarization_inc, 0,-1,0);
//init_polarization(&polarization_ref, 1,1,1);


int SIZE = Nx*Ny;


 
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



double *surface = (double*)malloc(Nx*Ny*sizeof(double));  
// Serni: Té sentit utilitzar "struct sea_surface" i afegir-hi els paràmetres que la defineixen: Nx, Ny, Dx, Dy, wind theta, Omega, Rand_mod? Potser en una versió futura que ja sigui operativa es poden fer aquests canvis.

generate_surface(surface,  Nx,  Ny, wind, theta,  Omega,  Rand_mod,  Dx,  Dy );


//BARREL OUTPUT DATA
double *barrel_result = (double*)malloc(azimut*inclination*sizeof(double)); 


//HAY QUE HACER BARRIDO, E_FIELD DEBE SER MATRIZ, Y LAS COORDENADAS RECEIVER DEBEN ITERANDO
//--HAY QUE CAMBIAR EL INIT SURFACE, CALCULAR EL VECTOR NORMAL DE LA SUPERFICIE, con elfouhaily solo se genera la altura de cada punto, data el normal

//--HAY QUE CAMBIAR LAS FUNCIONES, ALGUNAS USAN LAS COORDENADS X E Y, LA SUPERFICIE SE PASA SOLO CON LA COORDENADA Z (ALTURA), HAY QUE USAR LOS INDICES DEL BUCLE PARA GENERAR LAS COOR. X E Y
//--HAY QUE AÑADIR UNA VARIABLE DE DISTANCIA ENTRE PUNTOS, EN ELFOUHAILY SE GENERAN Y HAY UNA DISTANCIA EN METROS TOTAL, SE PODRIA USAR EL MISMO DELTAX PARA LA REFLEXION
//--HAY QUE AÑADIR VARIABLES PARA UBICAR LAS COORDENADAS INCIALES, DONDE ESTA EL 0,0,0, QUE HAY EN LA POSICION 0 DE LA MATRIZ....
//HAY QUE MIRAR Dx y Dy son int o REAL

compute_surface_reflection_Receiver_barrel(Nx, Ny, Dx, Dy, LANDA, Amp, Temp, Salinity,  &sat_coor_transmitter, &sat_coor_receiver, &polarization_inc, &polarization_ref, surface, barrel_result, azimut, inclination, RX_RADIO);



//plot_data_semisphere_sample(azimut, inclination);


maxfield = INT_MIN;
minfield = INT_MAX;
for (int i = 0; i < azimut*inclination; i++)
{
  
 if(maxfield<barrel_result[i]) maxfield=barrel_result[i];
 if(minfield>barrel_result[i]) minfield=barrel_result[i];

}


std::cout<<" max: "<<  maxfield<<std::endl;
std::cout<<" min: "<<  minfield<<std::endl;
plot_data_semisphere(azimut, inclination,  maxfield, minfield,  barrel_result, &sat_coor_transmitter);


        

 

free(surface);
free(barrel_result); // Serni

return  0;
}






 
