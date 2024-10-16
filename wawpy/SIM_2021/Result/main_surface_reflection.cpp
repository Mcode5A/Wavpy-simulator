#include "SurfaceReflection.hpp"


#include "math.h"
#include "stdio.h"
#include <iostream> 
#include <iomanip>
#include <complex>

using namespace std;




//****************************************************************************************************
// calculate_ponderacio  *2*PI/LANDA can be done later, 1 insteed of 3 times, can be merged with the division of module |xs -xt||xs-xr|
// so only on  div operation is computed,  
//ALSO CAN BE MOVED OUT OF INTEGRAL, SO ONLY IS CALCULATED ONE TIME INSTEAD OF N TIMES
//****************************************************************************************************


int main(){

cout << std::setprecision(15) << std::fixed;
int SIZE = 1000;
REAL Amp = 10.0;
REAL LANDA = 1.0;
REAL FREQ = 299792458.0/LANDA;

struct coor sat_coor_transmitter;
struct coor sat_coor_receiver; 
struct sea_surface surface;
struct coor polarization;
struct Vector_Fiel E_sum;

allocate_memory(&surface, SIZE);
    
init_Sat_ubi(&sat_coor_transmitter, 3, 36, 14);    //done!

init_Sat_ubi(&sat_coor_receiver,sqrt(SIZE),sqrt(SIZE), 10);  //done!

init_polarization(&polarization, 1,1,1);


std::complex<REAL> epsilon_permittivity_11(1.0, 0.0);
std::complex<REAL> epsilon_permittivity_22(1.0, 0.0);
std::complex<REAL> mu_permeability_11(1.0, 0.0);
std::complex<REAL> mu_permeability_22(1.0, 0.0);
REAL Temp = 20.0;
REAL salinity = 32.54;
epsilon_permittivity_22 = calculate_permittivity(Temp, salinity, FREQ);


std::complex<REAL> Eps_aux = (epsilon_permittivity_11*epsilon_permittivity_11)/(epsilon_permittivity_22*epsilon_permittivity_22);
std::complex<REAL> eta_SQ_11 = (mu_permeability_11/epsilon_permittivity_11); //eta^2
std::complex<REAL> eta_SQ_22 = (mu_permeability_22/epsilon_permittivity_22);
std::complex<REAL> A = eta_SQ_11*Eps_aux;
std::complex<REAL> B = eta_SQ_11*eta_SQ_11*Eps_aux/eta_SQ_22;
std::complex<REAL> eta_11=sqrt(eta_SQ_11);
std::complex<REAL> eta_22=sqrt(eta_SQ_22);
        
                                   
                
init_flat_surface(&surface, SIZE);  

calculate_incidence_unit_vec(&surface, &sat_coor_transmitter, SIZE);  

calculate_incidence_angle(&surface, SIZE);  

calculate_reflected_unit_vec(&surface, &sat_coor_receiver, SIZE); 

calculate_ponderacio(&surface, SIZE);  

calculate_propagation(&surface, SIZE, LANDA);  

calculate_polarization(&surface, &polarization, eta_11, eta_22, eta_SQ_11, eta_SQ_22, A, B, &E_sum, SIZE);
                   
      
    
std::cout<<Amp*E_sum.x.real()/(2.0*LANDA)<<" "<<Amp*E_sum.y.real()/(2.0*LANDA)<<" "<<Amp*E_sum.z.real()/(2.0*LANDA)<<" "<<std::endl;

 
 
 


}