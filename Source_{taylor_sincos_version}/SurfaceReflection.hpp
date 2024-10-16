#pragma once
#include <complex>
#define two_PI 2.0*M_PI
#define REAL double

struct coor{
    
    REAL x=0.0;
    REAL y=0.0;
    REAL z=0.0;

    };

struct complex_coor{
    
    std::complex<REAL> x=0.0;
    std::complex<REAL> y=0.0;
    std::complex<REAL> z=0.0;

    };

struct sea_surface {

    REAL * s_x;
    REAL * s_y;
    REAL * s_z; 

    //vecotr normal  !!UNITARIO!!, FACILITA CALCULOS POSTERIORES
    REAL * vec_L_x;  //
    REAL * vec_L_y;
    REAL * vec_L_z; 

    //vector incidente
    REAL * vec_inc_x;
    REAL * vec_inc_y;
    REAL * vec_inc_z; 

    //angulo incidente, && reflejado
    REAL * angle;
    
    //vector reflejado
    REAL * vec_ref_x;
    REAL * vec_ref_y;
    REAL * vec_ref_z; 


    //vector modulo incidente
    REAL * vec_module_inc;
     
    //vector modulo reflejado
    REAL * vec_module_ref;
     
    //vector ponderacio 
    REAL * vec_ponderacio;
 
    
    std::complex<REAL> * polarization_total_x;
    std::complex<REAL> * polarization_total_y;
    std::complex<REAL> * polarization_total_z;

 
    //Vector campo electrico
    std::complex<REAL> * E_field;
     
    };


void allocate_memory(struct sea_surface *s, int SIZE);
void free_memory(struct sea_surface *s, int SIZE);


void init_Sat_ubi(struct coor * xyz, REAL x, REAL y, REAL z);
void init_polarization(struct coor * polarization, REAL x, REAL y, REAL z);
void init_flat_surface(struct sea_surface *s, int SIZE_XX, int SIZE_YY);
void init_flat_surface(REAL *s_z , int SIZE_XX, int SIZE_YY);
void init_sin_cos_surface(struct sea_surface *s, int SIZE_XX, int SIZE_YY);
void init_normal_vector_surface(struct sea_surface *s, int SIZE_XX, int SIZE_YY, REAL deltaX, REAL deltaY);

void calculate_incidence_angle(struct sea_surface *s, int SIZE_XX, int SIZE_YY, REAL deltaX, REAL deltaY, REAL init_point_XX, REAL init_point_YY);
void calculate_incidence_unit_vec(struct sea_surface *s, struct coor *tx_xyz, int SIZE_XX, int SIZE_YY, REAL deltaX, REAL deltaY, REAL init_point_XX, REAL init_point_YY);
void calculate_reflected_unit_vec(struct sea_surface *s, struct coor *rx_xyz, int SIZE_XX, int SIZE_YY, REAL deltaX, REAL deltaY, REAL init_point_XX, REAL init_point_YY);
void calculate_module(struct sea_surface *s, bool inc, bool ref, int SIZE_XX, int SIZE_YY);



void calculate_ponderacio(struct sea_surface *s,  int SIZE_XX, int SIZE_YY);
void calculate_propagation(struct sea_surface *s, int SIZE_XX, int SIZE_YY, REAL LANDA);
void calculate_fresnell_coefficient_HV(std::complex<REAL> Ep11, std::complex<REAL> Ep22, std::complex<REAL> Mu11, std::complex<REAL> Mu22, REAL angle, std::complex<REAL> *RV, std::complex<REAL> *RH);
void calculate_polarization(struct sea_surface *s, struct coor * polarization_inc, REAL LANDA, std::complex<REAL> eta_11, std::complex<REAL> eta_22, std::complex<REAL> eta_SQ_11, std::complex<REAL> eta_SQ_22, std::complex<REAL> A, std::complex<REAL> B, int SIZE_XX, int SIZE_YY);


std::complex<REAL> calculate_permittivity(REAL temp, REAL salinity, REAL FREQ);
std::complex<REAL> calculate_refractive_index(std::complex<REAL> permittivity, std::complex<REAL> permeability);
std::complex<REAL> calculate_media_impedance(std::complex<REAL> permittivity, std::complex<REAL> permeability);



void compute_surface_reflection(int SIZE_XX, int SIZE_YY, REAL Dx, REAL Dy, REAL deltaX, REAL deltaY, REAL LANDA, REAL Amp,  REAL Temp, REAL Salinity, struct coor *sat_coor_transmitter, struct coor *sat_coor_receiver, struct coor *polarization_inc, struct coor * polarization_ref, REAL *surface_s_z);
void compute_surface_reflection_Receiver_barrel(int SIZE_XX, int SIZE_YY, REAL Dx, REAL Dy, REAL deltaX, REAL deltaY, REAL init_point_XX, REAL init_point_YY, REAL LANDA, REAL Amp,  REAL Temp, REAL Salinity, struct coor *sat_coor_transmitter, struct coor *sat_coor_receiver, struct coor *polarization_inc, struct coor *polarization_ref, char dir_pol_RX, REAL * surface_s_z, REAL * barrel_result, int azimut, int inclination, REAL RX_RADIO);