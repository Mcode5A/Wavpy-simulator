#define BOOST_TEST_MODULE ponderacio
#include <boost/test/unit_test.hpp>
#include "main.hpp"
#include <complex>
#include <iostream>
#define lindar 0.000000001

BOOST_AUTO_TEST_CASE(pond_prop_fresnell)
{

std::cout << std::setprecision(15) << std::fixed;
  //init data
  struct sea_surface s1;
  struct coor sat_coor_transmitter;
  struct coor sat_coor_receiver;
  allocate_memory(&s1);
  init_surface(&s1);
  int coor_xyz;
 
  complex<REAL>  Pinc_Vinc_RV, Pinc_Hinc_RH;
   
  /////TEST 1/////
  //redo with new values
  init_Sat_ubi(&sat_coor_transmitter, 4, 6, 9);     
  init_Sat_ubi(&sat_coor_receiver,9,9,9);  
  calculate_incidence_unit_vec(&s1, sat_coor_transmitter);  
  calculate_reflected_unit_vec(&s1, sat_coor_receiver);  
  calculate_incidence_angle(&s1);
  calculate_ponderacio(&s1);
  calculate_propagation(&s1);
   
  
  coor_xyz=6*DIM_SIZE+4; //surface(x,y)=(6,4)
  
  complex<REAL> eps11(1.0,0.0);
  complex<REAL> eps22(1.0,4.0);
  complex<REAL> mu11(1.0,0.0);
  complex<REAL> mu22(1.0,0.0);
  
  calculate_fresnell_coefficient_HV(eps11, eps22, mu11, mu22, s1.angle[coor_xyz],  &Pinc_Vinc_RV,  &Pinc_Hinc_RH);
  
  //angle 2.8371
  if ( fabs(s1.vec_ponderacio[coor_xyz]-(-11.2673372538)) > lindar  )
    BOOST_ERROR("ERROR FUNC_ponderacio_1_!");
   
  if ( fabs(s1.E_field[coor_xyz].real()-0.0054119563) > lindar)
    BOOST_ERROR("ERROR FUNC_propagation_1.1_!!"); //part real
  
  if ( fabs(s1.E_field[coor_xyz].imag()-0.0082713171) > lindar)
    BOOST_ERROR("ERROR FUNC_propagation_1.2_!!"); //part imag

  if ( fabs(Pinc_Vinc_RV.real() -0.3526951795) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnel_V_1.1_!!"); //part Vertical/parallel
  
  if ( (Pinc_Vinc_RV.imag() -0.3000209367) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnel_V_1.2_!!"); //part Vertical/parallel
    
  if ( fabs(Pinc_Hinc_RH.real() -(-0.3973820428)) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnelH_1.2_!!"); //part Horitzontal/normal
    
  if ( Pinc_Hinc_RH.imag() -(-0.3000759472) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnelH_1.2_!!"); //part Horitzontal/normal





 






   /////TEST 2///// 
   //redo with new values
  init_Sat_ubi(&sat_coor_transmitter, 2, 2, 7);     
  init_Sat_ubi(&sat_coor_receiver,15, 0, 8);  
  calculate_incidence_unit_vec(&s1, sat_coor_transmitter);  
  calculate_reflected_unit_vec(&s1, sat_coor_receiver);  
  calculate_incidence_angle(&s1);
  calculate_ponderacio(&s1);
  calculate_propagation(&s1);
  
  coor_xyz=20*DIM_SIZE+1; //surface(x,y)=(20,1)
  
  eps22.real(0.7);
  eps22.imag(8.7);
  mu22.real(9.0);
  mu22.imag(1.2);
  
  calculate_fresnell_coefficient_HV(eps11, eps22, mu11, mu22, s1.angle[coor_xyz],  &Pinc_Vinc_RV,  &Pinc_Hinc_RH);
  
  if ( fabs(s1.vec_ponderacio[coor_xyz]-(-7.5727175185 ))  > lindar )
    BOOST_ERROR("ERROR FUNC_ponderacio_2_!");  

  if ( fabs(s1.E_field[coor_xyz].real()-0.0025023189) > lindar)
    BOOST_ERROR("ERROR FUNC_propagation_2.1_!!"); //part real
    
  if ( fabs(s1.E_field[coor_xyz].imag()-(-0.0048422364)) > lindar)
    BOOST_ERROR("ERROR FUNC_propagation_2.2_!!"); //part imag
  
  if ( fabs(Pinc_Vinc_RV.real() -(-0.5201727043)) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnel_V_2.1_!!"); //part Vertical/parallel
     
   
  if ( (Pinc_Vinc_RV.imag() -0.3000209367) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnel_V_2.2_!!"); //part Vertical/parallel
  
  if ( fabs(Pinc_Hinc_RH.real() -(-0.5057101218)) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnelH_2.2_!!"); //part Horitzontal/normal
  
  if ( Pinc_Hinc_RH.imag() -(-0.273217251285163) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnelH_2.2_!!"); //part Horitzontal/normal

 
  
  
  
  
  
  /////TEST 3/////
  //redo with new values
  init_Sat_ubi(&sat_coor_transmitter, 3, 36, 14);     
  init_Sat_ubi(&sat_coor_receiver,34, 7, 4);  
  calculate_incidence_unit_vec(&s1, sat_coor_transmitter);  
  calculate_reflected_unit_vec(&s1, sat_coor_receiver);  
  calculate_incidence_angle(&s1);
  calculate_ponderacio(&s1);
  calculate_propagation(&s1);
  
  
  coor_xyz=16*DIM_SIZE+11; //surface(x,y)=(16,11)
  
  eps22.real(3.33);
  eps22.imag(4.9);
  mu22.real(1.0);
  mu22.imag(0.0);
  
  calculate_fresnell_coefficient_HV(eps11, eps22, mu11, mu22, s1.angle[coor_xyz],  &Pinc_Vinc_RV,  &Pinc_Hinc_RH);
  
 
  
  if ( fabs(s1.vec_ponderacio[coor_xyz]-(-4.1277309492)) > lindar) 
    BOOST_ERROR("ERROR FUNC_ponderacio_3_!"); 
  
  if ( fabs(s1.E_field[coor_xyz].real()-(-0.0008320695)) > lindar)
    BOOST_ERROR("ERROR FUNC_propagation_3.1_!!");  //part real
  
  if ( fabs(s1.E_field[coor_xyz].imag()-(0.0014645896)) > lindar)
    BOOST_ERROR("ERROR FUNC_propagation_3.2_!!"); //part imag
    
  if ( fabs(Pinc_Vinc_RV.real() -0.0604738133838627) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnel_V_3.1_!!"); //part Vertical/parallel
    
  if ( (Pinc_Vinc_RV.imag() -0.215863936227075) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnel_V_3.2_!!"); //part Vertical/parallel
  
  if ( fabs(Pinc_Hinc_RH.real() -(-0.709130255515308)) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnelH_3.1_!!"); //part Horitzontal/normal
  
  if ( Pinc_Hinc_RH.imag() -(-0.145048772833892) > lindar)
    BOOST_ERROR("ERROR FUNC__fresnelH_3.2_!!"); //part Horitzontal/normal

}
 