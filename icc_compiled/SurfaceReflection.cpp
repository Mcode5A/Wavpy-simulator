#include "SurfaceReflection.hpp"
#include "math.h"
#include <complex>
#include <stdio.h>
#include <iostream>
#include <omp.h>
 
void init_Sat_ubi(struct coor * xyz, REAL x, REAL y, REAL z){

        xyz->x=x;
        xyz->y=y;
        xyz->z=z;}
void init_polarization(struct coor * polarization, REAL x, REAL y, REAL z){
    
    REAL modulo = x*x+y*y+z*z;
    modulo = sqrt(modulo);

    polarization->x=x/modulo;
    polarization->y=y/modulo;
    polarization->z=z/modulo;}
void init_flat_surface(struct sea_surface *s, int SIZE_XX, int SIZE_YY){

    ////////////initialization of pixel coordenates /////////////
  
        for(int i=0; i < SIZE_XX*SIZE_YY; i++){
           
           // surface
            s->s_z[i]=0.0;
            
            //normal vec
           // s->vec_L_x[i]=0.0;
           // s->vec_L_y[i]=0.0;
           // s->vec_L_z[i]=1.0;  

        } }
void init_flat_surface(REAL * s_z , int SIZE_XX, int SIZE_YY){

    ////////////initialization of pixel coordenates /////////////
  
        for(int i=0; i < SIZE_XX*SIZE_YY; i++){
           
           // surface
            s_z[i]=0.0;
            
            //normal vec
           // s->vec_L_x[i]=0.0;
           // s->vec_L_y[i]=0.0;
           // s->vec_L_z[i]=1.0;  

        } }
void init_sin_cos_surface(struct sea_surface *s, int SIZE_XX, int SIZE_YY){
       

    // EM WAVE LENGHT 30mm
    REAL period = two_PI;
    REAL sampling_rate = 16;
    REAL increment = 1/sampling_rate;
         


    int DIM_SIZE;// =  sqrt(SIZE);

    ////////////initialization of pixel coordenates    ECUATION Z= SIN(2*PI*X/0.3)*COS(2*PI*Y/0.3) /////////////
     
       //NOT WORKING
        for(int i=0; i <0*SIZE_XX*SIZE_YY; i++){

            
            s->s_x[i] = increment*(i/DIM_SIZE);  

            s->s_y[i] = increment*(i%DIM_SIZE);  
            
            s->s_z[i] = sin(period*s->s_x[i])*cos(period*s->s_y[i]);


            // Calculation normal vec
            s->vec_L_x[i] =  period*cos(period*s->s_x[i])*cos(period*s->s_y[i]);
        
            s->vec_L_y[i] = -period*sin(period*s->s_x[i])*sin(period*s->s_y[i]);
        
            s->vec_L_z[i] = -1;


            //// make normal vactor unit 
            REAL modulo = s->vec_L_x[i]*s->vec_L_x[i] +    //// modulo^2 normal_vec
                            s->vec_L_y[i]*s->vec_L_y[i] +
                            s->vec_L_z[i]*s->vec_L_z[i];

            modulo=sqrt(modulo);

            s->vec_L_x[i]/=modulo;
            s->vec_L_y[i]/=modulo;
            s->vec_L_z[i]/=modulo;
          
            
        }
        

 
    }
 
std::complex<REAL> calculate_permittivity(REAL temp, REAL salinity, REAL FREQ){

    

    //E_sw0(T,S) = E_sw0(T,0)*a(T,S)  E.22

    //*********************************************************************************************************//
    //E_sw0(T,0)  E.23
    REAL E_sw0_T0 = 87.134 - 0.1949*temp - 0.01276*temp*temp + 0.0002491*temp*temp*temp;
    
    //a(T,S) E24
    REAL a_TS = 1.0 + 0.00001613*temp*salinity - 0.003656*salinity + 0.0000321*salinity*salinity - 0.0000004232*salinity*salinity*salinity;

    //---//E_sw0(T,S) E.22
    REAL E_sw0_TS = E_sw0_T0*a_TS;
    //*********************************************************************************************************//

    //tsw(T,S) =  tsw(T,0)*b(T,S) E.25

    //*********************************************************************************************************//
    //2*M_PI*tsw(T,0) E.17
    REAL tsw_T0 = 0.00000000011109 - 0.000000000003824*temp + 0.00000000000006938*temp*temp - 0.0000000000000005096*temp*temp*temp;

    //b(T,S) E.26
    REAL b_TS = 1.0 + 0.00002282*temp*salinity - 0.0007638*salinity + 0.00000776*salinity*salinity - 0.00000001105*salinity*salinity*salinity;

    //---//2*M_PI*tsw(T,S) E.25 
    REAL tsw_TS = tsw_T0*b_TS;
    //*********************************************************************************************************//

    //s(T,S) = s(25,S)*e^-F  E.27

    //*********************************************************************************************************//
    //conductivity(25,S) E.28
    REAL conductivity_S = salinity * (0.18252 - 0.0014619*salinity + 0.00002093*salinity*salinity - 0.0000001282*salinity*salinity*salinity);

    //? E.28
    REAL delta = 25 - temp;
    REAL exponent = delta*(0.02033 + 0.0001266*delta + 0.000002464*delta*delta -salinity*(0.00001849 - 0.0000002551*delta + 0.00000002551*delta*delta) );


    REAL conductivity_TS = conductivity_S * exp(-exponent);

    //*********************************************************************************************************//
    
    
    //E.21AB COMMON PART
    REAL aux = (E_sw0_TS - 4.9) / (1+(tsw_TS*FREQ)*(tsw_TS*FREQ));

    //E.21A
    REAL E_sw_REAL = 4.9 + aux;

    //E.21B
    REAL E_sw_IMAG = (FREQ*tsw_TS*aux) + conductivity_TS/(2.0*M_PI*FREQ*0.000000000008854);

    std::complex<REAL> ref_indx(E_sw_REAL,E_sw_IMAG);
    
    return ref_indx;
    }

std::complex<REAL> calculate_refractive_index(std::complex<REAL> permittivity, std::complex<REAL> permeability) {

    return sqrt(permittivity*permeability);

    }

std::complex<REAL> calculate_media_impedance(std::complex<REAL> permittivity, std::complex<REAL> permeability) {

    return sqrt(permeability/permittivity);

    } 



 

void init_inv_fact(REAL * vec, int n){

    vec[0]=1.0;
    
    for(int i=1; i<n; i++)
      vec[i]= (REAL)(vec[i-1]*(REAL)i);
    
    for(int i=1; i<n; i++)
      vec[i]= 1.0/vec[i];
      
     for(int i=2; i<n; i+=4){
      vec[i]= -vec[i];
      vec[i+1]= -vec[i+1];
      }
}

//MAIN FUNCTION WITH BARREL OVER PHI AND THETA, GROUPS FUNCTIONS ABOVE
void compute_surface_reflection_Receiver_barrel(int SIZE_XX, int SIZE_YY, REAL Dx, REAL Dy, REAL deltaX, REAL deltaY, REAL init_point_XX, REAL init_point_YY, REAL LANDA, REAL Amp,  REAL Temp, REAL Salinity, struct coor *sat_coor_transmitter, struct coor *sat_coor_receiver, struct coor *polarization_inc, struct coor *polarization_ref, char dir_pol_RX, REAL * surface_s_z, REAL * barrel_result, int azimut, int inclination, REAL RX_RADIO){



    REAL FREQ = 299792458.0/LANDA;


    struct sea_surface surface_ptr;
    surface_ptr.s_z = surface_s_z;

    struct sea_surface * __restrict__ surface = &surface_ptr;
 

    std::complex<REAL> polarization_total_sum;

 

    std::complex<REAL> epsilon_permittivity_11(1.0, 0.0);
    std::complex<REAL> epsilon_permittivity_22(1.0, 0.0);
    std::complex<REAL> mu_permeability_11(1.0, 0.0);
    std::complex<REAL> mu_permeability_22(1.0, 0.0);
    epsilon_permittivity_22 = calculate_permittivity(Temp, Salinity, FREQ);


    std::complex<REAL> Eps_aux = (epsilon_permittivity_11*epsilon_permittivity_11)/(epsilon_permittivity_22*epsilon_permittivity_22);
    std::complex<REAL> eta_SQ_11 = (mu_permeability_11/epsilon_permittivity_11); //eta^2
    std::complex<REAL> eta_SQ_22 = (mu_permeability_22/epsilon_permittivity_22);
    std::complex<REAL> A = eta_SQ_11*Eps_aux;
    std::complex<REAL> B = eta_SQ_11*eta_SQ_11*Eps_aux/eta_SQ_22;
    std::complex<REAL> eta_11=sqrt(eta_SQ_11);
    std::complex<REAL> eta_22=sqrt(eta_SQ_22);
            
                                       
   // init_flat_surface(surface, SIZE_XX, SIZE_YY);  ///??????? ESTO BORRARIA TODO EL INPUT RECIBIDO
    //init_normal_vector_surface(surface, SIZE_XX, SIZE_YY, deltaX, deltaY);


    REAL * __restrict__ RX_coor_xx = (REAL *)malloc(inclination*azimut*sizeof(REAL));  
    REAL * __restrict__ RX_coor_yy = (REAL *)malloc(inclination*azimut*sizeof(REAL));  
    REAL * __restrict__ RX_coor_zz = (REAL *)malloc(inclination*azimut*sizeof(REAL));  
    REAL * __restrict__ polarization_ref_xx = (REAL *)malloc(inclination*azimut*sizeof(REAL));  
    REAL * __restrict__ polarization_ref_yy = (REAL *)malloc(inclination*azimut*sizeof(REAL));  
    REAL * __restrict__ polarization_ref_zz = (REAL *)malloc(inclination*azimut*sizeof(REAL));  

 
    if(dir_pol_RX == 'H')
      for(int k =0; k< inclination; k++){ 
                     
          for(int l =0; l< azimut; l++){ 
           
              RX_coor_xx[(k*azimut + l)]    =   -(RX_RADIO * sin(((REAL)k / inclination) * M_PI/2.0)*cos(((REAL)l / azimut) * 2.0*M_PI)); 
              RX_coor_yy[(k*azimut + l)]    =   -(RX_RADIO * sin(((REAL)k / inclination) * M_PI/2.0)*sin(((REAL)l / azimut) * 2.0*M_PI)); 
              RX_coor_zz[(k*azimut + l)]    =   (RX_RADIO * cos(((REAL)k / inclination) * M_PI/2.0)); 
              
              polarization_ref_xx[(k*azimut + l)]    =   sin(((REAL)l / azimut) * 2.0*M_PI); 
              polarization_ref_yy[(k*azimut + l)]    =   -cos(((REAL)l / azimut) * 2.0*M_PI); 
              polarization_ref_zz[(k*azimut + l)]    =   0.0; 
   
          } 
      } 

      if(dir_pol_RX == 'V')
        for(int k =0; k< inclination; k++){ 
                       
            for(int l =0; l< azimut; l++){ 
             
                RX_coor_xx[(k*azimut + l)]    =   -(RX_RADIO * sin(((REAL)k / inclination) * M_PI/2.0)*cos(((REAL)l / azimut) * 2.0*M_PI)); 
                RX_coor_yy[(k*azimut + l)]    =   -(RX_RADIO * sin(((REAL)k / inclination) * M_PI/2.0)*sin(((REAL)l / azimut) * 2.0*M_PI)); 
                RX_coor_zz[(k*azimut + l)]    =   (RX_RADIO * cos(((REAL)k / inclination) * M_PI/2.0)); 
                
                polarization_ref_xx[(k*azimut + l)]    =   cos(((REAL)k / inclination) * M_PI/2.0)*cos(((REAL)l / azimut) * 2.0*M_PI); 
                polarization_ref_yy[(k*azimut + l)]    =   cos(((REAL)k / inclination) * M_PI/2.0)*sin(((REAL)l / azimut) * 2.0*M_PI); 
                polarization_ref_zz[(k*azimut + l)]    =   sin(((REAL)k / inclination) * M_PI/2.0); 
     
            } 
        } 


    int n= 30;
    REAL * Fact_inv_sign = (REAL *)malloc(n*sizeof(REAL));
    init_inv_fact(Fact_inv_sign, n);

    REAL two_pi_landa = 2.0*M_PI/LANDA; 
    REAL inv_two_pi = 1.0/ (2.0*M_PI);
    
    int num_threads = 4;


    REAL **  BR_aux_buff = (REAL ** )malloc(num_threads * sizeof(REAL*));
    for(int i = 0; i < num_threads; i++) BR_aux_buff[i] = (REAL * )malloc(inclination*azimut*sizeof(REAL));


    #pragma omp parallel  num_threads(4)
    {
        int idthread = omp_get_thread_num();
        REAL * __restrict__ barrel_result_priv =  BR_aux_buff[idthread];

        #pragma omp for 
        for(int i=0; i < SIZE_XX; i++){

            REAL coor_xx_aux =  init_point_XX + deltaX*i;
            REAL coor_yy_aux =  init_point_YY - deltaY;
            REAL coor_zz_aux;

            struct coor v1; 
            struct coor v2; 
            
            
            for(int j=0; j < SIZE_YY; j++){
            
            
    
               // v1.x = 0.0;         //((REAL)i*deltaX)-((REAL)i*deltaX); 
                v1.y = 2.0*deltaY;  //((REAL)(j+1)*deltaY)-((REAL)(j-1)*deltaY);  
                v1.z = surface->s_z[i*SIZE_YY + (j + 1)%SIZE_YY] - surface->s_z[i*SIZE_YY +(j - 1)%SIZE_YY]; 
    
    
    
                v2.x = 2.0*deltaX; //((REAL)(i-1)*deltaX)-((REAL)(i+1)*deltaX);  
                //v2.y = 0.0;         //((REAL)(j)*deltaY)-((REAL)(j)*deltaY); 
                v2.z = surface->s_z[((i-1))*SIZE_YY + j] - surface->s_z[((i+1))*SIZE_YY + j]; 
                 
                //v2xv1 
                
                REAL normal_xx_aux =  (- v2.z*v1.y); 
                REAL normal_yy_aux =  -(-v2.x*v1.z ); 
                REAL normal_zz_aux = (v2.x*v1.y ); 
    
 
 
                REAL modulo =  sqrt((normal_xx_aux*normal_xx_aux) + (normal_yy_aux*normal_yy_aux) + (normal_zz_aux*normal_zz_aux) ); 
                                 
                 normal_xx_aux /= modulo; 
                 normal_yy_aux /= modulo; 
                 normal_zz_aux /= modulo; 

 
                coor_zz_aux = surface->s_z[i*SIZE_YY+j]; 
                coor_yy_aux+=deltaY;


                REAL  vec_inc_x = coor_xx_aux - sat_coor_transmitter->x; 
           
                REAL vec_inc_y = coor_yy_aux - sat_coor_transmitter->y ; 
         
                REAL vec_inc_z = coor_zz_aux - sat_coor_transmitter->z;

                REAL vec_module_inc = vec_inc_x*vec_inc_x +    //// modulo^2 incidente
                                                  vec_inc_y*vec_inc_y +
                                                  vec_inc_z*vec_inc_z;

                vec_module_inc=sqrt(vec_module_inc);
         
                REAL INV_vec_module_inc = 1.0/vec_module_inc;
               
                vec_inc_x *= INV_vec_module_inc;
            
                vec_inc_y *= INV_vec_module_inc;
            
                vec_inc_z *= INV_vec_module_inc;


                ////____________________________________________________________


                REAL producto_escalar = vec_inc_x*normal_xx_aux +   // producto escalar  incidente con normal 
                                        vec_inc_y*normal_yy_aux + 
                                        vec_inc_z*normal_zz_aux; 
            
           
                REAL angle=(-producto_escalar); 

                ////____________________________________________________________

                REAL sin_SQ_angle = 1.0-angle*angle;
                REAL cos_angle = angle;


                std::complex<REAL> Pinc_Vinc_RV = (eta_11*cos_angle - sqrt(eta_SQ_22 - A*sin_SQ_angle))/(eta_11*cos_angle + sqrt(eta_SQ_22 - A*sin_SQ_angle));

                std::complex<REAL> Pinc_Hinc_RH = (eta_22*cos_angle - sqrt(eta_SQ_11 - B*sin_SQ_angle))/(eta_22*cos_angle + sqrt(eta_SQ_11 - B*sin_SQ_angle));


                struct coor H_inc_ref;
                ///calculo Hinc con n x K1
                H_inc_ref.x =   normal_yy_aux*vec_inc_z - normal_zz_aux*vec_inc_y;    
                H_inc_ref.y = -(normal_xx_aux*vec_inc_z - normal_zz_aux*vec_inc_x);  
                H_inc_ref.z =   normal_xx_aux*vec_inc_y - normal_yy_aux*vec_inc_x;    


                ///calculo Href con hinc x K2 
                // Serni: Aquest càlcul està equivocat perquè H_ref = H_ref
                //Mouad: He creat una sola varible, com son igual s'utilitza la mateixa 

                struct coor V_inc;
                ///calculo Vinc con K1 x hinc
                V_inc.x =   vec_inc_y*H_inc_ref.z - vec_inc_z*H_inc_ref.y;    
                V_inc.y = -(vec_inc_x*H_inc_ref.z - vec_inc_z*H_inc_ref.x);  
                V_inc.z =   vec_inc_x*H_inc_ref.y - vec_inc_y*H_inc_ref.x;



                // El vector vref = Kref x href, on Kref és el direcció de reflexió especular del trosset de superfície.
                // https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector
                // Kref = K1 - 2*(K1·n)*n
                // Mouad: Aquest calcul es el que es fa per calcular l'angle incident, es podet reaprofitar el rpducte escalar
                REAL K1n = producto_escalar;
                /* vec_inc_x*surface->vec_L_x[i*SIZE_YY+j] +
                           vec_inc_y*surface->vec_L_y[i*SIZE_YY+j] + 
                           vec_inc_z*surface->vec_L_z[i*SIZE_YY+j]; */  // es el producte escalar entre K1 i n

                struct coor K_ref;  // Cálcul de Kref
                K_ref.x = vec_inc_x - 2.0*K1n*normal_xx_aux;
                K_ref.y = vec_inc_y - 2.0*K1n*normal_yy_aux;
                K_ref.z = vec_inc_z - 2.0*K1n*normal_zz_aux;

                struct coor V_ref;
                // Serni: i ara podem calcular Vref = Kref x href
                V_ref.x =   K_ref.y*H_inc_ref.z - K_ref.z*H_inc_ref.y;
                V_ref.y = -(K_ref.x*H_inc_ref.z - K_ref.z*H_inc_ref.x);
                V_ref.z =   K_ref.x*H_inc_ref.y - K_ref.y*H_inc_ref.x;


                //prodcuto escalar --- projecting the [p^] onto the [h^inc] and [v^inc] vectors--- Multiplication between fresnell coef. and polarization vectors
                Pinc_Vinc_RV *= polarization_inc->x*V_inc.x + polarization_inc->y*V_inc.y + polarization_inc->z*V_inc.z;
                Pinc_Hinc_RH *= polarization_inc->x*H_inc_ref.x + polarization_inc->y*H_inc_ref.y + polarization_inc->z*H_inc_ref.z;


                // Serni: falta revisar codi des d'aquí  --------------------------------------------------------------------


                std::complex<REAL> polarization_total_x = (Pinc_Hinc_RH * H_inc_ref.x + Pinc_Vinc_RV * V_ref.x);
                std::complex<REAL> polarization_total_y = (Pinc_Hinc_RH * H_inc_ref.y + Pinc_Vinc_RV * V_ref.y);
                std::complex<REAL> polarization_total_z = (Pinc_Hinc_RH * H_inc_ref.z + Pinc_Vinc_RV * V_ref.z);


               

                for(int l =0; l< azimut*inclination; l++){


                    REAL vec_ref_x =  RX_coor_xx[l] - coor_xx_aux;   
                 
                    REAL vec_ref_y =  RX_coor_yy[l] - coor_yy_aux;     
             
                    REAL vec_ref_z =  RX_coor_zz[l] - coor_zz_aux;   
              

                    REAL vec_module_ref=   vec_ref_x*vec_ref_x +    //// modulo^2 refeljado
                                            vec_ref_y*vec_ref_y +
                                            vec_ref_z*vec_ref_z;


                    vec_module_ref=sqrt(vec_module_ref);
                    
                    REAL INV_vec_module_ref= 1.0/vec_module_ref; 
               
                    vec_ref_x *= INV_vec_module_ref;
                
                    vec_ref_y *= INV_vec_module_ref;
                
                    vec_ref_z *= INV_vec_module_ref;


                    REAL K1_K2_x;
                    REAL K1_K2_y;  
                    REAL K1_K2_z;
                
                    K1_K2_x = (vec_inc_x - vec_ref_x);
                    K1_K2_y = (vec_inc_y - vec_ref_y);  
                    K1_K2_z = (vec_inc_z - vec_ref_z);

                    REAL vec_ponderacio =   -( K1_K2_x*normal_xx_aux +
                                                K1_K2_y*normal_yy_aux +
                                                K1_K2_z*normal_zz_aux); 


                    vec_ponderacio /= sqrt(K1_K2_x*K1_K2_x + K1_K2_y*K1_K2_y + K1_K2_z*K1_K2_z);


                    REAL fase = (vec_module_inc)+(vec_module_ref);
                    REAL INV_modulos = INV_vec_module_inc*INV_vec_module_ref;
                    
                    fase = fase*two_pi_landa;

                    
                    std::complex<REAL> polarization_total_sum =  polarization_total_x * polarization_ref_xx[l] + polarization_total_y * polarization_ref_yy[l] + polarization_total_z *polarization_ref_zz[l];

                    
                     
              
                    REAL coseno_fase  = cos(fase); 
                    
                
                  
                    REAL seno_fase = sin(fase); 
                    


                    
                    REAL E_field_real;
                    REAL E_field_imag;
                    
                    E_field_real=coseno_fase*INV_modulos*polarization_total_sum.real() - seno_fase*INV_modulos*polarization_total_sum.imag();
                    E_field_imag=coseno_fase*INV_modulos*polarization_total_sum.imag() + seno_fase*INV_modulos*polarization_total_sum.real();


                    barrel_result_priv[l] += (E_field_real*E_field_real+E_field_imag*E_field_imag)*vec_ponderacio;

                    

                }
            }
        }

        #pragma omp for 
        for(int k =0; k< inclination*azimut; k++){ 
               
            for(int l =0; l<num_threads ; l++){ 
                    
                    barrel_result[k] += BR_aux_buff[l][k]; 
         
                } 

            barrel_result[k]*=two_pi_landa;
        } 
        
    }

    



    for(int i = 0; i < num_threads; i++) 
        free(BR_aux_buff[i]);

    free(BR_aux_buff);

}




