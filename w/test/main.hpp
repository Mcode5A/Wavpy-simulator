

// necesario punto del plano y vector normal

// COORDENADAS DE UN PUNTO DE LA SUPERFICIE, VECTOR NORMAL, VECTOR TX-SEA, VECTOR SEA-TR, VECTOR RELEXION

#define SIZE 1000
#define SIZE2 2000 
#define DIM_SIZE (int)sqrt(SIZE)
#define LANDA 1
#define PI 3.14159265358979323846  
#define two_PI 2*PI

#define REAL double
 

#include "math.h"
#include "stdio.h"
#include <iostream> 
#include <iomanip>
#include <complex>

using namespace std;


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
 
   
 
    //Vector campo electrico
    complex<REAL> * E_field;
     
    }
;
 

struct coor{
    
    REAL x=0;
    REAL y=0;
    REAL z=0;

    }
;

 


void calculate_incidence_angle(struct sea_surface *s){
    #pragma omp for
    for(long i= 0; i < SIZE; i++){
        
        REAL producto_escalar = s->vec_inc_x[i]*s->vec_L_x[i] +   // producto escalar  incidente con normal
                    s->vec_inc_y[i]*s->vec_L_y[i] +
                    s->vec_inc_z[i]*s->vec_L_z[i];
    
   
        s->angle[i]=acos(-producto_escalar); 
      

    }  

}

void init_surface(struct sea_surface *s){

    ////////////initialization of pixel coordenates /////////////
        #pragma omp for
        for(long i=0; i < SIZE; i++){
            s->s_z[i]=0;
            s->s_y[i]=(float)(i%DIM_SIZE);
            s->s_x[i]=i/DIM_SIZE; 
            
        
            s->vec_L_x[i]=0;
        
     
            s->vec_L_y[i]=0;
        
       
            s->vec_L_z[i]=1;   } 

    
     
    }

void init_surface2(struct sea_surface *s){
        // EM WAVE LENGHT 30mm
        REAL period = two_PI;
        REAL sampling_rate = 16;
        REAL increment = 1/sampling_rate;
         
    ////////////initialization of pixel coordenates    ECUATION Z= SIN(2*PI*X/0.3)*COS(2*PI*Y/0.3) /////////////
     
     
        for(long i=0; i < SIZE; i++){

            
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
 
 

//done
void init_Sat_ubi(struct coor * xyz, REAL x, REAL y, REAL z){

        
        xyz->x=x;
        xyz->y=y;
        xyz->z=z;
}
 

void calculate_incidence_unit_vec(struct sea_surface *s, struct coor tx_xyz){


        //CALCULATE VECTOR, DEST-ORIGIN
        #pragma omp  for  
        for(long i=0; i < SIZE; i++){
            s->vec_inc_x[i] = s->s_x[i] - tx_xyz.x;
       
       
            s->vec_inc_y[i] = s->s_y[i] - tx_xyz.y;
      
     
            s->vec_inc_z[i] = s->s_z[i] - tx_xyz.z;
      

            s->vec_module_inc[i]=   s->vec_inc_x[i]*s->vec_inc_x[i] +    //// modulo^2 incidente
                                    s->vec_inc_y[i]*s->vec_inc_y[i] +
                                    s->vec_inc_z[i]*s->vec_inc_z[i];

            s->vec_module_inc[i]=1.0/sqrt(s->vec_module_inc[i]);
     
            s->vec_inc_x[i] *= s->vec_module_inc[i];
       
            s->vec_inc_y[i] *= s->vec_module_inc[i];
      
            s->vec_inc_z[i] *= s->vec_module_inc[i];
        }


    }
void calculate_reflected_unit_vec(struct sea_surface *s, struct coor rx_xyz){
        
        //CALCULATE VECTOR, DEST-ORIGIN
         #pragma omp  for  
        for(long i=0; i < SIZE; i++){
            s->vec_ref_x[i] =  rx_xyz.x - s->s_x[i];     
        
            s->vec_ref_y[i] =  rx_xyz.y - s->s_y[i];    
     
            s->vec_ref_z[i] =  rx_xyz.z - s->s_z[i];     
      

            s->vec_module_ref[i]=   s->vec_ref_x[i]*s->vec_ref_x[i] +    //// modulo^2 refeljado
                                    s->vec_ref_y[i]*s->vec_ref_y[i] +
                                    s->vec_ref_z[i]*s->vec_ref_z[i];

            s->vec_module_ref[i]=1.0/sqrt(s->vec_module_ref[i]);
       
            s->vec_ref_x[i] *= s->vec_module_ref[i];
        
            s->vec_ref_y[i] *= s->vec_module_ref[i];
        
            s->vec_ref_z[i] *= s->vec_module_ref[i];
            
        }

        
    }
void calculate_ponderacio(struct sea_surface *s){

      #pragma omp  for  
        for(long i=0; i < SIZE; i++){
       
            s->vec_ponderacio[i] = (s->vec_inc_x[i] - s->vec_ref_x[i])*s->vec_L_x[i] + (s->vec_inc_y[i] - s->vec_ref_y[i])*s->vec_L_y[i] +(s->vec_inc_z[i] - s->vec_ref_z[i])*s->vec_L_z[i] ; 
            s->vec_ponderacio[i] = s->vec_ponderacio[i]*two_PI/LANDA;
          
            
            
      
       } 


    }
//end done


//////NEW GROUP NEED TO CHECK DEEPER

void calculate_propagation(struct sea_surface *s ){
    #pragma omp  for  
    for(long i=0; i < SIZE; i++) {

        REAL fase = (1.0/s->vec_module_inc[i])+(1.0/s->vec_module_ref[i]);
        REAL modulos = s->vec_module_inc[i]*s->vec_module_ref[i];

        
        fase = fase*2*PI/LANDA;


        s->E_field[i].real(cos(fase)*modulos);
        s->E_field[i].imag(sin(fase)*modulos);
          
    }



}

complex<REAL>  calculate_permittivity(REAL temp, REAL salinity){

    REAL freq=10000000000;

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
    //2*PI*tsw(T,0) E.17
    REAL tsw_T0 = 0.00000000011109 - 0.000000000003824*temp + 0.00000000000006938*temp*temp - 0.0000000000000005096*temp*temp*temp;

    //b(T,S) E.26
    REAL b_TS = 1.0 + 0.00002282*temp*salinity - 0.0007638*salinity + 0.00000776*salinity*salinity - 0.00000001105*salinity*salinity*salinity;

    //---//2*PI*tsw(T,S) E.25 
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
    REAL aux = (E_sw0_TS - 4.9) / (1+(tsw_TS*freq)*(tsw_TS*freq));

    //E.21A
    REAL E_sw_REAL = 4.9 + aux;

    //E.21B
    REAL E_sw_IMAG = (freq*tsw_TS*aux) + conductivity_TS/(two_PI*freq*0.000000000008854);

    complex<REAL> ref_indx(E_sw_REAL,E_sw_IMAG);
    
    return ref_indx;
}

complex<REAL> calculate_refractive_index(complex<REAL> permittivity, complex<REAL> permeability) {

    return sqrt(permittivity*permeability);

}

complex<REAL> calculate_media_impedance(complex<REAL> permittivity, complex<REAL> permeability) {

    return sqrt(permeability/permittivity);

}

void  calculate_fresnell_coefficient_HV(complex<REAL> Ep11, complex<REAL> Ep22, complex<REAL> Mu11, complex<REAL> Mu22, REAL angle, complex<REAL> *RV, complex<REAL> *RH)
{
    
   
    complex<REAL> aux = (Ep11*Ep11)/(Ep22*Ep22);

    complex<REAL> eta_SQ_11 = (Mu11/Ep11); //eta^2
    complex<REAL> eta_SQ_22 = (Mu22/Ep22);

    complex<REAL> A = eta_SQ_11*aux;
    complex<REAL> B = eta_SQ_11*eta_SQ_11*aux/eta_SQ_22;

    complex<REAL> eta_11=sqrt(eta_SQ_11);
    complex<REAL> eta_22=sqrt(eta_SQ_22);



    complex<REAL> sin_SQ_angle = sin((angle))*sin((angle));
    complex<REAL> cos_angle = cos(angle);
    

    *RV = (eta_11*cos_angle - sqrt(eta_SQ_22 - A*sin_SQ_angle))/(eta_11*cos_angle + sqrt(eta_SQ_22 - A*sin_SQ_angle));

    *RH = (eta_22*cos_angle - sqrt(eta_SQ_11 - B*sin_SQ_angle))/(eta_22*cos_angle + sqrt(eta_SQ_11 - B*sin_SQ_angle));

     
       
 }
 

void calculate_polarization(struct sea_surface *s, struct coor * polarization, complex<REAL> Ep11, complex<REAL> Ep22, complex<REAL> Mu11, complex<REAL> Mu22){

    struct coor H_inc;
    struct coor V_inc;
    struct coor H_ref;
    struct coor V_ref;
    complex<REAL> Pinc_Vinc_RV;
    complex<REAL> Pinc_Hinc_RH;
    
    #pragma omp  for  
    for(long i =0; i<SIZE; i++){
    	
     	
      
    	calculate_fresnell_coefficient_HV(Ep11, Ep22, Mu11, Mu22, s->angle[i],  &Pinc_Vinc_RV,  &Pinc_Hinc_RH);

        ///calculo Hinc con n x K1
        H_inc.x =   s->vec_L_y[i]*s->vec_inc_z[i] - s->vec_L_z[i]*s->vec_inc_y[i];    
        H_inc.y = -(s->vec_L_x[i]*s->vec_inc_z[i] - s->vec_L_z[i]*s->vec_inc_x[i]);  
        H_inc.z =   s->vec_L_x[i]*s->vec_inc_y[i] - s->vec_L_y[i]*s->vec_inc_x[i];    

        ///calculo Href con hinc x K2
        H_ref.x =  s->vec_L_y[i]*s->vec_ref_z[i] - s->vec_L_z[i]*s->vec_ref_y[i];    
        H_ref.y = -(s->vec_L_x[i]*s->vec_ref_z[i] - s->vec_L_z[i]*s->vec_ref_x[i]);  
        H_ref.z =  s->vec_L_x[i]*s->vec_ref_y[i] - s->vec_L_y[i]*s->vec_ref_x[i];   

        ///calculo Vinc con n x K1
        V_inc.x =  s->vec_inc_y[i]*H_inc.z - s->vec_inc_z[i]*H_inc.y;    
        V_inc.y = -(s->vec_inc_x[i]*H_inc.z - s->vec_inc_z[i]*H_inc.x);  
        V_inc.z =  s->vec_inc_x[i]*H_inc.y - s->vec_inc_y[i]*H_inc.x;

        
        ///calculo Vref con hinc x K2
        V_ref.x =  s->vec_ref_y[i]*H_ref.z - s->vec_ref_z[i]*H_ref.y;     
        V_ref.y = -(s->vec_ref_x[i]*H_ref.z - s->vec_ref_z[i]*H_ref.x);  
        V_ref.z =  s->vec_ref_x[i]*H_ref.y - s->vec_ref_y[i]*H_ref.x; 

        //prodcuto escalar
        Pinc_Vinc_RV *= polarization->x*V_inc.x + polarization->y*V_inc.y + polarization->z*V_inc.z;
		    Pinc_Hinc_RH *= polarization->x*H_inc.x + polarization->y*H_inc.y + polarization->z*H_inc.z;


        s->E_field[i] *= s->vec_ponderacio[i];

        s->E_field[i] * (Pinc_Hinc_RH * H_ref.x + Pinc_Vinc_RV * V_ref.x);
        s->E_field[i] * (Pinc_Hinc_RH * H_ref.y + Pinc_Vinc_RV * V_ref.y);
        s->E_field[i] * (Pinc_Hinc_RH * H_ref.z + Pinc_Vinc_RV * V_ref.z);

        
        
        
 
      //cout<< "E(x,y,z): "<< (Pinc_Hinc_RH * H_ref.x + Pinc_Vinc_RV * V_ref.x)<<"  ,  "<< (Pinc_Hinc_RH * H_ref.y + Pinc_Vinc_RV * V_ref.y)<<"  ,  "<< (Pinc_Hinc_RH * H_ref.z + Pinc_Vinc_RV * V_ref.z)<<endl;
    
    }



}
 

void init_polarization(struct coor * polarization, REAL x, REAL y, REAL z){
    
    REAL modulo = x*x+y*y+z*z;
    modulo = sqrt(modulo);

    polarization->x=x/modulo;
    polarization->y=y/modulo;
    polarization->z=z/modulo;

}
  

void allocate_memory(struct sea_surface *s){

     s->s_x = (REAL *)malloc(SIZE*sizeof(REAL));
     s->s_y = (REAL *)malloc(SIZE*sizeof(REAL));
     s->s_z = (REAL *)malloc(SIZE*sizeof(REAL)); 

    
     s->vec_L_x = (REAL *)malloc(SIZE*sizeof(REAL));  
     s->vec_L_y = (REAL *)malloc(SIZE*sizeof(REAL));
     s->vec_L_z = (REAL *)malloc(SIZE*sizeof(REAL)); 

    
    s->vec_inc_x = (REAL *)malloc(SIZE*sizeof(REAL));
    s->vec_inc_y = (REAL *)malloc(SIZE*sizeof(REAL));
    s->vec_inc_z = (REAL *)malloc(SIZE*sizeof(REAL)); 

    
     s->angle = (REAL *)malloc(SIZE*sizeof(REAL));
    
    
     s->vec_ref_x = (REAL *)malloc(SIZE*sizeof(REAL));
     s->vec_ref_y = (REAL *)malloc(SIZE*sizeof(REAL));
     s->vec_ref_z = (REAL *)malloc(SIZE*sizeof(REAL)); 


    
     s->vec_module_inc = (REAL *)malloc(SIZE*sizeof(REAL));
     
    
     s->vec_module_ref = (REAL *)malloc(SIZE*sizeof(REAL));
     
    
     s->vec_ponderacio = (REAL *)malloc(SIZE*sizeof(REAL));
 
    
     s->E_field = (complex<REAL> *)malloc(2*SIZE*sizeof(REAL));
     


} 


 



//****************************************************************************************************
// calculate_ponderacio  *2*PI/LANDA can be done later, 1 insteed of 3 times, can be merged with the division of module |xs -xt||xs-xr|
// so only on  div operation is computed,  
//ALSO CAN BE MOVED OUT OF INTEGRAL, SO ONLY IS CALCULATED ONE TIME INSTEAD OF N TIMES
//****************************************************************************************************

