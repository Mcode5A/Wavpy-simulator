

// necesario punto del plano y vector normal

// COORDENADAS DE UN PUNTO DE LA SUPERFICIE, VECTOR NORMAL, VECTOR TX-SEA, VECTOR SEA-TR, VECTOR RELEXION

#define SIZE 100000000
#define DIM_SIZE (long)sqrt(SIZE)
#define LANDA 1
#define PI 3.14159265358979323846  
#define two_PI 2*PI

#define REAL float


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

 


void calculate_incidence_angle(struct sea_surface *s, struct coor tx_xyz){
    #pragma omp for
    for(long i=0; i < SIZE; i++){
        
        REAL producto_escalar = s->vec_inc_x[i]*s->vec_L_x[i] +   // producto escalar  incidente con normal
                    s->vec_inc_y[i]*s->vec_L_y[i] +
                    s->vec_inc_z[i]*s->vec_L_z[i];
    
    //cout<< producto_escalar<<endl;
        /*REAL modulo= s->vec_inc_x[i]*s->vec_inc_x[i] +    //// modulo^2 incidente
                    s->vec_inc_y[i]*s->vec_inc_y[i] +
                    s->vec_inc_z[i]*s->vec_inc_z[i];

        modulo=sqrt(modulo); //modulo incidente
        cout<<modulo  <<" "<<endl;*/
    // cout<<"mod:"<<modulo<<endl;
        s->angle[i]=(-producto_escalar);///modulo);//* 180.0 / PI;

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
        
       
            s->vec_L_z[i]=1;   }// integer division, change in float coord... case 

    
    ////////NEED OTHER FUNCTION THAT CALCULATES NORMAL VECTOR BASED ON SURFACE CARACTERISITCS, DEPENDES OF VALUES OF DATA INPUT SURFACE (FILE, FUNCTIONS....)   
    ////////////////////END INIT OF NORMAL...////////////////// 

    }

void init_surface2(struct sea_surface *s){
        // EM WAVE LENGHT 30mm
        REAL period = two_PI;
        REAL sampling_rate = 16;
        REAL increment = 1/sampling_rate;
         
    ////////////initialization of pixel coordenates    ECUATION Z= SIN(2*PI*X/0.3)*COS(2*PI*Y/0.3) /////////////
     
     
        for(long i=0; i < SIZE; i++){
            
            s->s_x[i] = increment*(i/DIM_SIZE); //-1.0+0.125*(i&255);

            s->s_y[i] = increment*(i%DIM_SIZE); //-1.0+0.125*(i&((2^63)-1-255));
            
            s->s_z[i] = sin(period*s->s_x[i])*cos(period*s->s_y[i]);
            
        }
        
    ////////////////////END INIT OF COOR...//////////////////

    
    ///////////////// INIT OF NORMAL VECTOR/////////////////////

        // surface =>   Z= SIN(2*PI*X/0.3)*COS(2*PI*Y/0.3)   => SIN(2*PI*X/0.3)*COS(2*PI*Y/0.3) - Z = 0
        // ??F(x,y,z) = (?/?x,?/?y,?/?z) F(x,y,x) ==> (?/?x,?/?y,?/?z) * SIN(2*PI*X/0.3)*COS(2*PI*Y/0.3) - Z 
        //
        // ??F(x,y,z) = ( [2*PI/0.3]*COS(2*PI*X/0.3)*COS(2*PI*Y/0.3) , [-2*PI/0.3]SIN(2*PI*X/0.3)*SIN(2*PI*Y/0.3) , -1) 
        //

        for(long i=0; i < SIZE; i++){

            s->vec_L_x[i] =  period*cos(period*s->s_x[i])*cos(period*s->s_y[i]);
        
            s->vec_L_y[i] = -period*sin(period*s->s_x[i])*sin(period*s->s_y[i]);
        
            s->vec_L_z[i] = -1;
        }
         

                //// make normal vactor unit 

                for(long i=0; i < SIZE; i++){

                    REAL modulo = s->vec_L_x[i]*s->vec_L_x[i] +    //// modulo^2 incidente
                                  s->vec_L_y[i]*s->vec_L_y[i] +
                                  s->vec_L_z[i]*s->vec_L_z[i];

                    modulo=sqrt(modulo);

                    s->vec_L_x[i]/=modulo;
                    s->vec_L_y[i]/=modulo;
                    s->vec_L_z[i]/=modulo;
                    }


                /////end unit normal vector 

    /*      for(long i=0; i < SIZE; i++)
            s->vec_L_x[i]= period*cos(period*X)*cos(period*Y);

        for(long i=0; i < SIZE; i++)
            s->vec_L_y[i]= -period*sin(period*X)*sin(period*Y);
        
        for(long i=0; i < SIZE; i++)
            s->vec_L_z[i]=-1;   // integer division, change in float coord... case */
    ////////NEED OTHER FUNCTION THAT CALCULATES NORMAL VECTOR BASED ON SURFACE CARACTERISITCS, DEPENDES OF VALUES OF DATA INPUT SURFACE (FILE, FUNCTIONS....)   
    ////////////////////END INIT OF NORMAL...////////////////// 

    }

void print_angle(struct sea_surface *s){



    for(long i=0; i<SIZE; i++){
    //  cout<<"Angle coordenate (x,y,z): ("<<s->s_x[i]<<", "<< s->s_y[i] << ", "<< s->s_z[i] <<") : " << s->angle[i] <<" grados"<<endl;
   //  printf("%s%lf%s%lf%s%lf%s%lf%s","Angle coordenate   (x,y,z): (",s->s_x[i],", ",s->s_y[i],", ",s->s_z[i],"): ", s->angle[i], " grados\n");
    
    } 
     //cout<<endl<<endl;
   //   printf("\n\n");
}

 

//done
void init_Sat_ubi(struct coor * xyz, REAL x, REAL y, REAL z){

        
        xyz->x=x;
        xyz->y=y;
        xyz->z=z;
}
void calculate_module(struct sea_surface *s, bool inc, bool ref){

    if(inc)
        #pragma SIMD
        for(long i=0; i < SIZE; i++){

            s->vec_module_inc[i]=   s->vec_inc_x[i]*s->vec_inc_x[i] +    //// modulo^2 incidente
                                    s->vec_inc_y[i]*s->vec_inc_y[i] +
                                    s->vec_inc_z[i]*s->vec_inc_z[i];

            s->vec_module_inc[i]=sqrt(s->vec_module_inc[i]);
        }
    
    if(ref)
        for(long i=0; i < SIZE; i++){

            s->vec_module_ref[i]=   s->vec_ref_x[i]*s->vec_ref_x[i] +    //// modulo^2 refeljado
                                    s->vec_ref_y[i]*s->vec_ref_y[i] +
                                    s->vec_ref_z[i]*s->vec_ref_z[i];

            s->vec_module_ref[i]=sqrt(s->vec_module_ref[i]);
        }


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

            s->vec_module_inc[i]=sqrt(s->vec_module_inc[i]);
     
            s->vec_inc_x[i] /= s->vec_module_inc[i];
       
            s->vec_inc_y[i] /= s->vec_module_inc[i];
      
            s->vec_inc_z[i] /= s->vec_module_inc[i];
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

            s->vec_module_ref[i]=sqrt(s->vec_module_ref[i]);
       
            s->vec_ref_x[i] /= s->vec_module_ref[i];
        
            s->vec_ref_y[i] /= s->vec_module_ref[i];
        
            s->vec_ref_z[i] /= s->vec_module_ref[i];
            
        }

        
    }
void calculate_ponderacio(struct sea_surface *s){


        ///MULTIPLICATION BY NORMAL VECTOR OF SURFACE, PONDREATION OF POWER RECIVED AT RX

 

        /////// K1 - K2 ////////
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

        REAL fase = (s->vec_module_inc[i])+(s->vec_module_ref[i]);
        REAL modulos = s->vec_module_inc[i]*s->vec_module_ref[i];

        
        fase = fase*2*PI/LANDA;


        s->E_field[i].real(cos(fase)/modulos);
        s->E_field[i].imag(sin(fase)/modulos);
          
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



   complex<REAL> sin_SQ_angle = sin(acos(angle))*sin(acos(angle));
    complex<REAL> cos_angle = cos(angle);
    

    *RV = (eta_11*cos_angle - sqrt(eta_SQ_22 - A*sin_SQ_angle))/(eta_11*cos_angle + sqrt(eta_SQ_22 - A*sin_SQ_angle));

    *RH = (eta_22*cos_angle - sqrt(eta_SQ_11 - B*sin_SQ_angle))/(eta_22*cos_angle + sqrt(eta_SQ_11 - B*sin_SQ_angle));


   // cout<<*RH<<" "<<*RV<<endl;
     
 }
 

void calculate_polarization(struct sea_surface *s, struct coor * polarization, complex<REAL> Ep11, complex<REAL> Ep22, complex<REAL> Mu11, complex<REAL> Mu22){

    struct coor H_inc;
    struct coor V_inc;
    struct coor H_ref;
    struct coor V_ref;
    complex<REAL> Pinc_Vinc_RV;
    complex<REAL> Pinc_Hinc_RH;
    
    #pragma omp  for  
    for(long i = 0; i<SIZE; i++){
    	
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
// so only on  div operation is computed
//****************************************************************************************************


int main(){

cout << std::setprecision(15) << std::fixed;

REAL Amp=10.0;
REAL Temp = 20.0;
REAL salinity = 32.54;
struct coor sat_coor_transmitter;
struct coor sat_coor_receiver; 
struct sea_surface surface;

struct coor polarization;
complex<REAL> refractive_index_n1(1.0, 0.0);

complex<REAL> epsilon_permittivity_22(3.33, 4.9);

complex<REAL> mu_permeability_11(1.0, 0.0);
complex<REAL> epsilon_permittivity_11(1.0, 0.0);

allocate_memory(&surface); 

complex<REAL> mu_permeability_22(1.0, 0.0);
complex<REAL> refractive_index_n2(0.0,0.0);


for(long i=0; i<1; i++){ 
          
          init_Sat_ubi(&sat_coor_transmitter, 3, 36, 14);    //done!
          
          init_Sat_ubi(&sat_coor_receiver,DIM_SIZE-1, DIM_SIZE-1, 10);  //done!
          
          init_polarization(&polarization, 1,1,1);
           
          //init_surface(&surface); //done!
          
           #pragma omp parallel num_threads(6)
          {
          
          
          init_surface(&surface); //2s
          
         
          calculate_incidence_unit_vec(&surface, sat_coor_transmitter); //done! 2.8s
          
          calculate_incidence_angle(&surface, sat_coor_transmitter); //done!  1s    // need coordenates of a surface point and satellit point  (angulo critico y angulo Brewster)
          
          calculate_reflected_unit_vec(&surface, sat_coor_receiver); //done! 2.8s
          
     
          
          //for(long i = 0; i < SIZE; i++)
          //    surface.angle[i]=i*5*PI/180;
          
          ///////////////////////////////////////////////////
          
         
          
          
          
          ///*  test ponderacio
          
           calculate_ponderacio(&surface); //2s
            
           //*/ 
           
           ///* test propagation
           
           calculate_propagation(&surface); //2s
            
          // */ 
           
           // /* test epsilon
            
           epsilon_permittivity_22 = calculate_permittivity(Temp, salinity);
             // cout<<"ESPSILON: "<<epsilon_permittivity<<endl;
              
           //*/
            
          
           
           
           //  /* test FRESNEL COEF. and polarization
              
           calculate_polarization(&surface, &polarization, epsilon_permittivity_11, epsilon_permittivity_22, mu_permeability_11, mu_permeability_22);
             
            // */
            
            }
}
 



for(long i=0; i < 0; i++){
    for(long j=0; j < DIM_SIZE; j++){}
//long i=1, j=2;

    //cout<<"SIN*COS*(K1-K2)*normal_surface_Vector/(mod1*mod2) ("<<i<<","<<j<<",z) ("<< surface.s_x[i*DIM_SIZE+ j]<<","<<surface.s_y[i*DIM_SIZE+ j]<<","<<surface.s_z[i*DIM_SIZE+ j]<< "):"<<surface.vec_ponderacio[i*DIM_SIZE+ j]*surface.E_real[i*DIM_SIZE+ j] << " +j"<<surface.vec_ponderacio[i*DIM_SIZE+ j]*surface.E_img[i*DIM_SIZE+ j]<<endl;
}



std::cout<<"aaa"<<endl<<endl;
//std::cout<<"       coordenadas    TX   (x,y,z): ("<<sat_coor_transmitter.x<<", "<<sat_coor_transmitter.y<<", "<<sat_coor_transmitter.z<<")"<<endl;
//printf("\n\n");

//printf("%s%lf%s%lf%s%lf%s","coordenadas    TX  (x,y,z): (",sat_coor_transmitter.x,", ",sat_coor_transmitter.y,", ",sat_coor_transmitter.z,")\n");


 
//print_angle(&surface);
 


}