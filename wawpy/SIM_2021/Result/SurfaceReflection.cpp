#include "SurfaceReflection.hpp"
#include "math.h"
#include <complex>
#include <stdio.h>
#include <iostream>

void allocate_memory(struct sea_surface *s, int SIZE){
      
    // s->s_x = (REAL *)malloc(SIZE*sizeof(REAL));
    // s->s_y = (REAL *)malloc(SIZE*sizeof(REAL));
    // s->s_z = (REAL *)malloc(SIZE*sizeof(REAL)); 

    
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
 
    
     s->E_field = (std::complex<REAL> *)malloc(2*SIZE*sizeof(REAL));
 }

void free_memory(struct sea_surface *s, int SIZE){
      
    // s->s_x = (REAL *)malloc(SIZE*sizeof(REAL));
    // s->s_y = (REAL *)malloc(SIZE*sizeof(REAL));
    // s->s_z = (REAL *)malloc(SIZE*sizeof(REAL)); 

    
    free( s->vec_L_x);
    free( s->vec_L_y);
    free( s->vec_L_z);
    free( s->vec_inc_x);
    free(s->vec_inc_y);
    free(s->vec_inc_z);
    free(s->angle);
    free( s->vec_ref_x);
    free( s->vec_ref_y);
    free(s->vec_ref_z);
    free(s->vec_module_inc);
    free(s->vec_module_ref);
    free(s->vec_ponderacio);
    free(s->E_field);
 }


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

void init_normal_vector_surface(struct sea_surface *s, int SIZE_XX, int SIZE_YY, REAL deltaX, REAL deltaY){

    struct coor v1;
    struct coor v2;


    //inner matrix normal vector calculation
    for(int i=1; i<SIZE_XX-1; i++){
          
          for(int j=1; j<SIZE_YY-1; j++){


                v1.x = 0.0;         //((REAL)i*deltaX)-((REAL)i*deltaX);
                v1.y = 2.0*deltaY;  //((REAL)(j+1)*deltaY)-((REAL)(j-1)*deltaY); 
                v1.z = s->s_z[i*SIZE_YY + j + 1] - s->s_z[i*SIZE_YY + j - 1];



                v2.x = -2.0*deltaX; //((REAL)(i-1)*deltaX)-((REAL)(i+1)*deltaX); 
                v2.y = 0.0;         //((REAL)(j)*deltaY)-((REAL)(j)*deltaY);
                v2.z = s->s_z[(i-1)*SIZE_YY + j] - s->s_z[(i+1)*SIZE_YY + j];
                
                //v2xv1

                s->vec_L_x[i*SIZE_YY + j] =  -(v2.y*v1.z - v2.z*v1.y);
                s->vec_L_y[i*SIZE_YY + j] = -(-v2.x*v1.z + v2.z*v1.x);
                s->vec_L_z[i*SIZE_YY + j] =  -(v2.x*v1.y - v2.y*v1.x);

 

          }
    }


    //upper boundary  normal vector calculation
    for(int i=0; i<1; i++){

        for(int j=1; j<SIZE_YY-1; j++){


                v1.x = 0.0;         //((REAL)i*deltaX)-((REAL)i*deltaX);
                v1.y = 2.0*deltaY;  //((REAL)(j+1)*deltaY)-((REAL)(j-1)*deltaY); 
                v1.z = s->s_z[i*SIZE_YY + j + 1] - s->s_z[i*SIZE_YY + j - 1];



                v2.x = -2.0*deltaX; //((REAL)(i-1)*deltaX)-((REAL)(i+1)*deltaX); 
                v2.y = 0.0;         //((REAL)(j)*deltaY)-((REAL)(j)*deltaY);
                v2.z = s->s_z[(SIZE_XX-1)*SIZE_YY + j] - s->s_z[(i+1)*SIZE_YY + j];

                //v2xv1

                s->vec_L_x[i*SIZE_YY + j] =  v2.z*v1.y - v2.y*v1.z;
                s->vec_L_y[i*SIZE_YY + j] = -v2.z*v1.x + v2.x*v1.z;
                s->vec_L_z[i*SIZE_YY + j] =  v2.y*v1.x - v2.x*v1.y;

        }


    }


    // lower boundary  normal vector calculation
    for(int i=SIZE_XX-1; i<SIZE_XX; i++){

        for(int j=1; j<SIZE_YY-1; j++){


                v1.x = 0.0;         //((REAL)i*deltaX)-((REAL)i*deltaX);
                v1.y = 2.0*deltaY;  //((REAL)(j+1)*deltaY)-((REAL)(j-1)*deltaY); 
                v1.z = s->s_z[i*SIZE_YY + j + 1] - s->s_z[i*SIZE_YY + j - 1];



                v2.x = -2.0*deltaX; //((REAL)(i-1)*deltaX)-((REAL)(i+1)*deltaX); 
                v2.y = 0.0;         //((REAL)(j)*deltaY)-((REAL)(j)*deltaY);
                v2.z = s->s_z[(i-1)*SIZE_YY + j] - s->s_z[j]; //(0)*SIZE_YY + j

                //v2xv1

                s->vec_L_x[i*SIZE_YY + j] =  v2.z*v1.y - v2.y*v1.z;
                s->vec_L_y[i*SIZE_YY + j] = -v2.z*v1.x + v2.x*v1.z;
                s->vec_L_z[i*SIZE_YY + j] =  v2.y*v1.x - v2.x*v1.y;

        }


    }

    // left boundary  normal vector calculation
        for(int i=1; i<SIZE_XX-1; i++){
              
              for(int j=0; j<1; j++){


                    v1.x = 0.0;         //((REAL)i*deltaX)-((REAL)i*deltaX);
                    v1.y = 2.0*deltaY;  //((REAL)(j+1)*deltaY)-((REAL)(j-1)*deltaY); 
                    v1.z = s->s_z[i*SIZE_YY + j + 1] - s->s_z[i*SIZE_YY + SIZE_YY - 1];



                    v2.x = -2.0*deltaX; //((REAL)(i-1)*deltaX)-((REAL)(i+1)*deltaX); 
                    v2.y = 0.0;         //((REAL)(j)*deltaY)-((REAL)(j)*deltaY);
                    v2.z = s->s_z[(i-1)*SIZE_YY + j] - s->s_z[(i+1)*SIZE_YY + j];
                    
                    //v2xv1

                    s->vec_L_x[i*SIZE_YY + j] =  v2.z*v1.y - v2.y*v1.z;
                    s->vec_L_y[i*SIZE_YY + j] = -v2.z*v1.x + v2.x*v1.z;
                    s->vec_L_z[i*SIZE_YY + j] =  v2.y*v1.x - v2.x*v1.y;

     

              }
        }

    // right boundary  normal vector calculation
        for(int i=1; i<SIZE_XX-1; i++){
              
              for(int j=SIZE_YY-1; j<SIZE_YY; j++){


                    v1.x = 0.0;         //((REAL)i*deltaX)-((REAL)i*deltaX);
                    v1.y = 2.0*deltaY;  //((REAL)(j+1)*deltaY)-((REAL)(j-1)*deltaY); 
                    v1.z = s->s_z[i*SIZE_YY + 0] - s->s_z[i*SIZE_YY + j - 1];



                    v2.x = -2.0*deltaX; //((REAL)(i-1)*deltaX)-((REAL)(i+1)*deltaX); 
                    v2.y = 0.0;         //((REAL)(j)*deltaY)-((REAL)(j)*deltaY);
                    v2.z = s->s_z[(i-1)*SIZE_YY + j] - s->s_z[(i+1)*SIZE_YY + j];
                    
                    //v2xv1

                    s->vec_L_x[i*SIZE_YY + j] =  v2.z*v1.y - v2.y*v1.z;
                    s->vec_L_y[i*SIZE_YY + j] = -v2.z*v1.x + v2.x*v1.z;
                    s->vec_L_z[i*SIZE_YY + j] =  v2.y*v1.x - v2.x*v1.y;
     

              }
        }




    //4 corners normal vector calculation

            //SAME FOR ALL 4 CORNERS
            v1.x = 0.0;                          //((REAL)i*deltaX)-((REAL)i*deltaX);
            v1.y = -(REAL)(SIZE_YY-2)*deltaY;    //((REAL)(j+1)*deltaY)-((REAL)(j-1)*deltaY);
            
            v2.x = ((REAL)(SIZE_XX-2)*deltaX);   //((REAL)(i-1)*deltaX)-((REAL)(i+1)*deltaX); 
            v2.y = 0.0;                          //((REAL)(j)*deltaY)-((REAL)(j)*deltaY);
    
            //left-up corner

                v1.z = s->s_z[1] - s->s_z[SIZE_YY-1]; 
                v2.z = s->s_z[(SIZE_XX-1)*SIZE_YY] - s->s_z[SIZE_YY]; 

                s->vec_L_x[0] =  v2.z*v1.y - v2.y*v1.z;
                s->vec_L_y[0] = -v2.z*v1.x + v2.x*v1.z;
                s->vec_L_z[0] =  v2.y*v1.x - v2.x*v1.y;


            //right-up corner

                v1.z = s->s_z[0] - s->s_z[ SIZE_YY-2];
                v2.z = s->s_z[SIZE_XX*SIZE_YY-1] - s->s_z[2*SIZE_YY-1];

                s->vec_L_x[SIZE_YY-1] =  v2.z*v1.y - v2.y*v1.z;
                s->vec_L_y[SIZE_YY-1] = -v2.z*v1.x + v2.x*v1.z;
                s->vec_L_z[SIZE_YY-1] =  v2.y*v1.x - v2.x*v1.y;


        
            //left-down corner
                
                v1.z = s->s_z[(SIZE_XX-1)*SIZE_YY +  1] - s->s_z[SIZE_XX*SIZE_YY - 1];
                v2.z = s->s_z[(SIZE_XX-2)*SIZE_YY ] - s->s_z[0];

                s->vec_L_x[(SIZE_XX-1)*SIZE_YY] =  v2.z*v1.y - v2.y*v1.z;
                s->vec_L_y[(SIZE_XX-1)*SIZE_YY] = -v2.z*v1.x + v2.x*v1.z;
                s->vec_L_z[(SIZE_XX-1)*SIZE_YY] =  v2.y*v1.x - v2.x*v1.y;


            //right-down corner

                v1.z = s->s_z[(SIZE_XX-1)*SIZE_YY] - s->s_z[SIZE_XX*SIZE_YY - 2];
                v2.z = s->s_z[(SIZE_XX-1)*SIZE_YY -1] - s->s_z[SIZE_YY-1];

                s->vec_L_x[SIZE_XX*SIZE_YY-1] =  v2.z*v1.y - v2.y*v1.z;
                s->vec_L_y[SIZE_XX*SIZE_YY-1] = -v2.z*v1.x + v2.x*v1.z;
                s->vec_L_z[SIZE_XX*SIZE_YY-1] =  v2.y*v1.x - v2.x*v1.y;



            //unit vnormal vector
        for(int i=0; i<SIZE_XX*SIZE_YY; i++){

                    REAL modulo =  sqrt((s->vec_L_x[i]*s->vec_L_x[i]) + (s->vec_L_y[i]*s->vec_L_y[i]) + (s->vec_L_z[i]*s->vec_L_z[i]) );
                                    
                    s->vec_L_x[i] /= modulo;
                    s->vec_L_y[i] /= modulo;
                    s->vec_L_z[i] /= modulo;


    }
}

void calculate_incidence_angle(struct sea_surface *s, int SIZE_XX, int SIZE_YY, REAL deltaX, REAL deltaY, REAL init_point_XX, REAL init_point_YY){



        for(int i=0; i < SIZE_XX; i++){
            
            for(int j=0; j < SIZE_YY; j++){

                REAL producto_escalar = (((REAL)i*deltaX) + init_point_XX)*s->vec_L_x[i*SIZE_YY+j] +   // producto escalar  incidente con normal
                                        (((REAL)j*deltaY) + init_point_YY)*s->vec_L_y[i*SIZE_YY+j] +
                                        s->vec_inc_z[i*SIZE_YY+j]*s->vec_L_z[i*SIZE_YY+j];
            
           
                s->angle[i*SIZE_YY+j]=(-producto_escalar); 
            }

        } 
    }
void calculate_incidence_unit_vec(struct sea_surface *s, struct coor *tx_xyz, int SIZE_XX, int SIZE_YY, REAL deltaX, REAL deltaY, REAL init_point_XX, REAL init_point_YY){


        //CALCULATE VECTOR, DEST-ORIGIN
    
        for(int i=0; i < SIZE_XX; i++){
            
            for(int j=0; j < SIZE_YY; j++){

                s->vec_inc_x[i*SIZE_YY+j] = (((REAL)i*deltaX) + init_point_XX) - tx_xyz->x ;
           
           
                s->vec_inc_y[i*SIZE_YY+j] = (((REAL)j*deltaY) + init_point_YY) - tx_xyz->y ;
          
         
                s->vec_inc_z[i*SIZE_YY+j] = s->s_z[i*SIZE_YY+j] - tx_xyz->z;
          

                s->vec_module_inc[i*SIZE_YY+j]=   s->vec_inc_x[i*SIZE_YY+j]*s->vec_inc_x[i*SIZE_YY+j] +    //// modulo^2 incidente
                                                  s->vec_inc_y[i*SIZE_YY+j]*s->vec_inc_y[i*SIZE_YY+j] +
                                                  s->vec_inc_z[i*SIZE_YY+j]*s->vec_inc_z[i*SIZE_YY+j];

                s->vec_module_inc[i*SIZE_YY+j]=sqrt(s->vec_module_inc[i*SIZE_YY+j]);
         
                s->vec_inc_x[i*SIZE_YY+j] /= s->vec_module_inc[i*SIZE_YY+j];
           
                s->vec_inc_y[i*SIZE_YY+j] /= s->vec_module_inc[i*SIZE_YY+j];
          
                s->vec_inc_z[i*SIZE_YY+j] /= s->vec_module_inc[i*SIZE_YY+j];
            }
        }
    }
void calculate_reflected_unit_vec(struct sea_surface *s, struct coor *rx_xyz, int SIZE_XX, int SIZE_YY, REAL deltaX, REAL deltaY, REAL init_point_XX, REAL init_point_YY){
  
        //CALCULATE VECTOR, DEST-ORIGIN
 
        for(int i=0; i < SIZE_XX; i++){
            
            for(int j=0; j < SIZE_YY; j++){

                s->vec_ref_x[i*SIZE_YY+j] =  rx_xyz->x - (((REAL)i*deltaX) + init_point_XX);  
            
                s->vec_ref_y[i*SIZE_YY+j] =  rx_xyz->y - (((REAL)j*deltaY) + init_point_YY);    
         
                s->vec_ref_z[i*SIZE_YY+j] =  rx_xyz->z - s->s_z[i*SIZE_YY+j];     
          

                s->vec_module_ref[i*SIZE_YY+j]=   s->vec_ref_x[i*SIZE_YY+j]*s->vec_ref_x[i*SIZE_YY+j] +    //// modulo^2 refeljado
                                        s->vec_ref_y[i*SIZE_YY+j]*s->vec_ref_y[i*SIZE_YY+j] +
                                        s->vec_ref_z[i*SIZE_YY+j]*s->vec_ref_z[i*SIZE_YY+j];

                s->vec_module_ref[i*SIZE_YY+j]=sqrt(s->vec_module_ref[i*SIZE_YY+j]);
           
                s->vec_ref_x[i*SIZE_YY+j] /= s->vec_module_ref[i*SIZE_YY+j];
            
                s->vec_ref_y[i*SIZE_YY+j] /= s->vec_module_ref[i*SIZE_YY+j];
            
                s->vec_ref_z[i*SIZE_YY+j] /= s->vec_module_ref[i*SIZE_YY+j];

            }
            
        }
    }     
void calculate_module(struct sea_surface *s, bool inc, bool ref, int SIZE_XX, int SIZE_YY){

    if(inc)
        #pragma SIMD
        for(int i=0; i < SIZE_XX*SIZE_YY; i++){

            s->vec_module_inc[i]=   s->vec_inc_x[i]*s->vec_inc_x[i] +    //// modulo^2 incidente
                                    s->vec_inc_y[i]*s->vec_inc_y[i] +
                                    s->vec_inc_z[i]*s->vec_inc_z[i];

            s->vec_module_inc[i]=sqrt(s->vec_module_inc[i]);
        }
    
    if(ref)
       for(int i=0; i < SIZE_XX*SIZE_YY; i++){

            s->vec_module_ref[i]=   s->vec_ref_x[i]*s->vec_ref_x[i] +    //// modulo^2 refeljado
                                    s->vec_ref_y[i]*s->vec_ref_y[i] +
                                    s->vec_ref_z[i]*s->vec_ref_z[i];

            s->vec_module_ref[i]=sqrt(s->vec_module_ref[i]);
        }}



void calculate_ponderacio(struct sea_surface *s,  int SIZE_XX, int SIZE_YY){

        REAL K1_K2_x;
        REAL K1_K2_y;  
        REAL K1_K2_z;

         for(int i=0; i < SIZE_XX; i++){
            
            for(int j=0; j < SIZE_YY; j++){

            
                K1_K2_x = (s->vec_inc_x[i*SIZE_YY+j] - s->vec_ref_x[i*SIZE_YY+j]);
                K1_K2_y = (s->vec_inc_y[i*SIZE_YY+j] - s->vec_ref_y[i*SIZE_YY+j]);  
                K1_K2_z = (s->vec_inc_z[i*SIZE_YY+j] - s->vec_ref_z[i*SIZE_YY+j]);

                s->vec_ponderacio[i*SIZE_YY+j] =   -( K1_K2_x*s->vec_L_x[i*SIZE_YY+j] +
                                            K1_K2_y*s->vec_L_y[i*SIZE_YY+j] +
                                            K1_K2_z*s->vec_L_z[i*SIZE_YY+j]); 


                s->vec_ponderacio[i*SIZE_YY+j] /= sqrt(K1_K2_x*K1_K2_x + K1_K2_y*K1_K2_y + K1_K2_z*K1_K2_z);

                //LA PONDERACIO DONA NEGATIU, O ES CANVIA EL SIGNE AQUI, O EL VECTORL NORMAL HA DE SER NEGATIU

                //s->vec_ponderacio[i*SIZE_YY+j] = s->vec_ponderacio[i*SIZE_YY+j]; // Serni: Què fa aquesta línia de codi??  Mouad: Ara no fa res. Abans feia *2*pi/landa,
                                                                                                                // es va simplificar y treure fora de la integral
            }
        }
}

void calculate_propagation(struct sea_surface *s, int SIZE_XX, int SIZE_YY, REAL LANDA){
 
     for(int i=0; i < SIZE_XX; i++){
            
        for(int j=0; j < SIZE_YY; j++){


            REAL fase = (s->vec_module_inc[i*SIZE_YY+j])+(s->vec_module_ref[i*SIZE_YY+j]);
            REAL modulos = s->vec_module_inc[i*SIZE_YY+j]*s->vec_module_ref[i*SIZE_YY+j];
            
            fase = fase*2.0*M_PI/LANDA;

            s->E_field[i*SIZE_YY+j].real(cos(fase)/modulos);
            s->E_field[i*SIZE_YY+j].imag(sin(fase)/modulos);
          
        }
    }
}

void calculate_fresnell_coefficient_HV(std::complex<REAL> Ep11, std::complex<REAL> Ep22, std::complex<REAL> Mu11, std::complex<REAL> Mu22, REAL angle, std::complex<REAL> *RV, std::complex<REAL> *RH){
    
   
    std::complex<REAL> aux = (Ep11*Ep11)/(Ep22*Ep22);

    std::complex<REAL> eta_SQ_11 = (Mu11/Ep11); //eta^2
    std::complex<REAL> eta_SQ_22 = (Mu22/Ep22);

    std::complex<REAL> A = eta_SQ_11*aux;
    std::complex<REAL> B = eta_SQ_11*eta_SQ_11*aux/eta_SQ_22;

    std::complex<REAL> eta_11=sqrt(eta_SQ_11);
    std::complex<REAL> eta_22=sqrt(eta_SQ_22);



    REAL sin_SQ_angle = 1.0-angle*angle;
    REAL cos_angle = angle;
    

    *RV = (eta_11*cos_angle - sqrt(eta_SQ_22 - A*sin_SQ_angle))/(eta_11*cos_angle + sqrt(eta_SQ_22 - A*sin_SQ_angle));

    *RH = (eta_22*cos_angle - sqrt(eta_SQ_11 - B*sin_SQ_angle))/(eta_22*cos_angle + sqrt(eta_SQ_11 - B*sin_SQ_angle));

 
   /* cout<<*RH<<" "<<*RV<<endl;*/
}
void calculate_polarization(struct sea_surface *s, struct coor * polarization_inc, struct coor * polarization_ref, REAL LANDA, std::complex<REAL> eta_11, std::complex<REAL> eta_22, std::complex<REAL> eta_SQ_11, std::complex<REAL> eta_SQ_22, std::complex<REAL> A, std::complex<REAL> B, struct Vector_Fiel *E, int SIZE_XX, int SIZE_YY){

    struct coor H_inc_ref;  //Hinci Href son iguals
    struct coor V_inc;
   // struct coor H_ref;
    struct coor V_ref;
    std::complex<REAL> Pinc_Vinc_RV;
    std::complex<REAL> Pinc_Hinc_RH;
    
    std::complex<REAL> polarization_total;
     
    for(int i=0; i <  SIZE_XX*SIZE_YY; i++) {
        
        REAL sin_SQ_angle = 1.0-s->angle[i]*s->angle[i];
        REAL cos_angle = s->angle[i];
        
    
        Pinc_Vinc_RV = (eta_11*cos_angle - sqrt(eta_SQ_22 - A*sin_SQ_angle))/(eta_11*cos_angle + sqrt(eta_SQ_22 - A*sin_SQ_angle));
    
        Pinc_Hinc_RH = (eta_22*cos_angle - sqrt(eta_SQ_11 - B*sin_SQ_angle))/(eta_22*cos_angle + sqrt(eta_SQ_11 - B*sin_SQ_angle));

        ///calculo Hinc con n x K1
        H_inc_ref.x =   s->vec_L_y[i]*s->vec_inc_z[i] - s->vec_L_z[i]*s->vec_inc_y[i];    
        H_inc_ref.y = -(s->vec_L_x[i]*s->vec_inc_z[i] - s->vec_L_z[i]*s->vec_inc_x[i]);  
        H_inc_ref.z =   s->vec_L_x[i]*s->vec_inc_y[i] - s->vec_L_y[i]*s->vec_inc_x[i];    

       
        ///calculo Href con hinc x K2 
        // Serni: Aquest càlcul està equivocat perquè H_ref = H_ref
        //Mouad: He creat una sola varible, com son igual s'utilitza la mateixa 

      
        ///calculo Vinc con K1 x hinc
        V_inc.x =   s->vec_inc_y[i]*H_inc_ref.z - s->vec_inc_z[i]*H_inc_ref.y;    
        V_inc.y = -(s->vec_inc_x[i]*H_inc_ref.z - s->vec_inc_z[i]*H_inc_ref.x);  
        V_inc.z =   s->vec_inc_x[i]*H_inc_ref.y - s->vec_inc_y[i]*H_inc_ref.x;

        
 
        // El vector vref = Kref x href, on Kref és el direcció de reflexió especular del trosset de superfície.
        // https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector
        // Kref = K1 - 2*(K1·n)*n
        // Mouad: Aquest calcul es el que es fa per calcular l'angle incident, es podet reaprofitar el rpducte escalar
        REAL K1n = s->vec_inc_x[i]*s->vec_L_x[i] +
                   s->vec_inc_y[i]*s->vec_L_y[i] + 
                   s->vec_inc_z[i]*s->vec_L_z[i];   // es el producte escalar entre K1 i n

        struct coor K_ref;  // Cálcul de Kref
        K_ref.x = s->vec_inc_x[i] - 2*K1n*s->vec_L_x[i];
        K_ref.y = s->vec_inc_y[i] - 2*K1n*s->vec_L_y[i];
        K_ref.z = s->vec_inc_z[i] - 2*K1n*s->vec_L_z[i];
        
        // Serni: i ara podem calcular Vref = Kref x href
        V_ref.x =   K_ref.y*H_inc_ref.z - K_ref.z*H_inc_ref.y;
        V_ref.y = -(K_ref.x*H_inc_ref.z - K_ref.z*H_inc_ref.x);
        V_ref.z =   K_ref.x*H_inc_ref.y - K_ref.y*H_inc_ref.x;



        //prodcuto escalar --- projecting the [p^] onto the [h^inc] and [v^inc] vectors--- Multiplication between fresnell coef. and polarization vectors
        Pinc_Vinc_RV *= polarization_inc->x*V_inc.x + polarization_inc->y*V_inc.y + polarization_inc->z*V_inc.z;
        Pinc_Hinc_RH *= polarization_inc->x*H_inc_ref.x + polarization_inc->y*H_inc_ref.y + polarization_inc->z*H_inc_ref.z;


        // Serni: falta revisar codi des d'aquí  --------------------------------------------------------------------
       

        polarization_total =    (Pinc_Hinc_RH * H_inc_ref.x + Pinc_Vinc_RV * V_ref.x)*polarization_ref->x +
                                (Pinc_Hinc_RH * H_inc_ref.y + Pinc_Vinc_RV * V_ref.y)*polarization_ref->y +
                                (Pinc_Hinc_RH * H_inc_ref.z + Pinc_Vinc_RV * V_ref.z)*polarization_ref->z;


        s->E_field[i] = s->E_field[i]*polarization_total*s->vec_ponderacio[i]*2.0*M_PI/LANDA;
        //s->E_field[i] = s->vec_ponderacio[i]*2.0*M_PI/LANDA;



        // fins aquí. Cal plasmar-hi el que diu la nova versió del document de requeriments
        //--------------------------------------------------------------------

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



//MAIN FUNCTION, GROUPS FUNCTIONS ABOVE
void compute_surface_reflection(int SIZE_XX, int SIZE_YY, REAL Dx, REAL Dy, REAL LANDA, REAL Amp,  REAL Temp, REAL Salinity, struct coor *sat_coor_transmitter, struct coor *sat_coor_receiver, struct coor *polarization_inc, struct coor *polarization_ref, REAL *surface_s_z){



    REAL FREQ = 299792458.0/LANDA;
    REAL deltaX = (REAL)Dx/SIZE_XX;
    REAL deltaY = (REAL)Dy/SIZE_YY;

    REAL init_point_XX = -((REAL)Dx/2.0);
    REAL init_point_YY = -((REAL)Dy/2.0);

    struct sea_surface surface;
    surface.s_z = surface_s_z;


    struct Vector_Fiel E_sum;

    allocate_memory(&surface, SIZE_XX*SIZE_YY);
        


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
            
                                       
                    
    init_flat_surface(&surface, SIZE_XX, SIZE_YY);  ///???????

    calculate_incidence_unit_vec(&surface, sat_coor_transmitter, SIZE_XX, SIZE_YY, deltaX, deltaY, init_point_XX, init_point_YY);  

    calculate_incidence_angle(&surface, SIZE_XX, SIZE_YY,  deltaX, deltaY, init_point_XX, init_point_YY);  

    calculate_reflected_unit_vec(&surface, sat_coor_receiver, SIZE_XX, SIZE_YY,  deltaX, deltaY, init_point_XX, init_point_YY); 

    calculate_ponderacio(&surface, SIZE_XX, SIZE_YY);  

    calculate_propagation(&surface, SIZE_XX, SIZE_YY, LANDA);  

    calculate_polarization(&surface, polarization_inc, polarization_ref, LANDA, eta_11, eta_22, eta_SQ_11, eta_SQ_22, A, B, &E_sum, SIZE_XX, SIZE_YY);
                       
          
      
    std::cout<<Amp*E_sum.x.real()/(2.0*LANDA)<<" "<<Amp*E_sum.y.real()/(2.0*LANDA)<<" "<<Amp*E_sum.z.real()/(2.0*LANDA)<<" "<<std::endl;

    free_memory(&surface, SIZE_XX*SIZE_YY);

}


//MAIN FUNCTION WITH BARREL OVER PHI AND THETA, GROUPS FUNCTIONS ABOVE
void compute_surface_reflection_Receiver_barrel(int SIZE_XX, int SIZE_YY, REAL Dx, REAL Dy, REAL LANDA, REAL Amp,  REAL Temp, REAL Salinity, struct coor *sat_coor_transmitter, struct coor *sat_coor_receiver, struct coor *polarization_inc, struct coor *polarization_ref, REAL *surface_s_z, REAL *barrel_result, int azimut, int inclination, REAL RX_RADIO){



    REAL FREQ = 299792458.0/LANDA;
    REAL deltaX = (REAL)Dx/SIZE_XX;
    REAL deltaY = (REAL)Dy/SIZE_YY;

    REAL init_point_XX = -((REAL)Dx/2.0);
    REAL init_point_YY = -((REAL)Dy/2.0);

    struct sea_surface surface;
    surface.s_z = surface_s_z;


    struct Vector_Fiel E_sum;

    allocate_memory(&surface, SIZE_XX*SIZE_YY);
        


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
            
                                       
   // init_flat_surface(&surface, SIZE_XX, SIZE_YY);  ///??????? ESTO BORRARIA TODO EL INPUT RECIBIDO
    init_normal_vector_surface(&surface, SIZE_XX, SIZE_YY, deltaX, deltaY);

    calculate_incidence_unit_vec(&surface, sat_coor_transmitter, SIZE_XX, SIZE_YY, deltaX, deltaY, init_point_XX, init_point_YY);  

    calculate_incidence_angle(&surface, SIZE_XX, SIZE_YY,  deltaX, deltaY, init_point_XX, init_point_YY);  


    

    for(int i =0; i< inclination; i++){
      
        for(int j =0; j< azimut; j++){


            //CHECK!!!!
            sat_coor_receiver->x =  (RX_RADIO * sin(((REAL)i / inclination) * M_PI/2.0)*cos(((REAL)j / azimut) * 2.0*M_PI));
            sat_coor_receiver->y =  (RX_RADIO * sin(((REAL)i / inclination) * M_PI/2.0)*sin(((REAL)j / azimut) * 2.0*M_PI));
            sat_coor_receiver->z =  (RX_RADIO * cos(((REAL)i / inclination) * M_PI/2.0));

            init_polarization(polarization_ref, sin(((REAL)j / azimut) * 2.0*M_PI), -cos(((REAL)j / azimut) * 2.0*M_PI), 0);

         //   std::cout<< sat_coor_receiver->x<<", "<< sat_coor_receiver->y<<", "<< sat_coor_receiver->z<<"    "<<sqrt(( sat_coor_receiver->x* sat_coor_receiver->x)+( sat_coor_receiver->y* sat_coor_receiver->y)+( sat_coor_receiver->z* sat_coor_receiver->z));
            // init_Sat_ubi(&sat_coor_receiver,  , , 10);  //done!


            calculate_reflected_unit_vec(&surface, sat_coor_receiver, SIZE_XX, SIZE_YY,  deltaX, deltaY, init_point_XX, init_point_YY); 

            calculate_ponderacio(&surface, SIZE_XX, SIZE_YY);  

            calculate_propagation(&surface, SIZE_XX, SIZE_YY, LANDA);  

            calculate_polarization(&surface, polarization_inc, polarization_ref, LANDA, eta_11, eta_22, eta_SQ_11, eta_SQ_22, A, B, &E_sum, SIZE_XX, SIZE_YY);
            
            for(int p=0; p <  SIZE_XX*SIZE_YY; p++)
                barrel_result[i*azimut + j] += surface.E_field[p].real();//surface.vec_ponderacio[p];// (E_sum.x.real()*E_sum.x.real()+E_sum.y.real()*E_sum.y.real()+E_sum.z.real()*E_sum.z.real());
            
         //   std::cout<<" out: "<<  barrel_result[i*azimut + j]<<std::endl;
    
        }
        std::cout<< std::endl;
    }



    free_memory(&surface, SIZE_XX*SIZE_YY);

    //std::cout<<Amp*E_sum.x.real()/(2.0*LANDA)<<" "<<Amp*E_sum.y.real()/(2.0*LANDA)<<" "<<Amp*E_sum.z.real()/(2.0*LANDA)<<" "<<std::endl;


}
