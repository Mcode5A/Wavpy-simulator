#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define REAL double
double const pi = 3.14159265358979323846;

 
void init_fact(REAL * vec, int n){

    vec[0]=1.0;
    
    for(int i=1; i<n; i++)
      vec[i]= (REAL)(vec[i-1]*(REAL)i);
 
      
     for(int i=2; i<n; i+=4){
      vec[i]= -vec[i];
      vec[i+1]= -vec[i+1];
      }
}

 

  
REAL cos_new(REAL x, REAL *Fact_inv_sign, int n)
{
    //if number above 2pi
    REAL aux =  x / (2.0*pi);
    int int_part =(int)aux;
    x = (aux-int_part)*(2.0*pi);


    REAL Sum       = 1.0;
    REAL potencias_valor  = 1.0;
    REAL x2        = x * x;

       
    
    for (unsigned int i=2; i<n; i+=2)
    {
        potencias_valor   *= x2;
        Sum     += potencias_valor /Fact_inv_sign[i];
    }
    return Sum;
    
    
}
 
 
 
int main(){


    REAL val =45.0*pi/180.0;
    int N =1;
    int n = 1000; //par
    
    REAL res1 =0.0;
    REAL res2 =0.0;
   
    REAL * vec = (REAL *)malloc(N*sizeof(REAL));
    
    REAL * Fact_sign = (REAL *)malloc(n*sizeof(REAL));
    
    
     init_fact(Fact_sign, n);
    
    
    for(int i=0 ; i<N; i++){
    
        res1 +=cos(val);
        res2 +=cos_new(val, Fact_sign, n);
       
    
    }



printf("%0.50f, %0.50f, %0.50f\n", res1, res2, res1-res2);


return 0;}



 
 