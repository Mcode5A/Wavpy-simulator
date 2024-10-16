#include "../matplotlib-cpp-master/matplotlibcpp.h"
#include <vector>
#include <iostream>



void plot(int Nx,int Ny, double *PSI){
 
  std::vector<std::vector<double>> x, y, z;
  for (int i = 0; i < Nx;  i++) {
      std::vector<double> x_row, y_row, z_row;
      for (int j = 0; j < Ny; j++) {
          x_row.push_back(i);
          y_row.push_back(j);
          z_row.push_back(PSI[i*Ny+j]);
        
      }
      
      x.push_back(x_row);
      y.push_back(y_row);
      z.push_back(z_row);
  }
  
  plt::plot_surface(x, y, z);
  //plt::show();
  plt::title("PSI ");
   const char* filename = "./k.png";
  std::cout << "Saving result to " << filename << std::endl;;
  plt::save(filename);

}



 
void plot_complex(int Nx, int Ny, double *PSI, char *name,char *title, double wind, double Rand_mod, double deltakx, double deltaky){

  
  std::vector<std::vector<double>> x, y, z;
  for (int i = 0; i < Nx;  i++) {
      std::vector<double> x_row, y_row, z_row;
      for (int j = 0; j < Ny; j++) {
          x_row.push_back(2.0*M_PI*i/(deltakx*Nx));
          y_row.push_back(2.0*M_PI*j/(deltaky*Ny));
          z_row.push_back(PSI[i*2*Ny+2*j]);
          
      }
      
      x.push_back(x_row);
      y.push_back(y_row);
      z.push_back(z_row);
  }
  
  plt::plot_surface(x, y, z);
  plt::title(title);
  char* filename = (char *) malloc(20*sizeof(char));
  strcat(filename, "./"); 
  strcpy(filename, name);
  strcat(filename, ".png"); 
  std::cout << "Saving result to " << filename << std::endl;;
  plt::save(filename);


}


