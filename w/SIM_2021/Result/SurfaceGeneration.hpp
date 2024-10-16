#pragma once



int  elfouhaily_spectrum(double U10,double theta,double Omega,  double *PSI,  int Nx,int Ny, double deltakx, double deltaky);
void Generate_Random_fase(double *PSI_COMPLEX, double *PSI, int Nx, int Ny);
void matrix_transform(double * PSI, int Nx, int Ny);
void matrix_2d_FFT(double *PSI_COMPLEX, int Nx, int Ny);
void Complex_to_Real_Matrix(double *PSI_COMPLEX, double *PSI, int Nx, int Ny);
void  generate_surface(double * PSI, int Nx, int Ny, double wind, double theta, double Omega, double Rand_mod, int Dx, int Dy);

void check_MSS_Derivative(double *PSI_COMPLEX, int Nx, int Ny, double deltaX, double deltaY);
void check_MSS_surface(double *PSI_COMPLEX, int Nx, int Ny);
void check_MSS_spectrum(double * PSI, int Nx, int Ny,double deltakx,double deltaky);

