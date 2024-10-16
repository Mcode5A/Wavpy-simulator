#pragma once



int  elfouhaily_spectrum(double U10,double theta,double Omega,  double *PSI,  int Nx,int Ny, double deltakx, double deltak);
double  elfouhaily_spectrum_iter_calc(int i, int j, double U10,double theta,double Omega,   int Nx,int Ny, double deltakx, double deltaky);
void Generate_Random_fase(double *PSI_COMPLEX, double *PSI, int Nx, int Ny);
void matrix_transform(double * PSI, int Nx, int Ny);
void matrix_2d_FFT(double *PSI_COMPLEX, int Nx, int Ny);
void Complex_to_Real_Matrix(double *PSI_COMPLEX, double *PSI, int Nx, int Ny);
void generate_surface(double * PSI, int Nx, int Ny, double wind, double theta, double Omega, double Rand_mod, int Dx, int Dy, double deltaX, double deltaY);
void generate_surface_merge(double * PSI, int Nx, int Ny, int start_XX_rank, double wind, double theta, double Omega, double Rand_mod, int Dx, int Dy, double deltaX, double deltaY, int size_procs);

void  generate_surface_merge_mpi(double * PSI_COMPLEX, int Nx, int Ny, int start_XX_rank, double wind, double theta, double Omega, double Rand_mod, int Dx, int Dy, double deltaX, double deltaY, double randseed,  int size_procs);
void check_MSS_Derivative(double *PSI_COMPLEX, int Nx, int Ny, double deltaX, double deltaY);
void check_MSS_surface(double *PSI_COMPLEX, int Nx, int Ny);
void check_MSS_spectrum(double * PSI, int Nx, int Ny,double deltakx,double deltaky);

