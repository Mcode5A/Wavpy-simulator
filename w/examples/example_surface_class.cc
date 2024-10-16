#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include "reflecting_surface.h"
#include "ancillary_functions.h"

using namespace std;

//============================ PROGRAM =====================================
int main (int argc, char* argv[]) {
	Reflecting_surface mysurf;
	double rco[2], rcross[2], epsilon_air[2], kx_ky_spectrum[3];
	epsilon_air[0] = 1.0;
	epsilon_air[1] = 0.0;
	mysurf.dump_parameters();
	mysurf.c21_coeff = 1.0;
	mysurf.c03_coeff = 2.0;
	mysurf.wind_U10_speed = 14.0;
	mysurf.wind_U10_azimuth = 4.0;
	mysurf.epsilon_sea_ice(30.0);
	mysurf.dump_parameters();
	mysurf.compute_Rfresnel_circular( 30.0, epsilon_air, rco, rcross );
	printf("Fresnel CIRC: %f %f %f %f \n", rco[0], rco[1], rcross[0], rcross[1]);
	mysurf.compute_sea_spectrum(4, 100.5, 90.0, 1.0);
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			mysurf.get_surf_spectrum(i, j, kx_ky_spectrum);
			printf("SPEC %f %f %f\n", kx_ky_spectrum[0], kx_ky_spectrum[1], kx_ky_spectrum[2]);
		}
	}
	return 0;
}
