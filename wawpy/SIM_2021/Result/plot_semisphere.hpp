#pragma once


static int winW = 512;
static int winH = 512;



static float rotX = 4.0;
static float dX = 1.0;
static float updateRate = 30;




//COMING FROM OTHER CODE
static int azimut_max;
static int inc_max;
static double max_Efield;
static double* surface_barrel_output;

struct RGB {

    float R;  //max value 1.0
    float G;
    float B;
};




void reshape(int w, int h);



void timer(int);




void E_to_RGB(double * E_field, double max_Efield, RGB * color);


void plotsphere();
void drawTriangle(float* v1, float* v2, float* v3, RGB color);
// Main Function
void plot_data_semisphere(int az, int inc, double maxfield,double minfield, double* s, struct coor *tx_xyz);



void plotsphere_sample();
void drawTriangle_sample(float* v1, float* v2, float* v3, float low_high);
// Main Function sample
void plot_data_semisphere_sample(int az, int inc);
