/*#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/gl.h>  
#include <GL/glu.h>  */


 
#include <iostream>
#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <GL/glut.h>

//#include "plot_semisphere.hpp"

struct coor{
    
    double x=0.0;
    double y=0.0;
    double z=0.0;

    };
static int winW = 512;
static int winH = 512;



static float rotX = 4.0;
static float dX = 1.0;
static float updateRate = 50;


static struct coor tx_coor; 
static double* surface_barrel_output;


//COMING FROM OTHER CODE
static int azimut_max;
static int inc_max;
static double max_Efield;
static double min_Efield;

struct RGB {

    float R;  //max value 1.0
    float G;
    float B;
};




void reshape(int w, int h)
{
    winW = w;
    winH = h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, w, 0.0, h, -w, w);
    glMatrixMode(GL_MODELVIEW);
}

void timer(int)
{
    glutPostRedisplay();
    glutTimerFunc(1000.0 / updateRate, timer, 0);
}
 
void E_to_RGB(double * E_field, double max_Efield, double min_Efield, RGB * color) {

 double ratio = (E_field[0]-min_Efield) / (max_Efield-min_Efield); 
    double inv_ratio = 1.0 - ratio;

    color->G = ratio * 1.00000f;
    color->R = 1.0f;
    color->B = inv_ratio * 1.00000f;
    
    
    if (E_field[0]==max_Efield)
      color->G = 0.0f;
};


void drawTriangle(float* v1, float* v2, float* v3, RGB color)
{

    // Draw a triangle with specified vertices 
    glBegin(GL_TRIANGLES);
    //  glColor3f(0.0,0.0,0.8);
    glColor3f(color.R, color.G, color.B);
    //  glNormal3fv(v1); 
    //  glNormal3fv(v2); 
    //  glNormal3fv(v3);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glColor3f(color.R, color.G, color.B);
    //  glNormal3fv(v1); 
    //  glNormal3fv(v2); 
    //  glNormal3fv(v3);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
    glEnd();
}


void plotsphere()
{
    srand(5);

    //Enable Antialiasing
    glEnable(GL_LINE_SMOOTH);

    glClear(GL_COLOR_BUFFER_BIT);
    glClear(GL_DEPTH_BUFFER_BIT);


    glEnable(GL_DEPTH_TEST);

    glLoadIdentity();
    // Set the viewing transformation
    gluLookAt(0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    glPushMatrix();
    glTranslatef(250, 250, 250);
    glScalef(250, 250, 250);


    glRotatef(70, 2, 2.0, 0.0);


    glRotatef(rotX, 0.0, 1.0, 0.0);

    rotX += dX;




    float ptr1[3];
    float ptr2[3];
    float ptr3[3];
    float ptr4[3];

    float* first_low_point = ptr1;
    float* second_low_point = ptr2;
    float* first_top_point = ptr3;
    float* second_top_point = ptr4;




    float radio = 1.0f;

    float two_pi_sweep = 2.0f * 3.14159265359f;
    float half_pi_sweep = 3.14159265359f / 2.0f;

    struct RGB color;

    for (int j = 0; j < inc_max; j++)
    {


        float theta_low = half_pi_sweep * (float)j / inc_max;
        float theta_high = half_pi_sweep * (float)(j + 1.0) / inc_max;


        float sin_theta_low = sin(theta_low);
        float sin_theta_high = sin(theta_high);
        float cos_theta_low = cos(theta_low);
        float cos_theta_high = cos(theta_high);


        //INIT POINT --> end point
        float x_0 = radio * cos(((float)0.0 / azimut_max) * two_pi_sweep);
        float y_0 = radio * cos(theta_low);
        float z_0 = radio * sin(((float)0.0 / azimut_max) * two_pi_sweep);

        first_low_point[0] = x_0 * sin_theta_low;
        first_low_point[1] = y_0;
        first_low_point[2] = z_0 * sin_theta_low;

        first_top_point[0] = x_0 * sin_theta_high;
        first_top_point[1] = radio * cos_theta_high;
        first_top_point[2] = z_0 * sin_theta_high;


        E_to_RGB(&surface_barrel_output[j * azimut_max], max_Efield, min_Efield, &color);


        for (int i = 1; i < azimut_max; i++)
        {
            float x_1 = radio * cos(two_pi_sweep * ((float)i / azimut_max));
            float z_1 = radio * sin(two_pi_sweep * ((float)i / azimut_max));

            second_low_point[0] = x_1 * sin_theta_low;
            second_low_point[1] = radio * cos_theta_low;
            second_low_point[2] = z_1 * sin_theta_low;



            drawTriangle(first_low_point, second_low_point, first_top_point, color);

            second_top_point[0] = x_1 * sin_theta_high;
            second_top_point[1] = radio * cos_theta_high;
            second_top_point[2] = z_1 * sin_theta_high;


            drawTriangle(second_top_point, second_low_point, first_top_point, color);

            float* aux = first_low_point;
            first_low_point = second_low_point;
            second_low_point = aux;

            aux = first_top_point;
            first_top_point = second_top_point;
            second_top_point = aux;


           E_to_RGB(&surface_barrel_output[j * azimut_max + i], max_Efield, min_Efield, &color);

        }


        second_top_point[0] = x_0 * sin_theta_high;
        second_top_point[1] = radio * cos_theta_high;
        second_top_point[2] = z_0 * sin_theta_high;

        second_low_point[0] = x_0 * sin_theta_low;
        second_low_point[1] = radio * cos_theta_low;
        second_low_point[2] = z_0 * sin_theta_low;
     

        drawTriangle(first_low_point, second_low_point, first_top_point, color);
        drawTriangle(second_top_point, second_low_point, first_top_point, color);

    }

    float d = 0.1;
    float r = 1.1;
    first_low_point[0] = r*(tx_coor.x);
    first_low_point[1] = r*(tx_coor.y);
    first_low_point[2] = r*(tx_coor.z);
    second_low_point[0] = r*(tx_coor.x);
    second_low_point[1] = r*(tx_coor.y);
    second_low_point[2] = r*(tx_coor.z+d);
    first_top_point[0] = r*(tx_coor.x);
    first_top_point[1] = r*(tx_coor.y+d);
    first_top_point[2] = r*(tx_coor.z);

    /*second_top_point[0] = r*(tx_coor.x);
    second_top_point[1] = r*(tx_coor.y+d);
    second_top_point[2] = r*(tx_coor.z+d);*/

    color.R=0.0;
    color.G=1.0;
    color.B=0.0;
    drawTriangle(first_low_point, second_low_point, first_top_point, color);
  //  drawTriangle(second_top_point, second_low_point, first_top_point, color);

    first_low_point[0] = -r*(tx_coor.x);
    first_low_point[1] = r*(tx_coor.y);
    first_low_point[2] = -r*(tx_coor.z);
    second_low_point[0] = -r*(tx_coor.x);
    second_low_point[1] = r*(tx_coor.y);
    second_low_point[2] = -r*(tx_coor.z+d);
    first_top_point[0] = -r*(tx_coor.x);
    first_top_point[1] = r*(tx_coor.y+d);
    first_top_point[2] = -r*(tx_coor.z);

    /*second_top_point[0] = r*(tx_coor.x);
    second_top_point[1] = r*(tx_coor.y+d);
    second_top_point[2] = r*(tx_coor.z+d);*/

    color.R=0.0;
    color.G=0.0;
    color.B=1.0;
    drawTriangle(first_low_point, second_low_point, first_top_point, color);
  //  drawTriangle(second_top_point, second_low_point, first_top_point, color);












    glPopMatrix();
    glutSwapBuffers(); // If double buffering
}


void drawTriangle_sample(float* v1, float* v2, float* v3, float inc, float azimut, float low_high)
{


    low_high = 0.375;
    // Draw a triangle with specified vertices 
    glBegin(GL_TRIANGLES);
    //  glColor3f(0.0,0.0,0.8);
    glColor3f(inc, azimut, low_high );
    //  glNormal3fv(v1); 
    //  glNormal3fv(v2); 
    //  glNormal3fv(v3);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glColor3f(inc,azimut,low_high);
    //  glNormal3fv(v1); 
    //  glNormal3fv(v2); 
    //  glNormal3fv(v3);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
    glEnd();
}

void plotsphere_sample()
{
    srand(5);

    //Enable Antialiasing
    glEnable(GL_LINE_SMOOTH);
  
    glClear(GL_COLOR_BUFFER_BIT);
    glClear(GL_DEPTH_BUFFER_BIT);


    //Enable depth buffer
    glEnable(GL_DEPTH_TEST);

    glLoadIdentity();
    // Set the viewing transformation
    gluLookAt(0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    glPushMatrix();
    glTranslatef(250, 250, 250);
    glScalef(250, 250, 250);
    

    glRotatef(70,2, 2.0, 0.0);

    
    glRotatef(rotX, 1.0, 1.0, 0.0);

    rotX += dX;




    float ptr1[3];
    float ptr2[3];
    float ptr3[3];
    float ptr4[3];

    float* first_low_point = ptr1;
    float* second_low_point = ptr2;
    float* first_top_point = ptr3;
    float* second_top_point = ptr4;




    float radio = 1.0f;
  
    float two_pi_sweep =  2.0*3.14159265359;
    float half_pi_sweep = 3.14159265359/2.0 ;

   
    for (int j = 0; j < inc_max; j++)
    {


        float theta_low = half_pi_sweep * (float)j/ inc_max;
        float theta_high = half_pi_sweep * (float)(j+1.0) / inc_max;


        float sin_theta_low = sin(theta_low);
        float sin_theta_high = sin(theta_high);
        float cos_theta_low = cos(theta_low);
        float cos_theta_high = cos(theta_high);


        //INIT POINT --> end point
        float x_0 = radio * cos(((float)0.0 / azimut_max) * two_pi_sweep);
        float y_0 = radio * cos(theta_low);
        float z_0 = radio * sin(((float)0.0 / azimut_max) * two_pi_sweep);

        first_low_point[0] = x_0* sin_theta_low;
        first_low_point[1] = y_0;
        first_low_point[2] = z_0* sin_theta_low;

        first_top_point[0] = x_0* sin_theta_high;
        first_top_point[1] = radio * cos_theta_high;
        first_top_point[2] = z_0 * sin_theta_high;

  


        for (int i = 1; i < azimut_max; i++)
        {
            float x_1 = radio * cos(two_pi_sweep*((float)i / azimut_max) );
            float z_1 = radio * sin(two_pi_sweep * ((float)i / azimut_max));

            second_low_point[0] = x_1 * sin_theta_low;
            second_low_point[1] = radio * cos_theta_low;
            second_low_point[2] = z_1 * sin_theta_low;


            drawTriangle_sample(first_low_point, second_low_point, first_top_point, ((float)j/ inc_max), ((float)(i - 1)/ azimut_max), 1.0);
         
            second_top_point[0] = x_1 * sin_theta_high;
            second_top_point[1] = radio * cos_theta_high;
            second_top_point[2] = z_1 * sin_theta_high;
        

            drawTriangle_sample(second_top_point, second_low_point, first_top_point, ((float)j / inc_max), ((float)(i - 1) / azimut_max), 0.0);
          
            float* aux = first_low_point;
            first_low_point = second_low_point;
            second_low_point = aux;

            aux = first_top_point;
            first_top_point = second_top_point;
            second_top_point = aux;



        }


        second_top_point[0] = x_0 * sin_theta_high;
        second_top_point[1] = radio * cos_theta_high;
        second_top_point[2] = z_0 * sin_theta_high;

        second_low_point[0] = x_0 * sin_theta_low;
        second_low_point[1] = radio * cos_theta_low;
        second_low_point[2] = z_0 * sin_theta_low;


       drawTriangle_sample(first_low_point, second_low_point, first_top_point, ((float)j / inc_max), ((float)(azimut_max- 1) / azimut_max), 1.0);
       drawTriangle_sample(second_top_point, second_low_point, first_top_point, ((float)j / inc_max), ((float)(azimut_max - 1) / azimut_max), 0.0);

    }
    





    glPopMatrix();
    glutSwapBuffers(); // If double buffering
}



 

// Main Function
void plot_data_semisphere(int az, int inc, double maxfield, double minfield, double* s, struct coor *tx_xyz)
{

    int zero = 0;
    char** empty = new char* [1];
  
    float modulo_coor =  sqrt((tx_xyz->x*tx_xyz->x)+(tx_xyz->y*tx_xyz->y)+(tx_xyz->z*tx_xyz->z));

    tx_coor.x = tx_xyz->x / modulo_coor;
    tx_coor.y = tx_xyz->z / modulo_coor;
    tx_coor.z = tx_xyz->y / modulo_coor;

    azimut_max = az;
    inc_max = inc;
    max_Efield = maxfield;
    min_Efield = minfield;
    surface_barrel_output = s; // new  float[inc_max * azimut_max];
   



    // Initialize glut  and create your window here
    glutInit(&zero, empty);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000, 600);
    glutInitWindowPosition(200, 200);
    glutCreateWindow("myplot");
    // init();
     // Set your glut callbacks here
    glutDisplayFunc(plotsphere);
    glutReshapeFunc(reshape);
    glutTimerFunc(1000.0 / updateRate, timer, 0);
    // Enter the glut main loop here
    glutMainLoop();


}



void plot_data_semisphere_sample(int az, int inc)
{

    int zero = 0;
    char** empty = new char* [1];

    azimut_max = 36;//az;
    inc_max = 9;//inc;





    // Initialize glut  and create your window here
    glutInit(&zero, empty);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(200, 100);
    glutCreateWindow("myplot");
    // init();
     // Set your glut callbacks here
    glutDisplayFunc(plotsphere_sample);
    glutReshapeFunc(reshape);
    glutTimerFunc(1000.0 / updateRate, timer, 0);
    // Enter the glut main loop here
    glutMainLoop();


}