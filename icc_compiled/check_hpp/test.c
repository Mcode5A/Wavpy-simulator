//#include <windows.h> // Must have for Windows platform builds
//#include <GL\glew.h>
//#include <GL\freeglut.h>
#include <iostream>
//#include <gl\gl.h> // Microsoft OpenGL headers (version 1.1 by themselves)
//#include <gl\glu.h> // OpenGL Utilities



#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <GL/glut.h>




using namespace std;


 
int winW = 512;
int winH = 512;
// These are the 12 vertices for the icosahedron
 
 

 
static GLfloat rotX = 4.0;
static GLfloat dX = 1.0;
static GLfloat updateRate = 30;

//Helpers
 
float azimut_max = 36;
float inc_max = 9;

void drawTriangle_sample(float* v1, float* v2, float* v3, float inc, float azimut, float low_high)
{


    low_high = 0.375;
    // Draw a triangle with specified vertices 
    glBegin(GL_TRIANGLES);
    //  glColor3f(0.0,0.0,0.8);
    glColor3f(inc, azimut, low_high);
    //  glNormal3fv(v1); 
    //  glNormal3fv(v2); 
    //  glNormal3fv(v3);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glColor3f(inc, azimut, low_high);
    //  glNormal3fv(v1); 
    //  glNormal3fv(v2); 
    //  glNormal3fv(v3);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
    glEnd();
}



void plotsphere_sample() {
    
    
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


    glRotatef(70, 2, 2.0, 0.0);


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

    float two_pi_sweep = 2.0 * 3.14159265359;
    float half_pi_sweep = 3.14159265359 / 2.0;


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




        for (int i = 1; i < azimut_max; i++)
        {
            float x_1 = radio * cos(two_pi_sweep * ((float)i / azimut_max));
            float z_1 = radio * sin(two_pi_sweep * ((float)i / azimut_max));

            second_low_point[0] = x_1 * sin_theta_low;
            second_low_point[1] = radio * cos_theta_low;
            second_low_point[2] = z_1 * sin_theta_low;


            drawTriangle_sample(first_low_point, second_low_point, first_top_point, (float)(j / inc_max), (float)((i - 1) / azimut_max), 1.0);

            second_top_point[0] = x_1 * sin_theta_high;
            second_top_point[1] = radio * cos_theta_high;
            second_top_point[2] = z_1 * sin_theta_high;


            drawTriangle_sample(second_top_point, second_low_point, first_top_point, (float)(j / inc_max), (float)((i - 1) / azimut_max), 0.0);

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


        drawTriangle_sample(first_low_point, second_low_point, first_top_point, (float)(j / inc_max), (float)((azimut_max - 1) / azimut_max), 1.0);
        drawTriangle_sample(second_top_point, second_low_point, first_top_point, (float)(j / inc_max), (float)((azimut_max - 1) / azimut_max), 0.0);

    }






    glPopMatrix();
    glutSwapBuffers(); // If double buffering

}
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


 

//Timer function
void timer(int)
{
    glutPostRedisplay();
    glutTimerFunc(1000.0 / updateRate, timer, 0);
}





// Main Function
int main(int argc, char** argv)
{

    int zero = 0;
    char** empty = new char* [1];

    azimut_max = 36;//az;
    inc_max = 9;//inc;
  //  max_Efield = maxfield;
 //   surface_barrel_output = s; // new  float[inc_max * azimut_max];




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
