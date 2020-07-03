#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <map>
#include <time.h> 
#include <sys/time.h> 
#include <stdlib.h>
#include <string.h>
#include <GL/gl.h>
#include <GL/glu.h> /* fichiers d?entetes OpenGL, GLU et GLUT */
#include <GL/freeglut.h>
#include <vector>
#include <list>
#include <limits>

typedef double R;
const R Pi=3.14159265;
const GLint smooth=32;

GLint cH;
GLint cF;

GLint nQ;
GLint nI;

GLdouble theta_eye;
GLdouble phi_eye;
GLdouble rayon;
GLdouble x_eye;
GLdouble y_eye;
GLdouble z_eye;
GLdouble ee;
GLint x_g;
GLint x_d;
GLint y_g;
GLint y_d;
GLint click_g;
GLint click_d;

GLint nS;
GLdouble pS;
GLdouble rS;

GLdouble pI;

R    * LIST_R;
R    * LIST_X;
R    * LIST_Y;
R    * LIST_Z;

using namespace std;

inline R CPUtime(){
#ifdef SYSTIMES
  struct tms buf;
  if (times(&buf)!=-1)
    return ((R)buf.tms_utime+(R)buf.tms_stime)/(long) sysconf(_SC_CLK_TCK);
  else
#endif
    return ((R) clock())/CLOCKS_PER_SEC;
}

template<typename T>
std::string to_string( const T & Value )
{
    // utiliser un flux de sortie pour créer la chaîne
    ostringstream oss;
    // écrire la valeur dans le flux
    oss << Value;
    // renvoyer une string
    return oss.str();
}

template <class T>
void from_string(T& t, 
                 const string& s, 
                 ios_base& (*f)(ios_base&))
{
  istringstream iss(s);
  iss >> f >> t;
}

void setup_illumination(){
  // Intialise and set lighting parameters
  GLfloat light_pos[] = {0.,0.,0., 1.};
  GLfloat light_ka[] = {0., 0., 0.2, 1.};
  GLfloat light_kd[] = {1., 1., 1., 1.};
  GLfloat light_ks[] = {1., 1., 1., 1.};
  //cH/(2*cF)
  glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ka);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_kd);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_ks);

  
  // Initialise and set material parameters
  GLfloat material_ka[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat material_kd[] = {0.43, 0.47, 0.54, 1.0};
  GLfloat material_ks[] = {0.33, 0.33, 0.52, 1.0};
  GLfloat material_ke[] = {0.0, 0.0, 0.0, 0.0};
  GLfloat material_se[] = {10.0};

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  material_ka);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  material_kd);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  material_ks);
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,  material_ke);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, material_se);
  
// Smooth shading
glShadeModel(GL_SMOOTH);

//Enable lighting
glEnable (GL_LIGHTING);
glEnable (GL_LIGHT0);

}

static void InitGL(void)
{
glClearColor(0,0,0,0);
//enable z buffering
glEnable(GL_DEPTH_TEST);

/*definition de la couleur utilisee */
/*pour effacer */
//glColor3f(1.0,1.0,1.0);
/*couleur courante*/

theta_eye=0;
phi_eye=0;
x_eye = 0;
y_eye = 0;
z_eye = 2*cH/cF; 

x_g=280;
y_g=280;
x_d=280;
y_d=280;
rayon=2*cH/cF;
}

static void afficher()
{
  
//glClearDepth(1.);
//glClearColor(0.0,0.0,0.0,0.0);
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
//cout<<x_eye<<", "<<y_eye<<", "<<z_eye<<endl;

gluLookAt(x_eye,y_eye,z_eye,0.,0.,0.,0,1,0);

/*
//titre
glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,0.55*cH/cF);
const unsigned char myString [] = "Exemple de réseau d'hexaèdres aléatoires";
const unsigned char* title = &myString[0];
glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, title ); 

//infos
glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,-0.55*cH/cF);
string sstrnF = "Nombre d'hexaèdres : " + to_string(nHEX) + "";
const char * cstrnF = sstrnF.c_str();
const unsigned char *info1 = reinterpret_cast< const unsigned char * >(cstrnF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info1 ); 

glColor3f(255.,255.,255.); 
glRasterPos2f(-0.5*cH/cF,-0.6*cH/cF);
string sstrpF = "Pourcentage d'hexaèdres : " + to_string(pHEX) + "";
const char * cstrpF = sstrpF.c_str();
const unsigned char *info2 = reinterpret_cast< const unsigned char * >(cstrpF); 
glutBitmapString(GLUT_BITMAP_8_BY_13, info2 ); 
*/


int i,j,k,l;
i=0;
j=0;
k=0;



glColor3f(1.,1.,1.);
// Box
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();
	// BACK
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();      
	// LEFT
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();
	// RIGHT
	glBegin(GL_LINE_LOOP);
	glVertex3i((i+cH-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd();   	
	// BOTTOM
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k-cH/2)/cF);    
      	glEnd();
	// UP
	glBegin(GL_LINE_LOOP);
	glVertex3i((i-cH/2)/cF, (j-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i+cH-cH/2)/cF,  (j-cH/2)/cF,(k+cH-cH/2)/cF);      
	glVertex3i((i+cH-cH/2)/cF,  (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);     
	glVertex3i((i-cH/2)/cF, (j+cH-cH/2)/cF,(k+cH-cH/2)/cF);    
      	glEnd(); 



setup_illumination();
 
for(l=0;l<nS;l++){ 
  glPushMatrix();
  glTranslated((LIST_X[l]-cH/2)/cF,(LIST_Y[l]-cH/2)/cF,(LIST_Z[l]-cH/2)/cF);
  glutSolidSphere(LIST_R[l],50,50);
  glPopMatrix();
} 


glutSwapBuffers();
	//Attention : pas SwapBuffers(DC) !
glutPostRedisplay();
}

static void refenetrer(int w, int h)
{
glViewport(0, 0, (GLsizei) w, (GLsizei) h);
/*
modification des tailles
*/
/*du tampon d?affichage
*/
glMatrixMode(GL_PROJECTION);
/* pile courante = projection
*/
glLoadIdentity();
gluPerspective(
45,
float(w)/float(h),
0.1,
20*cH
); 	//Pour les explications, lire le tutorial sur OGL et win
  	glMatrixMode(GL_MODELVIEW); 	//Optionnel


}

static void clavier(unsigned char touche, int x, int y)
{

	switch (touche)
	{	
       
	case 27: {
	// sortie si touche ESC
	exit(0);break; }
	case 32: {
	// remise à 0 si space
	theta_eye=0;
	phi_eye=0;	
	rayon=2.*((GLdouble) cH)/((GLdouble) cF);
	x_eye=0;
	y_eye=0;
	z_eye=rayon;
	x_g=280;
	y_g=280;
	x_d=280;
	y_d=280;
	break;
	}
	case 81 : {
	// translation //z <0 si touche Q
        z_eye-=0.1*cH/cF;
	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}
	case 113 : {
	// translation //z >0 si touche q
        z_eye+=0.1*cH/cF; 	
	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	break;	  
	}
	
	}
}

static void processSpecialKeys(int key, int x, int y) {
  
        switch(key){
        case GLUT_KEY_LEFT: {
	// translation //X <0 si touche <-
        x_eye-=0.1*cH/cF;	

	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}
        case GLUT_KEY_RIGHT : {
	// translation //X >0 si touche ->
        x_eye+=0.1*cH/cF;	

	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}	
        case GLUT_KEY_DOWN : {
	// translation //Y <0 si touche v
        y_eye-=0.1*cH/cF;	

	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}	
        case GLUT_KEY_UP : {
	// translation //Y >0 si touche ^
        y_eye+=0.1*cH/cF;	

	phi_eye  = atan2(x_eye,z_eye);
	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye= asin((y_eye)/(rayon));
	
	break;	  
	}	
	}
}

static void gerer_souris(int bouton, int etat, int x, int y)
{
switch(bouton){
case GLUT_LEFT_BUTTON:
if(etat==GLUT_DOWN){
click_g = 1;
x_g     = x;
//cout<<"x_g : "<< x_g <<endl;
y_g     = 560-y;
//cout<<"y_g : "<< y_g <<endl;  
}
break;
case GLUT_MIDDLE_BUTTON:
break;
default:
if(etat==GLUT_DOWN){
click_d = 1;
x_d     = x;
//cout<<"x_d : "<< x_d <<endl;
y_d     = 560-y;
//cout<<"y_d : "<< y_d <<endl;    
theta_eye+= (3.14159/4.*(R(y_d)-280.)/280.);
phi_eye+= (3.14159/4.*(R(x_d)-280.)/280.);  

x_eye = rayon*sin(3.14159/2-theta_eye)*sin(phi_eye);
y_eye = rayon*cos(3.14159/2-theta_eye);
z_eye = rayon*sin(3.14159/2-theta_eye)*cos(phi_eye);
  
}  
break;
}
}

static void gerer_souris_mouvement(int x, int y)
{
/*
position courante (x,y) de la souris
*/
}

inline void ExpMeshMULTICOR(R SIZE, int NB_SPH, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R){

R H_MIN=0.;  
  
fstream g;
g.open("carte",fstream::out);
g.precision(4);

g<<"CONDLIMI"<<endl;
g<<setw(12)<<scientific<<left<<H_MIN<<setw(12)<<scientific<<left<<SIZE<<setw(12)<<scientific<<left<<H_MIN<<setw(12)<<scientific<<left<<SIZE<<setw(12)<<scientific<<left<<H_MIN<<setw(12)<<scientific<<left<<SIZE<<endl;
g<<"SNUMERI"<<endl;
g<<setw(12)<<scientific<<0.000001<<setw(12)<<scientific<<0.000001<<endl;
g<<"SCONTAC"<<endl;
g<<setw(12)<<scientific<<0.000005<<setw(14)<<left<<999999999<<endl;
g<<"MATER"<<endl;
g<<setw(12)<<scientific<<9.81<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<210e9<<setw(12)<<scientific<<0.02<<setw(12)<<scientific<<7.8e9<<endl;

g<<setw(10)<<"SPHERE"<<NB_SPH<<endl;
for(int i=0;i<NB_SPH;i++)
{
 g<<setw(10)<<(i+1)<<setw(12)<<scientific<<LIST_R[i]<<endl;
 g<<setw(14)<<scientific<<LIST_X[i]<<setw(14)<<scientific<<LIST_Y[i]<<setw(14)<<scientific<<LIST_Z[i]<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
 g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
 g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
 g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
}

g<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(10)<<1<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<endl;

g<<setw(10)<<"HSINUS"<<1<<endl;
g<<setw(10)<<1<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<endl;

g<<setw(10)<<"HOBLIQUE"<<1<<endl;
g<<setw(10)<<1<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<endl;
g<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<endl;

g<<setw(10)<<"PAROI"<<6<<endl;
g<<setw(10)<<1<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<2<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<3<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<4<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<5<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<6<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<SIZE<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g.close();
} 

int main(int argc, char *argv[])
{
R t1=CPUtime();

// TAILLE DU DOMAINE
string str_H_TOT = argv[1];
long H_TOT;
from_string<long>(H_TOT,str_H_TOT,dec);
long V_TOT = H_TOT*H_TOT;
long N_TOT = H_TOT*H_TOT*H_TOT;
cH=H_TOT;
cF  = 1;

// TAILLE REELLE DU DOMAINE
string str_SIZE = argv[2];
R SIZE;
from_string<R>(SIZE,str_SIZE,dec);

// RAPPORT D'ECHELLE
string str_R_E = argv[3];
long R_E;
from_string<long>(R_E,str_R_E,dec);
rS=R_E;

// EVALUATION DU NOMBRE DE SPHERES
int NB_SPH=(R_E+1)*(R_E+1)*(R_E+1)+3*(R_E+1)*R_E*R_E;
R R_SPH=R(H_TOT)/(2.+4./sqrt(2)*(R_E));
R L_E=(R(H_TOT)-2.*R_SPH)/R(R_E);
R_SPH=sqrt(2)*L_E/4.;

cout<<"L_E:"<<L_E<<", R_SPH:"<<R_SPH<<endl;
cout<<"H_TOT:"<<H_TOT<<", Heval:"<<(2*R_SPH+(R_E)*L_E)<<endl;
cout<<"Nombre de spheres : "<<NB_SPH<<endl;

// Paramètres des billes
LIST_X = new R[NB_SPH];
LIST_Y = new R[NB_SPH];
LIST_Z = new R[NB_SPH];
LIST_R = new R[NB_SPH];

NB_SPH=0;
for(int kte=0;kte<R_E+1;kte++){    
  for(int jte=0;jte<R_E+1;jte++){  
    for(int ite=0;ite<R_E+1;ite++){
    
    LIST_X[NB_SPH]=R_SPH+ite*L_E;LIST_Y[NB_SPH]=R_SPH+jte*L_E;LIST_Z[NB_SPH]=R_SPH+kte*L_E;LIST_R[NB_SPH]=R_SPH;NB_SPH++;
    
    }
  }
}

for(int kte=0;kte<R_E+1;kte++){    
  for(int jte=0;jte<R_E;jte++){  
    for(int ite=0;ite<R_E;ite++){
    
    LIST_X[NB_SPH]=R_SPH+L_E/2.+ite*L_E;LIST_Y[NB_SPH]=R_SPH+L_E/2.+jte*L_E;LIST_Z[NB_SPH]=R_SPH+kte*L_E;LIST_R[NB_SPH]=R_SPH;NB_SPH++;
    
    }
  }
}

for(int kte=0;kte<R_E;kte++){    
  for(int jte=0;jte<R_E;jte++){  
    for(int ite=0;ite<R_E+1;ite++){
    
    LIST_X[NB_SPH]=R_SPH+ite*L_E;LIST_Y[NB_SPH]=R_SPH+L_E/2.+jte*L_E;LIST_Z[NB_SPH]=R_SPH+L_E/2.+kte*L_E;LIST_R[NB_SPH]=R_SPH;NB_SPH++;
    
    }
  }
}

for(int kte=0;kte<R_E;kte++){    
  for(int jte=0;jte<R_E+1;jte++){  
    for(int ite=0;ite<R_E;ite++){
    
    LIST_X[NB_SPH]=R_SPH+L_E/2.+ite*L_E;LIST_Y[NB_SPH]=R_SPH+jte*L_E;LIST_Z[NB_SPH]=R_SPH+L_E/2.+kte*L_E;LIST_R[NB_SPH]=R_SPH;NB_SPH++;
    
    }
  }
}


cout<<"Nombre de spheres : "<<NB_SPH<<endl;
nS=NB_SPH;

// ExpMesh
for(int it=0;it<NB_SPH;it++){
LIST_X[it]=LIST_X[it]*SIZE/R(H_TOT);
LIST_Y[it]=LIST_Y[it]*SIZE/R(H_TOT);	
LIST_Z[it]=LIST_Z[it]*SIZE/R(H_TOT);	
LIST_R[it]=LIST_R[it]*SIZE/R(H_TOT);	
}

ExpMeshMULTICOR(SIZE,NB_SPH,LIST_X,LIST_Y,LIST_Z,LIST_R);


// AFFICHAGE

glutInit(&argc, argv);
glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE);
glutInitWindowSize(560,560);
glutInitWindowPosition(200,240);

InitGL();
glutCreateWindow("Empilement de billes");

glutReshapeFunc(refenetrer);
glutDisplayFunc(afficher);


glutKeyboardFunc(clavier);
glutSpecialFunc(processSpecialKeys);
glutMouseFunc(gerer_souris);
glutMotionFunc(gerer_souris_mouvement);

glutMainLoop();

char quit;
cin>> quit;
return 1;
}
