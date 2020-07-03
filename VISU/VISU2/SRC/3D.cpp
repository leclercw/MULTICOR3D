/********************************************************/
/*                     3D.c                             */
/********************************************************/
/* Visualisation en 3D de système multicontacts         */
/* Developpé par Fortin                                 */
/* Modifié le 6/03/05                                   */
/********************************************************/

/* inclusion des fichiers d'en-tete Glut */
#define GL_GLEXT_PROTOTYPES 1
#define GL3_PROTOTYPES 1

//#include <GL3/GL3.h>
#include <GL/freeglut.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <sys/time.h> 
#include <sys/resource.h> 
#include <string.h>
#include <unistd.h>
#include <time.h> 
#include <sys/time.h> 

#include "gl2ps.h"

using namespace std;

#define MAXSPHERE 1000000
#define MAXCONT 1000000
#define MAXSUBD 175000                       
#define MAXVOXEL 27000000
#define MAXELF 100000
#define MAXPAROI 400
#define MAX_PATH_LENGTH 1024
#define LDEBUG 0

// Prototype des fonctions 

void affichage();
void clavier(unsigned char key, int x, int y);
void processSpecialKeys(int key, int x, int y); 
void idle();
void menu(int choice);
void multi();
void reshape(int x,int y);
void mouse(int bouton,int etat,int x,int y);
void mousemotion(int x,int y);
void deffic();
void eval_sphere();
void setup_illumination();

void writefile(int format, int sort, int options, int nbcol,char *filename,const char *extension);
void sauvegarde(int pas);

// variables 
       
FILE *P_FICHIER; 
//ifstream g;
char NOM_FICHIER[MAX_PATH_LENGTH];
char NOM_FEPS[40];

int Nt=20;
int Np=20;
int NNt,NNp;

float views[MAXSUBD];

const  double  pi2 = 6.28318530718; 

// Data
GLuint vboId;

static float rotation = -60.;
static int nbsphere, nbparoi, nbver;
static int NX, NY, NZ, NTOT,NCONT;
static int jj;
static float qmaxv, qminv, qmaxh, qminh, qmaxz, qminz;
static float svmmin,svmmax,tracemax,tracemin,sig11max,sig11min,sig22max,sig22min,sig33max,sig33min;
static float sig12max,sig12min,sig13max,sig13min,sig23max,sig23min;
static float sig1max,sig1min,sig2max,sig2min,sig3max,sig3min;
static float dep1max,dep1min,dep2max,dep2min,dep3max,dep3min,maxt,mint,minfx,maxfx,minfy,maxfy,minfz,maxfz;
static float def11max,def11min,def22max,def22min,def33max,def33min;
static float defe11max,defe11min,defe22max,defe22min,defe33max,defe33min;
static float svmminv,svmmaxv,longmoy;
float xc[MAXSPHERE],  yc[MAXSPHERE], zc[MAXSPHERE], ray[MAXSPHERE],t;
int edge[MAXSPHERE];
int typs[MAXSPHERE];
float coul1[10], coul2[10], coul3[10], vcoulmax[10] ;
float SVM[MAXSPHERE],TRACE[MAXSPHERE],SIG11[MAXSPHERE],SIG22[MAXSPHERE],SIG33[MAXSPHERE];
float SIG12[MAXSPHERE],SIG13[MAXSPHERE],SIG23[MAXSPHERE];
float SIG1[MAXSPHERE],SIG2[MAXSPHERE],SIG3[MAXSPHERE];
float DEP1[MAXSPHERE],DEP2[MAXSPHERE],DEP3[MAXSPHERE],TEMP[MAXSPHERE],FLX[MAXSPHERE],FLY[MAXSPHERE],FLZ[MAXSPHERE];
float * vertices;

int PH[MAXSPHERE];
float DEF11[MAXSPHERE],DEF22[MAXSPHERE],DEF33[MAXSPHERE];
float DEFE11[MAXSPHERE],DEFE22[MAXSPHERE],DEFE33[MAXSPHERE];
float vali,tps,maxv,invmaxv,xfx1, xfx2, xfx3, xfx4, yfy1, yfy2, yfy3, yfy4;
float xmin1,xmin2,xmax1,xmax2,ymin1,ymin2,ymax1,ymax2;
float qcjx[MAXPAROI], qcjy[MAXPAROI], qcjz[MAXPAROI];
float qckx[MAXPAROI], qcky[MAXPAROI], qckz[MAXPAROI];
float qclx[MAXPAROI], qcly[MAXPAROI], qclz[MAXPAROI];
float qcmx[MAXPAROI], qcmy[MAXPAROI], qcmz[MAXPAROI];
float xcf[MAXELF][4], ycf[MAXELF][4], zcf[MAXELF][4], VONMISV[MAXELF];
float vcoh[MAXCONT][8];


static GLsizei window_w = 0; 
static GLsizei window_h = 0;
static int nfeps = 1 ;
static int incr = 0 ;
char nfeps_str[10];
static int lanime = 0 ;
static int ltype = 0 ;

bool bmove=1;

float X_TOT;
float Y_TOT;
float Z_TOT;

float X_TOT2;
float Y_TOT2;
float Z_TOT2;

static float pcut=0.;
static float angview=70.;
float posrel=0.;
float theta_eye;
float phi_eye;
float rayon;
float x_eye;
float y_eye;
float z_eye;
int x_g;
int x_d;
int y_g;
int y_d;
int click_g;
int click_d;

float dx=0.;
float dy=1.;
float dz=0.;

int type_ver;

char presse;
int anglex,angley,x,y,xold,yold,ntype,ind;


int main ( int argc, char **argv ) 
{ 

    const rlim_t kStackSize = 512 * 2048 * 2048;   // min stack size = 32 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
      //    cout << "StackLimit - expected : " << rl.rlim_cur << " - " << kStackSize << endl;

        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
	  cout << "Upgraded StackLimit soft - max : " << rl.rlim_cur << " - " << rl.rlim_max << endl;
        }
    }


  // Ouverture fichier 
  deffic();

  // Para sphères
  eval_sphere();

  // initialisation de glut et creation de la fenetre 
  glutInit(&argc, argv);  
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH);  
  glutInitWindowPosition(100,100);
  window_w = 900; 
  window_h = 900;
  glutInitWindowSize(window_w,window_h);
  glutCreateWindow("Visualisation multicor");

  // Initialisation d'OpenGL 
  glClearColor(1,1,1,1);
  glEnable(GL_DEPTH_TEST);
 
  // Parametres "joystick" 
	x_eye = -dx*2.*X_TOT;
	y_eye = -dy*2.*Y_TOT;
	z_eye = -dz*2.*Z_TOT;
	
	jj=1;

	x_g=window_w/2;
	y_g=window_h/2;
	x_d=window_w/2;
	y_d=window_h/2;

	rayon=sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye); 
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));

  glutDisplayFunc(affichage);
  glutIdleFunc(idle);
  glutReshapeFunc(reshape);

  glutKeyboardFunc(clavier);
  glutSpecialFunc(processSpecialKeys);
  glutMouseFunc(mouse); 
  glutMotionFunc(mousemotion);

   // enregistrement des fonctions de rappel 

  glutCreateMenu(menu);
  glutAddMenuEntry("[a] Stoppe animation", 1);
  glutAddMenuEntry("[r] Demarrer animation", 2);
  glutAddMenuEntry("[s] Pas a pas", 3);
  glutAddMenuEntry("[f] Format PS ou PDF", 4);
  glutAddMenuEntry("[i] Image", 5);
  glutAddMenuEntry("[v] Video", 6);
  glutAddMenuEntry("[esc] Quit", 7);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
     
    // Entree dans la boucle principale glut 
  glutMainLoop();
  return(0) ;
}  


void processSpecialKeys(int key, int x, int y) {
  
        switch(key){
        case GLUT_KEY_LEFT: {
	// translation //Y <0 si touche v
        x_eye-=0.1*dy*X_TOT;	
        y_eye-=0.1*dz*Y_TOT;	
        z_eye-=0.1*dx*Y_TOT;	

	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));	
	break;	  
	}
        case GLUT_KEY_RIGHT : {
	// translation //Y >0 si touche ^
        x_eye+=0.1*dy*X_TOT;	
        y_eye+=0.1*dz*Y_TOT;	
        z_eye+=0.1*dx*Y_TOT;	

	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));
	break;	  
	}	
        case GLUT_KEY_DOWN : {
        x_eye-=0.1*dz*X_TOT;	
        y_eye-=0.1*dx*Y_TOT;	
        z_eye-=0.1*dy*Y_TOT;	

	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));
	break;	  
	}	
        case GLUT_KEY_UP : {
        x_eye+=0.1*dz*X_TOT;	
        y_eye+=0.1*dx*Y_TOT;	
        z_eye+=0.1*dy*Y_TOT;		

	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));
	break;	  
	}
        case GLUT_KEY_PAGE_DOWN : {
        Nt=Nt-10;
        Np=Np-10;
        if(Nt<5) Nt=5;
        if(Np<5) Np=5;
  	eval_sphere();
        cout<<"Modif. Prec. :"<<Nt<<","<<Np<<endl;
	break;	  
        }
        case GLUT_KEY_PAGE_UP : {
        Nt=Nt+10;
        Np=Np+10;
 	eval_sphere();
        cout<<"Modif. Prec. :"<<Nt<<","<<Np<<endl;
	break;	  
        }

	
	}
}


void setup_illumination(float kx, float ky, float kz, float kt){
	
  // Intialise and set lighting parameters
  GLfloat light_pos[] = {300.*(dx+dy+dz),-300.*(dx+dy+dz),300.*(dx+dy+dz), 1.};
  GLfloat light_ka[] = {kx, ky, kz, kt};
  GLfloat light_kd[] = {1., 1., 1., 1.};
  GLfloat light_ks[] = {1., 1., 1., 1.};
  
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
  
  
void affichage(){

cout<<"Nt:"<<Nt<<endl;

  if(incr==0) {
glGenBuffers(1, &vboId);   
 }
 if(lanime==1||incr==0) {
 multi() ;
 }

float tcpu1=((float) clock())/CLOCKS_PER_SEC;

 float x,y,z;
  
glClear
(
GL_COLOR_BUFFER_BIT |
GL_DEPTH_BUFFER_BIT
); 	//Efface le frame buffer et le Z-buffer
glMatrixMode(GL_MODELVIEW); 	//Choisit la matrice MODELVIEW
glLoadIdentity(); 	//Réinitialise la matrice
gluLookAt(x_eye,y_eye,z_eye,0.,0.,0.,dz,dx,dy);
int i,j,k,l;
i=0;
j=0;
k=0;

if(bmove){

//Dir z
glColor3f(1.,0.,0.);   
GLUquadricObj *pObj;		//definition du pointeur d'objet
pObj = gluNewQuadric();		//creation d'un objet(Retourne 0 si No Memory)
gluQuadricNormals(pObj, GLU_SMOOTH);		//application des normals
glPushMatrix();	
glTranslated((-0.7*X_TOT),(-0.7*Y_TOT),(-0.7*Z_TOT));
gluCylinder(pObj, 1.,1., 0.15*Z_TOT, 10, 10);
glTranslated(0.,0.,0.15*Z_TOT);
glutSolidCone(0.02*Z_TOT, 0.05*Z_TOT, 10, 10);
glPopMatrix();	
gluDeleteQuadric (pObj);		//detruit l'objet


//Dir y
glColor3f(0.,1.,0.);   
GLUquadricObj *pObj2;		//definition du pointeur d'objet
pObj2 = gluNewQuadric();		//creation d'un objet(Retourne 0 si No Memory)
gluQuadricNormals(pObj2, GLU_SMOOTH);		//application des normals
glPushMatrix();	
glTranslated((-0.7*X_TOT),(-0.7*Y_TOT),(-0.7*Z_TOT));
glRotatef(-90,1,0,0);
gluCylinder(pObj2, 1.,1., 0.15*Z_TOT, 10, 10);
glTranslated(0.,0.,0.15*Z_TOT);
glutSolidCone(0.02*Z_TOT, 0.05*Z_TOT, 10, 10);
glPopMatrix();	
gluDeleteQuadric (pObj2);		//detruit l'objet


//Dir x
glColor3f(0.,0.,1.);   
GLUquadricObj *pObj3;		//definition du pointeur d'objet
pObj3 = gluNewQuadric();		//creation d'un objet(Retourne 0 si No Memory)
gluQuadricNormals(pObj3, GLU_SMOOTH);		//application des normals
glPushMatrix();	
glTranslated((-0.7*X_TOT),(-0.7*Y_TOT),(-0.7*Z_TOT));
glRotatef(90,0,1,0);
gluCylinder(pObj3, 1.,1., 0.15*Z_TOT, 10, 10);
glTranslated(0.,0.,0.15*Z_TOT);
glutSolidCone(0.02*Z_TOT, 0.05*Z_TOT, 10, 10);
glPopMatrix();	
gluDeleteQuadric (pObj3);	

}
else{

glColor3f(0.0f,0.0f,0.0f);

glMatrixMode(GL_PROJECTION);
glPushMatrix();
glLoadIdentity();
gluOrtho2D(0.0, window_w, 0.0, window_h);
glMatrixMode(GL_MODELVIEW);
glPushMatrix();
glLoadIdentity();

if(dx==1){
char txt[32] = "Plan yx" ;
glRasterPos2f(0.1*window_w,0.1*window_h);
gl2psText(txt, "Helvetica", 24);

int len = (int) strlen(txt);
for (int i = 0; i < len; i++){ glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, txt[i] ); }
}
else if(dy==1){
char txt[32] = "Plan xz" ;
glRasterPos2f(0.1*window_w,0.1*window_h);
gl2psText(txt, "Helvetica", 24);

int len = (int) strlen(txt);
for (int i = 0; i < len; i++){ glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, txt[i] ); }
}
else if(dz==1){
char txt[32] = "Plan yx" ;
glRasterPos2f(0.1*window_w,0.1*window_h);
gl2psText(txt, "Helvetica", 24);

int len = (int) strlen(txt);
for (int i = 0; i < len; i++){ glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, txt[i] ); }
}

char txt0[32] = "Pos. rel. du plan : " ;
char txt1[32];
sprintf(txt1,"%f",posrel);
strcat(txt0,txt1);

glRasterPos2f(0.1*window_w,0.06*window_h);
gl2psText(txt0, "Helvetica", 24);

int len0 = (int) strlen(txt0);
for (int i = 0; i < len0; i++){ glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, txt0[i] ); }

glMatrixMode(GL_MODELVIEW);
glPopMatrix();
glMatrixMode(GL_PROJECTION);
glPopMatrix();

}

  glColor3f(0.,0.,0.);  
    
  // Draw paroi
  //for (i=1;i<=nbparoi;i++ ) 
  //{ 
  // glBegin(GL_LINE_LOOP);
  //  x = qcjx[i]-X_TOT/2. ;
  //  y = qcjy[i]-Y_TOT/2.  ;
  //  z = qcjz[i]-Z_TOT/2. ;
  //  glVertex3f(x,y,z);
  //  x = qckx[i]-X_TOT/2.  ;
  //  y = qcky[i]-Y_TOT/2. ;
  //  z = qckz[i]-Z_TOT/2. ;
  //  glVertex3f(x,y,z);
  //  x = qclx[i]-X_TOT/2.  ;
  //  y = qcly[i]-Y_TOT/2. ;
  //  z = qclz[i]-Z_TOT/2. ;
  //  glVertex3f(x,y,z);
  //  x = qcmx[i]-X_TOT/2.  ;
  //  y = qcmy[i]-Y_TOT/2. ;
  //  z = qcmz[i]-Z_TOT/2. ;
  //  glVertex3f(x,y,z);
  //  glEnd();
 // }


if(ntype==3){ 

         setup_illumination(1.,0.,0.,1.);
  
	// Draw sphere 

	float xx1,yy1,zz1;
	float xx2,yy2,zz2;
	float angx,angy,angz;
	float longl,longp;

		for (i=1;i<=NCONT;i++ )
		{

				xx1=vcoh[i][1]-X_TOT/2.;
				yy1=vcoh[i][2]-Y_TOT/2.;
				zz1=vcoh[i][3]-Z_TOT/2.;
				xx2=vcoh[i][4]-X_TOT/2.;
				yy2=vcoh[i][5]-Y_TOT/2.;
				zz2=vcoh[i][6]-Z_TOT/2.;

			glPushMatrix();	
			glTranslated(xx1,yy1,zz1);
			glutSolidSphere(longmoy/8.,10, 10);
			glPopMatrix();	

			glPushMatrix();	
			glTranslated(xx2,yy2,zz2);
			glutSolidSphere(longmoy/8.,10, 10);
			glPopMatrix();	

			longl=sqrt((xx2-xx1)*(xx2-xx1)+(yy2-yy1)*(yy2-yy1)+(zz2-zz1)*(zz2-zz1));
			longp=sqrt((xx2-xx1)*(xx2-xx1)+(yy2-yy1)*(yy2-yy1));
			float bx,by,bz;
			
			if(longp>1e-10){
			bx=(yy1-yy2)/longp;
			by=(xx2-xx1)/longp;
			bz=0.;

			angx=acos((zz2-zz1)/longl)*180./3.14159265359;

			}else{
			bx=0.;
			by=1.;
			bz=0.;
                        angx=((zz2-zz1)>0.)?0.:-180.;
			}

			GLUquadricObj *pObjC;		//definition du pointeur d'objet
			pObjC = gluNewQuadric();		//creation d'un objet(Retourne 0 si No Memory)
			gluQuadricNormals(pObjC, GLU_SMOOTH);		//application des normals
			glPushMatrix();	

			glTranslated(xx1,yy1,zz1);
			glRotatef(angx,bx,by,bz);          		

			gluCylinder(pObjC, longmoy/10.,longmoy/10., longl, 10, 10);
			glPopMatrix();	
			gluDeleteQuadric (pObjC);


		}



	//Disable lighting
	glDisable (GL_LIGHTING);
	glDisable (GL_LIGHT0);	

}
else if(ntype==2){ 


	// Draw sphere 

	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(GLfloat)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

        // Draw wall

		GLfloat Transl[] = {0.5,0.5,0.5, 0.1};
		glColor4fv(Transl); 

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		for (i=1;i<=nbparoi;i++ )
		{
		glBegin(GL_QUADS);
			glVertex3f(qcjx[i]-X_TOT/2.,qcjy[i]-Y_TOT/2.,qcjz[i]-Z_TOT/2.); 
			glVertex3f(qckx[i]-X_TOT/2.,qcky[i]-Y_TOT/2.,qckz[i]-Z_TOT/2.); 
			glVertex3f(qclx[i]-X_TOT/2.,qcly[i]-Y_TOT/2.,qclz[i]-Z_TOT/2.); 
			glVertex3f(qcmx[i]-X_TOT/2.,qcmy[i]-Y_TOT/2.,qcmz[i]-Z_TOT/2.); 
		glEnd();
                }

		glDisable(GL_BLEND);

}
else if(ntype==1){


 if(ltype==0){
//cout<<"hello2"<<endl;

        setup_illumination(0.,0.,0.2,1.);
  
	// Draw sphere 

	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(float)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);


	//Disable lighting
	glDisable (GL_LIGHTING);
	glDisable (GL_LIGHT0);	

   }else
{
	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(GLfloat)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

		
        //Colorbar+title
		
	double coordx, coordy, deltax, deltay ;
	deltax = window_w/40. ; deltay = window_h/400. ;
	coordx = 4.27*window_w/5.; coordy = 0.65*window_h ;


	float fVal,fVal2 ;
        const char * unite;

	if(ltype==11)
	{			
	fVal = svmmin;
	fVal2 = svmmax;
	unite = "Contr. VMises (Pa)";	
	}
	else if(ltype==1)
	{			
	fVal = sig11min;
	fVal2 = sig11max;
	unite = "Contrainte Sig11 (Pa)";	
	}		
	else if(ltype==2)
	{			
	fVal = sig22min;
	fVal2 = sig22max;
	unite = "Contrainte Sig22 (Pa)";	
	}				        
	else if(ltype==3)
	{			
	fVal = sig33min;
	fVal2 = sig33max;
	unite = "Contrainte Sig33 (Pa)";
	}	
	else if(ltype==4)
	{			
	fVal = sig1min;
	fVal2 = sig1max;
	unite = "Contrainte Sig1 (Pa)";
	}		
	else if(ltype==5)
	{			
	fVal = sig2min;
	fVal2 = sig2max;
	unite = "Contrainte Sig2 (Pa)";	
	}				        
	else if(ltype==6)
	{			
	fVal = sig3min;
	fVal2 = sig3max;
	unite = "Contrainte Sig3 (Pa)";
	}		
	else if(ltype==7)
	{			
	fVal = sig12min;
	fVal2 = sig12max;
	unite = "Contrainte Sig12 (Pa)";
	}		
	else if(ltype==8)
	{			
	fVal = sig13min;
	fVal2 = sig13max;
	unite = "Contrainte Sig13 (Pa)";
	}				        
	else if(ltype==9)
	{			
	fVal = sig23min;
	fVal2 = sig23max;
	unite = "Contrainte Sig23 (Pa)";	
	}		
	else if(ltype==10)
	{			
	fVal = tracemin;
	fVal2 = tracemax;
	unite = "Contrainte Hydro. (Pa)";
	}			 
    
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, window_w, 0.0, window_w);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

                        //Colorbar
			for (k=1;k<=100;k++ ) 
			{
			glColor3f(1.,1.-k/100.,0.) ;
			x=coordx ;
			y=coordy + (k-1)*deltay ;

			glPolygonMode(GL_FRONT, GL_FILL);
			glBegin(GL_POLYGON);
			glVertex2f (x, y) ;
			glVertex2f (x+deltax, y) ;
			glVertex2f (x+deltax, y+deltay) ;
			glVertex2f (x, y+deltay) ;
			glVertex2f (x, y) ;
			glEnd();
			}

			glLineWidth(2.0);
			glColor3f(0.0f,0.0f,0.0f);
			glBegin(GL_LINE_LOOP); 
			glVertex2f(coordx,coordy) ;
			glVertex2f(coordx+deltax,coordy) ;
			glVertex2f(coordx+deltax,coordy+100.*deltay) ;
			glVertex2f(coordx,coordy+100.*deltay) ;
			glVertex2f(coordx,coordy) ;
			glEnd();
		   
			//Colorbar values
			char cVal[32] ;
			sprintf(cVal,"%.2E",fVal) ;
			cVal[9]='\0'  ;
			glRasterPos2f(coordx+deltax+deltax/6., coordy);
			gl2psText(cVal, "Helvetica", 24);
			for (k=0;k<=8;k++ ) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,cVal[k]);
                        }
			
			sprintf(cVal,"%.2E",fVal2) ;
			cVal[9]='\0'  ;
			glRasterPos2f(coordx+deltax+deltax/6., coordy+95.*deltay);
			gl2psText(cVal, "Helvetica", 24);
			for (k=0;k<=8;k++ ) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,cVal[k]);
                        }
		

		        //Title
			glRasterPos2f(4.7*coordx/5., 0.928*window_h);	

			glColor3f(0.f,0.f,0.f) ;
			gl2psText(unite, "Helvetica", 36);

			for (int i = 0; i < (int) strlen(unite); i++){
			glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, unite[i] );
			}			

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

}

  }
else if(ntype==4){


 if(ltype==0){


        setup_illumination(0.,0.,0.2,1.);
  
	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(GLfloat)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	//Disable lighting
	glDisable (GL_LIGHTING);
	glDisable (GL_LIGHT0);	

   }
   else{

	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(GLfloat)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
		
        //Colorbar+title
		
	double coordx, coordy, deltax, deltay ;
	deltax = window_w/40. ; deltay = window_h/400. ;
	coordx = 4.27*window_w/5.; coordy = 0.65*window_h ;


	float fVal,fVal2 ;
        const char * unite;

	if(ltype==1)
	{			
	fVal = def11min;
	fVal2 = def11max;
	unite = "Deformation Eps11";	
	}		
	else if(ltype==2)
	{			
	fVal = def22min;
	fVal2 = def22max;
	unite = "Deformation Eps22";	
	}				        
	else if(ltype==3)
	{			
	fVal = def33min;
	fVal2 = def33max;
	unite = "Deformation Eps33";
	}
	else if(ltype==4)
	{			
	fVal = defe11min;
	fVal2 = defe11max;
	unite = "Deformation Epse11";	
	}		
	else if(ltype==5)
	{			
	fVal = defe22min;
	fVal2 = defe22max;
	unite = "Deformation Epse22";	
	}				        
	else if(ltype>=6)
	{			
	fVal = defe33min;
	fVal2 = defe33max;
	unite = "Deformation Epse33";
	}

    
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, window_w, 0.0, window_w);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

                        //Colorbar
			for (k=1;k<=100;k++ ) 
			{
			glColor3f(k/100.,0.,1.-k/100.) ;
			x=coordx ;
			y=coordy + (k-1)*deltay ;

			glPolygonMode(GL_FRONT, GL_FILL);
			glBegin(GL_POLYGON);
			glVertex2f (x, y) ;
			glVertex2f (x+deltax, y) ;
			glVertex2f (x+deltax, y+deltay) ;
			glVertex2f (x, y+deltay) ;
			glVertex2f (x, y) ;
			glEnd();
			}

			glLineWidth(2.0);
			glColor3f(0.0f,0.0f,0.0f);
			glBegin(GL_LINE_LOOP); 
			glVertex2f(coordx,coordy) ;
			glVertex2f(coordx+deltax,coordy) ;
			glVertex2f(coordx+deltax,coordy+100.*deltay) ;
			glVertex2f(coordx,coordy+100.*deltay) ;
			glVertex2f(coordx,coordy) ;
			glEnd();
		   
			//Colorbar values
			char cVal[32] ;
			sprintf(cVal,"%.2E",fVal) ;
			cVal[9]='\0'  ;
			glRasterPos2f(coordx+deltax+deltax/6., coordy);
			gl2psText(cVal, "Helvetica", 24);
			for (k=0;k<=8;k++ ) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,cVal[k]);
                        }

			sprintf(cVal,"%.2E",fVal2) ;
			cVal[9]='\0'  ;
			glRasterPos2f(coordx+deltax+deltax/6., coordy+95.*deltay);
			gl2psText(cVal, "Helvetica", 24);
			for (k=0;k<=8;k++ ) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,cVal[k]);
                        }
		
		        //Title
			glRasterPos2f(4.7*coordx/5., 0.928*window_h);	

			glColor3f(0.f,0.f,0.f) ;
			gl2psText(unite, "Helvetica", 36);

			for (int i = 0; i < (int) strlen(unite); i++){
			glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, unite[i] );
			}			

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

}

  }else if(ntype==5){ 
  
  if(ltype==0){

       setup_illumination(0.,0.,0.2,1.);
  
	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(GLfloat)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	//Disable lighting
	glDisable (GL_LIGHTING);
	glDisable (GL_LIGHT0);	

   }
   else{
     
	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(GLfloat)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

			double coordx, coordy, deltax, deltay ;
			deltax = window_w/40. ; deltay = window_h/400. ;
			coordx = 4.27*window_w/5.; coordy = 0.65*window_h ;

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();
			gluOrtho2D(0.0, window_w, 0.0, window_w);
			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			glLoadIdentity();


float ka,kd,ks;
			for (k=1;k<=100;k++ ) {
////////////////////////////
			vali=4.*(1.-k/100.);
			if(vali<=1){ka=0.;kd=vali;ks=1.;}
			else if(vali<=2){ka=0.;kd=1.;ks=2.-vali;}
			else if(vali<=3){ka=vali-2.;kd=1.;ks=0.;}
			else if(vali<=4){ka=1.;kd=4.-vali;ks=0.;}
///////////////
			glColor3f(ka,kd,ks) ;
				x=coordx ;
		        
				y=coordy + (k-1)*deltay ;
				xfx2=x+deltax ;
				yfy2=y ;
				xfx3=xfx2 ;
				yfy3=y+deltay ;
				xfx4=x ;
				yfy4=yfy3 ;
				glPolygonMode(GL_FRONT, GL_FILL);
				glBegin(GL_POLYGON);
				glVertex2f (x, y) ;
				glVertex2f (xfx2, yfy2) ;
				glVertex2f (xfx3, yfy3) ;
				glVertex2f (xfx4, yfy4) ;
				glVertex2f (x, y) ;
				glEnd();
		        }
			
			xmin1=coordx ;
			ymin1=coordy ;

			xmin2=coordx+deltax ;
			ymin2=coordy ;

			xmax1=coordx+deltax ;
			ymax1=ymin1+100.*deltay ;

			xmax2=coordx ;
			ymax2=ymax1 ;

			glLineWidth(2.0);
			glColor3f(0.0f,0.0f,0.0f);
			glBegin(GL_LINES); 
			glVertex2f(xmin1,ymin1) ;
			glVertex2f(xmin2,ymin2) ;
			glVertex2f(xmin2,ymin2) ;
			glVertex2f(xmax1,ymax1) ;
			glVertex2f(xmax1,ymax1) ;
			glVertex2f(xmax2,ymax2) ;
			glVertex2f(xmax2,ymax2) ;
			glVertex2f(xmin1,ymin1) ;
			glEnd();

			float fVal,fVal2 ;
			const char * unite;

			if(ltype==1)
			{			
			fVal = mint;
			fVal2 = maxt;
			unite = "Temperature (degC)";	
			}
			else if(ltype==2)
			{			
			fVal = minfx;
			fVal2 = maxfx;
			unite = "FLUX X (W/m2)";	
			}		
			else if(ltype==3)
			{			
			fVal = minfy;
			fVal2 = maxfy;
			unite = "FLUX Y (W/m2)";		
			}
			else if(ltype>=4)
			{			
			fVal = minfz;
			fVal2 = maxfz;
			unite = "FLUX Z (W/m2)";		
			}
  
    	  
			//Colorbar values
			char cVal[32] ;
			sprintf(cVal,"%.2E",fVal) ;
			cVal[9]='\0'  ;
			glRasterPos2f(coordx+deltax+deltax/6., coordy);
			gl2psText(cVal, "Helvetica", 24);
			for (k=0;k<=8;k++ ) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,cVal[k]);
                        }
			
			sprintf(cVal,"%.2E",fVal2) ;
			cVal[9]='\0'  ;
			glRasterPos2f(coordx+deltax+deltax/6., coordy+95.*deltay);
			gl2psText(cVal, "Helvetica", 24);
			for (k=0;k<=8;k++ ) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,cVal[k]);
                        }
		

		        //Title
			glRasterPos2f(4.7*coordx/5., 0.928*window_h);	

			glColor3f(0.f,0.f,0.f) ;
			gl2psText(unite, "Helvetica", 36);

			for (int i = 0; i < (int) strlen(unite); i++){
			glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, unite[i] );
			}	



  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
	
    }

  }



  else if(ntype==0){ 

  if(ltype==0){

	setup_illumination(0.,0.,0.2,1.);

	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(GLfloat)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);


	//Disable lighting
	glDisable (GL_LIGHTING);
	glDisable (GL_LIGHT0);	


   }
   else{
     
	glUnmapBuffer(GL_ARRAY_BUFFER);

	//Render
	// Step 1
	glBindBuffer(GL_ARRAY_BUFFER,vboId);

	// Step 2
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	size_t nOffset =((2*NNp+2)*NNt)*nbsphere*sizeof(GLfloat)*3;
	size_t cOffset = 2*nOffset;

	// Step 3
	glVertexPointer(3, GL_FLOAT, 0, (void*)0);
	glNormalPointer(GL_FLOAT, 0, (void*)nOffset);
	glColorPointer(3, GL_FLOAT, 0, (void*)cOffset);

	// Step 4

	glDrawArrays(GL_TRIANGLE_STRIP,0,((2*NNp+2)*NNt)*nbsphere);

	// Step 5
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);


	float fVal,fVal2 ;
        const char * unite;

	if(ltype==1)
	{			
	fVal = dep1min;
	fVal2 = dep1max;
	unite = "Deplac. dep1 (m)";	
	}
	else if(ltype==2)
	{			
	fVal = dep2min;
	fVal2 = dep2max;
	unite = "Deplac. dep2 (m)";	
	}		
	else if(ltype>=3)
	{			
	fVal = dep3min;
	fVal2 = dep3max;
	unite = "Deplac. dep3 (m)";		
	}

	double coordx, coordy, deltax, deltay ;
	deltax = window_w/40. ; deltay = window_h/400. ;
	coordx = 4.27*window_w/5.; coordy = 0.65*window_h ;

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, window_w, 0.0, window_w);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
                        //Colorbar
			for (k=1;k<=100;k++ ) 
			{
			glColor3f(k/100.,0.,1.-k/100.) ;
			x=coordx ;
			y=coordy + (k-1)*deltay ;

			glPolygonMode(GL_FRONT, GL_FILL);
			glBegin(GL_POLYGON);
			glVertex2f (x, y) ;
			glVertex2f (x+deltax, y) ;
			glVertex2f (x+deltax, y+deltay) ;
			glVertex2f (x, y+deltay) ;
			glVertex2f (x, y) ;
			glEnd();
			}

			glLineWidth(2.0);
			glColor3f(0.0f,0.0f,0.0f);
			glBegin(GL_LINE_LOOP); 
			glVertex2f(coordx,coordy) ;
			glVertex2f(coordx+deltax,coordy) ;
			glVertex2f(coordx+deltax,coordy+100.*deltay) ;
			glVertex2f(coordx,coordy+100.*deltay) ;
			glVertex2f(coordx,coordy) ;
			glEnd();
		   
			//Colorbar values
			char cVal[32] ;
			sprintf(cVal,"%.2E",fVal) ;
			cVal[9]='\0'  ;
			glRasterPos2f(coordx+deltax+deltax/6., coordy);
			gl2psText(cVal, "Helvetica", 24);
			for (k=0;k<=8;k++ ) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,cVal[k]);
                        }
			
			sprintf(cVal,"%.2E",fVal2) ;
			cVal[9]='\0'  ;
			glRasterPos2f(coordx+deltax+deltax/6., coordy+95.*deltay);
			gl2psText(cVal, "Helvetica", 24);
			for (k=0;k<=8;k++ ) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,cVal[k]);
                        }
		

		        //Title
			glRasterPos2f(4.7*coordx/5., 0.928*window_h);	

			glColor3f(0.f,0.f,0.f) ;
			gl2psText(unite, "Helvetica", 36);

			for (int i = 0; i < (int) strlen(unite); i++){
			glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, unite[i] );
			}	
	  
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

    }

        // Draw wall

	GLfloat Transl[] = {0.5,0.5,0.5, 0.1};
	glColor4fv(Transl); 

        //setup_illumination(0.5,0.5,0.5, 0.1);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		for (i=1;i<=nbparoi;i++ )
		{
		glBegin(GL_QUADS);
			glVertex3f(qcjx[i]-X_TOT/2.,qcjy[i]-Y_TOT/2.,qcjz[i]-Z_TOT/2.); 
			glVertex3f(qckx[i]-X_TOT/2.,qcky[i]-Y_TOT/2.,qckz[i]-Z_TOT/2.); 
			glVertex3f(qclx[i]-X_TOT/2.,qcly[i]-Y_TOT/2.,qclz[i]-Z_TOT/2.); 
			glVertex3f(qcmx[i]-X_TOT/2.,qcmy[i]-Y_TOT/2.,qcmz[i]-Z_TOT/2.); 
		glEnd();
                }

		glDisable(GL_BLEND);

  }
 
 // Axes 
	if(bmove){
	glColor3f(0.,0.,1.);  
	glRasterPos3f(((0.22-0.7)*X_TOT),(-0.7*Y_TOT),(-0.7*Z_TOT));
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'x');
	glColor3f(0.,1.,0.);  
	glRasterPos3f(((-0.7)*X_TOT),((0.22-0.7)*Y_TOT),(-0.7*Z_TOT));
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'y');
	glColor3f(1.,0.,0.);  
	glRasterPos3f(((-0.7)*X_TOT),(-0.7*Y_TOT),((0.22-0.7)*Z_TOT));
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'z');
	}

 //On echange les buffers 
  glFlush();
  glutSwapBuffers();

float tcpu2=(((float) clock())/CLOCKS_PER_SEC)-tcpu1;
cout<<"tcpu_draw:"<<tcpu2<<endl;
}


void sauvegarde(int pas) {
  int opt,i ;
  const char *ext;
  static int format = GL2PS_EPS;
    lanime = 0 ;
//  glGenBuffers(1, &VertexVBOID);
//  glBindBuffer(GL_ARRAY_BUFFER,VertexVBOID); 
    multi() ;
    glutSetWindowTitle("Pas a pas");
    glutPostRedisplay();
    ext = (format == GL2PS_EPS) ? "eps" : "pdf";
    opt = GL2PS_NONE ;
    if (pas%1 == 0) {
    strcpy(NOM_FEPS,get_current_dir_name());
    if (nfeps<10) 
    { strcat(NOM_FEPS, "/IMAGE/image_visu00"); }
    else if (nfeps<100) 
    { strcat(NOM_FEPS, "/IMAGE/image_visu0"); }
    else
    { strcat(NOM_FEPS, "/IMAGE/image_visu"); }
    sprintf(nfeps_str,"%d",nfeps);
    strcat(NOM_FEPS, nfeps_str);
    writefile(format, GL2PS_SIMPLE_SORT, opt, 0, NOM_FEPS, ext);
    printf("GL2PS done... %s \n",NOM_FEPS);
    nfeps++;
    }
    pas=pas+1 ;
    glutTimerFunc(100,sauvegarde,pas++) ;
}

void clavier(unsigned char key, int x, int y){

  int opt;
  const char *ext;
  static int format = GL2PS_EPS;
  switch(key){
  case 'a':
    glutSetWindowTitle("Animation stoppe");
    lanime = 0 ;
    glutIdleFunc(NULL);
    break;
  case 'Q': {
if(bmove){
        x_eye-=0.1*dx*X_TOT;
        y_eye-=0.1*dy*Y_TOT;
        z_eye-=0.1*dz*Z_TOT;	

	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));
}
     break; }
  case 'q': { 
if(bmove){
        x_eye+=0.1*dx*X_TOT;
        y_eye+=0.1*dy*Y_TOT;
        z_eye+=0.1*dz*Z_TOT;	

	if(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye<2.*X_TOT*2.*X_TOT) {
        x_eye-=0.1*dx*X_TOT;
        y_eye-=0.1*dy*Y_TOT;
        z_eye-=0.1*dz*Z_TOT;	
        }


	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));
}
     break; }       
  case 'r':
    cout<<"On redemarre l'animation"<< endl;
    lanime = 1 ;
    glutSetWindowTitle("Animation redemarre");
    glutIdleFunc(idle);
    break;
  case 's':
    cout<<"Pas a pas ... "<< incr << ", " << tps <<endl;
    lanime = 0 ;
 // glGenBuffers(1, &VertexVBOID);
 // glBindBuffer(GL_ARRAY_BUFFER,VertexVBOID); 
    multi() ;
    glutSetWindowTitle("Pas a pas");
    glutPostRedisplay();
    break;
  case 'f':
    format = (format == GL2PS_EPS) ? GL2PS_PDF : GL2PS_EPS;
    printf("Print format changed to '%s'\n",
           (format == GL2PS_EPS) ? "EPS" : "PDF");
    break;  
  case '+':
    jj=jj+1;
    if(jj>NX) jj=NX;
   break;  
  case '-':
    jj=jj-1;
    if(jj<=0) jj=1;
   break;  
  case 'x':
   dx=1.;
   dy=0.;
   dz=0.;

// + remise à 0 
        pcut=0.;

	x_eye = -dx*2.*X_TOT;
	y_eye = -dy*2.*Y_TOT;
	z_eye = -dz*2.*Z_TOT; 

	x_g=window_w/2;
	y_g=window_h/2;
	x_d=window_w/2;
	y_d=window_h/2;

        angview=70.;
        glutReshapeWindow(window_w,window_h);
        cout<<"Plan zy"<<endl;

	rayon= sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));

   break;  
   case 'y':
   dx=0.;
   dy=1.;
   dz=0.;

// + remise à 0 
        pcut=0.;
	x_eye = -dx*2.*X_TOT;
	y_eye = -dy*2.*Y_TOT;
	z_eye = -dz*2.*Z_TOT; 

	x_g=window_w/2;
	y_g=window_h/2;
	x_d=window_w/2;
	y_d=window_h/2;

        angview=70.;
        glutReshapeWindow(window_w,window_h);
        cout<<"Plan xz"<<endl;

	rayon= sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));
   break;
   case 'z':
   dx=0.;
   dy=0.;
   dz=1.;

// + remise à 0 
        pcut=0.;
	x_eye = -dx*2.*X_TOT;
	y_eye = -dy*2.*Y_TOT;
	z_eye = -dz*2.*Z_TOT;

	x_g=window_w/2;
	y_g=window_h/2;
	x_d=window_w/2;
	y_d=window_h/2;

        angview=70.;
        glutReshapeWindow(window_w,window_h);
        cout<<"Plan yx"<<endl;

	rayon= sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));
   break;
  case 'i':

    ext = (format == GL2PS_EPS) ? "eps" : "pdf";
    opt = GL2PS_NONE ;
    strcpy(NOM_FEPS, get_current_dir_name());
    if (nfeps<10) 
    { strcat(NOM_FEPS, "/IMAGE/image_visu0"); }
    else
    { strcat(NOM_FEPS, "/IMAGE/image_visu"); }
    sprintf(nfeps_str,"%d",nfeps);
    strcat(NOM_FEPS, nfeps_str);
    writefile(format, GL2PS_SIMPLE_SORT, opt, 0, NOM_FEPS, ext);
    printf("GL2PS done... %s \n",NOM_FEPS);
    nfeps++;  
    break; 
  case 'v':
     glutTimerFunc(50,sauvegarde,0) ;
    break;     
  case 'p':
     bmove=0 ;
     cout<<"Mode plan actif"<<endl;
        
        pcut=0.;

	x_eye = -dx*1.5*X_TOT;
	y_eye = -dy*1.5*Y_TOT;
	z_eye = -dz*1.5*Z_TOT;

	x_g=window_w/2;
	y_g=window_h/2;
	x_d=window_w/2;
	y_d=window_h/2;

        angview=92.79;
        glutReshapeWindow(window_w,window_h);

	glutSpecialFunc(NULL);
	glutMouseFunc(NULL); 
	glutMotionFunc(NULL);

	rayon= sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));

    break;   
  case 'P':
     bmove=1 ;
     cout<<"Mode plan desactif"<<endl;
        pcut=0.;
	x_eye = -dx*2.*X_TOT;
	y_eye = -dy*2.*Y_TOT;
	z_eye = -dz*2.*Z_TOT;

	x_g=window_w/2;
	y_g=window_h/2;
	x_d=window_w/2;
	y_d=window_h/2;

        angview=70.;

        glutReshapeWindow(window_w,window_h);

	glutSpecialFunc(processSpecialKeys);
	glutMouseFunc(mouse); 
	glutMotionFunc(mousemotion);

	rayon= sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));


    break;     
  case 27: {
	// sortie si touche ESC
	exit(0);break; }
  case 32: {
	// remise à 0 si space
        bmove=1;
        pcut=0.;

	x_eye = -dx*2.*X_TOT;
	y_eye = -dy*2.*Y_TOT;
	z_eye = -dz*2.*Z_TOT;

	x_g=window_w/2;
	y_g=window_h/2;
	x_d=window_w/2;
	y_d=window_h/2;

        angview=70.;
        glutReshapeWindow(window_w,window_h);

        cout<<"Reinitialisation"<<endl;

	rayon= sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));
	
        break;
	}
	case 'C' : { 
        if(bmove){
	pcut=pcut-0.1;
	if(pcut<0.) {pcut=0.;}
        glutReshapeWindow(window_w,window_h);
        }
        else if(!bmove){
        x_eye-=0.1*dx*X_TOT;
        y_eye-=0.1*dy*Y_TOT;
        z_eye-=0.1*dz*Z_TOT;	

	if(x_eye<-dx*1.5*X_TOT) x_eye=-dx*1.5*X_TOT ;
	if(y_eye<-dy*1.5*Y_TOT) y_eye=-dy*1.5*Y_TOT ;
	if(z_eye<-dz*1.5*Z_TOT) z_eye=-dz*1.5*Z_TOT ;

	if(x_eye>-dx*0.5*X_TOT) x_eye=-dx*0.5*X_TOT;
	if(y_eye>-dy*0.5*Y_TOT) y_eye=-dy*0.5*Y_TOT;
	if(z_eye>-dz*0.5*Z_TOT) z_eye=-dz*0.5*Z_TOT;

        posrel=dx*(X_TOT/2.+X_TOT*(1+pcut)+x_eye)/X_TOT+dy*(Y_TOT/2.+Y_TOT*(1+pcut)+y_eye)/Y_TOT+dz*(Z_TOT/2.+Z_TOT*(1+pcut)+z_eye)/Z_TOT;
        angview=dx*180.*2.*atan2(1.05*X_TOT,X_TOT*max(1.,1.-posrel))/3.14159+dy*180.*2.*atan2(1.05*Y_TOT,Y_TOT*max(1.,1.-posrel))/3.14159+dz*180.*2.*atan2(1.05*Z_TOT,Z_TOT*max(1.,1.-posrel))/3.14159;
        glutReshapeWindow(window_w,window_h);

	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));

        cout<<"posi view:"<<dx*x_eye+dy*y_eye+dz*z_eye<<", posi plane:"<<posrel<<", angview:"<<angview<<endl;
        }
	break;	  
	}
	case 'c' : { 
        if(bmove){
        pcut=pcut+0.1;
        if(pcut>3.) {pcut=3.;}
        glutReshapeWindow(window_w,window_h);
        }
        else if(!bmove)
        {
        x_eye+=0.1*dx*X_TOT;
        y_eye+=0.1*dy*Y_TOT;
        z_eye+=0.1*dz*Z_TOT;	

	if(x_eye<-dx*1.5*X_TOT) x_eye=-dx*1.5*X_TOT ;
	if(y_eye<-dy*1.5*Y_TOT) y_eye=-dy*1.5*Y_TOT ;
	if(z_eye<-dz*1.5*Z_TOT) z_eye=-dz*1.5*Z_TOT ;

	if(x_eye>-dx*0.5*X_TOT) x_eye=-dx*0.5*X_TOT;
	if(y_eye>-dy*0.5*Y_TOT) y_eye=-dy*0.5*Y_TOT;
	if(z_eye>-dz*0.5*Z_TOT) z_eye=-dz*0.5*Z_TOT;

        posrel=dx*(X_TOT/2.+X_TOT*(1+pcut)+x_eye)/X_TOT+dy*(Y_TOT/2.+Y_TOT*(1+pcut)+y_eye)/Y_TOT+dz*(Z_TOT/2.+Z_TOT*(1+pcut)+z_eye)/Z_TOT;
        angview=dx*180.*2.*atan2(1.05*X_TOT,X_TOT*max(1.,1.-posrel))/3.14159+dy*180.*2.*atan2(1.05*Y_TOT,Y_TOT*max(1.,1.-posrel))/3.14159+dz*180.*2.*atan2(1.05*Z_TOT,Z_TOT*max(1.,1.-posrel))/3.14159;
        glutReshapeWindow(window_w,window_h);

	rayon    = sqrt(x_eye*x_eye+y_eye*y_eye+z_eye*z_eye);
	theta_eye  = dx*atan2(-z_eye,x_eye)+dy*atan2(-x_eye,y_eye)+dz*atan2(-y_eye,z_eye);
	phi_eye= dx*acos((y_eye)/(rayon))+dy*acos((z_eye)/(rayon))+dz*acos((x_eye)/(rayon));

        cout<<"posi view:"<<dx*x_eye+dy*y_eye+dz*z_eye<<", posi plane:"<<posrel<<", angview:"<<angview<<endl;
        }
	break;	  
	}    
	  case '0':
		ltype=0;
		break;     
	  case '1':
		ltype=1;
		break;       
	  case '2':
		ltype=2;
		break;      
	  case '3':
		ltype=3;
		break;         
	  case '4':
		ltype=4;
	        break;    	
	  case '5':
		ltype=5;
	        break;    
	  case '6':
		ltype=6;
	        break; 
	  case '7':
		ltype=7;
	        break;   
	  case '8':
		ltype=8;
	        break;   
	  case '9':
		ltype=9;
	        break;   
	  case '.':
		ltype=10;
	        break;   
	  case '*':
		ltype=11;
	        break;   
 	
  }

}

void eval_sphere(){

float thetap,stheta,sthetap,ctheta,cthetap;
float phip,sphip,cphip;

cthetap=1.;
sthetap=0.;
int nump=0;

	for (thetap = 2*3.141593/Nt; thetap < 2.00001*3.141593; thetap += 2*3.141593/Nt){
		ctheta=cthetap;
		stheta=sthetap;
		cthetap=cos(thetap);
		sthetap=sin(thetap);

	        views[3*nump]=0.;
	        views[3*nump+1]=0.;
	        views[3*nump+2]=1.;
		nump++;

	        views[3*nump]=0.;
	        views[3*nump+1]=0.;
	        views[3*nump+2]=1.; 
		nump++;

		for (phip = 3.14159/Np; phip < 1.00001*3.14159; phip += 3.14159/Np){
			       sphip=sin(phip);
			       cphip=cos(phip);

			        views[3*nump]=ctheta*sphip;
			        views[3*nump+1]=stheta*sphip;
			        views[3*nump+2]=cphip;
				nump++;

			        views[3*nump]=cthetap*sphip;
			        views[3*nump+1]=sthetap*sphip;
			        views[3*nump+2]=cphip; 
				nump++;

		}

	}
}


void deffic(){

int i ;
/* Première partie : Creer et remplir le fichier */

  cout<<"VISU deplac     num = 0  &  NOM fichier = vdepl "<<endl;
  cout<<"VISU contr      num = 1  &  NOM fichier = vcontr "<<endl;
  cout<<"VISU rupt       num = 2  &  NOM fichier = vrupt "<<endl;
  cout<<"VISU liens      num = 3  &  NOM fichier = vcoh "<<endl;
  cout<<"VISU def        num = 4  &  NOM fichier = vdef "<<endl;
  cout<<"VISU temp       num = 5  &  NOM fichier = vtemp "<<endl;
  cout<<"Entrez le num VISU :"<<endl;
  cin>>ntype;
  cout<<"Entrez le nom du fichier à lire : "<<endl;
  cin>>NOM_FICHIER;
  P_FICHIER = fopen(NOM_FICHIER, "r");  /* read */
//g.precision(4);
//g.open(NOM_FICHIER,fstream::in);	

  if (!P_FICHIER) 
//if(!g.is_open())
    {
     cout<<"Impossible d'ouvrir le fichier"<<endl;
     exit(-1);
    }  
  else
  {
    fscanf(P_FICHIER, "%f %f %f %f %f %f\n", &qmaxv,&qminv,&qmaxh,&qminh,&qmaxz,&qminz);
  //    g>>qmaxv>>qminv>>qmaxh>>qminh>>qmaxz>>qminz;
	Z_TOT2=fabs(qmaxz-qminz);
	Y_TOT2=fabs(qmaxv-qminv);
	X_TOT2=fabs(qmaxh-qminh);
	Z_TOT=100.;
	Y_TOT=100.;
	X_TOT=100.;
  
  }
/* Deuxieme partie :  appel de la fonction multi */
 // multi() ;
}


void idle(){
glutPostRedisplay();
}

void menu(int choice){

  int opt;
  const char *ext;
  static int format = GL2PS_EPS;
  
  switch (choice) {
  case 1:
    cout<<"On stoppe l'animation..."<<endl;
    lanime = 0 ;
    glutSetWindowTitle("Animation stoppe");
    glutIdleFunc(NULL);
    break;
  case 2:
    cout<<"On relance l'animation..."<<endl;
    lanime = 1 ;
    glutSetWindowTitle("Animation demarre");
    glutIdleFunc(idle);
    break;
  case 3:
    cout<<"On relance l'animation... "<<incr<<", "<<tps<<endl;
    glutSetWindowTitle("Pas a pas");
    lanime = 0 ;
 // glGenBuffers(1, &VertexVBOID);
 // glBindBuffer(GL_ARRAY_BUFFER,VertexVBOID); 
    multi() ;
    glutPostRedisplay();
    break;
  case 4:
   format = (format == GL2PS_EPS) ? GL2PS_PDF : GL2PS_EPS;
    printf("Print format changed to '%s'\n",
           (format == GL2PS_EPS) ? "EPS" : "PDF");
    break;    
  case 5:

    opt = GL2PS_OCCLUSION_CULL | GL2PS_DRAW_BACKGROUND;
    ext = (format == GL2PS_EPS) ? "eps" : "pdf";

    strcpy(NOM_FEPS,get_current_dir_name());
    if (nfeps<10) 
    { strcat(NOM_FEPS, "/IMAGE/image_visu00"); }
    else if (nfeps<100) 
    { strcat(NOM_FEPS, "/IMAGE/image_visu0"); }
    else
    { strcat(NOM_FEPS, "/IMAGE/image_visu"); }
    sprintf(nfeps_str,"%d",nfeps);
    strcat(NOM_FEPS, nfeps_str);

    writefile(format, GL2PS_SIMPLE_SORT, opt, 0, NOM_FEPS, ext);

    printf("GL2PS %d.%d.%d done with all images\n",
           GL2PS_MAJOR_VERSION, GL2PS_MINOR_VERSION, GL2PS_PATCH_VERSION);
    nfeps++;
    break;
  case 6:
    cout<<"On sort ..."<<endl;
    exit(1);
    break;
  }

}

void multi(){
float tcpu1=((float) clock())/CLOCKS_PER_SEC;

int i,ok ;

  incr++ ;
  if (!feof(P_FICHIER)) 
 // if (!g.eof()) 
  {
 
     if (ntype==0){

		fscanf(P_FICHIER, "%i %i %f \n", &nbsphere, &nbparoi, &tps);

		if(incr==1){
		Nt=int(200.*pow(nbsphere,-1./4));
		Np=int(200.*pow(nbsphere,-1./4));
		eval_sphere();
		}


			for (i=1;i<=nbparoi;i++)
			{
			fscanf(P_FICHIER, "%f %f %f\n", &qcjx[i], &qcjy[i], &qcjz[i]);
			fscanf(P_FICHIER, "%f %f %f\n", &qckx[i], &qcky[i], &qckz[i]);
			fscanf(P_FICHIER, "%f %f %f\n", &qclx[i], &qcly[i], &qclz[i]);
			fscanf(P_FICHIER, "%f %f %f\n", &qcmx[i], &qcmy[i], &qcmz[i]);

			qcjx[i]=qcjx[i]*X_TOT/X_TOT2;
			qcjy[i]=qcjy[i]*X_TOT/X_TOT2; 
			qcjz[i]=qcjz[i]*X_TOT/X_TOT2;

			qckx[i]=qckx[i]*X_TOT/X_TOT2;
			qcky[i]=qcky[i]*X_TOT/X_TOT2; 
			qckz[i]=qckz[i]*X_TOT/X_TOT2;

			qclx[i]=qclx[i]*X_TOT/X_TOT2;
			qcly[i]=qcly[i]*X_TOT/X_TOT2; 
			qclz[i]=qclz[i]*X_TOT/X_TOT2;

			qcmx[i]=qcmx[i]*X_TOT/X_TOT2;
			qcmy[i]=qcmy[i]*X_TOT/X_TOT2; 
			qcmz[i]=qcmz[i]*X_TOT/X_TOT2;
			}

		fscanf(P_FICHIER, "%f %f\n", &dep1max,&dep1min);
		fscanf(P_FICHIER, "%f %f\n", &dep2max,&dep2min);  
		fscanf(P_FICHIER, "%f %f\n", &dep3max,&dep3min);     

		glBindBuffer(GL_ARRAY_BUFFER, vboId);

		int nsize0=((2*Np+2)*Nt);
		int nsize=3*((2*Np+2)*Nt)*nbsphere;
		int nsize2=2*3*((2*Np+2)*Nt)*nbsphere;

		glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*nsize, NULL, GL_DYNAMIC_DRAW);
		vertices = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0,3*sizeof(float)*nsize,GL_MAP_WRITE_BIT));

		int nump=0;
		int edgei;
		float xxi,yyi,zzi,rays;
		float dpl1,dpl2,dpl3;
		int nump3,numo3,cnump3;

		if(ltype==0){
			for (i=1;i<=nbsphere;i++)
			{
			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &dpl1, &dpl2, &dpl3);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				//vertices[cnump3]=0;
				//vertices[cnump3+1]=0;
				//vertices[cnump3+2]=1;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==1){

		invmaxv = ((dep1max-dep1min)<=1.e-12)?0.:(1./(dep1max-dep1min)) ;

			for (i=1;i<=nbsphere;i++)
			{
			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &dpl1, &dpl2, &dpl3);

			vali=invmaxv*(dep1max-dpl1);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==2){

		invmaxv = ((dep2max-dep2min)<=1.e-12)?0.:(1./(dep2max-dep2min)) ;

			for (i=1;i<=nbsphere;i++)
			{
			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &dpl1, &dpl2, &dpl3);

			vali=invmaxv*(dep2max-dpl2);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==3){

		invmaxv = ((dep3max-dep3min)<=1.e-12)?0.:(1./(dep3max-dep3min)) ;

			for (i=1;i<=nbsphere;i++)
			{
			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &dpl1, &dpl2, &dpl3);

			vali=invmaxv*(dep3max-dpl3);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}


		//glBindBuffer(GL_ARRAY_BUFFER, 0);

		NNp=Np;
		NNt=Nt;

     }
     else if (ntype==1){
     fscanf(P_FICHIER, "%i %i %f \n", &nbsphere, &nbparoi, &tps);	
     
     if(incr==1){
	Nt=int(200.*pow(nbsphere,-1./4));
	Np=int(200.*pow(nbsphere,-1./4));
		eval_sphere();
     }


     fscanf(P_FICHIER, "%f %f\n", &svmmax,&svmmin);
     fscanf(P_FICHIER, "%f %f\n", &tracemax,&tracemin);
     fscanf(P_FICHIER, "%f %f\n", &sig11max,&sig11min);
     fscanf(P_FICHIER, "%f %f\n", &sig12max,&sig12min);
     fscanf(P_FICHIER, "%f %f\n", &sig13max,&sig13min);          
     fscanf(P_FICHIER, "%f %f\n", &sig22max,&sig22min);   
     fscanf(P_FICHIER, "%f %f\n", &sig23max,&sig23min);                    
     fscanf(P_FICHIER, "%f %f\n", &sig33max,&sig33min);     
     fscanf(P_FICHIER, "%f %f\n", &sig1max,&sig1min);
     fscanf(P_FICHIER, "%f %f\n", &sig2max,&sig2min);  
     fscanf(P_FICHIER, "%f %f\n", &sig3max,&sig3min);         
          
    tracemin=0.;

		glBindBuffer(GL_ARRAY_BUFFER, vboId);

		int nsize0=((2*Np+2)*Nt);
		int nsize=3*((2*Np+2)*Nt)*nbsphere;
		int nsize2=2*3*((2*Np+2)*Nt)*nbsphere;

		glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*nsize, NULL, GL_DYNAMIC_DRAW);
		vertices = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0,3*sizeof(float)*nsize,GL_MAP_WRITE_BIT));

		int nump=0;
		int edgei;
		float xxi,yyi,zzi,rays;
		int nump3,numo3,cnump3;
		float svm,trac,s11,s12,s13,s22,s23,s33,s1,s2,s3;

		if(ltype==0){
cout<<"hello! "<<nbsphere<<endl;
			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				//vertices[cnump3]=0;
				//vertices[cnump3+1]=0;
				//vertices[cnump3+2]=1;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==1){

		invmaxv = ((sig11max-sig11min)<=1.e-12)?0.:(1./(sig11max-sig11min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig11max-s11);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==2){

		invmaxv = ((sig22max-sig22min)<=1.e-12)?0.:(1./(sig22max-sig22min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig22max-s22);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==3){

		invmaxv = ((sig33max-sig33min)<=1.e-12)?0.:(1./(sig33max-sig33min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig33max-s33);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==4){

		invmaxv = ((sig1max-sig1min)<=1.e-12)?0.:(1./(sig1max-sig1min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig1max-s1);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==5){

		invmaxv = ((sig2max-sig2min)<=1.e-12)?0.:(1./(sig2max-sig2min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig2max-s2);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==6){

		invmaxv = ((sig3max-sig3min)<=1.e-12)?0.:(1./(sig3max-sig3min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig3max-s3);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==7){

		invmaxv = ((sig12max-sig12min)<=1.e-12)?0.:(1./(sig12max-sig12min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig12max-s12);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==8){

		invmaxv = ((sig13max-sig13min)<=1.e-12)?0.:(1./(sig13max-sig13min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig13max-s13);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==9){

		invmaxv = ((sig23max-sig23min)<=1.e-12)?0.:(1./(sig23max-sig23min)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(sig23max-s23);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==10){

	         invmaxv = ((tracemax-tracemin)<=1.e-12)?0.:(1./(tracemax-tracemin)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(tracemax-trac);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==11){

	         invmaxv = ((svmmax-svmmin)<=1.e-12)?0.:(1./(svmmax-svmmin)) ;

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &svm, &trac, &s11, &s12, &s13, &s22, &s23, &s33, &s1, &s2, &s3);

			vali=invmaxv*(svmmax-svm);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.;
				vertices[cnump3+1]=vali;
				vertices[cnump3+2]=0.;

				}
			nump+=3*nsize0;

			}
		}

		//glBindBuffer(GL_ARRAY_BUFFER, 0);

		NNp=Np;
		NNt=Nt;
     
     }    
     else if (ntype==2){


		fscanf(P_FICHIER, "%i %i %f \n", &nbsphere, &nbparoi, &tps);

		if(incr==1){
		Nt=int(200.*pow(nbsphere,-1./4));
		Np=int(200.*pow(nbsphere,-1./4));
		eval_sphere();
		}


			for (i=1;i<=nbparoi;i++)
			{
			fscanf(P_FICHIER, "%f %f %f\n", &qcjx[i], &qcjy[i], &qcjz[i]);
			fscanf(P_FICHIER, "%f %f %f\n", &qckx[i], &qcky[i], &qckz[i]);
			fscanf(P_FICHIER, "%f %f %f\n", &qclx[i], &qcly[i], &qclz[i]);
			fscanf(P_FICHIER, "%f %f %f\n", &qcmx[i], &qcmy[i], &qcmz[i]);

			qcjx[i]=qcjx[i]*X_TOT/X_TOT2;
			qcjy[i]=qcjy[i]*X_TOT/X_TOT2; 
			qcjz[i]=qcjz[i]*X_TOT/X_TOT2;

			qckx[i]=qckx[i]*X_TOT/X_TOT2;
			qcky[i]=qcky[i]*X_TOT/X_TOT2; 
			qckz[i]=qckz[i]*X_TOT/X_TOT2;

			qclx[i]=qclx[i]*X_TOT/X_TOT2;
			qcly[i]=qcly[i]*X_TOT/X_TOT2; 
			qclz[i]=qclz[i]*X_TOT/X_TOT2;

			qcmx[i]=qcmx[i]*X_TOT/X_TOT2;
			qcmy[i]=qcmy[i]*X_TOT/X_TOT2; 
			qcmz[i]=qcmz[i]*X_TOT/X_TOT2;
			}
 

		glBindBuffer(GL_ARRAY_BUFFER, vboId);

		int nsize0=((2*Np+2)*Nt);
		int nsize=3*((2*Np+2)*Nt)*nbsphere;
		int nsize2=2*3*((2*Np+2)*Nt)*nbsphere;

		glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*nsize, NULL, GL_DYNAMIC_DRAW);
		vertices = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0,3*sizeof(float)*nsize,GL_MAP_WRITE_BIT));

		int nump=0;
		int edgei;
		float xxi,yyi,zzi,rays;
		float typs,color1,color2,color3;
		int nump3,numo3,cnump3;

			for (i=1;i<=nbsphere;i++)
			{
			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f\n", &edgei, &xxi, &yyi, &zzi, &rays, &typs);

			if(typs<2) {color1=0.;color2=0.;color3=1.;}else{
			color1=0.5;color2=0.5;color3=0.5;
			}

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=color1;
				vertices[cnump3+1]=color2;
				vertices[cnump3+2]=color3;

				}
			nump+=3*nsize0;

			}

		//glBindBuffer(GL_ARRAY_BUFFER, 0);

		NNp=Np;
		NNt=Nt;
	
     
     }      
     else if (ntype==3){
	 fscanf(P_FICHIER, "%i %i %f \n", &NCONT, &nbparoi, &tps);
     
           longmoy=0.;

	     for (i=1;i<=NCONT;i++)
	     {

	      fscanf(P_FICHIER, "%f %f %f %f %f %f %f\n", &vcoh[i][1], &vcoh[i][2], &vcoh[i][3], &vcoh[i][4], &vcoh[i][5], &vcoh[i][6], &vcoh[i][7]);
	      vcoh[i][1]=vcoh[i][1]*X_TOT/X_TOT2;
	      vcoh[i][2]=vcoh[i][2]*X_TOT/X_TOT2;
	      vcoh[i][3]=vcoh[i][3]*X_TOT/X_TOT2;
	      vcoh[i][4]=vcoh[i][4]*X_TOT/X_TOT2;
	      vcoh[i][5]=vcoh[i][5]*X_TOT/X_TOT2;
	      vcoh[i][6]=vcoh[i][6]*X_TOT/X_TOT2;
	      vcoh[i][7]=vcoh[i][7]*X_TOT/X_TOT2;
		longmoy+=vcoh[i][7];
	     }

	longmoy/=NCONT;
     
     }else if (ntype==4){
     fscanf(P_FICHIER, "%i %i %f \n", &nbsphere, &nbparoi, &tps);	
     

     if(incr==1){
	Nt=int(200.*pow(nbsphere,-1./4));
	Np=int(200.*pow(nbsphere,-1./4));
		eval_sphere();
     }

     fscanf(P_FICHIER, "%f %f\n", &def11max,&def11min);        
     fscanf(P_FICHIER, "%f %f\n", &def22max,&def22min);                    
     fscanf(P_FICHIER, "%f %f\n", &def33max,&def33min);    
 
     fscanf(P_FICHIER, "%f %f\n", &defe11max,&defe11min);        
     fscanf(P_FICHIER, "%f %f\n", &defe22max,&defe22min);                    
     fscanf(P_FICHIER, "%f %f\n", &defe33max,&defe33min);  

		glBindBuffer(GL_ARRAY_BUFFER, vboId);

		int nsize0=((2*Np+2)*Nt);
		int nsize=3*((2*Np+2)*Nt)*nbsphere;
		int nsize2=2*3*((2*Np+2)*Nt)*nbsphere;

		glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*nsize, NULL, GL_DYNAMIC_DRAW);
		vertices = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0,3*sizeof(float)*nsize,GL_MAP_WRITE_BIT));

		int nump=0;
		int edgei;
		float xxi,yyi,zzi,rays;
		int nump3,numo3,cnump3;
		float def11,def22,def33,defe11,defe22,defe33;

		if(ltype==0){

			for (i=1;i<=nbsphere;i++)
			{
     fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f \n", &edgei, &xxi, &yyi, &zzi, &rays, &def11, &def22, &def33, &defe11, &defe22, &defe33);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				//vertices[cnump3]=0;
				//vertices[cnump3+1]=0;
				//vertices[cnump3+2]=1;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==1){

			invmaxv = ((def11max-def11min)<=1.e-12)?0.:(1./(def11max-def11min)) ;

			for (i=1;i<=nbsphere;i++)
			{

			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f \n", &edgei, &xxi, &yyi, &zzi, &rays, &def11, &def22, &def33, &defe11, &defe22, &defe33);
			
			vali=invmaxv*(def11max-def11);
			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==2){

			invmaxv = ((def22max-def22min)<=1.e-12)?0.:(1./(def22max-def22min)) ;

			for (i=1;i<=nbsphere;i++)
			{

			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f \n", &edgei, &xxi, &yyi, &zzi, &rays, &def11, &def22, &def33, &defe11, &defe22, &defe33);
			
			vali=invmaxv*(def22max-def22);
			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==3){

			invmaxv = ((def33max-def33min)<=1.e-12)?0.:(1./(def33max-def33min)) ;

			for (i=1;i<=nbsphere;i++)
			{

			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f \n", &edgei, &xxi, &yyi, &zzi, &rays, &def11, &def22, &def33, &defe11, &defe22, &defe33);
			
			vali=invmaxv*(def33max-def33);
			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==4){

			invmaxv = ((defe11max-defe11min)<=1.e-12)?0.:(1./(defe11max-defe11min)) ;

			for (i=1;i<=nbsphere;i++)
			{

			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f \n", &edgei, &xxi, &yyi, &zzi, &rays, &def11, &def22, &def33, &defe11, &defe22, &defe33);
			
			vali=invmaxv*(defe11max-defe11);
			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==5){

			invmaxv = ((defe22max-defe22min)<=1.e-12)?0.:(1./(defe22max-defe22min)) ;

			for (i=1;i<=nbsphere;i++)
			{

			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f \n", &edgei, &xxi, &yyi, &zzi, &rays, &def11, &def22, &def33, &defe11, &defe22, &defe33);
			
			vali=invmaxv*(defe22max-defe22);
			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==6){

			invmaxv = ((defe33max-defe33min)<=1.e-12)?0.:(1./(defe33max-defe33min)) ;

			for (i=1;i<=nbsphere;i++)
			{

			fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %f %f \n", &edgei, &xxi, &yyi, &zzi, &rays, &def11, &def22, &def33, &defe11, &defe22, &defe33);
			
			vali=invmaxv*(defe33max-defe33);
			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

		//	memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=1.-vali;
				vertices[cnump3+1]=0.;
				vertices[cnump3+2]=vali;

				}
			nump+=3*nsize0;

			}
		}

		//glBindBuffer(GL_ARRAY_BUFFER, 0);

		NNp=Np;
		NNt=Nt;


     }else if (ntype==5){
     fscanf(P_FICHIER, "%i %i %f \n", &nbsphere, &nbparoi, &tps);

     if(incr==1){
	Nt=int(200.*pow(nbsphere,-1./4));
	Np=int(200.*pow(nbsphere,-1./4));
		eval_sphere();
     }

     fscanf(P_FICHIER, "%f %f %f %f %f %f %f %f\n", &mint, &maxt, &minfx, &maxfx, &minfy, &maxfy, &minfz, &maxfz);
    

		glBindBuffer(GL_ARRAY_BUFFER, vboId);

		int nsize0=((2*Np+2)*Nt);
		int nsize=3*((2*Np+2)*Nt)*nbsphere;
		int nsize2=2*3*((2*Np+2)*Nt)*nbsphere;

		glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*nsize, NULL, GL_DYNAMIC_DRAW);
		vertices = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0,3*sizeof(float)*nsize,GL_MAP_WRITE_BIT));

		int nump=0;
		int edgei;
		float xxi,yyi,zzi,rays;
		int nump3,numo3,cnump3,ph;
		float temp,flx,fly,flz;
		float color1,color2,color3;

		if(ltype==0){

			for (i=1;i<=nbsphere;i++)
			{
      fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %i\n", &edgei, &xxi, &yyi, &zzi, &rays, &temp, &flx, &fly, &flz, &ph);

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				//vertices[cnump3]=0;
				//vertices[cnump3+1]=0;
				//vertices[cnump3+2]=1;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==1){

		invmaxv = ((maxt-mint)<=1.e-12)?0.:(1./(maxt-mint)) ;

			for (i=1;i<=nbsphere;i++)
			{
      fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %i\n", &edgei, &xxi, &yyi, &zzi, &rays, &temp, &flx, &fly, &flz, &ph);

			vali=4.*(1.-invmaxv*(maxt-temp));
			if(vali<=1){color1=0.;color2=vali;color3=1.;}
			else if(vali<=2){color1=0.;color2=1.;color3=2.-vali;}
			else if(vali<=3){color1=vali-2.;color2=1.;color3=0.;}
			else if(vali<=4){color1=1.;color2=4.-vali;color3=0.;}

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			//memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=color1;
				vertices[cnump3+1]=color2;
				vertices[cnump3+2]=color3;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==2){

 		 invmaxv = ((maxfx-minfx)<=1.e-12)?0.:(1./(maxfx-minfx)) ;

			for (i=1;i<=nbsphere;i++)
			{
      fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %i\n", &edgei, &xxi, &yyi, &zzi, &rays, &temp, &flx, &fly, &flz, &ph);

			vali=invmaxv*(maxfx-flx);
			if(vali<=1){color1=0.;color2=vali;color3=1.;}
			else if(vali<=2){color1=0.;color2=1.;color3=2.-vali;}
			else if(vali<=3){color1=vali-2.;color2=1.;color3=0.;}
			else if(vali<=4){color1=1.;color2=4.-vali;color3=0.;}

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			//memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=color1;
				vertices[cnump3+1]=color2;
				vertices[cnump3+2]=color3;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==3){

 		 invmaxv = ((maxfy-minfy)<=1.e-12)?0.:(1./(maxfy-minfy)) ;

			for (i=1;i<=nbsphere;i++)
			{
      fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %i\n", &edgei, &xxi, &yyi, &zzi, &rays, &temp, &flx, &fly, &flz, &ph);

			vali=invmaxv*(maxfy-fly);
			if(vali<=1){color1=0.;color2=vali;color3=1.;}
			else if(vali<=2){color1=0.;color2=1.;color3=2.-vali;}
			else if(vali<=3){color1=vali-2.;color2=1.;color3=0.;}
			else if(vali<=4){color1=1.;color2=4.-vali;color3=0.;}

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			//memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=color1;
				vertices[cnump3+1]=color2;
				vertices[cnump3+2]=color3;

				}
			nump+=3*nsize0;

			}
		}else if(ltype==4){

 		 invmaxv = ((maxfz-minfz)<=1.e-12)?0.:(1./(maxfz-minfz)) ;

			for (i=1;i<=nbsphere;i++)
			{
      fscanf(P_FICHIER, "%i %f %f %f %f %f %f %f %f %i\n", &edgei, &xxi, &yyi, &zzi, &rays, &temp, &flx, &fly, &flz, &ph);

			vali=invmaxv*(maxfz-flz);
			if(vali<=1){color1=0.;color2=vali;color3=1.;}
			else if(vali<=2){color1=0.;color2=1.;color3=2.-vali;}
			else if(vali<=3){color1=vali-2.;color2=1.;color3=0.;}
			else if(vali<=4){color1=1.;color2=4.-vali;color3=0.;}

			xxi=xxi*X_TOT/X_TOT2-X_TOT/2.;
			yyi=yyi*X_TOT/X_TOT2-Y_TOT/2.;
			zzi=zzi*X_TOT/X_TOT2-Z_TOT/2.;
			rays=rays*X_TOT/X_TOT2;

			//memcpy(vertices+(nsize+nump),views,3*nsize0*sizeof(float));  

				for (int numo = 0; numo<nsize0; numo++){
				nump3=nump+3*numo;numo3=3*numo;cnump3=nsize2+nump3;
				vertices[nump3]=xxi+rays*views[numo3];
				vertices[nump3+1]=yyi+rays*views[numo3+1];
				vertices[nump3+2]=zzi+rays*views[numo3+2]; 
				vertices[cnump3]=color1;
				vertices[cnump3+1]=color2;
				vertices[cnump3+2]=color3;

				}
			nump+=3*nsize0;

			}
		}






		//glBindBuffer(GL_ARRAY_BUFFER, 0);

		NNp=Np;
		NNt=Nt;

 }           
           
   
  } 
  else
  {
    fclose(P_FICHIER);
   //  g.close();
     incr = 0 ;
   
    P_FICHIER = fopen(NOM_FICHIER, "r");  
   //  g.open(NOM_FICHIER,fstream::in);  
   
   fscanf(P_FICHIER, "%f %f %f %f %f %f \n", &qmaxv,&qminv,&qmaxh,&qminh,&qmaxz,&qminz);
   //  g>>qmaxv>>qminv>>qmaxh>>qminh>>qmaxz>>qminz;
   nbsphere=0;
  glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*0, NULL, GL_DYNAMIC_DRAW);
  vertices = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0,3*sizeof(float)*0,GL_MAP_WRITE_BIT));

    lanime = 0 ;
    tps=0.;
  }



float tcpu2=(((float) clock())/CLOCKS_PER_SEC)-tcpu1;
cout<<"tcpu_lect:"<<tcpu2<<endl;
}


void reshape(int x, int y)
{

float w = (float) x;
float h = (float) y;	
	
glViewport(0, 0, (GLsizei) x, (GLsizei) y);
glMatrixMode(GL_PROJECTION);
glLoadIdentity();
     
gluPerspective(
angview,
1,
1.*(dx*X_TOT*(1.+pcut)+dy*Y_TOT*(1.+pcut)+dz*Z_TOT*(1.+pcut)),
10.*(dx*X_TOT*(1.+pcut)+dy*Y_TOT*(1.+pcut)+dz*Z_TOT*(1.+pcut))
);	//Pour les explications, lire le tutorial sur OGL et win
  
  	glMatrixMode(GL_MODELVIEW); 	//Optionnel


}


void writefile(int format, int sort, int options, int nbcol,
               char *filename, const char *extension){
       
FILE *fp;
  char file[256];
  int state = GL2PS_OVERFLOW, buffsize = 0;
  GLint viewport[4];

  strcpy(file, filename);
  strcat(file, ".");
  strcat(file, extension);

  viewport[0] = 0;
  viewport[1] = 0;
  viewport[2] = window_w;
  viewport[3] = window_h;
 
  fp = fopen(file, "wb");

  if(!fp){
    printf("Unable to open file %s for writing\n", file);
    exit(1);
  }

  printf("Saving image to file %s... ", file);
  fflush(stdout);

  while(state == GL2PS_OVERFLOW){
    buffsize += 1024*1024;
    gl2psBeginPage(file, "test", viewport, format, sort, options,
                   GL_RGBA, 0, NULL, nbcol, nbcol, nbcol, 
                   buffsize, fp, file);
    affichage();
    state = gl2psEndPage();
  }

  fclose(fp);

  printf("Done!\n");
  fflush(stdout);       
}
	
     
	void mouse(int button, int state,int x,int y)
	{
	switch(button){
	case GLUT_LEFT_BUTTON:
	if(state==GLUT_DOWN){
	click_g = 1;
	x_g     = x;
	y_g     = window_h-y;
	}
	break;
	case GLUT_MIDDLE_BUTTON:
	break;
	default:
	if(state==GLUT_DOWN){
	click_d = 1;
	x_d     = x;
	y_d     = window_h-y;
	
	theta_eye+= (3.14159/12.*(((float) x_d)-(((float) window_h)/2.))/(((float) window_h)/2.));
	phi_eye-= (3.14159/12.*(((float) y_d)-(((float) window_h)/2.))/(((float) window_h)/2.));
	if(theta_eye>3.14159) theta_eye-=(2.*3.14159);
	if(phi_eye>3.14159) phi_eye=(2.*3.14159)-phi_eye;
	if(theta_eye<-3.14159) theta_eye+=(2.*3.14159); 
	if(phi_eye<0) phi_eye=-phi_eye ;

 	x_eye = dx*rayon*sin(phi_eye)*cos(theta_eye)-dy*rayon*sin(phi_eye)*sin(theta_eye)+dz*rayon*cos(phi_eye);    
	y_eye = dx*rayon*cos(phi_eye)+dy*rayon*sin(phi_eye)*cos(theta_eye)-dz*rayon*sin(phi_eye)*sin(theta_eye);
	z_eye = -dx*rayon*sin(phi_eye)*sin(theta_eye)+dy*rayon*cos(phi_eye)+dz*rayon*sin(phi_eye)*cos(theta_eye);

	  
	}  
	break;
	}
	}       



  void mousemotion(int x,int y)
  {

  }

