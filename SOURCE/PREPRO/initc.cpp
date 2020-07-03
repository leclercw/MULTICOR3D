#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <time.h> 
#include <sys/time.h> 
#include <sys/resource.h> 
#include <string.h>
#include <map>
#include <cassert>
#include <vector>
#include <limits>

using namespace std;

#include "initc.h"
#include "init.h"
#include "../PROCESS/selco.h"
#include "omp.h"



void initc_sph(R coord, bool * TYPCO, int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc,R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

init_sph(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}

void initc_sphp(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

init_sphp(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}

}

}

void initc_sphp_cis(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi,int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc,R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,bool * EDGEC1, bool * EDGEC2,unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{
epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

init_sphp_cis(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_TX,LIST_TY,LIST_TZ,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,EDGEC1,EDGEC2,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);

coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}

}

}



void initc_sph_ind(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

init_sph_ind(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}

void initc_sphp_ind(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

init_sphp_ind(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}

void initc_sph_halo(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,  R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

coulh.clear();
for (int i = 0; i < vecsizeh; i++) {
    coulh.push_back(vector<int>()); // Add an empty row
}

init_sph_halo(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,vecsizexh,vecsizeyh,vecsizezh,coulh,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_H,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}

void initc_sph_halo_cyl(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,  R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,bool * EDGEC, unsigned int ** NOCONT,int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

coulh.clear();
for (int i = 0; i < vecsizeh; i++) {
    coulh.push_back(vector<int>()); // Add an empty row
}

init_sph_halo_cyl(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,vecsizexh,vecsizeyh,vecsizezh,coulh,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_H,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,EDGEC,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}



void initc_sph_halo2(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens1, R dens2, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO,bool * LIST_P, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

coulh.clear();
for (int i = 0; i < vecsizeh; i++) {
    coulh.push_back(vector<int>()); // Add an empty row
}

init_sph_halo2(TYPCO,NMAXCONT,Pi,dens1,dens2,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,vecsizexh,vecsizeyh,vecsizezh,coulh,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_H,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_P,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}

void initc_sphp_halo(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA,R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

coulh.clear();
for (int i = 0; i < vecsizeh; i++) {
    coulh.push_back(vector<int>()); // Add an empty row
}

init_sphp_halo(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,vecsizexh,vecsizeyh,vecsizezh,coulh,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_H,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}


void initc_sphp_halo_car(R coord,R EE, bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA,R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

coulh.clear();
for (int i = 0; i < vecsizeh; i++) {
    coulh.push_back(vector<int>()); // Add an empty row
}

init_sphp_halo_car(EE,TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,vecsizexh,vecsizeyh,vecsizezh,coulh,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_H,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}


void initc_sph2(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens1, R dens2, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,  R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO, bool * LIST_P, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

init_sph2(TYPCO,NMAXCONT,Pi,dens1,dens2,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBCONTCO,LIST_P,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}
}

}


void initc_sphp_cyl(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc, R & PR_SPH,R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, bool * EDGEC, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ)
{

epsi   = 0.1; //0.5 .. 1.75
R epsim  = 0.;
R epsip  = 1.;
R coord1 = 0.;

while(abs(coord1-coord)>0.01){

init_sphp_cyl(TYPCO,NMAXCONT,Pi,dens,epsi,vecsizex,vecsizey,vecsizez,coul,nocoul,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,EDGEC,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6);
selco_corci_ori(epsi,NB_SPH,NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,LIST_R,LIST_X,LIST_Y,LIST_Z,NOCONT,NBCONTCO,NMAXCONT,NMAXZ);
coord1= (R(NBCO)*2./R(NB_SPH));
if((coord1>coord)&&(abs(coord1-coord)>0.01)) {epsip=epsi;epsi=(epsi+epsim)/2.;}
if((coord1<coord)&&(abs(coord1-coord)>0.01)) {epsim=epsi;epsi=(epsi+epsip)/2.;}

}


}
