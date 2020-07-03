#include "../MAIN/def.h"

#ifndef __INITC__
#define __INITC__

using namespace std;

void initc_sph(R coord, bool * TYPCO, int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc, R & PR_SPH,  R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sphp_cis(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize,int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,bool * EDGEC1, bool * EDGEC2,unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sphp(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sph_ind(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sphp_ind(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sph_halo(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sph_halo_cyl(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,bool * EDGEC, unsigned int ** NOCONT,int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sph_halo2(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens1, R dens2, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO,bool * LIST_P, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sphp_halo(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA,R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sphp_halo(R coord,bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA,R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,bool * EDGEC, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sphp_halo_car(R coord,R EE, bool * TYPCO,int & NBCO,  int ** CONT, int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int ** numc, int * nonumc,int vecsizeh,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA,R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sph2(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens1, R dens2, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, unsigned int ** NOCONT,int * NBCONTCO, bool * LIST_P, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);

void initc_sphp_cyl(R coord,bool * TYPCO, int & NBCO,  int ** CONT,int NMAXCONT, R Pi, R dens, R & epsi, int vecsize, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, int ** numc, int * nonumc, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, bool * EDGEC, unsigned int ** NOCONT, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6, int NMAXZ);
#endif