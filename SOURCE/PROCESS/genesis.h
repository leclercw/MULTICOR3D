#include "../MAIN/def.h"

#ifndef _GENESIS_
#define _GENESIS_

bool bincl(R R_CYL,R AA,R BB,R CC, R X1, R Y1, R Z1, R X2, R Y2, R Z2, R XX, R YY, R ZZ);
void gener3_cyl_alea(int nC,R facf,R eche, R buffer, int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P);
void gener3_cyl_alea_per(int nC,R facf,R eche,  int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P) ;
void gener3_cyl_align(int nC,R facf,R eche, R buffer, int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P,int dir) ;
void gener3_cyl_align_per(int nC,R facf,R eche,  int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P,int dir) ;
void gener3_sphere_alea(int nS,R eche, R buffer, R tol,int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P) ;
void gener3_sphere_alea_per(int nS,R eche, R tol,int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P);
void gener3_cyl_inf(int nC,R eche, R buffer, R tol, int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P,int dir);

#endif

