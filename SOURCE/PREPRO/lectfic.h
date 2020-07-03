#include "../MAIN/def.h"

#ifndef __LECTFIC__
#define __LECTFIC__

void read_paroi(int & NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ);
void read_sph_ind(R & H_TOT, R & V_TOT, R & Z_TOT, R & H_POS, R & V_POS, R & Z_POS,int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX, R rind);
void read_sphp_ind(R & H_TOT, R & V_TOT, R & Z_TOT, R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX, R rind);
void read_sph(R & H_TOT, R & V_TOT, R & Z_TOT, R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX);
void read_sphp(R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX);
void read_sphp_car(R EE,R & H_TOT, R & V_TOT, R & Z_TOT, R & H_POS, R & V_POS, R & Z_POS,int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX);
void read_phase(R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, bool * LIST_P);
void read_phase_rupt(R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, bool * LIST_P, int * LIST_B);
void read_phase_compo(R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, bool * LIST_P);
void read_phase_dela(R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,  int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, bool * LIST_P, int * LIST_B);

#endif
