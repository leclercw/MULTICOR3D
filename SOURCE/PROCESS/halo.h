#include "../MAIN/def.h"

#ifndef __HALO__
#define __HALO__

void treat_halo(unsigned int ** NOHALO, int * NBHALO, R & RHALO, int SZHALO, R R_SPH, int * LIST_H,  int vecsizeh, int vecsizexh, int vecsizeyh, int vecsizezh, vector< vector<int> > & coulh,vector< vector<int> > & numch, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, int NB_SPH,R * LIST_V, R * VOLHALO, R coef1,int NMAXHALO);
void treat_halo2(bool * LIST_P, unsigned int ** NOHALO, int * NBHALO, R & RHALO, int SZHALO, R R_SPH, int * LIST_H, int vecsizeh, int vecsizexh, int vecsizeyh, int vecsizezh, vector< vector<int> > & coulh,vector< vector<int> > & numch, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, int NB_SPH,R * LIST_V, R * VOLHALO, R coef1,int NMAXHALO);
void explore_halo(int iglo, unsigned int ** NOCONT, int * NBCONTCO, int ** CONT, int * NOLOC, int * PHALO);
void search_halo(int nohalo, int NB_SPH, unsigned int ** NOHALO, int * NBHALO, int NMAXHALO, R * LIST_V, R * VOLHALO, R coef1, unsigned int ** NOCONT, int * NBCONTCO, int ** CONT);
void retreat_halo(int NB_SPH, unsigned int ** NOHALO, int * NBHALO, int NMAXHALO, R * LIST_V, R * VOLHALO, R coef1, unsigned int ** NOCONT, int * NBCONTCO, int ** CONT, vector<int> list_rupt);

#endif
