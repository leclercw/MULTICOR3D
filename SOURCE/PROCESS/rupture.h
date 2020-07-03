#include "../MAIN/def.h"

#ifndef __RUPTURE__
#define __RUPTURE__

void rupture1_bond_Rankine(int NBCO, int ** CONT, bool * TYPCO, R ** FOJI, R ** MTJI, R ** MTIJ, R * LIST_R, R rmu, R siglim, int * LIST_B, bool & brupt);
void rupture1_bond(int NBCO, int ** CONT, bool * TYPCO, R ** FOJI, R ** MTJI, R ** MTIJ, R * LIST_R, R rmu, R siglimt,R siglimcis, int * LIST_B, bool & brupt);
void rupture1a_trace(int NBCO, bool * TYPCO, R * TRACE, int ** CONT, R siglimt, R siglimc);
void rupture1b_trace(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * TRACE, R siglimt, R siglimc);
void rupture1c_trace(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * TRACE, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt, R siglimc, int * LIST_B, bool & brupt,vector<int> & list_rupt);
void rupture1_Griffith(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt, int * LIST_B, bool & brupt, int & nrt, int & nrcis, int & nrtot,vector<int> & list_rupt);
void rupture1_Griffith_incr(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt, int * LIST_B, bool & brupt, int & nrt, int & nrcis, int & nrtot,vector<int> & list_rupt,R ** FCJI);
void rupture1_Tresca(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimcis, int * LIST_B, bool & brupt,int &nrcis,int &nrtot,vector<int> & list_rupt);
void rupture1_Rankine(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt,R siglimc, int * LIST_B, bool & brupt,int &nrt,int &nrc,int &nrtot,vector<int> & list_rupt);
void rupture1_Mohr_Coulomb(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimcis, R phi, int * LIST_B, bool & brupt,int &nrcis,int &nrtot,vector<int> & list_rupt);
void rupture2a_trace(int NBCO, bool * TYPCO, R * TRACE, int ** CONT, R siglimt1, R siglimt2, R siglimti, bool * LIST_P);
void rupture2b_trace(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * TRACE,R siglimt1, R siglimt2, bool * LIST_P, int * LIST_B, bool & brupt,vector<int> & list_rupt);
void rupture2c_trace(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * TRACE, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt1, R siglimt2,bool * LIST_P, int * LIST_B, bool & brupt,vector<int> & list_rupt);
void rupture2_Rankine(int NB_SPH, int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt1, R siglimt2, R siglimc1, R siglimc2, int * LIST_B, bool & brupt,int &nrt,int &nrc,int &nrtot,vector<int> & list_rupt, bool * LIST_P);
#endif


