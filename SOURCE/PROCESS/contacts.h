#include "../MAIN/def.h"

#ifndef __CONTACTS__
#define __CONTACTS__

void force(R dt, int NBCO, int NB_SPH, int ** CONT,  R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONT,R * DCONTO, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI,R ** MTJI,R ** MTIJ, R E, R nu, R fric, R &fpar,R &fpar2, int NBCOP, int ** CONTP, R * DCONTP, R ** NCONTP, R * LIST_PM, R ** LIST_PVX, R ** LIST_PVY, R ** LIST_PVZ);

void forcecoh_ind(R & Epe, R & Epa,R dt, bool * TYPCO,int NBCO,  int ** CONT,  R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONT, R * DCONTO, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI,R ** MTJI,R ** MTIJ, R E1, R nu1,R E2, R nu2, R fric, int NB_SPH);

void forcep(int NB_SPH, R & fpar,R & fpar2,R dt, int NBCOP, int ** CONTP, R * DCONTP, R ** NCONTP, R * LIST_M,R * LIST_PM, R ** LIST_PVX, R ** LIST_PVY, R ** LIST_PVZ, R * LIST_R, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R E, R nu, R fric);

void forcecoh_qs(int ite, R & Epe, R & Epa,int NBCO, int NB_SPH, int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ, R amort);

void forcecoh_qs_test(int ite, R & Epe, R & Epa,int NBCO, int NB_SPH, int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ, R amort);

void forcecoh(R & Epe, R & Epa,R dt, bool * TYPCO,int NBCO,int NB_SPH,  int ** CONT,  R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONT, R * DCONTO, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ, R E, R nu, R fric);

void forcecoh2(R Kcon, R & Epe, R & Epa,R dt, bool * TYPCO,int NBCO, int NB_SPH,  int ** CONT,  R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONT, R * DCONTO, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI,R ** FOJI,R ** MTJI,R ** MTIJ, R E1, R nu1, R E2, R nu2, R fric, bool * LIST_P);

void forcecoh0(R & Epe, R & Epa,bool * TYPCO,int NBCO, int NB_SPH, int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ);

void forcecoh_qs_incr(int ite, R dt, R & Epe, R & Epa,int NBCO, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_XA, R * LIST_YA, R * LIST_ZA,R * LIST_VXA, R * LIST_VYA, R * LIST_VZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ,R * FIX, R * FIY, R * FIZ, R * MTIX, R * MTIY, R * MTIZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ);

void forcecoh0_incr(int ite, R dt, R & Epe, R & Epa, bool * TYPCO,int NBCO, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_XA, R * LIST_YA, R * LIST_ZA,R * LIST_VXA, R * LIST_VYA, R * LIST_VZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ,R * FIX, R * FIY, R * FIZ, R * MTIX, R * MTIY, R * MTIZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ);

void forcecoh2_qs_int_mI(int & nint,int & nsoft,R Cs, R & Epe, R & Epa,R dt, bool * TYPCO,int & nbco, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ,  R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI,R ** FOJI,R ** MTJI,R ** MTIJ, bool * LIST_P, int * LIST_B, R * vs, R * LIST_IND,R Pi,R Emu1, R rmu1, R Emu2, R rmu2, R amort, R GI,R GII, R snmax, R ssmax);

void forcecoh2_qs_int_mII(int & nint,int & nsoft,R Cs, R & Epe, R & Epa,R dt, bool * TYPCO,int & nbco, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ,  R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI,R ** FOJI,R ** MTJI,R ** MTIJ, bool * LIST_P, int * LIST_B, R * vs, R * LIST_IND,R Pi,R Emu1, R rmu1, R Emu2, R rmu2, R amort, R GI,R GII, R snmax, R ssmax);

void forcecoh2_qs_int_mm(int & nint,int & nsoft,R Cs, R & Epe, R & Epa,R dt, bool * TYPCO,int & nbco, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ,  R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI,R ** FOJI,R ** MTJI,R ** MTIJ, bool * LIST_P, int * LIST_B, R * vs, R * LIST_IND,R Pi,R Emu1, R rmu1, R Emu2, R rmu2, R amort, R GI,R GII, R snmax, R ssmax);
#endif

