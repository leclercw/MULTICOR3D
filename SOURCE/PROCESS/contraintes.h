#include "../MAIN/def.h"

#ifndef __CONTRAINTES__
#define __CONTRAINTES__

void contraintesloc(R Pi, R coef1, int NBENREG, R H_TOT, R V_TOT, R Z_TOT, R * LIST_R, R ** FCJI, R ** NCONT,int ite, int NB_SPH, R * VONMIS, R * TRACE, R * SIG11, R * SIG12, R * SIG13, R * SIG22, R * SIG23, R * SIG33, R * SIG1, R * SIG2, R * SIG3,R &minvm, R &maxvm,R &mintrac, R &maxtrac, R &minsig11, R &maxsig11, R &minsig12, R &maxsig12, R &minsig13, R &maxsig13, R &minsig22, R &maxsig22, R &minsig23, R &maxsig23, R &minsig33, R &maxsig33, R &minsig1, R &maxsig1, R &minsig2, R &maxsig2, R &minsig3, R &maxsig3, unsigned int ** NOCONT, int * NBCONTCO,bool * EDGE, R * LIST_V, R * LIST_IND);

void contrainteshalo(R Pi, R coef1, int ite,int NBENREG, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT,R * LIST_R, R ** FCJI, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * VONMIS, R * TRACE, R * SIG11, R * SIG12, R * SIG13, R * SIG22, R * SIG23, R * SIG33, R * SIG1, R * SIG2, R * SIG3,R &minvm, R &maxvm,R &mintrac, R &maxtrac, R &minsig11, R &maxsig11, R &minsig12, R &maxsig12, R &minsig13, R &maxsig13, R &minsig22, R &maxsig22, R &minsig23, R &maxsig23, R &minsig33, R &maxsig33, R &minsig1, R &maxsig1, R &minsig2, R &maxsig2, R &minsig3, R &maxsig3,bool * EDGE, unsigned int ** NOHALO, int * NBHALO, R * LIST_V, R * VOLHALO, R * LIST_IND);

void defoloc(R Pi, R coef1, R coef2, int NBENREG, R H_TOT, R V_TOT, R Z_TOT, int ite, int NB_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * DCONTXO, R * DCONTYO, R * DCONTZO, R * DCONTX, R * DCONTY, R * DCONTZ, int ** CONT, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * EPSI11, R * EPSI22, R * EPSI33, R * EPSE11, R * EPSE22, R * EPSE33, R & mindef11, R & maxdef11,R & mindef22, R & maxdef22,R & mindef33, R & maxdef33, R & mindefe11, R & maxdefe11,R & mindefe22, R & maxdefe22,R & mindefe33, R & maxdefe33, bool * EDGE, R * LIST_V, R * LIST_IND);

void defohalo(R Pi, R coef1, R coef2, int NBENREG, R H_TOT, R V_TOT, R Z_TOT, int ite, int NB_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * DCONTXO, R * DCONTYO, R * DCONTZO, R * DCONTX, R * DCONTY, R * DCONTZ, int ** CONT, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * EPSI11, R * EPSI22, R * EPSI33, R * EPSE11, R * EPSE22, R * EPSE33, R & mindef11, R & maxdef11,R & mindef22, R & maxdef22,R & mindef33, R & maxdef33, R & mindefe11, R & maxdefe11,R & mindefe22, R & maxdefe22,R & mindefe33, R & maxdefe33, bool * EDGE, unsigned int ** NOHALO, int * NBHALO, R * LIST_V, R * VOLHALO, R * LIST_IND);

void deploc(int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * DEP1, R * DEP2, R * DEP3, R & maxdep1, R & mindep1, R & maxdep2, R & mindep2, R & maxdep3, R & mindep3);

#endif
