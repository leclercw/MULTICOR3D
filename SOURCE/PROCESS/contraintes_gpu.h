#include "../MAIN/def.h"

#ifndef _CONTRAINTESGPU_
#define _CONTRAINTESGPU_

void contrainteshalo_gpu(R Pi, R coef1, int ite,int NBENREG, int NB_SPH, int NBCO, int NMAXZ,int NMAXHALO, int NMAXCONT, R H_TOT, R V_TOT, R Z_TOT,R * LIST_R, R ** FCJI, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * VONMIS, R * TRACE, R * SIG11, R * SIG12, R * SIG13, R * SIG22, R * SIG23, R * SIG33, R * SIG1, R * SIG2, R * SIG3,R &minvm, R &maxvm,R &mintrac, R &maxtrac, R &minsig11, R &maxsig11, R &minsig12, R &maxsig12, R &minsig13, R &maxsig13, R &minsig22, R &maxsig22, R &minsig23, R &maxsig23, R &minsig33, R &maxsig33, R &minsig1, R &maxsig1, R &minsig2, R &maxsig2, R &minsig3, R &maxsig3, bool * EDGE, unsigned int ** NOHALO, int * NBHALO, R * LIST_V, R * VOLHALO,R * LIST_IND);

#endif

