#include "../MAIN/def.h"

#ifndef _THERMIQ_
#define _THERMIQ_

void Exp_temp_hydriq(int ite, int NB_SPH, R dt, R p, R coef1, R rhos, R Mw, R * LIST_TEMP, R * LIST_V,R * LIST_X,R * LIST_Y,R * LIST_Z, R * LIST_CONW, R H_TOT, R V_TOT, R Z_TOT,bool bstart);

void temp_hydriq(int ite, R coef2,int NB_SPH, int NBCO, R dt, R Lambda_c,R p, R Cps, R Cpw, R coef1, int ** CONT, unsigned int ** NOCONT, R ** NCONT,int * NBCONTCO,int * NBHALO,R * VOLHALO, unsigned int **NOHALO, R *LIST_M, R rhos, R Mw,R Lw, R * LIST_R, R * LIST_TEMP, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4, bool * EDGE5, bool * EDGE6, R & mint, R & maxt, R * DCONTO, R * LIST_IND, R N1, R N2, R N3, R N4, R N5, R N6, R * DeltaT, R Dw, R * LIST_V, R * LIST_CONW, R * Cp_eff, R * LIST_GCX,R * LIST_GCY,R * LIST_GCZ, R & minc, R & maxc, R * DeltaC,R A, R B, R Rh,R Rg, R k, R h, R Tov, R Cov, R H_TOT, R V_TOT, R Z_TOT) ;

void ther_load_hygr(int ite, int nitd,bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO, bool * LIST_P, R alphag1, R alphag2);

void ther_load_hygrc(int NBENREG,int ite,R Mw, R rhow, R * LIST_CONW,R Cint, bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO);

void ther_load_dila(int ite, int nitd,bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO, bool * LIST_P,R alpha1, R alpha2);

void ther_load_dilat(int ite, R Tinit, R LIST_TEMP, bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO, bool * LIST_P,R alpha1, R alpha2);

void ther_load_hygrc_dilat(int NBENREG,int ite,R Mw,R rhow, R * LIST_CONW,R Cinit,R * LIST_TEMP, R Tinit,  bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO,bool * LIST_P, R alpha1, R alpha2);

void hydriq(R coef2,int NB_SPH, int NBCO, R dt, R Cov, R Clim,R k, R D_w, R coef1, int ** CONT, R *LIST_V, R * LIST_R, R * LIST_CONW, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, R & minc, R & maxc, R * DCONTO, R * LIST_IND, R N1, R N2, R N3, R N4, R N5, R N6, R * DeltaC, R H_TOT, R V_TOT, R Z_TOT, int dir, bool btype);  

void temp(R coef2,int NB_SPH, int NBCO, R dt, R Flimp, R Lambda_c, R Cp, R coef1, int ** CONT, R *LIST_M, R * LIST_R, R * LIST_TEMP,bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, R & mint, R & maxt, R * DCONTO, R * LIST_IND, R N1, R N2, R N3, R N4, R N5, R N6, R * DeltaT,int dir, bool btype);

void temp2(R coef2,int NB_SPH, int NBCO, R dt, R Flimp, R lambda1, R lambda2, R cp1, R cp2, R coef1, int ** CONT, R *LIST_M, R * LIST_R, R * LIST_TEMP,bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, R & mint, R & maxt, R * DCONTO, R * LIST_IND, bool * LIST_P, R N1, R N2, R N3, R N4, R N5, R N6,R * DeltaT, int dir, bool btype);

void flux(R coef2,int NB_SPH, int NBCO, R dt, R Lambda_eff, R coef1, int ** CONT, R ** NCONT, R *LIST_M, R * LIST_R, R * LIST_TEMP,bool * EDGE, R & minfx, R & maxfx, R & minfy, R & maxfy, R & minfz, R & maxfz, R * DCONTO, R H_TOT, R V_TOT, R Z_TOT, R * LIST_FLX,R * LIST_FLY,R * LIST_FLZ,int * NBHALO, unsigned int **NOHALO, int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO, R * LIST_V);

void flux2(R coef2,int NB_SPH, int NBCO, R dt, R lambda_eff1, R lambda_eff2, R coef1, int ** CONT, R ** NCONT, R *LIST_M, R * LIST_R, R * LIST_TEMP,bool * EDGE, R & minfx, R & maxfx, R & minfy, R & maxfy, R & minfz, R & maxfz, R * DCONTO,  bool * LIST_P, R H_TOT, R V_TOT, R Z_TOT, R * LIST_FLX,R * LIST_FLY,R * LIST_FLZ,int * NBHALO, unsigned int **NOHALO, int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO, R * LIST_V);

void grad_C(R coef2,int NB_SPH, int NBCO, R dt, R Dw, R coef1, int ** CONT, R ** NCONT, R *LIST_M, R * LIST_R, R * LIST_CONW,bool * EDGE, R & mingx, R & maxgx, R & mingy, R & maxgy, R & mingz, R & maxgz, R * DCONTO, R H_TOT, R V_TOT, R Z_TOT, R * LIST_GCX,R * LIST_GCY,R * LIST_GCZ,int * NBHALO,unsigned int **NOHALO, int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO, R * LIST_V); 

void grad_C_P_Halo(int jt,R coef2, R Dw,  int ** CONT, R ** NCONT, R * LIST_R, R * LIST_CONW, R * DCONTO, R * LIST_GCX,R * LIST_GCY,R * LIST_GCZ,int * NBHALO,unsigned int **NOHALO,int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO,R * LIST_V,R coef1);

void grad_C_P(int jt,R coef2, R Dw,  int ** CONT, R ** NCONT, R * LIST_R, R * LIST_CONW, R * DCONTO, R * LIST_GCX,R * LIST_GCY,R * LIST_GCZ,int * NBHALO, unsigned int **NOHALO, int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO,R * LIST_V,R coef1);   

R crit_time(R coord, R NBCO, R lambda, R cp, R dens, R * DCONTO, R Cs); 
R crit_time2(R coord, R NBCO, R lambda1, R lambda2, R cp1, R cp2, R dens1, R dens2, R * DCONTO, R Cs); 
R crit_time_couplage(R coord, R NBCO, R lambda, R cp, R dens, R * DCONTO, R Cs, R Dw);    

#endif



