#include "../MAIN/def.h"

#ifndef __CMDGKS__
#define __CMDGKS__

void ExpHalo(int NB_SPH, int * NBHALO, R * VOLHALO, int ite);

void ExpDEF(int NB_SPH, R * DEF, bool * EDGE);
void ExpSIG(int NB_SPH, R * SIG, bool * EDGE);

void ExpDEPL(R dplct, int ite,bool bstart);

void ExpREAC(R reac, int ite,bool bstart);

void ExpENERGY(R Ec, R Ep,R Epa, int ite,bool bstart);

void ExpVISU(int ite, R dt, int NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R,R * VONMIS,R * TRACE,R * SIG11,R * SIG12,R * SIG13,R * SIG22,R * SIG23,R * SIG33,R * SIG1,R * SIG2,R * SIG3,R * DEP1,R * DEP2,R * DEP3,R maxvm, R minvm, R maxtrac, R mintrac,R maxsig11, R minsig11,R maxsig12, R minsig12,R maxsig13, R minsig13,R maxsig22, R minsig22,R maxsig23, R minsig23,R maxsig33, R minsig33,R maxsig1, R minsig1,R maxsig2, R minsig2,R maxsig3, R minsig3,R maxdep1, R mindep1,R maxdep2, R mindep2,R maxdep3, R mindep3, bool bout, bool * EDGE, bool * LIST_P);

void ExpVISUT(int ite, R dt, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R, R* LIST_TEMP, R & maxt, R & mint, bool * EDGE, R & minfx, R & maxfx, R & minfy, R & maxfy, R & minfz, R & maxfz, R * LIST_FLX, R * LIST_FLY, R * LIST_FLZ, bool * LIST_P, bool bout);

void ExpVISUC(int ite, R dt, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R, R* LIST_CONW, R & maxc, R & minc, bool * EDGE, R & mingx, R & maxgx, R & mingy, R & maxgy, R & mingz, R & maxgz, R * LIST_GCX, R * LIST_GCY, R * LIST_GCZ, bool * LIST_P, bool bout);

void ExpVISU_COH(int ite, R dt, int NBCO, R H_TOT, R V_TOT, R Z_TOT, int ** CONT,R* LIST_X, R* LIST_Y, R* LIST_Z,R* LIST_R,bool bout);

void ExpVISU_rupt(int ite, R dt, int NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R, bool bout,int * LIST_B,  bool * EDGE);

void ExpVISU_def(int ite, R dt, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R, R * EPSI11, R * EPSI22, R * EPSI33, R * EPSE11, R * EPSE22, R * EPSE33, R mindef11, R maxdef11,R mindef22, R maxdef22,R mindef33, R maxdef33, R mindefe11, R maxdefe11,R mindefe22, R maxdefe22,R mindefe33, R maxdefe33, bool bout, bool * EDGE);

void ExpMeshHAMZA(R SIZEX, R SIZEY, R SIZEZ, int NB_SPH, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R);

void ExpMeshMULTICOR(R SIZEX, R SIZEY, R SIZEZ, int NB_SPH, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R);

void ExpDILA(int ite,R dpltx,R dplty,R dpltz,bool bstart);

void ExpPARAMETER(R Emoy,R nuu, int NB_SPH, int NBCO, R cpu, int ite, R dt, R coord1,bool bstart);

void ExpPARAMETER_rupt(int nrt, int nrc, int nrcis, int nrtot, int ite, bool bstart);

void ExpPARAMETER_ind(int ite, R ti, R force, int nbci, R pos, bool bstart);

void ExpPARAMETER_sig_dep(R ff, int ite, R ti, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R,R * VONMIS,R * TRACE,R * SIG11,R * SIG12,R * SIG13,R * SIG22,R * SIG23,R * SIG33,R * SIG1,R * SIG2,R * SIG3,R * DEP1,R * DEP2,R * DEP3,R maxvm, R minvm, R maxtrac, R mintrac,R maxsig11, R minsig11,R maxsig12, R minsig12,R maxsig13, R minsig13,R maxsig22, R minsig22,R maxsig23, R minsig23,R maxsig33, R minsig33,R maxsig1, R minsig1,R maxsig2, R minsig2,R maxsig3, R minsig3,R maxdep1, R mindep1,R maxdep2, R mindep2,R maxdep3, R mindep3,int * LIST_B, bool * EDGE);

void Exp_debonding(int ite, R dt, R ftot,R ftot2, int nint0, int nint, int npre, int nsoft,bool bstart);

#endif


