#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <time.h> 
#include <sys/time.h> 
#include <sys/resource.h> 
#include <string.h>
#include <map>
#include <cassert>
#include <vector>
#include <limits>
# include "omp.h"

#include "../POSTPRO/cmdgks.h"
#include "../PROCESS/contacts.h"
#include "../PROCESS/contraintes.h"
#include "../PROCESS/correct.h"
#include "../PROCESS/inicoh.h"
#include "../PREPRO/init.h"
#include "../PREPRO/initc.h"
#include "../PREPRO/lectfic.h"
#include "../PROCESS/selco.h"
#include "../PROCESS/ftstring.h"
#include "../PROCESS/ctime.h"
#include "../PROCESS/rupture.h"
#include "../PROCESS/halo.h"
#include "../PROCESS/heatmass.h"


#include "conf.h"

using namespace std;
typedef double R;

// Paramètres d'étude

const R Pi=3.14159265; //Nombre Pi
const R g=9.81;  //Constante de gravité

const R coord=7.5;  //Nombre de coordination
const R lmacro=0.02345; //Dimension réelle du domaine d'étude

const R coef1=1./0.64;
const R coef2=0.8058; /*
10000 ==> coef2=0.749;100000 ==> coef2=0.8058;150000 ==> coef2=0.8136;200000 ==> coef2=0.8185;
350000 ==> coef2=0.8244;500000 ==> coef2=0.8253;700000 ==> coef2=0.8293 ;8000000 ==> coef2=0.8496;
*/
const R Cse=0.3;  //Coefficient de sécurité (pas de temps élastique)
const R Cst=0.85;  //Coefficient de sécurité (pas de temps thermique)

const int nmax=100000000000; //Nombre d'itérations limite

const int NBENREG=500;
const int NB_THREADS=4;

const R lambda=7.;  // W/(mK) Thermal conductivity 
const R Tinit=20.; // °C Initial temperature 
const R Timp=50.;

const R Tov=50.; // °C Oven temperature
const R Po=exp(13.7-5120./(273.15+Tov));
const R Rg=8.314; //Ideal gas constant
const R Rh=0.75; //Relative humidity
const R Cov=101325.*Rh*Po/(Rg*(273.15+Tov)); // mol/m3 water concentration in the oven 

const R h=40.;  //(W/(m2K)  heat transfer coefficient
const R Dw=1.e-6;//33./(0.9*7800.); // m2/s water diffusion coefficient
const R Cinit=14217.16; // mol/m3 initial water concentration in the sample
const R Cimp=35.;

const R facm=1e12; //mass scaling
const R Mw=18.*1.e-3*facm;  // (kg/mol) water molar mass
const R Lw=2.3e6/facm;  // (J/kg) water vaporization latent heat
const R p=0.6; //sample porosity
const R Cps=830./facm; // J/(kgK) matter heat capacity
const R Cpw=4180/facm; // J/(kgK) water heat capacity 
const R Cpa=710./facm; // (J/(kg*K)) air heat capacity 
const R rhos=3580.*facm; //kg/m3 matter mass volume
const R rhow=1000.*facm; //kg/m3 matter mass volume
const R rhoa=1.2*facm;  // (kg/m3) air volume mass
const R k=h/(rhoa*Cpa);   // (m/s) water mass transfer coefficient
const R dens=rhos*(1-p);// Effective mass volume
const R A=0.008;
const R B=0.48; 

const R Emu=220.2e+9; //Module de Young micro échantillon
const R nuu=0.2; //Paramètre micro
const R rmu=0.722; //Paramètre micro

const R amort=0.; //Coefficient d'amortissement
const R alphag=1e-5;

const int SZHALO=1; //Taille halo

//////////////////////////////////////////////////////////////////////
///////// Début du programme principal
//////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

R t1=CPUtime();
R t1p=omp_get_wtime () ;

// Parallélisation
omp_set_num_threads(NB_THREADS);

// Allocation
allocat_memory();

 cout<<endl;
 cout<<endl;
 cout<<"           < DEBUT PRE-PROCESS >"<<endl;
 cout<<"________________________________________________"<<endl;
 cout<<endl;

// Initialisation
read_sph(H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH, R_SPH, LIST_R, LIST_X, LIST_Y, LIST_Z, RMAX);

// Allocation grille
allocat_grid(H_TOT,V_TOT,Z_TOT,RMAX,vecsizex,vecsizey, vecsizez, vecsize);

cout<<"Rayon par. moy/max :"<<R_SPH<<" - "<<RMAX<<endl;
init_sizeh(vecsizeh, vecsizexh,vecsizeyh,vecsizezh, coulh, numch, H_TOT, V_TOT, Z_TOT, RHALO,SZHALO,R_SPH);

// Init sph + Sélection des contacts 
initc_sph_halo(coord,TYPCO,NBCO,CONT,NMAXCONT,Pi,dens,epsi,vecsize,vecsizex,vecsizey,vecsizez,coul,nocoul,numc,nonumc,vecsizeh,vecsizexh,vecsizeyh,vecsizezh,coulh,PR_SPH,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,LIST_B,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_TX,LIST_TY,LIST_TZ,LIST_TXA,LIST_TYA,LIST_TZA,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,LIST_C,LIST_H,LIST_M,LIST_I,FX,FY,FZ,MTX,MTY,MTZ,EDGE,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NOCONT,NBCONTCO,LIST_V,LIST_IND,N1,N2,N3,N4,N5,N6,NMAXZ);

// Process
R t2=CPUtime();
R ti=0.;
int ite=1;

scale_sph(Pi,H_TOT,V_TOT,Z_TOT,H_POS,V_POS,Z_POS,NB_SPH,R_SPH,RMAX,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_XA,LIST_YA,LIST_ZA,LIST_M,LIST_I,epsi,lmacro,dens,LIST_V);

treat_halo(NOHALO,NBHALO,RHALO,SZHALO,R_SPH,LIST_H,vecsizeh,vecsizexh,vecsizeyh,vecsizezh,coulh,numch,LIST_R,LIST_X,LIST_Y,LIST_Z,NB_SPH,LIST_V,VOLHALO,coef1,NMAXHALO);

init_sphT(NB_SPH,Tinit,Timp,LIST_TEMP,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NCONT,NBCO,CONT,LIST_X,LIST_Y,LIST_Z,DCONTO,DeltaT,1,1);

init_sphC(NB_SPH,Cinit,Cimp,LIST_CONW,Cp_eff,LIST_M,LIST_V,p, Cps, Cpw, rhos, Mw,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NCONT,NBCO,CONT,LIST_X,LIST_Y,LIST_Z,DCONTO,DeltaC,1,1);

dt=crit_time_couplage(coord,NBCO, lambda, Cps, dens, DCONTO, Cst,Dw);

ExpVISUC(ite,dt,NB_SPH,H_TOT,V_TOT,Z_TOT,LIST_X,LIST_Y,LIST_Z,LIST_R,LIST_CONW,maxc,minc,EDGE,mingx,maxgx,mingy,maxgy,mingz,maxgz,LIST_GCX,LIST_GCY,LIST_GCZ,LIST_P,1);

ExpVISUT(ite,dt,NB_SPH,H_TOT,V_TOT,Z_TOT,LIST_X,LIST_Y,LIST_Z,LIST_R,LIST_TEMP,maxt,mint,EDGE,minfx,maxfx,minfy,maxfy,minfz,maxfz,LIST_FLX,LIST_FLY,LIST_FLZ,LIST_P,1);

Exp_temp_hydriq(ite,NB_SPH,dt,p,coef1,rhos,Mw,LIST_TEMP,LIST_V,LIST_X,LIST_Y,LIST_Z,LIST_CONW,H_TOT,V_TOT,Z_TOT,1);

R dte;
inicoh_coef(LIST_IND,TYPCO,dte,epsi,NBCO,CONT,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_M,LIST_I,DCONT,DCONTO,NCONT,VALCOH,VALAMO,NOCONT,NBCONTCO,DCONTX,DCONTY,DCONTZ,DCONTXO,DCONTYO,DCONTZO,Pi,Cse,Emu,rmu,nuu,amort);
dt=min(dt,dte);

 cout<<"Pas de temps de calcul           = "<<dt<<endl;   
 cout<<"Fraction volumique               = "<<PR_SPH<<endl;
 cout<<"Eps                              = "<<epsi<<endl;
 cout<<"Nombre de grilles                = "<<vecsizex*vecsizey*vecsizez<<endl;  	  
 cout<<"Nombre de contacts total         = "<<NBCO<<endl; 
 cout<<"Nombre de coordination           = "<<coord<<endl; 
 cout<<"________________________________________________"<<endl; 

 cout<<endl;
 cout<<endl;
 cout<<"             < DEBUT PROCESS >  "<<endl;
 cout<<"________________________________________________"<<endl; 
 cout<<endl;

R dpltx,dplty,dpltz;

while(ite<nmax)
     {

temp_hydriq(ite, coef2,NB_SPH,NBCO,dt,lambda,p,Cps,Cpw,coef1, CONT, NOCONT,NCONT,NBCONTCO,NBHALO,VOLHALO,NOHALO, LIST_M, rhos, Mw, Lw,  LIST_R, LIST_TEMP, EDGE, EDGE1, EDGE2, EDGE3, EDGE4,  EDGE5, EDGE6, mint, maxt, DCONTO, LIST_IND, N1, N2,  N3, N4, N5, N6, DeltaT,Dw,LIST_V, LIST_CONW,Cp_eff, LIST_GCX, LIST_GCY, LIST_GCZ,minc, maxc,DeltaC, A, B, Rh,Rg, k, h, Tov, Cov, H_TOT, V_TOT, Z_TOT); 
ther_load_hygrc(NBENREG,ite,Mw,rhow,LIST_CONW,Cinit,TYPCO,NBCO,CONT,DCONT,DCONTX,DCONTY,DCONTZ,DCONTO,DCONTXO,DCONTYO,DCONTZO);
forcecoh0(Epe,Epa,TYPCO,NBCO,NB_SPH,CONT,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_TX,LIST_TY,LIST_TZ,LIST_VX,LIST_VY,LIST_VZ,LIST_WX,LIST_WY,LIST_WZ,DCONTX,DCONTY,DCONTZ,NCONT,FX,FY,FZ,MTX,MTY,MTZ,VALCOH,VALAMO,FCJI,FOJI,MTJI,MTIJ);

      if(ite%NBENREG==0)
       {

        R t2p=omp_get_wtime () - t1p; 
 
        R  pcpuj=t2p/(24.*3600);
        int pnjcpu=int(pcpuj);
        R pcpuh=(pcpuj-pnjcpu)*24;
        int pnhcpu=int(pcpuh);
        R pcpum=(pcpuh-pnhcpu)*60.;
        int pnmcpu=int(pcpum);
        R pcpus=(pcpum-pnmcpu)*60.;
        int pnscpu=int(pcpus);                
 
        cout<<endl;
        cout<<"Pas:"<<ite<<"         PART:        "<<pnjcpu<<"j"<<" "<<pnhcpu<<"h"<<" "<<pnmcpu<<"mn"<<" "<<pnscpu<<"s"<<endl;

//grad_C(coef2,NB_SPH, NBCO, dt, Dw, coef1, CONT, NCONT, LIST_M, LIST_R, LIST_CONW, EDGE, mingx, maxgx, mingy, maxgy, mingz, maxgz, DCONTO, H_TOT,  V_TOT, Z_TOT, LIST_GCX, LIST_GCY, LIST_GCZ, NBHALO, NOHALO, NBCONTCO, NOCONT, VOLHALO, LIST_V); 

//R lambda_eff=lambda*(1-p);
//flux(coef2,NB_SPH,NBCO,dt,lambda_eff,coef1,CONT,NCONT,LIST_M,LIST_R,LIST_TEMP,EDGE,minfx,maxfx,minfy,maxfy,minfz,maxfz,DCONTO,H_TOT, V_TOT,Z_TOT,LIST_FLX,LIST_FLY,LIST_FLZ,NBHALO,NOHALO,NBCONTCO,NOCONT,VOLHALO,LIST_V);

        cout<<"min/max c:"<<minc<<"  "<<maxc<<endl;
        cout<<"min/max t:"<<mint<<"  "<<maxt<<endl;
        cout<<"min/max gz:"<<mingz<<"  "<<maxgz<<endl;
        cout<<"min/max fz:"<<minfz<<"  "<<maxfz<<endl;

	cout<<"________________________________________________"<<endl; 
	cout<<endl;

   // ExpVISUC(ite,dt,NB_SPH,H_TOT,V_TOT,Z_TOT,LIST_X,LIST_Y,LIST_Z,LIST_R,LIST_CONW,maxc,minc,EDGE,mingx,maxgx,mingy,maxgy,mingz,maxgz,LIST_GCX,LIST_GCY,LIST_GCZ,LIST_P,0);
    
   // ExpVISUT(ite,dt,NB_SPH,H_TOT,V_TOT,Z_TOT,LIST_X,LIST_Y,LIST_Z,LIST_R,LIST_TEMP,maxt,mint,EDGE,minfx,maxfx,minfy,maxfy,minfz,maxfz,LIST_FLX,LIST_FLY,LIST_FLZ,LIST_P,0);

	Exp_temp_hydriq(ite,NB_SPH,dt,p,coef1,rhos,Mw,LIST_TEMP,LIST_V,LIST_X,LIST_Y,LIST_Z,LIST_CONW,H_TOT,V_TOT,Z_TOT,0);
 
  }

verlet_qs_hygr(Ec,Epp,ite,NB_SPH,dt,LIST_M,LIST_I,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_XO,LIST_YO,LIST_ZO,LIST_TX,LIST_TY,LIST_TZ,LIST_VX,LIST_VY,LIST_VZ, LIST_WX,LIST_WY,LIST_WZ,LIST_AX,LIST_AY,LIST_AZ,LIST_AWX,LIST_AWY,LIST_AWZ,FX,FY,FZ,MTX,MTY,MTZ,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,NBENREG,dpltx,dplty,dpltz);
  
     ite++;
     ti+=dt;
    }

t2=CPUtime()-t1;
R t2p=omp_get_wtime () - t1p; 
cout<<"Iterations total number :"<<ite<<endl;


// Libération de la mémoire	
//free_memory();
cout<<"Process CPU time :"<<t2<<endl;
cout<<"Total CPU time :"<<CPUtime()-t1<<endl;
cout<<"Total Par. time :"<<omp_get_wtime () - t1p<<endl;

}

