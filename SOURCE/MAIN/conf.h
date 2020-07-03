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

using namespace std;

#ifndef __CONF__
#define __CONF__

////////////////////////////////////////////////////////////////////////////////
// Paramètres du système
int NB_SPH,NB_PAR,NBCO,NBCOP;
R N1,N2,N3,N4,N5,N6;
R R_SPH,H_TOT,V_TOT,Z_TOT,RMAX;
R H_POS,V_POS,Z_POS;
R dt,PR_SPH,epsi;
R Epe,Epa,Ec,Epp;
R RHALO;

R minvm,maxvm,sd;
R mintrac,maxtrac;
R minsig11,maxsig11;
R minsig12,maxsig12;
R minsig13,maxsig13;
R minsig22,maxsig22;
R minsig23,maxsig23;
R minsig33,maxsig33;
R minsig1,maxsig1;
R minsig2,maxsig2;
R minsig3,maxsig3;
R mindep1,maxdep1;
R mindep2,maxdep2;
R mindep3,maxdep3;
R mindef11,maxdef11;
R mindef22,maxdef22;
R mindef33,maxdef33;
R mindefe11,maxdefe11;
R mindefe22,maxdefe22;
R mindefe33,maxdefe33;
R mint,maxt,minc,maxc,minfx,maxfx,minfy,maxfy,minfz,maxfz,mingx,maxgx,mingy,maxgy,mingz,maxgz;

////////////////////////////////////////////////////////////////////////////////

// Grille
int ** coul;
int ** numc;

int * nocoul; 
int * nonumc;

int vecsizex;
int vecsizey;
int vecsizez;
int vecsize;

// Grille Halo
vector< vector<int> > coulh;
vector< vector<int> > numch;
int vecsizexh;
int vecsizeyh;
int vecsizezh;
int vecsizeh;

////////////////////////////////////////////////////////////////////////////////

// VAL MAX

long int NMAXMOYZ = 8;
long int NMAXZ    = 45;
long int NMAXHALO = 2000;
long int NMAXSPH;
long int NMAXCONT;
long int NMAXCONTP;
long int NMAXPAROI;

// OBJET SPHERE
R    * LIST_R;

R    * LIST_X;
R    * LIST_Y;
R    * LIST_Z;

R    * LIST_XA;
R    * LIST_YA;
R    * LIST_ZA;

R    * LIST_XO;
R    * LIST_YO;
R    * LIST_ZO;

R    * LIST_TX;
R    * LIST_TY;
R    * LIST_TZ;

R    * LIST_TXA;
R    * LIST_TYA;
R    * LIST_TZA;

R    * LIST_VX;
R    * LIST_VY;
R    * LIST_VZ;

R    * LIST_VXA;
R    * LIST_VYA;
R    * LIST_VZA;

R    * LIST_WX;
R    * LIST_WY;
R    * LIST_WZ;

R    * LIST_AX;
R    * LIST_AY;
R    * LIST_AZ;

R    * LIST_AWX;
R    * LIST_AWY;
R    * LIST_AWZ;

R    * LIST_TEMP;
R    * LIST_CONW;
R    * Cp_eff;
R    * DeltaT;
R    * DeltaC;

R    * LIST_FLX;
R    * LIST_FLY;
R    * LIST_FLZ;

R    * LIST_GCX;
R    * LIST_GCY;
R    * LIST_GCZ;

R    * LIST_M;
R    * LIST_V;
R    * LIST_I;
R    * LIST_IND;

R    * VOLHALO;	

int  * LIST_C;
int  * LIST_H;
bool * LIST_P;
int  * LIST_B;

R    * FX;
R    * FY;
R    * FZ;

R    * MTX;
R    * MTY;
R    * MTZ;

R    * FIX;
R    * FIY;
R    * FIZ;

R    * MTIX;
R    * MTIY;
R    * MTIZ;

bool * EDGE;
bool * EDGE1;
bool * EDGE2;
bool * EDGE3;
bool * EDGE4;
bool * EDGE5;
bool * EDGE6;
bool * EDGEC1;
bool * EDGEC2;
bool * EDGEC;
bool * EDGE1M1;
bool * EDGE1M2;
bool * EDGE13;
bool * EDGE43;
bool * EDGE6X;
bool * EDGE7X;

unsigned int ** NOCONT;
int * NBCONTCO;
unsigned int ** NOHALO;
int * NBHALO;

R * VONMIS;
R * TRACE;
R * SIG11;
R * SIG12;
R * SIG13;
R * SIG22;
R * SIG23;
R * SIG33;
R * SIG1;
R * SIG2;
R * SIG3;
R * DEP1;
R * DEP2;
R * DEP3;
R * EPSI11;
R * EPSI22;
R * EPSI33;
R * EPSE11;
R * EPSE22;
R * EPSE33;

// OBJET CONTACT
int ** CONT;
bool * TYPCO;

R * vs;
R * DCONTX;
R * DCONTY;
R * DCONTZ;
R * DCONTXO;
R * DCONTYO;
R * DCONTZO;
R * DCONT;
R * DCONTO;
R ** NCONT;
R ** FCJI;
R ** FOJI;
R ** MTJI;
R ** MTIJ;
R ** VALCOH;
R ** VALAMO;

// OBJET PAROI

R ** LIST_PX;
R ** LIST_PY;
R ** LIST_PZ;

R ** LIST_PVX;
R ** LIST_PVY;
R ** LIST_PVZ;

R ** LIST_PAX;
R ** LIST_PAY;
R ** LIST_PAZ;

R   * LIST_PM;
R  ** LIST_PN;

// OBJET CONTACT PAROI

int ** CONTP;
R * DCONTP;
R ** NCONTP;

////////////////////////////////////////////////////////////////////////////////

void allocat_grid(R H_TOT, R V_TOT, R Z_TOT, R RMAX, int & vecsizex, int & vecsizey, int & vecsizez, int & vecsize){

	vecsizex=(floor(H_TOT/(2*RMAX))>4)?(floor(H_TOT/(2*RMAX))-1):4;
	vecsizex-=2;
	vecsizey=(floor(V_TOT/(2*RMAX))>4)?(floor(V_TOT/(2*RMAX))-1):4;
	vecsizey-=2;
	vecsizez=(floor(Z_TOT/(2*RMAX))>4)?(floor(Z_TOT/(2*RMAX))-1):4;		
	vecsizez-=2;
	
if(vecsizex<3) vecsizex=3;
if(vecsizey<3) vecsizey=3;
if(vecsizez<3) vecsizez=3;

vecsize=vecsizex*vecsizey*vecsizez;

// grid

coul  = new int*[vecsize];
coul[0] = new int[15*vecsize];
for(int i=1;i<vecsize;i++){
coul[i]= coul[i-1] + 15;
}	

numc  = new int*[vecsize];
numc[0] = new int[14*vecsize];
for(int i=1;i<vecsize;i++){
numc[i]= numc[i-1] + 14;
}	

nocoul = new int[vecsize];
nonumc = new int[vecsize];

for (int ii = 0; ii < vecsize; ii++) {

	int NX  = ii%vecsizex;
	int NY  = (ii%(vecsizex*vecsizey))/vecsizex;
	int NZ  = ii/(vecsizex*vecsizey);
	
	int NXP1= (NX+vecsizex+1)%vecsizex;
	int NXM1= (NX+vecsizex-1)%vecsizex;
	int NYP1= (NY+vecsizey+1)%vecsizey;	
	int NYM1= (NY+vecsizey-1)%vecsizey;
	//int NZP1= (NZ+vecsizez+1)%vecsizez;	
	int NZM1= (NZ+vecsizez-1)%vecsizez;		

	numc[ii][0]=ii;
	nonumc[ii]=1;

	if(NZM1==(NZ-1)){
	 if((NYM1==(NY-1))&&(NXM1==(NX-1))){numc[ii][nonumc[ii]]=NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZM1;nonumc[ii]++;}

	 if(NYM1==(NY-1)){numc[ii][nonumc[ii]]=NX+vecsizex*NYM1+vecsizex*vecsizey*NZM1;nonumc[ii]++;}

	 if((NYM1==(NY-1))&&(NXP1==(NX+1))){numc[ii][nonumc[ii]]=NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZM1;nonumc[ii]++;}

	 if(NXM1==(NX-1)){numc[ii][nonumc[ii]]=NXM1+vecsizex*NY+vecsizex*vecsizey*NZM1;nonumc[ii]++;}	 
	
	 numc[ii][nonumc[ii]]=NX+vecsizex*NY+vecsizex*vecsizey*NZM1;nonumc[ii]++;

	 if(NXP1==(NX+1)){numc[ii][nonumc[ii]]=NXP1+vecsizex*NY+vecsizex*vecsizey*NZM1;nonumc[ii]++;}

	 if((NYP1==(NY+1))&&(NXM1==(NX-1))){numc[ii][nonumc[ii]]=NXM1+vecsizex*NYP1+vecsizex*vecsizey*NZM1;nonumc[ii]++;}

	 if(NYP1==(NY+1)){numc[ii][nonumc[ii]]=NX+vecsizex*NYP1+vecsizex*vecsizey*NZM1;nonumc[ii]++;}
	
	 if((NYP1==(NY+1))&&(NXP1==(NX+1))){numc[ii][nonumc[ii]]=NXP1+vecsizex*NYP1+vecsizex*vecsizey*NZM1;nonumc[ii]++;}	 	 
	}	
		
	if(NYM1==(NY-1)){
	
	 if(NXM1==(NX-1)){numc[ii][nonumc[ii]]=NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZ;nonumc[ii]++;}
	
	 numc[ii][nonumc[ii]]=NX+vecsizex*NYM1+vecsizex*vecsizey*NZ;nonumc[ii]++;
	}	

	if(NXP1==(NX+1)){
	
 	if(NYM1==(NY-1)){numc[ii][nonumc[ii]]=NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZ;nonumc[ii]++;}
	
 	numc[ii][nonumc[ii]]=NXP1+vecsizex*NY+vecsizex*vecsizey*NZ; nonumc[ii]++;

	}	
	
 
 }	 
 
 
 cout<<"NXYZ:"<<vecsizex<<", "<<vecsizey<<", "<<vecsizez<<endl;

}

////////////////////////////////////////////////////////////////////////////////

void allocat_memory(){

    const rlim_t kStackSize = 512 * 2048 * 2048;   // min stack size = 32 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
      //    cout << "StackLimit - expected : " << rl.rlim_cur << " - " << kStackSize << endl;

        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
	  cout << "Upgraded StackLimit soft - max : " << rl.rlim_cur << " - " << rl.rlim_max << endl;
        }
    }

// Lecture dimensions

ifstream f("APPLI/carttest/test",ios::in);
f >> NMAXSPH;
f.close();

ifstream g("APPLI/paroi/testp",ios::in);
g >> NMAXPAROI;
g.close();

NMAXCONT=NMAXSPH*NMAXMOYZ;
NMAXCONTP=NMAXSPH;
if(NMAXPAROI==0){
NMAXCONTP=0;
}


// OBJET SPHERE

LIST_R = new R[NMAXSPH];

LIST_X = new R[NMAXSPH];
LIST_Y = new R[NMAXSPH];
LIST_Z = new R[NMAXSPH];

LIST_XA = new R[NMAXSPH];
LIST_YA = new R[NMAXSPH];
LIST_ZA = new R[NMAXSPH];

LIST_XO = new R[NMAXSPH];
LIST_YO = new R[NMAXSPH];
LIST_ZO = new R[NMAXSPH];

LIST_TX = new R[NMAXSPH];
LIST_TY = new R[NMAXSPH];
LIST_TZ = new R[NMAXSPH];

LIST_TXA = new R[NMAXSPH];
LIST_TYA = new R[NMAXSPH];
LIST_TZA = new R[NMAXSPH];

LIST_VX = new R[NMAXSPH];
LIST_VY = new R[NMAXSPH];
LIST_VZ = new R[NMAXSPH];

LIST_VXA = new R[NMAXSPH];
LIST_VYA = new R[NMAXSPH];
LIST_VZA = new R[NMAXSPH];

LIST_WX  = new R[NMAXSPH];
LIST_WY  = new R[NMAXSPH];
LIST_WZ  = new R[NMAXSPH];

LIST_AX = new R[NMAXSPH];
LIST_AY = new R[NMAXSPH];
LIST_AZ = new R[NMAXSPH];

LIST_AWX = new R[NMAXSPH];
LIST_AWY = new R[NMAXSPH];
LIST_AWZ = new R[NMAXSPH];

LIST_TEMP = new R[NMAXSPH];
LIST_CONW = new R[NMAXSPH];
Cp_eff= new R[NMAXSPH];

DeltaT = new R[NMAXSPH];
DeltaC = new R[NMAXSPH];
LIST_FLX = new R[NMAXSPH];
LIST_FLY = new R[NMAXSPH];
LIST_FLZ = new R[NMAXSPH];
LIST_GCX = new R[NMAXSPH];
LIST_GCY = new R[NMAXSPH];
LIST_GCZ = new R[NMAXSPH];

LIST_M = new R[NMAXSPH];
LIST_V = new R[NMAXSPH];
LIST_IND = new R[NMAXSPH];
LIST_I = new R[NMAXSPH];
LIST_C = new int[NMAXSPH];
LIST_P = new bool[NMAXSPH];
LIST_B = new int[NMAXSPH];
LIST_H = new int[NMAXSPH];

VOLHALO = new R[NMAXSPH];

FX = new R[NMAXSPH];
FY = new R[NMAXSPH];
FZ = new R[NMAXSPH];

MTX = new R[NMAXSPH];
MTY = new R[NMAXSPH];
MTZ = new R[NMAXSPH];

FIX = new R[NMAXSPH];
FIY = new R[NMAXSPH];
FIZ = new R[NMAXSPH];

MTIX = new R[NMAXSPH];
MTIY = new R[NMAXSPH];
MTIZ = new R[NMAXSPH];

EDGE     = new bool[NMAXSPH];
EDGE1    = new bool[NMAXSPH];
EDGE2    = new bool[NMAXSPH];
EDGE3    = new bool[NMAXSPH];
EDGE4    = new bool[NMAXSPH];
EDGE5    = new bool[NMAXSPH];
EDGE6    = new bool[NMAXSPH];
EDGEC1    = new bool[NMAXSPH];
EDGEC2    = new bool[NMAXSPH];
EDGEC    = new bool[NMAXSPH];

EDGE1M1   = new bool[NMAXSPH];
EDGE1M2   = new bool[NMAXSPH];
EDGE13    = new bool[NMAXSPH];
EDGE43    = new bool[NMAXSPH];
EDGE6X    = new bool[NMAXSPH];
EDGE7X    = new bool[NMAXSPH];

NBCONTCO = new int[NMAXSPH];

NOCONT   = new unsigned int*[NMAXSPH];
NOCONT[0] = new unsigned int[NMAXSPH*NMAXZ];
for(int i=1;i<NMAXSPH;i++){
NOCONT[i]= NOCONT[i-1] + NMAXZ;
}


NBHALO   = new int[NMAXSPH];

NOHALO   = new unsigned int*[NMAXSPH];
NOHALO[0] = new unsigned int[NMAXSPH*NMAXHALO];
for(int i=1;i<NMAXSPH;i++){
NOHALO[i]= NOHALO[i-1] + NMAXHALO;
}


VONMIS = new R[NMAXSPH];
TRACE  = new R[NMAXSPH];
SIG11 = new R[NMAXSPH];
SIG12 = new R[NMAXSPH];
SIG13 = new R[NMAXSPH];
SIG22 = new R[NMAXSPH];
SIG23 = new R[NMAXSPH];
SIG33 = new R[NMAXSPH];
SIG1 = new R[NMAXSPH];
SIG2 = new R[NMAXSPH];
SIG3 = new R[NMAXSPH];
DEP1 = new R[NMAXSPH];
DEP2 = new R[NMAXSPH];
DEP3 = new R[NMAXSPH];

EPSI11  = new R[NMAXSPH];
EPSI22  = new R[NMAXSPH];
EPSI33  = new R[NMAXSPH];
EPSE11  = new R[NMAXSPH];
EPSE22  = new R[NMAXSPH];
EPSE33  = new R[NMAXSPH];

// OBJET CONTACT

CONT    = new int*[NMAXCONT];
CONT[0] = new int[NMAXCONT*2];
for(int i=1;i<NMAXCONT;i++){
CONT[i]= CONT[i-1] + 2;
}

TYPCO   = new bool[NMAXCONT];
DCONTX  = new R[NMAXCONT];
DCONTY  = new R[NMAXCONT];
DCONTZ  = new R[NMAXCONT];
DCONTXO = new R[NMAXCONT];
DCONTYO = new R[NMAXCONT];
DCONTZO = new R[NMAXCONT];
DCONT   = new R[NMAXCONT];
DCONTO  = new R[NMAXCONT];
vs	= new R[NMAXCONT];

NCONT   = new R *[NMAXCONT];
NCONT[0]= new R [NMAXCONT*9];
for(int i=1;i<NMAXCONT;i++){
NCONT[i]= NCONT[i-1] + 9;
}

FCJI    = new R *[NMAXCONT];
FCJI[0] = new R [NMAXCONT*3];
for(int i=1;i<NMAXCONT;i++){
FCJI[i]  = FCJI[i-1] + 3;
}

FOJI    = new R *[NMAXCONT];
FOJI[0] = new R [NMAXCONT*3];
for(int i=1;i<NMAXCONT;i++){
FOJI[i]  = FOJI[i-1] + 3;
}

MTJI    = new R *[NMAXCONT];
MTJI[0] = new R [NMAXCONT*3];
for(int i=1;i<NMAXCONT;i++){
MTJI[i]  = MTJI[i-1] + 3;
}

MTIJ    = new R *[NMAXCONT];
MTIJ[0] = new R [NMAXCONT*3];
for(int i=1;i<NMAXCONT;i++){
MTIJ[i]  = MTIJ[i-1] + 3;
}

VALCOH    = new R *[NMAXCONT];
VALCOH[0] = new R [NMAXCONT*6];
for(int i=1;i<NMAXCONT;i++){
VALCOH[i]  = VALCOH[i-1] + 6;
}

VALAMO    = new R *[NMAXCONT];
VALAMO[0] = new R [NMAXCONT*6];
for(int i=1;i<NMAXCONT;i++){
VALAMO[i]  = VALAMO[i-1] + 6;
}

// OBJET PAROI

LIST_PX    = new R *[NMAXPAROI];
LIST_PX[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PX[i]  = LIST_PX[i-1] + 4;	
}

LIST_PY    = new R *[NMAXPAROI];
LIST_PY[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PY[i]  = LIST_PY[i-1] + 4;	
}

LIST_PZ    = new R *[NMAXPAROI];
LIST_PZ[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PZ[i]  = LIST_PZ[i-1] + 4;	
}

LIST_PVX    = new R *[NMAXPAROI];
LIST_PVX[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PVX[i]  = LIST_PVX[i-1] + 4;	
}

LIST_PVY    = new R *[NMAXPAROI];
LIST_PVY[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PVY[i]  = LIST_PVY[i-1] + 4;	
}

LIST_PVZ    = new R *[NMAXPAROI];
LIST_PVZ[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PVZ[i]  = LIST_PVZ[i-1] + 4;	
}

LIST_PAX    = new R *[NMAXPAROI];
LIST_PAX[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PAX[i]  = LIST_PAX[i-1] + 4;	
}

LIST_PAY    = new R *[NMAXPAROI];
LIST_PAY[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PAY[i]  = LIST_PAY[i-1] + 4;	
}

LIST_PAZ    = new R *[NMAXPAROI];
LIST_PAZ[0] = new R [NMAXPAROI*4];
for(int i=1;i<NMAXPAROI;i++){
LIST_PAZ[i]  = LIST_PAZ[i-1] + 4;	
}

LIST_PM = new R[NMAXPAROI];

LIST_PN = new R *[NMAXCONT];
LIST_PN[0] = new R [NMAXCONT*9];
for(int i=1;i<NMAXCONT;i++){
LIST_PN[i]  = LIST_PN[i-1] + 9;	
}

// OBJET CONTACT PAROI

CONTP = new int*[NMAXCONTP];
CONTP[0] = new int[NMAXCONTP*2];
for(int i=1;i<NMAXCONTP;i++){
CONTP[i]  =CONTP[i-1] + 2;	
}

DCONTP  = new R[NMAXCONTP];

NCONTP  = new R *[NMAXCONTP];
NCONTP[0]  = new R [NMAXCONTP*9];
for(int i=1;i<NMAXCONTP;i++){
NCONTP[i] = NCONTP[i-1] + 9;	
}

}

////////////////////////////////////////////////////////////////////////////////

void free_memory(){

// OBJET CORCI 

delete [] LIST_R;

delete [] LIST_X;
delete [] LIST_Y;
delete [] LIST_Z;

delete [] LIST_XA;
delete [] LIST_YA;
delete [] LIST_ZA;

delete [] LIST_XO;
delete [] LIST_YO;
delete [] LIST_ZO;

delete [] LIST_TXA;
delete [] LIST_TYA;
delete [] LIST_TZA;

delete [] LIST_TX;
delete [] LIST_TY;
delete [] LIST_TZ;

delete [] LIST_VX;
delete [] LIST_VY;
delete [] LIST_VZ;

delete [] LIST_WX;
delete [] LIST_WY;
delete [] LIST_WZ;

delete [] LIST_AX;
delete [] LIST_AY;
delete [] LIST_AZ;

delete [] LIST_AWX;
delete [] LIST_AWY;
delete [] LIST_AWZ;

delete [] LIST_TEMP;
delete [] LIST_CONW;
delete [] Cp_eff;
delete [] DeltaT;
delete [] DeltaC;
delete [] LIST_FLX;
delete [] LIST_FLY;
delete [] LIST_FLZ;
delete [] LIST_GCX;
delete [] LIST_GCY;
delete [] LIST_GCZ;

delete [] LIST_M;
delete [] LIST_V;
delete [] LIST_IND;
delete [] LIST_I;
delete [] LIST_C;
delete [] LIST_P;
delete [] LIST_B;
delete [] LIST_H;

delete [] VOLHALO;

delete [] FX;
delete [] FY;
delete [] FZ;

delete [] MTX;
delete [] MTY;
delete [] MTZ;

delete [] EDGE;
delete [] EDGE1;
delete [] EDGE2;
delete [] EDGE3;
delete [] EDGE4;
delete [] EDGE5;
delete [] EDGE6;
delete [] EDGEC1;
delete [] EDGEC2;

delete [] EDGE1M1;
delete [] EDGE1M2;
delete [] EDGE13;
delete [] EDGE43;
delete [] EDGE6X;
delete [] EDGE7X;

for(int i=0;i<NMAXSPH;i++){
delete [] NOCONT[i];
}
delete [] NBCONTCO;

for(int i=0;i<NMAXSPH;i++){
delete [] NOHALO[i];
}
delete [] NBHALO;

delete [] VONMIS;
delete [] TRACE;
delete [] SIG11;
delete [] SIG12;
delete [] SIG13;
delete [] SIG22;
delete [] SIG23;
delete [] SIG33;
delete [] SIG1;
delete [] SIG2;
delete [] SIG3;
delete [] DEP1;
delete [] DEP2;
delete [] DEP3;
delete [] EPSI11;
delete [] EPSI22;
delete [] EPSI33;
delete [] EPSE11;
delete [] EPSE22;
delete [] EPSE33;

// OBJET CONTACT

for(int i=0;i<NMAXCONT;i++){
delete [] CONT[i];
}

delete [] TYPCO;
delete [] DCONTX;
delete [] DCONTY;
delete [] DCONTZ;
delete [] DCONTXO;
delete [] DCONTYO;
delete [] DCONTZO;
delete [] DCONT;
delete [] DCONTO;
delete [] vs;

for(int i=0;i<NMAXCONT;i++){
delete [] NCONT[i];
}

for(int i=0;i<NMAXCONT;i++){
delete [] FCJI[i];
}

for(int i=0;i<NMAXCONT;i++){
delete [] FOJI[i];
}

for(int i=0;i<NMAXCONT;i++){
delete [] MTJI[i];
}

for(int i=0;i<NMAXCONT;i++){
delete [] MTIJ[i];
}

for(int i=0;i<NMAXCONT;i++){
delete [] VALCOH[i];
}

for(int i=0;i<NMAXCONT;i++){
delete [] VALAMO[i];
}

// OBJET PAROI

for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PX[i];	
}
for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PY[i];	
}
for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PZ[i];	
}

for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PVX[i];		
}
for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PVY[i];	
}
for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PVZ[i];		
}

for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PAX[i];		
}
for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PAY[i];	
}
for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PAZ[i];		
}


delete [] LIST_PM;

for(int i=0;i<NMAXPAROI;i++){
delete [] LIST_PN[i];	
}


// OBJET CONTACT PAROI

for(int i=0;i<NMAXCONTP;i++){
delete [] CONTP[i];	
}

delete [] DCONTP;

for(int i=0;i<NMAXCONTP;i++){
delete [] NCONTP[i];	
}

// OBJET grid

for(int i=0;i<vecsize;i++){
delete [] coul[i];	
}

delete [] coul;
delete [] nocoul;

for(int i=0;i<vecsize;i++){
delete [] numc[i];	
}

delete [] numc;
delete [] nonumc;


}



#endif

