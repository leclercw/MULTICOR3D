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
#include "omp.h"

using namespace std;

typedef double R;

const R Pi=3.14159265;
const R g=9.81;
const R tol=1e-6;
const R dens=7800.;
const R cov=0.3;
const R eps=0.001;
const R dt=4e-4;
const R fric=0.;

const int NBCORCI=6000000;
const int NBCONT=32000000;
const int NB_THREADS=4;

R    * LIST_R = new R[NBCORCI];

R    * LIST_X = new R[NBCORCI];
R    * LIST_Y = new R[NBCORCI];
R    * LIST_Z = new R[NBCORCI];

R    * LIST_VX = new R[NBCORCI];
R    * LIST_VY = new R[NBCORCI];
R    * LIST_VZ = new R[NBCORCI];

R    * LIST_AX = new R[NBCORCI];
R    * LIST_AY = new R[NBCORCI];
R    * LIST_AZ = new R[NBCORCI];

R    * LIST_A = new R[NBCORCI];
R    * LIST_M = new R[NBCORCI];
int  * LIST_C = new int[NBCORCI];

R    * FX = new R[NBCORCI];
R    * FY = new R[NBCORCI];
R    * FZ = new R[NBCORCI];

int ** CONT = new int*[NBCONT];
R * DCONT  = new R[NBCONT];
R ** NCONT  = new R *[NBCONT];

vector< vector<int> > coul;
vector< vector<int> > numc;
vector< vector<int> > numc2;
int vecsizex;
int vecsizey;
int vecsizez;
int vecsize;

inline R CPUtime(){
#ifdef SYSTIMES
  struct tms buf;
  if (times(&buf)!=-1)
    return ((R)buf.tms_utime+(R)buf.tms_stime)/(long) sysconf(_SC_CLK_TCK);
  else
#endif
    return ((R) clock())/CLOCKS_PER_SEC;
}

template<typename T>
std::string to_string( const T & Value )
{
    // utiliser un flux de sortie pour créer la chaîne
    ostringstream oss;
    // écrire la valeur dans le flux
    oss << Value;
    // renvoyer une string
    return oss.str();
}

template <class T>
void from_string(T& t,
                 const string& s,
                 ios_base& (*f)(ios_base&))
{
  istringstream iss(s);
  iss >> f >> t;
}


inline void maxi(R &Reff,R * LIST_R,int nD)
{
Reff=0;

  for(int i=0;i<nD;i++){
  Reff=(Reff>LIST_R[i])?Reff:LIST_R[i];
  }

}

inline void init_size(int & vecsize, int & vecsizex,int & vecsizey, int & vecsizez, vector< vector<int>  > & coul,vector< vector<int>  > & numc,vector< vector<int>  > & numc2, R H_TOT, R V_TOT, R Z_TOT, int NB_SPH, R R_SPH)
{

// Vecteur des couleurs
/*
R Reff = ((2*R_SPH+1e-8)>(2*H_TOT/pow(NB_SPH,1./3)))?(2*R_SPH+1e-8):(2*H_TOT/pow(NB_SPH,1./3));
vecsizex=(floor(H_TOT/Reff)>4)?floor(H_TOT/Reff):4;
vecsizex+=2;
Reff = ((2*R_SPH+1e-8)>(2*V_TOT/pow(NB_SPH,1./3)))?(2*R_SPH+1e-8):(2*V_TOT/pow(NB_SPH,1./3));
vecsizey=(floor(V_TOT/Reff)>4)?floor(V_TOT/Reff):4;
vecsizey+=2;
Reff = ((2*R_SPH+1e-8)>(2*Z_TOT/pow(NB_SPH,1./3)))?(2*R_SPH+1e-8):(2*Z_TOT/pow(NB_SPH,1./3));
vecsizez=(floor(Z_TOT/Reff)>4)?floor(Z_TOT/Reff):4;
vecsizez+=2;*/

	R Reff = 2*R_SPH;
	vecsizex=(floor(H_TOT/Reff)>4)?(floor(H_TOT/Reff)-1):4;
	vecsizex-=2;
	Reff =   2*R_SPH;
	vecsizey=(floor(V_TOT/Reff)>4)?(floor(V_TOT/Reff)-1):4;
	vecsizey-=2;
	Reff =   2*R_SPH;
	vecsizez=(floor(Z_TOT/Reff)>4)?(floor(Z_TOT/Reff)-1):4;
	vecsizez-=2;

vecsize=vecsizex*vecsizey*vecsizez;

for (int i = 0; i < vecsize; i++) {
    coul.push_back(vector<int>()); // Add an empty row
}
for (int ii = 0; ii < vecsize; ii++) {
	numc.push_back(vector<int>()); // Add an empty row
}
for (int ii = 0; ii < vecsize; ii++) {
	numc2.push_back(vector<int>()); // Add an empty row
}

for (int ii = 0; ii < vecsize; ii++) {
 //cout<<"ii :"<<ii<<" - "<<vecsizex<<endl;
	int NX  = ii%vecsizex;
	int NY  = (ii%(vecsizex*vecsizey))/vecsizex;
	int NZ  = ii/(vecsizex*vecsizey);

	int NXP1= (NX+vecsizex+1)%vecsizex;
	int NXM1= (NX+vecsizex-1)%vecsizex;
	int NYP1= (NY+vecsizey+1)%vecsizey;
	int NYM1= (NY+vecsizey-1)%vecsizey;
	int NZP1= (NZ+vecsizez+1)%vecsizez;
	int NZM1= (NZ+vecsizez-1)%vecsizez;

	numc[ii].push_back(ii);
	if(NZM1==(NZ-1)){
	 if((NYM1==(NY-1))&&(NXM1==(NX-1))){numc[ii].push_back(NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
	 if(NYM1==(NY-1)){numc[ii].push_back(NX+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
	 if((NYM1==(NY-1))&&(NXP1==(NX+1))){numc[ii].push_back(NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
	 if(NXM1==(NX-1)){numc[ii].push_back(NXM1+vecsizex*NY+vecsizex*vecsizey*NZM1);}
	 numc[ii].push_back(NX+vecsizex*NY+vecsizex*vecsizey*NZM1);
	 if(NXP1==(NX+1)){numc[ii].push_back(NXP1+vecsizex*NY+vecsizex*vecsizey*NZM1);}
	 if((NYP1==(NY+1))&&(NXM1==(NX-1))){numc[ii].push_back(NXM1+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
	 if(NYP1==(NY+1)){numc[ii].push_back(NX+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
	 if((NYP1==(NY+1))&&(NXP1==(NX+1))){numc[ii].push_back(NXP1+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
	}

	if(NYM1==(NY-1)){
	 if(NXM1==(NX-1)){numc[ii].push_back((NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZ));}
	 numc[ii].push_back((NX+vecsizex*NYM1+vecsizex*vecsizey*NZ));
	}
	if(NXP1==(NX+1)){
	 if(NYM1==(NY-1)){numc[ii].push_back((NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZ));}
	 numc[ii].push_back((NXP1+vecsizex*NY+vecsizex*vecsizey*NZ));
	}

	numc2[ii].push_back(ii);

	if(NZM1==(NZ-1)){
	 if((NYM1==(NY-1))&&(NXM1==(NX-1))){numc2[ii].push_back(NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
	 if(NYM1==(NY-1)){numc2[ii].push_back(NX+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
	 if((NYM1==(NY-1))&&(NXP1==(NX+1))){numc2[ii].push_back(NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
	 if(NXM1==(NX-1)){numc2[ii].push_back(NXM1+vecsizex*NY+vecsizex*vecsizey*NZM1);}
	 numc2[ii].push_back(NX+vecsizex*NY+vecsizex*vecsizey*NZM1);
	 if(NXP1==(NX+1)){numc2[ii].push_back(NXP1+vecsizex*NY+vecsizex*vecsizey*NZM1);}
	 if((NYP1==(NY+1))&&(NXM1==(NX-1))){numc2[ii].push_back(NXM1+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
	 if(NYP1==(NY+1)){numc2[ii].push_back(NX+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
	 if((NYP1==(NY+1))&&(NXP1==(NX+1))){numc2[ii].push_back(NXP1+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
	}

	 if((NYM1==(NY-1))&&(NXM1==(NX-1))){numc2[ii].push_back(NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZ);}
	 if(NYM1==(NY-1)){numc2[ii].push_back(NX+vecsizex*NYM1+vecsizex*vecsizey*NZ);}
	 if((NYM1==(NY-1))&&(NXP1==(NX+1))){numc2[ii].push_back(NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZ);}
	 if(NXM1==(NX-1)){numc2[ii].push_back(NXM1+vecsizex*NY+vecsizex*vecsizey*NZ);}
	 if(NXP1==(NX+1)){numc2[ii].push_back(NXP1+vecsizex*NY+vecsizex*vecsizey*NZ);}
	 if((NYP1==(NY+1))&&(NXM1==(NX-1))){numc2[ii].push_back(NXM1+vecsizex*NYP1+vecsizex*vecsizey*NZ);}
	 if(NYP1==(NY+1)){numc2[ii].push_back(NX+vecsizex*NYP1+vecsizex*vecsizey*NZ);}
	 if((NYP1==(NY+1))&&(NXP1==(NX+1))){numc2[ii].push_back(NXP1+vecsizex*NYP1+vecsizex*vecsizey*NZ);}

	if(NZP1==(NZ+1)){
	 if((NYM1==(NY-1))&&(NXM1==(NX-1))){numc2[ii].push_back(NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZP1);}
	 if(NYM1==(NY-1)){numc2[ii].push_back(NX+vecsizex*NYM1+vecsizex*vecsizey*NZP1);}
	 if((NYM1==(NY-1))&&(NXP1==(NX+1))){numc2[ii].push_back(NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZP1);}
	 if(NXM1==(NX-1)){numc2[ii].push_back(NXM1+vecsizex*NY+vecsizex*vecsizey*NZP1);}
	 numc2[ii].push_back(NX+vecsizex*NY+vecsizex*vecsizey*NZP1);
	 if(NXP1==(NX+1)){numc2[ii].push_back(NXP1+vecsizex*NY+vecsizex*vecsizey*NZP1);}
	 if((NYP1==(NY+1))&&(NXM1==(NX-1))){numc2[ii].push_back(NXM1+vecsizex*NYP1+vecsizex*vecsizey*NZP1);}
	 if(NYP1==(NY+1)){numc2[ii].push_back(NX+vecsizex*NYP1+vecsizex*vecsizey*NZP1);}
	 if((NYP1==(NY+1))&&(NXP1==(NX+1))){numc2[ii].push_back(NXP1+vecsizex*NYP1+vecsizex*vecsizey*NZP1);}
	}



 }

cout<<"vecsize:"<<vecsizex<<", "<<vecsizey<<", "<<vecsizez<<endl;

}


inline void init_sph(R mu,int & nbd, R & pbd, int vecsizex, int vecsizey, int vecsizez, R PR_SPH, int H_TOT, int V_TOT, int Z_TOT, int NB_SPH, R R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, int * LIST_C, R * LIST_A, R * LIST_M, R * FX, R * FY, R * FZ)
{

R crit = pow(3.*H_TOT*V_TOT*Z_TOT*PR_SPH/(4.*Pi*NB_SPH),1./3);

cout<<"Nombre de spheres : "<<NB_SPH<<endl;

	// DISTRIBUTION INITIALE DE SPHERES (SUR LES BORDS ET LES DISQUES)
    nbd=0;
    pbd=0.;

    int nbbordx=int(R(H_TOT)/(2.*crit))+1;
    int nbbordy=int(R(V_TOT)/(2.*crit))+1;
    int nbbordz=int(R(Z_TOT)/(2.*crit))+1;

    R raybx,rayby,raybz,raybe;
    raybx=R(H_TOT)/(2.*nbbordx);
    rayby=R(V_TOT)/(2.*nbbordy);
    raybz=R(Z_TOT)/(2.*nbbordz);

    // Bords XY-/+
    raybe=(raybx<rayby)?raybx:rayby;

    for(int jt=0;jt<=nbbordy;jt++)
    {
		for(int it=0;it<=nbbordx;it++)
		{
				LIST_X[nbd] = 2.*it*raybx;
				LIST_Y[nbd] = 2.*jt*rayby;
				LIST_Z[nbd] = tol;
				LIST_R[nbd] = raybe	;
				pbd+=4.*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]/3.;

				int NX= ((int) LIST_X[nbd])*vecsizex/H_TOT;
				if(NX<0) NX=0;
				if(NX>=vecsizex) NX=vecsizex-1;

				int NY= ((int) LIST_Y[nbd])*vecsizey/V_TOT;
				if(NY<0) NY=0;
				if(NY>=vecsizey) NY=vecsizey-1;

				int NZ= ((int) LIST_Z[nbd])*vecsizez/Z_TOT;
				if(NZ<0) NZ=0;
				if(NZ>=vecsizez) NZ=vecsizez-1;

				coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(nbd);
				LIST_C[nbd] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;
				LIST_A[nbd] = 0.;

				LIST_VX[nbd] = 0.;
				LIST_VY[nbd] = 0.;
				LIST_VZ[nbd] = 0.;

				LIST_AX[nbd] = 0.;
				LIST_AY[nbd] = 0.;
				LIST_AZ[nbd] = 0.;

				FX[nbd] = 0.;
				FY[nbd] = 0.;
				FZ[nbd] = 0.;

				LIST_M[nbd]  = dens*(4./3)*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]*10.;

				nbd++;
		}
   }

    for(int jt=0;jt<=nbbordy;jt++)
    {
		for(int it=0;it<=nbbordx;it++)
		{
				LIST_X[nbd] = 2.*it*raybx;
				LIST_Y[nbd] = 2.*jt*rayby;
				LIST_Z[nbd] = R(Z_TOT)-tol;
				LIST_R[nbd] = raybe	;
				pbd+=4.*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]/3.;

				int NX= ((int) LIST_X[nbd])*vecsizex/H_TOT;
				if(NX<0) NX=0;
				if(NX>=vecsizex) NX=vecsizex-1;

				int NY= ((int) LIST_Y[nbd])*vecsizey/V_TOT;
				if(NY<0) NY=0;
				if(NY>=vecsizey) NY=vecsizey-1;

				int NZ= ((int) LIST_Z[nbd])*vecsizez/Z_TOT;
				if(NZ<0) NZ=0;
				if(NZ>=vecsizez) NZ=vecsizez-1;

				coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(nbd);
				LIST_C[nbd] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;
				LIST_A[nbd] = 0.;

				LIST_VX[nbd] = 0.;
				LIST_VY[nbd] = 0.;
				LIST_VZ[nbd] = 0.;

				LIST_AX[nbd] = 0.;
				LIST_AY[nbd] = 0.;
				LIST_AZ[nbd] = 0.;

				FX[nbd] = 0.;
				FY[nbd] = 0.;
				FZ[nbd] = 0.;

				LIST_M[nbd]  = dens*(4./3)*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]*10.;

				nbd++;
		}
   }

    // Bords YZ-/+
    raybe=(rayby<raybz)?rayby:raybz;

    for(int jt=1;jt<=nbbordz-1;jt++)
    {
		for(int it=0;it<=nbbordy;it++)
		{
				LIST_X[nbd] = tol	;
				LIST_Y[nbd] = 2.*it*rayby;
				LIST_Z[nbd] = 2.*jt*raybz;
				LIST_R[nbd] = raybe	;
				pbd+=4.*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]/3.;

				int NX= ((int) LIST_X[nbd])*vecsizex/H_TOT;
				if(NX<0) NX=0;
				if(NX>=vecsizex) NX=vecsizex-1;

				int NY= ((int) LIST_Y[nbd])*vecsizey/V_TOT;
				if(NY<0) NY=0;
				if(NY>=vecsizey) NY=vecsizey-1;

				int NZ= ((int) LIST_Z[nbd])*vecsizez/Z_TOT;
				if(NZ<0) NZ=0;
				if(NZ>=vecsizez) NZ=vecsizez-1;

				coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(nbd);
				LIST_C[nbd] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;
				LIST_A[nbd] = 0.;

				LIST_VX[nbd] = 0.;
				LIST_VY[nbd] = 0.;
				LIST_VZ[nbd] = 0.;

				LIST_AX[nbd] = 0.;
				LIST_AY[nbd] = 0.;
				LIST_AZ[nbd] = 0.;

				FX[nbd] = 0.;
				FY[nbd] = 0.;
				FZ[nbd] = 0.;

				LIST_M[nbd]  = dens*(4./3)*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]*10.;

				nbd++;
		}
   }

    for(int jt=1;jt<=nbbordz-1;jt++)
    {
		for(int it=0;it<=nbbordy;it++)
		{
				LIST_X[nbd] = R(H_TOT)-tol	;
				LIST_Y[nbd] = 2.*it*rayby;
				LIST_Z[nbd] = 2.*jt*raybz;
				LIST_R[nbd] = raybe	;
				pbd+=4.*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]/3.;

				int NX= ((int) LIST_X[nbd])*vecsizex/H_TOT;
				if(NX<0) NX=0;
				if(NX>=vecsizex) NX=vecsizex-1;

				int NY= ((int) LIST_Y[nbd])*vecsizey/V_TOT;
				if(NY<0) NY=0;
				if(NY>=vecsizey) NY=vecsizey-1;

				int NZ= ((int) LIST_Z[nbd])*vecsizez/Z_TOT;
				if(NZ<0) NZ=0;
				if(NZ>=vecsizez) NZ=vecsizez-1;

				coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(nbd);
				LIST_C[nbd] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;
				LIST_A[nbd] = 0.;

				LIST_VX[nbd] = 0.;
				LIST_VY[nbd] = 0.;
				LIST_VZ[nbd] = 0.;

				LIST_AX[nbd] = 0.;
				LIST_AY[nbd] = 0.;
				LIST_AZ[nbd] = 0.;

				FX[nbd] = 0.;
				FY[nbd] = 0.;
				FZ[nbd] = 0.;

				LIST_M[nbd]  = dens*(4./3)*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]*10.;

				nbd++;
		}
   }

    // Bords XZ-/+
    raybe=(raybx<raybz)?raybx:raybz;

    for(int jt=1;jt<=nbbordz-1;jt++)
    {
		for(int it=1;it<=nbbordx-1;it++)
		{
				LIST_X[nbd] = 2.*it*raybx;
				LIST_Y[nbd] = tol;
				LIST_Z[nbd] = 2.*jt*raybz;
				LIST_R[nbd] = raybe	;
				pbd+=4.*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]/3.;

				int NX= ((int) LIST_X[nbd])*vecsizex/H_TOT;
				if(NX<0) NX=0;
				if(NX>=vecsizex) NX=vecsizex-1;

				int NY= ((int) LIST_Y[nbd])*vecsizey/V_TOT;
				if(NY<0) NY=0;
				if(NY>=vecsizey) NY=vecsizey-1;

				int NZ= ((int) LIST_Z[nbd])*vecsizez/Z_TOT;
				if(NZ<0) NZ=0;
				if(NZ>=vecsizez) NZ=vecsizez-1;

				coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(nbd);
				LIST_C[nbd] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;
				LIST_A[nbd] = 0.;

				LIST_VX[nbd] = 0.;
				LIST_VY[nbd] = 0.;
				LIST_VZ[nbd] = 0.;

				LIST_AX[nbd] = 0.;
				LIST_AY[nbd] = 0.;
				LIST_AZ[nbd] = 0.;

				FX[nbd] = 0.;
				FY[nbd] = 0.;
				FZ[nbd] = 0.;

				LIST_M[nbd]  = dens*(4./3)*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]*10.;

				nbd++;
		}
   }

    for(int jt=1;jt<=nbbordz-1;jt++)
    {
		for(int it=1;it<=nbbordx-1;it++)
		{
				LIST_X[nbd] = 2.*it*raybx;
				LIST_Y[nbd] = R(V_TOT)-tol;
				LIST_Z[nbd] = 2.*jt*raybz;
				LIST_R[nbd] = raybe	;
				pbd+=4.*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]/3.;

				int NX= ((int) LIST_X[nbd])*vecsizex/H_TOT;
				if(NX<0) NX=0;
				if(NX>=vecsizex) NX=vecsizex-1;

				int NY= ((int) LIST_Y[nbd])*vecsizey/V_TOT;
				if(NY<0) NY=0;
				if(NY>=vecsizey) NY=vecsizey-1;

				int NZ= ((int) LIST_Z[nbd])*vecsizez/Z_TOT;
				if(NZ<0) NZ=0;
				if(NZ>=vecsizez) NZ=vecsizez-1;

				coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(nbd);
				LIST_C[nbd] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;
				LIST_A[nbd] = 0.;

				LIST_VX[nbd] = 0.;
				LIST_VY[nbd] = 0.;
				LIST_VZ[nbd] = 0.;

				LIST_AX[nbd] = 0.;
				LIST_AY[nbd] = 0.;
				LIST_AZ[nbd] = 0.;

				FX[nbd] = 0.;
				FY[nbd] = 0.;
				FZ[nbd] = 0.;

				LIST_M[nbd]  = dens*(4./3)*Pi*LIST_R[nbd]*LIST_R[nbd]*LIST_R[nbd]*10.;

				nbd++;
		}
   }


  pbd=pbd/2.;



    bool boolout;
    R dist,dx,dy,dz,erfm1,pg;
    int it,jt,kt,kjt,kk;
    int NX,NY,NZ,NN;

    R mas0=dens*(4./3)*Pi*R_SPH*R_SPH*R_SPH;
    int vecsizexy=vecsizex*vecsizey;

	 R coef1 =cov*mu*sqrt(2)*0.5*sqrt(Pi);
	 R coef3 =cov*mu*sqrt(2)*0.5*sqrt(Pi)*Pi/12;
	 R coef5 =cov*mu*sqrt(2)*0.5*sqrt(Pi)*7*pow(Pi,2)/480;
	 R coef7 =cov*mu*sqrt(2)*0.5*sqrt(Pi)*127*pow(Pi,3)/40320;
	 R coef9 =cov*mu*sqrt(2)*0.5*sqrt(Pi)*4369*pow(Pi,4)/5806080;
	 R coef11 =cov*mu*sqrt(2)*0.5*sqrt(Pi)*34807*pow(Pi,5)/182476800;

	// DISTRIBUTION INITIALE DE SPHERES
	for(it=nbd;it<NB_SPH;it++){
        if(it%50000==0){  cout<<"it :"<<it<<"/"<<NB_SPH<<endl;}

	  // PARAMETRES ALEATOIRES
		struct timeval tv ;
		gettimeofday(&tv, NULL) ;
		srand(tv.tv_usec) ;

		boolout=0;
		while(boolout==0){
		boolout=1;
		LIST_X[it] = (((R) (rand()%int((H_TOT-(2.*crit+2e-2))*1000000)))/1000000.)+crit+1e-2;
		LIST_Y[it] = (((R) (rand()%int((V_TOT-(2.*crit+2e-2))*1000000)))/1000000.)+crit+1e-2;
		LIST_Z[it] = (((R) (rand()%int((Z_TOT-(2.*crit+2e-2))*1000000)))/1000000.)+crit+1e-2;

		NX= ((int) LIST_X[it])*vecsizex/H_TOT;
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;
		NY= ((int) LIST_Y[it])*vecsizey/V_TOT;
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;
		NZ= ((int) LIST_Z[it])*vecsizez/Z_TOT;
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;

		NN=NX+NY*vecsizex+NZ*vecsizexy;

/////////////////////
/*
				jt=0;
				while((boolout==1)&&(jt<it)){

					    dx=(LIST_X[it]-LIST_X[jt]);
					    dy=(LIST_Y[it]-LIST_Y[jt]);
					    dz=(LIST_Z[it]-LIST_Z[jt]);

					dist=dx*dx+dy*dy+dz*dz;
					dist=sqrt(dist);
					if(dist<=(R_SPH+LIST_R[jt]+1e-2*R_SPH)){
					boolout=0;
					}
					jt++;
				}

*/
////////////////////////

			kt=0;
			while((boolout==1)&&(kt<numc2[NN].size())){

			         kk=numc2[NN][kt];
				 kjt=0;
				 while((boolout==1)&&(kjt<coul[kk].size())){
					 jt=coul[kk][kjt];

					  if(jt!=it){

						dx=(LIST_X[it]-LIST_X[jt]);
						dy=(LIST_Y[it]-LIST_Y[jt]);
						dz=(LIST_Z[it]-LIST_Z[jt]);

						dist=dx*dx+dy*dy+dz*dz;
						dist=sqrt(dist);
						if(dist<=(R_SPH+LIST_R[jt]+1e-2*R_SPH)) boolout=0;
					
					  }
					kjt++;
			          } // fin kjt
				kt++;
			 } //fin kt
 



////////////////


		}

		coul[NN].push_back(it);
		LIST_C[it] = NN;
		LIST_R[it] = R_SPH;

		pg=rand()%98;
                pg++;
		pg/=100;
		pg=2*pg-1;
		erfm1=coef1*pg+coef3*pow(pg,3)+coef5*pow(pg,5)+coef7*pow(pg,7)+coef9*pow(pg,9)+coef11*pow(pg,11);
		LIST_A[it] = mu+erfm1;
		if(LIST_A[it]<mu/2.){LIST_A[it]=mu/2.;}
		else if(LIST_A[it]>(150./100.*mu)){LIST_A[it]=(150./100.*mu);}

		LIST_VX[it] = ((R) (rand()%(100)))/50.;
		LIST_VY[it] = ((R) (rand()%(100)))/50.;
		LIST_VZ[it] = ((R) (rand()%(100)))/50.;
		LIST_VX[it]--;
		LIST_VY[it]--;
		LIST_VZ[it]--;

		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;

		FX[it] = 0.;
		FY[it] = 0.;
		FZ[it] = 0.;

		LIST_M[it] = mas0;

	}

numc2.clear();
}


inline void selco_corci(int & nbco,  int ** CONT, int vecsize, int vecsizex, int vecsizey, int vecsizez, int H_TOT, int V_TOT, int Z_TOT, vector< vector<int> > coul,vector< vector<int> > numc, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * DCONT, R ** NCONT)
{

nbco=0;

R dx,dy,dz,nd2,r2;
int ii,lt,it,kt,kjt,jt,nbco1;

 //   cout<<"vecsizex : "<< vecsizex <<endl;

# pragma omp parallel for schedule(dynamic,int(vecsize/NB_THREADS)) private (dx,dy,dz,nd2,r2,ii,lt,it,kt,kjt,jt,nbco1) 
	for(ii=0;ii<vecsize;ii++){
   //   cout<<ii<<", "<<vecsize<<", "<< vecsizex<<", "<< coul[ii].size()<<", "<< numc[ii].size()<<endl;

		for(lt=0;lt<coul[ii].size();lt++){

		it = coul[ii][lt];
		//    cout<<"it :"<< it<<" ("<<(lt+1)<<"/"<<coul[ii].size()<<")"<<endl;

			for(kt=0;kt<numc[ii].size();kt++){

				 for(kjt=0;kjt<coul[numc[ii][kt]].size();kjt++){
					 jt=coul[numc[ii][kt]][kjt];

			  //        cout<<"--jt :"<< jt<<" ("<<(kt+1)<<"/"<<numc.size()<<")"<<" ("<<(kjt+1)<<"/"<<coul[numc[kt]].size()<<")"<<endl;


					  if((kt>0)||((kt==0)&&(jt>it))){
					  dx=LIST_X[it]-LIST_X[jt];
					  dy=LIST_Y[it]-LIST_Y[jt];
					  dz=LIST_Z[it]-LIST_Z[jt];
					  nd2=dx*dx+dy*dy+dz*dz;
					  r2=(LIST_R[it]+LIST_R[jt])*(LIST_R[it]+LIST_R[jt]);

						  //collision
						  if(nd2<r2+tol)
						  {
		                                      # pragma omp critical
						       { 
						       nbco1=nbco;
						       nbco++;
						       }

						       CONT[nbco1][0]=it;
						       CONT[nbco1][1]=jt;
						       DCONT[nbco1]=sqrt(nd2);

						       // Normale

						       NCONT[nbco1][0]=dx/sqrt(nd2);
						       NCONT[nbco1][1]=dy/sqrt(nd2);
						       NCONT[nbco1][2]=dz/sqrt(nd2);

						  }

					  }

			     }

			}

		}

	}

}

inline void reset_size(int nbd, R pbd, int & vecsize, int & vecsizex, int & vecsizey, int & vecsizez, R & PR_EST,vector< vector<int>  > & coul,vector< vector<int>  > & numc,int * LIST_C,int H_TOT,int V_TOT,int Z_TOT,int nD,R dt,R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_A, R & raymax, R amax)
{
	R rmax=raymax;//+dt*amax;
  //

	/*R Reff = ((2*rmax+1e-8)>(2*H_TOT/pow(nD,1./3)))?(2*rmax+1e-8):(2*H_TOT/pow(nD,1./3));
	vecsizex=(floor(H_TOT/Reff)>4)?floor(H_TOT/Reff):4;
    Reff = ((2*rmax+1e-8)>(2*V_TOT/pow(nD,1./3)))?(2*rmax+1e-8):(2*V_TOT/pow(nD,1./3));
	vecsizey=(floor(V_TOT/Reff)>4)?floor(V_TOT/Reff):4;
    Reff = ((2*rmax+1e-8)>(2*Z_TOT/pow(nD,1./3)))?(2*rmax+1e-8):(2*Z_TOT/pow(nD,1./3));
	vecsizez=(floor(Z_TOT/Reff)>4)?floor(Z_TOT/Reff):4;*/

	R Reff = 2*rmax;
	vecsizex=(floor(H_TOT/Reff)>4)?(floor(H_TOT/Reff)-1):4;
	vecsizex-=2;
	Reff =   2*rmax;
	vecsizey=(floor(V_TOT/Reff)>4)?(floor(V_TOT/Reff)-1):4;
	vecsizey-=2;
	Reff =   2*rmax;
	vecsizez=(floor(Z_TOT/Reff)>4)?(floor(Z_TOT/Reff)-1):4;
	vecsizez-=2;

	vecsize=vecsizex*vecsizey*vecsizez;

	coul.clear();
	for (int i = 0; i < vecsize; i++) {
		coul.push_back(vector<int>()); // Add an empty row
	}
	numc.clear();
	for (int ii = 0; ii < vecsize; ii++) {
		numc.push_back(vector<int>()); // Add an empty row
    }

	for (int ii = 0; ii < vecsize; ii++) {
		int NX  = ii%vecsizex;
		int NY  = (ii%(vecsizex*vecsizey))/vecsizex;
		int NZ  = ii/(vecsizex*vecsizey);

		int NXP1= (NX+vecsizex+1)%vecsizex;
		int NXM1= (NX+vecsizex-1)%vecsizex;
		int NYP1= (NY+vecsizey+1)%vecsizey;
		int NYM1= (NY+vecsizey-1)%vecsizey;
		int NZP1= (NZ+vecsizez+1)%vecsizez;
		int NZM1= (NZ+vecsizez-1)%vecsizez;

		numc[ii].push_back(ii);
		if(NZM1==(NZ-1)){
		 if((NYM1==(NY-1))&&(NXM1==(NX-1))){numc[ii].push_back(NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
		 if(NYM1==(NY-1)){numc[ii].push_back(NX+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
		 if((NYM1==(NY-1))&&(NXP1==(NX+1))){numc[ii].push_back(NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZM1);}
		 if(NXM1==(NX-1)){numc[ii].push_back(NXM1+vecsizex*NY+vecsizex*vecsizey*NZM1);}
		 numc[ii].push_back(NX+vecsizex*NY+vecsizex*vecsizey*NZM1);
		 if(NXP1==(NX+1)){numc[ii].push_back(NXP1+vecsizex*NY+vecsizex*vecsizey*NZM1);}
		 if((NYP1==(NY+1))&&(NXM1==(NX-1))){numc[ii].push_back(NXM1+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
		 if(NYP1==(NY+1)){numc[ii].push_back(NX+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
		 if((NYP1==(NY+1))&&(NXP1==(NX+1))){numc[ii].push_back(NXP1+vecsizex*NYP1+vecsizex*vecsizey*NZM1);}
		}

		if(NYM1==(NY-1)){
		 if(NXM1==(NX-1)){numc[ii].push_back((NXM1+vecsizex*NYM1+vecsizex*vecsizey*NZ));}
		 numc[ii].push_back((NX+vecsizex*NYM1+vecsizex*vecsizey*NZ));
		}
		if(NXP1==(NX+1)){
		 if(NYM1==(NY-1)){numc[ii].push_back((NXP1+vecsizex*NYM1+vecsizex*vecsizey*NZ));}
		 numc[ii].push_back((NXP1+vecsizex*NY+vecsizex*vecsizey*NZ));
		}
	 }

	 PR_EST=pbd;

int it,NX,NY,NZ;
# pragma omp parallel for schedule(dynamic,int(nD/NB_THREADS)) private (it,NX,NY,NZ) reduction(+:PR_EST)
      for(it=0;it<nD;it++){
      //Volume fraction
	if(it>=nbd){
	PR_EST+=(4./3)*3.14159*LIST_R[it]*LIST_R[it]*LIST_R[it];
	}
     
      // update
      NX= floor((LIST_X[it]/H_TOT)*vecsizex);
      NY= floor((LIST_Y[it]/V_TOT)*vecsizey);
      NZ= floor((LIST_Z[it]/Z_TOT)*vecsizez);
      if(NX<0) NX=0;
      if(NY<0) NY=0;
      if(NZ<0) NZ=0;
      if(NX>=vecsizex) NX=vecsizex-1;
      if(NY>=vecsizey) NY=vecsizey-1;
      if(NZ>=vecsizez) NZ=vecsizez-1;

      if((NX+NY*vecsizex+NZ*vecsizex*vecsizey)>=(vecsizex*vecsizey*vecsizez)) cout<<"aie!"<<endl;
	# pragma omp critical
	{ 
	coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(it);
	}

	LIST_C[it]=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
      }
      PR_EST/=(Z_TOT*V_TOT*H_TOT);

cout<<"vecsize:"<<vecsizex<<", "<<vecsizey<<", "<<vecsizez<<endl;

}

inline void eval_v2(int nbd, R & v2, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_R, int NB_SPH){

fstream g;
g.open("eval_v2.txt",fstream::out);
g.precision(4);

g<<NB_SPH<<endl;
g<<nbd<<endl;
for(int i=0;i<NB_SPH;i++)
{
 g<<setw(10)<<(i+1)<<setw(18)<<scientific<<LIST_VX[i]<<setw(18)<<scientific<<LIST_VY[i]<<setw(18)<<scientific<<LIST_VZ[i]<<endl;
}

g.close();


v2=0.;

	for(int it=nbd;it<NB_SPH;it++){
    v2+=(LIST_VX[it]*LIST_VX[it])+(LIST_VY[it]*LIST_VY[it])+(LIST_VZ[it]*LIST_VZ[it]);
    }
    v2/=(NB_SPH-nbd);

}

inline void rayonint(R & rim, int NBCO, int ** CONT, R * DCONT, R * LIST_R){

rim=0.;
int dis1,dis2,nim;
R ri,r1,r2,d;
nim=0;

	for(int it=0;it<NBCO;it++){

	dis1=CONT[it][0];
	dis2=CONT[it][1];
	r1=LIST_R[dis1];
	r2=LIST_R[dis2];
	d=DCONT[it];

	ri=(r1*r1-((r1*r1-r2*r2+d*d)/(2*d))*((r1*r1-r2*r2+d*d)/(2*d)));
	if(ri>0){
        ri=sqrt(ri);
	ri=2*ri/(r1+r2);
	rim+=ri;
        nim++;
        }
    }
    rim/=nim;
}

inline void force( R dt, R kn, R kt, int NB_SPH, int NBCO,  int ** CONT, R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * DCONT, R ** NCONT,  R * FX, R * FY, R * FZ){

int numc1,numc2;
R dist,n1,n2,n3,rn,rnv,mas1,pred1,u1;
int it;

# pragma omp parallel for schedule(dynamic,int(NBCO/NB_THREADS)) private (it,numc1,numc2,dist,n1,n2,n3,rn,rnv,mas1,pred1,u1) reduction (+:FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH]) 
	for(it=0;it<NBCO;it++){

	     numc1=CONT[it][0];
	     numc2=CONT[it][1];
	     dist=DCONT[it];
	     mas1=LIST_M[numc1]*LIST_M[numc2]/(LIST_M[numc1]+LIST_M[numc2]);

	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];

	     u1=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);

	     rn=kn*(dist-(LIST_R[numc1]+LIST_R[numc2]));
	     rnv=sqrt(kn*mas1)*u1;
	     pred1=-rn-6.*rnv;


	     if(pred1<0){ //relachement
	     pred1=0.;
             }


	     FX[numc1]=FX[numc1]+n1*pred1;
	     FY[numc1]=FY[numc1]+n2*pred1;
	     FZ[numc1]=FZ[numc1]+n3*pred1;

	     FX[numc2]=FX[numc2]+(-(n1*pred1));
	     FY[numc2]=FY[numc2]+(-(n2*pred1));
	     FZ[numc2]=FZ[numc2]+(-(n3*pred1));

       }


}

inline void verlet(int nbd, int vecsize, int vecsizex, int vecsizey, int vecsizez, vector< vector<int>  > & coul,int * LIST_C,int H_TOT,int V_TOT,int Z_TOT, int nD,R dt,R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_A, R * FX, R * FY, R * FZ, R & raymax)
{

for(int it=0;it<vecsize;it++){
coul[it].clear();
}

raymax=0.;

int it,NX,NY,NZ;
# pragma omp parallel for schedule(dynamic,int(nbd/NB_THREADS)) private (it,NX,NY,NZ) 
for(it=0;it<nbd;it++){

	FX[it] = 0.;
	FY[it] = 0.;
	FZ[it] = 0.;

	NX= floor((LIST_X[it]/H_TOT)*vecsizex);
	NY= floor((LIST_Y[it]/V_TOT)*vecsizey);
	NZ= floor((LIST_Z[it]/Z_TOT)*vecsizez);
	if(NX<0) NX=0;
	if(NY<0) NY=0;
	if(NZ<0) NZ=0;
	if(NX>=vecsizex) NX=vecsizex-1;
	if(NY>=vecsizey) NY=vecsizey-1;
	if(NZ>=vecsizez) NZ=vecsizez-1;

	# pragma omp critical
	{
	coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(it);
	}

	LIST_C[it]=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
        if(raymax<LIST_R[it]) raymax=LIST_R[it];
}
# pragma omp parallel for schedule(dynamic,int((nD-nbd)/NB_THREADS)) private (it,NX,NY,NZ) 
for(it=nbd;it<nD;it++){

	 LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);
	 LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	 LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);

	 LIST_AX[it] = FX[it]/LIST_M[it];
	 LIST_AY[it] = FY[it]/LIST_M[it];
	 LIST_AZ[it] = FZ[it]/LIST_M[it];

	 LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
	 LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
	 LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];

	 FX[it] = 0.;
	 FY[it] = 0.;
	 FZ[it] = 0.;

      NX= floor((LIST_X[it]/H_TOT)*vecsizex);
      NY= floor((LIST_Y[it]/V_TOT)*vecsizey);
      NZ= floor((LIST_Z[it]/Z_TOT)*vecsizez);
      if(NX<0) NX=0;
      if(NY<0) NY=0;
      if(NZ<0) NZ=0;
      if(NX>=vecsizex) NX=vecsizex-1;
      if(NY>=vecsizey) NY=vecsizey-1;
      if(NZ>=vecsizez) NZ=vecsizez-1;

	# pragma omp critical
	{
	coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey].push_back(it);
	}
      LIST_C[it]=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
      LIST_R[it] = LIST_R[it]+LIST_A[it]*dt;
      if(raymax<LIST_R[it]) raymax=LIST_R[it];
}

}


inline void ExpMeshHAMZA(R SIZEX, R SIZEY, R SIZEZ, int NB_SPH, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R){

fstream g;
g.open("pristimap",fstream::out);
g.precision(4);

g<<NB_SPH<<endl;
for(int i=0;i<NB_SPH;i++)
{
 g<<setw(10)<<(i+1)<<setw(18)<<scientific<<LIST_X[i]<<setw(18)<<scientific<<LIST_Y[i]<<setw(18)<<scientific<<LIST_Z[i]<<setw(18)<<scientific<<LIST_R[i]<<endl;
}

g.close();
}



inline void ExpVISU(int NB_SPH, int H_TOT, int V_TOT, int Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R,bool bout){

fstream g;
g.precision(4);
if(!bout){
g.open("visu",fstream::out | fstream::app);
}
else{
g.open("visu",fstream::out);
}
g<<setw(12)<<NB_SPH<<setw(12)<<H_TOT<<endl;
for(int it=0;it<NB_SPH;it++){
g<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<endl;
}
g.close();
}

inline void ExpMeshMULTICOR(R SIZEX, R SIZEY, R SIZEZ, int NB_SPH, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R){

R H_MIN=0.;

fstream g;
g.open("carte",fstream::out);
g.precision(4);

g<<"CONDLIMI"<<endl;
g<<setw(12)<<scientific<<left<<H_MIN<<setw(12)<<scientific<<left<<SIZEX<<setw(12)<<scientific<<left<<H_MIN<<setw(12)<<scientific<<left<<SIZEY<<setw(12)<<scientific<<left<<H_MIN<<setw(12)<<scientific<<left<<SIZEZ<<endl;
g<<"SNUMERI"<<endl;
g<<setw(12)<<scientific<<0.000001<<setw(12)<<scientific<<0.000001<<endl;
g<<"SCONTAC"<<endl;
g<<setw(12)<<scientific<<0.000005<<setw(14)<<left<<999999999<<endl;
g<<"MATER"<<endl;
g<<setw(12)<<scientific<<9.81<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<210e9<<setw(12)<<scientific<<0.02<<setw(12)<<scientific<<7.8e9<<endl;

g<<setw(10)<<"SPHERE"<<NB_SPH<<endl;
for(int i=0;i<NB_SPH;i++)
{
 g<<setw(10)<<(i+1)<<setw(12)<<scientific<<LIST_R[i]<<endl;
 g<<setw(14)<<scientific<<LIST_X[i]<<setw(14)<<scientific<<LIST_Y[i]<<setw(14)<<scientific<<LIST_Z[i]<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
 g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
 g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
 g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
}

g<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(10)<<1<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<endl;

g<<setw(10)<<"HSINUS"<<1<<endl;
g<<setw(10)<<1<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<endl;

g<<setw(10)<<"HOBLIQUE"<<1<<endl;
g<<setw(10)<<1<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<endl;
g<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<setw(12)<<scientific<<0.<<endl;

g<<setw(10)<<"PAROI"<<6<<endl;
g<<setw(10)<<1<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<2<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<3<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<4<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<5<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(10)<<6<<setw(10)<<"HAFFINE"<<1<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<SIZEZ<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<SIZEX<<setw(14)<<scientific<<SIZEY<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<setw(14)<<scientific<<0.<<endl;
g.close();
}

int main(int argc, char *argv[])
{
    const rlim_t kStackSize = 128 * 2048 * 2048;   // min stack size = 32 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }

R t1=CPUtime();
R t1p=omp_get_wtime () ;

// Parallélisation
omp_set_num_threads(NB_THREADS);

// DISCRETISATIONS
string str_H_TOT = argv[1];
int H_TOT;
from_string<int>(H_TOT,str_H_TOT,dec);

string str_V_TOT = argv[2];
int V_TOT;
from_string<int>(V_TOT,str_V_TOT,dec);

string str_Z_TOT = argv[3];
int Z_TOT;
from_string<int>(Z_TOT,str_Z_TOT,dec);

// DIMENSIONS
string str_SIZEX = argv[4];
R SIZEX;
from_string<R>(SIZEX,str_SIZEX,dec);

string str_SIZEY = argv[5];
R SIZEY;
from_string<R>(SIZEY,str_SIZEY,dec);

string str_SIZEZ = argv[6];
R SIZEZ;
from_string<R>(SIZEZ,str_SIZEZ,dec);

// NOMBRE DE SPHERES
string str_NB_SPH = argv[7];
int NB_SPH;
from_string<int>(NB_SPH,str_NB_SPH,dec);

// FRACTION VOLUMIQUE DE SPHERES ATTENDUE
string str_PR_SPH = argv[8];
R PR_SPH;
from_string<R>(PR_SPH,str_PR_SPH,dec);
R pD=PR_SPH;

R mu=10./pow(NB_SPH,(1./3));

R R_SPH=pow((3.*H_TOT*V_TOT*Z_TOT*0.12/(4.*Pi*NB_SPH)),1./3);
R R_SPH2=pow((3.*H_TOT*V_TOT*Z_TOT/(4.*Pi*NB_SPH)),1./3);
R kn=1e+7;
R kt=1e+7;
int nbd;
R pbd;
int NBCO;

for(int i=0;i<NBCONT;i++){
CONT[i]= new int[2];
}
for(int i=0;i<NBCONT;i++){
NCONT[i]= new R[3];
}

// Initialisation
init_size(vecsize, vecsizex, vecsizey, vecsizez, coul, numc,numc2, H_TOT, V_TOT, Z_TOT, NB_SPH, R_SPH2);

init_sph(mu,nbd,pbd,vecsizex,vecsizey,vecsizez,PR_SPH,H_TOT,V_TOT,Z_TOT,NB_SPH,R_SPH,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_VX,LIST_VY,LIST_VZ,LIST_AX,LIST_AY,LIST_AZ,LIST_C,LIST_A,LIST_M,FX,FY,FZ);
//ExpVISU(NB_SPH,H_TOT,V_TOT,Z_TOT,LIST_X,LIST_Y,LIST_Z,LIST_R,1);

R gap=1.;
R PR_AV=0.;
R PR_EST=0.;
R ti=0.;
int ite=0;
R amax;
R raymax;
maxi(amax,LIST_A,NB_SPH);
maxi(raymax,LIST_R,NB_SPH);

cout<<"Début du process"<<endl;

while((PR_EST<pD)&&(ite<10000000)){
  selco_corci(NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,H_TOT,V_TOT,Z_TOT,coul,numc,LIST_R,LIST_X,LIST_Y,LIST_Z,DCONT,NCONT);
  force(dt,kn,kt,NB_SPH,NBCO,CONT,LIST_M,LIST_R,LIST_X,LIST_Y,LIST_Z,DCONT,NCONT,FX,FY,FZ);
  verlet(nbd,vecsize,vecsizex,vecsizey,vecsizez,coul,LIST_C,H_TOT,V_TOT,Z_TOT,NB_SPH,dt,LIST_M,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_VX,LIST_VY,LIST_VZ,LIST_AX,LIST_AY,LIST_AZ,LIST_A,FX,FY,FZ,raymax);

  if(ite%20==0){

	reset_size(nbd,pbd,vecsize,vecsizex,vecsizey,vecsizez,PR_EST,coul,numc,LIST_C,H_TOT,V_TOT,Z_TOT,NB_SPH,dt,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_A,raymax,amax); 

	PR_EST=pbd;
	for(int it=nbd;it<NB_SPH;it++){
	PR_EST+=(4./3)*3.14159*LIST_R[it]*LIST_R[it]*LIST_R[it];
	}
	PR_EST/=(H_TOT*V_TOT*Z_TOT);
	gap=abs(PR_AV-PR_EST)/PR_EST;
	PR_AV=PR_EST;

	R t2=CPUtime()-t1;	 
	R t2p=omp_get_wtime () - t1p; 
	cout<<"Iteration : "<<ite<<", Volume fraction : "<<PR_EST<<", err :"<<gap<<", CPU time :"<<t2<<", Par. time :"<<t2p<<endl;
	cout<<"NBCO : "<<NBCO<<", NBD : "<<nbd<<endl;

 // ExpVISU(NB_SPH,H_TOT,V_TOT,Z_TOT,LIST_X,LIST_Y,LIST_Z,LIST_R,0);
  }

ite++;
ti+=dt;
}

PR_EST=pbd;
for(int it=nbd;it<NB_SPH;it++){
      //Area fraction
      PR_EST+=(4./3)*3.14159*LIST_R[it]*LIST_R[it]*LIST_R[it];
}
PR_EST/=(H_TOT*V_TOT*Z_TOT);
cout<<"PR_EST av :"<<PR_EST <<endl;

R PR_OLD=PR_EST;
PR_EST=pbd;
for(int it=nbd;it<NB_SPH;it++){
LIST_A[it]=0.;

LIST_R[it]=pow((PR_SPH-pbd/(H_TOT*V_TOT*Z_TOT))/(PR_OLD-pbd/(H_TOT*V_TOT*Z_TOT)),1./3)*LIST_R[it];
      //Volume fraction
      PR_EST+=(4./3)*3.14159*LIST_R[it]*LIST_R[it]*LIST_R[it];
}
PR_EST/=(H_TOT*V_TOT*Z_TOT);
cout<<"PR_EST ap :"<<PR_EST <<endl;

reset_size(nbd,pbd,vecsize,vecsizex,vecsizey,vecsizez,PR_EST,coul,numc,LIST_C,H_TOT,V_TOT,Z_TOT,NB_SPH,dt,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_A,raymax,amax); 

ite=0;
R v2,rim;
eval_v2(nbd,v2,LIST_VX,LIST_VY,LIST_VZ,LIST_X,LIST_Y,LIST_Z,LIST_R,NB_SPH);
rayonint(rim,NBCO,CONT,DCONT,LIST_R);
kn=1e+7;
kt=1e+7;
cout<<"v2 : "<<v2<<", rim : "<<rim<<endl;

while(rim>2e-2){

  selco_corci(NBCO,CONT,vecsize,vecsizex,vecsizey,vecsizez,H_TOT,V_TOT,Z_TOT,coul,numc,LIST_R,LIST_X,LIST_Y,LIST_Z,DCONT,NCONT);
  force(dt,kn,kt,NB_SPH,NBCO,CONT,LIST_M,LIST_R,LIST_X,LIST_Y,LIST_Z, DCONT,NCONT,FX,FY,FZ);
  verlet(nbd,vecsize,vecsizex,vecsizey,vecsizez,coul,LIST_C,H_TOT,V_TOT,Z_TOT,NB_SPH,dt,LIST_M,LIST_R,LIST_X,LIST_Y,LIST_Z,LIST_VX,LIST_VY,LIST_VZ,LIST_AX,LIST_AY,LIST_AZ,LIST_A,FX,FY,FZ,raymax);

  if(ite%100==0){

  eval_v2(nbd,v2,LIST_VX,LIST_VY,LIST_VZ,LIST_X,LIST_Y,LIST_Z,LIST_R,NB_SPH);
  rayonint(rim,NBCO,CONT,DCONT,LIST_R);

  R t2=CPUtime()-t1;	 
  R t2p=omp_get_wtime () - t1p; 

  cout<<"Ite : "<<ite<<", area : "<<PR_EST<<", err :"<<gap<<", v2 :"<<v2<<", rim :"<<rim<<", CPU time :"<<t2<<", Par. time :"<<t2p<<endl;
cout<<"NBCO : "<<NBCO<<", NBD : "<<nbd<<endl;

  //ExpVISU(NB_SPH,H_TOT,V_TOT,Z_TOT,LIST_X,LIST_Y,LIST_Z,LIST_R,0);
  }


ite++;
ti+=dt;
}

cout<<"Iterations total number :"<<ite<<endl;
cout<<"Total area fraction :"<<PR_EST<<endl;
cout<<"Total CPU time :"<<CPUtime()-t1<<endl;
cout<<"Total par. time :"<<omp_get_wtime () - t1p<<endl;

 // Sortie data
fstream g;
g.open("data.csv",fstream::out);
g<<"RAY"<<endl;
for(int i=0;i<NB_SPH;i++){
g<<LIST_R[i]<<endl;
}
g.close();

// CONVERT
for(int i=0;i<NB_SPH;i++){
LIST_X[i]=LIST_X[i]*SIZEX/R(H_TOT);
LIST_Y[i]=LIST_Y[i]*SIZEY/R(V_TOT);
LIST_Z[i]=LIST_Z[i]*SIZEZ/R(Z_TOT);
LIST_R[i]=LIST_R[i]*SIZEX/R(H_TOT);
}

// Export MULTICOR
//ExpMeshMULTICOR(SIZEX,SIZEY,SIZEZ,NB_SPH,LIST_X,LIST_Y,LIST_Z,LIST_R);
ExpMeshHAMZA(SIZEX,SIZEY,SIZEZ,NB_SPH,LIST_X,LIST_Y,LIST_Z,LIST_R);
cout<<"Export MULTICOR ok"<<endl;

// free memory

delete [] CONT;
delete [] NCONT;

delete [] FX;
delete [] FY;
delete [] FZ;

delete [] DCONT;

delete [] LIST_M;
delete [] LIST_A;
delete [] LIST_VX;
delete [] LIST_VY;
delete [] LIST_VZ;
delete [] LIST_AX;
delete [] LIST_AY;
delete [] LIST_AZ;
delete [] LIST_C;
delete [] LIST_R;
delete [] LIST_X;
delete [] LIST_Y;
delete [] LIST_Z;

return 1;
}
