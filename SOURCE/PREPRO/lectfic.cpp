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

#include "lectfic.h"

void read_paroi(int & NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ)
{
ifstream f("APPLI/paroi/testp",ios::in);
f >> NB_PAR;

int tmp;

for (int it=0;it<NB_PAR;it++){
f >> tmp;	
	for (int jt=0;jt<4;jt++){
	f >> LIST_PX[it][jt];
	f >> LIST_PY[it][jt];
	f >> LIST_PZ[it][jt];
	}
}

f.close();
}


void read_sph(R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX)
{
R_SPH=0.;
RMAX=0.;  

H_TOT=-2e9;
V_TOT=-2e9;
Z_TOT=-2e9;

H_POS=2e9;
V_POS=2e9;
Z_POS=2e9;

ifstream f("APPLI/carttest/test",ios::in);
f >> NB_SPH;

int tmp;

for (int it=0;it<NB_SPH;it++){
f >> tmp;	
f >> LIST_X[it];
f >> LIST_Y[it];
f >> LIST_Z[it];
f >> LIST_R[it];
LIST_X[it]=LIST_X[it];
LIST_Y[it]=LIST_Y[it];
LIST_Z[it]=LIST_Z[it];
LIST_R[it]=LIST_R[it];
R_SPH+=LIST_R[it];

RMAX=(LIST_R[it]>RMAX)?LIST_R[it]:RMAX;
H_TOT=(LIST_X[it]>H_TOT)?LIST_X[it]:H_TOT;
V_TOT=(LIST_Y[it]>V_TOT)?LIST_Y[it]:V_TOT;
Z_TOT=(LIST_Z[it]>Z_TOT)?LIST_Z[it]:Z_TOT;
H_POS=(LIST_X[it]<H_POS)?LIST_X[it]:H_POS;
V_POS=(LIST_Y[it]<V_POS)?LIST_Y[it]:V_POS;
Z_POS=(LIST_Z[it]<Z_POS)?LIST_Z[it]:Z_POS;

}
R_SPH/=NB_SPH;

H_TOT=H_TOT-H_POS;
V_TOT=V_TOT-V_POS;
Z_TOT=Z_TOT-Z_POS;

f.close();
}

void read_sph_ind(R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX, R rind)
{
R_SPH=0.;
RMAX=0.;  

H_TOT=-2e9;
V_TOT=-2e9;
Z_TOT=-2e9;

H_POS=2e9;
V_POS=2e9;
Z_POS=2e9;

ifstream f("APPLI/carttest/test",ios::in);
f >> NB_SPH;

int tmp;

R vols=0.;

for (int it=0;it<NB_SPH;it++){
f >> tmp;	
f >> LIST_X[it];
f >> LIST_Y[it];
f >> LIST_Z[it];
f >> LIST_R[it];
LIST_X[it]=LIST_X[it];
LIST_Y[it]=LIST_Y[it];
LIST_Z[it]=LIST_Z[it];
LIST_R[it]=LIST_R[it];
R_SPH+=LIST_R[it];

vols=vols+4./3*3.14159*LIST_R[it]*LIST_R[it]*LIST_R[it];

RMAX=(LIST_R[it]>RMAX)?LIST_R[it]:RMAX;
H_TOT=((LIST_X[it])>H_TOT)?(LIST_X[it]):H_TOT;
V_TOT=((LIST_Y[it])>V_TOT)?(LIST_Y[it]):V_TOT;
Z_TOT=((LIST_Z[it])>Z_TOT)?(LIST_Z[it]):Z_TOT;
H_POS=(LIST_X[it]<H_POS)?LIST_X[it]:H_POS;
V_POS=(LIST_Y[it]<V_POS)?LIST_Y[it]:V_POS;
Z_POS=(LIST_Z[it]<Z_POS)?LIST_Z[it]:Z_POS;
}
R_SPH/=NB_SPH;

LIST_X[NB_SPH]=R(H_TOT)/2.;
LIST_Y[NB_SPH]=R(V_TOT)/2.;
LIST_Z[NB_SPH]=R(Z_TOT)+rind+R_SPH;
LIST_R[NB_SPH]=rind;
NB_SPH++;

f.close();

cout<<"Volume des particules :"<<vols<<endl;
}

void read_sphp_ind(R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX, R rind)
{
R_SPH=0.;
RMAX=0.;  

H_TOT=-2e9;
V_TOT=-2e9;
Z_TOT=-2e9;

H_POS=2e9;
V_POS=2e9;
Z_POS=2e9;

ifstream f("APPLI/carttest/test",ios::in);
f >> NB_SPH;

int tmp;

R vols=0.;

for (int it=0;it<NB_SPH;it++){
f >> tmp;	
f >> LIST_X[it];
f >> LIST_Y[it];
f >> LIST_Z[it];
f >> LIST_R[it];
LIST_X[it]=LIST_X[it];
LIST_Y[it]=LIST_Y[it];
LIST_Z[it]=LIST_Z[it];
LIST_R[it]=LIST_R[it];
R_SPH+=LIST_R[it];

vols=vols+4./3*3.14159*LIST_R[it]*LIST_R[it]*LIST_R[it];

RMAX=(LIST_R[it]>RMAX)?LIST_R[it]:RMAX;
H_TOT=((LIST_X[it]+LIST_R[it])>H_TOT)?(LIST_X[it]+LIST_R[it]):H_TOT;
V_TOT=((LIST_Y[it]+LIST_R[it])>V_TOT)?(LIST_Y[it]+LIST_R[it]):V_TOT;
Z_TOT=((LIST_Z[it]+LIST_R[it])>Z_TOT)?(LIST_Z[it]+LIST_R[it]):Z_TOT;
H_POS=((LIST_X[it]-LIST_R[it])<H_POS)?(LIST_X[it]-LIST_R[it]):H_POS;
V_POS=((LIST_Y[it]-LIST_R[it])<V_POS)?(LIST_Y[it]-LIST_R[it]):V_POS;
Z_POS=((LIST_Z[it]-LIST_R[it])<Z_POS)?(LIST_Z[it]-LIST_R[it]):Z_POS;
}
R_SPH/=NB_SPH;

LIST_X[NB_SPH]=R(H_TOT)/2.;
LIST_Y[NB_SPH]=R(V_TOT)/2.;
LIST_Z[NB_SPH]=R(Z_TOT)+rind;
LIST_R[NB_SPH]=rind;
NB_SPH++;

f.close();

cout<<"Volume des particules :"<<vols<<endl;
}

void read_sphp_car(R EE, R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX)
{
R_SPH=0.;
RMAX=0.;  

H_TOT=-2e9;
V_TOT=-2e9;
Z_TOT=-2e9;

H_POS=0.;
V_POS=0.;
Z_POS=0.;

ifstream f("APPLI/carttest/test",ios::in);
f >> NB_SPH;

int tmp;

R vols=0.;

for (int it=0;it<NB_SPH;it++){
f >> tmp;	
f >> LIST_X[it];
f >> LIST_Y[it];
f >> LIST_Z[it];
f >> LIST_R[it];
LIST_X[it]=LIST_X[it];
LIST_Y[it]=LIST_Y[it];
LIST_Z[it]=LIST_Z[it];
LIST_R[it]=LIST_R[it];
R_SPH+=LIST_R[it];

vols=vols+4./3*3.14159*LIST_R[it]*LIST_R[it]*LIST_R[it];

RMAX=(LIST_R[it]>RMAX)?LIST_R[it]:RMAX;
H_TOT=((LIST_X[it]+LIST_R[it]+EE)>H_TOT)?(LIST_X[it]+LIST_R[it]+EE):H_TOT;
V_TOT=((LIST_Y[it]+LIST_R[it])>V_TOT)?(LIST_Y[it]+LIST_R[it]):V_TOT;
Z_TOT=((LIST_Z[it]+LIST_R[it])>Z_TOT)?(LIST_Z[it]+LIST_R[it]):Z_TOT;
}
R_SPH/=NB_SPH;

f.close();

cout<<"Volume des particules :"<<vols<<endl;
}

void read_sphp(R & H_TOT, R & V_TOT, R & Z_TOT, R & H_POS, R & V_POS, R & Z_POS,int & NB_SPH, R & R_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R & RMAX)
{
R_SPH=0.;
RMAX=0.;  

H_TOT=-2e9;
V_TOT=-2e9;
Z_TOT=-2e9;

H_POS=2e9;
V_POS=2e9;
Z_POS=2e9;

ifstream f("APPLI/carttest/test",ios::in);
f >> NB_SPH;

int tmp;

R vols=0.;

for (int it=0;it<NB_SPH;it++){
f >> tmp;	
f >> LIST_X[it];
f >> LIST_Y[it];
f >> LIST_Z[it];
f >> LIST_R[it];
LIST_X[it]=LIST_X[it];
LIST_Y[it]=LIST_Y[it];
LIST_Z[it]=LIST_Z[it];
LIST_R[it]=LIST_R[it];
R_SPH+=LIST_R[it];

vols=vols+4./3*3.14159*LIST_R[it]*LIST_R[it]*LIST_R[it];

RMAX=(LIST_R[it]>RMAX)?LIST_R[it]:RMAX;
H_TOT=((LIST_X[it]+LIST_R[it])>H_TOT)?(LIST_X[it]+LIST_R[it]):H_TOT;
V_TOT=((LIST_Y[it]+LIST_R[it])>V_TOT)?(LIST_Y[it]+LIST_R[it]):V_TOT;
Z_TOT=((LIST_Z[it]+LIST_R[it])>Z_TOT)?(LIST_Z[it]+LIST_R[it]):Z_TOT;
H_POS=((LIST_X[it]-LIST_R[it])<H_POS)?(LIST_X[it]-LIST_R[it]):H_POS;
V_POS=((LIST_Y[it]-LIST_R[it])<V_POS)?(LIST_Y[it]-LIST_R[it]):V_POS;
Z_POS=((LIST_Z[it]-LIST_R[it])<Z_POS)?(LIST_Z[it]-LIST_R[it]):Z_POS;

}
R_SPH/=NB_SPH;

f.close();

H_TOT=H_TOT-H_POS;
V_TOT=V_TOT-V_POS;
Z_TOT=Z_TOT-Z_POS;

cout<<"POS:"<<H_POS<<", "<<V_POS<<", "<<Z_POS<<endl;

cout<<"Volume des particules :"<<vols<<endl;
}

void read_phase(R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, bool * LIST_P)
{

R rayc=R(V_TOT)/2.6;

R numi=0.;
R ntot=0.;

R coef,dist;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,coef,dist) reduction(+:numi,ntot)
	for(it=0;it<NB_SPH;it++){
        LIST_P[it]=0;
	
	coef=1.;
	
	    if(abs(LIST_X[it]-H_POS)<1e-5){
		coef/=2.;	
		}
        if(abs(LIST_X[it]-(H_TOT+H_POS))<1e-5){
		coef/=2.;	
		}		
        if(abs(LIST_Y[it]-V_POS)<1e-5){
		coef/=2.;	
		}		
        if(abs(LIST_Y[it]-(V_TOT+V_POS))<1e-5){
		coef/=2.;	
		}				
        if(abs(LIST_Z[it]-Z_POS)<1e-5){
		coef/=2.;	
		}		
        if(abs(LIST_Z[it]-(Z_TOT+Z_POS))<1e-5){
		coef/=2.;	
		}	 
	
	dist=sqrt((LIST_Z[it]-R(Z_TOT)/2.-(Z_POS))*(LIST_Z[it]-R(Z_TOT)/2.-(Z_POS))+(LIST_Y[it]-R(V_TOT)/2.-(V_POS))*(LIST_Y[it]-R(V_TOT)/2.-(V_POS)));	

	if(dist<rayc){
        numi=numi+coef;
        LIST_P[it]=1;}

	ntot=ntot+coef;
    } //fin for

R frac=numi/ntot;
cout<<"PR PART :"<<frac<<endl;

	
}

void read_phase_dela(R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,  int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, bool * LIST_P, int * LIST_B)
{

R numi=0.;
R ntot=0.;


R coef,dist;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,coef,dist) reduction(+:numi,ntot)
	for(it=0;it<NB_SPH;it++){
		LIST_P[it]=0;
		LIST_B[it]=0;

		coef=1.;
	
		if(abs(LIST_X[it]-H_POS)<1e-5){
		coef/=2.;	
		}
		if(abs(LIST_X[it]-(H_TOT+H_POS))<1e-5){
		coef/=2.;	
		}		
		if(abs(LIST_Y[it]-V_POS)<1e-5){
		coef/=2.;	
		}		
		if(abs(LIST_Y[it]-(V_TOT+V_POS))<1e-5){
		coef/=2.;	
		}				
		if(abs(LIST_Z[it]-Z_POS)<1e-5){
		coef/=2.;	
		}		
		if(abs(LIST_Z[it]-(Z_TOT+Z_POS))<1e-5){
		coef/=2.;	
		}	 

	     //DÉLAMINAGE
		if(LIST_Z[it]<Z_TOT/2.){
			numi=numi+coef;LIST_P[it]=1;LIST_B[it]=1;}
		ntot=ntot+coef;
	}

R frac=numi/ntot;
cout<<"PR PART :"<<frac<<endl;

	
}

void read_phase_rupt(R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,  int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, bool * LIST_P, int * LIST_B)
{

R rayc=R(V_TOT)/4.;
R numi=0.;
R ntot=0.;

R coef,dist;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,coef,dist) reduction(+:numi,ntot)
	for(it=0;it<NB_SPH;it++){
		LIST_P[it]=0;
		LIST_B[it]=0;
	
		coef=1.;
	
		if(abs(LIST_X[it]-H_POS)<1e-5){
		coef/=2.;	
		}
		if(abs(LIST_X[it]-(H_TOT+H_POS))<1e-5){
		coef/=2.;	
		}		
		if(abs(LIST_Y[it]-V_POS)<1e-5){
		coef/=2.;	
		}		
		if(abs(LIST_Y[it]-(V_TOT+V_POS))<1e-5){
		coef/=2.;	
		}				
		if(abs(LIST_Z[it]-Z_POS)<1e-5){
		coef/=2.;	
		}		
		if(abs(LIST_Z[it]-(Z_TOT+Z_POS))<1e-5){
		coef/=2.;	
		}	 
	
		dist=sqrt((LIST_Z[it]-R(Z_TOT)/2.-(Z_POS))*(LIST_Z[it]-R(Z_TOT)/2.-(Z_POS))+(LIST_X[it]-R(H_TOT)/2.-(H_POS))*(LIST_X[it]-R(H_TOT)/2.-(H_POS)));
		if(dist<rayc){numi=numi+coef;LIST_P[it]=1;LIST_B[it]=1;}
		ntot=ntot+coef;
	}

R frac=numi/ntot;
cout<<"PR PART :"<<frac<<endl;

	
}

void read_phase_compo(R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,  int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, bool * LIST_P)
{

int NB_INC;
R R_INC;

ifstream f("APPLI/compo/carte",ios::in);
f >> NB_INC;
f >> R_INC;
//R_INC/=2.;
R_INC*=H_TOT;
//R_INC*=1.0247;

cout<<"Nb inc:"<<NB_INC<<endl;

R * INC_X = new R[NB_INC];
R * INC_Y = new R[NB_INC];
R * INC_Z = new R[NB_INC];

R PR_INC=0.;

for(int it=0;it<NB_INC;it++){
f >> INC_X[it];
f >> INC_Y[it];
f >> INC_Z[it];

INC_X[it]*=H_TOT;
INC_Y[it]*=V_TOT;
INC_Z[it]*=Z_TOT;

INC_X[it]+=H_POS;
INC_Y[it]+=V_POS;
INC_Z[it]+=Z_POS;

PR_INC+=4./3*3.14159*R_INC*R_INC*R_INC;

}

f.close();

PR_INC/=(H_TOT*V_TOT*Z_TOT);

cout<<"Fraction volumique d'inc : "<<(PR_INC)<<endl;	
	
R dx,dy,dz,dxm,dym,dzm,dxp,dyp,dzp;

R r2=R_INC*R_INC;

int numi=0;	
	
for(int it=0;it<NB_SPH;it++){		
	LIST_P[it]=0;
  if(it%100000==0) cout<<"Particule :"<<it+1<<"/"<<NB_SPH<<endl;

	
  for(int jt=0;jt<NB_INC;jt++){	

  dx=LIST_X[it]-INC_X[jt];
  dxp=dx+H_TOT;
  dxm=dx-H_TOT;
  dx=dx*dx;
  dxm=dxm*dxm;
  dxp=dxp*dxp;
  dy=LIST_Y[it]-INC_Y[jt];    
  dyp=dy+V_TOT;
  dym=dy-V_TOT;    
  dy=dy*dy;
  dym=dym*dym;
  dyp=dyp*dyp;
  dz=LIST_Z[it]-INC_Z[jt];
  dzp=dz+Z_TOT;
  dzm=dz-Z_TOT;     
  dz=dz*dz;
  dzm=dzm*dzm;
  dzp=dzp*dzp;
    
  if(dx+dy+dz<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dy+dz<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dy+dz<=r2) {LIST_P[it]=1;numi++;}  
  else if(dx+dym+dz<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dym+dz<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dym+dz<=r2) {LIST_P[it]=1;numi++;}  
  else if(dx+dyp+dz<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dyp+dz<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dyp+dz<=r2) {LIST_P[it]=1;numi++;}  
  else if(dx+dy+dzm<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dy+dzm<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dy+dzm<=r2) {LIST_P[it]=1;numi++;}  
  else if(dx+dym+dzm<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dym+dzm<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dym+dzm<=r2) {LIST_P[it]=1;numi++;}  
  else if(dx+dyp+dzm<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dyp+dzm<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dyp+dzm<=r2) {LIST_P[it]=1;numi++;}    
  else if(dx+dy+dzp<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dy+dzp<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dy+dzp<=r2) {LIST_P[it]=1;numi++;}  
  else if(dx+dym+dzp<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dym+dzp<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dym+dzp<=r2) {LIST_P[it]=1;numi++;}  
  else if(dx+dyp+dzp<=r2) {LIST_P[it]=1;numi++;}
  else if(dxm+dyp+dzp<=r2) {LIST_P[it]=1;numi++;}  
  else if(dxp+dyp+dzp<=r2) {LIST_P[it]=1;numi++;}    	
  	
  	
  }
}  	

cout<<"PR PART :"<<(R(numi)/NB_SPH)<<endl;

	
}


