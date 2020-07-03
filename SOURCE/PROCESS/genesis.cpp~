#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <time.h> 
#include <sys/time.h> 
#include <list>

#include "genesis.h"

using namespace std;

bool bincl(R R_CYL,R AA,R BB,R CC, R X1, R Y1, R Z1, R X2, R Y2, R Z2, R XX, R YY, R ZZ){

bool boolinc=0;

R EQUAC1 = (XX-X1)*(XX-X1)+(YY-Y1)*(YY-Y1)+(ZZ-Z1)*(ZZ-Z1);
R EQUAC2 = (AA*(XX-X1)+BB*(YY-Y1)+CC*(ZZ-Z1))*(AA*(XX-X1)+BB*(YY-Y1)+CC*(ZZ-Z1));
R EQUAC  = EQUAC1-EQUAC2;

R EQUAP1 = -(AA*(XX-X1)+BB*(YY-Y1)+CC*(ZZ-Z1));
R EQUAP2 = (AA*(XX-X2)+BB*(YY-Y2)+CC*(ZZ-Z2));

if((EQUAC<=(R_CYL*R_CYL)+1e-8)&&(EQUAP1<=1e-8)&&(EQUAP2<=1e-8)){
boolinc=1;  
}
 
return boolinc;
}


void gener3_cyl_alea(int nC,R facf,R eche, R buffer, int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P)  {

buffer=buffer+100./(2.*eche);

if((100.-2.*buffer)<0){
cout<<"Fibres trop longues !!!"<<endl;
exit(0);
}

R lC=100./eche;  
R R_CYL=lC/(2.*facf);
 
R Pi=3.14159265;
int it;    
for(it=0;it<NB_SPH;it++){
LIST_P[it]=0;
}
   
// PARAMETRES ALEATOIRES
struct timeval tv ;
gettimeofday(&tv, NULL) ;
srand(tv.tv_usec) ;

it=0;
while(it<nC){
    
R X_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;
R Y_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;
R Z_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;

	R PH_REF  = rand()%360;    
	R GA_REF  = ((R) (rand()%100))/100.;
	GA_REF    = acos(GA_REF);
	R PS_REF  = rand()%360;  

	PH_REF    = PH_REF*Pi/180; 
	PS_REF    = PS_REF*Pi/180; 

	R AA = cos(PH_REF)*cos(PS_REF)-cos(GA_REF)*sin(PS_REF)*sin(PH_REF);
	R BB = sin(PH_REF)*cos(PS_REF)+cos(GA_REF)*sin(PS_REF)*cos(PH_REF);
	R CC = sin(GA_REF)*sin(PS_REF);

	R XX1=X_GRA-0.5*lC*AA;
	R YY1=Y_GRA-0.5*lC*BB;
	R ZZ1=Z_GRA-0.5*lC*CC;
	R XX2=X_GRA+0.5*lC*AA;
	R YY2=Y_GRA+0.5*lC*BB;
	R ZZ2=Z_GRA+0.5*lC*CC;

	bool boolc=0;
	list<int> list_tmp; 

        int jt=0;
	while((jt<NB_SPH)&&(boolc==0)){

		bool boolxa = bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt],LIST_Z[jt]);

		if((boolxa==1)&&(LIST_P[jt]==1))
		{
		boolc=1;  
		}
		else if((boolxa==1)&&(LIST_P[jt]==0)){
		list_tmp.push_back(jt);
		}

	jt++;    
	}

 //cout<<"it:"<<it<<", "<<boolc<<endl;

        if(boolc==0){
	       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
		 LIST_P[*kt]=1;
	       }	
	}
        else{
        it--;
        }
it++;
}

R phi_inc=0.;
for(it=0;it<NB_SPH;it++){
if(LIST_P[it]) phi_inc=phi_inc+1.;
}
phi_inc/=NB_SPH;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*R_CYL*R_CYL*lC)/1e6<<endl; 
cout<<"Nombre de particules par fibre :"<< int(phi_inc*NB_SPH/nC) <<endl; 
}

void gener3_cyl_alea_per(int nC,R facf,R eche,  int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P)  {

R lC=100./eche;  
R R_CYL=lC/(2.*facf);
 
R Pi=3.14159265;
int it;    
for(it=0;it<NB_SPH;it++){
LIST_P[it]=0;
}
   
// PARAMETRES ALEATOIRES
struct timeval tv ;
gettimeofday(&tv, NULL) ;
srand(tv.tv_usec) ;

it=0;
while(it<nC){
    
R X_GRA = ( rand()/(R)RAND_MAX ) * 100.;
R Y_GRA = ( rand()/(R)RAND_MAX ) * 100.;
R Z_GRA = ( rand()/(R)RAND_MAX ) * 100.;

	R PH_REF  = rand()%360;    
	R GA_REF  = ((R) (rand()%100))/100.;
	GA_REF    = acos(GA_REF);
	R PS_REF  = rand()%360;  

	PH_REF    = PH_REF*Pi/180; 
	PS_REF    = PS_REF*Pi/180; 

	R AA = cos(PH_REF)*cos(PS_REF)-cos(GA_REF)*sin(PS_REF)*sin(PH_REF);
	R BB = sin(PH_REF)*cos(PS_REF)+cos(GA_REF)*sin(PS_REF)*cos(PH_REF);
	R CC = sin(GA_REF)*sin(PS_REF);

	R XX1=X_GRA-0.5*lC*AA;
	R YY1=Y_GRA-0.5*lC*BB;
	R ZZ1=Z_GRA-0.5*lC*CC;
	R XX2=X_GRA+0.5*lC*AA;
	R YY2=Y_GRA+0.5*lC*BB;
	R ZZ2=Z_GRA+0.5*lC*CC;

	bool boolc=0;
	list<int> list_tmp; 

        int jt=0;
	while((jt<NB_SPH)&&(boolc==0)){

		bool boolxa = bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]-100.,LIST_Z[jt]-100.);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]-100.,LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]-100.,LIST_Z[jt]-100.),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt],LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt],LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt],LIST_Z[jt]-100.),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]+100.,LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]+100.,LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]+100.,LIST_Z[jt]-100.),boolxa);


		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]-100.,LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]-100.,LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]-100.,LIST_Z[jt]),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt],LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt],LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt],LIST_Z[jt]),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]+100.,LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]+100.,LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]+100.,LIST_Z[jt]),boolxa);


		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]-100.,LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]-100.,LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]-100.,LIST_Z[jt]+100.),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt],LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt],LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt],LIST_Z[jt]+100.),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]+100.,LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]+100.,LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]+100.,LIST_Z[jt]+100.),boolxa);


		if((boolxa==1)&&(LIST_P[jt]==1))
		{
		boolc=1;  
		}
		else if((boolxa==1)&&(LIST_P[jt]==0)){
		list_tmp.push_back(jt);
		}

	jt++;    
	}

 //cout<<"it:"<<it<<", "<<boolc<<endl;

        if(boolc==0){
	       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
		 LIST_P[*kt]=1;
	       }	
	}
        else{
        it--;
        }
it++;
}

R phi_inc=0.;
for(it=0;it<NB_SPH;it++){
if(LIST_P[it]) phi_inc=phi_inc+1.;
}
phi_inc/=NB_SPH;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*R_CYL*R_CYL*lC)/1e6<<endl; 
cout<<"Nombre de particules par fibre :"<< int(phi_inc*NB_SPH/nC) <<endl; 
}

void gener3_cyl_align(int nC,R facf,R eche, R buffer, int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P,int dir)  {

R bufferx,buffery,bufferz;
if(dir==1){
bufferx=buffer+100./(2.*eche);
buffery=buffer+100./(2.*eche*facf);
bufferz=buffer+100./(2.*eche*facf);
}else if(dir==2){
bufferx=buffer+100./(2.*eche*facf);
buffery=buffer+100./(2.*eche);
bufferz=buffer+100./(2.*eche*facf);
}else if(dir==3){
bufferx=buffer+100./(2.*eche*facf);
buffery=buffer+100./(2.*eche*facf);
bufferz=buffer+100./(2.*eche);
}else{
bufferx=buffer+100./(2.*eche);
buffery=buffer+100./(2.*eche*facf);
bufferz=buffer+100./(2.*eche*facf);
}

if((100.-2.*bufferx)<0||(100.-2.*buffery)<0||(100.-2.*bufferz)<0){
cout<<"Fibres trop longues !!!"<<endl;
exit(0);
}


R lC=100./eche;  
R R_CYL=lC/(2.*facf);
 
R Pi=3.14159265;
int it;    
for(it=0;it<NB_SPH;it++){
LIST_P[it]=0;
}
   
// PARAMETRES ALEATOIRES
struct timeval tv ;
gettimeofday(&tv, NULL) ;
srand(tv.tv_usec) ;

it=0;
while(it<nC){
    
R X_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*bufferx) + bufferx;
R Y_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffery) + buffery;
R Z_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*bufferz) + bufferz;

R XX1,YY1,ZZ1,XX2,YY2,ZZ2;
R AA,BB,CC;

if(dir==1){
XX1=X_GRA-0.5*lC;
YY1=Y_GRA;
ZZ1=Z_GRA;
XX2=X_GRA+0.5*lC;
YY2=Y_GRA;
ZZ2=Z_GRA;
AA=1.;
BB=0.;
CC=0.;
}else if(dir==2){
XX1=X_GRA;
YY1=Y_GRA-0.5*lC;
ZZ1=Z_GRA;
XX2=X_GRA;
YY2=Y_GRA+0.5*lC;
ZZ2=Z_GRA;
AA=0.;
BB=1.;
CC=0.;
}else if(dir==3){
XX1=X_GRA;
YY1=Y_GRA;
ZZ1=Z_GRA-0.5*lC;
XX2=X_GRA;
YY2=Y_GRA;
ZZ2=Z_GRA+0.5*lC;
AA=0.;
BB=0.;
CC=1.;
}else{
cout<<"Mauvaise direction !! Dir. x par défaut"<<endl;
XX1=X_GRA-0.5*lC;
YY1=Y_GRA;
ZZ1=Z_GRA;
XX2=X_GRA+0.5*lC;
YY2=Y_GRA;
ZZ2=Z_GRA;
AA=1.;
BB=0.;
CC=0.;
}

	bool boolc=0;
	list<int> list_tmp; 

        int jt=0;
	while((jt<NB_SPH)&&(boolc==0)){

		bool boolxa = bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt],LIST_Z[jt]);

		if((boolxa==1)&&(LIST_P[jt]==1))
		{
		boolc=1;  
		}
		else if((boolxa==1)&&(LIST_P[jt]==0)){
		list_tmp.push_back(jt);
		}

	jt++;    
	}

 //cout<<"it:"<<it<<", "<<boolc<<endl;

        if(boolc==0){
	       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
		 LIST_P[*kt]=1;
	       }	
	}
        else{
        it--;
        }
it++;
}

R phi_inc=0.;
for(it=0;it<NB_SPH;it++){
if(LIST_P[it]) phi_inc=phi_inc+1.;
}
phi_inc/=NB_SPH;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*R_CYL*R_CYL*lC)/1e6<<endl; 
cout<<"Nombre de particules par fibre :"<< int(phi_inc*NB_SPH/nC) <<endl; 
}

void gener3_cyl_align_per(int nC,R facf,R eche,  int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P, int dir)  {

R lC=100./eche;  
R R_CYL=lC/(2.*facf);
 
R Pi=3.14159265;
int it;    
for(it=0;it<NB_SPH;it++){
LIST_P[it]=0;
}
   
// PARAMETRES ALEATOIRES
struct timeval tv ;
gettimeofday(&tv, NULL) ;
srand(tv.tv_usec) ;

it=0;
while(it<nC){
    
R X_GRA = ( rand()/(R)RAND_MAX ) * 100.;
R Y_GRA = ( rand()/(R)RAND_MAX ) * 100.;
R Z_GRA = ( rand()/(R)RAND_MAX ) * 100.;

R XX1,YY1,ZZ1,XX2,YY2,ZZ2;
R AA,BB,CC;

if(dir==1){
XX1=X_GRA-0.5*lC;
YY1=Y_GRA;
ZZ1=Z_GRA;
XX2=X_GRA+0.5*lC;
YY2=Y_GRA;
ZZ2=Z_GRA;
AA=1.;
BB=0.;
CC=0.;
}else if(dir==2){
XX1=X_GRA;
YY1=Y_GRA-0.5*lC;
ZZ1=Z_GRA;
XX2=X_GRA;
YY2=Y_GRA+0.5*lC;
ZZ2=Z_GRA;
AA=0.;
BB=1.;
CC=0.;
}else if(dir==3){
XX1=X_GRA;
YY1=Y_GRA;
ZZ1=Z_GRA-0.5*lC;
XX2=X_GRA;
YY2=Y_GRA;
ZZ2=Z_GRA+0.5*lC;
AA=0.;
BB=0.;
CC=1.;
}else{
cout<<"Mauvaise direction !! Dir. x par défaut"<<endl;
XX1=X_GRA-0.5*lC;
YY1=Y_GRA;
ZZ1=Z_GRA;
XX2=X_GRA+0.5*lC;
YY2=Y_GRA;
ZZ2=Z_GRA;
AA=1.;
BB=0.;
CC=0.;
}

	bool boolc=0;
	list<int> list_tmp; 

        int jt=0;
	while((jt<NB_SPH)&&(boolc==0)){

		bool boolxa = bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]-100.,LIST_Z[jt]-100.);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]-100.,LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]-100.,LIST_Z[jt]-100.),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt],LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt],LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt],LIST_Z[jt]-100.),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]+100.,LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]+100.,LIST_Z[jt]-100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]+100.,LIST_Z[jt]-100.),boolxa);


		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]-100.,LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]-100.,LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]-100.,LIST_Z[jt]),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt],LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt],LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt],LIST_Z[jt]),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]+100.,LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]+100.,LIST_Z[jt]),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]+100.,LIST_Z[jt]),boolxa);


		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]-100.,LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]-100.,LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]-100.,LIST_Z[jt]+100.),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt],LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt],LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt],LIST_Z[jt]+100.),boolxa);

		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]-100.,LIST_Y[jt]+100.,LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt],LIST_Y[jt]+100.,LIST_Z[jt]+100.),boolxa);
		boolxa = max(bincl(R_CYL,AA,BB,CC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,LIST_X[jt]+100.,LIST_Y[jt]+100.,LIST_Z[jt]+100.),boolxa);


		if((boolxa==1)&&(LIST_P[jt]==1))
		{
		boolc=1;  
		}
		else if((boolxa==1)&&(LIST_P[jt]==0)){
		list_tmp.push_back(jt);
		}

	jt++;    
	}

 //cout<<"it:"<<it<<", "<<boolc<<endl;

        if(boolc==0){
	       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
		 LIST_P[*kt]=1;
	       }	
	}
        else{
        it--;
        }
it++;
}

R phi_inc=0.;
for(it=0;it<NB_SPH;it++){
if(LIST_P[it]) phi_inc=phi_inc+1.;
}
phi_inc/=NB_SPH;

cout<<"Fraction volumique de fibres :"<< phi_inc <<" / "<<(nC*Pi*R_CYL*R_CYL*lC)/1e6<<endl; 
cout<<"Nombre de particules par fibre :"<< int(phi_inc*NB_SPH/nC) <<endl; 
}

void gener3_sphere_alea(int nS,R eche, R buffer, R tol, int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P)  {

R Pi=3.14159265;
R rS=100./(2.*eche);
buffer=buffer+rS;

if((100.-2.*buffer)<0){
cout<<"Spheres trop grandes !!!"<<endl;
exit(0);
}

int it;    
for(it=0;it<NB_SPH;it++){
LIST_P[it]=0;
}
   
// PARAMETRES ALEATOIRES
struct timeval tv ;
gettimeofday(&tv, NULL) ;
srand(tv.tv_usec) ;

R XA[nS],YA[nS],ZA[nS];

it=0;
while(it<nS){
    
R X_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;
R Y_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;
R Z_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;

XA[it] = X_GRA;
YA[it] = Y_GRA;
ZA[it] = Z_GRA;

	bool boolc=0;
	list<int> list_tmp; 

        int jt=0;
	R distij;
	for(jt=0;jt<it;jt++){
	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;
	}

        jt=0;
	R dist;
	while((jt<NB_SPH)&&(boolc==0)){
                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);

		bool boolxa = (dist<rS);

		if((boolxa==1)&&(LIST_P[jt]==1))
		{
		boolc=1;  
		}
		else if((boolxa==1)&&(LIST_P[jt]==0)){
		list_tmp.push_back(jt);
		}

	jt++;    
	}

 //cout<<"it:"<<it<<", "<<boolc<<endl;

        if(boolc==0){
	       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
		 LIST_P[*kt]=1;
	       }	
	}
        else{
        it--;
        }
it++;
}

R phi_inc=0.;
for(it=0;it<NB_SPH;it++){
if(LIST_P[it]) phi_inc=phi_inc+1.;
}
phi_inc/=NB_SPH;

cout<<"Fraction volumique de spheres :"<< phi_inc <<" / "<<(nS*4./3*Pi*rS*rS*rS)/1e6<<endl; 
cout<<"Nombre de particules par sphere :"<< int(phi_inc*NB_SPH/nS) <<endl; 
}

void gener3_sphere_alea_per(int nS,R eche, R tol, int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P)  {

R Pi=3.14159265;
R rS=100./(2.*eche);

int it;    
for(it=0;it<NB_SPH;it++){
LIST_P[it]=0;
}
   
// PARAMETRES ALEATOIRES
struct timeval tv ;
gettimeofday(&tv, NULL) ;
srand(tv.tv_usec) ;

R XA[nS],YA[nS],ZA[nS];

it=0;
while(it<nS){
    
R X_GRA = ( rand()/(R)RAND_MAX ) * 100.;
R Y_GRA = ( rand()/(R)RAND_MAX ) * 100.;
R Z_GRA = ( rand()/(R)RAND_MAX ) * 100.;

XA[it] = X_GRA;
YA[it] = Y_GRA;
ZA[it] = Z_GRA;

	bool boolc=0;
	list<int> list_tmp; 

        int jt=0;
	R distij;
	for(jt=0;jt<it;jt++){

	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]-100.-ZA[jt])*(ZA[it]-100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;




	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;




	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]-100.-YA[jt])*(YA[it]-100.-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-100.-XA[jt])*(XA[it]-100.-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;

	distij=(XA[it]+100.-XA[jt])*(XA[it]+100.-XA[jt])+(YA[it]+100.-YA[jt])*(YA[it]+100.-YA[jt])+(ZA[it]+100.-ZA[jt])*(ZA[it]+100.-ZA[jt]);
        distij=sqrt(distij);
        if(distij<2.*rS*(1.+tol/100.)) boolc=1;



	}

        jt=0;
	R dist;
	while((jt<NB_SPH)&&(boolc==0)){

		bool boolxa=0;

                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]-100.-LIST_Z[jt])*(ZA[it]-100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);




                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);




                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]-100.-LIST_Y[jt])*(YA[it]-100.-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-100.-LIST_X[jt])*(XA[it]-100.-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);

                dist=(XA[it]+100.-LIST_X[jt])*(XA[it]+100.-LIST_X[jt])+(YA[it]+100.-LIST_Y[jt])*(YA[it]+100.-LIST_Y[jt])+(ZA[it]+100.-LIST_Z[jt])*(ZA[it]+100.-LIST_Z[jt]);
                dist=sqrt(dist);
                boolxa=max((dist<rS),boolxa);


		if((boolxa==1)&&(LIST_P[jt]==1))
		{
		boolc=1;  
		}
		else if((boolxa==1)&&(LIST_P[jt]==0)){
		list_tmp.push_back(jt);
		}

	jt++;    
	}

 //cout<<"it:"<<it<<", "<<boolc<<endl;

        if(boolc==0){
	       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
		 LIST_P[*kt]=1;
	       }	
	}
        else{
        it--;
        }
it++;
}

R phi_inc=0.;
for(it=0;it<NB_SPH;it++){
if(LIST_P[it]) phi_inc=phi_inc+1.;
}
phi_inc/=NB_SPH;

cout<<"Fraction volumique de spheres :"<< phi_inc <<" / "<<(nS*4./3*Pi*rS*rS*rS)/1e6<<endl; 
cout<<"Nombre de particules par sphere :"<< int(phi_inc*NB_SPH/nS) <<endl; 
}

void gener3_cyl_inf(int nC,R eche, R buffer, R tol, int NB_SPH, R * LIST_X, R * LIST_Y, R *LIST_Z, bool * LIST_P,int dir)  {

R Pi=3.14159265;
R rC=100./(2.*eche);
buffer=buffer+rC;

int it;    
for(it=0;it<NB_SPH;it++){
LIST_P[it]=0;
}
   
// PARAMETRES ALEATOIRES
struct timeval tv ;
gettimeofday(&tv, NULL) ;
srand(tv.tv_usec) ;

R XA[nC],YA[nC],ZA[nC];


if(dir==1){
	it=0;
	while(it<nC){
	    
	R Y_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;
	R Z_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;

	YA[it] = Y_GRA;
	ZA[it] = Z_GRA;

		bool boolc=0;
		list<int> list_tmp; 

		int jt=0;
		R distij;
		for(jt=0;jt<it;jt++){
		distij=(YA[it]-YA[jt])*(YA[it]-YA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
		distij=sqrt(distij);
		if(distij<2.*rC*(1.+tol/100.)) boolc=1;
		}

		jt=0;
		R dist;
		while((jt<NB_SPH)&&(boolc==0)){
		        dist=(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
		        dist=sqrt(dist);

			bool boolxa = (dist<rC);

			if((boolxa==1)&&(LIST_P[jt]==1))
			{
			boolc=1;  
			}
			else if((boolxa==1)&&(LIST_P[jt]==0)){
			list_tmp.push_back(jt);
			}

		jt++;    
		}

	 //cout<<"it:"<<it<<", "<<boolc<<endl;

		if(boolc==0){
		       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
			 LIST_P[*kt]=1;
		       }	
		}
		else{
		it--;
		}
	it++;
	}

}else if(dir==2){
	it=0;
	while(it<nC){
	    
	R X_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;
	R Z_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;

	XA[it] = X_GRA;
	ZA[it] = Z_GRA;

		bool boolc=0;
		list<int> list_tmp; 

		int jt=0;
		R distij;
		for(jt=0;jt<it;jt++){
		distij=(XA[it]-XA[jt])*(XA[it]-XA[jt])+(ZA[it]-ZA[jt])*(ZA[it]-ZA[jt]);
		distij=sqrt(distij);
		if(distij<2.*rC*(1.+tol/100.)) boolc=1;
		}

		jt=0;
		R dist;
		while((jt<NB_SPH)&&(boolc==0)){
		        dist=(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt])+(ZA[it]-LIST_Z[jt])*(ZA[it]-LIST_Z[jt]);
		        dist=sqrt(dist);

			bool boolxa = (dist<rC);

			if((boolxa==1)&&(LIST_P[jt]==1))
			{
			boolc=1;  
			}
			else if((boolxa==1)&&(LIST_P[jt]==0)){
			list_tmp.push_back(jt);
			}

		jt++;    
		}

	 //cout<<"it:"<<it<<", "<<boolc<<endl;

		if(boolc==0){
		       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
			 LIST_P[*kt]=1;
		       }	
		}
		else{
		it--;
		}
	it++;
	}

}else if(dir==3){
	it=0;
	while(it<nC){
	    
	R Y_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;
	R X_GRA = ( rand()/(R)RAND_MAX ) * (100.-2.*buffer) + buffer;

	YA[it] = Y_GRA;
	XA[it] = X_GRA;

		bool boolc=0;
		list<int> list_tmp; 

		int jt=0;
		R distij;
		for(jt=0;jt<it;jt++){
		distij=(YA[it]-YA[jt])*(YA[it]-YA[jt])+(XA[it]-XA[jt])*(XA[it]-XA[jt]);
		distij=sqrt(distij);
		if(distij<2.*rC*(1.+tol/100.)) boolc=1;
		}

		jt=0;
		R dist;
		while((jt<NB_SPH)&&(boolc==0)){
		        dist=(YA[it]-LIST_Y[jt])*(YA[it]-LIST_Y[jt])+(XA[it]-LIST_X[jt])*(XA[it]-LIST_X[jt]);
		        dist=sqrt(dist);

			bool boolxa = (dist<rC);

			if((boolxa==1)&&(LIST_P[jt]==1))
			{
			boolc=1;  
			}
			else if((boolxa==1)&&(LIST_P[jt]==0)){
			list_tmp.push_back(jt);
			}

		jt++;    
		}

	 //cout<<"it:"<<it<<", "<<boolc<<endl;

		if(boolc==0){
		       for (std::list<int>::iterator kt=list_tmp.begin(); kt!=list_tmp.end(); ++kt){	 
			 LIST_P[*kt]=1;
		       }	
		}
		else{
		it--;
		}
	it++;
	}

}else{
cout<<"Mauvaise direction !!!"<<endl; 
}




R phi_inc=0.;
for(it=0;it<NB_SPH;it++){
if(LIST_P[it]) phi_inc=phi_inc+1.;
}
phi_inc/=NB_SPH;

cout<<"Fraction volumique de spheres :"<< phi_inc <<" / "<<(nC*Pi*rC*rC)/1e4<<endl; 
cout<<"Nombre de particules par sphere :"<< int(phi_inc*NB_SPH/nC) <<endl; 
}



