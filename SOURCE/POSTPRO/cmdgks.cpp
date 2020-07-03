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

#include "cmdgks.h"

void ExpDEF(int NB_SPH, R * DEF, bool * EDGE){
   
int numc=1;

	fstream g;
	g.precision(4);
	g.open("DATA/datadef.txt",fstream::out);
	g<<setw(12)<<"NUM"<<setw(12)<<"DEF"<<endl;
	for(int it=0; it<NB_SPH; it++){
		if(!EDGE[it]){
		g<<setw(12)<<numc<<setw(12)<< DEF[it] <<endl;
		numc++;
		}
	}
	g.close();

}

void ExpSIG(int NB_SPH, R * SIG, bool * EDGE){
   
int numc=1;

	fstream g;
	g.precision(4);
	g.open("DATA/datasig.txt",fstream::out);
	g<<setw(12)<<"NUM"<<setw(12)<<"SIG"<<endl;
	for(int it=0; it<NB_SPH; it++){
		if(!EDGE[it]){
		g<<setw(12)<<numc<<setw(12)<< SIG[it] <<endl;
		numc++;
		}
	}
	g.close();

}

void ExpDEPL(R dplct, int ite,bool bstart){
  
fstream g;
g.precision(4);
if(bstart){
g.open("DATA/datadepl.txt",fstream::out );	
g<<setw(12)<<"ITE"<<setw(12)<<"U"<<endl;
}else{
g.open("DATA/datadepl.txt",fstream::out | fstream::app);	
}
g<<setw(12)<<ite<<setw(12)<<dplct<<endl;
g.close();
}

void ExpREAC(R reac, int ite,bool bstart){
  
fstream g;
g.precision(4);
if(bstart){
g.open("DATA/datareac.txt",fstream::out );	
g<<setw(12)<<"ITE"<<setw(12)<<"F"<<endl;
}else{
g.open("DATA/datareac.txt",fstream::out | fstream::app);	
}
g<<setw(12)<<ite<<setw(12)<<reac<<endl;
g.close();
}

void ExpHalo(int NB_SPH, int * NBHALO, R * VOLHALO, int ite){
  
ostringstream oss;
oss << ite;


string filename0 = "DATA/datahalo" + oss.str() + ".txt";
const char * filename1 = filename0.c_str();
  
fstream g;
g.precision(4);
g.open(filename1,fstream::out );	

g<<setw(12)<<"NUM"<<setw(12)<<"NBH"<<setw(12)<<"VOLH"<<endl;
for(int jte=0;jte<NB_SPH;jte++){
g<<setw(12)<<jte<<setw(12)<<NBHALO[jte]<<setw(12)<<VOLHALO[jte]<<endl;
}

g.close();

}

void ExpENERGY(R Ec, R Ep,R Epa, int ite,bool bstart){
  
fstream g;
g.precision(4);
if(bstart){
g.open("DATA/dataenergy.txt",fstream::out );	
g<<setw(12)<<"ITE"<<setw(12)<<"EC"<<setw(12)<<"EP"<<setw(12)<<"EPA"<<setw(12)<<"EM"<<endl;
}else{
g.open("DATA/dataenergy.txt",fstream::out | fstream::app);	
}
g<<setw(12)<<ite<<setw(12)<<Ec<<setw(12)<<Ep<<setw(12)<<Epa<<setw(12)<<(Ep+Epa+Ec)<<endl;
g.close();
}

void ExpDILA(int ite,R dpltx,R dplty,R dpltz,bool bstart){
  
fstream g;
g.precision(4);
if(bstart){
g.open("DATA/datadila.txt",fstream::out);
g<<setw(12)<<"ITE"<<setw(12)<<"DEPX"<<setw(12)<<"DEPY"<<setw(12)<<"DEPZ"<<endl;
}else{
g.open("DATA/datadila.txt",fstream::out | fstream::app);
}
g<<setw(12)<<ite<<setw(12)<< dpltx <<setw(12)<< dplty <<setw(12)<< dpltz <<endl;
g.close();
}

void ExpPARAMETER(R Emoy, R nuu, int NB_SPH, int NBCO, R cpu, int ite, R dt, R coord1,bool bstart){
  
fstream g;
g.precision(4);

if(bstart){
g.open("DATA/datacoh.txt",fstream::out);
g<<setw(12)<<"E"<<setw(12)<<"NU"<<setw(12)<<"NB_SPH"<<setw(12)<<"NBCO"<<setw(12)<<"CPU"<<setw(12)<<"ITE"<<setw(12)<<"TPS"<<setw(12)<<"COORD"<<endl;
}else{
g.open("DATA/datacoh.txt",fstream::out | fstream::app);
}

g<<setw(12)<<Emoy<<setw(12)<<nuu<<setw(12)<<NB_SPH<<setw(12)<<NBCO<<setw(12)<<cpu<<setw(12)<<ite<<setw(12)<<dt*ite<<setw(12)<<coord1<<endl;
g.close();
}

void ExpPARAMETER_rupt(int nrt, int nrc, int nrcis, int nrtot, int ite, bool bstart){
fstream g;
g.precision(4);

if(bstart){
g.open("DATA/datarupt.txt",fstream::out );
g<<setw(12)<<"ITE"<<setw(12)<<"NRT"<<setw(12)<<"NRC"<<setw(12)<<"NRCIS"<<setw(12)<<"NRTOT"<<endl;
}else{
g.open("DATA/datarupt.txt",fstream::out | fstream::app);
}

g<<setw(12)<<ite<<setw(12)<<nrt<<setw(12)<<nrc<<setw(12)<<nrcis<<setw(12)<<nrtot<<endl;
g.close();
}

void ExpPARAMETER_ind(int ite, R ti, R force, int nbci, R pos, bool bstart){
  
fstream g;
g.precision(4);

if(bstart){
g.open("DATA/dataind.txt",fstream::out);
g<<setw(12)<<"ITE"<<setw(12)<<"TPS"<<setw(12)<<"FORCE"<<setw(12)<<"NBCI"<<setw(12)<<"POS"<<endl;
}else{
g.open("DATA/dataind.txt",fstream::out | fstream::app);
}

g<<setw(12)<<ite<<setw(12)<<ti<<setw(12)<<force<<setw(12)<<nbci<<setw(12)<<pos<<endl;
g.close();
}

void ExpVISU(int ite, R dt, int NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R,R * VONMIS,R * TRACE,R * SIG11,R * SIG12,R * SIG13,R * SIG22,R * SIG23,R * SIG33,R * SIG1,R * SIG2,R * SIG3,R * DEP1,R * DEP2,R * DEP3,R maxvm, R minvm, R maxtrac, R mintrac,R maxsig11, R minsig11,R maxsig12, R minsig12,R maxsig13, R minsig13,R maxsig22, R minsig22,R maxsig23, R minsig23,R maxsig33, R minsig33,R maxsig1, R minsig1,R maxsig2, R minsig2,R maxsig3, R minsig3,R maxdep1, R mindep1,R maxdep2, R mindep2,R maxdep3, R mindep3, bool bout, bool * EDGE, bool * LIST_P){
  
  // fichier visu
fstream f1;
f1.precision(4);
if(!bout){
f1.open("VISU/vdepl",fstream::out | fstream::app);
}
else{
f1.open("VISU/vdepl",fstream::out); 
f1<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;
}

f1<<setw(12)<<NB_SPH<<setw(12)<<NB_PAR<<setw(12)<<scientific<<(dt*ite)<<endl;

for(int it=0;it<NB_PAR;it++){
f1<<setw(12)<<scientific<<LIST_PX[it][0]<<setw(12)<<scientific<<LIST_PY[it][0]<<setw(12)<<scientific<<LIST_PZ[it][0]<<endl;
f1<<setw(12)<<scientific<<LIST_PX[it][1]<<setw(12)<<scientific<<LIST_PY[it][1]<<setw(12)<<scientific<<LIST_PZ[it][1]<<endl;
f1<<setw(12)<<scientific<<LIST_PX[it][2]<<setw(12)<<scientific<<LIST_PY[it][2]<<setw(12)<<scientific<<LIST_PZ[it][2]<<endl;
f1<<setw(12)<<scientific<<LIST_PX[it][3]<<setw(12)<<scientific<<LIST_PY[it][3]<<setw(12)<<scientific<<LIST_PZ[it][3]<<endl;
}

f1<<setw(12)<<scientific<<maxdep1<<setw(12)<<scientific<<mindep1<<endl;
f1<<setw(12)<<scientific<<maxdep2<<setw(12)<<scientific<<mindep2<<endl;
f1<<setw(12)<<scientific<<maxdep3<<setw(12)<<scientific<<mindep3<<endl;

for(int it=0;it<NB_SPH;it++){
f1<<setw(12)<<scientific<<EDGE[it]<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<setw(12)<<scientific<<DEP1[it]<<setw(12)<<scientific<<DEP2[it]<<setw(12)<<scientific<<DEP3[it]<<endl; 
}

f1.close();

fstream f2;
f2.precision(4);
if(!bout){
f2.open("VISU/vcontr",fstream::out | fstream::app);
}
else{
f2.open("VISU/vcontr",fstream::out); 
f2<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;
}

f2<<setw(12)<<NB_SPH<<setw(12)<<0<<setw(12)<<scientific<<(dt*ite)<<endl;

f2<<setw(12)<<scientific<<maxvm<<setw(12)<<scientific<<minvm<<endl;
f2<<setw(12)<<scientific<<maxtrac<<setw(12)<<scientific<<mintrac<<endl;
f2<<setw(12)<<scientific<<maxsig11<<setw(12)<<scientific<<minsig11<<endl;
f2<<setw(12)<<scientific<<maxsig12<<setw(12)<<scientific<<minsig12<<endl;
f2<<setw(12)<<scientific<<maxsig13<<setw(12)<<scientific<<minsig13<<endl;
f2<<setw(12)<<scientific<<maxsig22<<setw(12)<<scientific<<minsig22<<endl;
f2<<setw(12)<<scientific<<maxsig23<<setw(12)<<scientific<<minsig23<<endl;
f2<<setw(12)<<scientific<<maxsig33<<setw(12)<<scientific<<minsig33<<endl;
f2<<setw(12)<<scientific<<maxsig1<<setw(12)<<scientific<<minsig1<<endl;
f2<<setw(12)<<scientific<<maxsig2<<setw(12)<<scientific<<minsig2<<endl;
f2<<setw(12)<<scientific<<maxsig3<<setw(12)<<scientific<<minsig3<<endl;

for(int it=0;it<NB_SPH;it++){
f2<<setw(12)<<scientific<<EDGE[it]<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<setw(12)<<scientific<<VONMIS[it]<<setw(12)<<scientific<<TRACE[it]<<setw(12)<<scientific<<SIG11[it]<<setw(12)<<scientific<<SIG12[it]<<setw(12)<<scientific<<SIG13[it]<<setw(12)<<scientific<<SIG22[it]<<setw(12)<<scientific<<SIG23[it]<<setw(12)<<scientific<<SIG33[it]<<setw(12)<<scientific<<SIG1[it]<<setw(12)<<scientific<<SIG2[it]<<setw(12)<<scientific<<SIG3[it]<<endl; 
}

f2.close();

//ajout vphase - modif Ahmed Ammar 28/05/19

fstream f3;
f3.precision(4);
if(!bout){
f3.open("VISU/vphase",fstream::out | fstream::app);
}
else{
f3.open("VISU/vphase",fstream::out); 
f3<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;
}

f3<<setw(12)<<NB_SPH<<setw(12)<<0<<setw(12)<<scientific<<(dt*ite)<<endl;
for(int it=0;it<NB_SPH;it++){	
f3<<setw(12)<<scientific<<EDGE[it]<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<setw(12)<<LIST_P[it]<<endl; 
}

f3.close();


}

void ExpVISUC(int ite, R dt, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R, R* LIST_CONW, R & maxc, R & minc, bool * EDGE, R & mingx, R & maxgx, R & mingy, R & maxgy, R & mingz, R & maxgz, R * LIST_GCX, R * LIST_GCY, R * LIST_GCZ, bool * LIST_P, bool bout){
  
  // fichier visu
fstream f1;
f1.precision(4);
if(!bout){
f1.open("VISU/vconw",fstream::out | fstream::app);
}
else{
f1.open("VISU/vconw",fstream::out); 
f1<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;
}

f1<<setw(12)<<NB_SPH<<setw(12)<<0<<setw(12)<<scientific<<(dt*ite)<<endl;

f1<<setw(12)<<scientific<<minc<<setw(12)<<scientific<<maxc<<setw(12)<<scientific<<mingx<<setw(12)<<scientific<<maxgx<<setw(12)<<scientific<<mingy<<setw(12)<<scientific<<maxgy<<setw(12)<<scientific<<mingz<<setw(12)<<scientific<<maxgz<<endl;

for(int it=0;it<NB_SPH;it++){
f1<<setw(12)<<scientific<<EDGE[it]<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<setw(12)<<scientific<<LIST_CONW[it]<<setw(12)<<scientific<<LIST_GCX[it]<<setw(12)<<scientific<<LIST_GCY[it]<<setw(12)<<scientific<<LIST_GCZ[it]<<setw(12)<<LIST_P[it]<<endl; 
}


f1.close();
}

void ExpVISUT(int ite, R dt, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R, R* LIST_TEMP, R & maxt, R & mint, bool * EDGE, R & minfx, R & maxfx, R & minfy, R & maxfy, R & minfz, R & maxfz, R * LIST_FLX, R * LIST_FLY, R * LIST_FLZ, bool * LIST_P, bool bout){
  
  // fichier visu
fstream f1;
f1.precision(4);
if(!bout){
f1.open("VISU/vtemp",fstream::out | fstream::app);
}
else{
f1.open("VISU/vtemp",fstream::out); 
f1<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;
}

f1<<setw(12)<<NB_SPH<<setw(12)<<0<<setw(12)<<scientific<<(dt*ite)<<endl;

f1<<setw(12)<<scientific<<mint<<setw(12)<<scientific<<maxt<<setw(12)<<scientific<<minfx<<setw(12)<<scientific<<maxfx<<setw(12)<<scientific<<minfy<<setw(12)<<scientific<<maxfy<<setw(12)<<scientific<<minfz<<setw(12)<<scientific<<maxfz<<endl;

for(int it=0;it<NB_SPH;it++){
f1<<setw(12)<<scientific<<EDGE[it]<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<setw(12)<<scientific<<LIST_TEMP[it]<<setw(12)<<scientific<<LIST_FLX[it]<<setw(12)<<scientific<<LIST_FLY[it]<<setw(12)<<scientific<<LIST_FLZ[it]<<setw(12)<<LIST_P[it]<<endl; 
}


f1.close();
}


void ExpVISU_def(int ite, R dt, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R, R * EPSI11, R * EPSI22, R * EPSI33, R * EPSE11, R * EPSE22, R * EPSE33, R mindef11, R maxdef11,R mindef22, R maxdef22,R mindef33, R maxdef33, R mindefe11, R maxdefe11,R mindefe22, R maxdefe22,R mindefe33, R maxdefe33, bool bout, bool * EDGE){

fstream f1;
f1.precision(4);
if(!bout){
f1.open("VISU/vdef",fstream::out | fstream::app);
}
else{
f1.open("VISU/vdef",fstream::out); 
f1<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;
}

f1<<setw(12)<<NB_SPH<<setw(12)<<0<<setw(12)<<scientific<<(dt*ite)<<endl;

f1<<setw(12)<<scientific<<maxdef11<<setw(12)<<scientific<<mindef11<<endl;
f1<<setw(12)<<scientific<<maxdef22<<setw(12)<<scientific<<mindef22<<endl;
f1<<setw(12)<<scientific<<maxdef33<<setw(12)<<scientific<<mindef33<<endl;
f1<<setw(12)<<scientific<<maxdefe11<<setw(12)<<scientific<<mindefe11<<endl;
f1<<setw(12)<<scientific<<maxdefe22<<setw(12)<<scientific<<mindefe22<<endl;
f1<<setw(12)<<scientific<<maxdefe33<<setw(12)<<scientific<<mindefe33<<endl;

for(int it=0;it<NB_SPH;it++){
f1<<setw(12)<<scientific<<EDGE[it]<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<setw(12)<<scientific<<EPSI11[it]<<setw(12)<<scientific<<EPSI22[it]<<setw(12)<<scientific<<EPSI33[it]<<setw(12)<<scientific<<EPSE11[it]<<setw(12)<<scientific<<EPSE22[it]<<setw(12)<<scientific<<EPSE33[it]<<endl; 
}

f1.close();

}


void ExpVISU_rupt(int ite, R dt, int NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R, bool bout,int * LIST_B,  bool * EDGE){

fstream f1;
f1.precision(4);
if(!bout){
f1.open("VISU/vrupt",fstream::out | fstream::app);
}
else{
f1.open("VISU/vrupt",fstream::out); 
f1<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;
}

f1<<setw(12)<<NB_SPH<<setw(12)<<NB_PAR<<setw(12)<<scientific<<(dt*ite)<<endl;

for(int it=0;it<NB_PAR;it++){
f1<<setw(12)<<scientific<<LIST_PX[it][0]<<setw(12)<<scientific<<LIST_PY[it][0]<<setw(12)<<scientific<<LIST_PZ[it][0]<<endl;
f1<<setw(12)<<scientific<<LIST_PX[it][1]<<setw(12)<<scientific<<LIST_PY[it][1]<<setw(12)<<scientific<<LIST_PZ[it][1]<<endl;
f1<<setw(12)<<scientific<<LIST_PX[it][2]<<setw(12)<<scientific<<LIST_PY[it][2]<<setw(12)<<scientific<<LIST_PZ[it][2]<<endl;
f1<<setw(12)<<scientific<<LIST_PX[it][3]<<setw(12)<<scientific<<LIST_PY[it][3]<<setw(12)<<scientific<<LIST_PZ[it][3]<<endl;
}

for(int it=0;it<NB_SPH;it++){
f1<<setw(12)<<scientific<<EDGE[it]<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<setw(12)<<LIST_B[it]<<endl; 
}
f1.close();

}


void ExpVISU_COH(int ite, R dt, int NBCO, R H_TOT, R V_TOT, R Z_TOT, int ** CONT,R* LIST_X, R* LIST_Y, R* LIST_Z,R* LIST_R,bool bout){
  
  // fichier visu
fstream f1;
f1.precision(4);
if(!bout){
f1.open("VISU/vcoh",fstream::out | fstream::app);
}
else{
f1.open("VISU/vcoh",fstream::out); 
f1<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;
}

f1<<setw(12)<<NBCO<<setw(12)<<0<<setw(12)<<scientific<<(dt*ite)<<endl;
for(int it=0;it<NBCO;it++){
int i1=CONT[it][0];
int i2=CONT[it][1];

f1<<setw(12)<<scientific<<LIST_X[i1]<<setw(12)<<scientific<<LIST_Y[i1]<<setw(12)<<scientific<<LIST_Z[i1]<<setw(12)<<scientific<<LIST_X[i2]<<setw(12)<<scientific<<LIST_Y[i2]<<setw(12)<<scientific<<LIST_Z[i2]<<setw(12)<<scientific<<(LIST_R[i1]+LIST_R[i2])<<endl; 
}
f1.close();

}


void ExpMeshHAMZA(R SIZEX, R SIZEY, R SIZEZ, int NB_SPH, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R){
  
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

void ExpMeshMULTICOR(R SIZEX, R SIZEY, R SIZEZ, int NB_SPH, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R){

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

void ExpPARAMETER_sig_dep(R ff, int ite, R ti, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT, R* LIST_X, R* LIST_Y, R* LIST_Z, R* LIST_R,R * VONMIS,R * TRACE,R * SIG11,R * SIG12,R * SIG13,R * SIG22,R * SIG23,R * SIG33,R * SIG1,R * SIG2,R * SIG3,R * DEP1,R * DEP2,R * DEP3,R maxvm, R minvm, R maxtrac, R mintrac,R maxsig11, R minsig11,R maxsig12, R minsig12,R maxsig13, R minsig13,R maxsig22, R minsig22,R maxsig23, R minsig23,R maxsig33, R minsig33,R maxsig1, R minsig1,R maxsig2, R minsig2,R maxsig3, R minsig3,R maxdep1, R mindep1,R maxdep2, R mindep2,R maxdep3, R mindep3, int * LIST_B, bool * EDGE){
 
ostringstream oss;
oss << ite;


string filename0 = "DATA/datasigdep" + oss.str() + ".txt";
const char * filename1 = filename0.c_str();

  
fstream g;
g.precision(4);
g.open(filename1,fstream::out);
//g<<setw(12)<<ite<<setw(12)<<scientific<<ti<<setw(12)<<NB_SPH<<setw(12)<<scientific<<ff<<endl;
//g<<setw(12)<<scientific<<H_TOT<<setw(12)<<0.<<setw(12)<<scientific<<V_TOT<<setw(12)<<0.<<setw(12)<<scientific<<Z_TOT<<setw(12)<<0.<<endl;


for(int it=0;it<NB_SPH;it++){
g<<setw(12)<<scientific<<ff<<setw(12)<<(it+1)<<setw(12)<<EDGE[it]<<setw(12)<<scientific<<LIST_X[it]<<setw(12)<<scientific<<LIST_Y[it]<<setw(12)<<scientific<<LIST_Z[it]<<setw(12)<<scientific<<LIST_R[it]<<setw(12)<<scientific<<VONMIS[it]<<setw(12)<<scientific<<TRACE[it]<<setw(12)<<scientific<<SIG11[it]<<setw(12)<<scientific<<SIG12[it]<<setw(12)<<scientific<<SIG13[it]<<setw(12)<<scientific<<SIG22[it]<<setw(12)<<scientific<<SIG23[it]<<setw(12)<<scientific<<SIG33[it]<<setw(12)<<scientific<<SIG1[it]<<setw(12)<<scientific<<SIG2[it]<<setw(12)<<scientific<<SIG3[it]<<setw(12)<<scientific<<DEP1[it]<<setw(12)<<scientific<<DEP2[it]<<setw(12)<<scientific<<DEP3[it]<<endl; 
}
/*
g<<setw(12)<<scientific<<maxvm<<setw(12)<<scientific<<minvm<<endl;
g<<setw(12)<<scientific<<maxtrac<<setw(12)<<scientific<<mintrac<<endl;
g<<setw(12)<<scientific<<maxsig11<<setw(12)<<scientific<<minsig11<<endl;
g<<setw(12)<<scientific<<maxsig12<<setw(12)<<scientific<<minsig12<<endl;
g<<setw(12)<<scientific<<maxsig13<<setw(12)<<scientific<<minsig13<<endl;
g<<setw(12)<<scientific<<maxsig22<<setw(12)<<scientific<<minsig22<<endl;
g<<setw(12)<<scientific<<maxsig23<<setw(12)<<scientific<<minsig23<<endl;
g<<setw(12)<<scientific<<maxsig33<<setw(12)<<scientific<<minsig33<<endl;
g<<setw(12)<<scientific<<maxsig1<<setw(12)<<scientific<<minsig1<<endl;
g<<setw(12)<<scientific<<maxsig2<<setw(12)<<scientific<<minsig2<<endl;
g<<setw(12)<<scientific<<maxsig3<<setw(12)<<scientific<<minsig3<<endl;
g<<setw(12)<<scientific<<maxdep1<<setw(12)<<scientific<<mindep1<<endl;
g<<setw(12)<<scientific<<maxdep2<<setw(12)<<scientific<<mindep2<<endl;
g<<setw(12)<<scientific<<maxdep3<<setw(12)<<scientific<<mindep3<<endl;*/

g.close();
}


void Exp_debonding(int ite, R dt, R ftot,R ftot2, int nint0, int nint, int npre, int nsoft,bool bstart)
{
	
fstream g;
g.precision(4);

	if(bstart)
	{
	g.open("CONV/result_R",fstream::out);
	g<<setw(9)<<"ite"<<setw(18)<<"dt"<<setw(18)<<"ftot"<<setw(18)<<"ftot2"<<setw(18)<<"nint0"<<setw(18)<<"nint"<<setw(18)<<"npre"<<setw(18)<<"nsoft"<<endl;
	}else{
	g.open("CONV/result_R",fstream::out | fstream::app);	
	g<<setw(9)<<ite<<setw(18)<<scientific<<ite*dt<<setw(18)<<scientific<<ftot<<setw(18)<<scientific<<ftot2<<setw(18)<<nint0<<setw(18)<<nint<<setw(18)<<npre<<setw(18)<<nsoft<<endl;
	}	



g.close();				

}	  
	

