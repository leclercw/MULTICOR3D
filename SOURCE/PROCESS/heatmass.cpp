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

#define PI 3.14159265

#include "heatmass.h"

void Exp_temp_hydriq(int ite, int NB_SPH, R dt, R p, R coef1, R rhos, R Mw, R * LIST_TEMP, R * LIST_V,R * LIST_X,R * LIST_Y,R * LIST_Z, R * LIST_CONW, R H_TOT, R V_TOT, R Z_TOT,bool bstart){ 
	

R epsd=H_TOT*5e-2;

int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

int nb1=0;

R Ctot=0.;
R Ctot1=0.;
R Ttot1=0.;


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(+:nb1,Ctot1,Ttot1,Ctot)
	for(it=0;it<NB_SPH;it++)
	   {
	
	   Ctot+= LIST_CONW[it]*LIST_V[it]*coef1; 
       if((fabs(LIST_X[it]-H_TOT/2.)<epsd)&&(fabs(LIST_Y[it]-V_TOT/2.)<epsd)&&(fabs(LIST_Z[it]-Z_TOT/2.)<epsd)) {Ctot1+=LIST_CONW[it];  nb1++;  Ttot1+=LIST_TEMP[it];}

	   }

if(nb1==0) {cout<<"!!!!!!!!!!!!!!!!!!! epsd petit !!!!!!!!!!!!!!!!!!!!!!"<<endl;} else {Ctot1/=nb1;  Ttot1/=nb1;}
	
	
     fstream f1;
     f1.precision(4);
     if(bstart)
      {
	   f1.open("DATA/result_thyd",fstream::out);
       f1<<setw(12)<<scientific<<"Temps"<<setw(12)<<scientific<<"X"<<setw(12)<<scientific<<"Cmilieu"<<setw(12)<<scientific<<"Tmilieu"<<endl;
      }else{
       f1.open("DATA/result_thyd",fstream::out | fstream::app);
      }	

    f1<<setw(12)<<ite*dt<<setw(12)<<scientific<<Ctot*Mw/(H_TOT*V_TOT*Z_TOT*rhos*(1-p))<<setw(12)<<scientific<<Ctot1<<setw(12)<<scientific<<Ttot1<<endl;
    f1.close();
  	
	}

void temp_hydriq(int ite, R coef2,int NB_SPH, int NBCO, R dt, R Lambda_c,R p, R Cps, R Cpw, R coef1, int ** CONT, unsigned int ** NOCONT, R ** NCONT,int * NBCONTCO,int * NBHALO,R * VOLHALO, unsigned int **NOHALO, R *LIST_M, R rhos, R Mw,R Lw, R * LIST_R, R * LIST_TEMP, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4, bool * EDGE5, bool * EDGE6, R & mint, R & maxt, R * DCONTO, R * LIST_IND, R N1, R N2, R N3, R N4, R N5, R N6, R * DeltaT, R Dw, R * LIST_V, R * LIST_CONW, R * Cp_eff, R * LIST_GCX,R * LIST_GCY,R * LIST_GCZ, R & minc, R & maxc, R * DeltaC,R A, R B, R Rh,R Rg, R k, R h, R Tov, R Cov, R H_TOT, R V_TOT, R Z_TOT) 
{

int i,j;
R Sij;

R MoyR,coeftot=0.;
R coefpi=PI*coef2;
int it,num_threads;
R Flimp;

   if((dt*ite)<21600) Tov=20.+60.*ite*dt/21600.; else Tov=80.;
   Cov= 101325.*Rh*exp(13.7-5120./(273.15+Tov))/(Rg*(Tov+273.15));



# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

mint=01.E+20;
maxt=-1.E+20;

minc=01.E+20;
maxc=-1.E+20;


R Lambda_eff=Lambda_c*(1-p);
R rho_eff, Clim,aw,X;
R rhos_eff=rhos*(1-p);    
R coeft=Lambda_eff*dt/coef1;
R coefw=Dw*dt/coef1;

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaT[0:NB_SPH],DeltaC[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
  
              DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*coeft*Sij/(Cp_eff[i]*LIST_M[i]*DCONTO[it]); 
              DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*coeft*Sij/(Cp_eff[j]*LIST_M[j]*DCONTO[it]);
              DeltaC[i]+= (LIST_CONW[j]-LIST_CONW[i])*coefw*Sij/(LIST_V[i]*DCONTO[it]); 
              DeltaC[j]+= (LIST_CONW[i]-LIST_CONW[j])*coefw*Sij/(LIST_V[j]*DCONTO[it]);
	   }



// Thermique
R DLMw=Dw*Lw*Mw;
R dcoef1=dt*V_TOT*Z_TOT/(N1*coef1);
R dcoef2=dt*H_TOT*Z_TOT/(N2*coef1);
R dcoef3=dt*H_TOT*V_TOT/(N3*coef1);
R dcoef4=dt*V_TOT*Z_TOT/(N4*coef1);
R dcoef5=dt*H_TOT*Z_TOT/(N5*coef1);
R dcoef6=dt*H_TOT*V_TOT/(N6*coef1);

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,Flimp) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {

		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;

		if(EDGE[it]==0){

        	maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	

		}
		else if(EDGE[it]==1){

			if(ite!=1) grad_C_P(it, coef2, Dw,  CONT,  NCONT, LIST_R, LIST_CONW, DCONTO, LIST_GCX, LIST_GCY, LIST_GCZ, NBHALO, NOHALO,  NBCONTCO, NOCONT, VOLHALO, LIST_V, coef1);


			if(EDGE1[it]==1) { 
				   Flimp=h*(Tov-LIST_TEMP[it])+DLMw*(-LIST_GCX[it]); LIST_TEMP[it]+=Flimp*dcoef1*LIST_IND[it]/(Cp_eff[it]*LIST_M[it]);}

			if(EDGE2[it]==1) { 
				   Flimp=h*(Tov-LIST_TEMP[it])+DLMw*(-LIST_GCY[it]);LIST_TEMP[it]+=Flimp*dcoef2*LIST_IND[it]/(Cp_eff[it]*LIST_M[it]);}


			if(EDGE3[it]==1) { 
				   Flimp=h*(Tov-LIST_TEMP[it])+DLMw*(-LIST_GCZ[it]);LIST_TEMP[it]+=Flimp*dcoef3*LIST_IND[it]/(Cp_eff[it]*LIST_M[it]);}


			if(EDGE4[it]==1) { 
				   Flimp=h*(Tov-LIST_TEMP[it])+DLMw*(LIST_GCX[it]);LIST_TEMP[it]+=Flimp*dcoef4*LIST_IND[it]/(Cp_eff[it]*LIST_M[it]);}


			if(EDGE5[it]==1) { 
				   Flimp=h*(Tov-LIST_TEMP[it])+DLMw*(LIST_GCY[it]);LIST_TEMP[it]+=Flimp*dcoef5*LIST_IND[it]/(Cp_eff[it]*LIST_M[it]);}

			if(EDGE6[it]==1) {
				   Flimp=h*(Tov-LIST_TEMP[it])+DLMw*(LIST_GCZ[it]); LIST_TEMP[it]+=Flimp*dcoef6*LIST_IND[it]/(Cp_eff[it]*LIST_M[it]);}


		}

         
	   }



// Hydrique
dcoef1=k*dcoef1;
dcoef2=k*dcoef2;
dcoef3=k*dcoef3;
dcoef4=k*dcoef4;
dcoef5=k*dcoef5;
dcoef6=k*dcoef6;

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,X,aw,Clim,rho_eff) reduction(min:minc) reduction(max:maxc)
	for(it=0;it<NB_SPH;it++)
	   {

	    LIST_CONW[it]+=DeltaC[it];
	    DeltaC[it]=0.;

		if(EDGE[it]==0){

		maxc=max(maxc,LIST_CONW[it]);
		minc=min(minc,LIST_CONW[it]);

		}
		else if(EDGE[it]==1){

		X = LIST_CONW[it]*Mw/rhos_eff;
		aw= pow(1.+pow(A/X,B),-1);
		Clim= 101325.*aw*exp(13.7-5120./(273.15+LIST_TEMP[it]))/(Rg*(LIST_TEMP[it]+273.15));
		if(EDGE1[it]==1)      LIST_CONW[it]+=dcoef1*LIST_IND[it]*(Cov-Clim)/LIST_V[it];
		if(EDGE2[it]==1)      LIST_CONW[it]+=dcoef2*LIST_IND[it]*(Cov-Clim)/LIST_V[it];
		if(EDGE3[it]==1)      LIST_CONW[it]+=dcoef3*LIST_IND[it]*(Cov-Clim)/LIST_V[it];
		if(EDGE4[it]==1)      LIST_CONW[it]+=dcoef4*LIST_IND[it]*(Cov-Clim)/LIST_V[it];
		if(EDGE5[it]==1)      LIST_CONW[it]+=dcoef5*LIST_IND[it]*(Cov-Clim)/LIST_V[it];
		if(EDGE6[it]==1)      LIST_CONW[it]+=dcoef6*LIST_IND[it]*(Cov-Clim)/LIST_V[it];

		}


                rho_eff=rhos_eff+Mw*LIST_CONW[it];
                Cp_eff[it]=(Cps*rhos_eff+Cpw*(rho_eff-rhos_eff))/rho_eff;
		LIST_M[it]=LIST_V[it]*rho_eff;

	   }



}








void hydriq(R coef2,int NB_SPH, int NBCO, R dt, R Cov,R Clim, R k, R D_w, R coef1, int ** CONT, R *LIST_V, R * LIST_R, R * LIST_CONW, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, R & minc, R & maxc, R * DCONTO, R * LIST_IND, R N1, R N2, R N3, R N4, R N5, R N6, R * DeltaC, R H_TOT, R V_TOT, R Z_TOT, int dir, bool btype)  
{

int i,j;
R Sij;

R MoyR;
R coefpi=PI*coef2;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

minc=01.E+20;
maxc=-1.E+20;

if(dir==1 && btype==0){
// cas dir 1 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaC[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;	    
 
	    if((not EDGE4[i]) and (not EDGE1[i])) DeltaC[i]+= (LIST_CONW[j]-LIST_CONW[i])*D_w*dt*Sij/(LIST_V[i]*coef1*DCONTO[it]); 
	    if((not EDGE4[j]) and (not EDGE1[j])) DeltaC[j]+= (LIST_CONW[i]-LIST_CONW[j])*D_w*dt*Sij/(LIST_V[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:minc) reduction(max:maxc)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_CONW[it]+=DeltaC[it];
		DeltaC[it]=0.;
		if(EDGE[it]==0){
		maxc=max(maxc,LIST_CONW[it]);
		minc=min(minc,LIST_CONW[it]);	
		}
	   }


}else if(dir==2 && btype==0){
// cas dir 2 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaC[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	    
	    if((not EDGE5[i]) and (not EDGE2[i])) DeltaC[i]+= (LIST_CONW[j]-LIST_CONW[i])*D_w*dt*Sij/(LIST_V[i]*coef1*DCONTO[it]); 
	    if((not EDGE5[j]) and (not EDGE2[j])) DeltaC[j]+= (LIST_CONW[i]-LIST_CONW[j])*D_w*dt*Sij/(LIST_V[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:minc) reduction(max:maxc)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_CONW[it]+=DeltaC[it];
		DeltaC[it]=0.;
		if(EDGE[it]==0){
		maxc=max(maxc,LIST_CONW[it]);
		minc=min(minc,LIST_CONW[it]);	
		}
	   }


}else if(dir==3 && btype==0){
// cas dir 3 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaC[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;

	    if((not EDGE6[i]) and (not EDGE3[i])) DeltaC[i]+= (LIST_CONW[j]-LIST_CONW[i])*D_w*dt*Sij/(LIST_V[i]*coef1*DCONTO[it]); 
	    if((not EDGE6[j]) and (not EDGE3[j])) DeltaC[j]+= (LIST_CONW[i]-LIST_CONW[j])*D_w*dt*Sij/(LIST_V[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:minc) reduction(max:maxc)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_CONW[it]+=DeltaC[it];
		DeltaC[it]=0.;
		if(EDGE[it]==0){
		maxc=max(maxc,LIST_CONW[it]);
		minc=min(minc,LIST_CONW[it]);	
		}
	   }

}

else if(dir==1 && btype==1){
// cas dir 1 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaC[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
  
	     if(not EDGE1[i]) DeltaC[i]+= (LIST_CONW[j]-LIST_CONW[i])*D_w*dt*Sij/(LIST_V[i]*coef1*DCONTO[it]); 
	     if(not EDGE1[j]) DeltaC[j]+= (LIST_CONW[i]-LIST_CONW[j])*D_w*dt*Sij/(LIST_V[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:minc) reduction(max:maxc)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_CONW[it]+=DeltaC[it];
		DeltaC[it]=0.;
      if(EDGE4[it])      LIST_CONW[it]+=dt*LIST_IND[it]/(N4*LIST_V[it]*coef1)*(k*(Cov-Clim)*V_TOT*Z_TOT);
		if(EDGE[it]==0){
		maxc=max(maxc,LIST_CONW[it]);
		minc=min(minc,LIST_CONW[it]);	
		}
	   }

}
else if(dir==2 && btype==1){
// cas dir 2 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij,D_w) reduction(+:DeltaC[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	     
	     if(not EDGE2[i]) DeltaC[i]+= (LIST_CONW[j]-LIST_CONW[i])*D_w*dt*Sij/(LIST_V[i]*coef1*DCONTO[it]); 
	     if(not EDGE2[j]) DeltaC[j]+= (LIST_CONW[i]-LIST_CONW[j])*D_w*dt*Sij/(LIST_V[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:minc) reduction(max:maxc)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_CONW[it]+=DeltaC[it];
		DeltaC[it]=0.;

      if(EDGE5[it])      LIST_CONW[it]+=dt*LIST_IND[it]/(N5*LIST_V[it]*coef1)*(k*(Cov-Clim)*H_TOT*Z_TOT);

		if(EDGE[it]==0){
		maxc=max(maxc,LIST_CONW[it]);
		minc=min(minc,LIST_CONW[it]);	
		}
	   }

}else if(dir==3 && btype==1){
// cas dir 3 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaC[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
  
	     if(not EDGE3[i]) DeltaC[i]+= (LIST_CONW[j]-LIST_CONW[i])*D_w*dt*Sij/(LIST_V[i]*coef1*DCONTO[it]); 
	     if(not EDGE3[j]) DeltaC[j]+= (LIST_CONW[i]-LIST_CONW[j])*D_w*dt*Sij/(LIST_V[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:minc) reduction(max:maxc)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_CONW[it]+=DeltaC[it];
		DeltaC[it]=0.;
      if(EDGE6[it])      LIST_CONW[it]+=dt*LIST_IND[it]/(N6*LIST_V[it]*coef1)*(k*(Cov-Clim)*H_TOT*V_TOT);
		if(EDGE[it]==0){
		maxc=max(maxc,LIST_CONW[it]);
		minc=min(minc,LIST_CONW[it]);	
		}
	   }

}




}


void ther_load_dila(int ite, int nitd,bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO, bool * LIST_P,R alpha1, R alpha2){

R pr,alpham;
if(ite<=nitd){
pr=R(ite)/nitd;
}else{
pr=1.;
}

int it1,it2;
int i,num_threads;


# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(i,it1,it2,alpham) 
	for (i=0;i<NBCO;i++){

		if(TYPCO[i]){
	
			it1=CONT[i][0];
			it2=CONT[i][1];

			if((!LIST_P[it1])&&(!LIST_P[it2])){
                        alpham=alpha1*pr;    
			}
			else if((LIST_P[it1])&&(LIST_P[it2])){
                        alpham=alpha2*pr;
			}
                        else{
                        alpham=sqrt(alpha1*alpha2)*pr;
                        }
             
			DCONT[i]=DCONTO[i]*(1.+alpham);
			DCONTX[i]=DCONTXO[i]*(1.+alpham);
			DCONTY[i]=DCONTYO[i]*(1.+alpham);
			DCONTZ[i]=DCONTZO[i]*(1.+alpham);

		}

	}


}

void ther_load_dilat(int ite, R Tinit, R * LIST_TEMP, bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO, bool * LIST_P,R alpha1, R alpha2){

int it1,it2;
int i,num_threads;
R alpham,alphami1,alphami2;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(i,it1,it2,alpham,alphami1,alphami2) 
	for (i=0;i<NBCO;i++){

		if(TYPCO[i]){
	
			it1=CONT[i][0];
			it2=CONT[i][1];

			if(!LIST_P[it1]){
			alphami1=alpha1*(LIST_TEMP[it1]-Tinit);
			}else{
			alphami1=alpha2*(LIST_TEMP[it1]-Tinit);
			}
		
			if(!LIST_P[it2]){
			alphami2=alpha1*(LIST_TEMP[it2]-Tinit);
			}else{
			alphami2=alpha2*(LIST_TEMP[it2]-Tinit);
			}
            
            alpham=(alphami1+alphami2)/2.; 
             
			DCONT[i]=DCONTO[i]*(1.+alpham);
			DCONTX[i]=DCONTXO[i]*(1.+alpham);
			DCONTY[i]=DCONTYO[i]*(1.+alpham);
			DCONTZ[i]=DCONTZO[i]*(1.+alpham);

		}

	}


}

void ther_load_hygr(int ite, int nitd,bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO, bool * LIST_P, R alphag1, R alphag2){

R pr,betam1,betam2,betam;
if(ite<=nitd){
pr=R(ite)/nitd;
}else{
pr=1.;
}
betam1=exp(-alphag1*pr);
betam2=exp(-alphag2*pr);

int it1,it2;
int i,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(i,it1,it2,betam) 
	for (i=0;i<NBCO;i++){

		if(TYPCO[i]){
	
			it1=CONT[i][0];
			it2=CONT[i][1];

			if((!LIST_P[it1])&&(!LIST_P[it2])){
			betam=betam1;    
			}
			else if((LIST_P[it1])&&(LIST_P[it2])){
			betam=betam2;
			}
			else{
			betam=sqrt(betam1*betam2);
			}             
                          
			DCONT[i]=DCONTO[i]*betam;
			DCONTX[i]=DCONTXO[i]*betam;
			DCONTY[i]=DCONTYO[i]*betam;
			DCONTZ[i]=DCONTZO[i]*betam;

		}

	}


}

void ther_load_hygrc(int NBENREG,int ite,R Mw,R rhow, R * LIST_CONW,R Cinit, bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO){

int it1,it2;
int i,num_threads;
R eps1,eps2,epsl,betal;
R epsm=0.;
R epsd=0.;
int ncom=0;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(i,it1,it2,eps1,eps2,epsl,betal) reduction(+:epsm,epsd,ncom)
	for (i=0;i<NBCO;i++){

		if(TYPCO[i]){
	
			ncom++;
	
			it1=CONT[i][0];
			it2=CONT[i][1];
       
			eps1=Mw*(LIST_CONW[it1]-Cinit)/rhow;
			eps2=Mw*(LIST_CONW[it2]-Cinit)/rhow;
			epsl=(eps1+eps2)/2.;
                        epsm+=epsl;
                        epsl=pow((1.+epsl),1./3)-1.;
                        epsd+=epsl;                        
			betal=1.+epsl;
             
		//	DCONT[i]=DCONTO[i]*betal;
			DCONTX[i]=DCONTXO[i]*betal;
			DCONTY[i]=DCONTYO[i]*betal;
            DCONTZ[i]=DCONTZO[i]*betal;

		}

	}

epsm/=ncom;
epsd/=ncom;

if(ite%NBENREG==0) {
cout<<"Retrait vol (%):"<<epsm*100.<<endl;
cout<<"Retrait lin (%):"<<epsd*100.<<endl;
}

}

void ther_load_hygrc_dilat(int NBENREG,int ite,R Mw,R rhow, R * LIST_CONW,R Cinit,R * LIST_TEMP, R Tinit,  bool * TYPCO, int NBCO, int ** CONT, R * DCONT, R * DCONTX,R * DCONTY,R * DCONTZ, R * DCONTO, R * DCONTXO,R * DCONTYO,R * DCONTZO,bool * LIST_P, R alpha1, R alpha2){

int it1,it2;
int i,num_threads;
R eps1,eps2,epsl,betal;
R epsm=0.;
R epsd=0.;
int ncom=0;
R alpham,alphami1,alphami2;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(i,it1,it2,eps1,eps2,epsl,betal,alpham,alphami1,alphami2) reduction(+:epsm,epsd,ncom)
	for (i=0;i<NBCO;i++){

		if(TYPCO[i]){
	
			ncom++;
	
			it1=CONT[i][0];
			it2=CONT[i][1];
       
			eps1=Mw*(LIST_CONW[it1]-Cinit)/rhow;
			eps2=Mw*(LIST_CONW[it2]-Cinit)/rhow;
			epsl=(eps1+eps2)/2.;
			epsm+=epsl;
			epsl=pow((1.+epsl),1./3)-1.;
			epsd+=epsl;                        
			betal=1.+epsl;

			if(!LIST_P[it1]){
			alphami1=alpha1*(LIST_TEMP[it1]-Tinit);
			}else{
			alphami1=alpha2*(LIST_TEMP[it1]-Tinit);
			}
		
			if(!LIST_P[it2]){
			alphami2=alpha1*(LIST_TEMP[it2]-Tinit);
			}else{
			alphami2=alpha2*(LIST_TEMP[it2]-Tinit);
			}
            
            alpham=(alphami1+alphami2)/2.; 
             
		//	DCONT[i]=DCONTO[i]*betal*(1.+alpham);
			DCONTX[i]=DCONTXO[i]*betal*(1.+alpham);
			DCONTY[i]=DCONTYO[i]*betal*(1.+alpham);
            DCONTZ[i]=DCONTZO[i]*betal*(1.+alpham);

		}

	}

epsm/=ncom;
epsd/=ncom;

if(ite%NBENREG==0) {
cout<<"Retrait vol (%):"<<epsm*100.<<endl;
cout<<"Retrait lin (%):"<<epsd*100.<<endl;
}

}

void temp(R coef2,int NB_SPH, int NBCO, R dt, R Flimp, R Lambda_c, R Cp, R coef1, int ** CONT, R *LIST_M, R * LIST_R, R * LIST_TEMP,bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, R & mint, R & maxt, R * DCONTO, R * LIST_IND, R N1, R N2, R N3, R N4, R N5, R N6, R * DeltaT, int dir, bool btype)  
{

int i,j;
R Sij;

R MoyR;
R coefpi=PI*coef2;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

mint=01.E+20;
maxt=-1.E+20;

if(dir==1 && btype==0){
// cas dir 1 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;	    
 
	    if((not EDGE4[i]) and (not EDGE1[i])) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	    if((not EDGE4[j]) and (not EDGE1[j])) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }


}else if(dir==2 && btype==0){
// cas dir 2 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	    
	    if((not EDGE5[i]) and (not EDGE2[i])) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	    if((not EDGE5[j]) and (not EDGE2[j])) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }


}else if(dir==3 && btype==0){
// cas dir 3 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;

	    if((not EDGE6[i]) and (not EDGE3[i])) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	    if((not EDGE6[j]) and (not EDGE3[j])) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }

}else if(dir==1 && btype==1){
// cas dir 1 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
  
	     if(not EDGE1[i]) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	     if(not EDGE1[j]) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
	   if(EDGE4[it]) LIST_TEMP[it]+=dt*Flimp*LIST_IND[it]/(N4*Cp*LIST_M[it]*coef1);
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }

}else if(dir==2 && btype==1){
// cas dir 2 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij,Lambda_c) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	     
	     if(not EDGE2[i]) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	     if(not EDGE2[j]) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
	   if(EDGE5[it]) LIST_TEMP[it]+=dt*Flimp*LIST_IND[it]/(N5*Cp*LIST_M[it]*coef1);
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }

}else if(dir==3 && btype==1){
// cas dir 3 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
  
	     if(not EDGE3[i]) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	     if(not EDGE3[j]) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
	   if(EDGE6[it]) LIST_TEMP[it]+=dt*Flimp*LIST_IND[it]/(N6*Cp*LIST_M[it]*coef1);
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }

}




}




void temp2(R coef2,int NB_SPH, int NBCO, R dt, R Flimp, R lambda1, R lambda2, R cp1, R cp2, R coef1, int ** CONT, R *LIST_M, R * LIST_R, R * LIST_TEMP,bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, R & mint, R & maxt, R * DCONTO, R * LIST_IND, bool * LIST_P, R N1, R N2, R N3, R N4, R N5, R N6, R * DeltaT, int dir, bool btype)  
{

int i,j;
R Sij;
R Lambda_c,Cp;

R MoyR;
R coefpi=PI*coef2;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

mint=01.E+20;
maxt=-1.E+20;

if(dir==1 && btype==0){
// cas dir 1 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij,Lambda_c,Cp) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	    
	    if (LIST_P[i]==1 && LIST_P[j]==1) {Lambda_c=lambda2; Cp=cp2; }
	    else if (LIST_P[i]==0 && LIST_P[j]==0) {Lambda_c=lambda1;  Cp=cp1;}
	    else {Lambda_c=2.*(lambda2*lambda1)/(lambda2+lambda1); Cp=2.*(cp1*cp2)/(cp1+cp2); }
  
	    if((not EDGE4[i]) and (not EDGE1[i])) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	    if((not EDGE4[j]) and (not EDGE1[j])) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,Cp) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}	
	   }


}else if(dir==2 && btype==0){
// cas dir 2 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij,Lambda_c,Cp) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	    
	    if (LIST_P[i]==1 && LIST_P[j]==1) {Lambda_c=lambda2; Cp=cp2; }
	    else if (LIST_P[i]==0 && LIST_P[j]==0) {Lambda_c=lambda1;  Cp=cp1;}
	    else {Lambda_c=2.*(lambda2*lambda1)/(lambda2+lambda1); Cp=2.*(cp1*cp2)/(cp1+cp2); }
  
	    if((not EDGE5[i]) and (not EDGE2[i])) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	    if((not EDGE5[j]) and (not EDGE2[j])) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,Cp) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }


}else if(dir==3 && btype==0){
// cas dir 3 & temp

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij,Lambda_c,Cp) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	    
	    if (LIST_P[i]==1 && LIST_P[j]==1) {Lambda_c=lambda2; Cp=cp2; }
	    else if (LIST_P[i]==0 && LIST_P[j]==0) {Lambda_c=lambda1;  Cp=cp1;}
	    else {Lambda_c=2.*(lambda2*lambda1)/(lambda2+lambda1); Cp=2.*(cp1*cp2)/(cp1+cp2); }
  
	    if((not EDGE6[i]) and (not EDGE3[i])) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	    if((not EDGE6[j]) and (not EDGE3[j])) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,Cp) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }

}else if(dir==1 && btype==1){
// cas dir 1 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij,Lambda_c,Cp) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	    
	    if (LIST_P[i]==1 && LIST_P[j]==1) {Lambda_c=lambda2; Cp=cp2; }
	    else if (LIST_P[i]==0 && LIST_P[j]==0) {Lambda_c=lambda1;  Cp=cp1;}
	    else {Lambda_c=2.*(lambda2*lambda1)/(lambda2+lambda1); Cp=2.*(cp1*cp2)/(cp1+cp2); }
  
	     if(not EDGE1[i]) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	     if(not EDGE1[j]) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,Cp) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
	   if (LIST_P[it]==1) { Cp=cp2; }
	    else { Cp=cp1; }
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
	   if(EDGE4[it]) LIST_TEMP[it]+=dt*Flimp*LIST_IND[it]/(N4*Cp*LIST_M[it]*coef1);
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }

}else if(dir==2 && btype==1){
// cas dir 2 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij,Lambda_c,Cp) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	    
	    if (LIST_P[i]==1 && LIST_P[j]==1) {Lambda_c=lambda2; Cp=cp2; }
	    else if (LIST_P[i]==0 && LIST_P[j]==0) {Lambda_c=lambda1;  Cp=cp1;}
	    else {Lambda_c=2.*(lambda2*lambda1)/(lambda2+lambda1); Cp=2.*(cp1*cp2)/(cp1+cp2); }
  
	     if(not EDGE2[i]) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	     if(not EDGE2[j]) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,Cp) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
	   if (LIST_P[it]==1) { Cp=cp2; }
	    else { Cp=cp1; }
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
	   if(EDGE5[it]) LIST_TEMP[it]+=dt*Flimp*LIST_IND[it]/(N5*Cp*LIST_M[it]*coef1);
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }

}else if(dir==3 && btype==1){
// cas dir 3 & flu

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,i,j,MoyR,Sij,Lambda_c,Cp) reduction(+:DeltaT[0:NB_SPH]) 
	for(it=0;it<NBCO;it++)
	   {

	    // Numeros des candidats 
	    i=CONT[it][0];
	    j=CONT[it][1]; 

	    MoyR=(LIST_R[i]+LIST_R[j])/2.;
	    Sij=MoyR*MoyR*coefpi;
	    
	    if (LIST_P[i]==1 && LIST_P[j]==1) {Lambda_c=lambda2; Cp=cp2; }
	    else if (LIST_P[i]==0 && LIST_P[j]==0) {Lambda_c=lambda1;  Cp=cp1;}
	    else {Lambda_c=2.*(lambda2*lambda1)/(lambda2+lambda1); Cp=2.*(cp1*cp2)/(cp1+cp2); }
  
	     if(not EDGE3[i]) DeltaT[i]+= (LIST_TEMP[j]-LIST_TEMP[i])*Lambda_c*dt*Sij/(Cp*LIST_M[i]*coef1*DCONTO[it]); 
	     if(not EDGE3[j]) DeltaT[j]+= (LIST_TEMP[i]-LIST_TEMP[j])*Lambda_c*dt*Sij/(Cp*LIST_M[j]*coef1*DCONTO[it]);  
	   }


# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,Cp) reduction(min:mint) reduction(max:maxt)
	for(it=0;it<NB_SPH;it++)
	   {
	   if (LIST_P[it]==1) { Cp=cp2; }
	    else { Cp=cp1; }
		LIST_TEMP[it]+=DeltaT[it];
		DeltaT[it]=0.;
	   if(EDGE6[it]) LIST_TEMP[it]+=dt*Flimp*LIST_IND[it]/(N6*Cp*LIST_M[it]*coef1);
		if(EDGE[it]==0){
		maxt=max(maxt,LIST_TEMP[it]);
		mint=min(mint,LIST_TEMP[it]);	
		}
	   }

}


}






void flux(R coef2,int NB_SPH, int NBCO, R dt, R Lambda_eff, R coef1, int ** CONT, R ** NCONT, R *LIST_M, R * LIST_R, R * LIST_TEMP,bool * EDGE, R & minfx, R & maxfx, R & minfy, R & maxfy, R & minfz, R & maxfz, R * DCONTO, R H_TOT, R V_TOT, R Z_TOT, R * LIST_FLX,R * LIST_FLY,R * LIST_FLZ,int * NBHALO, unsigned int **NOHALO, int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO, R * LIST_V)  
{

R DeltaF;
int jt,kt;
unsigned int numcor,numc;

int i,j;
R Sij;

R MoyR;
R coefpi=PI*coef2;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) 
	for(it=0;it<NB_SPH;it++)
	   {
	       LIST_FLX[it]=0.;
	       LIST_FLY[it]=0.;
	       LIST_FLZ[it]=0.;
	   }

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,kt,numc,numcor,i,j,MoyR,Sij,DeltaF) reduction(-:LIST_FLX[0:NB_SPH],LIST_FLY[0:NB_SPH],LIST_FLZ[0:NB_SPH])  
	for(jt=0;jt<NB_SPH;jt++)
	   {

		    for(it=0;it<NBHALO[jt];it++){
			numcor=NOHALO[jt][it];

				for(kt=0;kt<NBCONTCO[numcor];kt++){ 

					numc=NOCONT[numcor][kt];
					i=CONT[numc][0];
					j=CONT[numc][1];

					MoyR=(LIST_R[i]+LIST_R[j])/2.;
					Sij=MoyR*MoyR*coefpi;

					DeltaF=(LIST_TEMP[j]-LIST_TEMP[i])*Lambda_eff/DCONTO[numc];

					LIST_FLX[jt]-=DeltaF*NCONT[numc][0]*Sij*DCONTO[numc]/2;
					LIST_FLY[jt]-=DeltaF*NCONT[numc][1]*Sij*DCONTO[numc]/2;
					LIST_FLZ[jt]-=DeltaF*NCONT[numc][2]*Sij*DCONTO[numc]/2;

				}
		    }
	      
	     LIST_FLX[jt]/=VOLHALO[jt];      
	     LIST_FLY[jt]/=VOLHALO[jt];     
	     LIST_FLZ[jt]/=VOLHALO[jt];   
      		    
	   }


minfx=01.E+20;
maxfx=-1.E+20;

minfy=01.E+20;
maxfy=-1.E+20;

minfz=01.E+20;
maxfz=-1.E+20;



R FluTotx=0.;
R FluToty=0.;
R FluTotz=0.;

R volel;

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,volel) reduction(max:maxfx,maxfy,maxfz) reduction(min:minfx,minfy,minfz) reduction (+:FluTotx,FluToty,FluTotz)
	for(it=0;it<NB_SPH;it++)
	   {

		if(EDGE[it]==0){
		maxfx=(LIST_FLX[it]>maxfx)?LIST_FLX[it]:maxfx; 
		minfx=(LIST_FLX[it]<minfx)?LIST_FLX[it]:minfx; 

		maxfy=(LIST_FLY[it]>maxfy)?LIST_FLY[it]:maxfy; 
		minfy=(LIST_FLY[it]<minfy)?LIST_FLY[it]:minfy; 

		maxfz=(LIST_FLZ[it]>maxfz)?LIST_FLZ[it]:maxfz; 
		minfz=(LIST_FLZ[it]<minfz)?LIST_FLZ[it]:minfz; 
		}

		volel=coef1*LIST_V[it];

		FluTotx+=LIST_FLX[it]*volel;
		FluToty+=LIST_FLY[it]*volel;
		FluTotz+=LIST_FLZ[it]*volel;
	    
	      
	   }

cout<<"FLX "<<FluTotx/(H_TOT*V_TOT*Z_TOT)<<endl;
cout<<"FLY "<<FluToty/(H_TOT*V_TOT*Z_TOT)<<endl;
cout<<"FLZ "<<FluTotz/(H_TOT*V_TOT*Z_TOT)<<endl;

}







void flux2(R coef2,int NB_SPH, int NBCO, R dt, R lambda_eff1, R lambda_eff2, R coef1, int ** CONT, R ** NCONT, R *LIST_M, R * LIST_R, R * LIST_TEMP,bool * EDGE, R & minfx, R & maxfx, R & minfy, R & maxfy, R & minfz, R & maxfz, R * DCONTO,  bool * LIST_P, R H_TOT, R V_TOT, R Z_TOT, R * LIST_FLX,R * LIST_FLY,R * LIST_FLZ,int * NBHALO, unsigned int **NOHALO, int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO, R * LIST_V)  
{

R DeltaF;
int jt,kt;
unsigned int numcor,numc;

int i,j;
R Sij;
R Lambda_eff;

R MoyR;
R coefpi=PI*coef2;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) 
	for(it=0;it<NB_SPH;it++)
	   {
	       LIST_FLX[it]=0.;
	       LIST_FLY[it]=0.;
	       LIST_FLZ[it]=0.;
	   }

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,kt,numc,numcor,i,j,MoyR,Sij,DeltaF,Lambda_eff) reduction(-:LIST_FLX[0:NB_SPH],LIST_FLY[0:NB_SPH],LIST_FLZ[0:NB_SPH])  
	for(jt=0;jt<NB_SPH;jt++)
	   {

		    for(it=0;it<NBHALO[jt];it++){
			numcor=NOHALO[jt][it];

				for(kt=0;kt<NBCONTCO[numcor];kt++){ 

					numc=NOCONT[numcor][kt];
					i=CONT[numc][0];
					j=CONT[numc][1];

					if (LIST_P[i]==1 && LIST_P[j]==1) Lambda_eff=lambda_eff2;
					else if (LIST_P[i]==0 && LIST_P[j]==0) Lambda_eff=lambda_eff1;
					else Lambda_eff=2.*(lambda_eff2*lambda_eff1)/(lambda_eff2+lambda_eff1);

					MoyR=(LIST_R[i]+LIST_R[j])/2.;
					Sij=MoyR*MoyR*coefpi;

					DeltaF=(LIST_TEMP[j]-LIST_TEMP[i])*Lambda_eff/DCONTO[numc];

					LIST_FLX[jt]-=DeltaF*NCONT[numc][0]*Sij*DCONTO[numc]/2;
					LIST_FLY[jt]-=DeltaF*NCONT[numc][1]*Sij*DCONTO[numc]/2;
					LIST_FLZ[jt]-=DeltaF*NCONT[numc][2]*Sij*DCONTO[numc]/2;

				}
		    }
	      
	     LIST_FLX[jt]/=VOLHALO[jt];      
	     LIST_FLY[jt]/=VOLHALO[jt];     
	     LIST_FLZ[jt]/=VOLHALO[jt];   
      		    
	   }


minfx=01.E+20;
maxfx=-1.E+20;

minfy=01.E+20;
maxfy=-1.E+20;

minfz=01.E+20;
maxfz=-1.E+20;



R FluTotx=0.;
R FluToty=0.;
R FluTotz=0.;

R volel;

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,volel) reduction(max:maxfx,maxfy,maxfz) reduction(min:minfx,minfy,minfz) reduction (+:FluTotx,FluToty,FluTotz)
	for(it=0;it<NB_SPH;it++)
	   {

		if(EDGE[it]==0){
		maxfx=(LIST_FLX[it]>maxfx)?LIST_FLX[it]:maxfx; 
		minfx=(LIST_FLX[it]<minfx)?LIST_FLX[it]:minfx; 

		maxfy=(LIST_FLY[it]>maxfy)?LIST_FLY[it]:maxfy; 
		minfy=(LIST_FLY[it]<minfy)?LIST_FLY[it]:minfy; 

		maxfz=(LIST_FLZ[it]>maxfz)?LIST_FLZ[it]:maxfz; 
		minfz=(LIST_FLZ[it]<minfz)?LIST_FLZ[it]:minfz; 
		}

		volel=coef1*LIST_V[it];

		FluTotx+=LIST_FLX[it]*volel;
		FluToty+=LIST_FLY[it]*volel;
		FluTotz+=LIST_FLZ[it]*volel;
	    
	      
	   }

cout<<"FLX "<<FluTotx/(H_TOT*V_TOT*Z_TOT)<<endl;
cout<<"FLY "<<FluToty/(H_TOT*V_TOT*Z_TOT)<<endl;
cout<<"FLZ "<<FluTotz/(H_TOT*V_TOT*Z_TOT)<<endl;

}


void grad_C(R coef2,int NB_SPH, int NBCO, R dt, R Dw, R coef1, int ** CONT, R ** NCONT, R *LIST_M, R * LIST_R, R * LIST_CONW,bool * EDGE, R & mingx, R & maxgx, R & mingy, R & maxgy, R & mingz, R & maxgz, R * DCONTO, R H_TOT, R V_TOT, R Z_TOT, R * LIST_GCX,R * LIST_GCY,R * LIST_GCZ,int * NBHALO,unsigned int **NOHALO, int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO, R * LIST_V)  
{

R DeltaC;
int jt,kt,numcor,numc;

int i,j;
R Sij;

R MoyR;
R coefpi=PI*coef2;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) 
	for(it=0;it<NB_SPH;it++)
	   {
	       LIST_GCX[it]=0.;
	       LIST_GCY[it]=0.;
	       LIST_GCZ[it]=0.;
	   }

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,kt,numc,numcor,i,j,MoyR,Sij,DeltaC) reduction(-:LIST_GCX[0:NB_SPH],LIST_GCY[0:NB_SPH],LIST_GCZ[0:NB_SPH])  
	for(jt=0;jt<NB_SPH;jt++)
	   {

		    for(it=0;it<NBHALO[jt];it++){
			numcor=NOHALO[jt][it];

				for(kt=0;kt<NBCONTCO[numcor];kt++){ 

					numc=NOCONT[numcor][kt];
					i=CONT[numc][0];
					j=CONT[numc][1];

					MoyR=(LIST_R[i]+LIST_R[j])/2.;
					Sij=MoyR*MoyR*coefpi;

					DeltaC=(LIST_CONW[j]-LIST_CONW[i])/DCONTO[numc];

					LIST_GCX[jt]-=DeltaC*NCONT[numc][0]*Sij*DCONTO[numc]/2;
					LIST_GCY[jt]-=DeltaC*NCONT[numc][1]*Sij*DCONTO[numc]/2;
					LIST_GCZ[jt]-=DeltaC*NCONT[numc][2]*Sij*DCONTO[numc]/2;

				}
		    }
	      
	     LIST_GCX[jt]/=VOLHALO[jt];      
	     LIST_GCY[jt]/=VOLHALO[jt];     
	     LIST_GCZ[jt]/=VOLHALO[jt];   
      		    
	   }


mingx=01.E+20;
maxgx=-1.E+20;

mingy=01.E+20;
maxgy=-1.E+20;

mingz=01.E+20;
maxgz=-1.E+20;



R GCTotx=0.;
R GCToty=0.;
R GCTotz=0.;

R volel;

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,volel) reduction(max:maxgx,maxgy,maxgz) reduction(min:mingx,mingy,mingz) reduction (+:GCTotx,GCToty,GCTotz)
	for(it=0;it<NB_SPH;it++)
	   {

		if(EDGE[it]==0){
		maxgx=(LIST_GCX[it]>maxgx)?LIST_GCX[it]:maxgx; 
		mingx=(LIST_GCX[it]<mingx)?LIST_GCX[it]:mingx; 

		maxgy=(LIST_GCY[it]>maxgy)?LIST_GCY[it]:maxgy; 
		mingy=(LIST_GCY[it]<mingy)?LIST_GCY[it]:mingy; 

		maxgz=(LIST_GCZ[it]>maxgz)?LIST_GCZ[it]:maxgz; 
		mingz=(LIST_GCZ[it]<mingz)?LIST_GCZ[it]:mingz; 
		}

		volel=coef1*LIST_V[it];

		GCTotx+=LIST_GCX[it]*volel;
		GCToty+=LIST_GCY[it]*volel;
		GCTotz+=LIST_GCZ[it]*volel;
	    
	      
	   }

//cout<<"GCX "<<GCTotx/(H_TOT*V_TOT*Z_TOT)<<endl;
//cout<<"GCY "<<GCToty/(H_TOT*V_TOT*Z_TOT)<<endl;
//cout<<"GCZ "<<GCTotz/(H_TOT*V_TOT*Z_TOT)<<endl;

}



void grad_C_P_Halo(int jt,R coef2, R Dw,  int ** CONT, R ** NCONT, R * LIST_R, R * LIST_CONW, R * DCONTO, R * LIST_GCX,R * LIST_GCY,R * LIST_GCZ,int * NBHALO,unsigned int **NOHALO,int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO,R * LIST_V,R coef1)  
{

R DeltaC;
int it,kt,numcor,numc;

int i,j;
R Sij;

R MoyR;
R coefpi=PI*coef2;



	       LIST_GCX[jt]=0.;
	       LIST_GCY[jt]=0.;
	       LIST_GCZ[jt]=0.;



	
	   

		    for(it=0;it<NBHALO[jt];it++){
			numcor=NOHALO[jt][it];

				for(kt=0;kt<NBCONTCO[numcor];kt++){ 

					numc=NOCONT[numcor][kt];
					i=CONT[numc][0];
					j=CONT[numc][1];

					MoyR=(LIST_R[i]+LIST_R[j])/2.;
					Sij=MoyR*MoyR*coefpi;

					DeltaC=(LIST_CONW[j]-LIST_CONW[i])/DCONTO[numc];

					LIST_GCX[jt]-=DeltaC*NCONT[numc][0]*Sij*DCONTO[numc]/2;
					LIST_GCY[jt]-=DeltaC*NCONT[numc][1]*Sij*DCONTO[numc]/2;
					LIST_GCZ[jt]-=DeltaC*NCONT[numc][2]*Sij*DCONTO[numc]/2;

				}
		    }
	      
	     LIST_GCX[jt]/=VOLHALO[jt];      
	     LIST_GCY[jt]/=VOLHALO[jt];     
	     LIST_GCZ[jt]/=VOLHALO[jt];   
      		    
	   




}


void grad_C_P(int jt,R coef2, R Dw,  int ** CONT, R ** NCONT, R * LIST_R, R * LIST_CONW, R * DCONTO, R * LIST_GCX,R * LIST_GCY,R * LIST_GCZ,int * NBHALO, unsigned int **NOHALO, int * NBCONTCO, unsigned int ** NOCONT, R * VOLHALO,R * LIST_V,R coef1)  
{

R DeltaC;
int it,kt,numc;

int i,j;
R Sij;

R MoyR;
R coefpi=PI*coef2;



	       LIST_GCX[jt]=0.;
	       LIST_GCY[jt]=0.;
	       LIST_GCZ[jt]=0.;



	
	   

		    

				for(kt=0;kt<NBCONTCO[jt];kt++){ 

					numc=NOCONT[jt][kt];
					i=CONT[numc][0];
					j=CONT[numc][1];

					MoyR=(LIST_R[i]+LIST_R[j])/2.;
					Sij=MoyR*MoyR*coefpi;

					DeltaC=(LIST_CONW[j]-LIST_CONW[i])/DCONTO[numc];

					LIST_GCX[jt]-=DeltaC*NCONT[numc][0]*Sij*DCONTO[numc]/2;
					LIST_GCY[jt]-=DeltaC*NCONT[numc][1]*Sij*DCONTO[numc]/2;
					LIST_GCZ[jt]-=DeltaC*NCONT[numc][2]*Sij*DCONTO[numc]/2;

				}
		    
	      
	     LIST_GCX[jt]/=(coef1*LIST_V[jt]);      
	     LIST_GCY[jt]/=(coef1*LIST_V[jt]);     
	     LIST_GCZ[jt]/=(coef1*LIST_V[jt]);   
      		    
	   



}

R crit_time(R coord, R NBCO, R lambda, R cp, R dens, R * DCONTO, R Cs){
  R Length=1.E+09;
  for (int it=0; it<NBCO; it++) if(DCONTO[it]<Length) Length=DCONTO[it];
  return Cs*(Length*Length/coord)*(dens*cp/lambda);
}

R crit_time2(R coord, R NBCO, R lambda1, R lambda2, R cp1, R cp2, R dens1, R dens2, R * DCONTO, R Cs){
  R Length=1.E+09;
  for (int it=0; it<NBCO; it++) if(DCONTO[it]<Length) Length=DCONTO[it];

  if((dens2*cp2/lambda2)>=(dens1*cp1/lambda1)) 
    return Cs*(Length*Length/coord)*(dens1*cp1/lambda1);
  else
    return Cs*(Length*Length/coord)*(dens2*cp2/lambda2);

  
}

R crit_time_couplage(R coord, R NBCO, R lambda, R cp, R dens, R * DCONTO, R Cs, R Dw){
  R Length=1.E+09, dtt, dth;
  for (int it=0; it<NBCO; it++) if(DCONTO[it]<Length) Length=DCONTO[it];
   dtt=Cs*(Length*Length/coord)*(dens*cp/lambda);
   dth=Cs*(Length*Length/coord)*(1./Dw);

   if(dtt<dth) return dtt;
    else return dth;
}

