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

#include "init.h"
#include "omp.h"

void init_sphC(int NB_SPH, R Cinit, R Cimp, R * LIST_CONW,R * Cp_eff, R * LIST_M, R * LIST_V, R p, R Cps, R Cpw, R rhos, R Mw, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4, bool * EDGE5, bool * EDGE6, R ** NCONT,int NBCO, int ** CONT, R* LIST_X, R* LIST_Y, R* LIST_Z, R * DCONTO, R * DeltaC, int dir, bool btype)
{

R rho_eff;
R rhos_eff=rhos*(1-p);

// INITIATION DE LA CONCENTRATION
for (int it=0;it<NB_SPH;it++)
    {
     LIST_CONW[it] = Cinit;
     DeltaC[it]=0.;
  
     if((btype==0)&&(dir==3)&&(EDGE6[it])) LIST_CONW[it] = Cimp;
     if((btype==0)&&(dir==2)&&(EDGE5[it])) LIST_CONW[it] = Cimp;  
     if((btype==0)&&(dir==1)&&(EDGE4[it])) LIST_CONW[it] = Cimp;  

       rho_eff=rhos_eff+Mw*LIST_CONW[it];
       Cp_eff[it]=(Cps*rhos_eff+Cpw*(rho_eff-rhos_eff))/rho_eff;
       LIST_M[it]=LIST_V[it]*rho_eff;
 }


int i,j;
for(int it=0;it<NBCO;it++)
   {   

    i=CONT[it][0];
    j=CONT[it][1]; 

    DCONTO[it]=sqrt((LIST_X[j]- LIST_X[i])*(LIST_X[j]- LIST_X[i])+(LIST_Y[j]- LIST_Y[i])*(LIST_Y[j]- LIST_Y[i])+(LIST_Z[j]- LIST_Z[i])*(LIST_Z[j]- LIST_Z[i]));

    NCONT[it][0]= (LIST_X[i]- LIST_X[j])/DCONTO[it]; 
    NCONT[it][1]= (LIST_Y[i]- LIST_Y[j])/DCONTO[it]; 
    NCONT[it][2]= (LIST_Z[i]- LIST_Z[j])/DCONTO[it]; 

   }

	
}

void init_sphT(int NB_SPH, R Tinit, R Timp, R * LIST_TEMP, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,
R ** NCONT,int NBCO, int ** CONT, R* LIST_X, R* LIST_Y, R* LIST_Z, R * DCONTO, R * DeltaT, int dir, bool btype)
{

// INITIATION DE LA TEMPERATURE ET FLUX

for (int it=0;it<NB_SPH;it++)
    {
     LIST_TEMP[it] = Tinit;
     DeltaT[it]=0.;
  
     if((btype==0)&&(dir==3)&&(EDGE6[it])) LIST_TEMP[it] = Timp;
     if((btype==0)&&(dir==2)&&(EDGE5[it])) LIST_TEMP[it] = Timp;  
     if((btype==0)&&(dir==1)&&(EDGE4[it])) LIST_TEMP[it] = Timp;  
 }


int i,j;
for(int it=0;it<NBCO;it++)
   {   

    i=CONT[it][0];
    j=CONT[it][1]; 

    DCONTO[it]=sqrt((LIST_X[j]- LIST_X[i])*(LIST_X[j]- LIST_X[i])+(LIST_Y[j]- LIST_Y[i])*(LIST_Y[j]- LIST_Y[i])+(LIST_Z[j]- LIST_Z[i])*(LIST_Z[j]- LIST_Z[i]));

    NCONT[it][0]= (LIST_X[i]- LIST_X[j])/DCONTO[it]; 
    NCONT[it][1]= (LIST_Y[i]- LIST_Y[j])/DCONTO[it]; 
    NCONT[it][2]= (LIST_Z[i]- LIST_Z[j])/DCONTO[it]; 

   }

	
}




void init_sphp_cis(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,bool * EDGEC1, bool * EDGEC2, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT
    for(int it=0;it<NMAXCONT;it++){
		TYPCO[it]=1;		
	}
	
//cout<<"Nombre de spheres:"<<NB_SPH<<endl;

PR_SPH=0.;

int nb1=0;
int nb2=0;
int nb3=0;
int nb4=0;
int nb5=0;
int nb6=0;

int nbc1=0;
int nbc2=0;

for(int it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

	// DISTRIBUTION INITIALE DE SPHERES
	for(int it=0;it<NB_SPH;it++){

      LIST_IND[it]=1.;

      PR_SPH+=(4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it];
   
         LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;    
         
      NBCONTCO[it]= 0;
            
      if(fabs(LIST_X[it]-LIST_R[it]-H_POS)<epsi) { EDGE1[it]=1;  nb1++; }
      if(fabs(LIST_Y[it]-LIST_R[it]-V_POS)<epsi) { EDGE2[it]=1;  nb2++; }
      if(fabs(LIST_Z[it]-LIST_R[it]-Z_POS)<epsi) { EDGE3[it]=1;  nb3++; }
      if(fabs(LIST_X[it]+LIST_R[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1;  nb4++; }
      if(fabs(LIST_Y[it]+LIST_R[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;  nb5++; }
      if(fabs(LIST_Z[it]+LIST_R[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;  nb6++; }          	
     
     R dist=(LIST_X[it]-(H_TOT/2.+H_POS))*(LIST_X[it]-(H_TOT/2.+H_POS));
     dist=dist+(LIST_Y[it]-(V_TOT/2.+V_POS))*(LIST_Y[it]-(V_TOT/2.+V_POS));
     dist=sqrt(dist)+LIST_R[it];
     if((fabs(dist-(H_TOT/2.+H_POS))<5.*epsi)&&(LIST_Z[it]<=Z_TOT/2.+Z_POS)) { EDGEC1[it]=1; nbc1++;    }     
     if((fabs(dist-(H_TOT/2.+H_POS))<5.*epsi)&&(LIST_Z[it]>Z_TOT/2.+Z_POS)) { EDGEC2[it]=1;  nbc2++;   }     
     
    		
		int NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		int NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		int NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		
		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;
		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*1.;  
		LIST_M[it]  = dens*LIST_V[it];    
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

N1=R(nb1);
N2=R(nb2);
N3=R(nb3);
N4=R(nb4);
N5=R(nb5);
N6=R(nb6);


      cout<<"nbc:"<<nbc1<<", "<<nbc2<<endl;
      cout<<"nb36:"<<nb3<<", "<<nb6<<endl;
}

void init_sizeh(int & vecsizeh, int & vecsizexh,int & vecsizeyh, int & vecsizezh, vector< vector<int>  > & coulh,vector< vector<int>  > & numch, R H_TOT, R V_TOT, R Z_TOT, R & RHALO, int SZHALO, R R_SPH)
{
RHALO=SZHALO*R_SPH;

vecsizexh=floor(H_TOT/(RHALO))-1;
vecsizeyh=floor(V_TOT/(RHALO))-1;
vecsizezh=floor(Z_TOT/(RHALO))-1;

if(vecsizexh<=0) vecsizexh=1;
if(vecsizeyh<=0) vecsizeyh=1;
if(vecsizezh<=0) vecsizezh=1;

/*
	vecsizexh=1;
	vecsizeyh=1;
	vecsizezh=1;
*/	
vecsizeh=vecsizexh*vecsizeyh*vecsizezh;

for (int i = 0; i < vecsizeh; i++) {
    coulh.push_back(vector<int>()); // Add an empty row
}
for (int ii = 0; ii < vecsizeh; ii++) {
	numch.push_back(vector<int>()); // Add an empty row
}

for (int ii = 0; ii < vecsizeh; ii++) {

	int NX  = ii%vecsizexh;
	int NY  = (ii%(vecsizexh*vecsizeyh))/vecsizexh;
	int NZ  = ii/(vecsizexh*vecsizeyh);
	
	int NXP1= (NX+vecsizexh+1)%vecsizexh;
	int NXM1= (NX+vecsizexh-1)%vecsizexh;
	int NYP1= (NY+vecsizeyh+1)%vecsizeyh;	
	int NYM1= (NY+vecsizeyh-1)%vecsizeyh;
        int NZP1= (NZ+vecsizezh+1)%vecsizezh;	
	int NZM1= (NZ+vecsizezh-1)%vecsizezh;		

	
                if(NZM1==NZ-1){


		        if(NYM1==NY-1){

				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZM1);
				}

                                {
				numch[ii].push_back(NX+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZM1);
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZM1);			     
				}


		        }

                        {

				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NY+vecsizexh*vecsizeyh*NZM1);
				}

                                {
				numch[ii].push_back(NX+vecsizexh*NY+vecsizexh*vecsizeyh*NZM1);	
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NY+vecsizexh*vecsizeyh*NZM1);			     
				}

		        }

                        if(NYP1==NY+1){
	
				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZM1);
				}

                                {
				numch[ii].push_back(NX+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZM1);
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZM1);			     
				}
	     
		        }

                }
                
                {


		        if(NYM1==NY-1){

				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZ);
				}

                                {
				numch[ii].push_back(NX+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZ);
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZ);			     
				}


		        }

                        {

				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NY+vecsizexh*vecsizeyh*NZ);
				}

                                {
				numch[ii].push_back(ii);
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NY+vecsizexh*vecsizeyh*NZ);			     
				}

		        }

                        if(NYP1==NY+1){
	
				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZ);
				}

                                {
				numch[ii].push_back(NX+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZ);
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZ);			     
				}
	     
		        }

                }		

                if(NZP1==NZ+1){


		        if(NYM1==NY-1){

				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZP1);	
				}

                                {
				numch[ii].push_back(NX+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZP1);
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NYM1+vecsizexh*vecsizeyh*NZP1);			     
				}


		        }

                        {

				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NY+vecsizexh*vecsizeyh*NZP1);
				}

                                {
				numch[ii].push_back(NX+vecsizexh*NY+vecsizexh*vecsizeyh*NZP1);
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NY+vecsizexh*vecsizeyh*NZP1);			     
				}

		        }

                        if(NYP1==NY+1){
	
				if(NXM1==NX-1){
				numch[ii].push_back(NXM1+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZP1);
				}

                                {
				numch[ii].push_back(NX+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZP1);
				}

                                if(NXP1==NX+1){
				numch[ii].push_back(NXP1+vecsizexh*NYP1+vecsizexh*vecsizeyh*NZP1); 			     
				}
	     
		        }

                }	

//cout<<"ii:"<<ii<<"  "<<numch[ii].size()<<endl;

 }	 
  
 cout<<"NXYZ_halo:"<<vecsizexh<<", "<<vecsizeyh<<", "<<vecsizezh<<endl;
 	

}


void init_paroi(int NB_PAR, R densp, R epaisp, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ, R ** LIST_PVX, R ** LIST_PVY, R ** LIST_PVZ, R ** LIST_PAX, R ** LIST_PAY, R ** LIST_PAZ, R ** LIST_PN, R * LIST_PM)
{

R distij,distil;
R det1,det2,det3;
R fors,norms;
int it;

	// DISTRIBUTION INITIALE DE SPHERES
	for(it=0;it<NB_PAR;it++){
		
		for(int jt=0;jt<4;jt++){	
		
		LIST_PVX[it][jt]=0.;
		LIST_PVY[it][jt]=0.;		
		LIST_PVZ[it][jt]=0.;		

		LIST_PAX[it][jt]=0.;
		LIST_PAY[it][jt]=0.;		
		LIST_PAZ[it][jt]=0.;					
	        }
		
		distij=(LIST_PX[it][0]-LIST_PX[it][1])*(LIST_PX[it][0]-LIST_PX[it][1]);
		distij+=(LIST_PY[it][0]-LIST_PY[it][1])*(LIST_PY[it][0]-LIST_PY[it][1]);		
		distij+=(LIST_PZ[it][0]-LIST_PZ[it][1])*(LIST_PZ[it][0]-LIST_PZ[it][1]);			
		distij=sqrt(distij);
		
		distil=(LIST_PX[it][0]-LIST_PX[it][3])*(LIST_PX[it][0]-LIST_PX[it][3]);
		distil+=(LIST_PY[it][0]-LIST_PY[it][3])*(LIST_PY[it][0]-LIST_PY[it][3]);		
		distil+=(LIST_PZ[it][0]-LIST_PZ[it][3])*(LIST_PZ[it][0]-LIST_PZ[it][3]);			
		distil=sqrt(distil);	
		
		LIST_PM[it]=distij*distil*epaisp*densp;
	
		det1=((LIST_PY[it][1]-LIST_PY[it][0])*(LIST_PZ[it][3]-LIST_PZ[it][0])-(LIST_PY[it][3]-LIST_PY[it][0])*(LIST_PZ[it][1]-LIST_PZ[it][0]));
        det2=((LIST_PX[it][3]-LIST_PX[it][0])*(LIST_PZ[it][1]-LIST_PZ[it][0])-(LIST_PX[it][1]-LIST_PX[it][0])*(LIST_PZ[it][3]-LIST_PZ[it][0]));
        det3=((LIST_PX[it][1]-LIST_PX[it][0])*(LIST_PY[it][3]-LIST_PY[it][0])-(LIST_PX[it][3]-LIST_PX[it][0])*(LIST_PY[it][1]-LIST_PY[it][0]));  

        LIST_PN[it][0]=det1/sqrt(det1*det1+det2*det2+det3*det3);
        LIST_PN[it][1]=det2/sqrt(det1*det1+det2*det2+det3*det3);      
        LIST_PN[it][2]=det3/sqrt(det1*det1+det2*det2+det3*det3);  
        
	   // Zero numerique
	   
	   if(fabs(LIST_PN[it][0])<1e-18) LIST_PN[it][0]=0.;
	   if(fabs(LIST_PN[it][1])<1e-18) LIST_PN[it][1]=0.;
	   if(fabs(LIST_PN[it][2])<1e-18) LIST_PN[it][2]=0.;
	   
	   // Tangente s
	   
	   if((LIST_PN[it][0]*LIST_PN[it][1]==0.)&&(LIST_PN[it][1]*LIST_PN[it][2]==0.)&&(LIST_PN[it][2]*LIST_PN[it][0]==0.)){
		  
		  if(LIST_PN[it][0]!=0.){
		   LIST_PN[it][3] = 0.;
		   LIST_PN[it][4] = 0.;
		   LIST_PN[it][5] = 1.;						
		  }else if(LIST_PN[it][1]!=0.){
		   LIST_PN[it][3] = 1.;
		   LIST_PN[it][4] = 0.;
		   LIST_PN[it][5] = 0.;									  
		  }else if(LIST_PN[it][2]!=0.){
		   LIST_PN[it][3] = 0.;
		   LIST_PN[it][4] = 1.;
		   LIST_PN[it][5] = 0.;									  
		  }
		   
	   }else{
			
			if(LIST_PN[it][0]==0.){
			
			fors  = -(LIST_PN[it][1])/LIST_PN[it][2];
			norms = sqrt(fors*fors+2.);

			LIST_PN[it][3] = 1./norms;
			LIST_PN[it][4] = 1./norms;
			LIST_PN[it][5] = fors/norms;					
			
			}else{
				
		   fors  = -(LIST_PN[it][1]+LIST_PN[it][2])/LIST_PN[it][0];
		   norms = sqrt(fors*fors+2.);

		   LIST_PN[it][3] = fors/norms;
		   LIST_PN[it][4] = 1./norms;
		   LIST_PN[it][5] = 1./norms;							
			
			}
									
			
	   }
	   
	   
	   // Tangente t				   
	   
		LIST_PN[it][6] = LIST_PN[it][1]*LIST_PN[it][5] - LIST_PN[it][2]*LIST_PN[it][4];
		LIST_PN[it][7] = LIST_PN[it][2]*LIST_PN[it][3] - LIST_PN[it][0]*LIST_PN[it][5];
		LIST_PN[it][8] = LIST_PN[it][0]*LIST_PN[it][4] - LIST_PN[it][1]*LIST_PN[it][3];          
       
				
	}
	
	
} 


void init_sph(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT


for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}


PR_SPH=0.;
R coef;
int NX,NY,NZ;
int it;

N1=0.;
N2=0.;
N3=0.;
N4=0.;
N5=0.;
N6=0.;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}
     // DISTRIBUTION INITIALE DE SPHERES

  for(it=0;it<NB_SPH;it++){
      coef=1.;
      
      LIST_B[it]=0.;
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
         
      NBCONTCO[it]= 0;
      		
      if(fabs(LIST_X[it]-H_POS)<epsi) { EDGE1[it]=1;    EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Y[it]-V_POS)<epsi) { EDGE2[it]=1;   EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Z[it]-Z_POS)<epsi) { EDGE3[it]=1;   EDGE[it]=1; coef/=2.;}
      if(fabs(LIST_X[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1;   EDGE[it]=1; coef/=2.;}	
      if(fabs(LIST_Y[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;   EDGE[it]=1; coef/=2.;} 
      if(fabs(LIST_Z[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;    EDGE[it]=1;coef/=2.; } 

      if(EDGE1[it]) N1+=coef;
      if(EDGE2[it]) N2+=coef;
      if(EDGE3[it]) N3+=coef;
      if(EDGE4[it]) N4+=coef;
      if(EDGE5[it]) N5+=coef;
      if(EDGE6[it]) N6+=coef;

      PR_SPH+=(4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;
    		
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		
		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;
                LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];
		
		LIST_XA[it] = LIST_X[it];
		LIST_YA[it] = LIST_Y[it];
		LIST_ZA[it] = LIST_Z[it];	

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;
		
		LIST_TXA[it] = 0.;
		LIST_TYA[it] = 0.;
		LIST_TZA[it] = 0.;		
				
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;
		LIST_M[it]  = dens*LIST_V[it];		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

		LIST_IND[it]  = coef;

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

}


void init_sphp(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS,int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT

for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}

int it;
int NX,NY,NZ;

PR_SPH=0.;

int nb1=0;
int nb2=0;
int nb3=0;
int nb4=0;
int nb5=0;
int nb6=0;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

	// DISTRIBUTION INITIALE DE SPHERES

	for(it=0;it<NB_SPH;it++){

	LIST_IND[it]=1.;

      
      PR_SPH+=((4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]);
   
      LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
               
      NBCONTCO[it]= 0;
            
      if(fabs(LIST_X[it]-LIST_R[it]-H_POS)<epsi) { EDGE1[it]=1;  nb1++; }
      if(fabs(LIST_Y[it]-LIST_R[it]-V_POS)<epsi) { EDGE2[it]=1;  nb2++; }
      if(fabs(LIST_Z[it]-LIST_R[it]-Z_POS)<epsi) { EDGE3[it]=1;  nb3++; }
      if(fabs(LIST_X[it]+LIST_R[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1;  nb4++; }
      if(fabs(LIST_Y[it]+LIST_R[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;  nb5++; }
      if(fabs(LIST_Z[it]+LIST_R[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;  nb6++; }          	
  
    		
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		


		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;

		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_XA[it] = LIST_X[it];
		LIST_YA[it] = LIST_Y[it];
		LIST_ZA[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;

		LIST_TXA[it] = 0.;
		LIST_TYA[it] = 0.;
		LIST_TZA[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*1.;	
		LIST_M[it]  = dens*LIST_V[it];		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);


N1=R(nb1);
N2=R(nb2);
N3=R(nb3);
N4=R(nb4);
N5=R(nb5);
N6=R(nb6);



}

void init_sph_ind(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT
    for(int it=0;it<NMAXCONT;it++){
		TYPCO[it]=1;		
	}
	
//cout<<"Nombre de spheres:"<<NB_SPH<<endl;

PR_SPH=0.;

N1=0.;
N2=0.;
N3=0.;
N4=0.;
N5=0.;
N6=0.;

R coef;

for(int it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

	// DISTRIBUTION INITIALE DE SPHERES
	for(int it=0;it<NB_SPH;it++){

 	coef=1.;
      
      PR_SPH+=(4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;
   
      LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
         
      NBCONTCO[it]= 0;
            
      if(fabs(LIST_X[it]-H_POS)<epsi) { EDGE1[it]=1;    EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Y[it]-V_POS)<epsi) { EDGE2[it]=1;   EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Z[it]-Z_POS)<epsi) { EDGE3[it]=1;   EDGE[it]=1; coef/=2.;}
      if(fabs(LIST_X[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1;   EDGE[it]=1; coef/=2.;}	
      if(fabs(LIST_Y[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;   EDGE[it]=1; coef/=2.;} 
      if(fabs(LIST_Z[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;    EDGE[it]=1;coef/=2.; } 

      if(EDGE1[it]) N1+=coef;
      if(EDGE2[it]) N2+=coef;
      if(EDGE3[it]) N3+=coef;
      if(EDGE4[it]) N4+=coef;
      if(EDGE5[it]) N5+=coef;
      if(EDGE6[it]) N6+=coef;        	
  
    		
		int NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		int NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		int NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		
                if(it!=NB_SPH-1){
		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;

		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   
                }

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_XA[it] = LIST_X[it];
		LIST_YA[it] = LIST_Y[it];
		LIST_ZA[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;

		LIST_TXA[it] = 0.;
		LIST_TYA[it] = 0.;
		LIST_TZA[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;
		LIST_M[it]  = dens*LIST_V[it];	
		if(it==NB_SPH-1) {
		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it];
		LIST_M[it]  = dens*LIST_V[it]*1.3;		
		}
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	
		LIST_IND[it]=coef;

	}
	

	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

}

void init_sphp_ind(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT
    for(int it=0;it<NMAXCONT;it++){
		TYPCO[it]=1;		
	}
	
//cout<<"Nombre de spheres:"<<NB_SPH<<endl;

PR_SPH=0.;
int nb1=0;
int nb2=0;
int nb3=0;
int nb4=0;
int nb5=0;
int nb6=0;

for(int it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

	// DISTRIBUTION INITIALE DE SPHERES
	for(int it=0;it<NB_SPH;it++){
     
		LIST_IND[it]=1.;
 
      PR_SPH+=(4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it];
   
         LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
         
      NBCONTCO[it]= 0;
            
      if(fabs(LIST_X[it]-LIST_R[it]-H_POS)<epsi) { EDGE1[it]=1; nb1++;   }
      if(fabs(LIST_Y[it]-LIST_R[it]-V_POS)<epsi) { EDGE2[it]=1; nb2++;}
      if(fabs(LIST_Z[it]-LIST_R[it]-Z_POS)<epsi) { EDGE3[it]=1; nb3++; }
      if(fabs(LIST_X[it]+LIST_R[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1; nb4++; 	 }
      if(fabs(LIST_Y[it]+LIST_R[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;  nb5++; }
      if(fabs(LIST_Z[it]+LIST_R[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;  nb6++;  }          	
  
    		
		int NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		int NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		int NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		
                if(it!=NB_SPH-1){
		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;
		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   
                }

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_XA[it] = LIST_X[it];
		LIST_YA[it] = LIST_Y[it];
		LIST_ZA[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;

		LIST_TXA[it] = 0.;
		LIST_TYA[it] = 0.;
		LIST_TZA[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*1.;	
		LIST_M[it]  = dens*LIST_V[it]*1.;		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	
                 if(it==NB_SPH-1) {
		LIST_M[it]  = dens*LIST_V[it]*1.3;		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	
                 }


	}
	
N1=R(nb1);
N2=R(nb2);
N3=R(nb3);
N4=R(nb4);
N5=R(nb5);
N6=R(nb6);

	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);
}

void init_sph_halo(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul,int * nocoul,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT

for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}

PR_SPH=0.;

R coef;
int NX,NY,NZ;
int it;

N1=0.;
N2=0.;
N3=0.;
N4=0.;
N5=0.;
N6=0.;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}
	// DISTRIBUTION INITIALE DE SPHERES
	for(it=0;it<NB_SPH;it++){
        
      LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
         
      NBCONTCO[it]= 0;
      coef=1.;
      
      if(fabs(LIST_X[it]-H_POS)<epsi) { EDGE1[it]=1;    EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Y[it]-V_POS)<epsi) { EDGE2[it]=1;   EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Z[it]-Z_POS)<epsi) { EDGE3[it]=1;   EDGE[it]=1; coef/=2.;}
      if(fabs(LIST_X[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1;   EDGE[it]=1; coef/=2.;}	
      if(fabs(LIST_Y[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;   EDGE[it]=1; coef/=2.;} 
      if(fabs(LIST_Z[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;    EDGE[it]=1;coef/=2.; }    

      if(EDGE1[it]) N1+=coef;
      if(EDGE2[it]) N2+=coef;
      if(EDGE3[it]) N3+=coef;
      if(EDGE4[it]) N4+=coef;
      if(EDGE5[it]) N5+=coef;
      if(EDGE6[it]) N6+=coef;  	

	LIST_IND[it]=coef;
  
      PR_SPH+=((4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it])*coef;
    		
       // grille détection		
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;	
	

		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;
		
		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

       // grille halo	
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizexh/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizeyh/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizezh/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizexh) NX=vecsizexh-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizeyh) NY=vecsizeyh-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizezh) NZ=vecsizezh-1;		

		coulh[NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh].push_back(it);
		LIST_H[it] = NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;
		LIST_M[it]  = dens*LIST_V[it];		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

}

void init_sph_halo2(bool * TYPCO, int NMAXCONT, R Pi, R dens1,R dens2, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul,int * nocoul,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, bool * LIST_P, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT

for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}

PR_SPH=0.;

R coef;
int NX,NY,NZ;
int it;

N1=0.;
N2=0.;
N3=0.;
N4=0.;
N5=0.;
N6=0.;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

	// DISTRIBUTION INITIALE DE SPHERES
	for(it=0;it<NB_SPH;it++){
        
      LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
         
      NBCONTCO[it]= 0;
      coef=1.;
      
      if(fabs(LIST_X[it]-H_POS)<epsi) { EDGE1[it]=1;    EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Y[it]-V_POS)<epsi) { EDGE2[it]=1;   EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Z[it]-Z_POS)<epsi) { EDGE3[it]=1;   EDGE[it]=1; coef/=2.;}
      if(fabs(LIST_X[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1;   EDGE[it]=1; coef/=2.;}	
      if(fabs(LIST_Y[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;   EDGE[it]=1; coef/=2.;} 
      if(fabs(LIST_Z[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;    EDGE[it]=1;coef/=2.; }    	

      if(EDGE1[it]) N1+=coef;
      if(EDGE2[it]) N2+=coef;
      if(EDGE3[it]) N3+=coef;
      if(EDGE4[it]) N4+=coef;
      if(EDGE5[it]) N5+=coef;
      if(EDGE6[it]) N6+=coef;    

      PR_SPH+=((4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it])*coef;
    		
      LIST_IND[it]=coef;

       // grille détection		
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;	
	
		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;
	
		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

       // grille halo	
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizexh/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizeyh/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizezh/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizexh) NX=vecsizexh-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizeyh) NY=vecsizeyh-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizezh) NZ=vecsizezh-1;		

		coulh[NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh].push_back(it);
		LIST_H[it] = NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;
		if(LIST_P[it]==0){	
		LIST_M[it]  = dens1*LIST_V[it];		
		}
		else{
		LIST_M[it]  = dens2*LIST_V[it];			
		}			
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

}


void init_sphp_halo(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul,int * nocoul,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA,R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT

for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}

int NX,NY,NZ;
int it;

PR_SPH=0.;

int nb1=0;
int nb2=0;
int nb3=0;
int nb4=0;
int nb5=0;
int nb6=0;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}


	// DISTRIBUTION INITIALE DE SPHERES
	for(it=0;it<NB_SPH;it++){
          LIST_IND[it]=1.;
  
      PR_SPH+=((4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]);
   
      LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
         
      NBCONTCO[it]= 0;
            
      if(fabs(LIST_X[it]-LIST_R[it]-H_POS)<epsi) { EDGE1[it]=1; nb1++;   }
      if(fabs(LIST_Y[it]-LIST_R[it]-V_POS)<epsi) { EDGE2[it]=1; nb2++;}
      if(fabs(LIST_Z[it]-LIST_R[it]-Z_POS)<epsi) { EDGE3[it]=1; nb3++; }
      if(fabs(LIST_X[it]+LIST_R[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1; nb4++; 	 }
      if(fabs(LIST_Y[it]+LIST_R[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;  nb5++; }
      if(fabs(LIST_Z[it]+LIST_R[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;  nb6++;  }          	
  
    		
       // grille détection		
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		

		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;

		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

       // grille halo	
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizexh/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizeyh/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizezh/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizexh) NX=vecsizexh-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizeyh) NY=vecsizeyh-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizezh) NZ=vecsizezh-1;		

		coulh[NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh].push_back(it);
		LIST_H[it] = NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*1.;
		LIST_M[it]  = dens*LIST_V[it];		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

N1=R(nb1);
N2=R(nb2);
N3=R(nb3);
N4=R(nb4);
N5=R(nb5);
N6=R(nb6);

}


void init_sphp_halo_car(R EE, bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA,R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT

for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}


PR_SPH=0.;

int NX,NY,NZ;
int it;

int nb1=0;
int nb2=0;
int nb3=0;
int nb4=0;
int nb5=0;
int nb6=0;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

	// DISTRIBUTION INITIALE DE SPHERES
	for(it=0;it<NB_SPH;it++){
             LIST_IND[it]=1.;  
      PR_SPH+=((4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]);
   
      LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
               
      NBCONTCO[it]= 0;
            
      if(fabs(LIST_X[it]-EE-LIST_R[it])<epsi) { EDGE1[it]=1; nb1++;}
      if(fabs(LIST_Y[it]-LIST_R[it])<epsi) { EDGE2[it]=1; nb2++;}
      if(fabs(LIST_Z[it]-LIST_R[it])<epsi) { EDGE3[it]=1; nb3++; }
      if(fabs(LIST_X[it]+LIST_R[it]-H_TOT+EE)<epsi) { EDGE4[it]=1;nb4++; }
      if(fabs(LIST_Y[it]+LIST_R[it]-V_TOT)<epsi) { EDGE5[it]=1; nb5++;  }
      if(fabs(LIST_Z[it]+LIST_R[it]-Z_TOT)<epsi) { EDGE6[it]=1; nb6++;  }          	
  
    		
       // grille détection		
		NX= ((int) (LIST_X[it]*vecsizex/H_TOT));
		NY= ((int) (LIST_Y[it]*vecsizey/V_TOT));
		NZ= ((int) (LIST_Z[it]*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		
		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;
		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

       // grille halo	
		NX= ((int) (LIST_X[it]*vecsizexh/H_TOT));
		NY= ((int) (LIST_Y[it]*vecsizeyh/V_TOT));
		NZ= ((int) (LIST_Z[it]*vecsizezh/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizexh) NX=vecsizexh-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizeyh) NY=vecsizeyh-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizezh) NZ=vecsizezh-1;		
		
		coulh[NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh].push_back(it);
		LIST_H[it] = NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*1.;	
		LIST_M[it]  = dens*LIST_V[it];		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

N1=R(nb1);
N2=R(nb2);
N3=R(nb3);
N4=R(nb4);
N5=R(nb5);
N6=R(nb6);

}

void init_sph2(bool * TYPCO, int NMAXCONT, R Pi, R dens1, R dens2, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul,R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6, int * NBCONTCO, bool * LIST_P, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{


    // CONTACTS COHESIFS INITIALEMENT

for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}

int NX,NY,NZ;	
int it;
R coef;

PR_SPH=0.;


N1=0.;
N2=0.;
N3=0.;
N4=0.;
N5=0.;
N6=0.;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

	// DISTRIBUTION INITIALE DE SPHERES
	for(it=0;it<NB_SPH;it++){
      
      LIST_B[it]=0.;
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
         
      NBCONTCO[it]= 0;
      coef=1.;
      
      if(fabs(LIST_X[it]-H_POS)<epsi) { EDGE1[it]=1;    EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Y[it]-V_POS)<epsi) { EDGE2[it]=1;   EDGE[it]=1; coef/=2.; }
      if(fabs(LIST_Z[it]-Z_POS)<epsi) { EDGE3[it]=1;   EDGE[it]=1; coef/=2.;}
      if(fabs(LIST_X[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1;   EDGE[it]=1; coef/=2.;}	
      if(fabs(LIST_Y[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;   EDGE[it]=1; coef/=2.;} 
      if(fabs(LIST_Z[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;    EDGE[it]=1;coef/=2.; }    	

      if(EDGE1[it]) N1+=coef;
      if(EDGE2[it]) N2+=coef;
      if(EDGE3[it]) N3+=coef;
      if(EDGE4[it]) N4+=coef;
      if(EDGE5[it]) N5+=coef;
      if(EDGE6[it]) N6+=coef;  
	
  
      PR_SPH+=((4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it])*coef;
      LIST_IND[it]=coef;
    		
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		
		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;

		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_XA[it] = LIST_X[it];
		LIST_YA[it] = LIST_Y[it];
		LIST_ZA[it] = LIST_Z[it];
		
		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;

		LIST_TXA[it] = 0.;
		LIST_TYA[it] = 0.;
		LIST_TZA[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;

		if(LIST_P[it]==0){	
		LIST_M[it]  = dens1*LIST_V[it];		
		}
		else{
		LIST_M[it]  = dens2*LIST_V[it];			
		}	
	
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

}


void scale_paroi(int NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ, R * LIST_PM, int NBCOP, R * DCONTP,  R & epaisp, R H_TOT2, R lmacro)
{
	
R LG=H_TOT2;		
	
	for (int it=0;it<NB_PAR;it++){
	    for (int jt=0;jt<4;jt++){		
			
			LIST_PX[it][jt]=(LIST_PX[it][jt]/LG)*lmacro;
		    LIST_PY[it][jt]=(LIST_PY[it][jt]/LG)*lmacro;		
		    LIST_PZ[it][jt]=(LIST_PZ[it][jt]/LG)*lmacro;	
		
	    }
	    
	  LIST_PM[it]= (lmacro*lmacro*lmacro)*LIST_PM[it]/(LG*LG*LG) ;
	}
	
epaisp=(epaisp/LG)*lmacro;
	
	for (int it=0;it<NBCOP;it++){
		DCONTP[it]=DCONTP[it]*lmacro/LG;
	}	
	
}

void scale_sph(R Pi, R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R & RMAX, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_M, R * LIST_I, R & epsi, R lmacro, R dens, R * LIST_V)
{

R LG=H_TOT;	
R_SPH=0.;
RMAX=0.;  

R HG=H_POS;
R VG=V_POS;
R ZG=Z_POS;

int it;

for (it=0;it<NB_SPH;it++){

LIST_X[it]=(LIST_X[it]/LG)*lmacro;
LIST_Y[it]=(LIST_Y[it]/LG)*lmacro;
LIST_Z[it]=(LIST_Z[it]/LG)*lmacro;

LIST_XO[it] = LIST_X[it];
LIST_YO[it] = LIST_Y[it];
LIST_ZO[it] = LIST_Z[it];
	
LIST_XA[it] = LIST_X[it];
LIST_YA[it] = LIST_Y[it];
LIST_ZA[it] = LIST_Z[it];
		
LIST_R[it]=(LIST_R[it]/LG)*lmacro;

LIST_V[it]  = lmacro*lmacro*lmacro*LIST_V[it]/(LG*LG*LG);
LIST_M[it]  = dens*LIST_V[it];		
LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

R_SPH+=LIST_R[it];
RMAX=max(RMAX,LIST_R[it]);

}

H_TOT=(H_TOT/LG)*lmacro;
V_TOT=(V_TOT/LG)*lmacro;
Z_TOT=(Z_TOT/LG)*lmacro;

H_POS=(H_POS/LG)*lmacro;
V_POS=(V_POS/LG)*lmacro;
Z_POS=(Z_POS/LG)*lmacro;

R_SPH/=(NB_SPH);
epsi=(epsi/LG)*lmacro;
}

void scale_sphp(R Pi, R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R & RMAX, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_M, R * LIST_I, R & epsi, R lmacro, R dens, R * LIST_V)
{

R LG=H_TOT;	

R_SPH=0.;
RMAX=0.;  

R HG=H_POS;
R VG=V_POS;
R ZG=Z_POS;

int it;
for (it=0;it<NB_SPH;it++){

LIST_X[it]=(LIST_X[it]/LG)*lmacro;
LIST_Y[it]=(LIST_Y[it]/LG)*lmacro;
LIST_Z[it]=(LIST_Z[it]/LG)*lmacro;

LIST_XO[it] = LIST_X[it];
LIST_YO[it] = LIST_Y[it];
LIST_ZO[it] = LIST_Z[it];

LIST_XA[it] = LIST_X[it];
LIST_YA[it] = LIST_Y[it];
LIST_ZA[it] = LIST_Z[it];
		
LIST_R[it]=(LIST_R[it]/LG)*lmacro;

LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*1.;
LIST_M[it]  = dens*LIST_V[it];		
LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

R_SPH+=LIST_R[it];
RMAX=max(RMAX,LIST_R[it]);

}

H_TOT=(H_TOT/LG)*lmacro;
V_TOT=(V_TOT/LG)*lmacro;
Z_TOT=(Z_TOT/LG)*lmacro;

H_POS=(H_POS/LG)*lmacro;
V_POS=(V_POS/LG)*lmacro;
Z_POS=(Z_POS/LG)*lmacro;

R_SPH/=(NB_SPH);
epsi=(epsi/LG)*lmacro;
}

void scale_sph_ind(R Pi, R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS,  int & NB_SPH, R & R_SPH, R & RMAX, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_M, R * LIST_I, R & epsi, R lmacro, R dens, R rind, R * LIST_V)
{
R_SPH=0.;
RMAX=0.;  

R LG=H_TOT;	

R HG=H_POS;
R VG=V_POS;
R ZG=Z_POS;

R coef;
for (int it=0;it<NB_SPH-1;it++){

coef=1.;

if(fabs(LIST_X[it]-HG)<2.*epsi) {    coef/=2.; }
if(fabs(LIST_Y[it]-VG)<2.*epsi) {    coef/=2.; }
if(fabs(LIST_Z[it]-ZG)<2.*epsi) {   coef/=2.;}
if(fabs(LIST_X[it]-(LG+HG))<2.*epsi) {  coef/=2.;}	
if(fabs(LIST_Y[it]-(LG+VG))<2.*epsi) {  coef/=2.;} 
if(fabs(LIST_Z[it]-(LG+ZG))<2.*epsi) { coef/=2.;} 

LIST_X[it]=(LIST_X[it]/LG)*lmacro;
LIST_Y[it]=(LIST_Y[it]/LG)*lmacro;
LIST_Z[it]=(LIST_Z[it]/LG)*lmacro;

LIST_XO[it] = LIST_X[it];
LIST_YO[it] = LIST_Y[it];
LIST_ZO[it] = LIST_Z[it];

LIST_XA[it] = LIST_X[it];
LIST_YA[it] = LIST_Y[it];
LIST_ZA[it] = LIST_Z[it];
		
LIST_R[it]=(LIST_R[it]/LG)*lmacro;

LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;
LIST_M[it]  = dens*LIST_V[it];			
LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

R_SPH+=LIST_R[it];
RMAX=(RMAX>LIST_R[it])?RMAX:LIST_R[it];

}
R_SPH/=(NB_SPH-1);

H_TOT=(H_TOT/LG)*lmacro;
V_TOT=(V_TOT/LG)*lmacro;
Z_TOT=(Z_TOT/LG)*lmacro;

H_POS=(H_POS/LG)*lmacro;
V_POS=(V_POS/LG)*lmacro;
Z_POS=(Z_POS/LG)*lmacro;

LIST_X[NB_SPH-1]=R(H_TOT)/2.;
LIST_Y[NB_SPH-1]=R(V_TOT)/2.;
LIST_Z[NB_SPH-1]=R(Z_TOT)+rind+R_SPH;
LIST_R[NB_SPH-1]=rind;

LIST_V[NB_SPH-1]  = (4./3)*Pi*LIST_R[NB_SPH-1]*LIST_R[NB_SPH-1]*LIST_R[NB_SPH-1];
LIST_M[NB_SPH-1]  = 1.3*dens*LIST_V[NB_SPH-1];
LIST_I[NB_SPH-1]  = LIST_M[NB_SPH-1]*(2./5.)*LIST_R[NB_SPH-1]*LIST_R[NB_SPH-1];	

epsi=(epsi/LG)*lmacro;
}

void scale_sphp_ind(R Pi, R & H_TOT, R & V_TOT, R & Z_TOT,R & H_POS, R & V_POS, R & Z_POS,  int & NB_SPH, R & R_SPH, R & RMAX, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_M, R * LIST_I, R & epsi, R lmacro, R dens, R rind, R * LIST_V)
{
R_SPH=0.;
RMAX=0.;  

R LG=H_TOT;	

R HG=H_POS;
R VG=V_POS;
R ZG=Z_POS;

for (int it=0;it<NB_SPH-1;it++){

LIST_X[it]=(LIST_X[it]/LG)*lmacro;
LIST_Y[it]=(LIST_Y[it]/LG)*lmacro;
LIST_Z[it]=(LIST_Z[it]/LG)*lmacro;

LIST_XO[it] = LIST_X[it];
LIST_YO[it] = LIST_Y[it];
LIST_ZO[it] = LIST_Z[it];

LIST_XA[it] = LIST_X[it];
LIST_YA[it] = LIST_Y[it];
LIST_ZA[it] = LIST_Z[it];
		
LIST_R[it]=(LIST_R[it]/LG)*lmacro;

LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*1.;
LIST_M[it]  = dens*LIST_V[it];		
LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

R_SPH+=LIST_R[it];
RMAX=(RMAX>LIST_R[it])?RMAX:LIST_R[it];

}
R_SPH/=(NB_SPH-1);

H_TOT=(H_TOT/LG)*lmacro;
V_TOT=(V_TOT/LG)*lmacro;
Z_TOT=(Z_TOT/LG)*lmacro;

H_POS=(H_POS/LG)*lmacro;
V_POS=(V_POS/LG)*lmacro;
Z_POS=(Z_POS/LG)*lmacro;

LIST_X[NB_SPH-1]=R(H_TOT)/2.;
LIST_Y[NB_SPH-1]=R(V_TOT)/2.;
LIST_Z[NB_SPH-1]=R(Z_TOT)+rind;
LIST_R[NB_SPH-1]=rind;

LIST_V[NB_SPH-1]  = (4./3)*Pi*LIST_R[NB_SPH-1]*LIST_R[NB_SPH-1]*LIST_R[NB_SPH-1]*1.;
LIST_M[NB_SPH-1]  = 1.3*dens*LIST_V[NB_SPH-1];
LIST_I[NB_SPH-1]  = LIST_M[NB_SPH-1]*(2./5.)*LIST_R[NB_SPH-1]*LIST_R[NB_SPH-1];	

epsi=(epsi/LG)*lmacro;
}

void scale_sph2(R Pi, R & H_TOT, R & V_TOT, R & Z_TOT, R & H_POS, R & V_POS, R & Z_POS, int & NB_SPH, R & R_SPH, R & RMAX, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_M, R * LIST_I, R & epsi, R lmacro, R dens1, R dens2, bool * LIST_P, R * LIST_V)
{


R LG=H_TOT;	

R_SPH=0.;
RMAX=0.;  

R HG=H_POS;
R VG=V_POS;
R ZG=Z_POS;

R coef;

int it;
for (it=0;it<NB_SPH;it++){

coef=1.;

if(fabs(LIST_X[it]-HG)<2.*epsi) {    coef/=2.; }
if(fabs(LIST_Y[it]-VG)<2.*epsi) {    coef/=2.; }
if(fabs(LIST_Z[it]-ZG)<2.*epsi) {   coef/=2.;}
if(fabs(LIST_X[it]-(LG+HG))<2.*epsi) {  coef/=2.;}	
if(fabs(LIST_Y[it]-(LG+VG))<2.*epsi) {  coef/=2.;} 
if(fabs(LIST_Z[it]-(LG+ZG))<2.*epsi) { coef/=2.;} 

LIST_X[it]=(LIST_X[it]/LG)*lmacro;
LIST_Y[it]=(LIST_Y[it]/LG)*lmacro;
LIST_Z[it]=(LIST_Z[it]/LG)*lmacro;

LIST_XO[it] = LIST_X[it];
LIST_YO[it] = LIST_Y[it];
LIST_ZO[it] = LIST_Z[it];

LIST_XA[it] = LIST_X[it];
LIST_YA[it] = LIST_Y[it];
LIST_ZA[it] = LIST_Z[it];

LIST_R[it]=(LIST_R[it]/LG)*lmacro;
LIST_V[it]= (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;

if(LIST_P[it]==0){	
LIST_M[it]  = dens1*LIST_V[it];		
}
else{
LIST_M[it]  = dens2*LIST_V[it];			
}	
		
LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

R_SPH+=LIST_R[it];
RMAX=max(RMAX,LIST_R[it]);

}
R_SPH/=(NB_SPH);

H_TOT=(H_TOT/LG)*lmacro;
V_TOT=(V_TOT/LG)*lmacro;
Z_TOT=(Z_TOT/LG)*lmacro;

H_POS=(H_POS/LG)*lmacro;
V_POS=(V_POS/LG)*lmacro;
Z_POS=(Z_POS/LG)*lmacro;

epsi=(epsi/LG)*lmacro;
}


void init_sphp_cyl(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul, int * nocoul, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT, R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,bool * EDGEC, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT
for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}

int it,nbc;
int NX,NY,NZ;

PR_SPH=0.;

int nb1=0;
int nb2=0;
int nb3=0;
int nb4=0;
int nb5=0;
int nb6=0;
nbc=0;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

	// DISTRIBUTION INITIALE DE SPHERES
	for(it=0;it<NB_SPH;it++){
	LIST_IND[it]=1.;
      
      PR_SPH+=((4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]);
   
      LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;   
      EDGEC[it]=0;     
               
      NBCONTCO[it]= 0;

      if(fabs(LIST_X[it]-LIST_R[it]-H_POS)<epsi) { EDGE1[it]=1;  nb1++; }
      if(fabs(LIST_Y[it]-LIST_R[it]-V_POS)<epsi) { EDGE2[it]=1;  nb2++; }
      if(fabs(LIST_Z[it]-LIST_R[it]-Z_POS)<epsi) { EDGE3[it]=1;  nb3++; }
      if(fabs(LIST_X[it]+LIST_R[it]-(H_TOT+H_POS))<epsi) { EDGE4[it]=1;  nb4++; }
      if(fabs(LIST_Y[it]+LIST_R[it]-(V_TOT+V_POS))<epsi) { EDGE5[it]=1;  nb5++; }
      if(fabs(LIST_Z[it]+LIST_R[it]-(Z_TOT+Z_POS))<epsi) { EDGE6[it]=1;  nb6++; }   
  
     R dist=(LIST_X[it]-H_TOT/2.-H_POS)*(LIST_X[it]-H_TOT/2.-H_POS);
     dist=dist+(LIST_Y[it]-V_TOT/2.-V_POS)*(LIST_Y[it]-V_TOT/2.-V_POS);
     dist=sqrt(dist)+LIST_R[it];
     if(fabs(dist-H_TOT/2.-H_POS)<5.*epsi) { EDGEC[it]=1; nbc++;    }     

    		
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		

		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;

		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_XA[it] = LIST_X[it];
		LIST_YA[it] = LIST_Y[it];
		LIST_ZA[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;

		LIST_TXA[it] = 0.;
		LIST_TYA[it] = 0.;
		LIST_TZA[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*1.;
		LIST_M[it]  = dens*LIST_V[it];		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	

	}
	
	PR_SPH=PR_SPH/(H_TOT*V_TOT*Z_TOT);

N1=R(nb1);
N2=R(nb2);
N3=R(nb3);
N4=R(nb4);
N5=R(nb5);
N6=R(nb6);

cout<<"nbc:"<<nbc<<endl;
}

void init_sph_halo_cyl(bool * TYPCO, int NMAXCONT, R Pi, R dens, R epsi, int vecsizex, int vecsizey, int vecsizez,int ** coul,int * nocoul,int vecsizexh, int vecsizeyh, int vecsizezh,vector< vector<int>  > & coulh, R & PR_SPH, R H_TOT, R V_TOT, R Z_TOT,R H_POS, R V_POS, R Z_POS, int NB_SPH, R R_SPH, int * LIST_B, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, int * LIST_C, int * LIST_H, R * LIST_M, R * LIST_I, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE, bool * EDGE1, bool * EDGE2, bool * EDGE3, bool * EDGE4,  bool * EDGE5, bool * EDGE6,bool * EDGEC, int * NBCONTCO, R * LIST_V, R * LIST_IND, R & N1, R & N2, R & N3, R & N4, R & N5, R & N6)
{

    // CONTACTS COHESIFS INITIALEMENT

for(int it=0;it<NMAXCONT;it++){
TYPCO[it]=1;		
}

PR_SPH=0.;

R coef;
int NX,NY,NZ;
int it;
int n1,n2,n3,n4,n5,n6;
n1=0.;n2=0.;n3=0.;n4=0.;n5=0.;n6=0.;

N1=0.;
N2=0.;
N3=0.;
N4=0.;
N5=0.;
N6=0.;

for(it=0;it<vecsizex*vecsizey*vecsizez;it++){
nocoul[it]=0;
}

R volt=0.;

	// DISTRIBUTION INITIALE DE SPHERES
	for(it=0;it<NB_SPH;it++){
        
      LIST_B[it]=0.;   
      EDGE[it]=0;      
      EDGE1[it]=0;
      EDGE2[it]=0;
      EDGE3[it]=0;
      EDGE4[it]=0;
      EDGE5[it]=0;
      EDGE6[it]=0;     
      EDGEC[it]=0;     
         
      NBCONTCO[it]= 0;
      coef=1.;
      
      if(fabs(LIST_X[it]-H_POS)<epsi/5.) {n1++; EDGE1[it]=1;    }
      if(fabs(LIST_Y[it]-V_POS)<epsi/5.) {n2++; EDGE2[it]=1;    }
      if(fabs(LIST_Z[it]-Z_POS)<epsi/5.) {n3++; EDGE3[it]=1;   EDGE[it]=1; coef/=2.;}
      if(fabs(LIST_X[it]-(H_TOT+H_POS))<epsi/5.) {n4++; EDGE4[it]=1;   }	
      if(fabs(LIST_Y[it]-(V_TOT+V_POS))<epsi/5.) {n5++; EDGE5[it]=1;   } 
      if(fabs(LIST_Z[it]-(Z_TOT+Z_POS))<epsi/5.) {n6++; EDGE6[it]=1;    EDGE[it]=1;coef/=2.; }    

     R dist=(LIST_X[it]-H_TOT/2.-H_POS)*(LIST_X[it]-H_TOT/2.-H_POS);
     dist=dist+(LIST_Y[it]-V_TOT/2.-V_POS)*(LIST_Y[it]-V_TOT/2.-V_POS);
     dist=sqrt(dist);
     if(fabs(dist-H_TOT/2.-H_POS)<epsi/5.) { EDGEC[it]=1; EDGE[it]=1;coef/=2.; }    

      if(EDGE1[it]) N1+=coef;
      if(EDGE2[it]) N2+=coef;
      if(EDGE3[it]) N3+=coef;
      if(EDGE4[it]) N4+=coef;
      if(EDGE5[it]) N5+=coef;
      if(EDGE6[it]) N6+=coef;  	

	LIST_IND[it]=coef;
  
      PR_SPH+=((4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it])*coef;
    		
       // grille détection		
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;	
	

		coul[NX+NY*vecsizex+NZ*vecsizex*vecsizey][nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]]=it;
	nocoul[NX+NY*vecsizex+NZ*vecsizex*vecsizey]++;
		
		LIST_C[it] = NX+NY*vecsizex+NZ*vecsizex*vecsizey;   

       // grille halo	
		NX= ((int) ((LIST_X[it]-H_POS)*vecsizexh/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizeyh/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizezh/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizexh) NX=vecsizexh-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizeyh) NY=vecsizeyh-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizezh) NZ=vecsizezh-1;		

		coulh[NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh].push_back(it);
		LIST_H[it] = NX+NY*vecsizexh+NZ*vecsizexh*vecsizeyh;   

		LIST_XO[it] = LIST_X[it];
		LIST_YO[it] = LIST_Y[it];
		LIST_ZO[it] = LIST_Z[it];

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;
		LIST_TZ[it] = 0.;
		
		LIST_VX[it] = 0.;
		LIST_VY[it] = 0.; 
		LIST_VZ[it] = 0.;
					
		LIST_WX[it] = 0.;        
		LIST_WY[it] = 0.;    
		LIST_WZ[it] = 0.;		
		
		LIST_AX[it] = 0.;
		LIST_AY[it] = 0.;
		LIST_AZ[it] = 0.;	
			
		LIST_AWX[it] = 0.;  
		LIST_AWY[it] = 0.;  		 
		LIST_AWZ[it] = 0.;  
			    		
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;	
				
		MTX[it] = 0.;		
		MTY[it] = 0.;					
		MTZ[it] = 0.;		

		LIST_V[it]  = (4./3)*Pi*LIST_R[it]*LIST_R[it]*LIST_R[it]*coef;
		LIST_M[it]  = dens*LIST_V[it];		
		LIST_I[it]  = LIST_M[it]*(2./5.)*LIST_R[it]*LIST_R[it];	
               volt+=  LIST_V[it];
	}
	
	PR_SPH=PR_SPH/(3.14159*H_TOT*V_TOT*Z_TOT/4.);

cout<<"n1:"<<n1<<" -- n2:"<<n2<<" -- n3:"<<n3<<" -- n4:"<<n4<<" -- n5:"<<n5<<" -- n6:"<<n6<<endl;

}


