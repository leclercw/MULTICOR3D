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

#include "correct.h"
# include "omp.h"

void verlet_tracy(R & ftot, R & ftot2, R & Ec, R & Epp, R viti, int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul,int * LIST_C,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS,int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6)
{
ftot=0.;
ftot2=0.;
R Ec1=0.;
R Ep1=0.;
int it,jt,numxyz,num_threads,num_th;
int NX,NY,NZ;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}


if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}


	for(it=0;it<nD;it++){

	if((!EDGE2[it])&&(!EDGE5[it])) {	
	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
	}
	else if(EDGE5[it])
	{
	LIST_VY[it] = viti;
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it]  = LIST_Y[it]+dt*viti;	
	ftot+=FY[it];
	}
	else{
	LIST_VY[it] = -viti;
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it]  = LIST_Y[it]-dt*viti;	
	ftot2+=FY[it];
	}
	
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];    

	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	LIST_AZ[it] = FZ[it]/LIST_M[it];	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

	LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

	LIST_AWX[it] = MTX[it]/LIST_I[it];
	LIST_AWY[it] = MTY[it]/LIST_I[it];	 
	LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

	LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
	LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
	LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
  
	 Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	 Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
          
	 Ep1-=0.*(LIST_X[it]-LIST_XO[it]);
	 Ep1-=0.*(LIST_Y[it]-LIST_YO[it]);
 	 Ep1-=0.*(LIST_Z[it]-LIST_ZO[it]);         
          
	 FX[it] = 0.;
	 FY[it] = 0.;	
	 FZ[it] = 0.;
	 		 
         MTX[it] = 0.;
         MTY[it] = 0.;   
         MTZ[it] = 0.; 

	NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
	NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
	NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));			

	if(NX<0) NX=0;
	if(NX>=vecsizex) NX=vecsizex-1;		
	if(NY<0) NY=0;
	if(NY>=vecsizey) NY=vecsizey-1;		
	if(NZ<0) NZ=0;
	if(NZ>=vecsizez) NZ=vecsizez-1;		

        numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
	LIST_C[it] = numxyz;   
	coul[numxyz][nocoul[numxyz]]=it;
nocoul[numxyz]++;
	
	}

}else{

	int coul2[vecsize][num_threads-1][15]; 
        int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
		nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}


# pragma omp parallel for schedule(static) private(it,NX,NY,NZ,numxyz,num_th)  reduction(+:Ec1,ftot,ftot2) reduction(-:Ep1)
	for(it=0;it<nD;it++){

	num_th=omp_get_thread_num()-1;

	if((!EDGE2[it])&&(!EDGE5[it])) {	
	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
	}
	else if(EDGE5[it])
	{
	LIST_VY[it] = viti;
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it]  = LIST_Y[it]+dt*viti;	
	ftot+=FY[it];
	}
	else{
	LIST_VY[it] = -viti;
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it]  = LIST_Y[it]-dt*viti;	
	ftot2+=FY[it];
	}
	
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];    

	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	LIST_AZ[it] = FZ[it]/LIST_M[it];	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

	LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

	LIST_AWX[it] = MTX[it]/LIST_I[it];
	LIST_AWY[it] = MTY[it]/LIST_I[it];	 
	LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

	LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
	LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
	LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
  
	 Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	 Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
          
	 Ep1-=0.*(LIST_X[it]-LIST_XO[it]);
	 Ep1-=0.*(LIST_Y[it]-LIST_YO[it]);
 	 Ep1-=0.*(LIST_Z[it]-LIST_ZO[it]);         
          
	 FX[it] = 0.;
	 FY[it] = 0.;	
	 FZ[it] = 0.;
	 		 
         MTX[it] = 0.;
         MTY[it] = 0.;   
         MTZ[it] = 0.;  

	NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
	NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
	NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));			

	if(NX<0) NX=0;
	if(NX>=vecsizex) NX=vecsizex-1;		
	if(NY<0) NY=0;
	if(NY>=vecsizey) NY=vecsizey-1;		
	if(NZ<0) NZ=0;
	if(NZ>=vecsizez) NZ=vecsizez-1;		

        numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
	LIST_C[it] = numxyz;  

	if(num_th==-1){
	coul[numxyz][nocoul[numxyz]]=it;
	nocoul[numxyz]++;
	}else{
	coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
	nocoul2[numxyz][num_th]++;
	}


	}

	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}

}

 Ec=Ec1;
 Epp=Ep1;

}

void verlet_tracy_qs(R & ftot, R & ftot2, R & Ec, R & Epp, R viti, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6)
{
ftot=0.;
ftot2=0.;
R Ec1=0.;
int it;

# pragma omp parallel for schedule(static) private(it)  reduction(+:Ec1,ftot,ftot2) 
for(it=0;it<nD;it++){

	if((!EDGE2[it])&&(!EDGE5[it])) {	
	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
	}
	else if(EDGE5[it])
	{
	LIST_VY[it] = viti;
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it]  = LIST_Y[it]+dt*viti;	
	ftot+=FY[it];
	}
	else{
	LIST_VY[it] = -viti;
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it]  = LIST_Y[it]-dt*viti;	
	ftot2+=FY[it];
	}
	
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];    

	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	LIST_AZ[it] = FZ[it]/LIST_M[it];	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

	LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

	LIST_AWX[it] = MTX[it]/LIST_I[it];
	LIST_AWY[it] = MTY[it]/LIST_I[it];	 
	LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

	LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
	LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
	LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
  
	 Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	 Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
            
	 FX[it] = 0.;
	 FY[it] = 0.;	
	 FZ[it] = 0.;
	 		 
         MTX[it] = 0.;
         MTY[it] = 0.;   
         MTZ[it] = 0.;   


}

 Ec=Ec1;

}

void verlet_cis_cyl(int NBENREG, R & ftot,R & ftot2,R & ftot3, R epsi, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ,R ** LIST_PVX, R ** LIST_PVY, R ** LIST_PVZ,  R & Ec, R & Epp, int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul, int * LIST_C,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, bool * EDGEC1,bool * EDGEC2, R fimp, int nitf)
{

ftot=0.;
ftot2=0.;



int NX,NY,NZ;
int it,jt,num_th,numxyz,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}

	for(it=0;it<nD;it++){
		
		    if((!EDGEC1[it])&&(!EDGEC2[it])&&(!EDGE3[it])&&(!EDGE6[it])  ){
			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];   

                        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];			

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
		    }
			else if(EDGE6[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);	
/*
                        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];			

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; */
	
			}	    
		    else if(EDGEC2[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);	
/*
                        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];			

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; */
			}
			else if(EDGE3[it]){				
			ftot=ftot+FX[it];	
			ftot2=ftot2+FY[it];		
			ftot3=ftot3+FZ[it];					
			}	
			else if(EDGEC1[it]){
			ftot=ftot+FX[it];	
			ftot2=ftot2+FY[it];	
			ftot3=ftot3+FZ[it];						
			}
        
	
			Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
			Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
							
			FX[it] = 0.;
			FY[it] = 0.;	
			FZ[it] = 0.;
				 
			MTX[it] = 0.;
			MTY[it] = 0.;   
			MTZ[it] = 0.;   

			NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
			NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
			NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
			
			if(NX<0) NX=0;
			if(NX>=vecsizex) NX=vecsizex-1;		
			if(NY<0) NY=0;
			if(NY>=vecsizey) NY=vecsizey-1;		
			if(NZ<0) NZ=0;
			if(NZ>=vecsizez) NZ=vecsizez-1;		
			

		numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
		LIST_C[it] = numxyz;   
		coul[numxyz][nocoul[numxyz]]=it;
		nocoul[numxyz]++;


	}


}else{

	int coul2[vecsize][num_threads-1][15]; 
        int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
		nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}

# pragma omp parallel for schedule(static) private(it,NX,NY,NZ,num_th,numxyz) reduction(+:Ec,ftot,ftot2,ftot3) 
	for(it=0;it<nD;it++){
		
	num_th=omp_get_thread_num()-1;
		
		    if((!EDGEC1[it])&&(!EDGEC2[it])&&(!EDGE3[it])&&(!EDGE6[it])  ){
			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];   

                        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];			

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
		    }
			else if(EDGE6[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);	
/*
                        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];			

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; */
	
			}	    
		    else if(EDGEC2[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);	
/*
                        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];			

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; */
			}
			else if(EDGE3[it]){				
			ftot=ftot+FX[it];	
			ftot2=ftot2+FY[it];		
			ftot3=ftot3+FZ[it];					
			}	
			else if(EDGEC1[it]){
			ftot=ftot+FX[it];	
			ftot2=ftot2+FY[it];	
			ftot3=ftot3+FZ[it];						
			}
        
	
			Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
			Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
							
			FX[it] = 0.;
			FY[it] = 0.;	
			FZ[it] = 0.;
				 
			MTX[it] = 0.;
			MTY[it] = 0.;   
			MTZ[it] = 0.;   

			NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
			NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
			NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
			
			if(NX<0) NX=0;
			if(NX>=vecsizex) NX=vecsizex-1;		
			if(NY<0) NY=0;
			if(NY>=vecsizey) NY=vecsizey-1;		
			if(NZ<0) NZ=0;
			if(NZ>=vecsizez) NZ=vecsizez-1;		
			

		numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
		LIST_C[it] = numxyz;   
		if(num_th==-1){
		coul[numxyz][nocoul[numxyz]]=it;
		nocoul[numxyz]++;
		}else{
		coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
		nocoul2[numxyz][num_th]++;
		}


	}


	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}
}
   	    
   	  
}

void verlet_ind(R & Ec, R & Epp, int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul,int * LIST_C,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT, R H_POS,R V_POS,R Z_POS,int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, int NBENREG)
{

R ftot=0.;

R fite;
if(ite<50000){
fite=ite*0.01/50000.;	
}else{
fite=0.01;		
}

Ec=0.;
Epp=0.;

int NX,NY,NZ;
int it,jt,num_threads,numxyz,num_th;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}


	for(it=0;it<nD;it++){
	
	     if(it==nD-1){
			/*LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FZ[it]-fite)/LIST_M[it]);    
			LIST_AZ[it] = (FZ[it]-fite)/LIST_M[it];	
			LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	*/

			 
		LIST_VZ[it] =-viti;
		LIST_AZ[it] = FZ[it]/LIST_M[it];
		LIST_Z[it]  = LIST_Z[it]-dt*viti;	

		LIST_WX[it]  = 0.;
		LIST_WY[it]  = 0.;
		LIST_WZ[it]  = 0.;	

		LIST_AWX[it] = 0.;
		LIST_AWY[it] = 0.;
		LIST_AWZ[it] = 0.;   	 

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;    
		LIST_TZ[it] = 0.; 

		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   

	      }
	      else{
		    if(!EDGE3[it]) {	    

			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

			LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];

			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
		    }else{

		     ftot+=FZ[it];
		    }

		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
	    
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;  

		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		

		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		

		numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
		LIST_C[it] = numxyz;   
		coul[numxyz][nocoul[numxyz]]=it;
		nocoul[numxyz]++;


	      }

	}

}else{

	int coul2[vecsize][num_threads-1][15]; 
        int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
		nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}

	# pragma omp parallel for schedule(static) private(it,numxyz,num_th,NX,NY,NZ) reduction(+:Ec,ftot) 
	for(it=0;it<nD;it++){
	
	num_th=omp_get_thread_num()-1;

	     if(it==nD-1){
			/*LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FZ[it]-fite)/LIST_M[it]);    
			LIST_AZ[it] = (FZ[it]-fite)/LIST_M[it];	
			LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	*/

			 
		LIST_VZ[it] =-viti;
		LIST_AZ[it] = FZ[it]/LIST_M[it];
		LIST_Z[it]  = LIST_Z[it]-dt*viti;	

		LIST_WX[it]  = 0.;
		LIST_WY[it]  = 0.;
		LIST_WZ[it]  = 0.;	

		LIST_AWX[it] = 0.;
		LIST_AWY[it] = 0.;
		LIST_AWZ[it] = 0.;   	 

		LIST_TX[it] = 0.;
		LIST_TY[it] = 0.;    
		LIST_TZ[it] = 0.; 

		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   

	      }
	      else{
		    if(!EDGE3[it]) {	    

			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

			LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];

			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
		    }else{

		     ftot+=FZ[it];
		    }

		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
	    
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;  

		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		

		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		

		numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
		LIST_C[it] = numxyz;   

		if(num_th==-1){
		coul[numxyz][nocoul[numxyz]]=it;
		nocoul[numxyz]++;
		}else{
		coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
		nocoul2[numxyz][num_th]++;
		}

	      }

	}

	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}


}


}



void verlet(R & Ec, R & Epp,int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul,int * LIST_C,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS,int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ)
{
R Ec1=0.;
R Ep1=0.;
int it,jt,numxyz,num_threads,num_th;
int NX,NY,NZ;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}


if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}


	for(it=0;it<nD;it++){

	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];

	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]-1e7);    
	LIST_AZ[it] = (FZ[it]/LIST_M[it])-1e7;	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

	LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

	LIST_AWX[it] = MTX[it]/LIST_I[it];
	LIST_AWY[it] = MTY[it]/LIST_I[it];	 
	LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

	LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
	LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
	LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

	Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

	Ep1-=0.*(LIST_X[it]-LIST_XO[it]);
	Ep1-=0.*(LIST_Y[it]-LIST_YO[it]);
	Ep1-=0.*(LIST_Z[it]-LIST_ZO[it]);         

	FX[it] = 0.;
	FY[it] = 0.;	
	FZ[it] = 0.;
		 
	MTX[it] = 0.;
	MTY[it] = 0.;   
	MTZ[it] = 0.;   

	NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
	NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
	NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));			

	if(NX<0) NX=0;
	if(NX>=vecsizex) NX=vecsizex-1;		
	if(NY<0) NY=0;
	if(NY>=vecsizey) NY=vecsizey-1;		
	if(NZ<0) NZ=0;
	if(NZ>=vecsizez) NZ=vecsizez-1;		

        numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
	LIST_C[it] = numxyz;   
	coul[numxyz][nocoul[numxyz]]=it;
nocoul[numxyz]++;
	
	}

}else{

	int coul2[vecsize][num_threads-1][15]; 
        int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}

# pragma omp parallel for schedule(static) private(it,NX,NY,NZ,numxyz,num_th)  reduction(+:Ec1) reduction(-:Ep1)
	for(it=0;it<nD;it++){

	num_th=omp_get_thread_num()-1;

	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];

	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]-5e7);    
	LIST_AZ[it] = (FZ[it]/LIST_M[it])-5e7;	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

	LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

	LIST_AWX[it] = MTX[it]/LIST_I[it];
	LIST_AWY[it] = MTY[it]/LIST_I[it];	 
	LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

	LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
	LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
	LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

	Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

	Ep1-=0.*(LIST_X[it]-LIST_XO[it]);
	Ep1-=0.*(LIST_Y[it]-LIST_YO[it]);
	Ep1-=0.*(LIST_Z[it]-LIST_ZO[it]);         

	FX[it] = 0.;
	FY[it] = 0.;	
	FZ[it] = 0.;
		 
	MTX[it] = 0.;
	MTY[it] = 0.;   
	MTZ[it] = 0.;   

	NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
	NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
	NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));			

	if(NX<0) NX=0;
	if(NX>=vecsizex) NX=vecsizex-1;		
	if(NY<0) NY=0;
	if(NY>=vecsizey) NY=vecsizey-1;		
	if(NZ<0) NZ=0;
	if(NZ>=vecsizez) NZ=vecsizez-1;		

        numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
	LIST_C[it] = numxyz;   

if(num_th==-1){
coul[numxyz][nocoul[numxyz]]=it;
nocoul[numxyz]++;
}else{
coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
nocoul2[numxyz][num_th]++;
}



	}

	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}
}


}

void verlet_qs(R & Ec, R & Epp,int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ)
{

R Ec1=0.;
R Ep1=0.;
int it;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec1) reduction(-:Ep1)
	for(it=0;it<nD;it++){

	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];

	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	LIST_AZ[it] = FZ[it]/LIST_M[it];	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

	LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

	LIST_AWX[it] = MTX[it]/LIST_I[it];
	LIST_AWY[it] = MTY[it]/LIST_I[it];	 
	LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

	LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
	LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
	LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

	Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

	Ep1-=0.*(LIST_X[it]-LIST_XO[it]);
	Ep1-=0.*(LIST_Y[it]-LIST_YO[it]);
	Ep1-=0.*(LIST_Z[it]-LIST_ZO[it]);         

	FX[it] = 0.;
	FY[it] = 0.;	
	FZ[it] = 0.;
		 
	MTX[it] = 0.;
	MTY[it] = 0.;   
	MTZ[it] = 0.;   

	}
 

}


void verlet_symt_qs(R & Ec, R & Epp, bool &bout,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, int NBENREG)
{
R ftot=0.;
R dplty=0.; 
R dpltz=0.; 
int ndp=0;
Ec=0.;
Epp=0.;

int it,num_threads;
R dt2=dt/2.;
R ddt2=dt*dt/2.;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec,ndp,dpltz,dplty,ftot)  
	for(it=0;it<nD;it++){

	if((!EDGE1[it])&&(!EDGE4[it])) {	
	LIST_VX[it] = LIST_VX[it]+dt2*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+ddt2*LIST_AX[it];
	}
	else if(EDGE4[it])
	{
	LIST_VX[it] = viti;
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it]  = LIST_X[it]+dt*viti;	
	}
	else{
	ftot+=FX[it];
	}
	
    
	if(!EDGE2[it]) {	
	LIST_VY[it] = LIST_VY[it]+dt2*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+ddt2*LIST_AY[it];
	}

	if(!EDGE3[it]) {	    
	 LIST_VZ[it] = LIST_VZ[it]+dt2*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	 LIST_AZ[it] = FZ[it]/LIST_M[it];	
	 LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+ddt2*LIST_AZ[it];	
	}

	if(EDGE6[it]&&EDGE5[it]){
	dpltz+= (LIST_Z[it]-LIST_ZO[it]);
	dplty+= (LIST_Y[it]-LIST_YO[it]);
	ndp+=1;
	}
    
       
 	 LIST_WX[it]  = LIST_WX[it]+dt2*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	 LIST_WY[it]  = LIST_WY[it]+dt2*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	 LIST_WZ[it]  = LIST_WZ[it]+dt2*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	
	 
     LIST_AWX[it] = MTX[it]/LIST_I[it];
     LIST_AWY[it] = MTY[it]/LIST_I[it];	 
     LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

     LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+ddt2*LIST_AWX[it];
     LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+ddt2*LIST_AWY[it];     
     LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+ddt2*LIST_AWZ[it]; 
 
	 Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	 Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
          
	 FX[it] = 0.;
	 FY[it] = 0.;	
	 FZ[it] = 0.;
	 		 
     MTX[it] = 0.;
     MTY[it] = 0.;   
     MTZ[it] = 0.;   



}


 dplty/=ndp;
 dpltz/=ndp;

  // Eval nu
  
  R nuy=dplty/((ite-1)*dt*viti);
  R nuz=dpltz/((ite-1)*dt*viti);


  // Eval E

		 R EE,sigxx,epsx;	     
		 sigxx=ftot/(Z_TOT*V_TOT);
		 epsx=(ite-1)*dt*viti/H_TOT;

                 EE=sigxx/epsx; 
                 if(ite>int(nD/10.)){
					 Emoy=(Emoy*(ite-int(nD/10.)-1)+EE)/(ite-int(nD/10.));
				 } 
         	  		
	     if(ite%NBENREG==0){ 
			 
  				if((ite>int(nD/10.))&&(abs((Emoy1-Emoy)/Emoy)<1e-3)) bout=1;					 
				Emoy1=Emoy;

                cout<<"Defx "<<epsx<<", Sigxx :"<<sigxx<<", Module de Young inst :"<<EE<<", Module de Young moy :"<<Emoy<<endl;	              
                cout<<"nu : "<<nuy<<", "<<nuz<<endl;
                nuu=-(nuy+nuz)/2.;

			  	
		}      

}

void verlet_qs_dila(R & Ec, R & Epp, int ite,int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6,int NBENREG,R & dpltx,R & dplty,R & dpltz)
{

R ftotx=0.;
R ftoty=0.;
R ftotz=0.;

dpltx=0.; 
dplty=0.; 
dpltz=0.; 

int ndpx=0;
int ndpy=0;
int ndpz=0;

R Ec1=0.;
R Ep1=0.;

int it;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec1,ndpx,ndpy,ndpz,dpltx,dplty,dpltz,ftotx,ftoty,ftotz) reduction(-:Ep1)
	for(it=0;it<nD;it++){

	if(!EDGE1[it]) {	
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
	}
	else{
	ftotx+=FX[it];
	}	
    
	if(!EDGE2[it]) {	
	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
	}
	else{
	ftoty+=FY[it];
	}

	if(!EDGE3[it]) {	    
	 LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	 LIST_AZ[it] = FZ[it]/LIST_M[it];	
	 LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
	}
	else{
	ftotz+=FZ[it];
	}

	if(EDGE4[it]){
	dpltx+= (LIST_X[it]-LIST_XO[it]);
	ndpx+=1;
	}    
  	else if(EDGE5[it]){
	dplty+= (LIST_Y[it]-LIST_YO[it]);
	ndpy+=1;
	}         
  	else if(EDGE6[it]){
	dpltz+= (LIST_Z[it]-LIST_ZO[it]);
	ndpz+=1;
	}     

	LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

	LIST_AWX[it] = MTX[it]/LIST_I[it];
	LIST_AWY[it] = MTY[it]/LIST_I[it];	 
	LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

	LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
	LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
	LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

	Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

	Ep1-=0.*(LIST_X[it]-LIST_XO[it]);
	Ep1-=0.*(LIST_Y[it]-LIST_YO[it]);
	Ep1-=0.*(LIST_Z[it]-LIST_ZO[it]);         

	FX[it] = 0.;
	FY[it] = 0.;	
	FZ[it] = 0.;
		 
	MTX[it] = 0.;
	MTY[it] = 0.;   
	MTZ[it] = 0.;   



}

 Ec=Ec1;
 Epp=Ep1;

 dpltx/=ndpx;
 dplty/=ndpy;
 dpltz/=ndpz;

if(ite%NBENREG==0) {
cout<<"dpct:"<<dpltx<<", "<<dplty<<", "<<dpltz<<endl;
}
 
}

void verlet_qs_hygr(R & Ec, R & Epp, int ite,int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6,int NBENREG,R & dpltx,R & dplty,R & dpltz)
{

dpltx=0.; 
dplty=0.; 
dpltz=0.; 

int ndpx=0;
int ndpy=0;
int ndpz=0;

R Ec1=0.;
R Ep1=0.;

int it;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec1,ndpx,ndpy,ndpz,dpltx,dplty,dpltz) reduction(-:Ep1)
	for(it=0;it<nD;it++){
	
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
	 
	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    
	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	LIST_AZ[it] = FZ[it]/LIST_M[it];	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

	if(EDGE4[it]){
	dpltx+= (LIST_X[it]-LIST_XO[it]);
	ndpx+=1;
	}    
  	else if(EDGE5[it]){
	dplty+= (LIST_Y[it]-LIST_YO[it]);
	ndpy+=1;
	}         
  	else if(EDGE6[it]){
	dpltz+= (LIST_Z[it]-LIST_ZO[it]);
	ndpz+=1;
	}     

	LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

	LIST_AWX[it] = MTX[it]/LIST_I[it];
	LIST_AWY[it] = MTY[it]/LIST_I[it];	 
	LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

	LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
	LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
	LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

	Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

	Ep1-=0.*(LIST_X[it]-LIST_XO[it]);
	Ep1-=0.*(LIST_Y[it]-LIST_YO[it]);
	Ep1-=0.*(LIST_Z[it]-LIST_ZO[it]);         

	FX[it] = 0.;
	FY[it] = 0.;	
	FZ[it] = 0.;
		 
	MTX[it] = 0.;
	MTY[it] = 0.;   
	MTZ[it] = 0.;   



}

 Ec=Ec1;
 Epp=Ep1;

 dpltx/=ndpx;
 dplty/=ndpy;
 dpltz/=ndpz;

if(ite%NBENREG==0) {
cout<<"dpct:"<<dpltx<<", "<<dplty<<", "<<dpltz<<endl;
}
 
}


void verlet_symt(R & Ec, R & Epp, int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul, int * LIST_C,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO,R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, int NBENREG)
{

R ftot=0.;
R dplty=0.; 
R dpltz=0.; 
int ndp=0;
Ec=0.;
Epp=0.;

int it,jt;
int NX,NY,NZ,numxyz,num_threads,num_th;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}

	for(it=0;it<nD;it++){
	
		if((!EDGE1[it])&&(!EDGE4[it])) {	
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
	    }
	    else if(EDGE4[it])
	    {
		LIST_VX[it] = -viti;
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it]  = LIST_X[it]-dt*viti;	
		}
		else{
		ftot+=FX[it];
	    }
	
	    
	    if(!EDGE2[it]) {	
	    LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	    LIST_AY[it] = FY[it]/LIST_M[it];
	    LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
	    }
	    
	    if(!EDGE3[it]) {	    
		 LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
		 LIST_AZ[it] = FZ[it]/LIST_M[it];	
		 LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
	    }
	    
	    if(EDGE6[it]&&EDGE5[it]){
		dpltz	= dpltz+(LIST_Z[it]-LIST_ZO[it]);
		dplty	= dplty+(LIST_Y[it]-LIST_YO[it]);
		ndp++;
		}
	    
	       
		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

		Epp-=0.*(LIST_X[it]-LIST_XO[it]);
		Epp-=0.*(LIST_Y[it]-LIST_YO[it]);
		Epp-=0.*(LIST_Z[it]-LIST_ZO[it]);    	       

		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   

		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
			if(NX<0) NX=0;
			if(NX>=vecsizex) NX=vecsizex-1;		
			if(NY<0) NY=0;
			if(NY>=vecsizey) NY=vecsizey-1;		
			if(NZ<0) NZ=0;
			if(NZ>=vecsizez) NZ=vecsizez-1;		

		numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
		LIST_C[it] = numxyz;   
		coul[numxyz][nocoul[numxyz]]=it;
		nocoul[numxyz]++;

	}

}else{

	int coul2[vecsize][num_threads-1][15]; 
        int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}


# pragma omp parallel for schedule(static) private(it,NX,NY,NZ,numxyz,num_th) reduction(+:Ec,ndp,dpltz,dplty,ftot) reduction(-:Epp)
for(it=0;it<nD;it++){
	
	num_th=omp_get_thread_num()-1;


	if((!EDGE1[it])&&(!EDGE4[it])) {	
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
    }
    else if(EDGE4[it])
    {
	LIST_VX[it] = -viti;
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it]  = LIST_X[it]-dt*viti;	
	}
	else{
	ftot+=FX[it];
    }
	
    
    if(!EDGE2[it]) {	
    LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
    LIST_AY[it] = FY[it]/LIST_M[it];
    LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
    }
    
    if(!EDGE3[it]) {	    
	 LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	 LIST_AZ[it] = FZ[it]/LIST_M[it];	
	 LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
    }
    
    if(EDGE6[it]&&EDGE5[it]){
	dpltz	= dpltz+(LIST_Z[it]-LIST_ZO[it]);
	dplty	= dplty+(LIST_Y[it]-LIST_YO[it]);
	ndp++;
	}
    
       
 	 LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	 LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	 LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	
	 
     LIST_AWX[it] = MTX[it]/LIST_I[it];
     LIST_AWY[it] = MTY[it]/LIST_I[it];	 
     LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

     LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
     LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
     LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
 
 	 Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	 Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
          
	 Epp-=0.*(LIST_X[it]-LIST_XO[it]);
	 Epp-=0.*(LIST_Y[it]-LIST_YO[it]);
 	 Epp-=0.*(LIST_Z[it]-LIST_ZO[it]);    	       
          
	 FX[it] = 0.;
	 FY[it] = 0.;	
	 FZ[it] = 0.;
	 		 
     MTX[it] = 0.;
     MTY[it] = 0.;   
     MTZ[it] = 0.;   

	NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
	NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
	NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		

        numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
	LIST_C[it] = numxyz;   

	if(num_th==-1){
	coul[numxyz][nocoul[numxyz]]=it;
	nocoul[numxyz]++;
	}else{
	coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
	nocoul2[numxyz][num_th]++;
	} 

}

	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}


}





 dplty/=ndp;
 dpltz/=ndp;

  // Eval nu
  
  R nuy=dplty/((ite-1)*dt*viti);
  R nuz=dpltz/((ite-1)*dt*viti);


  // Eval E

		 R EE,sigxx,epsx;	     
		 sigxx=ftot/(Z_TOT*V_TOT);
		 epsx=(ite-1)*dt*viti/H_TOT;

                 EE=fabs(sigxx/epsx); 
                 if(ite>int(nD/10.)){
					 Emoy=(Emoy*(ite-int(nD/10.)-1)+EE)/(ite-int(nD/10.));
				 } 
         	  		
	     if(ite%NBENREG==0){ 
			 
  			//	if((ite>int(nD/10.))&&(abs((Emoy1-Emoy)/Emoy)<1e-3)) bout=1;					 
				Emoy1=Emoy;
                cout<<"ite : "<<ite<<", dt : "<<dt<<", viti : "<<viti<<", htot : "<<H_TOT<<endl;
                cout<<"Defx "<<epsx<<", Sigxx :"<<sigxx<<", Module de Young inst :"<<EE<<", Module de Young moy :"<<Emoy<<endl;	
              
                cout<<"nu : "<<nuy<<", "<<nuz<<endl;
                nuu=-(nuy+nuz)/2.;
			  	
		}      

}


void verlet_trac(R & Ec, R & Epp, int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul,int * LIST_C,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, int NBENREG, R fimp, int nitf)
{

R fimp1;
if(ite<nitf){
fimp1=(fimp*ite)/nitf;
}
else
{
fimp1=fimp;	
}

R ftot=0.;
R dtot=0.;
int ntot=0;
Ec=0.;
Epp=0.;

int it,jt,numxyz,num_threads,num_th;
int NX,NY,NZ;

if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}

for(it=0;it<nD;it++){
	
		if((!EDGE1[it])&&(!EDGE4[it])) {	
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
		}
		else if(EDGE4[it])
		{
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+(FX[it]+fimp1)/LIST_M[it]);		
		LIST_AX[it] = (FX[it]+fimp1)/LIST_M[it];
		LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
		}
		else{
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+(FX[it]-fimp1)/LIST_M[it]);		
		LIST_AX[it] = (FX[it]-fimp1)/LIST_M[it];
		LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
		}
	
    
		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
		LIST_AY[it] = FY[it]/LIST_M[it];
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];

		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
		LIST_AZ[it] = FZ[it]/LIST_M[it];	
		LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
		
		if(EDGE1[it]) {	
		Epp-=(fimp1)*(LIST_X[it]-LIST_XO[it]);
		ftot+=FX[it];
		dtot+=(LIST_X[it]-LIST_XO[it]);
		ntot++;
	    }
	    
		if(EDGE4[it]) {	
		Epp-=(-fimp1)*(LIST_X[it]-LIST_XO[it]);
	    }

	  

		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   

		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
				
		numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
		LIST_C[it] = numxyz;   
		coul[numxyz][nocoul[numxyz]]=it;
		nocoul[numxyz]++;

}



}else if(num_threads>1){

	int coul2[vecsize][num_threads-1][15]; 
        int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
		nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}

# pragma omp parallel for schedule(static) private(it,NX,NY,NZ,numxyz,num_th) reduction(+:Ec,ftot,dtot,ntot) reduction(-:Epp)
for(it=0;it<nD;it++){
	
	num_th=omp_get_thread_num()-1;

		if((!EDGE1[it])&&(!EDGE4[it])) {	
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
		}
		else if(EDGE4[it])
		{
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+(FX[it]+fimp1)/LIST_M[it]);		
		LIST_AX[it] = (FX[it]+fimp1)/LIST_M[it];
		LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
		}
		else{
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+(FX[it]-fimp1)/LIST_M[it]);		
		LIST_AX[it] = (FX[it]-fimp1)/LIST_M[it];
		LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
		}
	
    
		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
		LIST_AY[it] = FY[it]/LIST_M[it];
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];

		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
		LIST_AZ[it] = FZ[it]/LIST_M[it];	
		LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
		
		if(EDGE1[it]) {	
		Epp-=(fimp1)*(LIST_X[it]-LIST_XO[it]);
		ftot+=FX[it];
		dtot+=(LIST_X[it]-LIST_XO[it]);
		ntot++;
	    }
	    
		if(EDGE4[it]) {	
		Epp-=(-fimp1)*(LIST_X[it]-LIST_XO[it]);
	    }

	  

		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   

		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
				
        numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
	LIST_C[it] = numxyz;   

	if(num_th==-1){
	coul[numxyz][nocoul[numxyz]]=it;
	nocoul[numxyz]++;
	}else{
	coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
	nocoul2[numxyz][num_th]++;
	}

}

	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}

}




dtot=dtot/ntot;

}


void verlet_comp(int NBENREG, R & ftot,R & ftot2, R epsi, R & Ec, R & Epp, int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul,int * LIST_C,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, R fimp, int nitf){
					
Epp=0.;

ftot=0.;
ftot2=0.;

int it,jt,NX,NY,NZ,numxyz,num_threads,num_th;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}

	for(it=0;it<nD;it++){
				
		    if((!EDGE1[it])&&(!EDGE4[it])){
			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];   

		    }
		    else if(EDGE4[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);				
			ftot=ftot+FX[it];	
			}
			else{
			ftot2=ftot2+FX[it];	
			}
            

                        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	  
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
			
			Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
			Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
	
			FX[it] = 0.;
			FY[it] = 0.;	
			FZ[it] = 0.;
				 
			MTX[it] = 0.;
			MTY[it] = 0.;   
			MTZ[it] = 0.;   

			NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
			NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
			NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));			
			
			if(NX<0) NX=0;
			if(NX>=vecsizex) NX=vecsizex-1;		
			if(NY<0) NY=0;
			if(NY>=vecsizey) NY=vecsizey-1;		
			if(NZ<0) NZ=0;
			if(NZ>=vecsizez) NZ=vecsizez-1;		
			
			numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
			LIST_C[it] = numxyz;   
			coul[numxyz][nocoul[numxyz]]=it;
			nocoul[numxyz]++; 

	}


}else{

	int coul2[vecsize][num_threads-1][15]; 
        int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
		nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}

# pragma omp parallel for schedule(static) private(it,NX,NY,NZ,numxyz,num_th) reduction(+:Ec,ftot,ftot2)
	for(it=0;it<nD;it++){
	
	num_th=omp_get_thread_num()-1;
			
		    if((!EDGE1[it])&&(!EDGE4[it])){
			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];   

		    }
		    else if(EDGE4[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);				
			ftot=ftot+FX[it];	
			}
			else{
			ftot2=ftot2+FX[it];	
			}
            

                        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	  
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
			
			Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
			Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
	
			FX[it] = 0.;
			FY[it] = 0.;	
			FZ[it] = 0.;
				 
			MTX[it] = 0.;
			MTY[it] = 0.;   
			MTZ[it] = 0.;   

			NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
			NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
			NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));			
			
			if(NX<0) NX=0;
			if(NX>=vecsizex) NX=vecsizex-1;		
			if(NY<0) NY=0;
			if(NY>=vecsizey) NY=vecsizey-1;		
			if(NZ<0) NZ=0;
			if(NZ>=vecsizez) NZ=vecsizez-1;		
			
			numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
			LIST_C[it] = numxyz;   

			if(num_th==-1){
			coul[numxyz][nocoul[numxyz]]=it;
			nocoul[numxyz]++;
			}else{
			coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
			nocoul2[numxyz][num_th]++;
			}

	}

	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}

}

   	    
   	  
}

void verlet_torsion(R & Ec, R & Epp, int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul, int * LIST_C,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, int NBENREG, R fimp,R mimp, int nitf)
{

R fimp1;
if(ite<nitf){
fimp1=(fimp*R(ite))/nitf;
}
else
{
fimp1=fimp;	
}

R mimp1;
if(ite<nitf){
mimp1=(mimp*R(ite))/nitf;
}
else
{
mimp1=mimp;	
}

R ftot1=0.;
R dtot1=0.;
R z1=0.;
R zO1=0.;
int ntot1=0;
R dtot2=0.;
R z2=0.;
R zO2=0.;
R ftot2=0.;
int ntot2=0;
R mtot1=0.;
R mtot2=0.;

Ec=0.;
Epp=0.;

int it,jt,numxyz,num_threads,num_th;
int NX,NY,NZ;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}

for(it=0;it<nD;it++){
	
		if((!EDGE3[it])&&(!EDGE6[it])) {	
			
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it]  = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
		LIST_AY[it] = FY[it]/LIST_M[it];
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
			
		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
		LIST_AZ[it] = (FZ[it])/LIST_M[it];	
		LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];
			
		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
					
		}
		else if(EDGE6[it])
		{
			
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it]  = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+(FY[it])/LIST_M[it]);
		LIST_AY[it] = (FY[it])/LIST_M[it];
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
			
		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FZ[it]-fimp1)/LIST_M[it]);    
		LIST_AZ[it] = (FZ[it]-fimp1)/LIST_M[it];	
		LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];			
			
		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+(MTZ[it]+mimp1)/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = (MTZ[it]+mimp1)/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
		
		}
		else{
			
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it]  = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+(FY[it])/LIST_M[it]);
		LIST_AY[it] = (FY[it])/LIST_M[it];
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
			
		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FZ[it]+fimp1)/LIST_M[it]);    
		LIST_AZ[it] = (FZ[it]+fimp1)/LIST_M[it];	
		LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];
			
		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+(MTZ[it]-mimp1)/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = (MTZ[it]-mimp1)/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
		}
	


		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
		
		if(EDGE3[it]) {	
	//	Epp-=(fimp1)*(LIST_Z[it]-LIST_ZO[it]);
		dtot1+=fabs(LIST_Z[it]-LIST_ZO[it]);
		ftot1+=FZ[it];
		z1+=LIST_Z[it]-LIST_R[it];
		zO1+=LIST_ZO[it]-LIST_R[it];	
		mtot1+=MTZ[it];
		ntot1++;
	    }
	    
		if(EDGE6[it]) {	
	//	Epp-=(-fimp1)*(LIST_Z[it]-LIST_ZO[it]);
		dtot2+=fabs(LIST_Z[it]-LIST_ZO[it]);
		ftot2+=FZ[it];
		z2+=LIST_Z[it]+LIST_R[it];		
		zO2+=LIST_ZO[it]+LIST_R[it];		
		mtot2+=MTZ[it];	
		ntot2++;	
	    }
   
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   

		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		
		numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
		LIST_C[it] = numxyz;   
		coul[numxyz][nocoul[numxyz]]=it;
		nocoul[numxyz]++;
	 

}

}else{

	int coul2[vecsize][num_threads-1][15]; 
	int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}


# pragma omp parallel for schedule(static) private(it,NX,NY,NZ,numxyz,num_th) reduction(+:Ec,ftot1,dtot1,z1,zO1,ntot1,dtot2,z2,zO2,ftot2,ntot2,mtot1,mtot2) reduction(-:Epp)
for(it=0;it<nD;it++){
	
	num_th=omp_get_thread_num()-1;


		if((!EDGE3[it])&&(!EDGE6[it])) {	
			
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it]  = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
		LIST_AY[it] = FY[it]/LIST_M[it];
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
			
		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
		LIST_AZ[it] = (FZ[it])/LIST_M[it];	
		LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];
			
		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
					
		}
		else if(EDGE6[it])
		{
			
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it]  = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+(FY[it])/LIST_M[it]);
		LIST_AY[it] = (FY[it])/LIST_M[it];
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
			
		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FZ[it]-fimp1)/LIST_M[it]);    
		LIST_AZ[it] = (FZ[it]-fimp1)/LIST_M[it];	
		LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];			
			
		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+(MTZ[it]+mimp1)/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = (MTZ[it]+mimp1)/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
		
		}
		else{
			
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
		LIST_AX[it] = FX[it]/LIST_M[it];
		LIST_X[it]  = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+(FY[it])/LIST_M[it]);
		LIST_AY[it] = (FY[it])/LIST_M[it];
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
			
		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FZ[it]+fimp1)/LIST_M[it]);    
		LIST_AZ[it] = (FZ[it]+fimp1)/LIST_M[it];	
		LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];
			
		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+(MTZ[it]-mimp1)/LIST_I[it]);	

		LIST_AWX[it] = MTX[it]/LIST_I[it];
		LIST_AWY[it] = MTY[it]/LIST_I[it];	 
		LIST_AWZ[it] = (MTZ[it]-mimp1)/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
		}
	


		Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
		
		if(EDGE3[it]) {	
	//	Epp-=(fimp1)*(LIST_Z[it]-LIST_ZO[it]);
		dtot1+=fabs(LIST_Z[it]-LIST_ZO[it]);
		ftot1+=FZ[it];
		z1+=LIST_Z[it]-LIST_R[it];
		zO1+=LIST_ZO[it]-LIST_R[it];	
		mtot1+=MTZ[it];
		ntot1++;
	    }
	    
		if(EDGE6[it]) {	
	//	Epp-=(-fimp1)*(LIST_Z[it]-LIST_ZO[it]);
		dtot2+=fabs(LIST_Z[it]-LIST_ZO[it]);
		ftot2+=FZ[it];
		z2+=LIST_Z[it]+LIST_R[it];		
		zO2+=LIST_ZO[it]+LIST_R[it];		
		mtot2+=MTZ[it];	
		ntot2++;	
	    }
   
		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   

		NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
		NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
		NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
		
		if(NX<0) NX=0;
		if(NX>=vecsizex) NX=vecsizex-1;		
		if(NY<0) NY=0;
		if(NY>=vecsizey) NY=vecsizey-1;		
		if(NZ<0) NZ=0;
		if(NZ>=vecsizez) NZ=vecsizez-1;		
		
        numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
	LIST_C[it] = numxyz;   

	if(num_th==-1){
	coul[numxyz][nocoul[numxyz]]=it;
	nocoul[numxyz]++;
	}else{
	coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
	nocoul2[numxyz][num_th]++;
	} 

}

	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}

}



   	    
	    dtot1=dtot1/ntot1;
	    dtot2=dtot2/ntot2;	    
	    z1=z1/ntot1;
	    z2=z2/ntot2;	   
	    zO1=zO1/ntot1;
	    zO2=zO2/ntot2;	   
	    	     	    
			
		if(ite%NBENREG==0){

         cout<<"Moment théorique :" <<ntot2*1e-5*ite/100.<<", eval. :"<<mtot1<<"  "<<mtot2<<endl;
         cout<<"Forces :"<<ftot1<<"  "<<ftot2<<endl;
			  	
		}      
		
	

}

void verlet_brazil(int NBENREG, R & ftot,R & ftot2, R epsi, R & Ec, R & Epp, int vecsize, int vecsizex, int vecsizey,int vecsizez, int ** coul, int * nocoul, int * LIST_C,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, R fimp, int nitf)
{

					
Epp=0.;
		    
ftot=0.;
ftot2=0.;

int NX,NY,NZ,it,jt,numxyz,num_threads,num_th;


# pragma omp parallel
{
num_threads=omp_get_num_threads();
}


if(num_threads==1){

	for(it=0;it<vecsize;it++){
	nocoul[it]=0;
	}


	for(it=0;it<nD;it++){
				
		    if((!EDGE1[it])&&(!EDGE4[it])){
			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];   

		    }
		    else if(EDGE4[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);				
			ftot=ftot+FX[it];	
			}
			else{
			LIST_VX[it] = viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(viti);	
			ftot2=ftot2+FX[it];	
			}
            

           		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
			
			Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
			Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
		
			FX[it] = 0.;
			FY[it] = 0.;	
			FZ[it] = 0.;
				 
			MTX[it] = 0.;
			MTY[it] = 0.;   
			MTZ[it] = 0.;   

			NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
			NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
			NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
			
			if(NX<0) NX=0;
			if(NX>=vecsizex) NX=vecsizex-1;		
			if(NY<0) NY=0;
			if(NY>=vecsizey) NY=vecsizey-1;		
			if(NZ<0) NZ=0;
			if(NZ>=vecsizez) NZ=vecsizez-1;		
			
			numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
			LIST_C[it] = numxyz;   
			coul[numxyz][nocoul[numxyz]]=it;
			nocoul[numxyz]++;

 

	}
   	    


}else{

	int coul2[vecsize][num_threads-1][15]; 
        int nocoul2[vecsize][num_threads-1];


	for(it=0;it<vecsize;it++){
		nocoul[it]=0;
		for(jt=0;jt<num_threads-1;jt++){
		nocoul2[it][jt]=0;
		}
	}


# pragma omp parallel for schedule(static) private(it,NX,NY,NZ,numxyz,num_th) reduction(+:Ec,ftot,ftot2) 
	for(it=0;it<nD;it++){
				
	num_th=omp_get_thread_num()-1;

		    if((!EDGE1[it])&&(!EDGE4[it])){
			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];   

		    }
		    else if(EDGE4[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);				
			ftot=ftot+FX[it];	
			}
			else{
			LIST_VX[it] = viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(viti);	
			ftot2=ftot2+FX[it];	
			}
            

           		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
			
			Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
			Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
		
			FX[it] = 0.;
			FY[it] = 0.;	
			FZ[it] = 0.;
				 
			MTX[it] = 0.;
			MTY[it] = 0.;   
			MTZ[it] = 0.;   

			NX= ((int) ((LIST_X[it]-H_POS)*vecsizex/H_TOT));
			NY= ((int) ((LIST_Y[it]-V_POS)*vecsizey/V_TOT));
			NZ= ((int) ((LIST_Z[it]-Z_POS)*vecsizez/Z_TOT));		
			
			if(NX<0) NX=0;
			if(NX>=vecsizex) NX=vecsizex-1;		
			if(NY<0) NY=0;
			if(NY>=vecsizey) NY=vecsizey-1;		
			if(NZ<0) NZ=0;
			if(NZ>=vecsizez) NZ=vecsizez-1;		
			
			numxyz=NX+NY*vecsizex+NZ*vecsizex*vecsizey;
			LIST_C[it] = numxyz;   

			if(num_th==-1){
			coul[numxyz][nocoul[numxyz]]=it;
			nocoul[numxyz]++;
			}else{
			coul2[numxyz][num_th][nocoul2[numxyz][num_th]]=it;
			nocoul2[numxyz][num_th]++;
			}

	}
   	    

	for(it=0;it<vecsize;it++){
		for(jt=0;jt<num_threads-1;jt++){
			for(int kt=0;kt<nocoul2[it][jt];kt++){
			coul[it][nocoul[it]]=coul2[it][jt][kt];
			nocoul[it]++;
			}                   
		}

	}
}








   	  


}




void verlet_brazil_qs(int NBENREG, R & ftot,R & ftot2, R epsi, R & Ec, R & Epp,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, R fimp, int nitf)
{

ftot=0.;
ftot2=0.;
			
Epp=0.;
int it;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec,ftot,ftot2) 
	for(it=0;it<nD;it++){
				
		    if((!EDGE1[it])&&(!EDGE4[it])){
			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];   

		    }
		    else if(EDGE4[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);				
			ftot=ftot+FX[it];	
			}
			else{
			LIST_VX[it] = viti;	
			LIST_AX[it] = FX[it]/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(viti);	
			ftot2=ftot2+FX[it];	
			}
            

		        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
			LIST_AY[it] = FY[it]/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
			LIST_AZ[it] = FZ[it]/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	

			LIST_AWX[it] = MTX[it]/LIST_I[it];
			LIST_AWY[it] = MTY[it]/LIST_I[it];	 
			LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
			
			Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
			Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
	
			FX[it] = 0.;
			FY[it] = 0.;	
			FZ[it] = 0.;
				 
			MTX[it] = 0.;
			MTY[it] = 0.;   
			MTZ[it] = 0.;   

	}
   	    
   	  
}



void verlet_qs_incr(R & Ec, R & Epp, bool &bout,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XA, R * LIST_YA, R * LIST_ZA,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_VXA, R * LIST_VYA, R * LIST_VZA, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ,R * FIX, R * FIY, R * FIZ, R * MTIX, R * MTIY, R * MTIZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, int NBENREG)
{
R ftot=0.;
R dplty=0.; 
R dpltz=0.; 
int ndp=0;
R Ec1=0.;
R Ep1=0.;
int it;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec1,ndp,dpltz,dplty,ftot) reduction(-:Ep1)
	for(it=0;it<nD;it++){

		LIST_XA[it]=LIST_X[it];
		LIST_YA[it]=LIST_Y[it];
		LIST_ZA[it]=LIST_Z[it];  

		LIST_VXA[it]=LIST_VX[it];
		LIST_VYA[it]=LIST_VY[it];
		LIST_VZA[it]=LIST_VZ[it];  

		LIST_TXA[it]=LIST_TX[it];
		LIST_TYA[it]=LIST_TY[it];
		LIST_TZA[it]=LIST_TZ[it];   	


		if((!EDGE1[it])&&(!EDGE4[it])) {	
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+(FIX[it]+FX[it])/LIST_M[it]);		
		LIST_AX[it] = (FIX[it]+FX[it])/LIST_M[it];
		LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
		}
		else if(EDGE4[it])
		{
		LIST_VX[it] = viti;
		LIST_AX[it] = (FIX[it]+FX[it])/LIST_M[it];
		LIST_X[it]  = LIST_X[it]+dt*viti;	
		}
		else{
		ftot+=(FIX[it]+FX[it]);
		}


		if(!EDGE2[it]) {	
		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+(FIY[it]+FY[it])/LIST_M[it]);
		LIST_AY[it] = (FIY[it]+FY[it])/LIST_M[it];
		LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		}

		if(!EDGE3[it]) {	    
		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FIZ[it]+FZ[it])/LIST_M[it]);    
		LIST_AZ[it] = (FIZ[it]+FZ[it])/LIST_M[it];	
		LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
		}

		if(EDGE6[it]&&EDGE5[it]){
		dpltz+= (LIST_Z[it]-LIST_ZO[it]);
		dplty+= (LIST_Y[it]-LIST_YO[it]);
		ndp+=1;
		}
    
       
		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+(MTIX[it]+MTX[it])/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+(MTIY[it]+MTY[it])/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+(MTIZ[it]+MTZ[it])/LIST_I[it]);	

		LIST_AWX[it] = (MTIX[it]+MTX[it])/LIST_I[it];
		LIST_AWY[it] = (MTIY[it]+MTY[it])/LIST_I[it];	 
		LIST_AWZ[it] = (MTIZ[it]+MTZ[it])/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

		Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

		Ep1-=0.*(LIST_X[it]-LIST_XO[it]);
		Ep1-=0.*(LIST_Y[it]-LIST_YO[it]);
		Ep1-=0.*(LIST_Z[it]-LIST_ZO[it]);         

		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   
}

 Ec=Ec1;
 Epp=Ep1;

 dplty/=ndp;
 dpltz/=ndp;

  // Eval nu
  
  R nuy=dplty/((ite-1)*dt*viti);
  R nuz=dpltz/((ite-1)*dt*viti);


  // Eval E

		 R EE,sigxx,epsx;	     
		 sigxx=ftot/(Z_TOT*V_TOT);
		 epsx=(ite-1)*dt*viti/H_TOT;

                 EE=sigxx/epsx; 
                 if(ite>int(nD/10.)){
					 Emoy=(Emoy*(ite-int(nD/10.)-1)+EE)/(ite-int(nD/10.));
				 } 
         	  		
	     if(ite%NBENREG==0){ 
			 
  	if((ite>int(nD/10.))&&(abs((Emoy1-Emoy)/Emoy)<1e-3)) bout=1;					 
				Emoy1=Emoy;


              cout<<"Defx "<<epsx<<", Sigxx :"<<sigxx<<", Module de Young inst :"<<EE<<", Module de Young moy :"<<Emoy<<endl;	              
                cout<<"nu : "<<nuy<<", "<<nuz<<endl;
                nuu=-(nuy+nuz)/2.;

			  	
		}       

}



void verlet_incr(R & Ec, R & Epp, bool &bout,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT,R H_POS,R V_POS,R Z_POS, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XA, R * LIST_YA, R * LIST_ZA,R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_VXA, R * LIST_VYA, R * LIST_VZA, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ,R * FIX, R * FIY, R * FIZ, R * MTIX, R * MTIY, R * MTIZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6,bool * EDGEC, int NBENREG)
{

R dep1=0.;
R dep4=0.;
int nd1,nd4;
nd1=0;nd4=0;
R dx,dy;
R vx,vy;

int it;

# pragma omp parallel for schedule(static) private(it,dx,dy,vx,vy) reduction(+:nd4,nd1,dep4,dep1)
	for(it=0;it<nD;it++){

		LIST_XA[it]=LIST_X[it];
		LIST_YA[it]=LIST_Y[it];
		LIST_ZA[it]=LIST_Z[it];  

		LIST_VXA[it]=LIST_VX[it];
		LIST_VYA[it]=LIST_VY[it];
		LIST_VZA[it]=LIST_VZ[it];  

		LIST_TXA[it]=LIST_TX[it];
		LIST_TYA[it]=LIST_TY[it];
		LIST_TZA[it]=LIST_TZ[it];   	

		if(EDGEC[it]==0){
		LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+(FX[it]+FIX[it])/LIST_M[it]);		
		LIST_AX[it] = (FX[it]+FIX[it])/LIST_M[it];
		LIST_X[it]  = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];

		LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+(FY[it]+FIY[it])/LIST_M[it]-0.);
		LIST_AY[it] = (FY[it]+FIY[it])/LIST_M[it]-0.;
		LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		}else{
               

                dx=(LIST_X[it]-0.005);//-viti*dt*ite;
                dy=(LIST_Y[it]-0.005);    

                vx=viti*dy/0.005;//+viti;
		vy=-viti*dx/0.005;

		LIST_VX[it] = vx;
		LIST_X[it]  = LIST_X[it]+dt*vx;

		LIST_VY[it] = vy;
		LIST_Y[it]  = LIST_Y[it]+dt*vy;
		}


		dep1+=LIST_Y[it];
		dep4+=LIST_X[it];

		LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FZ[it]+FIZ[it])/LIST_M[it]);    
		LIST_AZ[it] = (FZ[it]+FIZ[it])/LIST_M[it];	
		LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

		LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+(MTX[it]+MTIX[it])/LIST_I[it]);
		LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+(MTY[it]+MTIY[it])/LIST_I[it]);	 
		LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+(MTZ[it]+MTIZ[it])/LIST_I[it]);	

		LIST_AWX[it] = (MTX[it]+MTIX[it])/LIST_I[it];
		LIST_AWY[it] = (MTY[it]+MTIY[it])/LIST_I[it];	 
		LIST_AWZ[it] = (MTZ[it]+MTIZ[it])/LIST_I[it];	   	 

		LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
		LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
		LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 

		FX[it] = 0.;
		FY[it] = 0.;	
		FZ[it] = 0.;
			 
		MTX[it] = 0.;
		MTY[it] = 0.;   
		MTZ[it] = 0.;   

	}




	dep1/=nD;
	dep4/=nD;

	if(ite%NBENREG==0){ 
		cout<<"--- ite :"<<ite<< "  --  " << (dep4)<< "  --  " << (dep1) <<endl; 
	}     


}


void verlet_brazil_qs_incr(int NBENREG, R & ftot,R & ftot2, R epsi, R & Ec, R & Epp,R & nuu, R &Emoy,R &Emoy1, int ite, R viti,R H_TOT,R V_TOT,R Z_TOT, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_XA, R * LIST_YA, R * LIST_ZA, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * LIST_TX, R * LIST_TY, R * LIST_TZ,R * LIST_TXA, R * LIST_TYA, R * LIST_TZA, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_VXA, R * LIST_VYA, R * LIST_VZA, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ,R * FIX, R * FIY, R * FIZ, R * MTIX, R * MTIY, R * MTIZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6, R fimp, int nitf)
{

			
int it;
ftot=0.;
ftot2=0.;

R dx=0.;
int n4=0;

R dxb=0.;
int n4b=0;

R Ec1=0.;
R Ep1=0.;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec1,ftot,ftot2,dx,n4,dxb,n4b) reduction(-:Ep1)
	for(it=0;it<nD;it++){
			
		LIST_XA[it]=LIST_X[it];
		LIST_YA[it]=LIST_Y[it];
		LIST_ZA[it]=LIST_Z[it];  

		LIST_VXA[it]=LIST_VX[it];
		LIST_VYA[it]=LIST_VY[it];
		LIST_VZA[it]=LIST_VZ[it];  

		LIST_TXA[it]=LIST_TX[it];
		LIST_TYA[it]=LIST_TY[it];
		LIST_TZA[it]=LIST_TZ[it];   	

	
		    if((!EDGE1[it])&&(!EDGE4[it])){
			LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+(FIX[it]+FX[it])/LIST_M[it]);		
			LIST_AX[it] = (FIX[it]+FX[it])/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];   

		    }
		    else if(EDGE4[it]){				
			LIST_VX[it] = -viti;	
			LIST_AX[it] =(FIX[it]+FX[it])/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(-viti);				
			ftot=ftot+(FIX[it]+FX[it]);	
	                dx+=(LIST_X[it]-LIST_XA[it]);
                        n4++;
			}
			else{
			LIST_VX[it] = viti;	
			LIST_AX[it] = (FIX[it]+FX[it])/LIST_M[it];
			LIST_X[it] = LIST_X[it]+dt*(viti);	
			ftot2=ftot2+(FIX[it]+FX[it]);	
	                dxb+=(LIST_X[it]-LIST_XA[it]);
                        n4b++;
			}
            

		        LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+(FIY[it]+FY[it])/LIST_M[it]);
			LIST_AY[it] = (FIY[it]+FY[it])/LIST_M[it];
			LIST_Y[it]  = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];
		    		
			LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+(FIZ[it]+FZ[it])/LIST_M[it]);    
			LIST_AZ[it] = (FIZ[it]+FZ[it])/LIST_M[it];	
			LIST_Z[it]  = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	

			LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+(MTIX[it]+MTX[it])/LIST_I[it]);
			LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+(MTIY[it]+MTY[it])/LIST_I[it]);	 
			LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+(MTIZ[it]+MTZ[it])/LIST_I[it]);	

			LIST_AWX[it] = (MTIX[it]+MTX[it])/LIST_I[it];
			LIST_AWY[it] = (MTIY[it]+MTY[it])/LIST_I[it];	 
			LIST_AWZ[it] = (MTIZ[it]+MTZ[it])/LIST_I[it];	   	 

			LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
			LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
			LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
			
		Ec1+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
		Ec1+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);

		Ep1-=(FIX[it]+FX[it])*(LIST_X[it]-LIST_XA[it]);
		Ep1-=(FIY[it]+FY[it])*(LIST_Y[it]-LIST_YA[it]);
	 	Ep1-=(FIZ[it]+FZ[it])*(LIST_Z[it]-LIST_ZA[it]);         
	  	Ep1-=(MTIX[it]+MTX[it])*(LIST_TX[it]-LIST_TXA[it])/2.;
		Ep1-=(MTIY[it]+MTY[it])*(LIST_TY[it]-LIST_TYA[it])/2.;
	 	Ep1-=(MTIZ[it]+MTZ[it])*(LIST_TZ[it]-LIST_TZA[it])/2.;    

			FX[it] = 0.;
			FY[it] = 0.;	
			FZ[it] = 0.;
				 
			MTX[it] = 0.;
			MTY[it] = 0.;   
			MTZ[it] = 0.;   

	}

 Ec=Ec1;
 Epp+=Ep1;

//Ec=(ite-1)*dt*viti*(fabs(ftot)+fabs(ftot2))/2.;

   	    
   	  
}


void verlet_symt_qs_DCB(R & ftot,R & ftot2,R & Ec, R viti, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6,bool * EDGE1M1,bool * EDGE1M2)
{
ftot=0.;
ftot2=0.;
Ec=0.;
int it,num_threads;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec,ftot,ftot2) 
	for(it=0;it<nD;it++){

	if((!EDGE1[it])&&(!EDGE4[it])) {	
		
	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	LIST_AZ[it] = FZ[it]/LIST_M[it];	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
		
	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
    
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
	
	}
    else if(EDGE1M1[it]) {		
	LIST_VZ[it] = -viti;
	LIST_AZ[it] = FZ[it]/LIST_M[it];
	LIST_Z[it]  = LIST_Z[it]-dt*viti;
	ftot+=FZ[it];	
	}
    else if(EDGE1M2[it]) {		
	LIST_VZ[it] = viti;
	LIST_AZ[it] = FZ[it]/LIST_M[it];
	LIST_Z[it]  = LIST_Z[it]+dt*viti;
	ftot2+=FZ[it];	
	}
    

  	
	 LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	 LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	 LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	
	 
     LIST_AWX[it] = MTX[it]/LIST_I[it];
     LIST_AWY[it] = MTY[it]/LIST_I[it];	 
     LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

     LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
     LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
     LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
     
	 Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	 Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
                    
	 FX[it] = 0.;
	 FY[it] = 0.;	
	 FZ[it] = 0.;
	 		 
     MTX[it] = 0.;
     MTY[it] = 0.;   
     MTZ[it] = 0.;   



}
	

}

void verlet_symt_qs_ELS(R & ftot,R & ftot2,R & Ec, R viti, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6,bool * EDGE1M1,bool * EDGE1M2)
{
ftot=0.;
ftot2=0.;
Ec=0.;
int it,num_threads;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec,ftot,ftot2) 
	for(it=0;it<nD;it++){

	if((!EDGE1M1[it])&&(!EDGE4[it])) {	
		
	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	LIST_AZ[it] = FZ[it]/LIST_M[it];	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
		
	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
    
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
	
	}
    else if(EDGE1M1[it]) {		
	LIST_VZ[it] = viti;
	LIST_AZ[it] = FZ[it]/LIST_M[it];
	LIST_Z[it]  = LIST_Z[it]+dt*viti;
		ftot2+=FZ[it];	
	}
	else  //encastrement
	{
	ftot+=FZ[it];	
    }
    

  	
	 LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	 LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	 LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	
	 
     LIST_AWX[it] = MTX[it]/LIST_I[it];
     LIST_AWY[it] = MTY[it]/LIST_I[it];	 
     LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

     LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
     LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
     LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
     
	 Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	 Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
                    
	 FX[it] = 0.;
	 FY[it] = 0.;	
	 FZ[it] = 0.;
	 		 
     MTX[it] = 0.;
     MTY[it] = 0.;   
     MTZ[it] = 0.;   



}
	

}

void verlet_symt_qs_MMB(R & ftot,R & ftot2,R & Ec, R viti, int nD,R dt,R * LIST_M,R * LIST_I, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ, R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * LIST_AX, R * LIST_AY, R * LIST_AZ, R * LIST_AWX, R * LIST_AWY, R * LIST_AWZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, bool * EDGE1,bool * EDGE2,bool * EDGE3,bool * EDGE4,bool * EDGE5,bool * EDGE6,bool * EDGE13,bool * EDGE43, bool * EDGE6X, bool * EDGE7X)
{
ftot=0.;
ftot2=0.;
Ec=0.;
int it,num_threads;

# pragma omp parallel for schedule(static) private(it) reduction(+:Ec,ftot,ftot2) 
	for(it=0;it<nD;it++){

	if((!EDGE13[it])&&(!EDGE43[it])&&(!EDGE6X[it])&&(!EDGE7X[it])) {	
		
	LIST_VZ[it] = LIST_VZ[it]+dt/2.*(LIST_AZ[it]+FZ[it]/LIST_M[it]);    
	LIST_AZ[it] = FZ[it]/LIST_M[it];	
	LIST_Z[it] = LIST_Z[it]+dt*LIST_VZ[it]+dt*dt/2.*LIST_AZ[it];	
		
	LIST_VY[it] = LIST_VY[it]+dt/2.*(LIST_AY[it]+FY[it]/LIST_M[it]);
	LIST_AY[it] = FY[it]/LIST_M[it];
	LIST_Y[it] = LIST_Y[it]+dt*LIST_VY[it]+dt*dt/2.*LIST_AY[it];	
    
	LIST_VX[it] = LIST_VX[it]+dt/2.*(LIST_AX[it]+FX[it]/LIST_M[it]);		
	LIST_AX[it] = FX[it]/LIST_M[it];
	LIST_X[it] = LIST_X[it]+dt*LIST_VX[it]+dt*dt/2.*LIST_AX[it];
	
	}
    else if(EDGE6X[it]) {		
	LIST_VZ[it] = -0.027*viti;
	LIST_AZ[it] = FZ[it]/LIST_M[it];
	LIST_Z[it]  = LIST_Z[it]-0.027*dt*viti;
	ftot2+=FZ[it];	
	}
    else if(EDGE7X[it]) {		
	LIST_VZ[it] = viti;
	LIST_AZ[it] = FZ[it]/LIST_M[it];
	LIST_Z[it]  = LIST_Z[it]+dt*viti;
	ftot+=FZ[it];	
	}
    

  	
	 LIST_WX[it]  = LIST_WX[it]+dt/2.*(LIST_AWX[it]+MTX[it]/LIST_I[it]);
	 LIST_WY[it]  = LIST_WY[it]+dt/2.*(LIST_AWY[it]+MTY[it]/LIST_I[it]);	 
	 LIST_WZ[it]  = LIST_WZ[it]+dt/2.*(LIST_AWZ[it]+MTZ[it]/LIST_I[it]);	
	 
     LIST_AWX[it] = MTX[it]/LIST_I[it];
     LIST_AWY[it] = MTY[it]/LIST_I[it];	 
     LIST_AWZ[it] = MTZ[it]/LIST_I[it];	   	 

     LIST_TX[it] = LIST_TX[it]+dt*LIST_WX[it]+dt*dt/2.*LIST_AWX[it];
     LIST_TY[it] = LIST_TY[it]+dt*LIST_WY[it]+dt*dt/2.*LIST_AWY[it];     
     LIST_TZ[it] = LIST_TZ[it]+dt*LIST_WZ[it]+dt*dt/2.*LIST_AWZ[it]; 
     
	 Ec+=0.5*LIST_M[it]*(LIST_VX[it]*LIST_VX[it]+LIST_VY[it]*LIST_VY[it]+LIST_VZ[it]*LIST_VZ[it]);
	 Ec+=0.5*LIST_I[it]*(LIST_WX[it]*LIST_WX[it]+LIST_WY[it]*LIST_WY[it]+LIST_WZ[it]*LIST_WZ[it]);
                    
	 FX[it] = 0.;
	 FY[it] = 0.;	
	 FZ[it] = 0.;
	 		 
     MTX[it] = 0.;
     MTY[it] = 0.;   
     MTZ[it] = 0.;   



}
	

}


