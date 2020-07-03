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

#include "rupture.h"

void rupture1_bond_Rankine(int NBCO, int ** CONT, bool * TYPCO, R ** FOJI, R ** MTJI, R ** MTIJ, R * LIST_R, R rmu, R siglim, int * LIST_B, bool & brupt){
	
R rankinemax=-1e8;
int numc1,numc2;
R rankine,sigfl,sigtr,sigmax,taumax;
R ray1,ray2,ray,S,I;	

int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,numc1,numc2,rankine,sigfl,sigtr,sigmax,taumax,ray1,ray2,ray,S,I) reduction(max:rankinemax)
	for(it=0;it<NBCO;it++){ 		
    
		if(TYPCO[it]==1){

		numc1=CONT[it][0];
		numc2=CONT[it][1];
    		
		ray1=LIST_R[numc1];
		ray2=LIST_R[numc2];   
		ray=rmu*(ray1+ray2)/2.;
		S=3.14159*ray*ray;
		I=3.14159*pow(ray,4)/4.;
		
		sigfl=sqrt((MTJI[it][1]+MTIJ[it][1])*(MTJI[it][1]+MTIJ[it][1])+(MTJI[it][2]+MTIJ[it][2])*(MTJI[it][2]+MTIJ[it][2]));
		sigfl=sigfl*ray/I;		
		sigtr=FOJI[it][0]/S;		
		sigmax=sigfl+sigtr;	 //> 0?		
	//	sigmax=(sigmax>0.)?sigmax:0.;		
		taumax=	4.*sqrt(FOJI[it][1]*FOJI[it][1]+FOJI[it][2]*FOJI[it][2])/(3.*S);
		
		rankine=0.5*(sigmax+sqrt(sigmax*sigmax+4.*taumax*taumax));					
		if(rankine>siglim) {brupt=1;TYPCO[it]=0; LIST_B[numc1]=2;LIST_B[numc2]=2; /*cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lien rompu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;*/	}
		rankinemax=max(rankine,rankinemax);		
		}   
    		
	}
	
	//cout<<"Contrainte max :"<<rankinemax<<endl;
	
}

void rupture1_bond(int NBCO, int ** CONT, bool * TYPCO, R ** FOJI, R ** MTJI, R ** MTIJ, R * LIST_R, R rmu, R siglimt,R siglimcis, int * LIST_B, bool & brupt){
	
int it,num_threads;
R sigfl,sigtr,sigmax,taumax;
int numc1,numc2;
R ray,ray1,ray2,S;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,sigfl,sigtr,sigmax,taumax,numc1,numc2,ray,ray1,ray2,S)
	for(it=0;it<NBCO;it++){ 		
    
		if(TYPCO[it]==1){

		numc1=CONT[it][0];
		numc2=CONT[it][1];
    		
		ray1=LIST_R[numc1];
		ray2=LIST_R[numc2];   
		ray=rmu*(ray1+ray2)/2.;
		S=3.14159*ray*ray;
		
		sigfl=0.;
		sigtr=FOJI[it][0]/S;		
		sigmax=sigfl+sigtr;				
		taumax=	sqrt(FOJI[it][1]*FOJI[it][1]+FOJI[it][2]*FOJI[it][2])/S;
							
		if((sigmax>siglimt)||(taumax>siglimcis)) { brupt=1;TYPCO[it]=0; LIST_B[numc1]=2;LIST_B[numc2]=2; 	/*	cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lien rompu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	*/}
		
		}   
    		
	}
	
}

 
void rupture1a_trace(int NBCO, bool * TYPCO, R * TRACE, int ** CONT, R siglimt, R siglimc){

int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
	
int numc1,numc2;
R trac;
# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,numc1,numc2,trac)
	for(it=0;it<NBCO;it++){ 
		
	    numc1=CONT[it][0];
	    numc2=CONT[it][1];
	    
	    trac=sqrt(TRACE[numc1]*TRACE[numc2]);		
	    		
	    if((trac>siglimt)||(trac<-siglimc)) {  TYPCO[it]=0; }		
    		
	}
	
}

void rupture1b_trace(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * TRACE, R siglimt, R siglimc){
	
int it,jt,num_threads;
R trac;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,trac)
	for(it=0;it<NB_SPH;it++){ 
		trac=TRACE[it];		
				
		if((trac>siglimt)||(trac<-siglimc)) { 
			for(jt=0;jt<NBCONTCO[it];jt++){ 
				TYPCO[NOCONT[it][jt]]=0;
			}
		}		
    		
	}
	
}

void rupture1c_trace(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * TRACE, R * LIST_X, R * LIST_Y, R * LIST_Z, R siglimt, R siglimc, int * LIST_B, bool & brupt,vector<int> & list_rupt){
	
int it,jt,num_threads;
R trac;
int nurupt=0;

	brupt=0;  
	list_rupt.clear();

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,trac)
	for(it=0;it<NB_SPH;it++){ 
		trac=TRACE[it];		
				
		if((trac>siglimt)||(trac<-siglimc)) { 
		   
		   LIST_X[it]=LIST_X[it]*1e10;
		   LIST_Y[it]=LIST_Y[it]*1e10;	
		   LIST_Z[it]=LIST_Z[it]*1e10;
		   LIST_B[it]=2;
		# pragma omp critical
		{
                           brupt=1;nurupt++;list_rupt.push_back(it);
                }	
		   	//	cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee Trac !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
			for(jt=0;jt<NBCONTCO[it];jt++){ 
				TYPCO[NOCONT[it][jt]]=0;
			}
		}		
    		
	}
	
} 


// Rupture Griffith
void rupture1_Griffith(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt, int * LIST_B, bool & brupt,int &nrt,int &nrcis,int &nrtot,vector<int> & list_rupt){

	brupt=0;  
	list_rupt.clear();

int nurupt=0;

R s1,s3;
int it,jt,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,s1,s3) reduction(+:nrt,nrtot,nrcis)

	for(it=0;it<NB_SPH;it++){ 

        s1=SIG1[it];
        s3=SIG3[it];            
        		
		if((3*s1+s3)>=0.){
				
			if(s1>siglimt) { 
				
				nrt++;
				nrtot++;
			   
			   LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;	
		# pragma omp critical
		{
                           brupt=1;nurupt++;list_rupt.push_back(it);
                }
				//   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
				}
				
			}	
    		
      }
      else if(((s1-s3)*(s1-s3)+(8.*siglimt*(s1+s3)))>=0){	
		  
		  		nrcis++;
				nrtot++;
				
		       LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;	

			# pragma omp critical
			{
		                   brupt=1;nurupt++;list_rupt.push_back(it);
		        }
				//	cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
				}
		  
	  }
           
      		
	}


}

// Rupture Griffith incr
void rupture1_Griffith_incr(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt, int * LIST_B, bool & brupt,int &nrt,int &nrcis,int &nrtot,vector<int> & list_rupt,R ** FCJI){

	brupt=0;  
	list_rupt.clear();

int nurupt=0;

R s1,s3;
int it,jt,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,s1,s3) reduction(+:nrt,nrtot,nrcis)

	for(it=0;it<NB_SPH;it++){ 

        s1=SIG1[it];
        s3=SIG3[it];            
        		
		if((3*s1+s3)>=0.){
				
			if(s1>siglimt) { 
				
				nrt++;
				nrtot++;
			   
			   LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;	
		# pragma omp critical
		{
                           brupt=1;nurupt++;list_rupt.push_back(it);
                }
				//   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
					FCJI[NOCONT[it][jt]][0]=0.;
					FCJI[NOCONT[it][jt]][1]=0.;
					FCJI[NOCONT[it][jt]][2]=0.;
				}
				
			}	
    		
      }
      else if(((s1-s3)*(s1-s3)+(8.*siglimt*(s1+s3)))>=0){	
		  
		  		nrcis++;
				nrtot++;
				
		       LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;	

			# pragma omp critical
			{
		                   brupt=1;nurupt++;list_rupt.push_back(it);
		        }
				//	cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
					FCJI[NOCONT[it][jt]][0]=0.;
					FCJI[NOCONT[it][jt]][1]=0.;
					FCJI[NOCONT[it][jt]][2]=0.;
				}
		  
	  }
           
      		
	}


}


// Rupture Tresca
void rupture1_Tresca(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimcis, int * LIST_B, bool & brupt,int &nrcis,int &nrtot,vector<int> & list_rupt){

	brupt=0;  
	list_rupt.clear();

R s1,s3;
int it,jt,num_threads;
int nurupt=0;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,s1,s3) reduction(+:nrtot,nrcis)

	for(it=0;it<NB_SPH;it++){ 

        s1=SIG1[it];
        s3=SIG3[it];            
        		
		if((s1-s3)/2.>=siglimcis){
				
				nrcis++;
				nrtot++;
			   
			   LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;
				# pragma omp critical
				{
					   brupt=1;nurupt++;list_rupt.push_back(it);
				}	
				//   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee cis !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
				}
	
	  }
           
      		
	}
}

// Rupture Mohr_Coulomb
void rupture1_Mohr_Coulomb(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimcis, R phi, int * LIST_B, bool & brupt,int &nrcis,int &nrtot,vector<int> & list_rupt){

	brupt=0;  
	list_rupt.clear();

R s1,s3,yield;
int it,jt,num_threads;
int nurupt=0;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,s1,s3,yield) reduction(+:nrtot,nrcis)
	for(it=0;it<NB_SPH;it++){ 

        s1=SIG1[it];
        s3=SIG3[it];            
        
        yield=sin(phi)*(s1+s3)+(s1-s3)-2*siglimcis*cos(phi);		
        		
		if(yield>=0){
				
				nrcis++;
				nrtot++;
						   
				   LIST_X[it]=LIST_X[it]*1e10;
				   LIST_Y[it]=LIST_Y[it]*1e10;	
				   LIST_Z[it]=LIST_Z[it]*1e10;
				   LIST_B[it]=2;	
				# pragma omp critical
				{
					   brupt=1;nurupt++;list_rupt.push_back(it);
				}

			//	   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee MohrC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
				}
	
	    }
           
      		
	}
}


// Rupture Rankine
void rupture1_Rankine(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt, R siglimc, int * LIST_B, bool & brupt,int &nrt,int &nrc,int &nrtot,vector<int> & list_rupt){

	brupt=0;  
	list_rupt.clear();

R s1,s3;
int it,jt,num_threads;
int nurupt=0;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,s1,s3) reduction(+:nrtot,nrc,nrt)
	for(it=0;it<NB_SPH;it++){ 

        s1=SIG1[it];   
        s3=SIG3[it];      
        		
		if(s1>=siglimt){
				
				nrt++;
				nrtot++;
						   
			   LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;	

				# pragma omp critical
				{
					   brupt=1;nurupt++;list_rupt.push_back(it);
				}

			//	   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee trac !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
				}
	
	    }
  		else if(s3<=-siglimc){
			
				nrc++;
				nrtot++;		
			   
			   LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;	

				# pragma omp critical
				{
				   brupt=1;nurupt++;list_rupt.push_back(it);
				}

			//	   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee comp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
				}
	
	    }         
      		
	}
}

void rupture2a_trace(int NBCO, bool * TYPCO, R * TRACE, int ** CONT, R siglimt1, R siglimt2, R siglimti, bool * LIST_P){
	
int numc1,numc2;
int it,num_threads;
R trac;
R siglimt;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)) private(it,numc1,numc2,trac,siglimt)
	for(it=0;it<NBCO;it++){ 
		
	    numc1=CONT[it][0];
	    numc2=CONT[it][1];
		
		if((LIST_P[numc1]==0)&&(LIST_P[numc2]==0)) siglimt=siglimt1;
		if((LIST_P[numc1]==0)&&(LIST_P[numc2]==1)) siglimt=siglimti;
		if((LIST_P[numc1]==1)&&(LIST_P[numc2]==0)) siglimt=siglimti;		
		if((LIST_P[numc1]==1)&&(LIST_P[numc2]==1)) siglimt=siglimt2;		 
    
    trac=sqrt(TRACE[numc1]*TRACE[numc2]);		
    		
    if(trac>siglimt) {  TYPCO[it]=0; }		
    		
	}
	
}

void rupture2b_trace(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * TRACE,R siglimt1, R siglimt2, bool * LIST_P,int * LIST_B, bool & brupt,vector<int> & list_rupt){
	
	brupt=0;  
	list_rupt.clear();

int it,jt,num_threads;
R trac;
R siglimt;
int nurupt=0;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,trac,siglimt)
	for(it=0;it<NB_SPH;it++){ 
		trac=TRACE[it];		
		siglimt;
		
		if(LIST_P[it]==0) siglimt=siglimt1;
		if(LIST_P[it]==1) siglimt=siglimt2;
						
		if(trac>siglimt) { 
			   LIST_B[it]=2;
				# pragma omp critical
				{
				   brupt=1;nurupt++;list_rupt.push_back(it);
				}			
			for(jt=0;jt<NBCONTCO[it];jt++){ 
				TYPCO[NOCONT[it][jt]]=0;
			}
		}		
    		
	}
	
}

void rupture2c_trace(int NB_SPH, unsigned int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * TRACE, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt1, R siglimt2,bool * LIST_P,int * LIST_B, bool & brupt,vector<int> & list_rupt){

	brupt=0;  
	list_rupt.clear();

int it,jt,num_threads;
R trac;
R siglimt;
int nurupt=0;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,trac,siglimt)
	for(it=0;it<NB_SPH;it++){ 
		trac=TRACE[it];	
		
		if(LIST_P[it]==0) siglimt=siglimt1;
		if(LIST_P[it]==1) siglimt=siglimt2;
		
				
		if(trac>siglimt) { 
		   
		   LIST_X[it]=LIST_X[it]*1e10;
		   LIST_Y[it]=LIST_Y[it]*1e10;	
		   LIST_Z[it]=LIST_Z[it]*1e10;
		   LIST_B[it]=2;
				# pragma omp critical
				{
				   brupt=1;nurupt++;list_rupt.push_back(it);
				}		   
		   	   	
			for(jt=0;jt<NBCONTCO[it];jt++){ 
				TYPCO[NOCONT[it][jt]]=0;
			}
		}		
    		
	}
	
} 








void rupture2_Rankine(int NB_SPH, int ** NOCONT, int * NBCONTCO, bool * TYPCO, R * SIG1, R * SIG3, R * LIST_X, R * LIST_Y,R * LIST_Z, R siglimt1, R siglimt2, R siglimc1, R siglimc2, int * LIST_B, bool & brupt,int &nrt,int &nrc,int &nrtot,vector<int> & list_rupt, bool * LIST_P){

	brupt=0;  
	list_rupt.clear();

R s1,s3;
int it,jt,num_threads;
int nurupt=0;
R siglimt,siglimc;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it,jt,s1,s3) reduction(+:nrtot,nrc,nrt)
	for(it=0;it<NB_SPH;it++){ 

		if(LIST_P[it]==0) {siglimt=siglimt1;siglimc=siglimc1;}
		if(LIST_P[it]==1) {siglimt=siglimt2;siglimc=siglimc2;}

        s1=SIG1[it];   
        s3=SIG3[it];      
        		
		if(s1>=siglimt){
				
				nrt++;
				nrtot++;
						   
			   LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;	

				# pragma omp critical
				{
					   brupt=1;nurupt++;list_rupt.push_back(it);
				}

		//		   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee trac !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
				}
	
	    }
  		else if(s3<=-siglimc){
			
				nrc++;
				nrtot++;		
			   
			   LIST_X[it]=LIST_X[it]*1e10;
			   LIST_Y[it]=LIST_Y[it]*1e10;	
			   LIST_Z[it]=LIST_Z[it]*1e10;
			   LIST_B[it]=2;	

				# pragma omp critical
				{
				   brupt=1;nurupt++;list_rupt.push_back(it);
				}

		//		   cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!! Particule ejectee comp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;	   		   	   	
				for(jt=0;jt<NBCONTCO[it];jt++){ 
					TYPCO[NOCONT[it][jt]]=0;
				}
	
	    }         
      		
	}
}




