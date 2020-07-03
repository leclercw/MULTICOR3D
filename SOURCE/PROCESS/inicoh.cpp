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

#include "inicoh.h"
#include "omp.h"


void inicoh2(bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, bool * LIST_P, R rmu1, R rmu2, R Emu1, R Emu2, R nuu1, R nuu2, R amort)
{

R dx,dy,dz,nd2,r2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,jt;
R Emui,nuui,rmui;      
          
dt=1e9;
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii,it,jt,fors,norms,kkn,kkt,dx,dy,dz,nd2,S,I,ray,ray1,ray2,Emui,nuui,rmui) reduction(min:dt)                
	for(ii=0;ii<nbco;ii++){

					   it=CONT[ii][0];
					   jt=CONT[ii][1];
					   
					   dx=LIST_X[it]-LIST_X[jt]; 
					   dy=LIST_Y[it]-LIST_Y[jt]; 
					   dz=LIST_Z[it]-LIST_Z[jt]; 					   
					   nd2=sqrt(dx*dx+dy*dy+dz*dz);	   
						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                         
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                      
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];      
						
						
					   if(TYPCO[ii]==1)
		               {
			        
                        
                       // Cohesion
 
						ray1=LIST_R[it];
						ray2=LIST_R[jt];   

               		     if((LIST_P[it]==0)&&(LIST_P[jt]==0)){	                                            
						 ray=rmu1*(ray1+ray2)/2.;
						 S=(Pi*ray*ray);
						 I=Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emu1*S/nd2;
						 kkt=12.*Emu1*I/(pow(nd2,3));   
					     
					 //    cout<<"knt:"<<kkn<<", "<<kkt<<endl;
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu1*I/((1+nuu1)*nd2) ;}
	               		 else if(((LIST_P[it]==0)&&(LIST_P[jt]==1))||((LIST_P[it]==1)&&(LIST_P[jt]==0))){	                                            
						 
				   	     rmui=(rmu1+rmu2)/2.;
                 		 Emui=(Emu1+Emu2)/2.; 
                 		 nuui=(nuu1+nuu2)/2.; 
                 		 					 
						 ray=rmui*(ray1+ray2)/2.;
						 S=(Pi*ray*ray);
						 I=Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emui*S/nd2;
						 kkt=12.*Emui*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emui*I/((1+nuui)*nd2) ;}				 
               		     else if((LIST_P[it]==1)&&(LIST_P[jt]==1)){	                                            
						 ray=rmu2*(ray1+ray2)/2.;
						 S=(Pi*ray*ray);
						 I=Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emu2*S/nd2;
						 kkt=12.*Emu2*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu2*I/((1+nuu2)*nd2) ;}
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt])));  
					 }
						 
                        
}
/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////
				
dt=Cs*dt;
cout<<"Pas de temps:"<<dt<<endl;
}



void inicoh2_int_mI(int & nint, int & nint0, R * vs, R * LIST_IND, R H_TOT,int * LIST_B,bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO,  unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, bool * LIST_P, R rmu1, R rmu2, R Emu1, R Emu2, R nuu1, R nuu2, R amort)
{

R dx,dy,dz,nd2,r2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,jt;
R Emui,nuui,rmui;  
R ind1,ind2,indm;

nint=0;
     
dt=1e9;
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
R Smoy=0.;
R indtot=0.;

# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii,dx,dy,dz,nd2,ind1,ind2,indm,r2,S,I,ray,ray1,ray2,kkn,kkt,fors,norms,it,jt,Emui,nuui,rmui) reduction(min:dt) reduction(+:nint,Smoy,indtot)                         
	for(ii=0;ii<nbco;ii++){
		
		vs[ii]=0.;
		
         //     cout<<"hello"<<endl;
					   it=CONT[ii][0];
					   jt=CONT[ii][1];
					   
					   dx=LIST_X[it]-LIST_X[jt]; 
					   dy=LIST_Y[it]-LIST_Y[jt]; 
					   dz=LIST_Z[it]-LIST_Z[jt]; 					   
					   nd2=sqrt(dx*dx+dy*dy+dz*dz);	   
						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                          
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                     
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];      
						
						
					   if(TYPCO[ii]==1)
		               {
			        
                        
                       // Cohesion
 
						ray1=LIST_R[it];
						ray2=LIST_R[jt];   

						ind1=LIST_IND[it];
						ind2=LIST_IND[jt];				                             
						indm=max(ind1,ind2);    

               		     if((LIST_P[it]==0)&&(LIST_P[jt]==0)){	            
							 							 					//	 cout<<"number1"<<endl;                               
						 ray=rmu1*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;

					     indtot+=indm;
						 Smoy+=S;
                         
						 kkn=Emu1*S/nd2;
						 kkt=12.*Emu1*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu1*I/((1+nuu1)*nd2) ;}
	               		 else if(((LIST_P[it]==0)&&(LIST_P[jt]==1))||((LIST_P[it]==1)&&(LIST_P[jt]==0))){	 
							 							 					//	 cout<<"number12"<<endl;                                           
						 
						 
						 //délaminage 25/02/2020
						     
					nint++;	 
						 
						 
				   	     rmui=(rmu1+rmu2)/2.;
                 		 Emui=(Emu1+Emu2)/2.; 
                 		 nuui=(nuu1+nuu2)/2.; 
                 		 					 
						 ray=rmui*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                        

			             indtot+=indm;
						 Smoy+=S;

						 //kkn=Emui*S/nd2;
						 kkn=Emui*S/nd2;
						 kkt=0.;//12.*Emui*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=0.;  
						 VALCOH[ii][2]=0. ;  	
						 VALCOH[ii][3]=0. ;  	
						 VALCOH[ii][4]=0. ;
						 VALCOH[ii][5]=0. ;  
						 
						 
						 }				 
               		     else if((LIST_P[it]==1)&&(LIST_P[jt]==1)){	         
							 						// cout<<"number2"<<endl;
						                                                  
						 ray=rmu2*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                         

						indtot+=indm;
						Smoy+=S;

						 kkn=Emu2*S/nd2;
						 kkt=12.*Emu2*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu2*I/((1+nuu2)*nd2) ;}
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt]))); 
						 
						 
						  
					 }
						 
                        
}
/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////	

dt=Cs*dt;
Smoy/=indtot;
nint0=nint;
cout<<"Surface moyenne:"<<Smoy<<endl;
cout<<"Pas de temps:"<<dt<<endl;
}




void inicoh2_int_mII(int & nint, int & nint0, R * vs, R * LIST_IND, R H_TOT,int * LIST_B,bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO,  unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, bool * LIST_P, R rmu1, R rmu2, R Emu1, R Emu2, R nuu1, R nuu2, R amort)
{

R dx,dy,dz,nd2,r2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,jt;
R Emui,nuui,rmui;  
R ind1,ind2,indm;

nint=0;
     
dt=1e9;
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
R Smoy=0.;
R indtot=0.;

# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii,dx,dy,dz,nd2,ind1,ind2,indm,r2,S,I,ray,ray1,ray2,kkn,kkt,fors,norms,it,jt,Emui,nuui,rmui) reduction(min:dt) reduction(+:nint,Smoy,indtot)                         
	for(ii=0;ii<nbco;ii++){
		
		vs[ii]=0.;
		
         //     cout<<"hello"<<endl;
					   it=CONT[ii][0];
					   jt=CONT[ii][1];
					   
					   dx=LIST_X[it]-LIST_X[jt]; 
					   dy=LIST_Y[it]-LIST_Y[jt]; 
					   dz=LIST_Z[it]-LIST_Z[jt]; 					   
					   nd2=sqrt(dx*dx+dy*dy+dz*dz);	   
						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                          
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                     
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];      
						
						
					   if(TYPCO[ii]==1)
		               {
			        
                        
                       // Cohesion
 
						ray1=LIST_R[it];
						ray2=LIST_R[jt];   

						ind1=LIST_IND[it];
						ind2=LIST_IND[jt];				                             
						indm=max(ind1,ind2);    

               		     if((LIST_P[it]==0)&&(LIST_P[jt]==0)){	            
							 							 					//	 cout<<"number1"<<endl;                               
						 ray=rmu1*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;

					     indtot+=indm;
						 Smoy+=S;
                         
						 kkn=Emu1*S/nd2;
						 kkt=12.*Emu1*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu1*I/((1+nuu1)*nd2) ;}
	               		 else if(((LIST_P[it]==0)&&(LIST_P[jt]==1))||((LIST_P[it]==1)&&(LIST_P[jt]==0))){	 
							 							 					//	 cout<<"number12"<<endl;                                           
						 
						 
						 //délaminage 25/02/2020
						     
					nint++;	 
						 
						 
				   	     rmui=(rmu1+rmu2)/2.;
                 		 Emui=(Emu1+Emu2)/2.; 
                 		 nuui=(nuu1+nuu2)/2.; 
                 		 					 
						 ray=rmui*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                        

			             indtot+=indm;
						 Smoy+=S;

						 kkn=0;
						 kkt=12.*Emui*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=0.;//kkn;
						 VALCOH[ii][1]=kkt;  
						 VALCOH[ii][2]=0. ;  	
						 VALCOH[ii][3]=0. ;  	
						 VALCOH[ii][4]=0. ;
						 VALCOH[ii][5]=0. ;  
						 
						 
						 }				 
               		     else if((LIST_P[it]==1)&&(LIST_P[jt]==1)){	         
							 						// cout<<"number2"<<endl;
						                                                  
						 ray=rmu2*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                         

						indtot+=indm;
						Smoy+=S;

						 kkn=Emu2*S/nd2;
						 kkt=12.*Emu2*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu2*I/((1+nuu2)*nd2) ;}
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt]))); 
						 
						 
						  
					 }
						 
                        
}
/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////	

dt=Cs*dt;
Smoy/=indtot;
nint0=nint;
cout<<"Surface moyenne:"<<Smoy<<endl;
cout<<"Pas de temps:"<<dt<<endl;
}





void inicoh2_int_mm(int & nint, int & nint0, R * vs, R * LIST_IND, R H_TOT,int * LIST_B,bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO,  unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, bool * LIST_P, R rmu1, R rmu2, R Emu1, R Emu2, R nuu1, R nuu2, R amort)
{

R dx,dy,dz,nd2,r2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,jt;
R Emui,nuui,rmui;  
R ind1,ind2,indm;

nint=0;
     
dt=1e9;
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
R Smoy=0.;
R indtot=0.;

# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii,dx,dy,dz,nd2,ind1,ind2,indm,r2,S,I,ray,ray1,ray2,kkn,kkt,fors,norms,it,jt,Emui,nuui,rmui) reduction(min:dt) reduction(+:nint,Smoy,indtot)                         
	for(ii=0;ii<nbco;ii++){
		
		vs[ii]=0.;
		
         //     cout<<"hello"<<endl;
					   it=CONT[ii][0];
					   jt=CONT[ii][1];
					   
					   dx=LIST_X[it]-LIST_X[jt]; 
					   dy=LIST_Y[it]-LIST_Y[jt]; 
					   dz=LIST_Z[it]-LIST_Z[jt]; 					   
					   nd2=sqrt(dx*dx+dy*dy+dz*dz);	   
						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                          
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                     
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];      
						
						
					   if(TYPCO[ii]==1)
		               {
			        
                        
                       // Cohesion
 
						ray1=LIST_R[it];
						ray2=LIST_R[jt];   

						ind1=LIST_IND[it];
						ind2=LIST_IND[jt];				                             
						indm=max(ind1,ind2);    

               		     if((LIST_P[it]==0)&&(LIST_P[jt]==0)){	            
							 							 					//	 cout<<"number1"<<endl;                               
						 ray=rmu1*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;

					     indtot+=indm;
						 Smoy+=S;
                         
						 kkn=Emu1*S/nd2;
						 kkt=12.*Emu1*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu1*I/((1+nuu1)*nd2) ;}
	               		 else if(((LIST_P[it]==0)&&(LIST_P[jt]==1))||((LIST_P[it]==1)&&(LIST_P[jt]==0))){	 
							 							 					//	 cout<<"number12"<<endl;                                           
						 
						 
						 //délaminage 25/02/2020
						     
					nint++;	 
						 
						 
				   	     rmui=(rmu1+rmu2)/2.;
                 		 Emui=(Emu1+Emu2)/2.; 
                 		 nuui=(nuu1+nuu2)/2.; 
                 		 					 
						 ray=rmui*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                        

			            		 indtot+=indm;
						 Smoy+=S;

						 kkn=Emui*S/nd2;
						 kkt=12.*Emui*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;  
						 VALCOH[ii][2]=0. ;  	
						 VALCOH[ii][3]=0. ;  	
						 VALCOH[ii][4]=0. ;
						 VALCOH[ii][5]=0. ;  
						 
						 
						 }				 
               		     else if((LIST_P[it]==1)&&(LIST_P[jt]==1)){	         
							 						// cout<<"number2"<<endl;
						                                                  
						 ray=rmu2*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                         

						indtot+=indm;
						Smoy+=S;

						 kkn=Emu2*S/nd2;
						 kkt=12.*Emu2*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu2*I/((1+nuu2)*nd2) ;}
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt]))); 
						 
						 
						  
					 }
						 
                        
}
/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////	

dt=Cs*dt;
Smoy/=indtot;
nint0=nint;
cout<<"Surface moyenne:"<<Smoy<<endl;
cout<<"Pas de temps:"<<dt<<endl;
}















void inicoh2_int_DCB(int & nint, int & nint0,int & npre, R * vs, R * LIST_IND, R H_TOT,int * LIST_B,bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO,  unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, bool * LIST_P, R rmu1, R rmu2, R Emu1, R Emu2, R nuu1, R nuu2, R amort)
{

R dx,dy,dz,nd2,r2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,jt;
R Emui,nuui,rmui;  
R ind1,ind2,indm;

nint=0,npre=0;

R Xa=0.0135;    

cout<<"Xa :"<<Xa<<endl; 
cout<<"% Xa:"<<Xa/H_TOT<<endl;
     
dt=1e9;
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
R Smoy=0.;
R indtot=0.;

# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii,dx,dy,dz,nd2,ind1,ind2,indm,r2,S,I,ray,ray1,ray2,kkn,kkt,fors,norms,it,jt,Emui,nuui,rmui) reduction(min:dt) reduction(+:npre,nint,Smoy,indtot)                         
	for(ii=0;ii<nbco;ii++){
		
		vs[ii]=0.;
		
         //     cout<<"hello"<<endl;
					   it=CONT[ii][0];
					   jt=CONT[ii][1];
					   
					   dx=LIST_X[it]-LIST_X[jt]; 
					   dy=LIST_Y[it]-LIST_Y[jt]; 
					   dz=LIST_Z[it]-LIST_Z[jt]; 					   
					   nd2=sqrt(dx*dx+dy*dy+dz*dz);	   
						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                          
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                     
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];      
						
						
					   if(TYPCO[ii]==1)
		               {
			        
                        
                       // Cohesion
 
						ray1=LIST_R[it];
						ray2=LIST_R[jt];   

						ind1=LIST_IND[it];
						ind2=LIST_IND[jt];				                             
						indm=max(ind1,ind2);    

               		     if((LIST_P[it]==0)&&(LIST_P[jt]==0)){	            
							 							 					//	 cout<<"number1"<<endl;                               
						 ray=rmu1*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;

					     indtot+=indm;
						 Smoy+=S;
                         
						 kkn=Emu1*S/nd2;
						 kkt=12.*Emu1*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu1*I/((1+nuu1)*nd2) ;}
	               		 else if(((LIST_P[it]==0)&&(LIST_P[jt]==1))||((LIST_P[it]==1)&&(LIST_P[jt]==0))){	 
							 							 					//	 cout<<"number12"<<endl;                                           
						 
						 
						 //délaminage 25/02/2020
						     
					nint++;	  
						 if ((LIST_X[it]<Xa)&&(LIST_X[jt]<Xa)){
							 
				            kkn=0.;
				            	
				            LIST_B[it]=2;LIST_B[jt]=2;TYPCO[ii]=0;
				            npre=npre+1;
				            }
						 
						 
				   	     rmui=(rmu1+rmu2)/2.;
                 		 Emui=(Emu1+Emu2)/2.; 
                 		 nuui=(nuu1+nuu2)/2.; 
                 		 					 
						 ray=rmui*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                        

			             indtot+=indm;
						 Smoy+=S;

						 //kkn=Emui*S/nd2;
						 kkn=Emui*S/nd2;
						 kkt=0.;//12.*Emui*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=0.;  
						 VALCOH[ii][2]=0. ;  	
						 VALCOH[ii][3]=0. ;  	
						 VALCOH[ii][4]=0. ;
						 VALCOH[ii][5]=0. ;  
						 
						 
						 }				 
               		     else if((LIST_P[it]==1)&&(LIST_P[jt]==1)){	         
							 						// cout<<"number2"<<endl;
						                                                  
						 ray=rmu2*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                         

						indtot+=indm;
						Smoy+=S;

						 kkn=Emu2*S/nd2;
						 kkt=12.*Emu2*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu2*I/((1+nuu2)*nd2) ;}
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt]))); 
						 
						 
						  
					 }
						 
                        
}
/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////	

dt=Cs*dt;
Smoy/=indtot;
nint0=nint;
nint=nint0-npre;
cout<<"% bonds:"<<npre/R(nint0)<<endl;
cout<<"Surface moyenne:"<<Smoy<<endl;
cout<<"Pas de temps:"<<dt<<endl;
}

void inicoh2_int_ELS(int & nint, int & nint0,int & npre, R * vs, R * LIST_IND, R H_TOT,int * LIST_B,bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO,  unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, bool * LIST_P, R rmu1, R rmu2, R Emu1, R Emu2, R nuu1, R nuu2, R amort)
{

R dx,dy,dz,nd2,r2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,jt;
R Emui,nuui,rmui;  
R ind1,ind2,indm;

nint=0,npre=0;

R Xa=0.06;    

cout<<"Xa :"<<Xa<<endl; 
cout<<"% Xa:"<<Xa/H_TOT<<endl;
     
dt=1e9;
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
R Smoy=0.;
R indtot=0.;

# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii,dx,dy,dz,nd2,ind1,ind2,indm,r2,S,I,ray,ray1,ray2,kkn,kkt,fors,norms,it,jt,Emui,nuui,rmui) reduction(min:dt) reduction(+:npre,nint,Smoy,indtot)                         
	for(ii=0;ii<nbco;ii++){
		
		vs[ii]=0.;
		
         //     cout<<"hello"<<endl;
					   it=CONT[ii][0];
					   jt=CONT[ii][1];
					   
					   dx=LIST_X[it]-LIST_X[jt]; 
					   dy=LIST_Y[it]-LIST_Y[jt]; 
					   dz=LIST_Z[it]-LIST_Z[jt]; 					   
					   nd2=sqrt(dx*dx+dy*dy+dz*dz);	   
						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                          
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                     
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];      
						
						
					   if(TYPCO[ii]==1)
		               {
			        
                        
                       // Cohesion
 
						ray1=LIST_R[it];
						ray2=LIST_R[jt];   

						ind1=LIST_IND[it];
						ind2=LIST_IND[jt];				                             
						indm=max(ind1,ind2);    

               		     if((LIST_P[it]==0)&&(LIST_P[jt]==0)){	            
							 							 					//	 cout<<"number1"<<endl;                               
						 ray=rmu1*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;

					     indtot+=indm;
						 Smoy+=S;
                         
						 kkn=Emu1*S/nd2;
						 kkt=12.*Emu1*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu1*I/((1+nuu1)*nd2) ;}
	               		 else if(((LIST_P[it]==0)&&(LIST_P[jt]==1))||((LIST_P[it]==1)&&(LIST_P[jt]==0))){	 
							 							 					//	 cout<<"number12"<<endl;                                           
						 
						 
						 //délaminage 25/02/2020
						     
					nint++;	  
						 if ((LIST_X[it]<Xa)&&(LIST_X[jt]<Xa)){
							 
				            kkn=0.;
				            	
				            LIST_B[it]=2;LIST_B[jt]=2;TYPCO[ii]=0;
				            npre=npre+1;
				            }
						 
						 
				   	     rmui=(rmu1+rmu2)/2.;
                 		 Emui=(Emu1+Emu2)/2.; 
                 		 nuui=(nuu1+nuu2)/2.; 
                 		 					 
						 ray=rmui*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                        

			             indtot+=indm;
						 Smoy+=S;

						 //kkn=Emui*S/nd2;
						 kkn=0.;
						 kkt=12.*Emui*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=0.;
						 VALCOH[ii][1]=kkt;  
						 VALCOH[ii][2]=0. ;  	
						 VALCOH[ii][3]=0. ;  	
						 VALCOH[ii][4]=0. ;
						 VALCOH[ii][5]=0. ;  
						 
						 
						 }				 
               		     else if((LIST_P[it]==1)&&(LIST_P[jt]==1)){	         
							 						// cout<<"number2"<<endl;
						                                                  
						 ray=rmu2*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                         

						indtot+=indm;
						Smoy+=S;

						 kkn=Emu2*S/nd2;
						 kkt=12.*Emu2*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu2*I/((1+nuu2)*nd2) ;}
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt]))); 
						 
						  
					 }
						 
                        
}
/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////	

dt=Cs*dt;
Smoy/=indtot;
nint0=nint;
nint=nint0-npre;
cout<<"% bonds:"<<npre/R(nint0)<<endl;
cout<<"Surface moyenne:"<<Smoy<<endl;
cout<<"Pas de temps:"<<dt<<endl;
}


void inicoh2_int_MMB(int & nint, int & nint0,int & npre, R * vs, R * LIST_IND, R H_TOT,int * LIST_B,bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO,  unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, bool * LIST_P, R rmu1, R rmu2, R Emu1, R Emu2, R nuu1, R nuu2, R amort)
{

R dx,dy,dz,nd2,r2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,jt;
R Emui,nuui,rmui;  
R ind1,ind2,indm;

nint=0,npre=0;

R Xa=0.03;    

cout<<"Xa :"<<Xa<<endl; 
cout<<"% Xa:"<<Xa/H_TOT<<endl;
     
dt=1e9;
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
R Smoy=0.;
R indtot=0.;

# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii,dx,dy,dz,nd2,ind1,ind2,indm,r2,S,I,ray,ray1,ray2,kkn,kkt,fors,norms,it,jt,Emui,nuui,rmui) reduction(min:dt) reduction(+:npre,nint,Smoy,indtot)                         
	for(ii=0;ii<nbco;ii++){
		
		vs[ii]=0.;
		
         //     cout<<"hello"<<endl;
					   it=CONT[ii][0];
					   jt=CONT[ii][1];
					   
					   dx=LIST_X[it]-LIST_X[jt]; 
					   dy=LIST_Y[it]-LIST_Y[jt]; 
					   dz=LIST_Z[it]-LIST_Z[jt]; 					   
					   nd2=sqrt(dx*dx+dy*dy+dz*dz);	   
						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                          
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                     
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];      
						
						
					   if(TYPCO[ii]==1)
		               {
			        
                        
                       // Cohesion
 
						ray1=LIST_R[it];
						ray2=LIST_R[jt];   

						ind1=LIST_IND[it];
						ind2=LIST_IND[jt];				                             
						indm=max(ind1,ind2);    

               		     if((LIST_P[it]==0)&&(LIST_P[jt]==0)){	            
							 							 					//	 cout<<"number1"<<endl;                               
						 ray=rmu1*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;

					     indtot+=indm;
						 Smoy+=S;
                         
						 kkn=Emu1*S/nd2;
						 kkt=12.*Emu1*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu1*I/((1+nuu1)*nd2) ;}
	               		 else if(((LIST_P[it]==0)&&(LIST_P[jt]==1))||((LIST_P[it]==1)&&(LIST_P[jt]==0))){	 
							 							 					//	 cout<<"number12"<<endl;                                           
						 
						 
						 //délaminage 25/02/2020
						     
					nint++;	  
					 
				   	     rmui=(rmu1+rmu2)/2.;
                 		 Emui=(Emu1+Emu2)/2.; 
                 		 nuui=(nuu1+nuu2)/2.; 
                 		 					 
						 ray=rmui*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                        

			             indtot+=indm;
						 Smoy+=S;

						 kkn=Emui*S/nd2;
						 kkt=12.*Emui*I/(pow(nd2,3));   

						 if ((LIST_X[it]<Xa)&&(LIST_X[jt]<Xa)){
							 
				            kkn=0.;
				            kkt=0.;
				            	
				            LIST_B[it]=2;LIST_B[jt]=2;TYPCO[ii]=0;
				            npre=npre+1;
				          }
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;  
						 VALCOH[ii][2]=0. ;  	
						 VALCOH[ii][3]=0. ;  	
						 VALCOH[ii][4]=0. ;
						 VALCOH[ii][5]=0. ;  
					 
						 
						 }				 
               		     else if((LIST_P[it]==1)&&(LIST_P[jt]==1)){	         
							 						// cout<<"number2"<<endl;
						                                                  
						 ray=rmu2*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                         

						indtot+=indm;
						Smoy+=S;

						 kkn=Emu2*S/nd2;
						 kkt=12.*Emu2*I/(pow(nd2,3));   
					     
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu2*I/((1+nuu2)*nd2) ;}
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt]))); 
						 
						 
						  
					 }
						 
                        
}
/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////	

dt=Cs*dt;
Smoy/=indtot;
nint0=nint;
nint=nint0-npre;
cout<<"% bonds:"<<npre/R(nint0)<<endl;
cout<<"Surface moyenne:"<<Smoy<<endl;
cout<<"Pas de temps:"<<dt<<endl;
}







void inicoh(bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, R Emu, R rmu, R nuu, R amort)
{

R dx,dy,dz,nd2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,jt;
          
R mkn,mkt;
int nbcoc=0;          
mkn=0.;
mkt=0.;     
dt=1e9;
     
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(it,jt,fors,norms,kkn,kkt,dx,dy,dz,nd2,S,I,ray,ray1,ray2) reduction(+:mkn,mkt,nbcoc) reduction(min:dt)    
	for(ii=0;ii<nbco;ii++){
         
			it=CONT[ii][0];
			jt=CONT[ii][1];

			dx=LIST_X[it]-LIST_X[jt]; 
			dy=LIST_Y[it]-LIST_Y[jt]; 
			dz=LIST_Z[it]-LIST_Z[jt]; 					   
			nd2=sqrt(dx*dx+dy*dy+dz*dz);	   

						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                         
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                    
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }
                       
                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];              
                        
                                               
                       // Cohesion
 
 					   if(TYPCO[ii]==1)
		               {
						   

						ray1=LIST_R[it];
						ray2=LIST_R[jt];   
				                                                           
						 ray=rmu*(ray1+ray2)/2.;
						 S=(Pi*ray*ray);
						 I=Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emu*S/nd2;
						 kkt=12.*Emu*I/(pow(nd2,3)); 
						 
						  nbcoc++;
						  
						 mkn=mkn+kkn;
						 mkt=mkt+kkt;    
						                   
						 //                         cout<<ray<<", "<<S<<", "<<I<<", "<<kkn<<", "<<kkt<<endl;            
         //    if(ii==0) cout<<"kkn:"<<kkn<<endl;  
                               
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu*I/((1+nuu)*nd2) ;
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt])));  
                        }
}

mkt=mkt/nbcoc;
mkn=mkn/nbcoc;
cout<<"mkkt:"<<mkt<<", mkkn:"<<mkn<<endl;  

/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////	
		
dt=Cs*dt;
cout<<"Pas de temps:"<<dt<<endl;
}


void inicoh_coef(R * LIST_IND, bool * TYPCO, R & dt, R epsi, int  nbco,  int ** CONT, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, R Emu, R rmu, R nuu, R amort)
{

R dx,dy,dz,nd2,S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
R ind1,ind2,indm;

int it,jt;
          
R mkn,mkt;
int nbcoc=0;          
mkn=0.;
mkt=0.;     
dt=1e9;
     
int ii,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}
      
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(it,jt,ind1,ind2,indm,fors,norms,kkn,kkt,dx,dy,dz,nd2,S,I,ray,ray1,ray2) reduction(+:mkn,mkt,nbcoc) reduction(min:dt)    
	for(ii=0;ii<nbco;ii++){
         
			it=CONT[ii][0];
			jt=CONT[ii][1];

			dx=LIST_X[it]-LIST_X[jt]; 
			dy=LIST_Y[it]-LIST_Y[jt]; 
			dz=LIST_Z[it]-LIST_Z[jt]; 					   
			nd2=sqrt(dx*dx+dy*dy+dz*dz);	   

						   
                       DCONTX[ii]=dx;
                       DCONTY[ii]=dy;
                       DCONTZ[ii]=dz;  
                       DCONTXO[ii]=DCONTX[ii];
                       DCONTYO[ii]=DCONTY[ii];
                       DCONTZO[ii]=DCONTZ[ii];                         
                       DCONT[ii]=nd2;
                       DCONTO[ii]=nd2;                    
                       
                       // Normale
                       
                       NCONT[ii][0]=dx/nd2;
                       NCONT[ii][1]=dy/nd2; 
                       NCONT[ii][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[ii][0])<1e-18) NCONT[ii][0]=0.;
                       if(fabs(NCONT[ii][1])<1e-18) NCONT[ii][1]=0.;
                       if(fabs(NCONT[ii][2])<1e-18) NCONT[ii][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[ii][0]*NCONT[ii][1]==0.)&&(NCONT[ii][1]*NCONT[ii][2]==0.)&&(NCONT[ii][2]*NCONT[ii][0]==0.)){
						  
						  if(NCONT[ii][0]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 1.;						
					      }else if(NCONT[ii][1]!=0.){
						   NCONT[ii][3] = 1.;
						   NCONT[ii][4] = 0.;
						   NCONT[ii][5] = 0.;									  
						  }else if(NCONT[ii][2]!=0.){
						   NCONT[ii][3] = 0.;
						   NCONT[ii][4] = 1.;
						   NCONT[ii][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[ii][0]==0.){
							
							fors  = -(NCONT[ii][1])/NCONT[ii][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[ii][3] = 1./norms;
						    NCONT[ii][4] = 1./norms;
						    NCONT[ii][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[ii][1]+NCONT[ii][2])/NCONT[ii][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[ii][3] = fors/norms;
						   NCONT[ii][4] = 1./norms;
						   NCONT[ii][5] = 1./norms;							
							
						    }
													
							
					   }
                       
                       
					   			   
					   // Tangente t				   
					   
						NCONT[ii][6] = NCONT[ii][1]*NCONT[ii][5] - NCONT[ii][2]*NCONT[ii][4];
						NCONT[ii][7] = NCONT[ii][2]*NCONT[ii][3] - NCONT[ii][0]*NCONT[ii][5];
						NCONT[ii][8] = NCONT[ii][0]*NCONT[ii][4] - NCONT[ii][1]*NCONT[ii][3];              
                        
                                               
                       // Cohesion
 
 					   if(TYPCO[ii]==1)
		               {
						   

						ray1=LIST_R[it];
						ray2=LIST_R[jt];   

						ind1=LIST_IND[it];
						ind2=LIST_IND[jt];				                             
						indm=max(ind1,ind2);                              
						
						 ray=rmu*(ray1+ray2)/2.;
						 S=indm*(Pi*ray*ray);
						 I=indm*Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emu*S/nd2;
						 kkt=12.*Emu*I/(pow(nd2,3)); 
						 
						  nbcoc++;
						  
						 mkn=mkn+kkn;
						 mkt=mkt+kkt;    
						                   
						 //                         cout<<ray<<", "<<S<<", "<<I<<", "<<kkn<<", "<<kkt<<endl;            
         //    if(ii==0) cout<<"kkn:"<<kkn<<endl;  
                               
 						 VALCOH[ii][0]=kkn;
						 VALCOH[ii][1]=kkt;
						 VALCOH[ii][2]=nd2*kkt/2.;	
						 VALCOH[ii][3]=nd2*nd2*kkt/3.;	
						 VALCOH[ii][4]=nd2*nd2*kkt/6.;
						 VALCOH[ii][5]=Emu*I/((1+nuu)*nd2) ;
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt])));  
                        }
}

mkt=mkt/nbcoc;
mkn=mkn/nbcoc;
cout<<"mkkt:"<<mkt<<", mkkn:"<<mkn<<endl;  

/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)) private(ii) 
	for(ii=0;ii<nbco;ii++){
  						 VALAMO[ii][0]=amort*dt*VALCOH[ii][0];
						 VALAMO[ii][1]=amort*dt*VALCOH[ii][1];
						 VALAMO[ii][2]=amort*dt*VALCOH[ii][2];	
						 VALAMO[ii][3]=amort*dt*VALCOH[ii][3];	
						 VALAMO[ii][4]=amort*dt*VALCOH[ii][4];                       
						 VALAMO[ii][5]=amort*dt*VALCOH[ii][5];  
        }
if(amort>1) dt=dt/amort;
}
/////////////////////////////////////////////////////////////////////////	
		
dt=Cs*dt;
cout<<"Pas de temps:"<<dt<<endl;
}
