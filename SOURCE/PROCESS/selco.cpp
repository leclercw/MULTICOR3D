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
# include "omp.h"

using namespace std;

#include "selco.h"


void selco(R & dt, R epsi, int & nbco, int NB_SPH, int ** CONT, int vecsize, int vecsizex, int vecsizey, int vecsizez, int ** coul, int * nocoul, int ** numc, int * nonumc, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT, R * DCONTO, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R Pi, R E, R nu, int & nbcop,int ** CONTP, R * DCONTP, R ** NCONTP, int NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ, R ** LIST_PN, int NMAXCONT, int NMAXCONTP, int NMAXZ)
{

R dx,dy,dz,nd2,r2;
R fors,norms;
unsigned int nbco1;
int ii,nbcoi,nbcoj;
int nbcop1,num_threads,numth;
R p,dist;
nbcop=0;
nbco=0;
int numax=0;

int kjt;
long int it,jt;

int kt,lt;


  for(it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  } 

// Recherche nouveaux contacts sphere-sphere


# pragma omp parallel 
{
num_threads=omp_get_num_threads();
}


if(num_threads==1){


	for(ii=0;ii<vecsize;ii++){
		
			for(lt=0;lt<nocoul[ii];lt++){
			it = coul[ii][lt];

				for(kt=0;kt<nonumc[ii];kt++){

					 for(kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
						 
						  jt=coul[numc[ii][kt]][kjt];

						  if((kt>0)||((kt==0)&&(jt>it))){
						  dx=LIST_X[it]-LIST_X[jt]; 
						  dy=LIST_Y[it]-LIST_Y[jt]; 
						  dz=LIST_Z[it]-LIST_Z[jt]; 					  
						  nd2=sqrt(dx*dx+dy*dy+dz*dz);
						  r2=(LIST_R[it]+LIST_R[jt]);
					

	/////////////////////////////////////////////////////////////////////////	//Contact
							  if(nd2<r2*(1.+epsi))
							  {  


				                                nbco1=nbco;
				                                nbco++;
								nbcoi=NBCONTCO[it];
								NBCONTCO[it]++;
								nbcoj=NBCONTCO[jt];
								NBCONTCO[jt]++;  
								if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
								if(NBCONTCO[it]>numax) numax=NBCONTCO[it];

								CONT[nbco1][0]=it;
								CONT[nbco1][1]=jt;

								NOCONT[it][nbcoi]=nbco1;
								NOCONT[jt][nbcoj]=nbco1;                 
									   
									   DCONTX[nbco1]=dx;
									   DCONTY[nbco1]=dy;
									   DCONTZ[nbco1]=dz;	   
									   DCONT[nbco1] =nd2;
									   DCONTO[nbco1]=r2;
								   
						         // Normale
									   
									   NCONT[nbco1][0]=dx/nd2;
									   NCONT[nbco1][1]=dy/nd2; 
									   NCONT[nbco1][2]=dz/nd2;       
									   
									   // Zero numerique
									   
									   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
									   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
									   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
									   
									   // Tangente s
									   
									   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){

										   NCONT[nbco1][3] = NCONT[nbco1][1];
										   NCONT[nbco1][4] = NCONT[nbco1][2];
										   NCONT[nbco1][5] = NCONT[nbco1][0];
										   
									   }else{

									
											if(NCONT[nbco1][0]==0.){
									
											fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
											norms = sqrt(fors*fors+2.);

											NCONT[nbco1][3] = 1./norms;
											NCONT[nbco1][4] = 1./norms;
											NCONT[nbco1][5] = fors/norms;					
									
											}else{
										
											fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
											norms = sqrt(fors*fors+2.);

											NCONT[nbco1][3] = fors/norms;
											NCONT[nbco1][4] = 1./norms;
											NCONT[nbco1][5] = 1./norms;							
									
											}
															
									
									   }

									   
									   // Tangente t				   
									   
										NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
										NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
										NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];  




		                                                

							  }
	/////////////////////////////////////////////////////////////////////////	

						}




		                   }
			       }







			}


	}


}else if(num_threads>1){

int pnbco[num_threads-1];

int PNBCONTCO[num_threads-1][NB_SPH];
int PCONT[num_threads-1][8*NB_SPH][2];
unsigned int PNOCONT[num_threads-1][NB_SPH][NMAXZ];
R PDCONTX[num_threads-1][8*NB_SPH];
R PDCONTY[num_threads-1][8*NB_SPH];
R PDCONTZ[num_threads-1][8*NB_SPH];
R PDCONT[num_threads-1][8*NB_SPH];
R PDCONTO[num_threads-1][8*NB_SPH];
R PNCONT[num_threads-1][8*NB_SPH][9];

for(it=0;it<num_threads-1;it++){

pnbco[it]=0;

	for(jt=0;jt<NB_SPH;jt++){
	PNBCONTCO[it][jt]=0;
	}

}


nbco=0;
# pragma omp parallel for schedule(static) private(ii,it,jt,kt,kjt,lt,dx,dy,dz,nd2,r2,numth,nbco1,nbcoi,nbcoj,fors,norms)
	for(ii=0;ii<vecsize;ii++){

	numth=omp_get_thread_num()-1;

		for(lt=0;lt<nocoul[ii];lt++){
			it = coul[ii][lt];

				for(kt=0;kt<nonumc[ii];kt++){

					 for(kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
						 
						  jt=coul[numc[ii][kt]][kjt];

 						 if((kt>0)||((kt==0)&&(jt>it))){

						  dx=LIST_X[it]-LIST_X[jt]; 
						  dy=LIST_Y[it]-LIST_Y[jt]; 
						  dz=LIST_Z[it]-LIST_Z[jt]; 					  
						  nd2=sqrt(dx*dx+dy*dy+dz*dz);
						  r2=(LIST_R[it]+LIST_R[jt]);
					


	/////////////////////////////////////////////////////////////////////////	//Contact
							  if(nd2<r2*(1.+epsi))
							  {  




if(numth==-1){


				                                nbco1=nbco;
				                                nbco++;

								nbcoi=NBCONTCO[it];
								NBCONTCO[it]++;
								nbcoj=NBCONTCO[jt];
								NBCONTCO[jt]++;  

								CONT[nbco1][0]=it;
								CONT[nbco1][1]=jt;

								NOCONT[it][nbcoi]=nbco1;
								NOCONT[jt][nbcoj]=nbco1;                 
									   
									   DCONTX[nbco1]=dx;
									   DCONTY[nbco1]=dy;
									   DCONTZ[nbco1]=dz;	   
									   DCONT[nbco1] =nd2;
									   DCONTO[nbco1]=r2;
								   
						         // Normale
									   
									   NCONT[nbco1][0]=dx/nd2;
									   NCONT[nbco1][1]=dy/nd2; 
									   NCONT[nbco1][2]=dz/nd2;       
									   
									   // Zero numerique
									   
									   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
									   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
									   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
									   
									   // Tangente s
									   
									   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){

										   NCONT[nbco1][3] = NCONT[nbco1][1];
										   NCONT[nbco1][4] = NCONT[nbco1][2];
										   NCONT[nbco1][5] = NCONT[nbco1][0];
										   
									   }else{

									
											if(NCONT[nbco1][0]==0.){
									
											fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
											norms = sqrt(fors*fors+2.);

											NCONT[nbco1][3] = 1./norms;
											NCONT[nbco1][4] = 1./norms;
											NCONT[nbco1][5] = fors/norms;					
									
											}else{
										
											fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
											norms = sqrt(fors*fors+2.);

											NCONT[nbco1][3] = fors/norms;
											NCONT[nbco1][4] = 1./norms;
											NCONT[nbco1][5] = 1./norms;							
									
											}
															
									
									   }

									   
									   // Tangente t				   
									   
										NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
										NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
										NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];  

}
else{



		                                        nbco1=pnbco[numth];
		                                        pnbco[numth]++;
							nbcoi=PNBCONTCO[numth][it];
							PNBCONTCO[numth][it]++;
							nbcoj=PNBCONTCO[numth][jt];
							PNBCONTCO[numth][jt]++;  

							PCONT[numth][nbco1][0]=it;
							PCONT[numth][nbco1][1]=jt;

							PNOCONT[numth][it][nbcoi]=nbco1;
							PNOCONT[numth][jt][nbcoj]=nbco1;                 
								   
								   PDCONTX[numth][nbco1]=dx;
								   PDCONTY[numth][nbco1]=dy;
								   PDCONTZ[numth][nbco1]=dz;	   
								   PDCONT[numth][nbco1] =nd2;
								   PDCONTO[numth][nbco1]=r2;
							   
				                 // Normale
								   
								   PNCONT[numth][nbco1][0]=dx/nd2;
								   PNCONT[numth][nbco1][1]=dy/nd2; 
								   PNCONT[numth][nbco1][2]=dz/nd2;       
								   
								   // Zero numerique
								   
								   if(fabs(PNCONT[numth][nbco1][0])<1e-18) PNCONT[numth][nbco1][0]=0.;
								   if(fabs(PNCONT[numth][nbco1][1])<1e-18) PNCONT[numth][nbco1][1]=0.;
								   if(fabs(PNCONT[numth][nbco1][2])<1e-18) PNCONT[numth][nbco1][2]=0.;
								   
								   // Tangente s
								   
								   if((PNCONT[numth][nbco1][0]*PNCONT[numth][nbco1][1]==0.)&&(PNCONT[numth][nbco1][1]*PNCONT[numth][nbco1][2]==0.)&&(PNCONT[numth][nbco1][2]*PNCONT[numth][nbco1][0]==0.)){
									  
				
									   PNCONT[numth][nbco1][3] = PNCONT[numth][nbco1][1] ;
									   PNCONT[numth][nbco1][4] = PNCONT[numth][nbco1][2];
									   PNCONT[numth][nbco1][5] = PNCONT[numth][nbco1][0] ;
									   
								   }else{
									
										if(PNCONT[numth][nbco1][0]==0.){
									
										fors  = -(PNCONT[numth][nbco1][1])/PNCONT[numth][nbco1][2];
										norms = sqrt(fors*fors+2.);

										PNCONT[numth][nbco1][3] = 1./norms;
										PNCONT[numth][nbco1][4] = 1./norms;
										PNCONT[numth][nbco1][5] = fors/norms;					
									
										}else{
										
									   fors  = -(PNCONT[numth][nbco1][1]+PNCONT[numth][nbco1][2])/PNCONT[numth][nbco1][0];
									   norms = sqrt(fors*fors+2.);

									   PNCONT[numth][nbco1][3] = fors/norms;
									   PNCONT[numth][nbco1][4] = 1./norms;
									   PNCONT[numth][nbco1][5] = 1./norms;							
									
										}
															
									
								   }
								   
								   
								   // Tangente t				   
								   
									PNCONT[numth][nbco1][6] = PNCONT[numth][nbco1][1]*PNCONT[numth][nbco1][5] - PNCONT[numth][nbco1][2]*PNCONT[numth][nbco1][4];
									PNCONT[numth][nbco1][7] = PNCONT[numth][nbco1][2]*PNCONT[numth][nbco1][3] - PNCONT[numth][nbco1][0]*PNCONT[numth][nbco1][5];
									PNCONT[numth][nbco1][8] = PNCONT[numth][nbco1][0]*PNCONT[numth][nbco1][4] - PNCONT[numth][nbco1][1]*PNCONT[numth][nbco1][3];  

                                                   
}

							    }

	/////////////////////////////////////////////////////////////////////////	
 							

						}



		                   }
			       }







			}


	}






for(it=0;it<num_threads-1;it++){ 

	 for(jt=0;jt<pnbco[it];jt++){
	 CONT[nbco+jt][0]=PCONT[it][jt][0];
	 CONT[nbco+jt][1]=PCONT[it][jt][1];

	 NCONT[nbco+jt][0]=PNCONT[it][jt][0];
	 NCONT[nbco+jt][1]=PNCONT[it][jt][1];
	 NCONT[nbco+jt][2]=PNCONT[it][jt][2];
	 NCONT[nbco+jt][3]=PNCONT[it][jt][3];
	 NCONT[nbco+jt][4]=PNCONT[it][jt][4];
	 NCONT[nbco+jt][5]=PNCONT[it][jt][5];
	 NCONT[nbco+jt][6]=PNCONT[it][jt][6];
	 NCONT[nbco+jt][7]=PNCONT[it][jt][7];
	 NCONT[nbco+jt][8]=PNCONT[it][jt][8];

	 DCONT[nbco+jt]=PDCONT[it][jt];
	 DCONTO[nbco+jt]=PDCONTO[it][jt];
	 DCONTX[nbco+jt]=PDCONTX[it][jt];
	 DCONTY[nbco+jt]=PDCONTY[it][jt];
	 DCONTZ[nbco+jt]=PDCONTZ[it][jt];
	 }

	for(jt=0;jt<NB_SPH;jt++){
		for(kjt=0;kjt<PNBCONTCO[it][jt];kjt++){ 
		NOCONT[jt][kjt+NBCONTCO[jt]]=PNOCONT[it][jt][kjt]+nbco;
		}
	NBCONTCO[jt]+=PNBCONTCO[it][jt];
        if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];	
        }

nbco+=pnbco[it];
}









}

// Recherche nouveaux contacts sphere-paroi

	for(jt=0;jt<NB_SPH;jt++){

		for(it=0;it<NB_PAR;it++){

          
		p= (LIST_PN[it][0]*(LIST_X[jt]-LIST_PX[it][0])+LIST_PN[it][1]*(LIST_Y[jt]-LIST_PY[it][0])+LIST_PN[it][2]*(LIST_Z[jt]-LIST_PZ[it][0]));
                dist=fabs(p);
     
		  if(dist<LIST_R[jt]*(1.+epsi)){
				nbcop1=nbcop;
				nbcop++;

				CONTP[nbcop1][0]=it;
	          		CONTP[nbcop1][1]=jt;
				DCONTP[nbcop1]=dist;
				NCONTP[nbcop1][0]=-p*LIST_PN[it][0]/dist;
				NCONTP[nbcop1][1]=-p*LIST_PN[it][1]/dist; 
				NCONTP[nbcop1][2]=-p*LIST_PN[it][2]/dist;   
			  
		  }
		  
		
		}
		
	}

if(numax>NMAXZ){
cout<<"Nombre maximal de contacts par sph atteint: "<<numax<<"/"<<NMAXZ<<" - stop "<<endl;
exit(0);
}

if(nbco>NMAXCONT){
cout<<"Nombre maximal de contacts sph-sph atteint:"<<nbco<<"/"<<NMAXCONT<<" - stop "<<endl;
exit(0);
}

if(nbcop>NMAXCONTP){
cout<<"Nombre maximal de contacts sph-paroi atteint:"<<nbcop<<"/"<<NMAXCONTP<<" - stop "<<endl;
exit(0);
}

}

void selco_paroi(R epsi, int & nbcop,int ** CONTP, R * DCONTP, R ** NCONTP, int NB_PAR, R ** LIST_PX, R ** LIST_PY, R ** LIST_PZ, R ** LIST_PN, int NB_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, int NMAXCONTP) 
{
	
R p,dist;	
nbcop=0;
int it,jt,nbcop1;

	for(jt=0;jt<NB_SPH;jt++){
		
		for(it=0;it<NB_PAR;it++){

          
		p= (LIST_PN[it][0]*(LIST_X[jt]-LIST_PX[it][0])+LIST_PN[it][1]*(LIST_Y[jt]-LIST_PY[it][0])+LIST_PN[it][2]*(LIST_Z[jt]-LIST_PZ[it][0]));
                dist=fabs(p);
     
		  if(dist<LIST_R[jt]*(1.+epsi)){
				nbcop1=nbcop;
				nbcop++;

				CONTP[nbcop1][0]=it;
	          		CONTP[nbcop1][1]=jt;
				DCONTP[nbcop1]=dist;
				NCONTP[nbcop1][0]=-p*LIST_PN[it][0]/dist;
				NCONTP[nbcop1][1]=-p*LIST_PN[it][1]/dist; 
				NCONTP[nbcop1][2]=-p*LIST_PN[it][2]/dist;   
			  
		  }
		  
		
		}
		
	}	


if(nbcop>NMAXCONTP){
cout<<"Nombre maximal de contacts sph-paroi atteint:"<<nbcop<<"/"<<NMAXCONTP<<" - stop "<<endl;
exit(0);
}
	
}


void selco_corci_qs(R & dt, R epsi, unsigned int & nbco,  int ** CONT, int vecsize, int vecsizex, int vecsizey, int vecsizez, R H_TOT, R V_TOT, R Z_TOT, int ** coul, int * nocoul, int ** numc, int * nonumc, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R Pi, R Cs, R Emu, R rmu, R nuu, R amort, int NMAXCONT, int NMAXZ)
{

nbco=0;

R dx,dy,dz,nd2,r2, S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
int it,ii,lt,kt,kjt,jt;
int numax=0;


                      
	for(ii=0;ii<vecsize;ii++){
		
			for(lt=0;lt<nocoul[ii];lt++){
			it = coul[ii][lt];

				for(kt=0;kt<nonumc[ii];kt++){

					 for(kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
					 
					  jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Collision
					  if(nd2<r2*(1.+epsi))
					  {
					   CONT[nbco][0]=it;
                       CONT[nbco][1]=jt;
                       
			NOCONT[it][NBCONTCO[it]]=nbco;
			NOCONT[jt][NBCONTCO[jt]]=nbco;
			NBCONTCO[it]++;
			NBCONTCO[jt]++;   
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];					
	   
                       DCONTX[nbco]=dx;
                       DCONTY[nbco]=dy;
                       DCONTZ[nbco]=dz;  
                       DCONT[nbco]=nd2;
                       
                       // Normale
                       
                       NCONT[nbco][0]=dx/nd2;
                       NCONT[nbco][1]=dy/nd2; 
                       NCONT[nbco][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
                       if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
                       if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
						  
						  if(NCONT[nbco][0]!=0.){
						   NCONT[nbco][3] = 0.;
						   NCONT[nbco][4] = 0.;
						   NCONT[nbco][5] = 1.;						
					      }else if(NCONT[nbco][1]!=0.){
						   NCONT[nbco][3] = 1.;
						   NCONT[nbco][4] = 0.;
						   NCONT[nbco][5] = 0.;									  
						  }else if(NCONT[nbco][2]!=0.){
						   NCONT[nbco][3] = 0.;
						   NCONT[nbco][4] = 1.;
						   NCONT[nbco][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[nbco][0]==0.){
							
							fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[nbco][3] = 1./norms;
						    NCONT[nbco][4] = 1./norms;
						    NCONT[nbco][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[nbco][3] = fors/norms;
						   NCONT[nbco][4] = 1./norms;
						   NCONT[nbco][5] = 1./norms;							
							
						    }
													
							
					   }
                       
                  				   			   
					   // Tangente t				   
					   
						NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
						NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
						NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];              
                       
                       // Cohesion
					 
			                                                           
						 ray=rmu*r2/2.;
						 S=(Pi*ray*ray);
						 I=Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emu*S/nd2;
						 kkt=12.*Emu*I/(pow(nd2,3));                                   
                         
 						 VALCOH[nbco][0]=kkn;
						 VALCOH[nbco][1]=kkt;
						 VALCOH[nbco][2]=nd2*kkt/2.;	
						 VALCOH[nbco][3]=nd2*nd2*kkt/3.;	
						 VALCOH[nbco][4]=nd2*nd2*kkt/6.;
						 VALCOH[nbco][5]=Emu*I/((1+nuu)*nd2) ;
                        
  						 VALAMO[nbco][0]=amort*VALCOH[nbco][0];
						 VALAMO[nbco][1]=amort*VALCOH[nbco][1];
						 VALAMO[nbco][2]=amort*VALCOH[nbco][2];	
						 VALAMO[nbco][3]=amort*VALCOH[nbco][3];	
						 VALAMO[nbco][4]=amort*VALCOH[nbco][4];                       
						 VALAMO[nbco][5]=amort*VALCOH[nbco][5];                        
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt])));  

                                                nbco++;      
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

if(numax>NMAXZ){
cout<<"Nombre maximal de contacts par sph atteint: "<<numax<<"/"<<NMAXZ<<" - stop "<<endl;
exit(0);
}

if(nbco>NMAXCONT){
cout<<"Nombre maximal de contacts sph-sph atteint:"<<nbco<<"/"<<NMAXCONT<<" - stop "<<endl;
exit(0);
}


dt=Cs*dt;
cout<<"Pas de temps:"<<dt<<endl;
}


void selco_corci2_qs(R & dt, R epsi, unsigned int & nbco,  int ** CONT, int vecsize, int vecsizex, int vecsizey, int vecsizez, R H_TOT, R V_TOT, R Z_TOT, int ** coul, int * nocoul, int ** numc, int * nonumc, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R Pi, R Cs, bool * LIST_P, R rmu1, R rmu2, R Emu1, R Emu2, R nuu1, R nuu2, R amort, int NMAXCONT, int NMAXZ)
{

nbco=0;

R dx,dy,dz,nd2,r2, S,I,ray,ray1,ray2;
R kkn,kkt;
R fors,norms;
R Emui,nuui,rmui;    
int numax=0;      
                        
	for(int ii=0;ii<vecsize;ii++){
		
			for(int lt=0;lt<nocoul[ii];lt++){
			int it = coul[ii][lt];

				for(int kt=0;kt<nonumc[ii];kt++){

					 for(int kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
					 
					  int jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Collision
					  if(nd2<r2*(1.+epsi))
					  {
					   CONT[nbco][0]=it;
                       CONT[nbco][1]=jt;
                       
                       NOCONT[it][NBCONTCO[it]]=nbco;
					   NOCONT[jt][NBCONTCO[jt]]=nbco;
					   NBCONTCO[it]++;
					   NBCONTCO[jt]++;   
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];
						   
                       DCONTX[nbco]=dx;
                       DCONTY[nbco]=dy;
                       DCONTZ[nbco]=dz;  
                       DCONT[nbco]=nd2;
                       
                       // Normale
                       
                       NCONT[nbco][0]=dx/nd2;
                       NCONT[nbco][1]=dy/nd2; 
                       NCONT[nbco][2]=dz/nd2;       
                       
                       // Zero numerique
                       
                       if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
                       if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
                       if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
                       
					   // Tangente s
                       
                       if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
						  
						  if(NCONT[nbco][0]!=0.){
						   NCONT[nbco][3] = 0.;
						   NCONT[nbco][4] = 0.;
						   NCONT[nbco][5] = 1.;						
					      }else if(NCONT[nbco][1]!=0.){
						   NCONT[nbco][3] = 1.;
						   NCONT[nbco][4] = 0.;
						   NCONT[nbco][5] = 0.;									  
						  }else if(NCONT[nbco][2]!=0.){
						   NCONT[nbco][3] = 0.;
						   NCONT[nbco][4] = 1.;
						   NCONT[nbco][5] = 0.;									  
						  }
						   
					   }else{
							
							if(NCONT[nbco][0]==0.){
							
							fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
						    norms = sqrt(fors*fors+2.);

						    NCONT[nbco][3] = 1./norms;
						    NCONT[nbco][4] = 1./norms;
						    NCONT[nbco][5] = fors/norms;					
							
					    	}else{
								
						   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
						   norms = sqrt(fors*fors+2.);

						   NCONT[nbco][3] = fors/norms;
						   NCONT[nbco][4] = 1./norms;
						   NCONT[nbco][5] = 1./norms;							
							
						    }
													
							
					   }
                       
				   			   
					   // Tangente t				   
					   
						NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
						NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
						NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];              
                       
                       // Cohesion
					 
				
               		     if((LIST_P[it]==0)&&(LIST_P[jt]==0)){	                                            

						 ray=rmu1*r2/2.;
						 S=(Pi*ray*ray);
						 I=Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emu1*S/nd2;
						 kkt=12.*Emu1*I/(pow(nd2,3));   
					     
 						 VALCOH[nbco][0]=kkn;
						 VALCOH[nbco][1]=kkt;
						 VALCOH[nbco][2]=nd2*kkt/2.;	
						 VALCOH[nbco][3]=nd2*nd2*kkt/3.;	
						 VALCOH[nbco][4]=nd2*nd2*kkt/6.;
						 VALCOH[nbco][5]=Emu1*I/((1+nuu1)*nd2) ;}
	               		 else if(((LIST_P[it]==0)&&(LIST_P[jt]==1))||((LIST_P[it]==1)&&(LIST_P[jt]==0))){	                                            
						 					 
				   	     rmui=(rmu1+rmu2)/2.;
                 		 Emui=(Emu1+Emu2)/2.; 
                 		 nuui=(nuu1+nuu2)/2.; 
                 		 					 
						 ray=rmui*r2/2.;
						 S=(Pi*ray*ray);
						 I=Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emui*S/nd2;
						 kkt=12.*Emui*I/(pow(nd2,3));   
					     
 						 VALCOH[nbco][0]=kkn;
						 VALCOH[nbco][1]=kkt;
						 VALCOH[nbco][2]=nd2*kkt/2.;	
						 VALCOH[nbco][3]=nd2*nd2*kkt/3.;	
						 VALCOH[nbco][4]=nd2*nd2*kkt/6.;
						 VALCOH[nbco][5]=Emui*I/((1+nuui)*nd2) ;}					 
               		     else if((LIST_P[it]==1)&&(LIST_P[jt]==1)){	    
                         
						 ray=rmu2*r2/2.;
						 S=(Pi*ray*ray);
						 I=Pi*pow((2.*ray),4)/64.;
                         
						 kkn=Emu2*S/nd2;
						 kkt=12.*Emu2*I/(pow(nd2,3));   
					     
 						 VALCOH[nbco][0]=kkn;
						 VALCOH[nbco][1]=kkt;
						 VALCOH[nbco][2]=nd2*kkt/2.;	
						 VALCOH[nbco][3]=nd2*nd2*kkt/3.;	
						 VALCOH[nbco][4]=nd2*nd2*kkt/6.;
						 VALCOH[nbco][5]=Emu2*I/((1+nuu2)*nd2) ;}						 
						 
						 
                        
  						 VALAMO[nbco][0]=amort*VALCOH[nbco][0];
						 VALAMO[nbco][1]=amort*VALCOH[nbco][1];
						 VALAMO[nbco][2]=amort*VALCOH[nbco][2];	
						 VALAMO[nbco][3]=amort*VALCOH[nbco][3];	
						 VALAMO[nbco][4]=amort*VALCOH[nbco][4];                       
						 VALAMO[nbco][5]=amort*VALCOH[nbco][5];                        
                       
						 dt=min(dt,sqrt(LIST_M[it]/kkn));
						 dt=min(dt,sqrt(LIST_M[jt]/kkn));
						 dt=min(dt,sqrt(LIST_I[it]/(kkt*LIST_R[it]*LIST_R[it])));
						 dt=min(dt,sqrt(LIST_I[jt]/(kkt*LIST_R[jt]*LIST_R[jt])));  

                       nbco++;      
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

if(numax>NMAXZ){
cout<<"Nombre maximal de contacts par sph atteint: "<<numax<<"/"<<NMAXZ<<" - stop "<<endl;
exit(0);
}

if(nbco>NMAXCONT){
cout<<"Nombre maximal de contacts sph-sph atteint:"<<nbco<<"/"<<NMAXCONT<<" - stop "<<endl;
exit(0);
}


dt=Cs*dt;
cout<<"Pas de temps:"<<dt<<endl;
}


void selco_corci_ori(R epsi, int NB_SPH, int & nbco,  int ** CONT, int vecsize, int vecsizex, int vecsizey, int vecsizez, int ** coul, int * nocoul, int ** numc, int * nonumc, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, unsigned int ** NOCONT, int * NBCONTCO, int NMAXCONT, int NMAXZ)
{

R dx,dy,dz,nd2,r2;

unsigned int nbco1;
int ii,nbcoi,nbcoj;
int num_threads,numth;

nbco=0;
int numax=0;

int kjt;
long int it,jt;

int kt,lt;


  for(it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  } 

// Recherche nouveaux contacts sphere-sphere


# pragma omp parallel 
{
num_threads=omp_get_num_threads();
}


if(num_threads==1){


	for(ii=0;ii<vecsize;ii++){
		
			for(lt=0;lt<nocoul[ii];lt++){
			it = coul[ii][lt];

				for(kt=0;kt<nonumc[ii];kt++){

					 for(kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
						 
						  jt=coul[numc[ii][kt]][kjt];

						  if((kt>0)||((kt==0)&&(jt>it))){
						  dx=LIST_X[it]-LIST_X[jt]; 
						  dy=LIST_Y[it]-LIST_Y[jt]; 
						  dz=LIST_Z[it]-LIST_Z[jt]; 					  
						  nd2=sqrt(dx*dx+dy*dy+dz*dz);
						  r2=(LIST_R[it]+LIST_R[jt]);
					

	/////////////////////////////////////////////////////////////////////////	//Contact
							  if(nd2<r2*(1.+epsi))
							  {  


				                                nbco1=nbco;
				                                nbco++;
								nbcoi=NBCONTCO[it];
								NBCONTCO[it]++;
								nbcoj=NBCONTCO[jt];
								NBCONTCO[jt]++;  
								if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
								if(NBCONTCO[it]>numax) numax=NBCONTCO[it];

								CONT[nbco1][0]=it;
								CONT[nbco1][1]=jt;

								NOCONT[it][nbcoi]=nbco1;
								NOCONT[jt][nbcoj]=nbco1;                 
									   


							  }
	/////////////////////////////////////////////////////////////////////////	

						}




		                   }
			       }







			}


	}


}else{

int pnbco[num_threads-1];

int PNBCONTCO[num_threads-1][NB_SPH];
int PCONT[num_threads-1][8*NB_SPH][2];
unsigned int PNOCONT[num_threads-1][NB_SPH][NMAXZ];


for(it=0;it<num_threads-1;it++){

pnbco[it]=0;

	for(jt=0;jt<NB_SPH;jt++){
	PNBCONTCO[it][jt]=0;
	}

}


nbco=0;
# pragma omp parallel for schedule(static) private(ii,it,jt,kt,kjt,lt,dx,dy,dz,nd2,r2,numth,nbco1,nbcoi,nbcoj)
	for(ii=0;ii<vecsize;ii++){

	numth=omp_get_thread_num()-1;

		for(lt=0;lt<nocoul[ii];lt++){
			it = coul[ii][lt];

				for(kt=0;kt<nonumc[ii];kt++){

					 for(kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
						 
						  jt=coul[numc[ii][kt]][kjt];

 						 if((kt>0)||((kt==0)&&(jt>it))){

						  dx=LIST_X[it]-LIST_X[jt]; 
						  dy=LIST_Y[it]-LIST_Y[jt]; 
						  dz=LIST_Z[it]-LIST_Z[jt]; 					  
						  nd2=sqrt(dx*dx+dy*dy+dz*dz);
						  r2=(LIST_R[it]+LIST_R[jt]);
					


	/////////////////////////////////////////////////////////////////////////	//Contact
							  if(nd2<r2*(1.+epsi))
							  {  

								if(numth==-1){


								nbco1=nbco;
								nbco++;

								nbcoi=NBCONTCO[it];
								NBCONTCO[it]++;
								nbcoj=NBCONTCO[jt];
								NBCONTCO[jt]++;  

								CONT[nbco1][0]=it;
								CONT[nbco1][1]=jt;

								NOCONT[it][nbcoi]=nbco1;
								NOCONT[jt][nbcoj]=nbco1;                 
									   

								}
								else{



								nbco1=pnbco[numth];
								pnbco[numth]++;
								nbcoi=PNBCONTCO[numth][it];
								PNBCONTCO[numth][it]++;
								nbcoj=PNBCONTCO[numth][jt];
								PNBCONTCO[numth][jt]++;  

								PCONT[numth][nbco1][0]=it;
								PCONT[numth][nbco1][1]=jt;

								PNOCONT[numth][it][nbcoi]=nbco1;
								PNOCONT[numth][jt][nbcoj]=nbco1;                 
								   

								}

							    }

	/////////////////////////////////////////////////////////////////////////	
 							

						}



		                   }
			       }







			}


	}


for(it=0;it<num_threads-1;it++){ 

	 for(jt=0;jt<pnbco[it];jt++){
	 CONT[nbco+jt][0]=PCONT[it][jt][0];
	 CONT[nbco+jt][1]=PCONT[it][jt][1];
	 }

	for(jt=0;jt<NB_SPH;jt++){
		for(kjt=0;kjt<PNBCONTCO[it][jt];kjt++){ 
		NOCONT[jt][kjt+NBCONTCO[jt]]=PNOCONT[it][jt][kjt]+nbco;
		}
	NBCONTCO[jt]+=PNBCONTCO[it][jt];
        if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];	
        }

nbco+=pnbco[it];
}



}

cout<<"nbco:"<<nbco<<endl;

if(numax>NMAXZ){
cout<<"Nombre maximal de contacts par sph atteint: "<<numax<<"/"<<NMAXZ<<" - stop "<<endl;
exit(0);
}

if(nbco>NMAXCONT){
cout<<"Nombre maximal de contacts sph-sph atteint:"<<nbco<<"/"<<NMAXCONT<<" - stop "<<endl;
exit(0);
}


}


void selco_corci_ind(bool * TYPCO, R & dt, R epsi, int & nbco, int & nbcoh, int & nbcoe, int & nbci, int NB_SPH, int ** CONT, int vecsize, int vecsizex, int vecsizey, int vecsizez, R H_TOT, R V_TOT, R Z_TOT, int ** coul, int * nocoul, int ** numc, int * nonumc,R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT, R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R Pi, R Cs, R E1, R nu1, R E2, R nu2, int NMAXCONT, int NMAXZ)
{

R dx,dy,dz,nd2,r2,rij,hij,iEij,Eij;
R Kn,Kt;
R fors,norms;

// Mise à blanc des contacts par corci
  for(int it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  }

dt=1e9;

// Traitement anciens contacts
unsigned int nbco1=0;
nbcoh=0;
nbcoe=0;
int numax=0;

  for(int ii=0;ii<nbco;ii++){
	  
	  if(TYPCO[ii]){
	       nbcoh++;
		   TYPCO[nbco1]=1;  	
		   
                  CONT[nbco1][0]=CONT[ii][0];
		   CONT[nbco1][1]=CONT[ii][1];                       
		   
		   int it=CONT[nbco1][0];
		   int jt=CONT[nbco1][1];     
		    
   		   DCONTX[nbco1]=DCONTX[ii];
		   DCONTY[nbco1]=DCONTY[ii];
		   DCONTZ[nbco1]=DCONTZ[ii];		   
		   DCONT[nbco1] =DCONT[ii];		
		   DCONTO[nbco1]=DCONTO[ii];	

		   NCONT[nbco1][0]=NCONT[ii][0];
		   NCONT[nbco1][1]=NCONT[ii][1]; 
		   NCONT[nbco1][2]=NCONT[ii][2];
		   NCONT[nbco1][3]=NCONT[ii][3];
		   NCONT[nbco1][4]=NCONT[ii][4]; 
		   NCONT[nbco1][5]=NCONT[ii][5];     
		   NCONT[nbco1][6]=NCONT[ii][6];
		   NCONT[nbco1][7]=NCONT[ii][7]; 
		   NCONT[nbco1][8]=NCONT[ii][8];     		                     		    
                       		    
		   NOCONT[it][NBCONTCO[it]]=nbco1;
		   NOCONT[jt][NBCONTCO[jt]]=nbco1;
		   NBCONTCO[it]++;
		   NBCONTCO[jt]++; 		                       		    
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];

			 VALCOH[nbco1][0]=VALCOH[ii][0];
			 VALCOH[nbco1][1]=VALCOH[ii][1];
			 VALCOH[nbco1][2]=VALCOH[ii][2];	
			 VALCOH[nbco1][3]=VALCOH[ii][3];	
			 VALCOH[nbco1][4]=VALCOH[ii][4];
			 VALCOH[nbco1][5]=VALCOH[ii][5];
			 VALCOH[nbco1][6]=VALCOH[ii][6];	
			
			 VALAMO[nbco1][0]=VALAMO[ii][0];
			 VALAMO[nbco1][1]=VALAMO[ii][1];
			 VALAMO[nbco1][2]=VALAMO[ii][2];
			 VALAMO[nbco1][3]=VALAMO[ii][3];
			 VALAMO[nbco1][4]=VALAMO[ii][4];                      
			 VALAMO[nbco1][5]=VALAMO[ii][5];
			 VALAMO[nbco1][6]=VALAMO[ii][6];	                       
		   
			 dt=min(dt,sqrt(LIST_M[it]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_M[jt]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_I[it]/(VALCOH[nbco1][1]*LIST_R[it]*LIST_R[it])));
			 dt=min(dt,sqrt(LIST_I[jt]/(VALCOH[nbco1][1]*LIST_R[jt]*LIST_R[jt])));  

		  
		  nbco1++;
	  }
	  else{/*
	
					   int it=CONT[ii][0];
					   int jt=CONT[ii][1];
		          
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 	
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);  
				/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
						
						   TYPCO[nbco1]=0;  		
						 	  nbcoe++; 
						       CONT[nbco1][0]=it;
							   CONT[nbco1][1]=jt;                       
							   
							   DCONTX[nbco1]=dx;
							   DCONTY[nbco1]=dy;
							   DCONTZ[nbco1]=dz;	   
							   DCONT[nbco1] =nd2;
				               DCONTO[nbco1]=DCONTO[ii];	
  
		                         // Normale
							   
							   NCONT[nbco1][0]=dx/nd2;
							   NCONT[nbco1][1]=dy/nd2; 
							   NCONT[nbco1][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
							   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
							   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){
								  
								  if(NCONT[nbco1][0]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 1.;						
								  }else if(NCONT[nbco1][1]!=0.){
								   NCONT[nbco1][3] = 1.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 0.;									  
								  }else if(NCONT[nbco1][2]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 1.;
								   NCONT[nbco1][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco1][0]==0.){
									
									fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = 1./norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco1][3] = fors/norms;
								   NCONT[nbco1][4] = 1./norms;
								   NCONT[nbco1][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
								NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
								NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco1;
								NOCONT[jt][NBCONTCO[jt]]=nbco1;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
								if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
								if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	
														   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco1]-(LIST_R[it]+LIST_R[jt]));
								 iEij=2.*(1.-nu2*nu2)/E2;
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=10.*4.*Eij*sqrt(rij*hij)/3.;						
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  				
						
						
						
						  nbco1++;		  
					  }	  
					  	  
		*/  
	  }
	  
  }	  

nbco=nbco1;

// Recherche nouveaux contacts
                        /*
	for(int ii=0;ii<vecsize;ii++){
		
			for(int lt=0;lt<nocoul[ii];lt++){
			int it = coul[ii][lt];

				for(int kt=0;kt<nonumc[ii];kt++){

					 for(int kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
					 
					  int jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
                           
							   bool booloc=0 ; 
								  
							   if(NBCONTCO[it]>0){
								 int ite=0;
								 while((ite<NBCONTCO[it])&&(!booloc)){
								   unsigned int numco=NOCONT[it][ite];
								 
								   if((numco<nbco1)&&(((CONT[numco][0]==it)&&(CONT[numco][1]==jt))||((CONT[numco][0]==jt)&&(CONT[numco][1]==it)))) booloc=1;			   	
								
								  ite++;						   
								 } 
							   }	
									
							   if(!booloc){  
                             //      cout<<"hello!"<<endl;
							   TYPCO[nbco]=0;
										  nbcoe++;
							   CONT[nbco][0]=it;
							   CONT[nbco][1]=jt;                       
							   
							   DCONTX[nbco]=dx;
							   DCONTY[nbco]=dy;
							   DCONTZ[nbco]=dz;	   
							   DCONT[nbco] =nd2;
							   DCONTO[nbco1]=r2;
							   
		                         // Normale
							   
							   NCONT[nbco][0]=dx/nd2;
							   NCONT[nbco][1]=dy/nd2; 
							   NCONT[nbco][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
							   if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
							   if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
								  
								  if(NCONT[nbco][0]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 1.;						
								  }else if(NCONT[nbco][1]!=0.){
								   NCONT[nbco][3] = 1.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 0.;									  
								  }else if(NCONT[nbco][2]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 1.;
								   NCONT[nbco][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco][0]==0.){
									
									fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco][3] = 1./norms;
									NCONT[nbco][4] = 1./norms;
									NCONT[nbco][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco][3] = fors/norms;
								   NCONT[nbco][4] = 1./norms;
								   NCONT[nbco][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
								NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
								NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco;
								NOCONT[jt][NBCONTCO[jt]]=nbco;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 	
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco]-(LIST_R[it]+LIST_R[jt]));
								 iEij=2.*(1.-nu2*nu2)/E2;
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=10.*4.*Eij*sqrt(rij*hij)/3.;	 			
								 				                                        Kn=Kn/10000.;		
		
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  		
															   
							   nbco++;						   
							   }														   
						   
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

*/

// Contacts indenteur
								
int jt=NB_SPH-1;

nbci=0.;

int it,nbcoi,nbcoj;

		for(it=0;it<NB_SPH-1;it++){

			
			dx=LIST_X[it]-LIST_X[jt]; 
			dy=LIST_Y[it]-LIST_Y[jt]; 
			dz=LIST_Z[it]-LIST_Z[jt]; 					  
			nd2=sqrt(dx*dx+dy*dy+dz*dz);
			r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {

		                                        nbco1=nbco;
		                                        nbco++;
								nbcoi=NBCONTCO[it];
								NBCONTCO[it]++;
								nbcoj=NBCONTCO[jt];
								NBCONTCO[jt]++;  
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];

                                                           // cout<<"contact détecté"<<endl;
							   TYPCO[nbco1]=0;
									
							   CONT[nbco1][0]=it;
							   CONT[nbco1][1]=jt;                       
							   
							   DCONTX[nbco1]=dx;
							   DCONTY[nbco1]=dy;
							   DCONTZ[nbco1]=dz;	   
							   DCONT[nbco1] =nd2;
							   DCONTO[nbco1] =r2;

		                                           // Normale							   
							   NCONT[nbco1][0]=dx/nd2;
							   NCONT[nbco1][1]=dy/nd2; 
							   NCONT[nbco1][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
							   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
							   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){
								  
								  if(NCONT[nbco1][0]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 0.;
					            		   NCONT[nbco1][5] = 1.;						
								  }else if(NCONT[nbco1][1]!=0.){
								   NCONT[nbco1][3] = 1.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 0.;									  
								  }else if(NCONT[nbco1][2]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 1.;
								   NCONT[nbco1][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco1][0]==0.){
									
									fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = 1./norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = fors/norms;					
									
									}else{
										
									fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = fors/norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
								NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
								NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];    
								
 															   
								NOCONT[it][nbcoi]=nbco1;
								NOCONT[jt][nbcoj]=nbco1;
									   
								rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								hij=fabs(DCONT[nbco1]-(LIST_R[it]+LIST_R[jt]));
								iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;								
								Eij=1./iEij;	
								Kn=10.*4.*Eij*sqrt(rij*hij)/3.;							
								Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  		

								nbci++;							   
													   
						   
					  }

 }

if(numax>NMAXZ){
cout<<"Nombre maximal de contacts par sph atteint: "<<numax<<"/"<<NMAXZ<<" - stop "<<endl;
exit(0);
}

if(nbco>NMAXCONT){
cout<<"Nombre maximal de contacts sph-sph atteint:"<<nbco<<"/"<<NMAXCONT<<" - stop "<<endl;
exit(0);
}


dt=Cs*dt;
//cout<<"Pas de temps:"<<dt<<endl;
}

void selco_corci(bool * TYPCO, R & dt, R epsi, unsigned int & nbco,int & nbcoh,int & nbcoe, int NB_SPH, int ** CONT, int vecsize, int vecsizex, int vecsizey, int vecsizez, R H_TOT, R V_TOT, R Z_TOT, int ** coul, int * nocoul, int ** numc, int * nonumc,R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT, R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R Pi, R Cs, R E, R nu, int NMAXCONT, int NMAXZ)
{

R dx,dy,dz,nd2,r2,rij,hij,iEij,Eij;
R Kn,Kt;
R fors,norms;

// Mise à blanc des contacts par corci
  for(int it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  }

dt=1e9;

// Traitement anciens contacts
unsigned int nbco1=0;
nbcoh=0;
nbcoe=0;
int numax=0;

  for(int ii=0;ii<nbco;ii++){
	  
	  if(TYPCO[ii]){
	       nbcoh++;
		   TYPCO[nbco1]=1;  	
		   
           CONT[nbco1][0]=CONT[ii][0];
		   CONT[nbco1][1]=CONT[ii][1];                       
		   
		   int it=CONT[nbco1][0];
		   int jt=CONT[nbco1][1];     
		    
   		   DCONTX[nbco1]=DCONTX[ii];
		   DCONTY[nbco1]=DCONTY[ii];
		   DCONTZ[nbco1]=DCONTZ[ii];		   
		   DCONT[nbco1] =DCONT[ii];		
		   DCONTO[nbco1]=DCONTO[ii];	

		   NCONT[nbco1][0]=NCONT[ii][0];
		   NCONT[nbco1][1]=NCONT[ii][1]; 
		   NCONT[nbco1][2]=NCONT[ii][2];
		   NCONT[nbco1][3]=NCONT[ii][3];
		   NCONT[nbco1][4]=NCONT[ii][4]; 
		   NCONT[nbco1][5]=NCONT[ii][5];     
		   NCONT[nbco1][6]=NCONT[ii][6];
		   NCONT[nbco1][7]=NCONT[ii][7]; 
		   NCONT[nbco1][8]=NCONT[ii][8];     		                     		    
                       		    
		   NOCONT[it][NBCONTCO[it]]=nbco1;
		   NOCONT[jt][NBCONTCO[jt]]=nbco1;
		   NBCONTCO[it]++;
		   NBCONTCO[jt]++; 		                       		    
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];

			 VALCOH[nbco1][0]=VALCOH[ii][0];
			 VALCOH[nbco1][1]=VALCOH[ii][1];
			 VALCOH[nbco1][2]=VALCOH[ii][2];	
			 VALCOH[nbco1][3]=VALCOH[ii][3];	
			 VALCOH[nbco1][4]=VALCOH[ii][4];
			 VALCOH[nbco1][5]=VALCOH[ii][5];
			 VALCOH[nbco1][6]=VALCOH[ii][6];	
			
			 VALAMO[nbco1][0]=VALAMO[ii][0];
			 VALAMO[nbco1][1]=VALAMO[ii][1];
			 VALAMO[nbco1][2]=VALAMO[ii][2];
			 VALAMO[nbco1][3]=VALAMO[ii][3];
			 VALAMO[nbco1][4]=VALAMO[ii][4];                      
			 VALAMO[nbco1][5]=VALAMO[ii][5];
			 VALAMO[nbco1][6]=VALAMO[ii][6];	                       
		   
			 dt=min(dt,sqrt(LIST_M[it]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_M[jt]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_I[it]/(VALCOH[nbco1][1]*LIST_R[it]*LIST_R[it])));
			 dt=min(dt,sqrt(LIST_I[jt]/(VALCOH[nbco1][1]*LIST_R[jt]*LIST_R[jt])));  

		  
		  nbco1++;
	  }
	  else{
	
					   int it=CONT[ii][0];
					   int jt=CONT[ii][1];
		          
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 	
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);  
				/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
						
						   TYPCO[nbco1]=0;  		
						 	  nbcoe++; 
						       CONT[nbco1][0]=it;
							   CONT[nbco1][1]=jt;                       
							   
							   DCONTX[nbco1]=dx;
							   DCONTY[nbco1]=dy;
							   DCONTZ[nbco1]=dz;	   
							   DCONT[nbco1] =nd2;
				               DCONTO[nbco1]=DCONTO[ii];	
  
		                         // Normale
							   
							   NCONT[nbco1][0]=dx/nd2;
							   NCONT[nbco1][1]=dy/nd2; 
							   NCONT[nbco1][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
							   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
							   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){
								  
								  if(NCONT[nbco1][0]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 1.;						
								  }else if(NCONT[nbco1][1]!=0.){
								   NCONT[nbco1][3] = 1.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 0.;									  
								  }else if(NCONT[nbco1][2]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 1.;
								   NCONT[nbco1][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco1][0]==0.){
									
									fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = 1./norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco1][3] = fors/norms;
								   NCONT[nbco1][4] = 1./norms;
								   NCONT[nbco1][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
								NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
								NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco1;
								NOCONT[jt][NBCONTCO[jt]]=nbco1;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco1]-(LIST_R[it]+LIST_R[jt]));
								 iEij=2.*(1.-nu*nu)/E;
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=2.5e8;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  				
						
						
						
						  nbco1++;		  
					  }	  
					  	  
		  
	  }
	  
  }	  

nbco=nbco1;

// Recherche nouveaux contacts
                        
	for(int ii=0;ii<vecsize;ii++){
		
			for(int lt=0;lt<nocoul[ii];lt++){
			int it = coul[ii][lt];

				for(int kt=0;kt<nonumc[ii];kt++){

					 for(int kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
					 
					  int jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
                           
							   bool booloc=0 ; 
								  
							   if(NBCONTCO[it]>0){
								 int ite=0;
								 while((ite<NBCONTCO[it])&&(!booloc)){
								   unsigned int numco=NOCONT[it][ite];
								 
								   if((numco<nbco1)&&(((CONT[numco][0]==it)&&(CONT[numco][1]==jt))||((CONT[numco][0]==jt)&&(CONT[numco][1]==it)))) booloc=1;			   	
								
								  ite++;						   
								 } 
							   }	
									
							   if(!booloc){  
                             //      cout<<"hello!"<<endl;
							   TYPCO[nbco]=0;
										  nbcoe++;
							   CONT[nbco][0]=it;
							   CONT[nbco][1]=jt;                       
							   
							   DCONTX[nbco]=dx;
							   DCONTY[nbco]=dy;
							   DCONTZ[nbco]=dz;	   
							   DCONT[nbco] =nd2;
							   DCONTO[nbco1]=r2;
							   
		                         // Normale
							   
							   NCONT[nbco][0]=dx/nd2;
							   NCONT[nbco][1]=dy/nd2; 
							   NCONT[nbco][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
							   if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
							   if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
								  
								  if(NCONT[nbco][0]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 1.;						
								  }else if(NCONT[nbco][1]!=0.){
								   NCONT[nbco][3] = 1.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 0.;									  
								  }else if(NCONT[nbco][2]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 1.;
								   NCONT[nbco][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco][0]==0.){
									
									fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco][3] = 1./norms;
									NCONT[nbco][4] = 1./norms;
									NCONT[nbco][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco][3] = fors/norms;
								   NCONT[nbco][4] = 1./norms;
								   NCONT[nbco][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
								NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
								NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco;
								NOCONT[jt][NBCONTCO[jt]]=nbco;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco]-(LIST_R[it]+LIST_R[jt]));
								 iEij=2.*(1.-nu*nu)/E;
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=2.5e8;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  		
															   
							   nbco++;						   
							   }														   
						   
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

if(numax>NMAXZ){
cout<<"Nombre maximal de contacts par sph atteint: "<<numax<<"/"<<NMAXZ<<" - stop "<<endl;
exit(0);
}

if(nbco>NMAXCONT){
cout<<"Nombre maximal de contacts sph-sph atteint:"<<nbco<<"/"<<NMAXCONT<<" - stop "<<endl;
exit(0);
}


dt=Cs*dt;
//cout<<"Pas de temps:"<<dt<<endl;
}

void selco_corci2(bool * TYPCO, R & dt, R epsi, unsigned int & nbco,int & nbcoh,int & nbcoe,int NB_SPH, int ** CONT, int vecsize, int vecsizex, int vecsizey, int vecsizez, R H_TOT, R V_TOT, R Z_TOT, int **  coul, int * nocoul, int ** numc, int * nonumc, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT, R * DCONTO,  R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R Pi, R Cs, R E1, R nu1, R E2, R nu2, bool * LIST_P, int NMAXCONT, int NMAXZ)
{

// Mise à blanc des contacts par corci
  for(int it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  }

R dx,dy,dz,nd2,r2,rij,hij,iEij,Eij;
R Kn,Kt;
R fors,norms;

dt=1e9;

// Traitement anciens contacts
unsigned int nbco1=0;
nbcoh=0;
nbcoe=0;
int numax=0;

  for(int ii=0;ii<nbco;ii++){
	  
	  if(TYPCO[ii]){
	       nbcoh++;
		   TYPCO[nbco1]=1;  	
		   
           CONT[nbco1][0]=CONT[ii][0];
		   CONT[nbco1][1]=CONT[ii][1];                       
		   
		   int it=CONT[nbco1][0];
		   int jt=CONT[nbco1][1];     
		    
   		   DCONTX[nbco1]=DCONTX[ii];
		   DCONTY[nbco1]=DCONTY[ii];
		   DCONTZ[nbco1]=DCONTZ[ii];		   
		   DCONT[nbco1] =DCONT[ii];		
		   DCONTO[nbco1]=DCONTO[ii];	

		   NCONT[nbco1][0]=NCONT[ii][0];
		   NCONT[nbco1][1]=NCONT[ii][1]; 
		   NCONT[nbco1][2]=NCONT[ii][2];
		   NCONT[nbco1][3]=NCONT[ii][3];
		   NCONT[nbco1][4]=NCONT[ii][4]; 
		   NCONT[nbco1][5]=NCONT[ii][5];     
		   NCONT[nbco1][6]=NCONT[ii][6];
		   NCONT[nbco1][7]=NCONT[ii][7]; 
		   NCONT[nbco1][8]=NCONT[ii][8];     		                     		    
                       		    
		   NOCONT[it][NBCONTCO[it]]=nbco1;
		   NOCONT[jt][NBCONTCO[jt]]=nbco1;
		   NBCONTCO[it]++;
		   NBCONTCO[jt]++; 	
			if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
			if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	                       		    

			 VALCOH[nbco1][0]=VALCOH[ii][0];
			 VALCOH[nbco1][1]=VALCOH[ii][1];
			 VALCOH[nbco1][2]=VALCOH[ii][2];	
			 VALCOH[nbco1][3]=VALCOH[ii][3];	
			 VALCOH[nbco1][4]=VALCOH[ii][4];
			 VALCOH[nbco1][5]=VALCOH[ii][5];
			 VALCOH[nbco1][6]=VALCOH[ii][6];	
			
			 VALAMO[nbco1][0]=VALAMO[ii][0];
			 VALAMO[nbco1][1]=VALAMO[ii][1];
			 VALAMO[nbco1][2]=VALAMO[ii][2];
			 VALAMO[nbco1][3]=VALAMO[ii][3];
			 VALAMO[nbco1][4]=VALAMO[ii][4];                      
			 VALAMO[nbco1][5]=VALAMO[ii][5];
			 VALAMO[nbco1][6]=VALAMO[ii][6];	                       
		   
			 dt=min(dt,sqrt(LIST_M[it]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_M[jt]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_I[it]/(VALCOH[nbco1][1]*LIST_R[it]*LIST_R[it])));
			 dt=min(dt,sqrt(LIST_I[jt]/(VALCOH[nbco1][1]*LIST_R[jt]*LIST_R[jt])));  

		  
		  nbco1++;
	  }
	  else{
	
					   int it=CONT[ii][0];
					   int jt=CONT[ii][1];
		          
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 	
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);  
				/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
						
						   TYPCO[nbco1]=0;  		
						 	  nbcoe++; 
						       CONT[nbco1][0]=it;
							   CONT[nbco1][1]=jt;                       
							   
							   DCONTX[nbco1]=dx;
							   DCONTY[nbco1]=dy;
							   DCONTZ[nbco1]=dz;	   
							   DCONT[nbco1] =nd2;
				               DCONTO[nbco1]=DCONTO[ii];	
  
		                         // Normale
							   
							   NCONT[nbco1][0]=dx/nd2;
							   NCONT[nbco1][1]=dy/nd2; 
							   NCONT[nbco1][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
							   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
							   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){
								  
								  if(NCONT[nbco1][0]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 1.;						
								  }else if(NCONT[nbco1][1]!=0.){
								   NCONT[nbco1][3] = 1.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 0.;									  
								  }else if(NCONT[nbco1][2]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 1.;
								   NCONT[nbco1][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco1][0]==0.){
									
									fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = 1./norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco1][3] = fors/norms;
								   NCONT[nbco1][4] = 1./norms;
								   NCONT[nbco1][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
								NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
								NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco1;
								NOCONT[jt][NBCONTCO[jt]]=nbco1;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
								if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
								if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco1]-(LIST_R[it]+LIST_R[jt]));
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
									
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=2.5e8;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  				
						
						
						
						  nbco1++;		  
					  }	  
					  	  
		  
	  }
	  
	  
  }	  

nbco=nbco1;
// Recherche nouveaux contacts

                       
	for(int ii=0;ii<vecsize;ii++){
		
			for(int lt=0;lt<nocoul[ii];lt++){
			int it = coul[ii][lt];

				for(int kt=0;kt<nonumc[ii];kt++){

					 for(int kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 

					 
					  int jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
							   bool booloc=0 ; 
								  
							   if(NBCONTCO[it]>0){
								 int ite=0;
								 while((ite<NBCONTCO[it])&&(!booloc)){
								   unsigned int numco=NOCONT[it][ite];
								 
								   if((numco<nbco1)&&(((CONT[numco][0]==it)&&(CONT[numco][1]==jt))||((CONT[numco][0]==jt)&&(CONT[numco][1]==it)))) booloc=1;			   	
								
								  ite++;						   
								 } 
							   }	
									
							   if(!booloc){  
                            
							   TYPCO[nbco]=0;
						       nbcoe++;
							   CONT[nbco][0]=it;
							   CONT[nbco][1]=jt;                       
							   
							   DCONTX[nbco]=dx;
							   DCONTY[nbco]=dy;
							   DCONTZ[nbco]=dz;	   
							   DCONT[nbco] =nd2;
							   DCONTO[nbco1]=r2;
							   
		                         // Normale
							   
							   NCONT[nbco][0]=dx/nd2;
							   NCONT[nbco][1]=dy/nd2; 
							   NCONT[nbco][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
							   if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
							   if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
								  
								  if(NCONT[nbco][0]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 1.;						
								  }else if(NCONT[nbco][1]!=0.){
								   NCONT[nbco][3] = 1.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 0.;									  
								  }else if(NCONT[nbco][2]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 1.;
								   NCONT[nbco][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco][0]==0.){
									
									fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco][3] = 1./norms;
									NCONT[nbco][4] = 1./norms;
									NCONT[nbco][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco][3] = fors/norms;
								   NCONT[nbco][4] = 1./norms;
								   NCONT[nbco][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
								NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
								NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco;
								NOCONT[jt][NBCONTCO[jt]]=nbco;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
								if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
								if(NBCONTCO[it]>numax) numax=NBCONTCO[it];		
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco]-(LIST_R[it]+LIST_R[jt]));
								 
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
								 
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=2.5e8;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  		
															   
							   nbco++;						   
							   }	
							   
							   
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

if(numax>NMAXZ){
cout<<"Nombre maximal de contacts par sph atteint: "<<numax<<"/"<<NMAXZ<<" - stop "<<endl;
exit(0);
}

if(nbco>NMAXCONT){
cout<<"Nombre maximal de contacts sph-sph atteint:"<<nbco<<"/"<<NMAXCONT<<" - stop "<<endl;
exit(0);
}


dt=Cs*dt;
//cout<<"Pas de temps:"<<dt<<endl;
//cout<<"nbco:"<<nbco<<endl;
}

void selco_corci2_mIexp(int & nint,int & nsoft,R Kcon, R * vs,R * LIST_VX, R * LIST_VY, R * LIST_VZ,bool * TYPCO, R & dt, R epsi, int & nbco,int & nbcoh,int & nbcoe, int NB_SPH, int ** CONT, int vecsize, int ** coul, int * nocoul, int ** numc, int * nonumc, R * LIST_IND, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ,R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, R E1, R nu1, R E2, R nu2, R Emu1, R rmu1, R Emu2, R rmu2, R amort, bool * LIST_P, int * LIST_B, int ite,  R GI, R GII, R snmax, R ssmax)
{

R ind1,ind2,indm;
R dx,dy,dz,nd2,r2,rij,hij,iEij,Eij;
R Kn,Kt;
R fors,norms;
R Emui,rmui;
R S,kkn,kkn0,un,unc,unm,ray;

int ii;
R dxij,dyij,dzij;
R I,kkt,kkt0,us,usc,usm,D;

unm=2.*GI/snmax;
usm=2.*GII/ssmax;

// Mise à blanc des contacts par corci
  for(int it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  }

R v2,vn;
dt=1e9;

// Traitement anciens contacts
int nbco1=0;
nbcoh=0;
nbcoe=0;
nint=0;
nsoft=0;

  for(ii=0;ii<nbco;ii++){
	 
	  if(TYPCO[ii]){
	
            nbcoh++;
		   TYPCO[nbco1]=1;  	
		   
           CONT[nbco1][0]=CONT[ii][0];
		   CONT[nbco1][1]=CONT[ii][1];                       
		   
		   int it=CONT[nbco1][0];
		   int jt=CONT[nbco1][1];     
		    
   		   DCONTX[nbco1]=DCONTX[ii];
		   DCONTY[nbco1]=DCONTY[ii];
		   DCONTZ[nbco1]=DCONTZ[ii];		   
		   DCONT[nbco1] =DCONT[ii];		
		   DCONTO[nbco1]=DCONTO[ii];	

		   NCONT[nbco1][0]=NCONT[ii][0];
		   NCONT[nbco1][1]=NCONT[ii][1]; 
		   NCONT[nbco1][2]=NCONT[ii][2];
		   NCONT[nbco1][3]=NCONT[ii][3];
		   NCONT[nbco1][4]=NCONT[ii][4]; 
		   NCONT[nbco1][5]=NCONT[ii][5];     
		   NCONT[nbco1][6]=NCONT[ii][6];
		   NCONT[nbco1][7]=NCONT[ii][7]; 
		   NCONT[nbco1][8]=NCONT[ii][8];     		                     		    
                       		    
		   NOCONT[it][NBCONTCO[it]]=nbco1;
		   NOCONT[jt][NBCONTCO[jt]]=nbco1;
		   NBCONTCO[it]++;
		   NBCONTCO[jt]++; 		                       		    

			 VALCOH[nbco1][0]=VALCOH[ii][0];
			 VALCOH[nbco1][1]=VALCOH[ii][1];
			 VALCOH[nbco1][2]=VALCOH[ii][2];	
			 VALCOH[nbco1][3]=VALCOH[ii][3];	
			 VALCOH[nbco1][4]=VALCOH[ii][4];
			 VALCOH[nbco1][5]=VALCOH[ii][5];
			 VALCOH[nbco1][6]=VALCOH[ii][6];	
			

 // Cas interfaciel		
              
         
           
 			if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
            
			rmui=(rmu1+rmu2)/2.;
            Emui=(Emu1+Emu2)/2.;
            
		    ray=rmui*(LIST_R[it]+LIST_R[jt])/2.;

			ind1=LIST_IND[it];
			ind2=LIST_IND[jt];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;
            
			dxij=(LIST_VX[it]-LIST_VX[jt]);        	           
		    dyij=(LIST_VY[it]-LIST_VY[jt]);
		    dzij=(LIST_VZ[it]-LIST_VZ[jt]);

		    dx=LIST_X[it]-LIST_X[jt]; 
		    dy=LIST_Y[it]-LIST_Y[jt]; 			
		    dz=LIST_Z[it]-LIST_Z[jt]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

			un=(dx/nd2)*(dx-DCONTX[nbco1]);
			un+=(dy/nd2)*(dy-DCONTY[nbco1]);	
			un+=(dz/nd2)*(dz-DCONTZ[nbco1]);           

     	    kkn0=Emui*S/nd2;
		    unc=snmax*S/kkn0;
		    
		    nint++;
			
			if(un<=unc) {D=0.;}
			else if((un>unc)&&(un<15.*unc)) {D=1.-exp(-(un-unc)/unc);nsoft++;}
			else if(un>=15.*unc) {D=1.;TYPCO[nbco1]=0;
			
		    LIST_B[it]=2;
		    LIST_B[jt]=2;  
		    
		    } 			
	
			kkn=(1.-D)*kkn0;					
		        VALCOH[nbco1][0]=kkn;
 		
            }

      
             dt=min(dt,sqrt(LIST_M[it]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_M[jt]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_I[it]/(VALCOH[nbco1][1]*LIST_R[it]*LIST_R[it])));
		 	 dt=min(dt,sqrt(LIST_I[jt]/(VALCOH[nbco1][1]*LIST_R[jt]*LIST_R[jt])));  
		    
		  nbco1++;
	  
      }	 
 	  /* else{
	
					   int it=CONT[ii][0];
					   int jt=CONT[ii][1];
		          
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 	
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);  
				/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
						
						   TYPCO[nbco1]=0;  		
						 	  nbcoe++; 
						       CONT[nbco1][0]=it;
							   CONT[nbco1][1]=jt;                       
							   
							   DCONTX[nbco1]=dx;
							   DCONTY[nbco1]=dy;
							   DCONTZ[nbco1]=dz;	   
							   DCONT[nbco1] =nd2;
				               DCONTO[nbco1]=DCONTO[ii];	
  
		                         // Normale
							   
							   NCONT[nbco1][0]=dx/nd2;
							   NCONT[nbco1][1]=dy/nd2; 
							   NCONT[nbco1][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
							   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
							   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){
								  
								  if(NCONT[nbco1][0]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 1.;						
								  }else if(NCONT[nbco1][1]!=0.){
								   NCONT[nbco1][3] = 1.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 0.;									  
								  }else if(NCONT[nbco1][2]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 1.;
								   NCONT[nbco1][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco1][0]==0.){
									
									fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = 1./norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco1][3] = fors/norms;
								   NCONT[nbco1][4] = 1./norms;
								   NCONT[nbco1][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
								NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
								NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco1;
								NOCONT[jt][NBCONTCO[jt]]=nbco1;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
	//							if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
		//						if(NBCONTCO[it]>numax) numax=NBCONTCO[it];		
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco1]-(LIST_R[it]+LIST_R[jt]));
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
									
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=Kcon;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  				
						
						
						
						  nbco1++;		  
					  }	  
					  	  
		  
	  }*/
      }	   

nbco=nbco1;

/*
// Nouveau contacts
                      
	for(int ii=0;ii<vecsize;ii++){
		
			for(int lt=0;lt<nocoul[ii];lt++){
			int it = coul[ii][lt];

				for(int kt=0;kt<nonumc[ii];kt++){

					 for(int kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
					 
					  int jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
                            
							   bool booloc=0 ; 
								  
							   if(NBCONTCO[it]>0){
								 int ite=0;
								 while((ite<NBCONTCO[it])&&(!booloc)){
								   unsigned int numco=NOCONT[it][ite];
								 
								   if((numco<nbco1)&&(((CONT[numco][0]==it)&&(CONT[numco][1]==jt))||((CONT[numco][0]==jt)&&(CONT[numco][1]==it)))) booloc=1;			   	
								
								  ite++;						   
								 } 
							   }	
							   
							   if(!booloc){  
                            
							   TYPCO[nbco]=0;
						       nbcoe++;
							   CONT[nbco][0]=it;
							   CONT[nbco][1]=jt;                       
							   
							   DCONTX[nbco]=dx;
							   DCONTY[nbco]=dy;
							   DCONTZ[nbco]=dz;	   
							   DCONT[nbco] =nd2;
							   DCONTO[nbco1]=r2;
							   
		                         // Normale
							   
							   NCONT[nbco][0]=dx/nd2;
							   NCONT[nbco][1]=dy/nd2; 
							   NCONT[nbco][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
							   if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
							   if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
								  
								  if(NCONT[nbco][0]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 1.;						
								  }else if(NCONT[nbco][1]!=0.){
								   NCONT[nbco][3] = 1.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 0.;									  
								  }else if(NCONT[nbco][2]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 1.;
								   NCONT[nbco][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco][0]==0.){
									
									fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco][3] = 1./norms;
									NCONT[nbco][4] = 1./norms;
									NCONT[nbco][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco][3] = fors/norms;
								   NCONT[nbco][4] = 1./norms;
								   NCONT[nbco][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
								NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
								NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco;
								NOCONT[jt][NBCONTCO[jt]]=nbco;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
			//					if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
				//				if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco]-(LIST_R[it]+LIST_R[jt]));
								 
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
								 
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=Kcon;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  		
															   
							   nbco++;						   
							   }														   
						
						  
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

*/

/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
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

//cout<<"Pas de temps:"<<dt<<endl;
//cout<<"nbco:"<<nbco<<endl;
}









void selco_corci2_mI(int & nint,int & nsoft,R Kcon, R * vs,R * LIST_VX, R * LIST_VY, R * LIST_VZ,bool * TYPCO, R & dt, R epsi, int & nbco,int & nbcoh,int & nbcoe, int NB_SPH, int ** CONT, int vecsize, int ** coul, int * nocoul, int ** numc, int * nonumc, R * LIST_IND, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ,R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, R E1, R nu1, R E2, R nu2, R Emu1, R rmu1, R Emu2, R rmu2, R amort, bool * LIST_P, int * LIST_B, int ite,  R GI, R GII, R snmax, R ssmax)
{

R ind1,ind2,indm;
R dx,dy,dz,nd2,r2,rij,hij,iEij,Eij;
R Kn,Kt;
R fors,norms;
R Emui,rmui;
R S,kkn,kkn0,un,unc,unm,ray;

int ii;
R dxij,dyij,dzij;
R I,kkt,kkt0,us,usc,usm,D;

unm=2.*GI/snmax;
usm=2.*GII/ssmax;

// Mise à blanc des contacts par corci
  for(int it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  }

R v2,vn;
dt=1e9;

// Traitement anciens contacts
int nbco1=0;
nbcoh=0;
nbcoe=0;
nint=0;
nsoft=0;

  for(ii=0;ii<nbco;ii++){
	 
	  if(TYPCO[ii]){
	
            nbcoh++;
		   TYPCO[nbco1]=1;  	
		   
           CONT[nbco1][0]=CONT[ii][0];
		   CONT[nbco1][1]=CONT[ii][1];                       
		   
		   int it=CONT[nbco1][0];
		   int jt=CONT[nbco1][1];     
		    
   		   DCONTX[nbco1]=DCONTX[ii];
		   DCONTY[nbco1]=DCONTY[ii];
		   DCONTZ[nbco1]=DCONTZ[ii];		   
		   DCONT[nbco1] =DCONT[ii];		
		   DCONTO[nbco1]=DCONTO[ii];	

		   NCONT[nbco1][0]=NCONT[ii][0];
		   NCONT[nbco1][1]=NCONT[ii][1]; 
		   NCONT[nbco1][2]=NCONT[ii][2];
		   NCONT[nbco1][3]=NCONT[ii][3];
		   NCONT[nbco1][4]=NCONT[ii][4]; 
		   NCONT[nbco1][5]=NCONT[ii][5];     
		   NCONT[nbco1][6]=NCONT[ii][6];
		   NCONT[nbco1][7]=NCONT[ii][7]; 
		   NCONT[nbco1][8]=NCONT[ii][8];     		                     		    
                       		    
		   NOCONT[it][NBCONTCO[it]]=nbco1;
		   NOCONT[jt][NBCONTCO[jt]]=nbco1;
		   NBCONTCO[it]++;
		   NBCONTCO[jt]++; 		                       		    

			 VALCOH[nbco1][0]=VALCOH[ii][0];
			 VALCOH[nbco1][1]=VALCOH[ii][1];
			 VALCOH[nbco1][2]=VALCOH[ii][2];	
			 VALCOH[nbco1][3]=VALCOH[ii][3];	
			 VALCOH[nbco1][4]=VALCOH[ii][4];
			 VALCOH[nbco1][5]=VALCOH[ii][5];
			 VALCOH[nbco1][6]=VALCOH[ii][6];	
			

 // Cas interfaciel		
              
         
           
 			if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
            
			rmui=(rmu1+rmu2)/2.;
            Emui=(Emu1+Emu2)/2.;
            
		    ray=rmui*(LIST_R[it]+LIST_R[jt])/2.;

			ind1=LIST_IND[it];
			ind2=LIST_IND[jt];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;
            
			dxij=(LIST_VX[it]-LIST_VX[jt]);        	           
		    dyij=(LIST_VY[it]-LIST_VY[jt]);
		    dzij=(LIST_VZ[it]-LIST_VZ[jt]);

		    dx=LIST_X[it]-LIST_X[jt]; 
		    dy=LIST_Y[it]-LIST_Y[jt]; 			
		    dz=LIST_Z[it]-LIST_Z[jt]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

			un=(dx/nd2)*(dx-DCONTX[nbco1]);
			un+=(dy/nd2)*(dy-DCONTY[nbco1]);	
			un+=(dz/nd2)*(dz-DCONTZ[nbco1]);           

     	    kkn0=Emui*S/nd2;
		    unc=snmax*S/kkn0;
		    
		    nint++;
			
			if(un<=unc) {D=0.;}
			else if((un>unc)&&(un<unm)) {D=unm*(un-unc)/(un*(unm-unc));nsoft++;}
			else if(un>=unm) {D=1.;TYPCO[nbco1]=0;
				
		    LIST_B[it]=2;
		    LIST_B[jt]=2;  
		    
		    } 			
	
			kkn=(1.-D)*kkn0;					
		        VALCOH[nbco1][0]=kkn;
 		
            }

      
             dt=min(dt,sqrt(LIST_M[it]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_M[jt]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_I[it]/(VALCOH[nbco1][1]*LIST_R[it]*LIST_R[it])));
		 	 dt=min(dt,sqrt(LIST_I[jt]/(VALCOH[nbco1][1]*LIST_R[jt]*LIST_R[jt])));  
		    
		  nbco1++;
	  
      }	 
 	  /* else{
	
					   int it=CONT[ii][0];
					   int jt=CONT[ii][1];
		          
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 	
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);  
				/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
						
						   TYPCO[nbco1]=0;  		
						 	  nbcoe++; 
						       CONT[nbco1][0]=it;
							   CONT[nbco1][1]=jt;                       
							   
							   DCONTX[nbco1]=dx;
							   DCONTY[nbco1]=dy;
							   DCONTZ[nbco1]=dz;	   
							   DCONT[nbco1] =nd2;
				               DCONTO[nbco1]=DCONTO[ii];	
  
		                         // Normale
							   
							   NCONT[nbco1][0]=dx/nd2;
							   NCONT[nbco1][1]=dy/nd2; 
							   NCONT[nbco1][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
							   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
							   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){
								  
								  if(NCONT[nbco1][0]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 1.;						
								  }else if(NCONT[nbco1][1]!=0.){
								   NCONT[nbco1][3] = 1.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 0.;									  
								  }else if(NCONT[nbco1][2]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 1.;
								   NCONT[nbco1][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco1][0]==0.){
									
									fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = 1./norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco1][3] = fors/norms;
								   NCONT[nbco1][4] = 1./norms;
								   NCONT[nbco1][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
								NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
								NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco1;
								NOCONT[jt][NBCONTCO[jt]]=nbco1;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
	//							if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
		//						if(NBCONTCO[it]>numax) numax=NBCONTCO[it];		
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco1]-(LIST_R[it]+LIST_R[jt]));
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
									
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=Kcon;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  				
						
						
						
						  nbco1++;		  
					  }	  
					  	  
		  
	  }*/
      }	   

nbco=nbco1;

/*
// Nouveau contacts
                      
	for(int ii=0;ii<vecsize;ii++){
		
			for(int lt=0;lt<nocoul[ii];lt++){
			int it = coul[ii][lt];

				for(int kt=0;kt<nonumc[ii];kt++){

					 for(int kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
					 
					  int jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
                            
							   bool booloc=0 ; 
								  
							   if(NBCONTCO[it]>0){
								 int ite=0;
								 while((ite<NBCONTCO[it])&&(!booloc)){
								   unsigned int numco=NOCONT[it][ite];
								 
								   if((numco<nbco1)&&(((CONT[numco][0]==it)&&(CONT[numco][1]==jt))||((CONT[numco][0]==jt)&&(CONT[numco][1]==it)))) booloc=1;			   	
								
								  ite++;						   
								 } 
							   }	
							   
							   if(!booloc){  
                            
							   TYPCO[nbco]=0;
						       nbcoe++;
							   CONT[nbco][0]=it;
							   CONT[nbco][1]=jt;                       
							   
							   DCONTX[nbco]=dx;
							   DCONTY[nbco]=dy;
							   DCONTZ[nbco]=dz;	   
							   DCONT[nbco] =nd2;
							   DCONTO[nbco1]=r2;
							   
		                         // Normale
							   
							   NCONT[nbco][0]=dx/nd2;
							   NCONT[nbco][1]=dy/nd2; 
							   NCONT[nbco][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
							   if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
							   if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
								  
								  if(NCONT[nbco][0]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 1.;						
								  }else if(NCONT[nbco][1]!=0.){
								   NCONT[nbco][3] = 1.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 0.;									  
								  }else if(NCONT[nbco][2]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 1.;
								   NCONT[nbco][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco][0]==0.){
									
									fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco][3] = 1./norms;
									NCONT[nbco][4] = 1./norms;
									NCONT[nbco][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco][3] = fors/norms;
								   NCONT[nbco][4] = 1./norms;
								   NCONT[nbco][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
								NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
								NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco;
								NOCONT[jt][NBCONTCO[jt]]=nbco;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
			//					if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
				//				if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco]-(LIST_R[it]+LIST_R[jt]));
								 
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
								 
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=Kcon;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  		
															   
							   nbco++;						   
							   }														   
						
						  
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

*/

/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
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

//cout<<"Pas de temps:"<<dt<<endl;
//cout<<"nbco:"<<nbco<<endl;
}










void selco_corci2_mII(int & nint,int & nsoft,R Kcon, R * vs,R * LIST_VX, R * LIST_VY, R * LIST_VZ,bool * TYPCO, R & dt, R epsi, int & nbco,int & nbcoh,int & nbcoe, int NB_SPH, int ** CONT, int vecsize, int ** coul, int * nocoul, int ** numc, int * nonumc, R * LIST_IND, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ,R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, R E1, R nu1, R E2, R nu2, R Emu1, R rmu1, R Emu2, R rmu2, R amort, bool * LIST_P, int * LIST_B, int ite,  R GI, R GII, R snmax, R ssmax)
{

R ind1,ind2,indm;
R dx,dy,dz,nd2,r2,rij,hij,iEij,Eij;
R Kn,Kt;
R fors,norms;
R Emui,rmui;
R S,kkn,kkn0,un,unc,unm,ray;

int ii;
R dxij,dyij,dzij;
R I,kkt,kkt0,us,usc,usm,D;

// Mise à blanc des contacts par corci
  for(int it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  }

R v2,vn;
R dta=dt;
dt=1e9;

// Traitement anciens contacts
int nbco1=0;
nbcoh=0;
nbcoe=0;
nint=0;
nsoft=0;

  for(ii=0;ii<nbco;ii++){
	 
	  if(TYPCO[ii]){
	
            nbcoh++;
		   TYPCO[nbco1]=1;  	
		   
           CONT[nbco1][0]=CONT[ii][0];
		   CONT[nbco1][1]=CONT[ii][1];                       
		   
		   int it=CONT[nbco1][0];
		   int jt=CONT[nbco1][1];     
		    
   		   DCONTX[nbco1]=DCONTX[ii];
		   DCONTY[nbco1]=DCONTY[ii];
		   DCONTZ[nbco1]=DCONTZ[ii];		   
		   DCONT[nbco1] =DCONT[ii];		
		   DCONTO[nbco1]=DCONTO[ii];	

		   NCONT[nbco1][0]=NCONT[ii][0];
		   NCONT[nbco1][1]=NCONT[ii][1]; 
		   NCONT[nbco1][2]=NCONT[ii][2];
		   NCONT[nbco1][3]=NCONT[ii][3];
		   NCONT[nbco1][4]=NCONT[ii][4]; 
		   NCONT[nbco1][5]=NCONT[ii][5];     
		   NCONT[nbco1][6]=NCONT[ii][6];
		   NCONT[nbco1][7]=NCONT[ii][7]; 
		   NCONT[nbco1][8]=NCONT[ii][8];     		                     		    
                       		    
		   NOCONT[it][NBCONTCO[it]]=nbco1;
		   NOCONT[jt][NBCONTCO[jt]]=nbco1;
		   NBCONTCO[it]++;
		   NBCONTCO[jt]++; 		                       		    

			 VALCOH[nbco1][0]=VALCOH[ii][0];
			 VALCOH[nbco1][1]=VALCOH[ii][1];
			 VALCOH[nbco1][2]=VALCOH[ii][2];	
			 VALCOH[nbco1][3]=VALCOH[ii][3];	
			 VALCOH[nbco1][4]=VALCOH[ii][4];
			 VALCOH[nbco1][5]=VALCOH[ii][5];
			 VALCOH[nbco1][6]=VALCOH[ii][6];	
			

 // Cas interfaciel		
              
         
           
 			if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
            
			rmui=(rmu1+rmu2)/2.;
            Emui=(Emu1+Emu2)/2.;
            
		    ray=rmui*(LIST_R[it]+LIST_R[jt])/2.;

			ind1=LIST_IND[it];
			ind2=LIST_IND[jt];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;
            
			dxij=(LIST_VX[it]-LIST_VX[jt]);        	           
		    dyij=(LIST_VY[it]-LIST_VY[jt]);
		    dzij=(LIST_VZ[it]-LIST_VZ[jt]);

		    dx=LIST_X[it]-LIST_X[jt]; 
		    dy=LIST_Y[it]-LIST_Y[jt]; 			
		    dz=LIST_Z[it]-LIST_Z[jt]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

		    vn=(dx/nd2)*dxij+(dy/nd2)*dyij+(dz/nd2)*dzij;
			v2=(dxij*dxij+dyij*dyij+dzij*dzij-vn*vn);	
		    if(v2>0) {vs[nbco1]+=sqrt(v2);}	
		    
		    us=vs[nbco1]*dta;
			
     		kkt0=12.*Emui*I/(pow(nd2,3));
			
		    usc=ssmax*S/kkt0;
			usm=2.*GII/ssmax;
			nint++;
			
			if(us<=usc) {D=0.;}
			else if((us>usc)&&(us<usm)) {D=usm*(us-usc)/(us*(usm-usc));nsoft++;}
			else if(us>=usm) {D=1.;TYPCO[nbco1]=0;LIST_B[it]=2;LIST_B[jt]=2;}	
	
			kkt=(1.-D)*kkt0;					
			VALCOH[nbco1][1]=kkt;
 		
            }

      
             dt=min(dt,sqrt(LIST_M[it]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_M[jt]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_I[it]/(VALCOH[nbco1][1]*LIST_R[it]*LIST_R[it])));
		 	 dt=min(dt,sqrt(LIST_I[jt]/(VALCOH[nbco1][1]*LIST_R[jt]*LIST_R[jt])));  
		    
		  nbco1++;
	  
      }	 
 	  /* else{
	
					   int it=CONT[ii][0];
					   int jt=CONT[ii][1];
		          
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 	
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);  
				/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
						
						   TYPCO[nbco1]=0;  		
						 	  nbcoe++; 
						       CONT[nbco1][0]=it;
							   CONT[nbco1][1]=jt;                       
							   
							   DCONTX[nbco1]=dx;
							   DCONTY[nbco1]=dy;
							   DCONTZ[nbco1]=dz;	   
							   DCONT[nbco1] =nd2;
				               DCONTO[nbco1]=DCONTO[ii];	
  
		                         // Normale
							   
							   NCONT[nbco1][0]=dx/nd2;
							   NCONT[nbco1][1]=dy/nd2; 
							   NCONT[nbco1][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
							   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
							   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){
								  
								  if(NCONT[nbco1][0]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 1.;						
								  }else if(NCONT[nbco1][1]!=0.){
								   NCONT[nbco1][3] = 1.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 0.;									  
								  }else if(NCONT[nbco1][2]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 1.;
								   NCONT[nbco1][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco1][0]==0.){
									
									fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = 1./norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco1][3] = fors/norms;
								   NCONT[nbco1][4] = 1./norms;
								   NCONT[nbco1][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
								NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
								NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco1;
								NOCONT[jt][NBCONTCO[jt]]=nbco1;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
	//							if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
		//						if(NBCONTCO[it]>numax) numax=NBCONTCO[it];		
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco1]-(LIST_R[it]+LIST_R[jt]));
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
									
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=Kcon;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  				
						
						
						
						  nbco1++;		  
					  }	  
					  	  
		  
	  }*/
      }	   

nbco=nbco1;

/*
// Nouveau contacts
                      
	for(int ii=0;ii<vecsize;ii++){
		
			for(int lt=0;lt<nocoul[ii];lt++){
			int it = coul[ii][lt];

				for(int kt=0;kt<nonumc[ii];kt++){

					 for(int kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
					 
					  int jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
                            
							   bool booloc=0 ; 
								  
							   if(NBCONTCO[it]>0){
								 int ite=0;
								 while((ite<NBCONTCO[it])&&(!booloc)){
								   unsigned int numco=NOCONT[it][ite];
								 
								   if((numco<nbco1)&&(((CONT[numco][0]==it)&&(CONT[numco][1]==jt))||((CONT[numco][0]==jt)&&(CONT[numco][1]==it)))) booloc=1;			   	
								
								  ite++;						   
								 } 
							   }	
							   
							   if(!booloc){  
                            
							   TYPCO[nbco]=0;
						       nbcoe++;
							   CONT[nbco][0]=it;
							   CONT[nbco][1]=jt;                       
							   
							   DCONTX[nbco]=dx;
							   DCONTY[nbco]=dy;
							   DCONTZ[nbco]=dz;	   
							   DCONT[nbco] =nd2;
							   DCONTO[nbco1]=r2;
							   
		                         // Normale
							   
							   NCONT[nbco][0]=dx/nd2;
							   NCONT[nbco][1]=dy/nd2; 
							   NCONT[nbco][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
							   if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
							   if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
								  
								  if(NCONT[nbco][0]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 1.;						
								  }else if(NCONT[nbco][1]!=0.){
								   NCONT[nbco][3] = 1.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 0.;									  
								  }else if(NCONT[nbco][2]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 1.;
								   NCONT[nbco][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco][0]==0.){
									
									fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco][3] = 1./norms;
									NCONT[nbco][4] = 1./norms;
									NCONT[nbco][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco][3] = fors/norms;
								   NCONT[nbco][4] = 1./norms;
								   NCONT[nbco][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
								NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
								NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco;
								NOCONT[jt][NBCONTCO[jt]]=nbco;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
			//					if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
				//				if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco]-(LIST_R[it]+LIST_R[jt]));
								 
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
								 
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=Kcon;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  		
															   
							   nbco++;						   
							   }														   
						
						  
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

*/

/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
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

//cout<<"Pas de temps:"<<dt<<endl;
//cout<<"nbco:"<<nbco<<endl;
}




void selco_corci2_mm(int & nint,int & nsoft,R Kcon, R * vs,R * LIST_VX, R * LIST_VY, R * LIST_VZ,bool * TYPCO, R & dt, R epsi, int & nbco,int & nbcoh,int & nbcoe, int NB_SPH, int ** CONT, int vecsize, int ** coul, int * nocoul, int ** numc, int * nonumc, R * LIST_IND, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_M, R * LIST_I, R * DCONT,R * DCONTO, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ,R * DCONTXO, R * DCONTYO, R * DCONTZO, R Pi, R Cs, R E1, R nu1, R E2, R nu2, R Emu1, R rmu1, R Emu2, R rmu2, R amort, bool * LIST_P, int * LIST_B, int ite,  R GI, R GII, R snmax, R ssmax)
{

R ind1,ind2,indm;
R dx,dy,dz,nd2,r2,rij,hij,iEij,Eij;
R Kn,Kt;
R fors,norms;
R Emui,rmui;
R S,kkn,kkn0,un,unc,unm,ray;

int ii;
R dxij,dyij,dzij;
R I,kkt,kkt0,us,usc,usm,D;

R ue,uem,uec,eta,alpha;

rmui=(rmu1+rmu2)/2.;
Emui=(Emu1+Emu2)/2.;
unm=2.*GI/snmax;
usm=2.*GII/ssmax;
alpha=1.;

// Mise à blanc des contacts par corci
  for(int it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  }

R v2,vn;
R dta=dt;
dt=1e9;

// Traitement anciens contacts
int nbco1=0;
nbcoh=0;
nbcoe=0;
nint=0;
nsoft=0;

  for(ii=0;ii<nbco;ii++){
	 
	  if(TYPCO[ii]){
	
            nbcoh++;
		   TYPCO[nbco1]=1;  	
		   
           CONT[nbco1][0]=CONT[ii][0];
		   CONT[nbco1][1]=CONT[ii][1];                       
		   
		   int it=CONT[nbco1][0];
		   int jt=CONT[nbco1][1];     
		    
   		   DCONTX[nbco1]=DCONTX[ii];
		   DCONTY[nbco1]=DCONTY[ii];
		   DCONTZ[nbco1]=DCONTZ[ii];		   
		   DCONT[nbco1] =DCONT[ii];		
		   DCONTO[nbco1]=DCONTO[ii];	

		   NCONT[nbco1][0]=NCONT[ii][0];
		   NCONT[nbco1][1]=NCONT[ii][1]; 
		   NCONT[nbco1][2]=NCONT[ii][2];
		   NCONT[nbco1][3]=NCONT[ii][3];
		   NCONT[nbco1][4]=NCONT[ii][4]; 
		   NCONT[nbco1][5]=NCONT[ii][5];     
		   NCONT[nbco1][6]=NCONT[ii][6];
		   NCONT[nbco1][7]=NCONT[ii][7]; 
		   NCONT[nbco1][8]=NCONT[ii][8];     		                     		    
                       		    
		   NOCONT[it][NBCONTCO[it]]=nbco1;
		   NOCONT[jt][NBCONTCO[jt]]=nbco1;
		   NBCONTCO[it]++;
		   NBCONTCO[jt]++; 		                       		    

			 VALCOH[nbco1][0]=VALCOH[ii][0];
			 VALCOH[nbco1][1]=VALCOH[ii][1];
			 VALCOH[nbco1][2]=VALCOH[ii][2];	
			 VALCOH[nbco1][3]=VALCOH[ii][3];	
			 VALCOH[nbco1][4]=VALCOH[ii][4];
			 VALCOH[nbco1][5]=VALCOH[ii][5];
			 VALCOH[nbco1][6]=VALCOH[ii][6];	
			

 // Cas interfaciel		
              
         
           
 			if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
                        
		    ray=rmui*(LIST_R[it]+LIST_R[jt])/2.;

			ind1=LIST_IND[it];
			ind2=LIST_IND[jt];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;
            
			dxij=(LIST_VX[it]-LIST_VX[jt]);        	           
		    dyij=(LIST_VY[it]-LIST_VY[jt]);
		    dzij=(LIST_VZ[it]-LIST_VZ[jt]);

		    dx=LIST_X[it]-LIST_X[jt]; 
		    dy=LIST_Y[it]-LIST_Y[jt]; 			
		    dz=LIST_Z[it]-LIST_Z[jt]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

		    vn=(dx/nd2)*dxij+(dy/nd2)*dyij+(dz/nd2)*dzij;
			v2=(dxij*dxij+dyij*dyij+dzij*dzij-vn*vn);	
		    if(v2>0) {vs[nbco1]+=sqrt(v2);}	
		    
		    us=vs[nbco1]*dta;			

			un=(dx/nd2)*(dx-DCONTX[nbco1]);
			un+=(dy/nd2)*(dy-DCONTY[nbco1]);	
			un+=(dz/nd2)*(dz-DCONTZ[nbco1]);			

			ue=sqrt(un*un+us*us);			
			eta=(un>0)?(us/un):0.;	
			
     		kkt0=12.*Emui*I/(pow(nd2,3));
    		kkn0=Emui*S/nd2;

            unc=snmax*S/kkn0;
		    usc=ssmax*S/kkt0;
		    
			uec=unc*usc*sqrt((1+eta*eta)/(usc*usc+eta*eta*unc*unc));
			uem=(2.*(1+eta*eta)/uec)*pow(pow(kkn0/GI,alpha)+pow(eta*eta*kkt0/GII,alpha),-1/alpha);			    
		    
			nint++;
			if(ue<=uec) {D=0.;}
			else if((ue>uec)&&(ue<uem)) {D=uem*(ue-uec)/(ue*(uem-uec));nsoft++;}
			else if(ue>=uem) {D=1.;TYPCO[nbco1]=0;
								
		    # pragma omp critical
		    { 	//dangereux
		    LIST_B[it]=2;
		    LIST_B[jt]=2;  
		    }			
		    
		    } 	
	
			kkn=(1.-D)*kkn0;	
			kkt=(1.-D)*kkt0;	
					
			VALCOH[nbco1][0]=kkn;
			VALCOH[nbco1][1]=kkt;
 		
            }

      
             dt=min(dt,sqrt(LIST_M[it]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_M[jt]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_I[it]/(VALCOH[nbco1][1]*LIST_R[it]*LIST_R[it])));
		 	 dt=min(dt,sqrt(LIST_I[jt]/(VALCOH[nbco1][1]*LIST_R[jt]*LIST_R[jt])));  
		    
		  nbco1++;
	  
      }	 
 	  /* else{
	
					   int it=CONT[ii][0];
					   int jt=CONT[ii][1];
		          
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 	
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);  
				/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
						
						   TYPCO[nbco1]=0;  		
						 	  nbcoe++; 
						       CONT[nbco1][0]=it;
							   CONT[nbco1][1]=jt;                       
							   
							   DCONTX[nbco1]=dx;
							   DCONTY[nbco1]=dy;
							   DCONTZ[nbco1]=dz;	   
							   DCONT[nbco1] =nd2;
				               DCONTO[nbco1]=DCONTO[ii];	
  
		                         // Normale
							   
							   NCONT[nbco1][0]=dx/nd2;
							   NCONT[nbco1][1]=dy/nd2; 
							   NCONT[nbco1][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco1][0])<1e-18) NCONT[nbco1][0]=0.;
							   if(fabs(NCONT[nbco1][1])<1e-18) NCONT[nbco1][1]=0.;
							   if(fabs(NCONT[nbco1][2])<1e-18) NCONT[nbco1][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco1][0]*NCONT[nbco1][1]==0.)&&(NCONT[nbco1][1]*NCONT[nbco1][2]==0.)&&(NCONT[nbco1][2]*NCONT[nbco1][0]==0.)){
								  
								  if(NCONT[nbco1][0]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 1.;						
								  }else if(NCONT[nbco1][1]!=0.){
								   NCONT[nbco1][3] = 1.;
								   NCONT[nbco1][4] = 0.;
								   NCONT[nbco1][5] = 0.;									  
								  }else if(NCONT[nbco1][2]!=0.){
								   NCONT[nbco1][3] = 0.;
								   NCONT[nbco1][4] = 1.;
								   NCONT[nbco1][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco1][0]==0.){
									
									fors  = -(NCONT[nbco1][1])/NCONT[nbco1][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco1][3] = 1./norms;
									NCONT[nbco1][4] = 1./norms;
									NCONT[nbco1][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco1][1]+NCONT[nbco1][2])/NCONT[nbco1][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco1][3] = fors/norms;
								   NCONT[nbco1][4] = 1./norms;
								   NCONT[nbco1][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco1][6] = NCONT[nbco1][1]*NCONT[nbco1][5] - NCONT[nbco1][2]*NCONT[nbco1][4];
								NCONT[nbco1][7] = NCONT[nbco1][2]*NCONT[nbco1][3] - NCONT[nbco1][0]*NCONT[nbco1][5];
								NCONT[nbco1][8] = NCONT[nbco1][0]*NCONT[nbco1][4] - NCONT[nbco1][1]*NCONT[nbco1][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco1;
								NOCONT[jt][NBCONTCO[jt]]=nbco1;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
	//							if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
		//						if(NBCONTCO[it]>numax) numax=NBCONTCO[it];		
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco1]-(LIST_R[it]+LIST_R[jt]));
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
									
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=Kcon;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  				
						
						
						
						  nbco1++;		  
					  }	  
					  	  
		  
	  }*/
      }	   

nbco=nbco1;

/*
// Nouveau contacts
                      
	for(int ii=0;ii<vecsize;ii++){
		
			for(int lt=0;lt<nocoul[ii];lt++){
			int it = coul[ii][lt];

				for(int kt=0;kt<nonumc[ii];kt++){

					 for(int kjt=0;kjt<nocoul[numc[ii][kt]];kjt++){ 
					 
					  int jt=coul[numc[ii][kt]][kjt];
		
						
					  if((kt>0)||((kt==0)&&(jt>it))){	      
					  dx=LIST_X[it]-LIST_X[jt]; 
					  dy=LIST_Y[it]-LIST_Y[jt]; 
					  dz=LIST_Z[it]-LIST_Z[jt]; 					  
					  nd2=sqrt(dx*dx+dy*dy+dz*dz);
					  r2=(LIST_R[it]+LIST_R[jt]);
					
/////////////////////////////////////////////////////////////////////////	//Contact
					  if(nd2<r2*(1.+epsi))
					  {
                            
							   bool booloc=0 ; 
								  
							   if(NBCONTCO[it]>0){
								 int ite=0;
								 while((ite<NBCONTCO[it])&&(!booloc)){
								   unsigned int numco=NOCONT[it][ite];
								 
								   if((numco<nbco1)&&(((CONT[numco][0]==it)&&(CONT[numco][1]==jt))||((CONT[numco][0]==jt)&&(CONT[numco][1]==it)))) booloc=1;			   	
								
								  ite++;						   
								 } 
							   }	
							   
							   if(!booloc){  
                            
							   TYPCO[nbco]=0;
						       nbcoe++;
							   CONT[nbco][0]=it;
							   CONT[nbco][1]=jt;                       
							   
							   DCONTX[nbco]=dx;
							   DCONTY[nbco]=dy;
							   DCONTZ[nbco]=dz;	   
							   DCONT[nbco] =nd2;
							   DCONTO[nbco1]=r2;
							   
		                         // Normale
							   
							   NCONT[nbco][0]=dx/nd2;
							   NCONT[nbco][1]=dy/nd2; 
							   NCONT[nbco][2]=dz/nd2;       
							   
							   // Zero numerique
							   
							   if(fabs(NCONT[nbco][0])<1e-18) NCONT[nbco][0]=0.;
							   if(fabs(NCONT[nbco][1])<1e-18) NCONT[nbco][1]=0.;
							   if(fabs(NCONT[nbco][2])<1e-18) NCONT[nbco][2]=0.;
							   
							   // Tangente s
							   
							   if((NCONT[nbco][0]*NCONT[nbco][1]==0.)&&(NCONT[nbco][1]*NCONT[nbco][2]==0.)&&(NCONT[nbco][2]*NCONT[nbco][0]==0.)){
								  
								  if(NCONT[nbco][0]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 1.;						
								  }else if(NCONT[nbco][1]!=0.){
								   NCONT[nbco][3] = 1.;
								   NCONT[nbco][4] = 0.;
								   NCONT[nbco][5] = 0.;									  
								  }else if(NCONT[nbco][2]!=0.){
								   NCONT[nbco][3] = 0.;
								   NCONT[nbco][4] = 1.;
								   NCONT[nbco][5] = 0.;									  
								  }
								   
							   }else{
									
									if(NCONT[nbco][0]==0.){
									
									fors  = -(NCONT[nbco][1])/NCONT[nbco][2];
									norms = sqrt(fors*fors+2.);

									NCONT[nbco][3] = 1./norms;
									NCONT[nbco][4] = 1./norms;
									NCONT[nbco][5] = fors/norms;					
									
									}else{
										
								   fors  = -(NCONT[nbco][1]+NCONT[nbco][2])/NCONT[nbco][0];
								   norms = sqrt(fors*fors+2.);

								   NCONT[nbco][3] = fors/norms;
								   NCONT[nbco][4] = 1./norms;
								   NCONT[nbco][5] = 1./norms;							
									
									}
															
									
							   }
							   
							   
							   // Tangente t				   
							   
								NCONT[nbco][6] = NCONT[nbco][1]*NCONT[nbco][5] - NCONT[nbco][2]*NCONT[nbco][4];
								NCONT[nbco][7] = NCONT[nbco][2]*NCONT[nbco][3] - NCONT[nbco][0]*NCONT[nbco][5];
								NCONT[nbco][8] = NCONT[nbco][0]*NCONT[nbco][4] - NCONT[nbco][1]*NCONT[nbco][3];    
								
															   
								NOCONT[it][NBCONTCO[it]]=nbco;
								NOCONT[jt][NBCONTCO[jt]]=nbco;
								NBCONTCO[it]++;
								NBCONTCO[jt]++; 
			//					if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
				//				if(NBCONTCO[it]>numax) numax=NBCONTCO[it];	
									   
								 rij=(LIST_R[it]*LIST_R[jt])/(LIST_R[it]+LIST_R[jt]);								
								 hij=fabs(DCONT[nbco]-(LIST_R[it]+LIST_R[jt]));
								 
								 
									if((LIST_P[it]==0)&&(LIST_P[jt]==0)){				  
									iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
									}
									if(((LIST_P[it]==1)&&(LIST_P[jt]==0))||((LIST_P[it]==0)&&(LIST_P[jt]==1))) {	
									iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
									}
									if((LIST_P[it]==1)&&(LIST_P[jt]==1)){				  
									iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
									}
								 
								 Eij=1./iEij;	
								// Kn=4.*Eij*rij;
								 Kn=4.*Eij*sqrt(rij*hij)/3.;	
								 	     Kn=Kn/100.;	
									     Kn=Kcon;	 							
								 Kt=Kn;
								 
								dt=min(dt,sqrt(LIST_M[it]/Kn));
								dt=min(dt,sqrt(LIST_M[jt]/Kn));
								dt=min(dt,sqrt(LIST_I[it]/(Kt*LIST_R[it]*LIST_R[it])));
								dt=min(dt,sqrt(LIST_I[jt]/(Kt*LIST_R[jt]*LIST_R[jt])));  		
															   
							   nbco++;						   
							   }														   
						
						  
					  }
/////////////////////////////////////////////////////////////////////////	
				
					  } 

			     }

			}	

		}
		
	}

*/

/////////////////////////////////////////////////////////////////////////
if(amort!=0.){
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

//cout<<"Pas de temps:"<<dt<<endl;
//cout<<"nbco:"<<nbco<<endl;
}






void selco_corci0(bool * TYPCO, R & dt, int & nbco, int NB_SPH, int ** CONT, R * LIST_R, R * LIST_M, R * LIST_I, R * DCONT, R ** NCONT, R ** VALCOH, R ** VALAMO, unsigned int ** NOCONT, int * NBCONTCO, R * DCONTX, R * DCONTY, R * DCONTZ, R Cs, int NMAXCONT, int NMAXZ)
{

// Mise à blanc des contacts par corci
  for(int it=0;it<NB_SPH;it++){
     NBCONTCO[it]=0;
  }

dt=1e9;
unsigned int nbco1=0;
int ii,it,jt;
int numax=0;

  for(ii=0;ii<nbco;ii++){
	  
	  if(TYPCO[ii]){
	
		   TYPCO[nbco1]=1;  	
		   
                   CONT[nbco1][0]=CONT[ii][0];
		   CONT[nbco1][1]=CONT[ii][1];                       
		   
		   it=CONT[nbco1][0];
		   jt=CONT[nbco1][1];     
		    
   		   DCONTX[nbco1]=DCONTX[ii];
		   DCONTY[nbco1]=DCONTY[ii];
		   DCONTZ[nbco1]=DCONTZ[ii];		   
		   DCONT[nbco1] =DCONT[ii];		 

		   NCONT[nbco1][0]=NCONT[ii][0];
		   NCONT[nbco1][1]=NCONT[ii][1]; 
		   NCONT[nbco1][2]=NCONT[ii][2];
		   NCONT[nbco1][3]=NCONT[ii][3];
		   NCONT[nbco1][4]=NCONT[ii][4]; 
		   NCONT[nbco1][5]=NCONT[ii][5];     
		   NCONT[nbco1][6]=NCONT[ii][6];
		   NCONT[nbco1][7]=NCONT[ii][7]; 
		   NCONT[nbco1][8]=NCONT[ii][8];     		                     		    
                       		    
		   NOCONT[it][NBCONTCO[it]]=nbco1;
		   NOCONT[jt][NBCONTCO[jt]]=nbco1;
		   NBCONTCO[it]++;
		   NBCONTCO[jt]++; 		                       		    
		if(NBCONTCO[jt]>numax) numax=NBCONTCO[jt];
		if(NBCONTCO[it]>numax) numax=NBCONTCO[it];

			 VALCOH[nbco1][0]=VALCOH[ii][0];
			 VALCOH[nbco1][1]=VALCOH[ii][1];
			 VALCOH[nbco1][2]=VALCOH[ii][2];	
			 VALCOH[nbco1][3]=VALCOH[ii][3];	
			 VALCOH[nbco1][4]=VALCOH[ii][4];
			 VALCOH[nbco1][5]=VALCOH[ii][5];
			 VALCOH[nbco1][6]=VALCOH[ii][6];	
			
			 VALAMO[nbco1][0]=VALAMO[ii][0];
			 VALAMO[nbco1][1]=VALAMO[ii][1];
			 VALAMO[nbco1][2]=VALAMO[ii][2];
			 VALAMO[nbco1][3]=VALAMO[ii][3];
			 VALAMO[nbco1][4]=VALAMO[ii][4];                      
			 VALAMO[nbco1][5]=VALAMO[ii][5];
			 VALAMO[nbco1][6]=VALAMO[ii][6];	                       
		   
			 dt=min(dt,sqrt(LIST_M[it]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_M[jt]/VALCOH[nbco1][0]));
			 dt=min(dt,sqrt(LIST_I[it]/(VALCOH[nbco1][1]*LIST_R[it]*LIST_R[it])));
			 dt=min(dt,sqrt(LIST_I[jt]/(VALCOH[nbco1][1]*LIST_R[jt]*LIST_R[jt])));  

		  
		  nbco1++;
	  }
	  
  }	  

nbco=nbco1;

if(numax>NMAXZ){
cout<<"Nombre maximal de contacts par sph atteint: "<<numax<<"/"<<NMAXZ<<" - stop "<<endl;
exit(0);
}

if(nbco>NMAXCONT){
cout<<"Nombre maximal de contacts sph-sph atteint:"<<nbco<<"/"<<NMAXCONT<<" - stop "<<endl;
exit(0);
}


dt=Cs*dt;
//cout<<"Pas de temps:"<<dt<<endl;
}





