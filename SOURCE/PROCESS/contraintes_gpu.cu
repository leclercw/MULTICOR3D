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

#include "contraintes_gpu.h"

__global__
void slocal(int NB_SPH,R * LIST_IND, int * NBCONTCO, unsigned int * NOCONT, R * LIST_R, R * NCONT, R * FCJI, R * PSIG11, R * PSIG12, R * PSIG13, R * PSIG22, R * PSIG23, R * PSIG33, int NMAXZ) {

	int numcor = (blockIdx.x*blockDim.x)+threadIdx.x;

	if (numcor<NB_SPH) {

		R sig11=0.;
		R sig12=0.;
		R sig13=0.;
		R sig21=0.;
		R sig22=0.;
		R sig23=0.;
		R sig31=0.;
		R sig32=0.;
		R sig33=0.;  			   
		  
		R coefij=LIST_IND[numcor];

					 for(int kt=0;kt<NBCONTCO[numcor];kt++){ 

				                unsigned int numc=NOCONT[numcor*NMAXZ+kt]; 
//  if(numcor==0){ printf("NUMC: (%i,%i)\n",kt,numc);}
						R ray=LIST_R[numcor];
						R n1=-NCONT[numc*9+0];
						R n2=-NCONT[numc*9+1];
						R n3=-NCONT[numc*9+2];
/*
						sig11=numc;
						sig12=numc;
						sig13=numc;
						sig21=numc;
						sig22=numc;
						sig23=numc;
						sig31=numc;
						sig32=numc;
						sig33=numc;*/

						sig11=sig11+coefij*ray*n1*FCJI[numc*3+0];
						sig12=sig12+coefij*ray*n1*FCJI[numc*3+1];
						sig13=sig13+coefij*ray*n1*FCJI[numc*3+2];
						sig21=sig21+coefij*ray*n2*FCJI[numc*3+0];
						sig22=sig22+coefij*ray*n2*FCJI[numc*3+1];
						sig23=sig23+coefij*ray*n2*FCJI[numc*3+2];
						sig31=sig31+coefij*ray*n3*FCJI[numc*3+0];
						sig32=sig32+coefij*ray*n3*FCJI[numc*3+1];
						sig33=sig33+coefij*ray*n3*FCJI[numc*3+2];

					 }
	

				sig12=(sig12+sig21)/2.;
				sig13=(sig13+sig31)/2.;
				sig23=(sig23+sig32)/2.;                
						 
                               
				PSIG11[numcor]=sig11;
				PSIG12[numcor]=sig12;
				PSIG13[numcor]=sig13;                              
				PSIG22[numcor]=sig22;
				PSIG23[numcor]=sig23;                                      
				PSIG33[numcor]=sig33; 

	  	}
}

__global__
void shalo(int NB_SPH, int * NBHALO, unsigned int * NOHALO, R * VOLHALO, R * PSIG11, R * PSIG12, R * PSIG13, R * PSIG22, R * PSIG23, R * PSIG33, R * SIG11, R * SIG12, R * SIG13, R * SIG22, R * SIG23, R * SIG33, R * VONMIS, R * TRACE, R * SIG1, R *SIG2, R * SIG3, R minsig11, R maxsig11, R minsig12, R maxsig12, R minsig13, R maxsig13, R minsig22, R maxsig22, R minsig23, R maxsig23, R minsig33, R maxsig33, R minvm, R maxvm, R mintrac, R maxtrac, R minsig1, R maxsig1, R minsig2, R maxsig2, R minsig3, R maxsig3, int NMAXHALO) {
	
	int jt = (blockIdx.x*blockDim.x)+threadIdx.x;

	if (jt<NB_SPH) {

		R sig11=PSIG11[jt];
		R sig12=PSIG12[jt];
		R sig13=PSIG13[jt];
		R sig22=PSIG22[jt];
		R sig23=PSIG23[jt];
		R sig33=PSIG33[jt];


         	  for(int kt=1;kt<NBHALO[jt];kt++){ 
			//printf("NBHALOav: %i\n",long(700000)*NMAXHALO+kt);
			unsigned int lt=NOHALO[long(jt)*NMAXHALO+kt];			
			//printf("NBHALOap: %i\n",lt);
			sig11=sig11+PSIG11[lt];
			sig12=sig12+PSIG12[lt];
			sig13=sig13+PSIG13[lt];                             
			sig22=sig22+PSIG22[lt];
			sig23=sig23+PSIG23[lt];                                      
			sig33=sig33+PSIG33[lt];  
		  }


				sig11/=VOLHALO[jt];
				sig12/=VOLHALO[jt];
				sig13/=VOLHALO[jt];
				sig22/=VOLHALO[jt];
				sig23/=VOLHALO[jt];
				sig33/=VOLHALO[jt];	

				SIG11[jt]=sig11;
				SIG12[jt]=sig12;
				SIG13[jt]=sig13;
				SIG22[jt]=sig22;
				SIG23[jt]=sig23;
				SIG33[jt]=sig33;

				R i1=sig11+sig22+sig33;
				R i2=sig11*sig22+sig22*sig33+sig33*sig11-sig12*sig12-sig23*sig23-sig13*sig13;
				R i3=sig11*(sig22*sig33-sig23*sig23)-sig12*(sig12*sig33-sig13*sig23)+sig13*(sig12*sig23-sig22*sig13);
					
				R b=-i1;
				R c=i2;
				R d=-i3;
				R p=c-b*b/3.;
				R q=d-b*c/3.+2*b*b*b/27.;
				R detd=4*c*c*c+27*d*d+4*d*b*b*b-b*b*c*c-18*b*c*d; 


				R s1=0.;
				R s2=0.;
				R s3=0.;

				if (fabs(detd)<=1e-30){
				double t=-q/2.;
				s1=2.*pow(t,1./3)-b/3.;
				s2=-pow(t,1./3)-b/3.;
				s3=s2;  			
				}
				else{
				R r=sqrt(-p*p*p/27.);
				R theta=acos(-q/(2*r));
				s1=2.*sqrt(-p/3.)*cos(theta/3.)-b/3.;
				s2=2.*sqrt(-p/3.)*cos((theta+2.*3.14159265358979323846)/3.)-b/3.;
				s3=2.*sqrt(-p/3.)*cos((theta+4.*3.14159265358979323846)/3.)-b/3.;

				}
			    
			    
				R smax=max(max(s1,s2),s3);
				R smin=min(min(s1,s2),s3);
				if(s1==smax){s2=max(s2,s3);}
				else if(s2==smax){s2=max(s1,s3);}
				else if(s3==smax){s2=max(s1,s2);}              

				s1=smax;
				s3=smin;
				R trac=s1+s2+s3;                
				R vmis=sqrt((sig11-sig22)*(sig11-sig22)+(sig33-sig22)*(sig33-sig22)+(sig11-sig33)*(sig11-sig33)+6.*(sig12*sig12+sig13*sig13+sig23*sig23))/sqrt(2.);
				if(vmis!=vmis) {vmis=0.;}    

				VONMIS[jt]=vmis;
				TRACE[jt]=trac;
				SIG1[jt]=s1;
				SIG2[jt]=s2;
				SIG3[jt]=s3;
/*
				minvm=fmin(vmis,minvm);
				maxvm=fmax(vmis,maxvm);
				
				mintrac=fmin(trac,mintrac);
				maxtrac=fmax(trac,maxtrac);			
				
				minsig11=fmin(sig11,minsig11);
				maxsig11=fmax(sig11,maxsig11);
				minsig12=fmin(sig12,minsig12);
				maxsig12=fmax(sig12,maxsig12);
				minsig13=fmin(sig13,minsig13);
				maxsig13=fmax(sig13,maxsig13);				
				minsig22=fmin(sig22,minsig22);
				maxsig22=fmax(sig22,maxsig22);					
				minsig23=fmin(sig23,minsig23);
				maxsig23=fmax(sig23,maxsig23);	
				minsig33=fmin(sig33,minsig33);
				maxsig33=fmax(sig33,maxsig33);					
								
				minsig1=fmin(s1,minsig1);
				maxsig1=fmax(s1,maxsig1);
				minsig2=fmin(s2,minsig2);
				maxsig2=fmax(s2,maxsig2);	
				minsig3=fmin(s3,minsig3);
				maxsig3=fmax(s3,maxsig3);	*/
		}				
	}



void contrainteshalo_gpu(R Pi, R coef1, int ite,int NBENREG, int NB_SPH, int NBCO, int NMAXZ,int NMAXHALO, int NMAXCONT, R H_TOT, R V_TOT, R Z_TOT,R * LIST_R, R ** FCJI, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * VONMIS, R * TRACE, R * SIG11, R * SIG12, R * SIG13, R * SIG22, R * SIG23, R * SIG33, R * SIG1, R * SIG2, R * SIG3,R &minvm, R &maxvm,R &mintrac, R &maxtrac, R &minsig11, R &maxsig11, R &minsig12, R &maxsig12, R &minsig13, R &maxsig13, R &minsig22, R &maxsig22, R &minsig23, R &maxsig23, R &minsig33, R &maxsig33, R &minsig1, R &maxsig1, R &minsig2, R &maxsig2, R &minsig3, R &maxsig3, bool * EDGE, unsigned int ** NOHALO, int * NBHALO, R * LIST_V, R * VOLHALO,R * LIST_IND) {

int it;

mintrac = 1e12;
maxtrac =  -1e12;
minvm = 1e12;
maxvm =  -1e12;
minsig11 = 1e12;
maxsig11 = -1e12;
minsig12 = 1e12;
maxsig12 = -1e12;
minsig13 = 1e12;
maxsig13 = -1e12;
minsig22 = 1e12;
maxsig22 = -1e12;
minsig23 = 1e12;
maxsig23 = -1e12;
minsig33 = 1e12;
maxsig33 = -1e12;
minsig1 = 1e12;
maxsig1 = -1e12;
minsig2 = 1e12;
maxsig2 = -1e12;
minsig3 = 1e12;
maxsig3 = -1e12;

/////////////////////////////////
// Vecteurs/matrices device
R * dLIST_IND;
R * dLIST_R;
R * dPSIG11;
R * dPSIG12;
R * dPSIG13;
R * dPSIG22;
R * dPSIG23;
R * dPSIG33;
R * dSIG11;
R * dSIG12;
R * dSIG13;
R * dSIG22;
R * dSIG23;
R * dSIG33;
R * dVONMIS;
R * dTRACE;
R * dSIG1;
R * dSIG2;
R * dSIG3;
R * dVOLHALO;
int * dNBHALO;
int * dNBCONTCO;
unsigned int * dNOCONT;
R * dNCONT;
R * dFCJI;
unsigned int * dNOHALO;

/////////////////////////////////
// Allocation mémoire

cudaMalloc((void **)&dLIST_IND, NB_SPH*sizeof(R));
cudaMemcpy(dLIST_IND, LIST_IND, NB_SPH*sizeof(R), cudaMemcpyHostToDevice);

cudaMalloc((void **)&dLIST_R, NB_SPH*sizeof(R));
cudaMemcpy(dLIST_R, LIST_R, NB_SPH*sizeof(R), cudaMemcpyHostToDevice);

cudaMalloc((void **)&dPSIG11, NB_SPH*sizeof(R));
cudaMalloc((void **)&dPSIG12, NB_SPH*sizeof(R));
cudaMalloc((void **)&dPSIG13, NB_SPH*sizeof(R));
cudaMalloc((void **)&dPSIG22, NB_SPH*sizeof(R));
cudaMalloc((void **)&dPSIG23, NB_SPH*sizeof(R));
cudaMalloc((void **)&dPSIG33, NB_SPH*sizeof(R));

cudaMalloc((void **)&dNBCONTCO, NB_SPH*sizeof(int));
cudaMemcpy(dNBCONTCO, NBCONTCO, NB_SPH*sizeof(int), cudaMemcpyHostToDevice);

cudaMalloc((void **)&dNOCONT, NB_SPH*NMAXZ*sizeof(unsigned int));
cudaMemcpy(dNOCONT, NOCONT[0], NB_SPH*NMAXZ*sizeof(unsigned int), cudaMemcpyHostToDevice);

cudaMalloc((void **)&dNCONT, NMAXCONT*9*sizeof(R));
cudaMemcpy(dNCONT, NCONT[0], NMAXCONT*9*sizeof(R), cudaMemcpyHostToDevice);

cudaMalloc((void **)&dFCJI, NMAXCONT*3*sizeof(R));
cudaMemcpy(dFCJI, FCJI[0], NMAXCONT*3*sizeof(R), cudaMemcpyHostToDevice);

/////////////////////////////////
// Contraintes à l'échelle de la particule

dim3 DimGrid ((NB_SPH-1)/256+1,1,1) ;
dim3 DimBlock (256,1,1) ;
slocal<<<DimGrid, DimBlock>>>(NB_SPH,dLIST_IND,dNBCONTCO,dNOCONT,dLIST_R,dNCONT,dFCJI,dPSIG11,dPSIG12, dPSIG13, dPSIG22,dPSIG23,dPSIG33,NMAXZ);
	
/////////////////////////////////
// Libération mémoire
	cudaFree(dLIST_IND);
	cudaFree(dLIST_R);
	cudaFree(dNBCONTCO);
	cudaFree(dNOCONT);
	cudaFree(dNCONT);
	cudaFree(dFCJI);

/////////////////////////////////
// Allocation mémoire

cudaMalloc((void **)&dSIG11, NB_SPH*sizeof(R));
cudaMalloc((void **)&dSIG12, NB_SPH*sizeof(R));
cudaMalloc((void **)&dSIG13, NB_SPH*sizeof(R));
cudaMalloc((void **)&dSIG22, NB_SPH*sizeof(R));
cudaMalloc((void **)&dSIG23, NB_SPH*sizeof(R));
cudaMalloc((void **)&dSIG33, NB_SPH*sizeof(R));
cudaMalloc((void **)&dVONMIS, NB_SPH*sizeof(R));
cudaMalloc((void **)&dTRACE, NB_SPH*sizeof(R));
cudaMalloc((void **)&dSIG1, NB_SPH*sizeof(R));
cudaMalloc((void **)&dSIG2, NB_SPH*sizeof(R));
cudaMalloc((void **)&dSIG3, NB_SPH*sizeof(R));

cudaMalloc((void **)&dVOLHALO, NB_SPH*sizeof(R));
cudaMemcpy(dVOLHALO, VOLHALO, NB_SPH*sizeof(R), cudaMemcpyHostToDevice);

cudaMalloc((void **)&dNBHALO, NB_SPH*sizeof(int));
cudaMemcpy(dNBHALO, NBHALO, NB_SPH*sizeof(int), cudaMemcpyHostToDevice);

cudaMalloc((void **)&dNOHALO, long(NB_SPH)*NMAXHALO*sizeof(unsigned int));
cudaMemcpy(dNOHALO, NOHALO[0], long(NB_SPH)*NMAXHALO*sizeof(unsigned int), cudaMemcpyHostToDevice);

size_t free, total;

printf("\n");

cudaMemGetInfo(&free,&total);

printf("%d KB free of total %d KB\n",free/1024,total/1024);

/////////////////////////////////
// Contraintes à l'échelle du halo
shalo<<<DimGrid, DimBlock>>>(NB_SPH, dNBHALO,dNOHALO,dVOLHALO,dPSIG11,dPSIG12,dPSIG13,dPSIG22,dPSIG23,dPSIG33,dSIG11,dSIG12,dSIG13,dSIG22,dSIG23,dSIG33,dVONMIS,dTRACE,dSIG1,dSIG2,dSIG3,minsig11, maxsig11,minsig12,maxsig12,minsig13,maxsig13,minsig22,maxsig22,minsig23,maxsig23,minsig33,maxsig33,minvm,maxvm,mintrac,maxtrac,minsig1,maxsig1,minsig2,maxsig2,minsig3,maxsig3,NMAXHALO);

/////////////////////////////////
// Copies des vecteurs/matrices utiles vers l'host
    cudaMemcpy(SIG11, dSIG11, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(SIG12, dSIG12, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(SIG13, dSIG13, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(SIG22, dSIG22, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(SIG23, dSIG23, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(SIG33, dSIG33, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(VONMIS, dVONMIS, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(TRACE, dTRACE, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);	
    cudaMemcpy(SIG1, dSIG1, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(SIG2, dSIG2, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);
    cudaMemcpy(SIG3, dSIG3, NB_SPH*sizeof(R), cudaMemcpyDeviceToHost);

/////////////////////////////////
// Libération mémoire

	cudaFree(dPSIG11);
	cudaFree(dPSIG12);
	cudaFree(dPSIG13);
	cudaFree(dPSIG22);
	cudaFree(dPSIG23);
	cudaFree(dPSIG33);
	cudaFree(dSIG11);
	cudaFree(dSIG12);
	cudaFree(dSIG13);
	cudaFree(dSIG22);
	cudaFree(dSIG23);
	cudaFree(dSIG33);
	cudaFree(dVONMIS);
	cudaFree(dTRACE);
	cudaFree(dSIG1);
	cudaFree(dSIG2);
	cudaFree(dSIG3);
	cudaFree(dNBHALO);
	cudaFree(dVOLHALO);
	cudaFree(dNOHALO);

/////////////////////////////////
// Traitement Min/max


for(it=0;it<NB_SPH;it++){
				minvm=min(VONMIS[it],minvm);
				maxvm=max(VONMIS[it],maxvm);
				
				mintrac=min(TRACE[it],mintrac);
				maxtrac=max(TRACE[it],maxtrac);			
				
				minsig11=min(SIG11[it],minsig11);
				maxsig11=max(SIG11[it],maxsig11);
				minsig12=min(SIG12[it],minsig12);
				maxsig12=max(SIG12[it],maxsig12);
				minsig13=min(SIG13[it],minsig13);
				maxsig13=max(SIG13[it],maxsig13);				
				minsig22=min(SIG22[it],minsig22);
				maxsig22=max(SIG22[it],maxsig22);					
				minsig23=min(SIG23[it],minsig23);
				maxsig23=max(SIG23[it],maxsig23);	
				minsig33=min(SIG33[it],minsig33);
				maxsig33=max(SIG33[it],maxsig33);						
								
				minsig1=min(SIG1[it],minsig1);
				maxsig1=max(SIG1[it],maxsig1);
				minsig2=min(SIG2[it],minsig2);
				maxsig2=max(SIG2[it],maxsig2);	
				minsig3=min(SIG3[it],minsig3);
				maxsig3=max(SIG3[it],maxsig3);	

}	

if(ite%NBENREG==0){ 
cout<<"Maxsig11:"<<maxsig11<<endl;
cout<<"Minsig11:"<<minsig11<<endl;
cout<<"Maxsig12:"<<maxsig12<<endl;
cout<<"Minsig12:"<<minsig12<<endl;
cout<<"Maxsig13:"<<maxsig13<<endl;
cout<<"Minsig13:"<<minsig13<<endl;
cout<<"Maxsig22:"<<maxsig22<<endl;
cout<<"Minsig22:"<<minsig22<<endl;
cout<<"Maxsig23:"<<maxsig23<<endl;
cout<<"Minsig23:"<<minsig23<<endl;
cout<<"Maxsig33:"<<maxsig33<<endl;
cout<<"Minsig33:"<<minsig33<<endl;
cout<<"Maxsig1:"<<maxsig1<<endl;
cout<<"Minsig1:"<<minsig1<<endl;
cout<<"Maxsig2:"<<maxsig2<<endl;
cout<<"Minsig2:"<<minsig2<<endl;
cout<<"Maxsig3:"<<maxsig3<<endl;
cout<<"Minsig3:"<<minsig3<<endl;
cout<<"Maxvm:"<<maxvm<<endl;
cout<<"Minvm:"<<minvm<<endl;
cout<<"Maxtrac:"<<maxtrac<<endl;
cout<<"Mintrac:"<<mintrac<<endl;
}



}
