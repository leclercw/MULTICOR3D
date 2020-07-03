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

#include "contraintes.h"

void contrainteshalo(R Pi, R coef1, int ite,int NBENREG, int NB_SPH, R H_TOT, R V_TOT, R Z_TOT,R * LIST_R, R ** FCJI, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * VONMIS, R * TRACE, R * SIG11, R * SIG12, R * SIG13, R * SIG22, R * SIG23, R * SIG33, R * SIG1, R * SIG2, R * SIG3,R &minvm, R &maxvm,R &mintrac, R &maxtrac, R &minsig11, R &maxsig11, R &minsig12, R &maxsig12, R &minsig13, R &maxsig13, R &minsig22, R &maxsig22, R &minsig23, R &maxsig23, R &minsig33, R &maxsig33, R &minsig1, R &maxsig1, R &minsig2, R &maxsig2, R &minsig3, R &maxsig3, bool * EDGE, unsigned int ** NOHALO, int * NBHALO, R * LIST_V, R * VOLHALO,R * LIST_IND){

unsigned int numc;
R ray;
R vmis;
R sig11,sig12,sig13,sig21,sig22,sig23,sig31,sig32,sig33;
R sig1,sig2,sig3,trac;
R i1,i2,i3;
R detd,b,c,d,t,p,q,r,theta;
R s1,s2,s3,smin,smax;
R n1,n2,n3;
int jt,kt,numcor;
unsigned it;

R coefij;

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

R PSIG11[NB_SPH];
R PSIG12[NB_SPH];
R PSIG13[NB_SPH];
R PSIG22[NB_SPH];
R PSIG23[NB_SPH];
R PSIG33[NB_SPH];

int num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(coefij,numc,ray,n1,n2,n3,sig11,sig12,sig13,sig21,sig22,sig23,sig31,sig32,sig33,kt,numcor)

	for(numcor=0;numcor<NB_SPH;numcor++){
	
		 sig11=0.;
		 sig12=0.;
		 sig13=0.;
		 sig21=0.;
		 sig22=0.;
		 sig23=0.;
		 sig31=0.;
		 sig32=0.;
		 sig33=0.;  			   
		  
		   coefij=LIST_IND[numcor];

					 for(kt=0;kt<NBCONTCO[numcor];kt++){ 
					 
						numc=NOCONT[numcor][kt];
						ray=LIST_R[numcor];
						n1=-NCONT[numc][0];
						n2=-NCONT[numc][1];
						n3=-NCONT[numc][2];
		
						sig11=sig11+coefij*ray*n1*FCJI[numc][0];
						sig12=sig12+coefij*ray*n1*FCJI[numc][1];
						sig13=sig13+coefij*ray*n1*FCJI[numc][2];
						sig21=sig21+coefij*ray*n2*FCJI[numc][0];
						sig22=sig22+coefij*ray*n2*FCJI[numc][1];
						sig23=sig23+coefij*ray*n2*FCJI[numc][2];
						sig31=sig31+coefij*ray*n3*FCJI[numc][0];
						sig32=sig32+coefij*ray*n3*FCJI[numc][1];
						sig33=sig33+coefij*ray*n3*FCJI[numc][2]; 

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


R sumsig=0.;

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(vmis,sig11,sig12,sig13,sig22,sig23,sig33,sig1,sig2,sig3,trac,i1,i2,i3,detd,b,c,d,t,p,q,r,theta,s1,s2,s3,smin,smax,it,jt,kt)  reduction(min:minvm,mintrac,minsig11,minsig12,minsig13,minsig22,minsig23,minsig33,minsig1,minsig2,minsig3) reduction(max:maxvm,maxtrac,maxsig11,maxsig12,maxsig13,maxsig22,maxsig23,maxsig33,maxsig1,maxsig2,maxsig3) reduction(+:sumsig)
	for(jt=0;jt<NB_SPH;jt++){

		sig11=PSIG11[jt];
		sig12=PSIG12[jt];
		sig13=PSIG13[jt];
		sig22=PSIG22[jt];
		sig23=PSIG23[jt];
		sig33=PSIG33[jt];


         	  for(kt=1;kt<NBHALO[jt];kt++){ 
			it=NOHALO[jt][kt];
			sig11+=PSIG11[it];
			sig12+=PSIG12[it];
			sig13+=PSIG13[it];                             
			sig22+=PSIG22[it];
			sig23+=PSIG23[it];                                      
			sig33+=PSIG33[it];  
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

				i1=sig11+sig22+sig33;
				i2=sig11*sig22+sig22*sig33+sig33*sig11-sig12*sig12-sig23*sig23-sig13*sig13;
				i3=sig11*(sig22*sig33-sig23*sig23)-sig12*(sig12*sig33-sig13*sig23)+sig13*(sig12*sig23-sig22*sig13);
					
				b=-i1;
				c=i2;
				d=-i3;
				p=c-b*b/3.;
				q=d-b*c/3.+2*b*b*b/27.;
				detd=4*c*c*c+27*d*d+4*d*b*b*b-b*b*c*c-18*b*c*d; 


				s1=0.;
				s2=0.;
				s3=0.;

				if (fabs(detd)<=1e-30){
				t=-q/2.;
				s1=2.*pow(t,1./3)-b/3.;
				s2=-pow(t,1./3)-b/3.;
				s3=s2;  			
				}
				else{
				r=sqrt(-p*p*p/27.);
				theta=acos(-q/(2*r));
				s1=2.*sqrt(-p/3.)*cos(theta/3.)-b/3.;
				s2=2.*sqrt(-p/3.)*cos((theta+2.*3.14159265358979323846)/3.)-b/3.;
				s3=2.*sqrt(-p/3.)*cos((theta+4.*3.14159265358979323846)/3.)-b/3.;

				}
			    
			    
				smax=max(max(s1,s2),s3);
				smin=min(min(s1,s2),s3);
				if(s1==smax){s2=max(s2,s3);}
				else if(s2==smax){s2=max(s1,s3);}
				else if(s3==smax){s2=max(s1,s2);}              

				s1=smax;
				s3=smin;
				trac=s1+s2+s3;                
				vmis=sqrt((sig11-sig22)*(sig11-sig22)+(sig33-sig22)*(sig33-sig22)+(sig11-sig33)*(sig11-sig33)+6.*(sig12*sig12+sig13*sig13+sig23*sig23))/sqrt(2.);
				if(vmis!=vmis) {vmis=0.;}    

				VONMIS[jt]=vmis;
				TRACE[jt]=trac;
				SIG1[jt]=s1;
				SIG2[jt]=s2;
				SIG3[jt]=s3;

				minvm=min(vmis,minvm);
				maxvm=max(vmis,maxvm);
				
				mintrac=min(trac,mintrac);
				maxtrac=max(trac,maxtrac);			
				
				minsig11=min(sig11,minsig11);
				maxsig11=max(sig11,maxsig11);
				minsig12=min(sig12,minsig12);
				maxsig12=max(sig12,maxsig12);
				minsig13=min(sig13,minsig13);
				maxsig13=max(sig13,maxsig13);				
				minsig22=min(sig22,minsig22);
				maxsig22=max(sig22,maxsig22);							
				minsig23=min(sig23,minsig23);
				maxsig23=max(sig23,maxsig23);	
				minsig33=min(sig33,minsig33);
				maxsig33=max(sig33,maxsig33);						
								
				minsig1=min(s1,minsig1);
				maxsig1=max(s1,maxsig1);
				minsig2=min(s2,minsig2);
				maxsig2=max(s2,maxsig2);	
				minsig3=min(s3,minsig3);
				maxsig3=max(s3,maxsig3);	
	
sumsig+=(s1+s2+s3);
					
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
cout<<"sumsig:"<<sumsig<<endl;
}

}



void contraintesloc(R Pi, R coef1, int NBENREG, R H_TOT, R V_TOT, R Z_TOT, R * LIST_R, R ** FCJI, R ** NCONT,int ite, int NB_SPH, R * VONMIS, R * TRACE, R * SIG11, R * SIG12, R * SIG13, R * SIG22, R * SIG23, R * SIG33, R * SIG1, R * SIG2, R * SIG3,R &minvm, R &maxvm,R &mintrac, R &maxtrac, R &minsig11, R &maxsig11, R &minsig12, R &maxsig12, R &minsig13, R &maxsig13, R &minsig22, R &maxsig22, R &minsig23, R &maxsig23, R &minsig33, R &maxsig33, R &minsig1, R &maxsig1, R &minsig2, R &maxsig2, R &minsig3, R &maxsig3, unsigned int ** NOCONT, int * NBCONTCO,bool * EDGE, R * LIST_V, R * LIST_IND){

unsigned int numc;
int numcor;
R ray;
R vmis;
R sig11,sig12,sig13,sig21,sig22,sig23,sig31,sig32,sig33;
R sig1,sig2,sig3,trac;
R i1,i2,i3;
R detd,b,c,d,t,p,q,r,theta;
R s1,s2,s3,smin,smax;
R n1,n2,n3;
int kt;

R mmintrac = 1e12;
R mmaxtrac =  -1e12;
R mminvm = 1e12;
R mmaxvm =  -1e12;
R mminsig11 = 1e12;
R mmaxsig11 = -1e12;
R mminsig12 = 1e12;
R mmaxsig12 = -1e12;
R mminsig13 = 1e12;
R mmaxsig13 = -1e12;
R mminsig22 = 1e12;
R mmaxsig22 = -1e12;
R mminsig23 = 1e12;
R mmaxsig23 = -1e12;
R mminsig33 = 1e12;
R mmaxsig33 = -1e12;
R mminsig1 = 1e12;
R mmaxsig1 = -1e12;
R mminsig2 = 1e12;
R mmaxsig2 = -1e12;
R mminsig3 = 1e12;
R mmaxsig3 = -1e12;

int num_dis=0;
R volt=0.;
R volt2=0.;

int num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(numc,numcor,ray,vmis,sig11,sig12,sig13,sig21,sig22,sig23,sig31,sig32,sig33,sig1,sig2,sig3,trac,i1,i2,i3,detd,b,c,d,t,p,q,r,theta,s1,s2,s3,smin,smax,n1,n2,n3,kt) reduction(+:num_dis,volt,volt2) reduction(min:mminvm,mmintrac,mminsig11,mminsig12,mminsig13,mminsig22,mminsig23,mminsig33,mminsig1,mminsig2,mminsig3) reduction(max:mmaxvm,mmaxtrac,mmaxsig11,mmaxsig12,mmaxsig13,mmaxsig22,mmaxsig23,mmaxsig33,mmaxsig1,mmaxsig2,mmaxsig3)
	for(numcor=0;numcor<NB_SPH;numcor++){ 
      
        num_dis++;
  
        sig11=0.;
        sig12=0.;
        sig13=0.;
        sig21=0.;
        sig22=0.;
        sig23=0.;
        sig31=0.;
        sig32=0.;
        sig33=0.;        
   
        volt=volt+coef1*LIST_V[numcor];
      volt2=volt2+LIST_V[numcor];

             for(kt=0;kt<NBCONTCO[numcor];kt++){ 
            
		numc=NOCONT[numcor][kt];
		ray=LIST_R[numcor];
		n1=-NCONT[numc][0];
		n2=-NCONT[numc][1];
		n3=-NCONT[numc][2];

             sig11=sig11+ray*n1*FCJI[numc][0];
             sig12=sig12+ray*n1*FCJI[numc][1];
             sig13=sig13+ray*n1*FCJI[numc][2];
             sig21=sig21+ray*n2*FCJI[numc][0];
             sig22=sig22+ray*n2*FCJI[numc][1];
             sig23=sig23+ray*n2*FCJI[numc][2];
             sig31=sig31+ray*n3*FCJI[numc][0];
             sig32=sig32+ray*n3*FCJI[numc][1];
             sig33=sig33+ray*n3*FCJI[numc][2];  
             }
             
             
			 sig11=LIST_IND[numcor]*sig11/(coef1*LIST_V[numcor]);
			 sig12=LIST_IND[numcor]*sig12/(coef1*LIST_V[numcor]);
			 sig13=LIST_IND[numcor]*sig13/(coef1*LIST_V[numcor]);
			 sig21=LIST_IND[numcor]*sig21/(coef1*LIST_V[numcor]);
			 sig22=LIST_IND[numcor]*sig22/(coef1*LIST_V[numcor]);
			 sig23=LIST_IND[numcor]*sig23/(coef1*LIST_V[numcor]);	
	                 sig31=LIST_IND[numcor]*sig31/(coef1*LIST_V[numcor]);
			 sig32=LIST_IND[numcor]*sig32/(coef1*LIST_V[numcor]);
			 sig33=LIST_IND[numcor]*sig33/(coef1*LIST_V[numcor]);	

				sig12=(sig12+sig21)/2.;
				sig13=(sig13+sig31)/2.;
				sig23=(sig23+sig32)/2.;                  
						 
				i1=sig11+sig22+sig33;
				i2=sig11*sig22+sig22*sig33+sig33*sig11-sig12*sig12-sig23*sig23-sig13*sig13;
				i3=sig11*(sig22*sig33-sig23*sig23)-sig12*(sig12*sig33-sig13*sig23)+sig13*(sig12*sig23-sig22*sig13);
					
				b=-i1;
				c=i2;
				d=-i3;
				p=c-b*b/3.;
				q=d-b*c/3.+2*b*b*b/27.;
				detd=4*c*c*c+27*d*d+4*d*b*b*b-b*b*c*c-18*b*c*d; 


				s1=0.;
				s2=0.;
				s3=0.;

				if (fabs(detd)<=1e-30){
				t=-q/2.;
				s1=2.*pow(t,1./3)-b/3.;
				s2=-pow(t,1./3)-b/3.;
				s3=s2;  			
				}
				else{
				r=sqrt(-p*p*p/27.);
				//			if((-q/(2*r)<-1)||(-q/(2*r)>1)){cout<<"Hamza est passé par là" <<(-q/(2*r))<<", "<<p<<", "<<detd<<endl;char quit;cin>>quit;}
				theta=acos(-q/(2*r));
				s1=2.*sqrt(-p/3.)*cos(theta/3.)-b/3.;
				s2=2.*sqrt(-p/3.)*cos((theta+2.*3.14159265358979323846)/3.)-b/3.;
				s3=2.*sqrt(-p/3.)*cos((theta+4.*3.14159265358979323846)/3.)-b/3.;

				}
			    
			    
                smax=max(max(s1,s2),s3);
                smin=min(min(s1,s2),s3);
                if(s1==smax){s2=max(s2,s3);}
                else if(s2==smax){s2=max(s1,s3);}
                else if(s3==smax){s2=max(s1,s2);}              
                
                s1=smax;
                s3=smin;
                trac=s1+s2+s3;                
				vmis=sqrt((sig11-sig22)*(sig11-sig22)+(sig33-sig22)*(sig33-sig22)+(sig11-sig33)*(sig11-sig33)+6.*(sig12*sig21+sig13*sig31+sig23*sig32))/sqrt(2.);
				if(vmis!=vmis) {vmis=0.;}    
                
                            
				mminvm=min(vmis,mminvm);
				mmaxvm=max(vmis,mmaxvm);
				
				mmintrac=min(trac,mmintrac);
				mmaxtrac=max(trac,mmaxtrac);			
				
				mminsig11=min(sig11,mminsig11);
				mmaxsig11=max(sig11,mmaxsig11);
				mminsig12=min(sig12,mminsig12);
				mmaxsig12=max(sig12,mmaxsig12);
				mminsig13=min(sig13,mminsig13);
				mmaxsig13=max(sig13,mmaxsig13);				
				mminsig22=min(sig22,mminsig22);
				mmaxsig22=max(sig22,mmaxsig22);							
				mminsig23=min(sig23,mminsig23);
				mmaxsig23=max(sig23,mmaxsig23);	
				mminsig33=min(sig33,mminsig33);
				mmaxsig33=max(sig33,mmaxsig33);						
								
				mminsig1=min(s1,mminsig1);
				mmaxsig1=max(s1,mmaxsig1);
				mminsig2=min(s2,mminsig2);
				mmaxsig2=max(s2,mmaxsig2);	
				mminsig3=min(s3,mminsig3);
				mmaxsig3=max(s3,mmaxsig3);		
							
				VONMIS[numcor]=vmis;  
				TRACE[numcor]=trac;                                    
				SIG11[numcor]=sig11;
				SIG12[numcor]=sig12;
				SIG13[numcor]=sig13;                                               
				SIG22[numcor]=sig22;
				SIG23[numcor]=sig23;                                      
				SIG33[numcor]=sig33;              
				SIG1[numcor]=s1;
				SIG2[numcor]=s2;           
				SIG3[numcor]=s3;    

	    
      }

minvm=mminvm;
maxvm=mmaxvm;
mintrac=mmintrac;
maxtrac=mmaxtrac;
minsig11=mminsig11;
maxsig11=mmaxsig11;
minsig12=mminsig12;
maxsig12=mmaxsig12;
minsig13=mminsig13;
maxsig13=mmaxsig13;
minsig22=mminsig22;
maxsig22=mmaxsig22;
minsig23=mminsig23;
maxsig23=mmaxsig23;
minsig33=mminsig33;
maxsig33=mmaxsig33;
minsig1=mminsig1;
maxsig1=mmaxsig1;
minsig2=mminsig2;
maxsig2=mmaxsig2;
minsig3=mminsig3;
maxsig3=mmaxsig3;

if(ite%NBENREG==0){ 
cout<<"Maxtrac:"<<maxtrac<<endl;
cout<<"Mintrac:"<<mintrac<<endl;
cout<<"voltot:"<<volt<<", "<<volt2<<endl;
}

}

void defoloc(R Pi, R coef1, R coef2, int NBENREG, R H_TOT, R V_TOT, R Z_TOT, int ite, int NB_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * DCONTXO, R * DCONTYO, R * DCONTZO, R * DCONTX, R * DCONTY, R * DCONTZ, int ** CONT, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * EPSI11, R * EPSI22, R * EPSI33, R * EPSE11, R * EPSE22, R * EPSE33, R & mindef11, R & maxdef11,R & mindef22, R & maxdef22,R & mindef33, R & maxdef33, R & mindefe11, R & maxdefe11,R & mindefe22, R & maxdefe22,R & mindefe33, R & maxdefe33, bool * EDGE, R * LIST_V, R * LIST_IND){

unsigned int numc;
int numcor,kt;
R ray;

R eps11,eps12,eps13,eps21,eps22,eps23,eps31,eps32,eps33;
R eps11e,eps12e,eps13e,eps21e,eps22e,eps23e,eps31e,eps32e,eps33e;
R n1,n2,n3;
R dxoij,dyoij,dzoij,dxoije,dyoije,dzoije;
int numcor1,numcor2;  

R moy=0.;
R epst11=0.;
R epst12=0.;
R epst13=0.;
R epst21=0.;
R epst22=0.;
R epst23=0.;
R epst31=0.;
R epst32=0.;
R epst33=0.; 
int num_dis=0;
R volt=0.;
R coefi,coefj,coefij;

maxdef11=-1e9;
mindef11=1e9;
maxdef22=-1e9;
mindef22=1e9;
maxdef33=-1e9;
mindef33=1e9;

maxdefe11=-1e9;
mindefe11=1e9;
maxdefe22=-1e9;
mindefe22=1e9;
maxdefe33=-1e9;
mindefe33=1e9;

int num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(numc,numcor,kt,ray,eps11,eps12,eps13,eps21,eps22,eps23,eps31,eps32,eps33,eps11e,eps12e,eps13e,eps21e,eps22e,eps23e,eps31e,eps32e,eps33e,n1,n2,n3,dxoij,dyoij,dzoij,dxoije,dyoije,dzoije,numcor1,numcor2) reduction(+:moy,epst11,epst12,epst13,epst21,epst22,epst23,epst31,epst32,epst33,num_dis,volt) reduction(min:mindef11,mindef22,mindef33,mindefe11,mindefe22,mindefe33) reduction(max:maxdef11,maxdef22,maxdef33,maxdefe11,maxdefe22,maxdefe33)
	for(numcor=0;numcor<NB_SPH;numcor++){ 
           
        num_dis++;
        
        eps11=0.;
        eps12=0.;
        eps13=0.;
        eps21=0.;
        eps22=0.;
        eps23=0.;
        eps31=0.;
        eps32=0.;
        eps33=0.;        
        
        eps11e=0.;
        eps12e=0.;
        eps13e=0.;
        eps21e=0.;
        eps22e=0.;
        eps23e=0.;
        eps31e=0.;
        eps32e=0.;
        eps33e=0.;   
      

        volt=volt+coef1*LIST_V[numcor];
             
             for(kt=0;kt<NBCONTCO[numcor];kt++){ 
             
		numc=NOCONT[numcor][kt];
		ray=LIST_R[numcor];
		n1=NCONT[numc][0];
		n2=NCONT[numc][1];
		n3=NCONT[numc][2];

		numcor1=CONT[numc][0];
		numcor2=CONT[numc][1];  

		dxoij=(LIST_X[numcor1]-LIST_X[numcor2])-DCONTXO[numc];        	           
		dyoij=(LIST_Y[numcor1]-LIST_Y[numcor2])-DCONTYO[numc];
		dzoij=(LIST_Z[numcor1]-LIST_Z[numcor2])-DCONTZO[numc];

  		dxoije=(LIST_X[numcor1]-LIST_X[numcor2])-DCONTX[numc];        	           
		dyoije=(LIST_Y[numcor1]-LIST_Y[numcor2])-DCONTY[numc];
		dzoije=(LIST_Z[numcor1]-LIST_Z[numcor2])-DCONTZ[numc];  
                        
             eps11e=eps11e+Pi*coef2*ray*ray*n1*dxoije/2.;
             eps12e=eps12e+Pi*coef2*ray*ray*n1*dyoije/2.;
             eps13e=eps13e+Pi*coef2*ray*ray*n1*dzoije/2.;
             eps21e=eps21e+Pi*coef2*ray*ray*n2*dxoije/2.;
             eps22e=eps22e+Pi*coef2*ray*ray*n2*dyoije/2.;
             eps23e=eps23e+Pi*coef2*ray*ray*n2*dzoije/2.;
             eps31e=eps31e+Pi*coef2*ray*ray*n3*dxoije/2.;
             eps32e=eps32e+Pi*coef2*ray*ray*n3*dyoije/2.;
             eps33e=eps33e+Pi*coef2*ray*ray*n3*dzoije/2.;  
  
             eps11=eps11+Pi*coef2*ray*ray*n1*dxoij/2.;
             eps12=eps12+Pi*coef2*ray*ray*n1*dyoij/2.;
             eps13=eps13+Pi*coef2*ray*ray*n1*dzoij/2.;
             eps21=eps21+Pi*coef2*ray*ray*n2*dxoij/2.;
             eps22=eps22+Pi*coef2*ray*ray*n2*dyoij/2.;
             eps23=eps23+Pi*coef2*ray*ray*n2*dzoij/2.;
             eps31=eps31+Pi*coef2*ray*ray*n3*dxoij/2.;
             eps32=eps32+Pi*coef2*ray*ray*n3*dyoij/2.;
             eps33=eps33+Pi*coef2*ray*ray*n3*dzoij/2.;    
                                        
             epst11=epst11+(Pi*coef2*ray*ray*n1*dxoij/2.);   
             epst12=epst12+(Pi*coef2*ray*ray*n1*dyoij/2.);
             epst13=epst13+(Pi*coef2*ray*ray*n1*dzoij/2.);
             epst21=epst21+(Pi*coef2*ray*ray*n2*dxoij/2.);
             epst22=epst22+(Pi*coef2*ray*ray*n2*dyoij/2.);
             epst23=epst23+(Pi*coef2*ray*ray*n2*dzoij/2.);
             epst31=epst31+(Pi*coef2*ray*ray*n3*dxoij/2.);
             epst32=epst32+(Pi*coef2*ray*ray*n3*dyoij/2.);
             epst33=epst33+(Pi*coef2*ray*ray*n3*dzoij/2.);                 
             }

			eps11=LIST_IND[numcor]*eps11/(coef1*LIST_V[numcor]);
			eps12=LIST_IND[numcor]*eps12/(coef1*LIST_V[numcor]);
			eps13=LIST_IND[numcor]*eps13/(coef1*LIST_V[numcor]);
			eps21=LIST_IND[numcor]*eps21/(coef1*LIST_V[numcor]);
			eps22=LIST_IND[numcor]*eps22/(coef1*LIST_V[numcor]);
			eps23=LIST_IND[numcor]*eps23/(coef1*LIST_V[numcor]);
			eps31=LIST_IND[numcor]*eps31/(coef1*LIST_V[numcor]);
			eps32=LIST_IND[numcor]*eps32/(coef1*LIST_V[numcor]);
			eps33=LIST_IND[numcor]*eps33/(coef1*LIST_V[numcor]);	

			eps11e=LIST_IND[numcor]*eps11e/(coef1*LIST_V[numcor]);
			eps12e=LIST_IND[numcor]*eps12e/(coef1*LIST_V[numcor]);
			eps13e=LIST_IND[numcor]*eps13e/(coef1*LIST_V[numcor]);
			eps21e=LIST_IND[numcor]*eps21e/(coef1*LIST_V[numcor]);
			eps22e=LIST_IND[numcor]*eps22e/(coef1*LIST_V[numcor]);
			eps23e=LIST_IND[numcor]*eps23e/(coef1*LIST_V[numcor]);
			eps31e=LIST_IND[numcor]*eps31e/(coef1*LIST_V[numcor]);
			eps32e=LIST_IND[numcor]*eps32e/(coef1*LIST_V[numcor]);
			eps33e=LIST_IND[numcor]*eps33e/(coef1*LIST_V[numcor]);

			EPSI11[numcor]=eps11;
			EPSI22[numcor]=eps22;
			EPSI33[numcor]=eps33;

                        if(maxdef11<eps11) maxdef11=eps11;
                        if(maxdef22<eps22) maxdef22=eps22;
                        if(maxdef33<eps33) maxdef33=eps33;
                        if(mindef11>eps11) mindef11=eps11;
                        if(mindef22>eps22) mindef22=eps22;
                        if(mindef33>eps33) mindef33=eps33;

			EPSE11[numcor]=eps11e;
			EPSE22[numcor]=eps22e;
			EPSE33[numcor]=eps33e;

                        if(maxdefe11<eps11e) maxdefe11=eps11e;
                        if(maxdefe22<eps22e) maxdefe22=eps22e;
                        if(maxdefe33<eps33e) maxdefe33=eps33e;
                        if(mindefe11>eps11e) mindefe11=eps11e;
                        if(mindefe22>eps22e) mindefe22=eps22e;
                        if(mindefe33>eps33e) mindefe33=eps33e;

			moy+=eps11;
	    
      }

   moy/=num_dis;    

			 epst11=epst11/(H_TOT*V_TOT*Z_TOT);
			 epst12=epst12/(H_TOT*V_TOT*Z_TOT);
			 epst13=epst13/(H_TOT*V_TOT*Z_TOT);
			 epst21=epst21/(H_TOT*V_TOT*Z_TOT);
			 epst22=epst22/(H_TOT*V_TOT*Z_TOT);
			 epst23=epst23/(H_TOT*V_TOT*Z_TOT);		
			 epst31=epst31/(H_TOT*V_TOT*Z_TOT);
			 epst32=epst32/(H_TOT*V_TOT*Z_TOT);
			 epst33=epst33/(H_TOT*V_TOT*Z_TOT);	

if(ite%NBENREG==0){ 
cout<<"maxdef11 : "<<maxdef11<<endl;
cout<<"mindef11 : "<<mindef11<<endl;
cout<<"voltot:"<<volt<<endl;
}


}



void defohalo(R Pi, R coef1, R coef2, int NBENREG, R H_TOT, R V_TOT, R Z_TOT, int ite, int NB_SPH, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * DCONTXO, R * DCONTYO, R * DCONTZO, R * DCONTX, R * DCONTY, R * DCONTZ, int ** CONT, R ** NCONT, unsigned int ** NOCONT, int * NBCONTCO, R * EPSI11, R * EPSI22, R * EPSI33, R * EPSE11, R * EPSE22, R * EPSE33, R & mindef11, R & maxdef11,R & mindef22, R & maxdef22,R & mindef33, R & maxdef33, R & mindefe11, R & maxdefe11,R & mindefe22, R & maxdefe22,R & mindefe33, R & maxdefe33, bool * EDGE, unsigned int ** NOHALO, int * NBHALO, R * LIST_V, R * VOLHALO, R * LIST_IND){

unsigned int numc;
int it,jt,kt;
unsigned int numcor; 
R ray;

R eps11,eps12,eps13,eps21,eps22,eps23,eps31,eps32,eps33;
R eps11e,eps12e,eps13e,eps21e,eps22e,eps23e,eps31e,eps32e,eps33e;
R n1,n2,n3;
R dxoij,dyoij,dzoij,dxoije,dyoije,dzoije;
int numcor1,numcor2;  

R coefij;

R moy=0.;
R moyt=0.;
int num_dis=0;
R volt=0.;

maxdef11=-1e9;
mindef11=1e9;
maxdef22=-1e9;
mindef22=1e9;
maxdef33=-1e9;
mindef33=1e9;

maxdefe11=-1e9;
mindefe11=1e9;
maxdefe22=-1e9;
mindefe22=1e9;
maxdefe33=-1e9;
mindefe33=1e9;

int num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(coefij,numc,numcor,it,jt,kt,ray,eps11,eps12,eps13,eps21,eps22,eps23,eps31,eps32,eps33,eps11e,eps12e,eps13e,eps21e,eps22e,eps23e,eps31e,eps32e,eps33e,n1,n2,n3,dxoij,dyoij,dzoij,dxoije,dyoije,dzoije,numcor1,numcor2) reduction(+:moy,moyt,num_dis,volt) reduction(min:mindef11,mindef22,mindef33,mindefe11,mindefe22,mindefe33) reduction(max:maxdef11,maxdef22,maxdef33,maxdefe11,maxdefe22,maxdefe33)
	for(jt=0;jt<NB_SPH;jt++){ 

        num_dis++;
        
        eps11=0.;
        eps12=0.;
        eps13=0.;
        eps21=0.;
        eps22=0.;
        eps23=0.;
        eps31=0.;
        eps32=0.;
        eps33=0.;        
        
        eps11e=0.;
        eps12e=0.;
        eps13e=0.;
        eps21e=0.;
        eps22e=0.;
        eps23e=0.;
        eps31e=0.;
        eps32e=0.;
        eps33e=0.;   
      
	for(it=0;it<NBHALO[jt];it++){
	numcor=NOHALO[jt][it];
	coefij=LIST_IND[numcor];

		for(kt=0;kt<NBCONTCO[numcor];kt++){ 

		numc=NOCONT[numcor][kt];
		ray=LIST_R[numcor];
		n1=NCONT[numc][0];
		n2=NCONT[numc][1];
		n3=NCONT[numc][2];

		numcor1=CONT[numc][0];
		numcor2=CONT[numc][1];  

		dxoij=(LIST_X[numcor1]-LIST_X[numcor2])-DCONTXO[numc];        	           
		dyoij=(LIST_Y[numcor1]-LIST_Y[numcor2])-DCONTYO[numc];
		dzoij=(LIST_Z[numcor1]-LIST_Z[numcor2])-DCONTZO[numc];

		dxoije=(LIST_X[numcor1]-LIST_X[numcor2])-DCONTX[numc];        	           
		dyoije=(LIST_Y[numcor1]-LIST_Y[numcor2])-DCONTY[numc];
		dzoije=(LIST_Z[numcor1]-LIST_Z[numcor2])-DCONTZ[numc];  

		eps11e=eps11e+Pi*coef2*coefij*ray*ray*n1*dxoije/2.;
		eps12e=eps12e+Pi*coef2*coefij*ray*ray*n1*dyoije/2.;
		eps13e=eps13e+Pi*coef2*coefij*ray*ray*n1*dzoije/2.;
		eps21e=eps21e+Pi*coef2*coefij*ray*ray*n2*dxoije/2.;
		eps22e=eps22e+Pi*coef2*coefij*ray*ray*n2*dyoije/2.;
		eps23e=eps23e+Pi*coef2*coefij*ray*ray*n2*dzoije/2.;
		eps31e=eps31e+Pi*coef2*coefij*ray*ray*n3*dxoije/2.;
		eps32e=eps32e+Pi*coef2*coefij*ray*ray*n3*dyoije/2.;
		eps33e=eps33e+Pi*coef2*coefij*ray*ray*n3*dzoije/2.;  

		eps11=eps11+Pi*coef2*coefij*ray*ray*n1*dxoij/2.;
		eps12=eps12+Pi*coef2*coefij*ray*ray*n1*dyoij/2.;
		eps13=eps13+Pi*coef2*coefij*ray*ray*n1*dzoij/2.;
		eps21=eps21+Pi*coef2*coefij*ray*ray*n2*dxoij/2.;
		eps22=eps22+Pi*coef2*coefij*ray*ray*n2*dyoij/2.;
		eps23=eps23+Pi*coef2*coefij*ray*ray*n2*dzoij/2.;
		eps31=eps31+Pi*coef2*coefij*ray*ray*n3*dxoij/2.;
		eps32=eps32+Pi*coef2*coefij*ray*ray*n3*dyoij/2.;
		eps33=eps33+Pi*coef2*coefij*ray*ray*n3*dzoij/2.;    
           
		}

	}

			eps11=eps11/VOLHALO[jt];
			eps12=eps12/VOLHALO[jt];
			eps13=eps13/VOLHALO[jt];
			eps21=eps21/VOLHALO[jt];
			eps22=eps22/VOLHALO[jt];
			eps23=eps23/VOLHALO[jt];
			eps31=eps31/VOLHALO[jt];
			eps32=eps32/VOLHALO[jt];
			eps33=eps33/VOLHALO[jt];	

			eps11e=eps11e/VOLHALO[jt];
			eps12e=eps12e/VOLHALO[jt];
			eps13e=eps13e/VOLHALO[jt];
			eps21e=eps21e/VOLHALO[jt];
			eps22e=eps22e/VOLHALO[jt];
			eps23e=eps23e/VOLHALO[jt];
			eps31e=eps31e/VOLHALO[jt];
			eps32e=eps32e/VOLHALO[jt];
			eps33e=eps33e/VOLHALO[jt];

			EPSI11[jt]=eps11;
			EPSI22[jt]=eps22;
			EPSI33[jt]=eps33;

                         maxdef11=max(maxdef11,eps11);
                         maxdef22=max(maxdef22,eps22);
                         maxdef33=max(maxdef33,eps33);
                         mindef11=min(mindef11,eps11);
                         mindef22=min(mindef22,eps22);
                         mindef33=min(mindef33,eps33);

			EPSE11[jt]=eps11e;
			EPSE22[jt]=eps22e;
			EPSE33[jt]=eps33e;

                         maxdefe11=max(maxdefe11,eps11e);
                         maxdefe22=max(maxdefe22,eps22e);
                         maxdefe33=max(maxdefe33,eps33e);
                         mindefe11=min(mindefe11,eps11e);
                         mindefe22=min(mindefe22,eps22e);
                         mindefe33=min(mindefe33,eps33e);

			moy+=eps11;
			moyt+=eps11*coef1*LIST_V[jt];
	    	        volt+=coef1*LIST_V[jt];    
      }

   moy/=num_dis;    
   moyt/=(H_TOT*V_TOT*Z_TOT);

if(ite%NBENREG==0){ 
cout<<"maxdef11 : "<<maxdef11<<endl;
cout<<"mindef11 : "<<mindef11<<endl;
cout<<"voltot:"<<volt<<endl;
}

}




void deploc(int NB_SPH, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_XO, R * LIST_YO, R * LIST_ZO, R * DEP1, R * DEP2, R * DEP3, R & maxdep1, R & mindep1, R & maxdep2, R & mindep2, R & maxdep3, R & mindep3){

R mmaxdep1=-1e12;
R mmindep1=1e12;
R mmaxdep2=-1e12;
R mmindep2=1e12;
R mmaxdep3=-1e12;
R mmindep3=1e12;

int num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

int it;

# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)) private(it) reduction(min:mmindep1,mmindep2,mmindep3) reduction(max:mmaxdep1,mmaxdep2,mmaxdep3)
	for(it=0;it<NB_SPH;it++){ 
		
		DEP1[it]=LIST_X[it]-LIST_XO[it];
		DEP2[it]=LIST_Y[it]-LIST_YO[it];
		DEP3[it]=LIST_Z[it]-LIST_ZO[it];	
		
		mmindep1=min(DEP1[it],mmindep1);
		mmaxdep1=max(DEP1[it],mmaxdep1);
		mmindep2=min(DEP2[it],mmindep2);
		mmaxdep2=max(DEP2[it],mmaxdep2);	
		mmindep3=min(DEP3[it],mmindep3);
		mmaxdep3=max(DEP3[it],mmaxdep3);		
		
	}

mindep1=mmindep1;
maxdep1=mmaxdep1;
mindep2=mmindep2;
maxdep2=mmaxdep2;	
mindep3=mmindep3;
maxdep3=mmaxdep3;	

}

