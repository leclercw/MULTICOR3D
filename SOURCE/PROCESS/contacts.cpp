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

#include "contacts.h"
#include "omp.h"

void force(R dt, int NBCO, int NB_SPH, int ** CONT,  R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONT,R * DCONTO, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI,R ** MTJI,R ** MTIJ, R E, R nu, R fric, R &fpar,R &fpar2, int NBCOP, int ** CONTP, R * DCONTP, R ** NCONTP, R * LIST_PM, R ** LIST_PVX, R ** LIST_PVY, R ** LIST_PVZ){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,mxij,mxji,myij,myji,mzij,mzji;
R fxx,fyy,fzz;
R abspred2;
R dist,rn,rnv,crn,crnv,drn,drnv,mas1,pred1,pred2,pred3,u1,u2,u3,Kn,Kt;

R Eij,iEij,rij,hij;

int numc,nump;
R pred1t;
R vpx,vpy,vpz;	
fpar=0.;
fpar2=0.;

int it,jt,kt,lt,num_threads,numth;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}


if(num_threads==1){

# pragma omp parallel for schedule(static) private(it,kt,lt,numth,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,dist,mas1,u1,u2,u3,rij,hij,Eij,iEij,Kn,Kt,rn,rnv,crn,crnv,drn,drnv,abspred2,pred1,pred2,pred3,mxij,mxji,myij,myji,mzij,mzji,fxx,fyy,fzz) reduction(+:FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH])
	for(it=0;it<NBCO;it++){ 

		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     
		
	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
	     
	  	 dist=DCONT[it];     
	     mas1=LIST_M[numc1]*LIST_M[numc2]/(LIST_M[numc1]+LIST_M[numc2]);      
 
        // Contact sph-sph (Non-cohésif)                   
	    
         u1=0.;
	     u2=-LIST_R[numc1]*t1*LIST_WX[numc1]-LIST_R[numc1]*t2*LIST_WY[numc1]-LIST_R[numc1]*t3*LIST_WZ[numc1];
	     u2=u2-LIST_R[numc2]*t1*LIST_WX[numc2]-LIST_R[numc2]*t2*LIST_WY[numc2]-LIST_R[numc2]*t3*LIST_WZ[numc2];
	     u3=LIST_R[numc1]*s1*LIST_WX[numc1]+LIST_R[numc1]*s2*LIST_WY[numc1]+LIST_R[numc1]*s3*LIST_WZ[numc1];
	     u3=u3+LIST_R[numc2]*s1*LIST_WX[numc2]+LIST_R[numc2]*s2*LIST_WY[numc2]+LIST_R[numc2]*s3*LIST_WZ[numc2];	  
	     
	     u1=u1+n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u2=u2+s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u3=u3+t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	 

	     rij=(LIST_R[numc1]*LIST_R[numc2])/(LIST_R[numc1]+LIST_R[numc2]);
	     hij=abs(dist-(LIST_R[numc1]+LIST_R[numc2]));

		 iEij=2.*(1.-nu*nu)/E;		 
		 Eij=1./iEij;	
	     Kn=4.*Eij*sqrt(rij*hij)/3.;
	     Kn=Kn/100.;
	     Kn=1e6;
	     Kt=Kn;
	     	     	 	     
	     rn=Kn*(dist-(LIST_R[numc1]+LIST_R[numc2]));
	     crn=Kt*u2*dt;
	     drn=Kt*u3*dt;
		  
	     rnv=sqrt(Kn*mas1)*u1;   
	     crnv=sqrt(Kn*mas1)*u2;  
	     drnv=sqrt(Kn*mas1)*u3;  

	     pred1=-rn-rnv;
	     pred2=-crn-crnv;    
	     pred3=-drn-drnv;    

	     abspred2=sqrt(pred2*pred2+pred3*pred3);
	       
	     if(fric==0){pred2=0.;pred3=0.;}

	     if(pred1<0){ //relachement
		 pred1=0.;
		 pred2=0.;	 
		 pred3=0.;	 
		 }
		 else if(abspred2<=fric*pred1){ //adherence
		 pred1=pred1;
		 pred2=pred2;	 
		 pred3=pred3;	 
		 }
		 else if(fric>0){ // Glissement avec frottement
		 pred1=pred1;
		 pred2=-abs(fric*pred1)*u2/abs(u2); 
		 pred3=-abs(fric*pred1)*u3/abs(u3); 
		 }	 
		 else // Glissement sans frottement
		 {
		 pred1=pred1;
		 pred2=0.;
		 pred3=0.;
	         }
 
	          
          fxx=n1*pred1+s1*pred2+t1*pred3;
          fyy=n2*pred1+s2*pred2+t2*pred3;
          fzz=n3*pred1+s3*pred2+t3*pred3;
		  
		  mxji = -LIST_R[numc1]*t1*pred2+LIST_R[numc1]*s1*pred3;
		  mxij = -LIST_R[numc2]*t1*pred2+LIST_R[numc2]*s1*pred3;
		  myji = -LIST_R[numc1]*t2*pred2+LIST_R[numc1]*s2*pred3;
		  myij = -LIST_R[numc2]*t2*pred2+LIST_R[numc2]*s2*pred3;
		  mzji = -LIST_R[numc1]*t3*pred2+LIST_R[numc1]*s3*pred3; 
		  mzij = -LIST_R[numc2]*t3*pred2+LIST_R[numc2]*s3*pred3;
			 
			FCJI[it][0]=fxx;
			FCJI[it][1]=fyy;
			FCJI[it][2]=fzz; 

			  FX[numc1]+=fxx;
			  FY[numc1]+=fyy;
			  FZ[numc1]+=fzz;
			  FX[numc2]+=(-fxx);
			  FY[numc2]+=(-fyy);
			  FZ[numc2]+=(-fzz);
			  MTX[numc1]+=mxji;
			  MTX[numc2]+=mxij;
			  MTY[numc1]+=myji;
			  MTY[numc2]+=myij;	
			  MTZ[numc1]+=mzji;
			  MTZ[numc2]+=mzij;

			          
   }



    


}
else{

R FX2[NB_SPH][num_threads-1];
R FY2[NB_SPH][num_threads-1];
R FZ2[NB_SPH][num_threads-1];
R MTX2[NB_SPH][num_threads-1];
R MTY2[NB_SPH][num_threads-1];
R MTZ2[NB_SPH][num_threads-1];

# pragma omp parallel for collapse(2)
	for(it=0;it<NB_SPH;it++){ 

		for(jt=0;jt<num_threads-1;jt++){ 

		 FX2[it][jt]=0;
		 FY2[it][jt]=0;
		 FZ2[it][jt]=0;
		 MTX2[it][jt]=0;
		 MTY2[it][jt]=0;
		 MTZ2[it][jt]=0;

		}
	}



# pragma omp parallel for schedule(static) private(it,kt,lt,numth,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,dist,mas1,u1,u2,u3,rij,hij,Eij,iEij,Kn,Kt,rn,rnv,crn,crnv,drn,drnv,abspred2,pred1,pred2,pred3,mxij,mxji,myij,myji,mzij,mzji,fxx,fyy,fzz) 
	for(it=0;it<NBCO;it++){ 

		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     
		
	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
	     
	  	 dist=DCONT[it];     
	     mas1=LIST_M[numc1]*LIST_M[numc2]/(LIST_M[numc1]+LIST_M[numc2]);      
 
        // Contact sph-sph (Non-cohésif)                   
	    
         u1=0.;
	     u2=-LIST_R[numc1]*t1*LIST_WX[numc1]-LIST_R[numc1]*t2*LIST_WY[numc1]-LIST_R[numc1]*t3*LIST_WZ[numc1];
	     u2=u2-LIST_R[numc2]*t1*LIST_WX[numc2]-LIST_R[numc2]*t2*LIST_WY[numc2]-LIST_R[numc2]*t3*LIST_WZ[numc2];
	     u3=LIST_R[numc1]*s1*LIST_WX[numc1]+LIST_R[numc1]*s2*LIST_WY[numc1]+LIST_R[numc1]*s3*LIST_WZ[numc1];
	     u3=u3+LIST_R[numc2]*s1*LIST_WX[numc2]+LIST_R[numc2]*s2*LIST_WY[numc2]+LIST_R[numc2]*s3*LIST_WZ[numc2];	  
	     
	     u1=u1+n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u2=u2+s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u3=u3+t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	 

	     rij=(LIST_R[numc1]*LIST_R[numc2])/(LIST_R[numc1]+LIST_R[numc2]);
	     hij=abs(dist-(LIST_R[numc1]+LIST_R[numc2]));

		 iEij=2.*(1.-nu*nu)/E;		 
		 Eij=1./iEij;	
	     Kn=4.*Eij*sqrt(rij*hij)/3.;
	     Kn=Kn/100.;
	     Kn=1e6;
	     Kt=Kn;
	     	     	 	     
	     rn=Kn*(dist-(LIST_R[numc1]+LIST_R[numc2]));
	     crn=Kt*u2*dt;
	     drn=Kt*u3*dt;
		  
	     rnv=sqrt(Kn*mas1)*u1;   
	     crnv=sqrt(Kn*mas1)*u2;  
	     drnv=sqrt(Kn*mas1)*u3;  

	     pred1=-rn-rnv;
	     pred2=-crn-crnv;    
	     pred3=-drn-drnv;    

	     abspred2=sqrt(pred2*pred2+pred3*pred3);
	       
	     if(fric==0){pred2=0.;pred3=0.;}

	     if(pred1<0){ //relachement
		 pred1=0.;
		 pred2=0.;	 
		 pred3=0.;	 
		 }
		 else if(abspred2<=fric*pred1){ //adherence
		 pred1=pred1;
		 pred2=pred2;	 
		 pred3=pred3;	 
		 }
		 else if(fric>0){ // Glissement avec frottement
		 pred1=pred1;
		 pred2=-abs(fric*pred1)*u2/abs(u2); 
		 pred3=-abs(fric*pred1)*u3/abs(u3); 
		 }	 
		 else // Glissement sans frottement
		 {
		 pred1=pred1;
		 pred2=0.;
		 pred3=0.;
	         }
 
	          
          fxx=n1*pred1+s1*pred2+t1*pred3;
          fyy=n2*pred1+s2*pred2+t2*pred3;
          fzz=n3*pred1+s3*pred2+t3*pred3;
		  
		  mxji = -LIST_R[numc1]*t1*pred2+LIST_R[numc1]*s1*pred3;
		  mxij = -LIST_R[numc2]*t1*pred2+LIST_R[numc2]*s1*pred3;
		  myji = -LIST_R[numc1]*t2*pred2+LIST_R[numc1]*s2*pred3;
		  myij = -LIST_R[numc2]*t2*pred2+LIST_R[numc2]*s2*pred3;
		  mzji = -LIST_R[numc1]*t3*pred2+LIST_R[numc1]*s3*pred3; 
		  mzij = -LIST_R[numc2]*t3*pred2+LIST_R[numc2]*s3*pred3;
			 
			FCJI[it][0]=fxx;
			FCJI[it][1]=fyy;
			FCJI[it][2]=fzz; 

			  numth=omp_get_thread_num()-1;

		if(numth==-1){
			  FX[numc1]+=fxx;
			  FY[numc1]+=fyy;
			  FZ[numc1]+=fzz;
			  FX[numc2]+=(-fxx);
			  FY[numc2]+=(-fyy);
			  FZ[numc2]+=(-fzz);
			  MTX[numc1]+=mxji;
			  MTX[numc2]+=mxij;
			  MTY[numc1]+=myji;
			  MTY[numc2]+=myij;	
			  MTZ[numc1]+=mzji;
			  MTZ[numc2]+=mzij;
		}else{


			  FX2[numc1][numth]+=fxx;
			  FY2[numc1][numth]+=fyy;
			  FZ2[numc1][numth]+=fzz;
			  FX2[numc2][numth]+=(-fxx);
			  FY2[numc2][numth]+=(-fyy);
			  FZ2[numc2][numth]+=(-fzz);
			  MTX2[numc1][numth]+=mxji;
			  MTX2[numc2][numth]+=mxij;
			  MTY2[numc1][numth]+=myji;
			  MTY2[numc2][numth]+=myij;	
			  MTZ2[numc1][numth]+=mzji;
			  MTZ2[numc2][numth]+=mzij;


		}
			          
   }

# pragma omp parallel for collapse(2) 
	for(it=0;it<NB_SPH;it++){
 
		for(jt=0;jt<num_threads-1;jt++){ 

			FX[it]+=FX2[it][jt];
			FY[it]+=FY2[it][jt];
			FZ[it]+=FZ2[it][jt];
			MTX[it]+=MTX2[it][jt];
			MTY[it]+=MTY2[it][jt];
			MTZ[it]+=MTZ2[it][jt];
		}


        }





}






	for(it=0;it<NBCOP;it++){ 
		
		 // Numeros des candidats 
	     nump=CONTP[it][0];
	     numc=CONTP[it][1];     
	     
	     // normal ext 1->2
	     n1=NCONTP[it][0];
	     n2=NCONTP[it][1];
	     n3=NCONTP[it][2];
	     s1=NCONTP[it][3];
	     s2=NCONTP[it][4];
	     s3=NCONTP[it][5];
	     t1=NCONTP[it][6];
	     t2=NCONTP[it][7];
	     t3=NCONTP[it][8];  		     

	  	 dist=DCONTP[it];     
	     //mas1=LIST_M[numc]*LIST_M[numc]/(LIST_M[numc]+LIST_M[numc]);  
	     mas1=LIST_PM[nump]*1.;
	     
	     vpx=(LIST_PVX[nump][0]+LIST_PVX[nump][1]+LIST_PVX[nump][2]+LIST_PVX[nump][3])/4.;
	     vpy=(LIST_PVY[nump][0]+LIST_PVY[nump][1]+LIST_PVY[nump][2]+LIST_PVY[nump][3])/4.;
	     vpz=(LIST_PVZ[nump][0]+LIST_PVZ[nump][1]+LIST_PVZ[nump][2]+LIST_PVZ[nump][3])/4.;	     	     
	     
        // Contact sph-par (Non-cohésif)                   
	     
         u1=0.;
	     u2=u2-LIST_R[numc]*t1*LIST_WX[numc]-LIST_R[numc]*t2*LIST_WY[numc]-LIST_R[numc]*t3*LIST_WZ[numc];
	     u3=u3+LIST_R[numc]*s1*LIST_WX[numc]+LIST_R[numc]*s2*LIST_WY[numc]+LIST_R[numc]*s3*LIST_WZ[numc];	  
	     
	     u1=u1+n1*(vpx-LIST_VX[numc])+n2*(vpy-LIST_VY[numc])+n3*(vpz-LIST_VZ[numc]);
	     u2=u2+s1*(vpx-LIST_VX[numc])+s2*(vpy-LIST_VY[numc])+s3*(vpz-LIST_VZ[numc]);
	     u3=u3+t1*(vpx-LIST_VX[numc])+t2*(vpy-LIST_VY[numc])+t3*(vpz-LIST_VZ[numc]);	 

	     rij=LIST_R[numc];
	     hij=fabs(dist-LIST_R[numc]);
		 iEij=(1.-nu*nu)/E;		 
		 Eij=1./iEij;	
		 Kn=4.*Eij*sqrt(rij*hij)/3.;
		 //Kn=4.*Eij*rij;   
	//     Kn=3.35e7; 1 PAROI
	//     Kn=1.155e7; 2 parois
	         Kn=1e7;
	   	 Kt=Kn;
	     	     	 	     
	     rn=Kn*(dist-LIST_R[numc]);
	     crn=Kt*u2*dt;
	     drn=Kt*u3*dt;
		  
	     rnv=sqrt(Kn*mas1)*u1;   
	     crnv=sqrt(Kn*mas1)*u2;  
	     drnv=sqrt(Kn*mas1)*u3;  

	     pred1=-rn-rnv;
	     pred1t=pred1;
	     pred2=-crn-crnv;    
	     pred3=-drn-drnv;    

	     abspred2=sqrt(pred2*pred2+pred3*pred3);

	     if(fric==0){pred2=0.;pred3=0.;}

	     if(pred1<0){ //relachement
		 pred1=0.;
		 pred2=0.;	 
		 pred3=0.;	 
		 }
		 else if(abspred2<=fric*pred1){ //adherence
		 pred1=pred1;
		 pred2=pred2;	 
		 pred3=pred3;	 
		 }
		 else if(fric>0){ // Glissement avec frottement
		 pred1=pred1;
		 pred2=-abs(fric*pred1)*u2/abs(u2); 
		 pred3=-abs(fric*pred1)*u3/abs(u3); 
		 }	 
		 else // Glissement sans frottement
		 {
		 pred1=pred1;
		 pred2=0.;
		 pred3=0.;
	     }
		          
		fxx=n1*pred1+s1*pred2+t1*pred3;
		fyy=n2*pred1+s2*pred2+t2*pred3;
		fzz=n3*pred1+s3*pred2+t3*pred3;		  

		mxij = -LIST_R[numc]*t1*pred2+LIST_R[numc]*s1*pred3;
		myij = -LIST_R[numc]*t2*pred2+LIST_R[numc]*s2*pred3;
		mzij = -LIST_R[numc]*t3*pred2+LIST_R[numc]*s3*pred3;	     

		FX[numc]+=(-fxx);
		FY[numc]+=(-fyy);
		FZ[numc]+=(-fzz);			    
			 
		if(nump==0){fpar=fpar+fxx;}	 
		if(nump==1){fpar2=fpar2+fxx;}	 
					 
		MTX[numc]+=mxij;
		MTY[numc]+=myij;	
		MTZ[numc]+=mzij;	
	     	
    }




}


void forcecoh_ind(R & Epe, R & Epa,R dt, bool * TYPCO,int NBCO,  int ** CONT,  R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONT, R * DCONTO, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI,R ** MTJI,R ** MTIJ, R E1, R nu1,R E2, R nu2, R fric, int NB_SPH){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij,abspred2;
R mnji,mtji,mbji;
R rij,hij;
R dist,rn,rnv,crn,crnv,drn,drnv,mas1,pred1,pred2,pred3,u1,u2,u3,Kn,Kt;
R Eij,iEij;
Epe=0.;
Epa=0.;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)+1) private(it,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij,dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,abspred2,mnji,mtji,mbji,rij,hij, Eij,iEij) reduction(+:Epe,Epa,FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH])
	for(it=0;it<NBCO;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     
		
	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
	     
	  	 dist=DCONT[it];     
	     mas1=LIST_M[numc1]*LIST_M[numc2]/(LIST_M[numc1]+LIST_M[numc2]);      
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			/* vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
             */        
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			fnji=-VALCOH[it][0]*un_ij;
			ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  

			mnji=-VALCOH[it][5]*(thetani-thetanj);
			mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
					 
			mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
			mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;


// Modèle Kn/Kt
/*
			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij;
			  fbji=-VALCOH[it][1]*ub_ij;  
			  
			  mnji=0.;
			  mtji=0.;
			  mbji=0.;
			 			 
			  mtij=0.;  
              mbij=0.;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;	*/

// Amortissement 
/*
			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  
			  
			  mnji=mnji-(VALAMO[it][5]*(LIST_WX[numc2]-LIST_WX[numc1]));
			  mtji=mtji-(VALAMO[it][3]*LIST_WY[numc1]+VALAMO[it][4]*LIST_WY[numc2]);
			  mbji=mbji-(VALAMO[it][3]*LIST_WZ[numc1]+VALAMO[it][4]*LIST_WZ[numc2]);
			  
			  mtij=mtij-(VALAMO[it][3]*LIST_WY[numc2]+VALAMO[it][4]*LIST_WY[numc1]);
			  mbij=mbij-(VALAMO[it][3]*LIST_WZ[numc2]+VALAMO[it][4]*LIST_WZ[numc1]);*/

/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/
		
              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  		
		
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
		 	  

	
	      }
          else{
			  
        // Contact sph-sph (Non-cohésif)                   
	    fxx=0.;
		fyy=0.;
		fzz=0.;

		mxji=0.;
		myji=0.;
		mzji=0.;

		mxij=0.;
		myij=0.;
		mzij=0.;
		
       	     u1=0.;
	     u2=-LIST_R[numc1]*t1*LIST_WX[numc1]-LIST_R[numc1]*t2*LIST_WY[numc1]-LIST_R[numc1]*t3*LIST_WZ[numc1];
	     u2=u2-LIST_R[numc2]*t1*LIST_WX[numc2]-LIST_R[numc2]*t2*LIST_WY[numc2]-LIST_R[numc2]*t3*LIST_WZ[numc2];
	     u3=LIST_R[numc1]*s1*LIST_WX[numc1]+LIST_R[numc1]*s2*LIST_WY[numc1]+LIST_R[numc1]*s3*LIST_WZ[numc1];
	     u3=u3+LIST_R[numc2]*s1*LIST_WX[numc2]+LIST_R[numc2]*s2*LIST_WY[numc2]+LIST_R[numc2]*s3*LIST_WZ[numc2];	  
	     
	     u1=u1+n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u2=u2+s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u3=u3+t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	 

	     rij=(LIST_R[numc1]*LIST_R[numc2])/(LIST_R[numc1]+LIST_R[numc2]);
	     hij=fabs(dist-DCONTO[it]);

             if((numc1==NB_SPH-1)||(numc2==NB_SPH-1)){
             //contacts indenteur
	     iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;
	     Eij=1./iEij;	
	     Kn=10.*4.*Eij*sqrt(rij*hij)/3.;
		//cout<<"Kn:"<<Kn<<endl;	     

             }else 
             {
             //contacts dans l'échantillon
             //  cout<<"contacts dans l'échantillon"<<endl;

	     iEij=2.*(1.-nu2*nu2)/E2;
	     Eij=1./iEij;	
	     Kn=4.*Eij*sqrt(rij*hij)/3.; 
             Kn=0.;
             }	
             Kt=Kn;
	     	     	 	     
	     rn=Kn*(dist-DCONTO[it]);
	     crn=Kt*u2*dt;
	     drn=Kt*u3*dt;
		  
	     rnv=sqrt(Kn*mas1)*u1;   
	     crnv=sqrt(Kn*mas1)*u2;  
	     drnv=sqrt(Kn*mas1)*u3;  



	     pred1=-rn-rnv;
	     pred2=-crn-crnv;    
	     pred3=-drn-drnv;    

	     abspred2=sqrt(pred2*pred2+pred3*pred3);

	     if(fric==0){pred2=0.;pred3=0.;}

	     if(pred1<0){ //relachement
		 pred1=0.;
		 pred2=0.;	 
		 pred3=0.;	 
		 }
		 else if(abspred2<=fric*pred1){ //adherence
		 pred1=pred1;
		 pred2=pred2;	 
		 pred3=pred3;	 
		 }
		 else if(fric>0){ // Glissement avec frottement
		 pred1=pred1;
		 pred2=-abs(fric*pred1)*u2/abs(u2); 
		 pred3=-abs(fric*pred1)*u3/abs(u3); 
		 }	 
		 else // Glissement sans frottement
		 {
		 pred1=pred1;
		 pred2=0.;
		 pred3=0.;
	     }
	          
          fxx=n1*pred1+s1*pred2+t1*pred3;
          fyy=n2*pred1+s2*pred2+t2*pred3;
          fzz=n3*pred1+s3*pred2+t3*pred3;
		  
		  mxji = -LIST_R[numc1]*t1*pred2+LIST_R[numc1]*s1*pred3;
		  mxij = -LIST_R[numc2]*t1*pred2+LIST_R[numc2]*s1*pred3;
		  myji = -LIST_R[numc1]*t2*pred2+LIST_R[numc1]*s2*pred3;
		  myij = -LIST_R[numc2]*t2*pred2+LIST_R[numc2]*s2*pred3;
		  mzji = -LIST_R[numc1]*t3*pred2+LIST_R[numc1]*s3*pred3; 
		  mzij = -LIST_R[numc2]*t3*pred2+LIST_R[numc2]*s3*pred3;
			 
	      }   


			  FCJI[it][0]=fxx;
			  FCJI[it][1]=fyy;
			  FCJI[it][2]=fzz; 

			  FX[numc1]=FX[numc1]+fxx;			 
			  FY[numc1]=FY[numc1]+fyy;			  
			  FZ[numc1]=FZ[numc1]+fzz;

			  FX[numc2]=FX[numc2]-fxx;
			  FY[numc2]=FY[numc2]-fyy;
			  FZ[numc2]=FZ[numc2]-fzz;			    
					 
			  MTX[numc1]=MTX[numc1]+mxji;
			  MTX[numc2]=MTX[numc2]+mxij;

			  MTY[numc1]=MTY[numc1]+myji;
			  MTY[numc2]=MTY[numc2]+myij;			  

			  MTZ[numc1]=MTZ[numc1]+mzji;
			  MTZ[numc2]=MTZ[numc2]+mzij;	
			           
        }	
  
}

void forcecoh_qs(int ite, R & Epe, R & Epa,int NBCO, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ, R amort){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij;
R mnji,mtji,mbji;
R Epe1=0.;
R Epa1=0.;
int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)+1) private(it,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij,dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,mnji,mtji,mbji) reduction(+:Epe1,Epa1,FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH])
	for(it=0;it<NBCO;it++){ 
		
	  // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];   
		
	  // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];


         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
                      
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe1+=-fnji*un_ij/2.;
			Epe1+=-ftji*ut_ij/2.;                      
			Epe1+=-fbji*ub_ij/2.;      

			Epe1+=-mnji*(thetani)/2.;    
			Epe1+=-mtji*(thetati)/2.;                   
			Epe1+=-mbji*(thetabi)/2.;                     

			Epe1+=-(-mnji)*(thetanj)/2.;                      
			Epe1+=-mtij*(thetatj)/2.;                   
			Epe1+=-mbij*(thetabj)/2.;


if(amort!=0.){
	// Vitesses relatives dans le repère local
		vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
		vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
		vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	

	// Amortissement 

		fnji=fnji-(VALAMO[it][0]*vn_ij);
		ftji=ftji-(VALAMO[it][1]*vt_ij);
		fbji=fbji-(VALAMO[it][1]*vb_ij);		  

		Epa1+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
		Epa1+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
		Epa1+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
}
				

			FOJI[it][0]=fnji;
			FOJI[it][1]=ftji;
			FOJI[it][2]=fbji;

			MTJI[it][0]=mnji;
			MTJI[it][1]=mtji;
			MTJI[it][2]=mbji;                                      

			MTIJ[it][0]=-mnji;
			MTIJ[it][1]=mtij;
			MTIJ[it][2]=mbij;  

			fxx=n1*fnji+s1*ftji+t1*fbji;
			fyy=n2*fnji+s2*ftji+t2*fbji;
			fzz=n3*fnji+s3*ftji+t3*fbji;

			FCJI[it][0]=fxx;
			FCJI[it][1]=fyy;
			FCJI[it][2]=fzz;  

			mxji=n1*mnji+s1*mtji+t1*mbji; 
			myji=n2*mnji+s2*mtji+t2*mbji; 
			mzji=n3*mnji+s3*mtji+t3*mbji; 

			mxij=-n1*mnji+s1*mtij+t1*mbij; 
			myij=-n2*mnji+s2*mtij+t2*mbij; 
			mzij=-n3*mnji+s3*mtij+t3*mbij; 

		//	  # pragma omp atomic
			  FX[numc1]+=fxx;
		//	  # pragma omp atomic
			  FY[numc1]+=fyy;
		//	  # pragma omp atomic
			  FZ[numc1]+=fzz;
		//	  # pragma omp atomic
			  FX[numc2]+=(-fxx);
		//	  # pragma omp atomic
			  FY[numc2]+=(-fyy);
		//	  # pragma omp atomic
			  FZ[numc2]+=(-fzz);
		//	  # pragma omp atomic
			  MTX[numc1]+=mxji;
		//	  # pragma omp atomic
			  MTX[numc2]+=mxij;
		//	  # pragma omp atomic
			  MTY[numc1]+=myji;
		//	  # pragma omp atomic
			  MTY[numc2]+=myij;	
		//	  # pragma omp atomic
			  MTZ[numc1]+=mzji;
		//	  # pragma omp atomic
			  MTZ[numc2]+=mzij;
    
        }	
  

Epe=Epe1;
Epa=Epa1;

}

void forcecoh_qs_test(int ite, R & Epe, R & Epa,int NBCO, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ, R amort){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij;
R mnji,mtji,mbji;
R Epe1=0.;
R Epa1=0.;
int it,jt,numth,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){
	for(it=0;it<NBCO;it++){ 
		
	  // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];   
		
	  // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
                      
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
            
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe1+=-fnji*un_ij/2.;
			Epe1+=-ftji*ut_ij/2.;                      
			Epe1+=-fbji*ub_ij/2.;      

			Epe1+=-mnji*(thetani)/2.;    
			Epe1+=-mtji*(thetati)/2.;                   
			Epe1+=-mbji*(thetabi)/2.;                     

			Epe1+=-(-mnji)*(thetanj)/2.;                      
			Epe1+=-mtij*(thetatj)/2.;                   
			Epe1+=-mbij*(thetabj)/2.;


if(amort!=0.){
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	

// Amortissement 
			fnji=fnji-(VALAMO[it][0]*vn_ij);
			ftji=ftji-(VALAMO[it][1]*vt_ij);
			fbji=fbji-(VALAMO[it][1]*vb_ij);		  

			Epa1+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
			Epa1+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
			Epa1+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
}

				

			FOJI[it][0]=fnji;
			FOJI[it][1]=ftji;
			FOJI[it][2]=fbji;

			MTJI[it][0]=mnji;
			MTJI[it][1]=mtji;
			MTJI[it][2]=mbji;                                      

			MTIJ[it][0]=-mnji;
			MTIJ[it][1]=mtij;
			MTIJ[it][2]=mbij;  

			fxx=n1*fnji+s1*ftji+t1*fbji;
			fyy=n2*fnji+s2*ftji+t2*fbji;
			fzz=n3*fnji+s3*ftji+t3*fbji;

			FCJI[it][0]=fxx;
			FCJI[it][1]=fyy;
			FCJI[it][2]=fzz;  

			mxji=n1*mnji+s1*mtji+t1*mbji; 
			myji=n2*mnji+s2*mtji+t2*mbji; 
			mzji=n3*mnji+s3*mtji+t3*mbji; 

			mxij=-n1*mnji+s1*mtij+t1*mbij; 
			myij=-n2*mnji+s2*mtij+t2*mbij; 
			mzij=-n3*mnji+s3*mtij+t3*mbij; 

			FX[numc1]+=fxx;
			FY[numc1]+=fyy;
			FZ[numc1]+=fzz;
			FX[numc2]+=(-fxx);
			FY[numc2]+=(-fyy);
			FZ[numc2]+=(-fzz);
			MTX[numc1]+=mxji;
			MTX[numc2]+=mxij;
			MTY[numc1]+=myji;
			MTY[numc2]+=myij;
			MTZ[numc1]+=mzji;
			MTZ[numc2]+=mzij;
    
       }	

}else{

R FXP[num_threads-1][NB_SPH];
R FYP[num_threads-1][NB_SPH];
R FZP[num_threads-1][NB_SPH];

R MTXP[num_threads-1][NB_SPH];
R MTYP[num_threads-1][NB_SPH];
R MTZP[num_threads-1][NB_SPH];

# pragma omp parallel for schedule(static,1) private(it,jt)
for(it=0;it<num_threads-1;it++){
	for(jt=0;jt<NB_SPH;jt++){
	FXP[it][jt]=0;
	FYP[it][jt]=0;
	FZP[it][jt]=0;
	MTXP[it][jt]=0;
	MTYP[it][jt]=0;
	MTZP[it][jt]=0;
	}
}

//cout<<int(NBCO/num_threads)+1<< " " <<NBCO<<endl;
# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)+1) private(it,numth,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij,dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,mnji,mtji,mbji) reduction(+:Epe1,Epa1)
	for(it=0;it<NBCO;it++){ 
		
	numth=omp_get_thread_num()-1;

	  // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];   
		
	  // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
                      
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
            
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe1+=-fnji*un_ij/2.;
			Epe1+=-ftji*ut_ij/2.;                      
			Epe1+=-fbji*ub_ij/2.;      

			Epe1+=-mnji*(thetani)/2.;    
			Epe1+=-mtji*(thetati)/2.;                   
			Epe1+=-mbji*(thetabi)/2.;                     

			Epe1+=-(-mnji)*(thetanj)/2.;                      
			Epe1+=-mtij*(thetatj)/2.;                   
			Epe1+=-mbij*(thetabj)/2.;


// Amortissement 
if(amort!=0.){
        
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	


			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  

			//	Epa1+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
			//	Epa1+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
			//	Epa1+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 			
}
		
			FOJI[it][0]=fnji;
			FOJI[it][1]=ftji;
			FOJI[it][2]=fbji;

			MTJI[it][0]=mnji;
			MTJI[it][1]=mtji;
			MTJI[it][2]=mbji;                                      

			MTIJ[it][0]=-mnji;
			MTIJ[it][1]=mtij;
			MTIJ[it][2]=mbij;  

			fxx=n1*fnji+s1*ftji+t1*fbji;
			fyy=n2*fnji+s2*ftji+t2*fbji;
			fzz=n3*fnji+s3*ftji+t3*fbji;

			FCJI[it][0]=fxx;
			FCJI[it][1]=fyy;
			FCJI[it][2]=fzz;  

			mxji=n1*mnji+s1*mtji+t1*mbji; 
			myji=n2*mnji+s2*mtji+t2*mbji; 
			mzji=n3*mnji+s3*mtji+t3*mbji; 

			mxij=-n1*mnji+s1*mtij+t1*mbij; 
			myij=-n2*mnji+s2*mtij+t2*mbij; 
			mzij=-n3*mnji+s3*mtij+t3*mbij; 

			if(numth==-1){
			FX[numc1]+=fxx;
			FY[numc1]+=fyy;
			FZ[numc1]+=fzz;
			FX[numc2]+=(-fxx);
			FY[numc2]+=(-fyy);
			FZ[numc2]+=(-fzz);
			MTX[numc1]+=mxji;
			MTX[numc2]+=mxij;
			MTY[numc1]+=myji;
			MTY[numc2]+=myij;
			MTZ[numc1]+=mzji;
			MTZ[numc2]+=mzij;
			}else{
			FXP[numth][numc1]+=fxx;
			FYP[numth][numc1]+=fyy;
			FZP[numth][numc1]+=fzz;
			FXP[numth][numc2]+=(-fxx);
			FYP[numth][numc2]+=(-fyy);
			FZP[numth][numc2]+=(-fzz);
			MTXP[numth][numc1]+=mxji;
			MTXP[numth][numc2]+=mxij;
			MTYP[numth][numc1]+=myji;
			MTYP[numth][numc2]+=myij;
			MTZP[numth][numc1]+=mzji;
			MTZP[numth][numc2]+=mzij;

			}


    
        }	
  



# pragma omp parallel for schedule(dynamic,int(NB_SPH/num_threads)+1) private(it,jt)
for(jt=0;jt<NB_SPH;jt++){
	for(it=0;it<num_threads-1;it++){ 
	FX[jt]+=FXP[it][jt];
	FY[jt]+=FYP[it][jt];
	FZ[jt]+=FZP[it][jt];
	MTX[jt]+=MTXP[it][jt];
	MTY[jt]+=MTYP[it][jt];
	MTZ[jt]+=MTZP[it][jt];
        }
}



}

Epe=Epe1;
Epa=Epa1;

}

void forcep(int NB_SPH, R &fpar,R &fpar2,R dt, int NBCOP, int ** CONTP, R * DCONTP, R ** NCONTP, R * LIST_M,  R * LIST_PM, R ** LIST_PVX, R ** LIST_PVY, R ** LIST_PVZ, R * LIST_R, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R E, R nu, R fric)
{
int numc,nump;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,mxij,myij,mzij;
R fxx,fyy,fzz;
R abspred2,pred1t;
R Eij,iEij,rij,hij;
R dist,rn,rnv,crn,crnv,drn,drnv,mas1,pred1,pred2,pred3,u1,u2,u3,Kn,Kt;	
R vpx,vpy,vpz;	
fpar=0.;
fpar2=0.;

int it,num_threads;
/*
# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCOP/num_threads)) private(it,nump,numc,n1,n2,n3,s1,s2,s3,t1,t2,t3,dist,mas1,vpx,vpy,vpz,u1,u2,u3,rij,hij,iEij,Eij,Kn,Kt,rn,crn,drn,rnv,crnv,drnv,pred1,pred2,pred3,pred1t,abspred2,fxx,fyy,fzz,mxij,myij,mzij) reduction(+:fpar,fpar2,FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH]) 


*/
	for(it=0;it<NBCOP;it++){ 
		
		 // Numeros des candidats 
	     nump=CONTP[it][0];
	     numc=CONTP[it][1];     
	     
	     // normal ext 1->2
	     n1=NCONTP[it][0];
	     n2=NCONTP[it][1];
	     n3=NCONTP[it][2];
	     s1=NCONTP[it][3];
	     s2=NCONTP[it][4];
	     s3=NCONTP[it][5];
	     t1=NCONTP[it][6];
	     t2=NCONTP[it][7];
	     t3=NCONTP[it][8];  		     

	  	 dist=DCONTP[it];     
	     //mas1=LIST_M[numc]*LIST_M[numc]/(LIST_M[numc]+LIST_M[numc]);  
	     mas1=LIST_PM[nump]*1.;
	     
	     vpx=(LIST_PVX[nump][0]+LIST_PVX[nump][1]+LIST_PVX[nump][2]+LIST_PVX[nump][3])/4.;
	     vpy=(LIST_PVY[nump][0]+LIST_PVY[nump][1]+LIST_PVY[nump][2]+LIST_PVY[nump][3])/4.;
	     vpz=(LIST_PVZ[nump][0]+LIST_PVZ[nump][1]+LIST_PVZ[nump][2]+LIST_PVZ[nump][3])/4.;	     	     
	     
        // Contact sph-par (Non-cohésif)                   
	     
         u1=0.;
	     u2=u2-LIST_R[numc]*t1*LIST_WX[numc]-LIST_R[numc]*t2*LIST_WY[numc]-LIST_R[numc]*t3*LIST_WZ[numc];
	     u3=u3+LIST_R[numc]*s1*LIST_WX[numc]+LIST_R[numc]*s2*LIST_WY[numc]+LIST_R[numc]*s3*LIST_WZ[numc];	  
	     
	     u1=u1+n1*(vpx-LIST_VX[numc])+n2*(vpy-LIST_VY[numc])+n3*(vpz-LIST_VZ[numc]);
	     u2=u2+s1*(vpx-LIST_VX[numc])+s2*(vpy-LIST_VY[numc])+s3*(vpz-LIST_VZ[numc]);
	     u3=u3+t1*(vpx-LIST_VX[numc])+t2*(vpy-LIST_VY[numc])+t3*(vpz-LIST_VZ[numc]);	 

	     rij=LIST_R[numc];
	     hij=fabs(dist-LIST_R[numc]);
		 iEij=(1.-nu*nu)/E;		 
		 Eij=1./iEij;	
		 Kn=4.*Eij*sqrt(rij*hij)/3.;
		 //Kn=4.*Eij*rij;   
	//     Kn=3.35e7; 1 PAROI
	//     Kn=1.155e7; 2 parois
	         Kn=1e7;
	   	 Kt=Kn;
	     	     	 	     
	     rn=Kn*(dist-LIST_R[numc]);
	     crn=Kt*u2*dt;
	     drn=Kt*u3*dt;
		  
	     rnv=sqrt(Kn*mas1)*u1;   
	     crnv=sqrt(Kn*mas1)*u2;  
	     drnv=sqrt(Kn*mas1)*u3;  

	     pred1=-rn-rnv;
	     pred1t=pred1;
	     pred2=-crn-crnv;    
	     pred3=-drn-drnv;    

	     abspred2=sqrt(pred2*pred2+pred3*pred3);

	     if(fric==0){pred2=0.;pred3=0.;}

	     if(pred1<0){ //relachement
		 pred1=0.;
		 pred2=0.;	 
		 pred3=0.;	 
		 }
		 else if(abspred2<=fric*pred1){ //adherence
		 pred1=pred1;
		 pred2=pred2;	 
		 pred3=pred3;	 
		 }
		 else if(fric>0){ // Glissement avec frottement
		 pred1=pred1;
		 pred2=-abs(fric*pred1)*u2/abs(u2); 
		 pred3=-abs(fric*pred1)*u3/abs(u3); 
		 }	 
		 else // Glissement sans frottement
		 {
		 pred1=pred1;
		 pred2=0.;
		 pred3=0.;
	     }
		          
		fxx=n1*pred1+s1*pred2+t1*pred3;
		fyy=n2*pred1+s2*pred2+t2*pred3;
		fzz=n3*pred1+s3*pred2+t3*pred3;		  

		mxij = -LIST_R[numc]*t1*pred2+LIST_R[numc]*s1*pred3;
		myij = -LIST_R[numc]*t2*pred2+LIST_R[numc]*s2*pred3;
		mzij = -LIST_R[numc]*t3*pred2+LIST_R[numc]*s3*pred3;	     

		FX[numc]+=(-fxx);
		FY[numc]+=(-fyy);
		FZ[numc]+=(-fzz);			    
			 
		if(nump==0){fpar=fpar+fxx;}	 
		if(nump==1){fpar2=fpar2+fxx;}	 
					 
		MTX[numc]+=mxij;
		MTY[numc]+=myij;	
		MTZ[numc]+=mzij;	
	     	
    }
    

			  
	
}	

	      
void forcecoh(R & Epe, R & Epa,R dt, bool * TYPCO,int NBCO, int NB_SPH, int ** CONT,  R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONT,R * DCONTO, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI,R ** MTJI,R ** MTIJ, R E, R nu, R fric){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij,abspred2;
R mnji,mtji,mbji;
R rij,hij;
R dist,rn,rnv,crn,crnv,drn,drnv,mas1,pred1,pred2,pred3,u1,u2,u3,Kn,Kt;
R Eij,iEij;
Epe=0.;
Epa=0.;

int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)+1) private(it,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,numc1,numc2,vn_ij,vt_ij,vb_ij,dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,abspred2,mnji,mtji,mbji,rij,hij,dist,rn,rnv,crn,crnv,drn,drnv,mas1,pred1,pred2,pred3,u1,u2,u3,Kn,Kt,Eij,iEij) reduction(+:Epa,Epe,FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH]) 
	for(it=0;it<NBCO;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     
		
	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
	     
	  	 dist=DCONT[it];     
	     mas1=LIST_M[numc1]*LIST_M[numc2]/(LIST_M[numc1]+LIST_M[numc2]);      
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	

             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

		fnji=-VALCOH[it][0]*un_ij;
		ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
		fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  

		mnji=-VALCOH[it][5]*(thetani-thetanj);
		mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
		mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
				 
		mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
		mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;

				
// Modèle Kn/Kt
/*
				fnji=-VALCOH[it][0]*un_ij;
				ftji=-VALCOH[it][1]*ut_ij;
				fbji=-VALCOH[it][1]*ub_ij;  

				mnji=0.;
				mtji=0.;
				mbji=0.;
						 
				mtij=0.;  
				mbij=0.;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;  
				
				*/
				

// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  


/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/
		
              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  		
		
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
		 	  

	
	      }
          else{

		fxx=0.;
		fyy=0.;
		fzz=0.;

		mxji=0.;
		myji=0.;
		mzji=0.;

		mxij=0.;
		myij=0.;
		mzij=0.;


        // Contact sph-sph (Non-cohésif)                   
	    
         u1=0.;
	     u2=-LIST_R[numc1]*t1*LIST_WX[numc1]-LIST_R[numc1]*t2*LIST_WY[numc1]-LIST_R[numc1]*t3*LIST_WZ[numc1];
	     u2=u2-LIST_R[numc2]*t1*LIST_WX[numc2]-LIST_R[numc2]*t2*LIST_WY[numc2]-LIST_R[numc2]*t3*LIST_WZ[numc2];
	     u3=LIST_R[numc1]*s1*LIST_WX[numc1]+LIST_R[numc1]*s2*LIST_WY[numc1]+LIST_R[numc1]*s3*LIST_WZ[numc1];
	     u3=u3+LIST_R[numc2]*s1*LIST_WX[numc2]+LIST_R[numc2]*s2*LIST_WY[numc2]+LIST_R[numc2]*s3*LIST_WZ[numc2];	  
	     
	     u1=u1+n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u2=u2+s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u3=u3+t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	 

	     rij=(LIST_R[numc1]*LIST_R[numc2])/(LIST_R[numc1]+LIST_R[numc2]);
	     hij=abs(dist-DCONTO[it]);

		 iEij=2.*(1.-nu*nu)/E;		 
		 Eij=1./iEij;	
	     Kn=4.*Eij*sqrt(rij*hij)/3.;
	     Kn=Kn/100.;
	     Kn=2.5e8;
	     Kt=Kn;
	     	     	 	     
	 //    rn=Kn*(dist-(LIST_R[numc1]+LIST_R[numc2]));
	     rn=Kn*(dist-DCONTO[it]);	     
	     crn=Kt*u2*dt;
	     drn=Kt*u3*dt;
		  
	     rnv=sqrt(Kn*mas1)*u1;   
	     crnv=sqrt(Kn*mas1)*u2;  
	     drnv=sqrt(Kn*mas1)*u3;  

	     pred1=-rn-rnv;
	     pred2=-crn-crnv;    
	     pred3=-drn-drnv;    

	     abspred2=sqrt(pred2*pred2+pred3*pred3);
	       
	     if(fric==0){pred2=0.;pred3=0.;}

	     if(pred1<0){ //relachement
		 pred1=0.;
		 pred2=0.;	 
		 pred3=0.;	 
		 }
		 else if(abspred2<=fric*pred1){ //adherence
		 pred1=pred1;
		 pred2=pred2;	 
		 pred3=pred3;	 
		 }
		 else if(fric>0){ // Glissement avec frottement
		 pred1=pred1;
		 pred2=-abs(fric*pred1)*u2/abs(u2); 
		 pred3=-abs(fric*pred1)*u3/abs(u3); 
		 }	 
		 else // Glissement sans frottement
		 {
		 pred1=pred1;
		 pred2=0.;
		 pred3=0.;
	         }
 
	          
          fxx=n1*pred1+s1*pred2+t1*pred3;
          fyy=n2*pred1+s2*pred2+t2*pred3;
          fzz=n3*pred1+s3*pred2+t3*pred3;
		  
		  mxji = -LIST_R[numc1]*t1*pred2+LIST_R[numc1]*s1*pred3;
		  mxij = -LIST_R[numc2]*t1*pred2+LIST_R[numc2]*s1*pred3;
		  myji = -LIST_R[numc1]*t2*pred2+LIST_R[numc1]*s2*pred3;
		  myij = -LIST_R[numc2]*t2*pred2+LIST_R[numc2]*s2*pred3;
		  mzji = -LIST_R[numc1]*t3*pred2+LIST_R[numc1]*s3*pred3; 
		  mzij = -LIST_R[numc2]*t3*pred2+LIST_R[numc2]*s3*pred3;
			 
	      }   


			FCJI[it][0]=fxx;
			FCJI[it][1]=fyy;
			FCJI[it][2]=fzz; 


		//	  # pragma omp atomic
			  FX[numc1]+=fxx;
		//	  # pragma omp atomic
			  FY[numc1]+=fyy;
		//	  # pragma omp atomic
			  FZ[numc1]+=fzz;
		//	  # pragma omp atomic
			  FX[numc2]+=(-fxx);
		//	  # pragma omp atomic
			  FY[numc2]+=(-fyy);
		//	  # pragma omp atomic
			  FZ[numc2]+=(-fzz);
		//	  # pragma omp atomic
			  MTX[numc1]+=mxji;
		//	  # pragma omp atomic
			  MTX[numc2]+=mxij;
		//	  # pragma omp atomic
			  MTY[numc1]+=myji;
		//	  # pragma omp atomic
			  MTY[numc2]+=myij;	
		//	  # pragma omp atomic
			  MTZ[numc1]+=mzji;
		//	  # pragma omp atomic
			  MTZ[numc2]+=mzij;	
			           
        }	
       

}


void forcecoh2(R Kcon, R & Epe, R & Epa,R dt, bool * TYPCO,int NBCO, int NB_SPH,  int ** CONT,  R * LIST_M, R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONT, R * DCONTO, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI,R ** FOJI,R ** MTJI,R ** MTIJ, R E1, R nu1, R E2, R nu2, R fric, bool * LIST_P){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij,abspred2;
R mnji,mtji,mbji;
R rij,hij;
R dist,rn,rnv,crn,crnv,drn,drnv,mas1,pred1,pred2,pred3,u1,u2,u3,Kn,Kt;
R Eij,iEij;

Epe=0.;
Epa=0.;


int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)+1) private(it,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij, dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,abspred2,mnji,mtji,mbji,rij,hij,dist,rn,rnv,crn,crnv,drn,drnv,mas1,pred1,pred2,pred3,u1,u2,u3,Kn,Kt,Eij,iEij) reduction(+:Epe,Epa,FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH])
	for(it=0;it<NBCO;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     
		
	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
	     
	  	 dist=DCONT[it];     
	     mas1=LIST_M[numc1]*LIST_M[numc2]/(LIST_M[numc1]+LIST_M[numc2]);      
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                 
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
              mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;
		
// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  


/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/	

              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  
              
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
		 	  

	
	      }
          else{
			  
			fxx=0.;
			fyy=0.;
			fzz=0.;

			mxji=0.;
			myji=0.;
			mzji=0.;

			mxij=0.;
			myij=0.;
			mzij=0.;
			  
/*			  
        // Contact sph-sph (Non-cohésif)                   
	    
         u1=0.;
	     u2=-LIST_R[numc1]*t1*LIST_WX[numc1]-LIST_R[numc1]*t2*LIST_WY[numc1]-LIST_R[numc1]*t3*LIST_WZ[numc1];
	     u2=u2-LIST_R[numc2]*t1*LIST_WX[numc2]-LIST_R[numc2]*t2*LIST_WY[numc2]-LIST_R[numc2]*t3*LIST_WZ[numc2];
	     u3=LIST_R[numc1]*s1*LIST_WX[numc1]+LIST_R[numc1]*s2*LIST_WY[numc1]+LIST_R[numc1]*s3*LIST_WZ[numc1];
	     u3=u3+LIST_R[numc2]*s1*LIST_WX[numc2]+LIST_R[numc2]*s2*LIST_WY[numc2]+LIST_R[numc2]*s3*LIST_WZ[numc2];	  
	     
	     u1=u1+n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u2=u2+s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
	     u3=u3+t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	 

	     rij=(LIST_R[numc1]*LIST_R[numc2])/(LIST_R[numc1]+LIST_R[numc2]);
	     hij=abs(dist-DCONTO[it]);
	     
		if((LIST_P[numc1]==0)&&(LIST_P[numc2]==0)){				  
		iEij=(1.-nu1*nu1)/E1+(1.-nu1*nu1)/E1;		  
		}
		if(((LIST_P[numc1]==1)&&(LIST_P[numc2]==0))||((LIST_P[numc1]==0)&&(LIST_P[numc2]==1))) {	
		iEij=(1.-nu1*nu1)/E1+(1.-nu2*nu2)/E2;	
		}
		if((LIST_P[numc1]==1)&&(LIST_P[numc2]==1)){				  
		iEij=(1.-nu2*nu2)/E2+(1.-nu2*nu2)/E2;	     
		}
			 
		 Eij=1./iEij;	
		 
	     Kn=4.*Eij*sqrt(rij*hij)/3.;
	     Kn=Kn/100.;
	     Kn=Kcon;
	     Kt=Kn;
	     	     	 	     
	 //    rn=Kn*(dist-(LIST_R[numc1]+LIST_R[numc2]));
	     rn=Kn*(dist-DCONTO[it]);	 
	     crn=Kt*u2*dt;
	     drn=Kt*u3*dt;
		  
	     rnv=sqrt(Kn*mas1)*u1;   
	     crnv=sqrt(Kn*mas1)*u2;  
	     drnv=sqrt(Kn*mas1)*u3;  

	     pred1=-rn-rnv;
	     pred2=-crn-crnv;    
	     pred3=-drn-drnv;    

	     abspred2=sqrt(pred2*pred2+pred3*pred3);

	     if(fric==0){pred2=0.;pred3=0.;}

	     if(pred1<0){ //relachement
		 pred1=0.;
		 pred2=0.;	 
		 pred3=0.;	 
		 }
		 else if(abspred2<=fric*pred1){ //adherence
		 pred1=pred1;
		 pred2=pred2;	 
		 pred3=pred3;	 
		 }
		 else if(fric>0){ // Glissement avec frottement
		 pred1=pred1;
		 pred2=-abs(fric*pred1)*u2/abs(u2); 
		 pred3=-abs(fric*pred1)*u3/abs(u3); 
		 }	 
		 else // Glissement sans frottement
		 {
		 pred1=pred1;
		 pred2=0.;
		 pred3=0.;
	     }
		          
          fxx=n1*pred1+s1*pred2+t1*pred3;
          fyy=n2*pred1+s2*pred2+t2*pred3;
          fzz=n3*pred1+s3*pred2+t3*pred3;
		  
		  mxji = -LIST_R[numc1]*t1*pred2+LIST_R[numc1]*s1*pred3;
		  mxij = -LIST_R[numc2]*t1*pred2+LIST_R[numc2]*s1*pred3;
		  myji = -LIST_R[numc1]*t2*pred2+LIST_R[numc1]*s2*pred3;
		  myij = -LIST_R[numc2]*t2*pred2+LIST_R[numc2]*s2*pred3;
		  mzji = -LIST_R[numc1]*t3*pred2+LIST_R[numc1]*s3*pred3; 
		  mzij = -LIST_R[numc2]*t3*pred2+LIST_R[numc2]*s3*pred3;
			  */
	      }   



              FCJI[it][0]=fxx;
              FCJI[it][1]=fyy;
              FCJI[it][2]=fzz; 

		//	  # pragma omp atomic
			  FX[numc1]+=fxx;
		//	  # pragma omp atomic
			  FY[numc1]+=fyy;
		//	  # pragma omp atomic
			  FZ[numc1]+=fzz;
		//	  # pragma omp atomic
			  FX[numc2]+=(-fxx);
		//	  # pragma omp atomic
			  FY[numc2]+=(-fyy);
		//	  # pragma omp atomic
			  FZ[numc2]+=(-fzz);
		//	  # pragma omp atomic
			  MTX[numc1]+=mxji;
		//	  # pragma omp atomic
			  MTX[numc2]+=mxij;
		//	  # pragma omp atomic
			  MTY[numc1]+=myji;
		//	  # pragma omp atomic
			  MTY[numc2]+=myij;	
		//	  # pragma omp atomic
			  MTZ[numc1]+=mzji;
		//	  # pragma omp atomic
			  MTZ[numc2]+=mzij;
			           
        }	
  
 }

void forcecoh0(R & Epe, R & Epa,bool * TYPCO,int NBCO,int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij;
R mnji,mtji,mbji;
Epe=0.;
Epa=0.;

int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)+1) private(it,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij,dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,mnji,mtji,mbji) reduction(+:Epe,Epa,FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH])
	for(it=0;it<NBCO;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     
		
	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	  
	      
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	

             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
              mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;
			
// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  
		/*	  
			  mnji=mnji-(VALAMO[it][5]*(LIST_WX[numc2]-LIST_WX[numc1]));
			  mtji=mtji-(VALAMO[it][3]*LIST_WY[numc1]+VALAMO[it][4]*LIST_WY[numc2]);
			  mbji=mbji-(VALAMO[it][3]*LIST_WZ[numc1]+VALAMO[it][4]*LIST_WZ[numc2]);
			  
			  mtij=mtij-(VALAMO[it][3]*LIST_WY[numc2]+VALAMO[it][4]*LIST_WY[numc1]);
			  mbij=mbij-(VALAMO[it][3]*LIST_WZ[numc2]+VALAMO[it][4]*LIST_WZ[numc1]);

/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/
		
              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  		
		
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
		 	  

	
	      }
          else{
	
			  fxx=0.;
			  fyy=0.;
			  fzz=0.;
	
	          mxji=0.;
	          myji=0.;
	          mzji=0.;

	          mxij=0.;
	          myij=0.; 
	          mzij=0.;
		 
	      }   


              FCJI[it][0]=fxx;
              FCJI[it][1]=fyy;
              FCJI[it][2]=fzz; 

			  FX[numc1]+=fxx;
			  FY[numc1]+=fyy;
			  FZ[numc1]+=fzz;
			  FX[numc2]+=(-fxx);
			  FY[numc2]+=(-fyy);
			  FZ[numc2]+=(-fzz);
			  MTX[numc1]+=mxji;
			  MTX[numc2]+=mxij;
			  MTY[numc1]+=myji;
			  MTY[numc2]+=myij;	
			  MTZ[numc1]+=mzji;
			  MTZ[numc2]+=mzij;
			           
        }	
       

}


void forcecoh_qs_incr(int ite, R dt, R & Epe, R & Epa,int NBCO, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_XA, R * LIST_YA, R * LIST_ZA,R * LIST_VXA, R * LIST_VYA, R * LIST_VZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ,R * FIX, R * FIY, R * FIZ, R * MTIX, R * MTIY, R * MTIZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,n1_t,n2_t,n3_t,s1_t,s2_t,s3_t,t1_t,t2_t,t3_t,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij_t,ut_ij_t,ub_ij_t,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetaxi_t,thetaxj_t,thetayi_t,thetayj_t,thetazi_t,thetazj_t;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R thetani_t,thetanj_t,thetati_t,thetatj_t,thetabi_t,thetabj_t,dcont;
R vn_ij,vt_ij,vb_ij,vn,vb;
R dxij,dyij,dzij; 
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij;
R mnji,mtji,mbji;
R b1,b2,b3;
R cx,cy,cz,sx,sy,sz;
R n1d,n2d,n3d,s1d,s2d,s3d,t1d,t2d,t3d;

R fnji_tot,ftji_tot,fbji_tot,mnji_tot,mtji_tot,mbji_tot,mtij_tot,mbij_tot;

R fnjia,ftjia,fbjia;
R mtija,mbija;
R mnjia,mtjia,mbjia;
R fxxa,fyya,fzza;
R mxija,mxjia,myija,myjia,mzija,mzjia;

int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)+1) private(fxxa,fyya,fzza,mxija,mxjia,myija,myjia,mzija,mzjia,fnjia,ftjia,fbjia,mtija,mbija,mnjia,mtjia,mbjia,it,numc1,numc2,n1d,n2d,n3d,s1d,s2d,s3d,t1d,t2d,t3d,n1,n2,n3,s1,s2,s3,t1,t2,t3,n1_t,n2_t,n3_t,s1_t,s2_t,s3_t,t1_t,t2_t,t3_t,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij,un_ij,ut_ij,ub_ij,un_ij_t,ut_ij_t,ub_ij_t,fxx,fyy,fzz,mtij,mbij,mnji,mtji,mbji,thetaxi_t,thetaxj_t,thetayi_t,thetayj_t,thetazi_t,thetazj_t,thetani_t,thetanj_t,thetati_t,thetatj_t,thetabi_t,thetabj_t,dcont,dxij,dyij,dzij,vn,vb,b1,b2,b3,cx,cy,cz,sx,sy,sz) reduction(+:Epa,FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH],FIX[0:NB_SPH],FIY[0:NB_SPH],FIZ[0:NB_SPH],MTIX[0:NB_SPH],MTIY[0:NB_SPH],MTIZ[0:NB_SPH]) reduction(-:Epe)
	for(it=0;it<NBCO;it++){ 
		
		// Numeros des candidats 
	        numc1=CONT[it][0];
	        numc2=CONT[it][1];   
	
//////////////////////////////////////// normales obtenues à partir de la vitesse //////////////////////////////////////

		// normal au pas t précédent
		n1_t=NCONT[it][0];
		n2_t=NCONT[it][1];
		n3_t=NCONT[it][2];
		s1_t=NCONT[it][3];
		s2_t=NCONT[it][4];
		s3_t=NCONT[it][5];
		t1_t=NCONT[it][6];
		t2_t=NCONT[it][7];
		t3_t=NCONT[it][8]; 

		// normal au pas actuel

		dcont=sqrt((LIST_XA[numc1]-LIST_XA[numc2])*(LIST_XA[numc1]-LIST_XA[numc2])+(LIST_YA[numc1]-LIST_YA[numc2])*(LIST_YA[numc1]-LIST_YA[numc2])+(LIST_ZA[numc1]-LIST_ZA[numc2])*(LIST_ZA[numc1]-LIST_ZA[numc2]));

		n1=(LIST_XA[numc1]-LIST_XA[numc2])/(dcont);
		n2=(LIST_YA[numc1]-LIST_YA[numc2])/(dcont);
		n3=(LIST_ZA[numc1]-LIST_ZA[numc2])/(dcont);

		dxij=(LIST_VXA[numc1]-LIST_VXA[numc2]);        	           
		dyij=(LIST_VYA[numc1]-LIST_VYA[numc2]);
		dzij=(LIST_VZA[numc1]-LIST_VZA[numc2]);

		vn=n1*dxij+n2*dyij+n3*dzij;
		vb=sqrt(dxij*dxij+dyij*dyij+dzij*dzij-vn*vn);

		if(fabs(vb)>1e-10){

		b1=(dxij-vn*n1)/vb;
		b2=(dyij-vn*n2)/vb;
		b3=(dzij-vn*n3)/vb;

		s1=sqrt(2)/2.*(b1*(n1*n1+1.)+b2*(n1*n2+n3)+b3*(n1*n3-n2));
		s2=sqrt(2)/2.*(b1*(n1*n2-n3)+b2*(n2*n2+1.)+b3*(n2*n3+n1));
		s3=sqrt(2)/2.*(b1*(n1*n3+n2)+b2*(n2*n3-n1)+b3*(n3*n3+1.));

		t1=n2*s3-n3*s2;
		t2=n3*s1-n1*s3;
		t3=n1*s2-n2*s1; 
		
		}else{

		s1=s1_t;
		s2=s2_t;
		s3=s3_t;
		t1=t1_t;
		t2=t2_t;
		t3=t3_t;

		}


		NCONT[it][0]=n1;
		NCONT[it][1]=n2;
		NCONT[it][2]=n3;
		NCONT[it][3]=s1;
		NCONT[it][4]=s2;
		NCONT[it][5]=s3;
		NCONT[it][6]=t1;
		NCONT[it][7]=t2;
		NCONT[it][8]=t3; 

/// normales egales aux normales initiales ///

/*
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8]; 

	     n1_t=NCONT[it][0];
	     n2_t=NCONT[it][1];
	     n3_t=NCONT[it][2];
	     s1_t=NCONT[it][3];
	     s2_t=NCONT[it][4];
	     s3_t=NCONT[it][5];
	     t1_t=NCONT[it][6];
	     t2_t=NCONT[it][7];
	     t3_t=NCONT[it][8]; 
*/

////////////////////////////////////////////////


          // Déplacements globaux dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements globaux dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
   
         // Déplacements relatifs dans le repère global
             dxij=(LIST_X[numc1]-LIST_X[numc2])-(LIST_XA[numc1]-LIST_XA[numc2]); 
             dyij=(LIST_Y[numc1]-LIST_Y[numc2])-(LIST_YA[numc1]-LIST_YA[numc2]);
             dzij=(LIST_Z[numc1]-LIST_Z[numc2])-(LIST_ZA[numc1]-LIST_ZA[numc2]);


         // Déplacements relatifs dans le repère local             
             un_ij_t=n1*dxij+n2*dyij+n3*dzij;
             ut_ij_t=s1*dxij+s2*dyij+s3*dzij;
             ub_ij_t=t1*dxij+t2*dyij+t3*dzij;       
  
           // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];             
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
             thetaxi_t=LIST_TXA[numc1];
             thetaxj_t=LIST_TXA[numc2];
             thetayi_t=LIST_TYA[numc1];             
             thetayj_t=LIST_TYA[numc2];
             thetazi_t=LIST_TZA[numc1];
             thetazj_t=LIST_TZA[numc2];     
                        
             thetani=n1*(thetaxi-thetaxi_t)+n2*(thetayi-thetayi_t)+n3*(thetazi-thetazi_t);   
             thetati=s1*(thetaxi-thetaxi_t)+s2*(thetayi-thetayi_t)+s3*(thetazi-thetazi_t);     
             thetabi=t1*(thetaxi-thetaxi_t)+t2*(thetayi-thetayi_t)+t3*(thetazi-thetazi_t);  
             thetanj=n1*(thetaxj-thetaxj_t)+n2*(thetayj-thetayj_t)+n3*(thetazj-thetazj_t);   
             thetatj=s1*(thetaxj-thetaxj_t)+s2*(thetayj-thetayj_t)+s3*(thetazj-thetazj_t);     
             thetabj=t1*(thetaxj-thetaxj_t)+t2*(thetayj-thetayj_t)+t3*(thetazj-thetazj_t);    
        
          // Vitesses relatives dans le repère local
		vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
		vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
		vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                    
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 



			fnji=-VALCOH[it][0]*un_ij_t;
			ftji=-VALCOH[it][1]*ut_ij_t-VALCOH[it][2]*((thetabi)+(thetabj));
			fbji=-VALCOH[it][1]*ub_ij_t+VALCOH[it][2]*((thetati)+(thetatj));  

			mnji=-VALCOH[it][5]*((thetani)-(thetanj));
			mtji=VALCOH[it][2]*ub_ij_t-VALCOH[it][3]*(thetati)-VALCOH[it][4]*(thetatj);
			mbji=-VALCOH[it][2]*ut_ij_t-VALCOH[it][3]*(thetabi)-VALCOH[it][4]*(thetabj);
					 
			mtij=VALCOH[it][2]*ub_ij_t-VALCOH[it][3]*(thetatj)-VALCOH[it][4]*(thetati);  
			mbij=-VALCOH[it][2]*ut_ij_t-VALCOH[it][3]*(thetabj)-VALCOH[it][4]*(thetabi);

			fnji_tot=-VALCOH[it][0]*un_ij;
			ftji_tot=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*((thetabi)+(thetabj));
			fbji_tot=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*((thetati)+(thetatj));  

			mnji_tot=-VALCOH[it][5]*((thetani)-(thetanj));
			mtji_tot=VALCOH[it][2]*ub_ij-VALCOH[it][3]*(thetati)-VALCOH[it][4]*(thetatj);
			mbji_tot=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*(thetabi)-VALCOH[it][4]*(thetabj);
					 
			mtij_tot=VALCOH[it][2]*ub_ij-VALCOH[it][3]*(thetatj)-VALCOH[it][4]*(thetati);  
			mbij_tot=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*(thetabj)-VALCOH[it][4]*(thetabi);

			Epe-=fnji_tot*un_ij_t;
			Epe-=ftji_tot*ut_ij_t;                      
			Epe-=fbji_tot*ub_ij_t;      

			Epe-=mnji_tot*(thetani);            
			Epe-=mtji_tot*(thetati);                   
			Epe-=mbji_tot*(thetabi);                     

			Epe-=(-mnji_tot)*(thetanj);                      
			Epe-=mtij_tot*(thetatj);                   
			Epe-=mbij_tot*(thetabj);      


// Amortissement 

			fnjia=-(VALAMO[it][0]*vn_ij);
			ftjia=-(VALAMO[it][1]*vt_ij);
			fbjia=-(VALAMO[it][1]*vb_ij);	
/*	  
			Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
			Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
			Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 	
*/	
						
              
				FOJI[it][0]=fnji;
				FOJI[it][1]=ftji;
				FOJI[it][2]=fbji;

				MTJI[it][0]=mnji;
				MTJI[it][1]=mtji;
				MTJI[it][2]=mbji;                                      

				MTIJ[it][0]=-mnji;
				MTIJ[it][1]=mtij;
				MTIJ[it][2]=mbij;  

				fxx=n1*fnji+s1*ftji+t1*fbji;
				fyy=n2*fnji+s2*ftji+t2*fbji;
				fzz=n3*fnji+s3*ftji+t3*fbji;

				fxxa=n1*fnjia+s1*ftjia+t1*fbjia;
				fyya=n2*fnjia+s2*ftjia+t2*fbjia;
				fzza=n3*fnjia+s3*ftjia+t3*fbjia;

                                if(ite!=1){
				FCJI[it][0]+=fxx;
				FCJI[it][1]+=fyy;
				FCJI[it][2]+=fzz;  
				}else{
				FCJI[it][0]=fxx;
				FCJI[it][1]=fyy;
				FCJI[it][2]=fzz;  
				}


				mxji=n1*mnji+s1*mtji+t1*mbji; 
				myji=n2*mnji+s2*mtji+t2*mbji; 
				mzji=n3*mnji+s3*mtji+t3*mbji; 

				mxij=-n1*mnji+s1*mtij+t1*mbij; 
				myij=-n2*mnji+s2*mtij+t2*mbij; 
				mzij=-n3*mnji+s3*mtij+t3*mbij; 

				mxjia=n1*mnjia+s1*mtjia+t1*mbjia; 
				myjia=n2*mnjia+s2*mtjia+t2*mbjia; 
				mzjia=n3*mnjia+s3*mtjia+t3*mbjia; 

				mxija=-n1*mnjia+s1*mtija+t1*mbija; 
				myija=-n2*mnjia+s2*mtija+t2*mbija; 
				mzija=-n3*mnjia+s3*mtija+t3*mbija; 

				FIX[numc1]=FIX[numc1]+fxx;			 
				FIY[numc1]=FIY[numc1]+fyy;			  
				FIZ[numc1]=FIZ[numc1]+fzz;

				FIX[numc2]=FIX[numc2]-fxx;
				FIY[numc2]=FIY[numc2]-fyy;
				FIZ[numc2]=FIZ[numc2]-fzz;			    
					 
				MTIX[numc1]=MTIX[numc1]+mxji;
				MTIX[numc2]=MTIX[numc2]+mxij;

				MTIY[numc1]=MTIY[numc1]+myji;
				MTIY[numc2]=MTIY[numc2]+myij;			  

				MTIZ[numc1]=MTIZ[numc1]+mzji;
				MTIZ[numc2]=MTIZ[numc2]+mzij;	

				FX[numc1]=FX[numc1]+fxxa;			 
				FY[numc1]=FY[numc1]+fyya;			  
				FZ[numc1]=FZ[numc1]+fzza;

				FX[numc2]=FX[numc2]-fxxa;
				FY[numc2]=FY[numc2]-fyya;
				FZ[numc2]=FZ[numc2]-fzza;			    
				 
				MTX[numc1]=MTX[numc1]+mxjia;
				MTX[numc2]=MTX[numc2]+mxija;

				MTY[numc1]=MTY[numc1]+myjia;
				MTY[numc2]=MTY[numc2]+myija;			  

				MTZ[numc1]=MTZ[numc1]+mzjia;
				MTZ[numc2]=MTZ[numc2]+mzija;		
    
        }	
       

}


void forcecoh0_incr(int ite, R dt, R & Epe, R & Epa, bool * TYPCO,int NBCO, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z,R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_XA, R * LIST_YA, R * LIST_ZA,R * LIST_VXA, R * LIST_VYA, R * LIST_VZA, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_TXA, R * LIST_TYA, R * LIST_TZA,R * LIST_WX, R * LIST_WY, R * LIST_WZ, R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT, R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ,R * FIX, R * FIY, R * FIZ, R * MTIX, R * MTIY, R * MTIZ, R ** VALCOH, R ** VALAMO, R ** FCJI, R ** FOJI, R ** MTJI, R ** MTIJ){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,n1_t,n2_t,n3_t,s1_t,s2_t,s3_t,t1_t,t2_t,t3_t,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij_t,ut_ij_t,ub_ij_t,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetaxi_t,thetaxj_t,thetayi_t,thetayj_t,thetazi_t,thetazj_t;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R thetani_t,thetanj_t,thetati_t,thetatj_t,thetabi_t,thetabj_t,dcont;
R vn_ij,vt_ij,vb_ij,vn,vb;
R dxij,dyij,dzij; 
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij;
R mnji,mtji,mbji;
R b1,b2,b3;
R cx,cy,cz,sx,sy,sz;
R n1d,n2d,n3d,s1d,s2d,s3d,t1d,t2d,t3d;

R fnji_tot,ftji_tot,fbji_tot,mnji_tot,mtji_tot,mbji_tot,mtij_tot,mbij_tot;

R fnjia,ftjia,fbjia;
R mtija,mbija;
R mnjia,mtjia,mbjia;
R fxxa,fyya,fzza;
R mxija,mxjia,myija,myjia,mzija,mzjia;

int it,num_threads;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

# pragma omp parallel for schedule(dynamic,int(NBCO/num_threads)+1) private(fnji_tot,ftji_tot,fbji_tot,mnji_tot,mtji_tot,mbji_tot,mtij_tot,mbij_tot,fxxa,fyya,fzza,mxija,mxjia,myija,myjia,mzija,mzjia,fnjia,ftjia,fbjia,mtija,mbija,mnjia,mtjia,mbjia,it,numc1,numc2,n1d,n2d,n3d,s1d,s2d,s3d,t1d,t2d,t3d,n1,n2,n3,s1,s2,s3,t1,t2,t3,n1_t,n2_t,n3_t,s1_t,s2_t,s3_t,t1_t,t2_t,t3_t,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij,un_ij,ut_ij,ub_ij,un_ij_t,ut_ij_t,ub_ij_t,fxx,fyy,fzz,mtij,mbij,mnji,mtji,mbji,thetaxi_t,thetaxj_t,thetayi_t,thetayj_t,thetazi_t,thetazj_t,thetani_t,thetanj_t,thetati_t,thetatj_t,thetabi_t,thetabj_t,dcont,dxij,dyij,dzij,dxoij,dyoij,dzoij,vn,vb,b1,b2,b3,cx,cy,cz,sx,sy,sz) reduction(+:Epa,FX[0:NB_SPH],FY[0:NB_SPH],FZ[0:NB_SPH],MTX[0:NB_SPH],MTY[0:NB_SPH],MTZ[0:NB_SPH],FIX[0:NB_SPH],FIY[0:NB_SPH],FIZ[0:NB_SPH],MTIX[0:NB_SPH],MTIY[0:NB_SPH],MTIZ[0:NB_SPH]) reduction(-:Epe)
	for(it=0;it<NBCO;it++){ 
		
		// Numeros des candidats 
	        numc1=CONT[it][0];
	        numc2=CONT[it][1];   
	
//////////////////////////////////////// normales obtenues à partir de la vitesse //////////////////////////////////////

if(TYPCO[it]==1){      

		// normal au pas t précédent
		n1_t=NCONT[it][0];
		n2_t=NCONT[it][1];
		n3_t=NCONT[it][2];
		s1_t=NCONT[it][3];
		s2_t=NCONT[it][4];
		s3_t=NCONT[it][5];
		t1_t=NCONT[it][6];
		t2_t=NCONT[it][7];
		t3_t=NCONT[it][8]; 

		// normal au pas actuel

		dcont=sqrt((LIST_XA[numc1]-LIST_XA[numc2])*(LIST_XA[numc1]-LIST_XA[numc2])+(LIST_YA[numc1]-LIST_YA[numc2])*(LIST_YA[numc1]-LIST_YA[numc2])+(LIST_ZA[numc1]-LIST_ZA[numc2])*(LIST_ZA[numc1]-LIST_ZA[numc2]));

		n1=(LIST_XA[numc1]-LIST_XA[numc2])/(dcont);
		n2=(LIST_YA[numc1]-LIST_YA[numc2])/(dcont);
		n3=(LIST_ZA[numc1]-LIST_ZA[numc2])/(dcont);

		dxij=(LIST_VXA[numc1]-LIST_VXA[numc2]);        	           
		dyij=(LIST_VYA[numc1]-LIST_VYA[numc2]);
		dzij=(LIST_VZA[numc1]-LIST_VZA[numc2]);

		vn=n1*dxij+n2*dyij+n3*dzij;
		vb=sqrt(dxij*dxij+dyij*dyij+dzij*dzij-vn*vn);

		if(fabs(vb)>1e-10){

		b1=(dxij-vn*n1)/vb;
		b2=(dyij-vn*n2)/vb;
		b3=(dzij-vn*n3)/vb;

		s1=sqrt(2)/2.*(b1*(n1*n1+1.)+b2*(n1*n2+n3)+b3*(n1*n3-n2));
		s2=sqrt(2)/2.*(b1*(n1*n2-n3)+b2*(n2*n2+1.)+b3*(n2*n3+n1));
		s3=sqrt(2)/2.*(b1*(n1*n3+n2)+b2*(n2*n3-n1)+b3*(n3*n3+1.));

		t1=n2*s3-n3*s2;
		t2=n3*s1-n1*s3;
		t3=n1*s2-n2*s1; 
		
		}else{

		s1=s1_t;
		s2=s2_t;
		s3=s3_t;
		t1=t1_t;
		t2=t2_t;
		t3=t3_t;

		}


		NCONT[it][0]=n1;
		NCONT[it][1]=n2;
		NCONT[it][2]=n3;
		NCONT[it][3]=s1;
		NCONT[it][4]=s2;
		NCONT[it][5]=s3;
		NCONT[it][6]=t1;
		NCONT[it][7]=t2;
		NCONT[it][8]=t3; 

/// normales egales aux normales initiales ///

/*
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8]; 

	     n1_t=NCONT[it][0];
	     n2_t=NCONT[it][1];
	     n3_t=NCONT[it][2];
	     s1_t=NCONT[it][3];
	     s2_t=NCONT[it][4];
	     s3_t=NCONT[it][5];
	     t1_t=NCONT[it][6];
	     t2_t=NCONT[it][7];
	     t3_t=NCONT[it][8]; 
*/

////////////////////////////////////////////////

          // Déplacements globaux dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements globaux dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
         
         // Déplacements relatifs dans le repère global

             dxij=(LIST_X[numc1]-LIST_X[numc2])-(LIST_XA[numc1]-LIST_XA[numc2]); 
             dyij=(LIST_Y[numc1]-LIST_Y[numc2])-(LIST_YA[numc1]-LIST_YA[numc2]);
             dzij=(LIST_Z[numc1]-LIST_Z[numc2])-(LIST_ZA[numc1]-LIST_ZA[numc2]);


         // Déplacements relatifs dans le repère local
             
             un_ij_t=n1*dxij+n2*dyij+n3*dzij;
             ut_ij_t=s1*dxij+s2*dyij+s3*dzij;
             ub_ij_t=t1*dxij+t2*dyij+t3*dzij;       
  
           // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];             
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
             thetaxi_t=LIST_TXA[numc1];
             thetaxj_t=LIST_TXA[numc2];
             thetayi_t=LIST_TYA[numc1];             
             thetayj_t=LIST_TYA[numc2];
             thetazi_t=LIST_TZA[numc1];
             thetazj_t=LIST_TZA[numc2];     
                        
             thetani=n1*(thetaxi-thetaxi_t)+n2*(thetayi-thetayi_t)+n3*(thetazi-thetazi_t);   
             thetati=s1*(thetaxi-thetaxi_t)+s2*(thetayi-thetayi_t)+s3*(thetazi-thetazi_t);     
             thetabi=t1*(thetaxi-thetaxi_t)+t2*(thetayi-thetayi_t)+t3*(thetazi-thetazi_t);  
             thetanj=n1*(thetaxj-thetaxj_t)+n2*(thetayj-thetayj_t)+n3*(thetazj-thetazj_t);   
             thetatj=s1*(thetaxj-thetaxj_t)+s2*(thetayj-thetayj_t)+s3*(thetazj-thetazj_t);     
             thetabj=t1*(thetaxj-thetaxj_t)+t2*(thetayj-thetayj_t)+t3*(thetazj-thetazj_t);    
        
          // Vitesses relatives dans le repère local
		vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
		vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
		vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                    
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 



			fnji=-VALCOH[it][0]*un_ij_t;
			ftji=-VALCOH[it][1]*ut_ij_t-VALCOH[it][2]*((thetabi)+(thetabj));
			fbji=-VALCOH[it][1]*ub_ij_t+VALCOH[it][2]*((thetati)+(thetatj));  

			mnji=-VALCOH[it][5]*((thetani)-(thetanj));
			mtji=VALCOH[it][2]*ub_ij_t-VALCOH[it][3]*(thetati)-VALCOH[it][4]*(thetatj);
			mbji=-VALCOH[it][2]*ut_ij_t-VALCOH[it][3]*(thetabi)-VALCOH[it][4]*(thetabj);
					 
			mtij=VALCOH[it][2]*ub_ij_t-VALCOH[it][3]*(thetatj)-VALCOH[it][4]*(thetati);  
			mbij=-VALCOH[it][2]*ut_ij_t-VALCOH[it][3]*(thetabj)-VALCOH[it][4]*(thetabi);

			fnji_tot=-VALCOH[it][0]*un_ij;
			ftji_tot=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*((thetabi)+(thetabj));
			fbji_tot=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*((thetati)+(thetatj));  

			mnji_tot=-VALCOH[it][5]*((thetani)-(thetanj));
			mtji_tot=VALCOH[it][2]*ub_ij-VALCOH[it][3]*(thetati)-VALCOH[it][4]*(thetatj);
			mbji_tot=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*(thetabi)-VALCOH[it][4]*(thetabj);
					 
			mtij_tot=VALCOH[it][2]*ub_ij-VALCOH[it][3]*(thetatj)-VALCOH[it][4]*(thetati);  
			mbij_tot=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*(thetabj)-VALCOH[it][4]*(thetabi);

			Epe-=fnji_tot*un_ij_t;
			Epe-=ftji_tot*ut_ij_t;                      
			Epe-=fbji_tot*ub_ij_t;      

			Epe-=mnji_tot*(thetani);            
			Epe-=mtji_tot*(thetati);                   
			Epe-=-mbji_tot*(thetabi);                     

			Epe-=(-mnji_tot)*(thetanj);                      
			Epe-=mtij_tot*(thetatj);                   
			Epe-=-mbij_tot*(thetabj);   


     

// Amortissement 

			fnjia=-(VALAMO[it][0]*vn_ij);
			ftjia=-(VALAMO[it][1]*vt_ij);
			fbjia=-(VALAMO[it][1]*vb_ij);	
/*	  
	
			Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
			Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
			Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 	
*/	
						
              
				FOJI[it][0]=fnji;
				FOJI[it][1]=ftji;
				FOJI[it][2]=fbji;

				MTJI[it][0]=mnji;
				MTJI[it][1]=mtji;
				MTJI[it][2]=mbji;                                      

				MTIJ[it][0]=-mnji;
				MTIJ[it][1]=mtij;
				MTIJ[it][2]=mbij;  

				fxx=n1*fnji+s1*ftji+t1*fbji;
				fyy=n2*fnji+s2*ftji+t2*fbji;
				fzz=n3*fnji+s3*ftji+t3*fbji;

				fxxa=n1*fnjia+s1*ftjia+t1*fbjia;
				fyya=n2*fnjia+s2*ftjia+t2*fbjia;
				fzza=n3*fnjia+s3*ftjia+t3*fbjia;

				mxji=n1*mnji+s1*mtji+t1*mbji; 
				myji=n2*mnji+s2*mtji+t2*mbji; 
				mzji=n3*mnji+s3*mtji+t3*mbji; 

				mxij=-n1*mnji+s1*mtij+t1*mbij; 
				myij=-n2*mnji+s2*mtij+t2*mbij; 
				mzij=-n3*mnji+s3*mtij+t3*mbij; 

				mxjia=n1*mnjia+s1*mtjia+t1*mbjia; 
				myjia=n2*mnjia+s2*mtjia+t2*mbjia; 
				mzjia=n3*mnjia+s3*mtjia+t3*mbjia; 

				mxija=-n1*mnjia+s1*mtija+t1*mbija; 
				myija=-n2*mnjia+s2*mtija+t2*mbija; 
				mzija=-n3*mnjia+s3*mtija+t3*mbija; 


	      }
          else{
	
			  fxx=0.;
			  fyy=0.;
			  fzz=0.;
	
	          mxji=0.;
	          myji=0.;
	          mzji=0.;

	          mxij=0.;
	          myij=0.; 
	          mzij=0.;

		  fxxa=0.;
		  fyya=0.;
		  fzza=0.;
	
	          mxjia=0.;
	          myjia=0.;
	          mzjia=0.;

	          mxija=0.;
	          myija=0.; 
	          mzija=0.;		 
	      }   

                                if(ite!=1){
				FCJI[it][0]+=fxx;
				FCJI[it][1]+=fyy;
				FCJI[it][2]+=fzz;  
				}else{
				FCJI[it][0]=fxx;
				FCJI[it][1]=fyy;
				FCJI[it][2]=fzz;  
				}


				FIX[numc1]=FIX[numc1]+fxx;			 
				FIY[numc1]=FIY[numc1]+fyy;			  
				FIZ[numc1]=FIZ[numc1]+fzz;

				FIX[numc2]=FIX[numc2]-fxx;
				FIY[numc2]=FIY[numc2]-fyy;
				FIZ[numc2]=FIZ[numc2]-fzz;			    
					 
				MTIX[numc1]=MTIX[numc1]+mxji;
				MTIX[numc2]=MTIX[numc2]+mxij;

				MTIY[numc1]=MTIY[numc1]+myji;
				MTIY[numc2]=MTIY[numc2]+myij;			  

				MTIZ[numc1]=MTIZ[numc1]+mzji;
				MTIZ[numc2]=MTIZ[numc2]+mzij;	

				FX[numc1]=FX[numc1]+fxxa;			 
				FY[numc1]=FY[numc1]+fyya;			  
				FZ[numc1]=FZ[numc1]+fzza;

				FX[numc2]=FX[numc2]-fxxa;
				FY[numc2]=FY[numc2]-fyya;
				FZ[numc2]=FZ[numc2]-fzza;			    
				 
				MTX[numc1]=MTX[numc1]+mxjia;
				MTX[numc2]=MTX[numc2]+mxija;

				MTY[numc1]=MTY[numc1]+myjia;
				MTY[numc2]=MTY[numc2]+myija;			  

				MTZ[numc1]=MTZ[numc1]+mzjia;
				MTZ[numc2]=MTZ[numc2]+mzija;		
    
        }	
       

}


void forcecoh2_qs_int_mI(int & nint,int & nsoft,R Cs, R & Epe, R & Epa,R dt, bool * TYPCO,int & nbco, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ,  R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI,R ** FOJI,R ** MTJI,R ** MTIJ, bool * LIST_P, int * LIST_B, R * vs, R * LIST_IND,R Pi,R Emu1, R rmu1, R Emu2, R rmu2, R amort, R GI,R GII, R snmax, R ssmax){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij;
R mnji,mtji,mbji;

R dxij,dyij,dzij; 

R ind1,ind2,indm;
R dx,dy,dz,nd2;
R Emui,rmui,ray,S,I;
R kkn,kkn0,kkt,kkt0,D;
R un,unc,unm;
R us,usc,usm;
R vs2,vn;

nint=0;
nsoft=0;

Epe=0.;
Epa=0.;

rmui=(rmu1+rmu2)/2.;
Emui=(Emu1+Emu2)/2.;
usm=2.*GII/ssmax;
unm=2.*GI/snmax;

int it,num_threads;
int jt,numth;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){

   for(it=0;it<nbco;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     

	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	 
	     
 
         if((LIST_P[numc1]!=LIST_P[numc2])&&TYPCO[it]){
			                        
		    ray=rmui*(LIST_R[numc1]+LIST_R[numc2])/2.;

			ind1=LIST_IND[numc1];
			ind2=LIST_IND[numc2];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;			
			
			dxij=(LIST_VX[numc1]-LIST_VX[numc2]);        	           
		    dyij=(LIST_VY[numc1]-LIST_VY[numc2]);
		    dzij=(LIST_VZ[numc1]-LIST_VZ[numc2]);

		    dx=LIST_X[numc1]-LIST_X[numc2]; 
		    dy=LIST_Y[numc1]-LIST_Y[numc2]; 			
		    dz=LIST_Z[numc1]-LIST_Z[numc2]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

			un=(dx/nd2)*(dx-DCONTX[it]);
			un+=(dy/nd2)*(dy-DCONTY[it]);	
			un+=(dz/nd2)*(dz-DCONTZ[it]);           

     	    kkn0=Emui*S/nd2;
		    unc=snmax*S/kkn0;
		    
		    nint++;
			
			if(un<=unc) {D=0.;}
			else if((un>unc)&&(un<unm)) {D=unm*(un-unc)/(un*(unm-unc));nsoft++;}
			else if(un>=unm) {D=1.;TYPCO[it]=0;
				
		    LIST_B[numc1]=2;
		    LIST_B[numc2]=2;  
		    
		    } 			
	
			kkn=(1.-D)*kkn0;	
				
		    VALCOH[it][0]=kkn;
		//	VALCOH[it][1]=kkt;
	
			indm=max(amort,1.);
		    VALAMO[it][0]=(amort*indm*dt/Cs)*kkn;
		//	VALAMO[it][1]=(amort*indm*dt/Cs)*kkt;
			 
	     } 
 
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                 
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;
                        
                   			
// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  


/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/	

              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  
              
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
	
	      }
          else{
			  
			fxx=0.;
			fyy=0.;
			fzz=0.;

			mxji=0.;
			myji=0.;
			mzji=0.;

			mxij=0.;
			myij=0.;
			mzij=0.;
			  
	      }   



              FCJI[it][0]=fxx;
              FCJI[it][1]=fyy;
              FCJI[it][2]=fzz; 

			  FX[numc1]+=fxx;
			  FY[numc1]+=fyy;
			  FZ[numc1]+=fzz;
			  FX[numc2]+=(-fxx);
			  FY[numc2]+=(-fyy);
			  FZ[numc2]+=(-fzz);
			  MTX[numc1]+=mxji;
			  MTX[numc2]+=mxij;
			  MTY[numc1]+=myji;
			  MTY[numc2]+=myij;	
			  MTZ[numc1]+=mzji;
			  MTZ[numc2]+=mzij;
			           
        }	
        
        
 
} 
 
 
else{

R FX2[NB_SPH][num_threads-1];
R FY2[NB_SPH][num_threads-1];
R FZ2[NB_SPH][num_threads-1];
R MTX2[NB_SPH][num_threads-1];
R MTY2[NB_SPH][num_threads-1];
R MTZ2[NB_SPH][num_threads-1];

# pragma omp parallel for collapse(2)
	for(it=0;it<NB_SPH;it++){ 

		for(jt=0;jt<num_threads-1;jt++){ 

		 FX2[it][jt]=0;
		 FY2[it][jt]=0;
		 FZ2[it][jt]=0;
		 MTX2[it][jt]=0;
		 MTY2[it][jt]=0;
		 MTZ2[it][jt]=0;

		}
	}
 


# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)+1) private(numth,jt,vs2,vn,us,usc,D,kkt,kkt0,dxij,dyij,dzij,dx,dy,dz,nd2,ray,ind1,ind2,indm,S,I,it,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij, dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,mnji,mtji,mbji) reduction(+:nint,nsoft,Epe,Epa)
	for(it=0;it<nbco;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     

	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	 
	     
 
         if((LIST_P[numc1]!=LIST_P[numc2])&&TYPCO[it]){
			                        
		    ray=rmui*(LIST_R[numc1]+LIST_R[numc2])/2.;

			ind1=LIST_IND[numc1];
			ind2=LIST_IND[numc2];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;			
			
			dxij=(LIST_VX[numc1]-LIST_VX[numc2]);        	           
		    dyij=(LIST_VY[numc1]-LIST_VY[numc2]);
		    dzij=(LIST_VZ[numc1]-LIST_VZ[numc2]);

		    dx=LIST_X[numc1]-LIST_X[numc2]; 
		    dy=LIST_Y[numc1]-LIST_Y[numc2]; 			
		    dz=LIST_Z[numc1]-LIST_Z[numc2]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

     	    kkn0=Emui*S/nd2;
		    unc=snmax*S/kkn0;
		    
		    nint++;
			
			if(un<=unc) {D=0.;}
			else if((un>unc)&&(un<unm)) {D=unm*(un-unc)/(un*(unm-unc));nsoft++;}
			else if(un>=unm) {D=1.;TYPCO[it]=0;
				
		    # pragma omp critical
		    { 	//dangereux
		    LIST_B[numc1]=2;
		    LIST_B[numc2]=2;  
		    } 
		    
		    } 			
	
			kkn=(1.-D)*kkn0;		 
			
		    VALCOH[it][0]=kkn;
		//	VALCOH[it][1]=kkt;

			indm=max(amort,1.);
		    VALAMO[it][0]=(amort*indm*dt/Cs)*kkn;
		//	VALAMO[it][1]=(amort*indm*dt/Cs)*kkt;
			 
	     } 
 
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                 
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;
                        
                   			
// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  


/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/	

              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  
              
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
	
	      }
          else{
			  
			fxx=0.;
			fyy=0.;
			fzz=0.;

			mxji=0.;
			myji=0.;
			mzji=0.;

			mxij=0.;
			myij=0.;
			mzij=0.;
			  
	      }   



              FCJI[it][0]=fxx;
              FCJI[it][1]=fyy;
              FCJI[it][2]=fzz; 

			  numth=omp_get_thread_num()-1;

				if(numth==-1){
					  FX[numc1]+=fxx;
					  FY[numc1]+=fyy;
					  FZ[numc1]+=fzz;
					  FX[numc2]+=(-fxx);
					  FY[numc2]+=(-fyy);
					  FZ[numc2]+=(-fzz);
					  MTX[numc1]+=mxji;
					  MTX[numc2]+=mxij;
					  MTY[numc1]+=myji;
					  MTY[numc2]+=myij;	
					  MTZ[numc1]+=mzji;
					  MTZ[numc2]+=mzij;
				}else{


					  FX2[numc1][numth]+=fxx;
					  FY2[numc1][numth]+=fyy;
					  FZ2[numc1][numth]+=fzz;
					  FX2[numc2][numth]+=(-fxx);
					  FY2[numc2][numth]+=(-fyy);
					  FZ2[numc2][numth]+=(-fzz);
					  MTX2[numc1][numth]+=mxji;
					  MTX2[numc2][numth]+=mxij;
					  MTY2[numc1][numth]+=myji;
					  MTY2[numc2][numth]+=myij;	
					  MTZ2[numc1][numth]+=mzji;
					  MTZ2[numc2][numth]+=mzij;


				}
			           
        }	




# pragma omp parallel for collapse(2) 
	for(it=0;it<NB_SPH;it++){
 
		for(jt=0;jt<num_threads-1;jt++){ 

			FX[it]+=FX2[it][jt];
			FY[it]+=FY2[it][jt];
			FZ[it]+=FZ2[it][jt];
			MTX[it]+=MTX2[it][jt];
			MTY[it]+=MTY2[it][jt];
			MTZ[it]+=MTZ2[it][jt];
		}


        }





}

      

}



void forcecoh2_qs_int_mII(int & nint,int & nsoft,R Cs, R & Epe, R & Epa,R dt, bool * TYPCO,int & nbco, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ,  R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI,R ** FOJI,R ** MTJI,R ** MTIJ, bool * LIST_P, int * LIST_B, R * vs, R * LIST_IND,R Pi,R Emu1, R rmu1, R Emu2, R rmu2, R amort, R GI,R GII, R snmax, R ssmax){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij;
R mnji,mtji,mbji;

R dxij,dyij,dzij; 

R ind1,ind2,indm;
R dx,dy,dz,nd2;
R Emui,rmui,ray,S,I;
R kkn,kkn0,kkt,kkt0,D;
R un,unc,unm;
R us,usc,usm;
R vs2,vn;

nint=0;
nsoft=0;

Epe=0.;
Epa=0.;

rmui=(rmu1+rmu2)/2.;
Emui=(Emu1+Emu2)/2.;
usm=2.*GII/ssmax;

int it,num_threads;
int jt,numth;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){

   for(it=0;it<nbco;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     

	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	 
	     
 
         if((LIST_P[numc1]!=LIST_P[numc2])&&TYPCO[it]){
			                        
		    ray=rmui*(LIST_R[numc1]+LIST_R[numc2])/2.;

			ind1=LIST_IND[numc1];
			ind2=LIST_IND[numc2];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;			
			
			dxij=(LIST_VX[numc1]-LIST_VX[numc2]);        	           
		    dyij=(LIST_VY[numc1]-LIST_VY[numc2]);
		    dzij=(LIST_VZ[numc1]-LIST_VZ[numc2]);

		    dx=LIST_X[numc1]-LIST_X[numc2]; 
		    dy=LIST_Y[numc1]-LIST_Y[numc2]; 			
		    dz=LIST_Z[numc1]-LIST_Z[numc2]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

		    vn=(dx/nd2)*dxij+(dy/nd2)*dyij+(dz/nd2)*dzij;
		    vs2=(dxij*dxij+dyij*dyij+dzij*dzij-vn*vn);
		    if(vs2>0){vs[it]+=sqrt(vs2);}
		    	 
		    us=vs[it]*dt;
			
     		kkt0=12.*Emui*I/(pow(nd2,3)); 
	
		
		    usc=ssmax*S/kkt0;
			nint++;
			if(us<=usc) {D=0.;}
			else if((us>usc)&&(us<usm)) {D=usm*(us-usc)/(us*(usm-usc));nsoft++;}
			else if(us>=usm) {D=1.;TYPCO[it]=0;
				
		    LIST_B[numc1]=2;
		    LIST_B[numc2]=2;  
		    
		    } 			
	
			kkt=(1.-D)*kkt0;

	
				
		//    VALCOH[it][0]=kkn;
			VALCOH[it][1]=kkt;
	
			indm=max(amort,1.);
		  //  VALAMO[it][0]=(amort*indm*dt/Cs)*kkn;
			VALAMO[it][1]=(amort*indm*dt/Cs)*kkt;
			 
	     } 
 
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                 
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;
                        
                   			
// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  


/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/	

              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  
              
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
	
	      }
          else{
			  
			fxx=0.;
			fyy=0.;
			fzz=0.;

			mxji=0.;
			myji=0.;
			mzji=0.;

			mxij=0.;
			myij=0.;
			mzij=0.;
			  
	      }   



              FCJI[it][0]=fxx;
              FCJI[it][1]=fyy;
              FCJI[it][2]=fzz; 

			  FX[numc1]+=fxx;
			  FY[numc1]+=fyy;
			  FZ[numc1]+=fzz;
			  FX[numc2]+=(-fxx);
			  FY[numc2]+=(-fyy);
			  FZ[numc2]+=(-fzz);
			  MTX[numc1]+=mxji;
			  MTX[numc2]+=mxij;
			  MTY[numc1]+=myji;
			  MTY[numc2]+=myij;	
			  MTZ[numc1]+=mzji;
			  MTZ[numc2]+=mzij;
			           
        }	
        
        
 
} 
 
 
else{

R FX2[NB_SPH][num_threads-1];
R FY2[NB_SPH][num_threads-1];
R FZ2[NB_SPH][num_threads-1];
R MTX2[NB_SPH][num_threads-1];
R MTY2[NB_SPH][num_threads-1];
R MTZ2[NB_SPH][num_threads-1];

# pragma omp parallel for collapse(2)
	for(it=0;it<NB_SPH;it++){ 

		for(jt=0;jt<num_threads-1;jt++){ 

		 FX2[it][jt]=0;
		 FY2[it][jt]=0;
		 FZ2[it][jt]=0;
		 MTX2[it][jt]=0;
		 MTY2[it][jt]=0;
		 MTZ2[it][jt]=0;

		}
	}
 


# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)+1) private(numth,jt,vs2,vn,us,usc,D,kkt,kkt0,dxij,dyij,dzij,dx,dy,dz,nd2,ray,ind1,ind2,indm,S,I,it,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij, dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,mnji,mtji,mbji) reduction(+:nint,nsoft,Epe,Epa)
	for(it=0;it<nbco;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     

	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	 
	     
 
         if((LIST_P[numc1]!=LIST_P[numc2])&&TYPCO[it]){
			                        
		    ray=rmui*(LIST_R[numc1]+LIST_R[numc2])/2.;

			ind1=LIST_IND[numc1];
			ind2=LIST_IND[numc2];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;			
			
			dxij=(LIST_VX[numc1]-LIST_VX[numc2]);        	           
		    dyij=(LIST_VY[numc1]-LIST_VY[numc2]);
		    dzij=(LIST_VZ[numc1]-LIST_VZ[numc2]);

		    dx=LIST_X[numc1]-LIST_X[numc2]; 
		    dy=LIST_Y[numc1]-LIST_Y[numc2]; 			
		    dz=LIST_Z[numc1]-LIST_Z[numc2]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

		    vn=(dx/nd2)*dxij+(dy/nd2)*dyij+(dz/nd2)*dzij;
		    vs2=(dxij*dxij+dyij*dyij+dzij*dzij-vn*vn);
		    if(vs2>0){vs[it]+=sqrt(vs2);}
		    	 
		    us=vs[it]*dt;
			
     		kkt0=12.*Emui*I/(pow(nd2,3));	 
	
		
		    usc=ssmax*S/kkt0;
			nint++;
			if(us<=usc) {D=0.;}
			else if((us>usc)&&(us<usm)) {D=usm*(us-usc)/(us*(usm-usc));nsoft++;}
			else if(us>=usm) {D=1.;TYPCO[it]=0;
				
		    # pragma omp critical
		    { 	//dangereux
		    LIST_B[numc1]=2;
		    LIST_B[numc2]=2;  
		    }			
		    
		    } 			
	
			kkt=(1.-D)*kkt0;

	
				
		//    VALCOH[it][0]=kkn;
			VALCOH[it][1]=kkt;
	
			indm=max(amort,1.);
		  //  VALAMO[it][0]=(amort*indm*dt/Cs)*kkn;
			VALAMO[it][1]=(amort*indm*dt/Cs)*kkt;
			 
	     } 
 
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                 
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;
                        
                   			
// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  


/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/	

              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  
              
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
	
	      }
          else{
			  
			fxx=0.;
			fyy=0.;
			fzz=0.;

			mxji=0.;
			myji=0.;
			mzji=0.;

			mxij=0.;
			myij=0.;
			mzij=0.;
			  
	      }   



              FCJI[it][0]=fxx;
              FCJI[it][1]=fyy;
              FCJI[it][2]=fzz; 

			  numth=omp_get_thread_num()-1;

				if(numth==-1){
					  FX[numc1]+=fxx;
					  FY[numc1]+=fyy;
					  FZ[numc1]+=fzz;
					  FX[numc2]+=(-fxx);
					  FY[numc2]+=(-fyy);
					  FZ[numc2]+=(-fzz);
					  MTX[numc1]+=mxji;
					  MTX[numc2]+=mxij;
					  MTY[numc1]+=myji;
					  MTY[numc2]+=myij;	
					  MTZ[numc1]+=mzji;
					  MTZ[numc2]+=mzij;
				}else{


					  FX2[numc1][numth]+=fxx;
					  FY2[numc1][numth]+=fyy;
					  FZ2[numc1][numth]+=fzz;
					  FX2[numc2][numth]+=(-fxx);
					  FY2[numc2][numth]+=(-fyy);
					  FZ2[numc2][numth]+=(-fzz);
					  MTX2[numc1][numth]+=mxji;
					  MTX2[numc2][numth]+=mxij;
					  MTY2[numc1][numth]+=myji;
					  MTY2[numc2][numth]+=myij;	
					  MTZ2[numc1][numth]+=mzji;
					  MTZ2[numc2][numth]+=mzij;


				}
			           
        }	




# pragma omp parallel for collapse(2) 
	for(it=0;it<NB_SPH;it++){
 
		for(jt=0;jt<num_threads-1;jt++){ 

			FX[it]+=FX2[it][jt];
			FY[it]+=FY2[it][jt];
			FZ[it]+=FZ2[it][jt];
			MTX[it]+=MTX2[it][jt];
			MTY[it]+=MTY2[it][jt];
			MTZ[it]+=MTZ2[it][jt];
		}


        }





}
      
}







void forcecoh2_qs_int_mm(int & nint,int & nsoft,R Cs, R & Epe, R & Epa,R dt, bool * TYPCO,int & nbco, int NB_SPH,  int ** CONT,  R * LIST_R, R * LIST_X, R * LIST_Y, R * LIST_Z, R * LIST_TX, R * LIST_TY, R * LIST_TZ, R * LIST_VX, R * LIST_VY, R * LIST_VZ,R * LIST_WX, R * LIST_WY, R * LIST_WZ,  R * DCONTX, R * DCONTY, R * DCONTZ, R ** NCONT,  R * FX, R * FY, R * FZ, R * MTX, R * MTY, R * MTZ, R ** VALCOH, R ** VALAMO, R ** FCJI,R ** FOJI,R ** MTJI,R ** MTIJ, bool * LIST_P, int * LIST_B, R * vs, R * LIST_IND,R Pi,R Emu1, R rmu1, R Emu2, R rmu2, R amort, R GI,R GII, R snmax, R ssmax){

int numc1,numc2;
R n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj;
R thetani,thetanj,thetati,thetatj,thetabi,thetabj;
R vn_ij,vt_ij,vb_ij;
R dxoij,dyoij,dzoij; 
R fxx,fyy,fzz;
R mtij,mbij;
R mnji,mtji,mbji;

R dxij,dyij,dzij; 

R ind1,ind2,indm;
R dx,dy,dz,nd2;
R Emui,rmui,ray,S,I;
R kkn,kkn0,kkt,kkt0,D;
R un,unc,unm;
R us,usc,usm;
R vs2,vn;
R ue,uem,uec,eta,alpha;

nint=0;
nsoft=0;

Epe=0.;
Epa=0.;

rmui=(rmu1+rmu2)/2.;
Emui=(Emu1+Emu2)/2.;
unm=2.*GI/snmax;
usm=2.*GII/ssmax;
alpha=1.;

int it,num_threads;
int jt,numth;

# pragma omp parallel
{
num_threads=omp_get_num_threads();
}

if(num_threads==1){

   for(it=0;it<nbco;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     

	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	 
	     
          if((LIST_P[numc1]!=LIST_P[numc2])&&TYPCO[it]){
			                        
		    ray=rmui*(LIST_R[numc1]+LIST_R[numc2])/2.;

			ind1=LIST_IND[numc1];
			ind2=LIST_IND[numc2];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;			
			
			dxij=(LIST_VX[numc1]-LIST_VX[numc2]);        	           
		    dyij=(LIST_VY[numc1]-LIST_VY[numc2]);
		    dzij=(LIST_VZ[numc1]-LIST_VZ[numc2]);

		    dx=LIST_X[numc1]-LIST_X[numc2]; 
		    dy=LIST_Y[numc1]-LIST_Y[numc2]; 			
		    dz=LIST_Z[numc1]-LIST_Z[numc2]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

		    vn=(dx/nd2)*dxij+(dy/nd2)*dyij+(dz/nd2)*dzij;
		    vs2=(dxij*dxij+dyij*dyij+dzij*dzij-vn*vn);
		    if(vs2>0){vs[it]+=sqrt(vs2);}
		    	 
		    us=vs[it]*dt;
			
			un=(dx/nd2)*(dx-DCONTX[it]);
			un+=(dy/nd2)*(dy-DCONTY[it]);	
			un+=(dz/nd2)*(dz-DCONTZ[it]);			
			
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
			else if(ue>=uem) {D=1.;TYPCO[it]=0;
								
		    # pragma omp critical
		    { 	//dangereux
		    LIST_B[numc1]=2;
		    LIST_B[numc2]=2;  
		    }			
		    
		    } 			
	
			kkn=(1.-D)*kkn0;	
			kkt=(1.-D)*kkt0;
				
		    VALCOH[it][0]=kkn;
			VALCOH[it][1]=kkt;
	
indm=max(amort,1.);
		    VALAMO[it][0]=(amort*indm*dt/Cs)*kkn;
			VALAMO[it][1]=(amort*indm*dt/Cs)*kkt;
			 
	     } 
 
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                 
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;
                        
                   			
// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  


/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/	

              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  
              
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
	
	      }
          else{
			  
			fxx=0.;
			fyy=0.;
			fzz=0.;

			mxji=0.;
			myji=0.;
			mzji=0.;

			mxij=0.;
			myij=0.;
			mzij=0.;
			  
	      }   



              FCJI[it][0]=fxx;
              FCJI[it][1]=fyy;
              FCJI[it][2]=fzz; 

			  FX[numc1]+=fxx;
			  FY[numc1]+=fyy;
			  FZ[numc1]+=fzz;
			  FX[numc2]+=(-fxx);
			  FY[numc2]+=(-fyy);
			  FZ[numc2]+=(-fzz);
			  MTX[numc1]+=mxji;
			  MTX[numc2]+=mxij;
			  MTY[numc1]+=myji;
			  MTY[numc2]+=myij;	
			  MTZ[numc1]+=mzji;
			  MTZ[numc2]+=mzij;
			           
        }	
        
        
 
} 
 
 
else{

R FX2[NB_SPH][num_threads-1];
R FY2[NB_SPH][num_threads-1];
R FZ2[NB_SPH][num_threads-1];
R MTX2[NB_SPH][num_threads-1];
R MTY2[NB_SPH][num_threads-1];
R MTZ2[NB_SPH][num_threads-1];

# pragma omp parallel for collapse(2)
	for(it=0;it<NB_SPH;it++){ 

		for(jt=0;jt<num_threads-1;jt++){ 

		 FX2[it][jt]=0;
		 FY2[it][jt]=0;
		 FZ2[it][jt]=0;
		 MTX2[it][jt]=0;
		 MTY2[it][jt]=0;
		 MTZ2[it][jt]=0;

		}
	}
 


# pragma omp parallel for schedule(dynamic,int(nbco/num_threads)+1) private(numth,jt,vs2,eta,vn,un,unc,us,usc,ue,uec,uem,D,kkn,kkn0,kkt,kkt0,dxij,dyij,dzij,dx,dy,dz,nd2,ray,ind1,ind2,indm,S,I,it,numc1,numc2,n1,n2,n3,s1,s2,s3,t1,t2,t3,fnji,ftji,fbji,mxij,mxji,myij,myji,mzij,mzji,un_ij,ut_ij,ub_ij,thetaxi,thetaxj,thetayi,thetayj,thetazi,thetazj,thetani,thetanj,thetati,thetatj,thetabi,thetabj,vn_ij,vt_ij,vb_ij, dxoij,dyoij,dzoij,fxx,fyy,fzz,mtij,mbij,mnji,mtji,mbji) reduction(+:nint,nsoft,Epe,Epa)
	for(it=0;it<nbco;it++){ 
		
		 // Numeros des candidats 
	     numc1=CONT[it][0];
	     numc2=CONT[it][1];     

	     // normal ext 1->2
	     n1=NCONT[it][0];
	     n2=NCONT[it][1];
	     n3=NCONT[it][2];
	     s1=NCONT[it][3];
	     s2=NCONT[it][4];
	     s3=NCONT[it][5];
	     t1=NCONT[it][6];
	     t2=NCONT[it][7];
	     t3=NCONT[it][8];  	 
	     

         if((LIST_P[numc1]!=LIST_P[numc2])&&TYPCO[it]){
			
                        
		    ray=rmui*(LIST_R[numc1]+LIST_R[numc2])/2.;

			ind1=LIST_IND[numc1];
			ind2=LIST_IND[numc2];				                             
			indm=max(ind1,ind2); 
			S=indm*(Pi*ray*ray);
			I=indm*Pi*pow((2.*ray),4)/64.;			
			
			dxij=(LIST_VX[numc1]-LIST_VX[numc2]);        	           
		    dyij=(LIST_VY[numc1]-LIST_VY[numc2]);
		    dzij=(LIST_VZ[numc1]-LIST_VZ[numc2]);

		    dx=LIST_X[numc1]-LIST_X[numc2]; 
		    dy=LIST_Y[numc1]-LIST_Y[numc2]; 			
		    dz=LIST_Z[numc1]-LIST_Z[numc2]; 					
			nd2=sqrt(dx*dx+dy*dy+dz*dz);

		    vn=(dx/nd2)*dxij+(dy/nd2)*dyij+(dz/nd2)*dzij;
		    vs2=(dxij*dxij+dyij*dyij+dzij*dzij-vn*vn);
		    if(vs2>0){vs[it]+=sqrt(vs2);}
		    	 
		    us=vs[it]*dt;
			
			un=(dx/nd2)*(dx-DCONTX[it]);
			un+=(dy/nd2)*(dy-DCONTY[it]);	
			un+=(dz/nd2)*(dz-DCONTZ[it]);			
			
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
			else if(ue>=uem) {D=1.;TYPCO[it]=0;
								
		    # pragma omp critical
		    { 	//dangereux
		    LIST_B[numc1]=2;
		    LIST_B[numc2]=2;  
		    }			
		    
		    } 			
	
			kkn=(1.-D)*kkn0;	
			kkt=(1.-D)*kkt0;
				
		    VALCOH[it][0]=kkn;
			VALCOH[it][1]=kkt;
	
			indm=max(amort,1.);
		    VALAMO[it][0]=(amort*indm*dt/Cs)*kkn;
			VALAMO[it][1]=(amort*indm*dt/Cs)*kkt;
			 
	     } 
 
 
         if(TYPCO[it]==1){      
            
         // Déplacements relatifs dans le repère global
             dxoij=(LIST_X[numc1]-LIST_X[numc2])-DCONTX[it];        	           
             dyoij=(LIST_Y[numc1]-LIST_Y[numc2])-DCONTY[it];
             dzoij=(LIST_Z[numc1]-LIST_Z[numc2])-DCONTZ[it];

         // Déplacements relatifs dans le repère local
             un_ij=n1*dxoij+n2*dyoij+n3*dzoij;
             ut_ij=s1*dxoij+s2*dyoij+s3*dzoij;
             ub_ij=t1*dxoij+t2*dyoij+t3*dzoij; 
             
          // Rotations      
             thetaxi=LIST_TX[numc1];
             thetaxj=LIST_TX[numc2];
             thetayi=LIST_TY[numc1];
             thetayj=LIST_TY[numc2];
             thetazi=LIST_TZ[numc1];
             thetazj=LIST_TZ[numc2];
             
         
             thetani=n1*thetaxi+n2*thetayi+n3*thetazi;   
             thetati=s1*thetaxi+s2*thetayi+s3*thetazi;                       
             thetabi=t1*thetaxi+t2*thetayi+t3*thetazi;            

             thetanj=n1*thetaxj+n2*thetayj+n3*thetazj;   
             thetatj=s1*thetaxj+s2*thetayj+s3*thetazj;                       
             thetabj=t1*thetaxj+t2*thetayj+t3*thetazj;  
            
        
          // Vitesses relatives dans le repère local
			 vn_ij=n1*(LIST_VX[numc1]-LIST_VX[numc2])+n2*(LIST_VY[numc1]-LIST_VY[numc2])+n3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vt_ij=s1*(LIST_VX[numc1]-LIST_VX[numc2])+s2*(LIST_VY[numc1]-LIST_VY[numc2])+s3*(LIST_VZ[numc1]-LIST_VZ[numc2]);
			 vb_ij=t1*(LIST_VX[numc1]-LIST_VX[numc2])+t2*(LIST_VY[numc1]-LIST_VY[numc2])+t3*(LIST_VZ[numc1]-LIST_VZ[numc2]);	
                 
             
             
// Modèle Poutre 3D
//	ES/L	 VALCOH[nbco][0]=kn;
//  12EI/L3	 VALCOH[nbco][1]=kt;
//	6EI/L2	 VALCOH[nbco][2]=nd2*kt/2.;	
//	4EI/L2	 VALCOH[nbco][3]=nd2*nd2*kt/3.;	
//	2EI/L2	 VALCOH[nbco][4]=nd2*nd2*kt/6.;
//	2GI/L	 VALCOH[nbco][5]		 

			  fnji=-VALCOH[it][0]*un_ij;
			  ftji=-VALCOH[it][1]*ut_ij-VALCOH[it][2]*(thetabi+thetabj);
			  fbji=-VALCOH[it][1]*ub_ij+VALCOH[it][2]*(thetati+thetatj);  
			  
			  mnji=-VALCOH[it][5]*(thetani-thetanj);
			  mtji=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetati-VALCOH[it][4]*thetatj;
			  mbji=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabi-VALCOH[it][4]*thetabj;
			 			 
			  mtij=VALCOH[it][2]*ub_ij-VALCOH[it][3]*thetatj-VALCOH[it][4]*thetati;  
                          mbij=-VALCOH[it][2]*ut_ij-VALCOH[it][3]*thetabj-VALCOH[it][4]*thetabi;

			Epe+=-fnji*un_ij/2.;
			Epe+=-ftji*ut_ij/2.;                      
			Epe+=-fbji*ub_ij/2.;      

			Epe+=-mnji*(thetani)/2.;    
			Epe+=-mtji*(thetati)/2.;                   
			Epe+=-mbji*(thetabi)/2.;                     

			Epe+=-(-mnji)*(thetanj)/2.;                      
			Epe+=-mtij*(thetatj)/2.;                   
			Epe+=-mbij*(thetabj)/2.;
                        
                   			
// Amortissement 

			  fnji=fnji-(VALAMO[it][0]*vn_ij);
			  ftji=ftji-(VALAMO[it][1]*vt_ij);
			  fbji=fbji-(VALAMO[it][1]*vb_ij);		  


/*
				Epa+=(0.5*VALAMO[it][0]*vn_ij*un_ij);
				Epa+=(0.5*VALAMO[it][1]*ut_ij*vt_ij); 
				Epa+=(0.5*VALAMO[it][1]*ub_ij*vb_ij); 
	
		*/	

              FOJI[it][0]=fnji;
              FOJI[it][1]=ftji;
              FOJI[it][2]=fbji;
              
              MTJI[it][0]=mnji;
              MTJI[it][1]=mtji;
              MTJI[it][2]=mbji;                                      
	
              MTIJ[it][0]=-mnji;
              MTIJ[it][1]=mtij;
              MTIJ[it][2]=mbij;  
              
			  fxx=n1*fnji+s1*ftji+t1*fbji;
			  fyy=n2*fnji+s2*ftji+t2*fbji;
			  fzz=n3*fnji+s3*ftji+t3*fbji;
	
	          mxji=n1*mnji+s1*mtji+t1*mbji; 
	          myji=n2*mnji+s2*mtji+t2*mbji; 
	          mzji=n3*mnji+s3*mtji+t3*mbji; 

	          mxij=-n1*mnji+s1*mtij+t1*mbij; 
	          myij=-n2*mnji+s2*mtij+t2*mbij; 
	          mzij=-n3*mnji+s3*mtij+t3*mbij; 
	
	      }
          else{
			  
			fxx=0.;
			fyy=0.;
			fzz=0.;

			mxji=0.;
			myji=0.;
			mzji=0.;

			mxij=0.;
			myij=0.;
			mzij=0.;
			  
	      }   



              FCJI[it][0]=fxx;
              FCJI[it][1]=fyy;
              FCJI[it][2]=fzz; 

			  numth=omp_get_thread_num()-1;

				if(numth==-1){
					  FX[numc1]+=fxx;
					  FY[numc1]+=fyy;
					  FZ[numc1]+=fzz;
					  FX[numc2]+=(-fxx);
					  FY[numc2]+=(-fyy);
					  FZ[numc2]+=(-fzz);
					  MTX[numc1]+=mxji;
					  MTX[numc2]+=mxij;
					  MTY[numc1]+=myji;
					  MTY[numc2]+=myij;	
					  MTZ[numc1]+=mzji;
					  MTZ[numc2]+=mzij;
				}else{


					  FX2[numc1][numth]+=fxx;
					  FY2[numc1][numth]+=fyy;
					  FZ2[numc1][numth]+=fzz;
					  FX2[numc2][numth]+=(-fxx);
					  FY2[numc2][numth]+=(-fyy);
					  FZ2[numc2][numth]+=(-fzz);
					  MTX2[numc1][numth]+=mxji;
					  MTX2[numc2][numth]+=mxij;
					  MTY2[numc1][numth]+=myji;
					  MTY2[numc2][numth]+=myij;	
					  MTZ2[numc1][numth]+=mzji;
					  MTZ2[numc2][numth]+=mzij;


				}
			           
        }	




# pragma omp parallel for collapse(2) 
	for(it=0;it<NB_SPH;it++){
 
		for(jt=0;jt<num_threads-1;jt++){ 

			FX[it]+=FX2[it][jt];
			FY[it]+=FY2[it][jt];
			FZ[it]+=FZ2[it][jt];
			MTX[it]+=MTX2[it][jt];
			MTY[it]+=MTY2[it][jt];
			MTZ[it]+=MTZ2[it][jt];
		}


        }





}

      

}



