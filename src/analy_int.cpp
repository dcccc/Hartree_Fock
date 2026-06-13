#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include "math.h"
using namespace std;
using namespace boost;


extern "C" double boys(double n, double x);
extern "C" double s_ab(double x1,double i1,double a1,double x2,double i2,double a2);
extern "C" double sab(double* ga,double* gb);
extern "C" double t_ab(double x1,double i1,double a1,double x2,double i2,double a2);
extern "C" double tab(double* ga,double* gb);
extern "C" double hermite(double x1,int i1,double a1,double x2,int i2,double a2,int t);
extern "C" double ruvt(int u,int v,int t,double* ga,double* gb,double n,double* nu);
extern "C" double vab(double* ga,double* gb,double* nu);
extern "C" double ruvt2(int u,int v,int t,double* ga,double* gb,double* gc,double* gd,double n);
extern "C" double gabcd(double* ga,double* gb,double* gc,double* gd);
extern "C" void smat_pri(double* ba,double* s_mat,double* s_mat_dia,int len);
extern "C" void smat_con(double* ba,double* s_mat,double* s_mat_dia,int* contract_list,int co_len,int ba_len);
extern "C" void tmat_pri(double* ba,double* s_mat,double* s_mat_dia,int len);
extern "C" void tmat_con(double* ba,double* t_mat,double* s_mat_dia,int* contract_list,int len);
extern "C" void vmat_con(double* ba,double* v_mat,double* s_mat_dia,int* contract_list,double* nu,int con_len,int nu_len);
extern "C" void vmat_pri(double* ba,double* v_mat,double* s_mat_dia,double* nu,int ba_len,int nu_len);
extern "C" void gab_pri(double* ba,double* g_mat,double* s_mat_dia,int len);
extern "C" void gab_con(double* ba,double* g_mat,double* s_mat_dia,int* contract_list,int con_len);



double boys(double n, double x){
	double b=0.0;
	if ( x==0.0 ){
		b=1.0/(2.0*n+1.0);
	}
	else if (n==0.0 and x !=0.0 ){
		b=pow(3.1415926,0.5)*erf(pow(x,0.5))/pow(x,0.5)*0.5;
	}
	else {
		b=0.5/pow(x,n+0.5)*boost::math::tgamma(n+0.5)*boost::math::gamma_p(n+0.5,x);
	}
	return b;
}


double s_ab(double x1,double i1,double a1,double x2,double i2,double a2){
	double p=a1+a2;
	double px=(a1*x1+a2*x2)/p;
	double kab=exp(-a1*a2/(a1+a2)*pow(x1-x2,2));
	double s=0.0;
	
	if (i1<0.0 or i2<0.0){
		s=0.0;
	}
	else if (i1==0.0 and i2==0.0){
		s=pow(3.1415926/p,0.5)*kab;
	}
	else if (i1>0.0 and i2==0.0){
		s=(px-x1)*s_ab(x1,i1-1,a1,x2,i2,a2)+0.5/p*((i1-1.0)*s_ab(x1,i1-2.0,a1,x2,i2,a2));
	}
	else{
		s=(px-x2)*s_ab(x1,i1,a1,x2,i2-1.0,a2)+0.5/p*(i1*s_ab(x1,i1-1.0,a1,x2,i2-1.0,a2)+
		  (i2-1.0)*s_ab(x1,i1,a1,x2,i2-2.0,a2));
	}
	return s;
}


double sab(double* ga,double* gb){
	double s=s_ab(ga[0],ga[3],ga[6],gb[0],gb[3],gb[6])*
	         s_ab(ga[1],ga[4],ga[6],gb[1],gb[4],gb[6])*
	         s_ab(ga[2],ga[5],ga[6],gb[2],gb[5],gb[6]);
	return s;
}

double t_ab(double x1,double i1,double a1,double x2,double i2,double a2){
	double p=a1+a2;
	double px=(a1*x1+a2*x2)/p;
	double xa=px-x1;
	double xb=px-x2;
	double t=0.0;
	
	if (i1<0.0 or i2<0.0){
		t=0.0;
	}
	else if (i1==0.0 and i2==0.0){
		t=(a1-2.0*a1*a1*(xa*xa+0.5/p))*s_ab(x1,0,a1,x2,0,a2);
	}
	else if (i1>0.0 and i2==0.0){
		t=xa*t_ab(x1,i1-1.0,a1,x2,0.0,a2)+0.5/p*((i1-1.0)*t_ab(x1,i1-2.0,a1,x2,0.0,a2))+
		a2/p*(2.0*a1*s_ab(x1,i1,a1,x2,0.0,a2)-(i1-1.0)*s_ab(x1,i1-2.0,a1,x2,0.0,a2));
	}
	else{
		t=xb*t_ab(x1,i1,a1,x2,i2-1.0,a2)+
		0.5/p*(i1*t_ab(x1,i1-1.0,a1,x2,i2-1.0,a2)+(i2-1.0)*t_ab(x1,i1,a1,x2,i2-2.0,a2))+
		a1/p*(2.0*a2*s_ab(x1,i1,a1,x2,i2,a2)-(i2-1.0)*s_ab(x1,i1,a1,x2,i2-2.0,a2));
	}
	return t;
}



double tab(double* ga,double* gb){
	double t=t_ab(ga[0],ga[3],ga[6],gb[0],gb[3],gb[6])*
	         s_ab(ga[1],ga[4],ga[6],gb[1],gb[4],gb[6])*
	         s_ab(ga[2],ga[5],ga[6],gb[2],gb[5],gb[6])+
	         s_ab(ga[0],ga[3],ga[6],gb[0],gb[3],gb[6])*
	         t_ab(ga[1],ga[4],ga[6],gb[1],gb[4],gb[6])*
	         s_ab(ga[2],ga[5],ga[6],gb[2],gb[5],gb[6])+
	         s_ab(ga[0],ga[3],ga[6],gb[0],gb[3],gb[6])*
	         s_ab(ga[1],ga[4],ga[6],gb[1],gb[4],gb[6])*
	         t_ab(ga[2],ga[5],ga[6],gb[2],gb[5],gb[6]);
	return t;
}


double hermite(double x1,int i1,double a1,double x2, int i2,double a2, int t){
	double p=a1+a2;
	double e=0.0;
	
	if (t<0 or t > i1+i2 or i1<0 or i2<0){
		e=0.0;
	}
	else if (t==0){
		e=pow(p/3.1415926,0.5)*s_ab(x1,i1,a1,x2,i2,a2);
	}
	else{
		e=(i1*hermite(x1,i1-1,a1,x2,i2,a2,t-1)+i2*hermite(x1,i1,a1,x2,i2-1,a2,t-1))*0.5/p/t;
	}
	return e;
}



double ruvt(int u,int v,int t,double* ga,double* gb,double n,double* nu){
	double p=ga[6]+gb[6];
	double px=(ga[6]*ga[0]+gb[6]*gb[0])/p;
	double py=(ga[6]*ga[1]+gb[6]*gb[1])/p;
	double pz=(ga[6]*ga[2]+gb[6]*gb[2])/p;
	double xa=px-nu[0];
	double ya=py-nu[1];
	double za=pz-nu[2];
	double r=0.0;
	
	if (u<0 or v<0 or t<0 or n<0.0){
		r=0.0;
	}
	else if (u==0 and v==0 and t==0){
		r=pow(-2.0*p,n)*boys(n,p*(xa*xa+ya*ya+za*za));

	}
	else if (u>0 and v==0 and t==0){
		r=(u-1)*ruvt(u-2,0,0,ga,gb,n+1.0,nu)+xa*ruvt(u-1,0,0,ga,gb,n+1.0,nu);
	}
	else if (v>0 and t==0){
		r=(v-1)*ruvt(u,v-2,0.0,ga,gb,n+1.0,nu)+ya*ruvt(u,v-1,0,ga,gb,n+1.0,nu);
	}
	else{
		r=(t-1)*ruvt(u,v,t-2,ga,gb,n+1.0,nu)+za*ruvt(u,v,t-1,ga,gb,n+1.0,nu);
	}
	return r;
}


double vab(double* ga,double* gb,double* nu){
	double p=ga[6]+gb[6];
	double vv=0.0;
	
	for (int u=0;u<(int) ga[3]+gb[3]+1.0;u+=1){
	for (int v=0;v<(int) ga[4]+gb[4]+1.0;v+=1){
	for (int t=0;t<(int) ga[5]+gb[5]+1.0;t+=1){
			vv+=hermite(ga[0], ga[3],ga[6],gb[0],gb[3],gb[6],u)*
			    hermite(ga[1], ga[4],ga[6],gb[1],gb[4],gb[6],v)*
				hermite(ga[2], ga[5],ga[6],gb[2],gb[5],gb[6],t)*
				ruvt(u,v,t,ga,gb,0.0,nu);
	}}}
	return vv*2.0*3.1415926*nu[3]/p;
}
	
	

double ruvt2(int u,int v,int t,double* ga,double* gb,double* gc,double* gd,double n){
	double p1=ga[6]+gb[6];
	double p2=gc[6]+gd[6];
	double pp=p1*p2/(p1+p2);
	double px=(ga[6]*ga[0]+gb[6]*gb[0])/p1-(gc[6]*gc[0]+gd[6]*gd[0])/p2;
	double py=(ga[6]*ga[1]+gb[6]*gb[1])/p1-(gc[6]*gc[1]+gd[6]*gd[1])/p2;
	double pz=(ga[6]*ga[2]+gb[6]*gb[2])/p1-(gc[6]*gc[2]+gd[6]*gd[2])/p2;
	
	double r=0.0;
	
	if (u<0 or v<0 or t<0 or n<0){
		r=0.0;
	}
	else if (u==0 and v==0 and t==0){
		r=pow(-2.0*pp,n)*boys(n,pp*(px*px+py*py+pz*pz));
	}
	else if (u>0){
		r=(u-1)*ruvt2(u-2,v,t,ga,gb,gc,gd,n+1.0)+px*ruvt2(u-1,v,t,ga,gb,gc,gd,n+1.0);
	}
	else if (v>0){
		r=(v-1)*ruvt2(u,v-2,t,ga,gb,gc,gd,n+1.0)+py*ruvt2(u,v-1,t,ga,gb,gc,gd,n+1.0);
	}
	else if (t>0){
		r=(t-1)*ruvt2(u,v,t-2,ga,gb,gc,gd,n+1.0)+pz*ruvt2(u,v,t-1,ga,gb,gc,gd,n+1.0);
	}
	return r;
}
	
	
double gabcd(double* ga,double* gb,double* gc,double* gd){
	double p1=ga[6]+gb[6];
	double p2=gc[6]+gd[6];
	double pp=p1*p2/(p1+p2);
	double gv=0.0;
	double e1=0.0,e2=0.0;
	for (int u=0;u<(int) ga[3]+gb[3]+1.0;u+=1){
	for (int v=0;v<(int) ga[4]+gb[4]+1.0;v+=1){
	for (int t=0;t<(int) ga[5]+gb[5]+1.0;t+=1){
		e1=hermite(ga[0],(int) ga[3],ga[6],gb[0],(int) gb[3],gb[6],u)*
		   hermite(ga[1],(int) ga[4],ga[6],gb[1],(int) gb[4],gb[6],v)*
		   hermite(ga[2],(int) ga[5],ga[6],gb[2],(int) gb[5],gb[6],t);
		for (int u1=0;u1<(int) gc[3]+gd[3]+1.0;u1+=1){
		for (int v1=0;v1<(int) gc[4]+gd[4]+1.0;v1+=1){
		for (int t1=0;t1<(int) gc[5]+gd[5]+1.0;t1+=1){
			e2=hermite(gc[0],(int) gc[3],gc[6],gd[0],(int) gd[3],gd[6],u1)*
			   hermite(gc[1],(int) gc[4],gc[6],gd[1],(int) gd[4],gd[6],v1)*
		       hermite(gc[2],(int) gc[5],gc[6],gd[2],(int) gd[5],gd[6],t1);
			gv+=e1*e2*ruvt2(u+u1,v+v1,t+t1,ga,gb,gc,gd,0)*pow(-1,u1+v1+t1);
		}}}
	}}}
	return 2.0*pow(3.1415926,2.5)*gv/(p1*p2*pow(p1+p2,0.5));
}
	
void smat_pri(double* ba,double* s_mat,double* s_mat_dia,int len){
	
	for (int i=0;i<len;i++){
		s_mat_dia[i]=sab(ba+i*8,ba+i*8);
	}
	
	for (int i=0;i<len;i++){
	for (int j=0;j<len;j++){
		if (i>j){
			s_mat[i*len+j]=sab(ba+i*8,ba+j*8)/pow(s_mat_dia[i]*s_mat_dia[j],0.5);
			s_mat[j*len+i]=s_mat[i*len+j];
		}
		else if (i==j){
			s_mat[i*len+j]=1.0;
		}
	}}
}
	
void smat_con(double* ba,double* s_mat,double* s_mat_dia,int* contract_list,int co_len,int ba_len){
	
	for (int i=0;i<ba_len;i++){
		s_mat_dia[i]=sab(ba+i*8,ba+i*8);
	}
	
	double s=0.0;
	
	for (int i=0;i<co_len;i++){
	for (int j=0;j<co_len;j++){
		if (i>=j){
			s=0.0;
			for (int k=contract_list[i];k<contract_list[i+1];k++){
			for (int l=contract_list[j];l<contract_list[j+1];l++){
				s+=sab(ba+k*8,ba+l*8)/pow(s_mat_dia[k]*s_mat_dia[l],0.5)*ba[k*8+7]*ba[l*8+7];
			}}		
		s_mat[i*co_len+j]=s;
		s_mat[j*co_len+i]=s;
		}
	}}
}


void tmat_pri(double* ba,double* t_mat,double* s_mat_dia,int len){
	
	for (int i=0;i<len;i++){
	for (int j=0;j<len;j++){
		if (i>=j){
			t_mat[i*len+j]=tab(ba+i*8,ba+j*8)/pow(s_mat_dia[i]*s_mat_dia[j],0.5);
			t_mat[j*len+i]=t_mat[i*len+j];
		}
	}}
}
	
void tmat_con(double* ba,double* t_mat,double* s_mat_dia,int* contract_list,int len){
	double t=0.0;
	for (int i=0;i<len;i++){
	for (int j=0;j<len;j++){
		if (i>=j){
			t=0.0;
			for (int k=contract_list[i];k<contract_list[i+1];k++){
			for (int l=contract_list[j];l<contract_list[j+1];l++){
				t+=tab(ba+k*8,ba+l*8)/pow(s_mat_dia[k]*s_mat_dia[l],0.5)*ba[k*8+7]*ba[l*8+7];
			}}		
			t_mat[i*len+j]=t;
			t_mat[j*len+i]=t;
		}
	}}
}

void vmat_pri(double* ba,double* v_mat,double* s_mat_dia,double* nu,int ba_len,int nu_len){
	double v=0.0;
	for (int i=0;i<ba_len;i++){
	for (int j=0;j<ba_len;j++){
		if (i>=j){
			v=0.0;
			for (int k=0;k<nu_len;k++){
				v+=vab(ba+i*8,ba+j*8,nu+4*k);
			}

			v_mat[i*ba_len+j]=v/pow(s_mat_dia[i]*s_mat_dia[j],0.5);
			v_mat[j*ba_len+i]=v_mat[i*ba_len+j];
		}
	}}
}


void vmat_con(double* ba,double* v_mat,double* s_mat_dia,int* contract_list,double* nu,int con_len,int nu_len){
	double v1=0.0;
	double v2=0.0;
	for (int i=0;i<con_len;i++){
	for (int j=0;j<con_len;j++){
		if (i>=j){
			v1=0.0;
			for (int k=contract_list[i];k<contract_list[i+1];k++){
			for (int l=contract_list[j];l<contract_list[j+1];l++){
				v2=0.0;
				for (int kk=0;kk<nu_len;kk++){
					v2+=vab(ba+k*8,ba+l*8,nu+4*kk);
				}
				v1+=v2/pow(s_mat_dia[k]*s_mat_dia[l],0.5)*ba[k*8+7]*ba[l*8+7];
			}}
			v_mat[i*con_len+j]=v1;
			v_mat[j*con_len+i]=v1;
		}
	}}
}

void gab_pri(double* ba,double* g_mat,double* s_mat_dia,int len){
	double g=0.0;
	for (int i=0;i<len;i++){
	for (int j=0;j<len;j++){
	for (int k=0;k<len;k++){
	for (int l=0;l<len;l++){
	if (i>=j && k>=l && i>=k){
		g=gabcd(ba+i*8,ba+j*8,ba+k*8,ba+l*8)/
		  pow(s_mat_dia[i]*s_mat_dia[j]*s_mat_dia[k]*s_mat_dia[l],0.5);
		g_mat[i*len*len*len+j*len*len+k*len+l]=g;
		g_mat[j*len*len*len+i*len*len+k*len+l]=g;
		g_mat[j*len*len*len+i*len*len+l*len+k]=g;
		g_mat[i*len*len*len+j*len*len+l*len+k]=g;
		g_mat[k*len*len*len+l*len*len+i*len+j]=g;
		g_mat[l*len*len*len+k*len*len+i*len+j]=g;
		g_mat[l*len*len*len+k*len*len+j*len+i]=g;
		g_mat[k*len*len*len+l*len*len+j*len+i]=g;
	}}}}}
}

void gab_con(double* ba,double* g_mat,double* s_mat_dia,int* contract_list,int con_len){
	double g=0.0;
	for (int i=0;i<con_len;i++){
	for (int j=0;j<con_len;j++){
	for (int k=0;k<con_len;k++){
	for (int l=0;l<con_len;l++){
	if (i>=j && k>=l && i>=k){
		g=0.0;
		for (int ii=contract_list[i];ii<contract_list[i+1];ii++){
		for (int jj=contract_list[j];jj<contract_list[j+1];jj++){
		for (int kk=contract_list[k];kk<contract_list[k+1];kk++){
		for (int ll=contract_list[l];ll<contract_list[l+1];ll++){
			g+=gabcd(ba+ii*8,ba+jj*8,ba+kk*8,ba+ll*8)/
		  pow(s_mat_dia[ii]*s_mat_dia[jj]*s_mat_dia[kk]*s_mat_dia[ll],0.5)*
		  ba[ii*8+7]*ba[jj*8+7]*ba[kk*8+7]*ba[ll*8+7];
		}}}}
		g_mat[i*con_len*con_len*con_len+j*con_len*con_len+k*con_len+l]=g;
		g_mat[j*con_len*con_len*con_len+i*con_len*con_len+k*con_len+l]=g;
		g_mat[j*con_len*con_len*con_len+i*con_len*con_len+l*con_len+k]=g;
		g_mat[i*con_len*con_len*con_len+j*con_len*con_len+l*con_len+k]=g;
		g_mat[k*con_len*con_len*con_len+l*con_len*con_len+i*con_len+j]=g;
		g_mat[l*con_len*con_len*con_len+k*con_len*con_len+i*con_len+j]=g;
		g_mat[l*con_len*con_len*con_len+k*con_len*con_len+j*con_len+i]=g;
		g_mat[k*con_len*con_len*con_len+l*con_len*con_len+j*con_len+i]=g;
	}
	}}}}
}








