// #include "aomo.h"
#include "math.h"
// #include <iostream>
using namespace std;

extern "C" void aomo1(double* ao,double* mo,double* c,int n);
extern "C" double density_pot(double* ga, double* gb, double x,double y,double z );
extern "C" double density_pot_t(double* ba, double* p_mat,double* s_mat_dia,int basis_num, double x,double y,double z );
extern "C" double density_pot_t_contra(double* ba, double* p_mat,double* s_mat_dia,int basis_num, double x,double y,double z );
extern "C" double grid_int(double* ga,double* gb, double* xyzw, double* density_list,int len);
extern "C" void density_pot_list(double* ba,double* p_mat, double* s_mat_dia, double* xyzw, double* density_list,int basis_num,int len);
extern "C" void density_pot_list_contra(double* ba,double* p_mat, double* s_mat_dia, double* xyzw, double* density_list,int basis_num,int len);



void aomo1(double* ao,double* mo,double* c,int n){
    int ao1;
    int c1;
    int c2;
    int c3;
    int c4;
    int mo1;
    double temp;
	for (int i=0;i<n;i++){
	for (int j=0;j<n;j++){
	for (int k=0;k<n;k++){
	for (int l=0;l<n;l++){
	if (i>=j && k>=l && i>=k){
		temp=0.0;
		for (int ii=0;ii<n;ii++){
		for (int jj=0;jj<n;jj++){
		for (int kk=0;kk<n;kk++){
		for (int ll=0;ll<n;ll++){
	        ao1=ii*n*n*n+jj*n*n+kk*n+ll;
	        c1=ii*n+i;
	        c2=jj*n+j;
	        c3=kk*n+k;
	        c4=ll*n+l;
			temp+=ao[ao1]*c[c1]*c[c2]*c[c3]*c[c4];}}}}
										

		mo1=i*n*n*n+j*n*n+k*n+l;
		mo[mo1]=temp;
		mo1=j*n*n*n+i*n*n+k*n+l;
		mo[mo1]=temp;
		mo1=j*n*n*n+i*n*n+l*n+k;
		mo[mo1]=temp;
		mo1=i*n*n*n+j*n*n+l*n+k;
		mo[mo1]=temp;
		mo1=k*n*n*n+l*n*n+i*n+j;
		mo[mo1]=temp;
		mo1=l*n*n*n+k*n*n+i*n+j;
		mo[mo1]=temp;
		mo1=l*n*n*n+k*n*n+j*n+i;
		mo[mo1]=temp;
		mo1=k*n*n*n+l*n*n+j*n+i;
		mo[mo1]=temp;


	}}}}}
}



double density_pot(double* ga, double* gb, double x,double y,double z ){
	double density=0.0;
	double x1=x-ga[0];
	double y1=y-ga[1];
	double z1=z-ga[2];
	double x2=x-gb[0];
	double y2=y-gb[1];
	double z2=z-gb[2];
	density =pow(x1,ga[3])*pow(y1,ga[4])*pow(z1,ga[5])*exp(-(x1*x1+y1*y1+z1*z1)*ga[6]);
	density*=pow(x2,gb[3])*pow(y2,gb[4])*pow(z2,gb[5])*exp(-(x2*x2+y2*y2+z2*z2)*gb[6]);
	return density;
}


double density_pot_t(double* ba, double* p_mat,double* s_mat_dia,int basis_num, double x,double y,double z ){
	double density=0.0;
	for (int i=0;i<basis_num;i++){
	for (int j=0;j<basis_num;j++){
		double* ga=ba+i*8;
		double* gb=ba+j*8;
		if (i>j){
			density += 2.0*density_pot(ga, gb,  x, y, z )/pow(s_mat_dia[i]*s_mat_dia[j],0.5)*p_mat[i*basis_num+j];
		}
		if (i==j){
			density += density_pot(ga, gb, x, y, z )/s_mat_dia[i]*p_mat[i*basis_num+j];
		}}}
	return density;
}

double density_pot_t_contra(double* ba, double* p_mat,double* s_mat_dia,int basis_num, double x,double y,double z ){
	double density=0.0;
	for (int i=0;i<basis_num;i++){
	for (int j=0;j<basis_num;j++){
		double* ga=ba+i*8;
		double* gb=ba+j*8;
		if (i>j){
			density += 2.0*density_pot(ga, gb,  x, y, z )/pow(s_mat_dia[i]*s_mat_dia[j],0.5)*p_mat[i*basis_num+j]*ba[i*8+7]*ba[j*8+7];
		}
		if (i==j){
			density += density_pot(ga, gb, x, y, z )/s_mat_dia[i]*p_mat[i*basis_num+j]*ba[i*8+7]*ba[i*8+7];
		}}}
	return density;
}


double grid_int(double* ga,double* gb, double* xyzw, double* density_list,int len){
	double density=0.0;
	double x,y,z;
	for (int i=0;i<len;i++){
		x=xyzw[i*4];
		y=xyzw[i*4+1];
		z=xyzw[i*4+2];
		density+=density_pot(ga,gb,x,y,z)*pow(density_list[i],1.0/3.0)*xyzw[i*4+3];
		}
	return density;
}



void density_pot_list(double* ba,double* p_mat, double* s_mat_dia, double* xyzw, double* density_list,int basis_num,int len){
	double x,y,z;
	for (int i=0;i<len;i++){
		x=xyzw[i*4];
		y=xyzw[i*4+1];
		z=xyzw[i*4+2];
		density_list[i]=density_pot_t(ba, p_mat,s_mat_dia,basis_num, x, y, z );
		}
}


void density_pot_list_contra(double* ba,double* p_mat, double* s_mat_dia, double* xyzw, double* density_list,int basis_num,int len){
	double x,y,z;
	for (int i=0;i<len;i++){
		x=xyzw[i*4];
		y=xyzw[i*4+1];
		z=xyzw[i*4+2];
		density_list[i]=density_pot_t_contra(ba, p_mat,s_mat_dia,basis_num, x, y, z );
		}
}
