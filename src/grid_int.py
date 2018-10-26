#coding:utf-8
import numpy as np
import math,copy
import  lebedev_laikov
from ctypes import cdll,c_double,POINTER,c_int


cfun=cdll.LoadLibrary("./grid_int.so")
cfun.grid_int.restype=c_double



def niuw(ra,rb,rab):
	temp=(ra-rb)/rab
	temp=temp*1.5-temp**3*0.5
	temp=temp*1.5-temp**3*0.5
	temp=temp*1.5-temp**3*0.5
	temp=(1-temp)*0.5
	return(temp)


def beck_w(point,ri,atom_i,atom_xyz,atom_distan):
	

	atom_i_w1=np.where(atom_distan[atom_i]<=20)[0]
	atom_i_w=np.where(atom_i_w1==atom_i)[0][0]
	atom_xyz_w=[atom_xyz[i]  for i in atom_i_w1]
	atom_num_w=len(atom_xyz_w)

	pot_distan=[]
	xa=point[0]*ri+atom_xyz_w[atom_i_w][1]
	ya=point[1]*ri+atom_xyz_w[atom_i_w][2]
	za=point[2]*ri+atom_xyz_w[atom_i_w][3]
	for i in range(atom_num_w):
		temp=((xa-atom_xyz_w[i][1])**2+(ya-atom_xyz_w[i][2])**2+(za-atom_xyz_w[i][3])**2)**0.5
		pot_distan.append(temp)
	
	atom_w=np.ones((atom_num_w,atom_num_w))
	for i in range(atom_num_w):
		for j in range(atom_num_w):
			if i>j:
				ra=pot_distan[i]
				rb=pot_distan[j]
				rab=atom_distan[atom_i_w1[i],atom_i_w1[j]]
				niu_a=niuw(ra,rb,rab)
				atom_w[i,j]=niu_a
				atom_w[j,i]=1-niu_a

	atom_w=np.prod(atom_w,axis=1).T
	w_sum=np.sum(atom_w)

	return(atom_w[atom_i_w]/w_sum)




#                a1     a2     a3      a4
#        H-He  0.2500 0.5000 1.0000  4.5000
#       Li-Ne  0.1667 0.5000 0.9000  3.5000
#       Na-Ar  0.1000 0.4000 0.8000  2.5000

sg0_a_list=[[0.2500,0.5000,1.0000,4.5000],\
            [0.1667,0.5000,0.9000,3.5000],\
            [0.1000,0.4000,0.8000,2.5000]]


radii={"H" :[1.0000, 0.6614 ],"He":[0.5882, 0.5291 ],\
       "Li":[3.0769, 2.4188 ],"Be":[2.0513, 1.8141 ],\
       "B" :[1.5385, 1.5874 ],"C" :[1.2308, 1.3795 ],\
       "N" :[1.0256, 1.3417 ],"O" :[0.8791, 1.2472 ],\
       "F" :[0.7692, 1.0771 ],"Ne":[0.6838, 1.0960 ],\
       "Na":[4.0909, 3.1369 ],"Mg":[3.1579, 2.6645 ],\
       "Al":[2.5714, 2.2866 ],"Si":[2.1687, 2.0976 ],\
       "P" :[1.8750, 2.0220 ],"S" :[1.6514, 1.9842 ],\
       "Cl":[1.4754, 1.9275 ],"Ar":[1.3333, 2.0031 ]}


def xyzw2(atom_i,atom_xyz,atom_distan,w=1):
	xyzw_list=[]
	radial_n=50

	atom_xyz_i=atom_xyz[atom_i]
	atom_r=radii[atom_xyz_i[0]][0]
	p=radii[atom_xyz_i[0]][1]
	# print(p)

	if atom_xyz_i[0] in ["H","He"]:
		a_list=sg0_a_list[0]
	elif atom_xyz_i[0] in ["Li","Be","B","C","N","O","F","Ne"]:
		a_list=sg0_a_list[1]
	else:
		a_list=sg0_a_list[2]


	for i in range(1,radial_n+1):
		xi=i/(radial_n+1.0)
		ri=p*xi**2/(1-xi)**2
		wi=2*p**3/(radial_n+1.0)*xi**5/(1-xi)**7


		if ri<a_list[0]*atom_r:
			pot_num=6
		elif ri>=a_list[0]*atom_r and ri<a_list[1]*atom_r:
			pot_num=38
		elif ri>=a_list[1]*atom_r and ri<a_list[2]*atom_r:
			pot_num=86
		elif ri>=a_list[2]*atom_r and ri<a_list[3]*atom_r:
			pot_num=194
		else:
			pot_num=86

		if w==1 and ri<20:
			temp=lebedev_laikov.ld(pot_num).T
			for j in range(pot_num):
				temp[j,3]*=beck_w(temp[j,:],ri,atom_i,atom_xyz,atom_distan)*wi
				temp[j,0]=atom_xyz_i[1]+temp[j,0]*ri
				temp[j,1]=atom_xyz_i[2]+temp[j,1]*ri
				temp[j,2]=atom_xyz_i[3]+temp[j,2]*ri

			xyzw_list+=temp.tolist()


	return(np.array(xyzw_list))








def xyzw_list_gen(atom_xyz):
	atom_num=len(atom_xyz)
	atom_distan=np.zeros((atom_num,atom_num))
	for i in range(atom_num):
		for j in range(atom_num):
			if i>j:
				temp=(atom_xyz[i][1]-atom_xyz[j][1])**2+(atom_xyz[i][2]-atom_xyz[j][2])**2+\
				(atom_xyz[i][3]-atom_xyz[j][3])**2
				atom_distan[i,j]=temp**0.5
				atom_distan[j,i]=temp**0.5

	xyzw_list=[]
	for i in range(atom_num):
		xyzw_list.append(xyzw2(i,atom_xyz,atom_distan))

	return(xyzw_list)
	

def density_pot_list_c(ba,xyzw_list,s_mat_dia,p_mat):
	density_list=[]
	ba=np.array(ba)

	basis_num=c_int(len(ba))

	s_mat_dia1=s_mat_dia.ctypes.data_as(POINTER(c_double))
	p_mat1=p_mat.ctypes.data_as(POINTER(c_double))

	ba1=ba.ctypes.data_as(POINTER(c_double))


	for i in range(len(xyzw_list)):
		len1=len(xyzw_list[i])
		temp1=np.zeros((len1))
		xyzw1=xyzw_list[i].ctypes.data_as(POINTER(c_double))
		temp21=temp1.ctypes.data_as(POINTER(c_double))
		cfun.density_pot_list(ba1,p_mat1,s_mat_dia1,xyzw1,temp21,basis_num,c_int(len1))

		density_list.append(temp1)
	return(density_list)



def density_pot_list_c_contra(ba,xyzw_list,s_mat_dia,p_mat):
	density_list=[]
	ba=np.array(ba)

	basis_num=c_int(len(ba))

	s_mat_dia1=s_mat_dia.ctypes.data_as(POINTER(c_double))
	p_mat1=p_mat.ctypes.data_as(POINTER(c_double))

	ba1=ba.ctypes.data_as(POINTER(c_double))


	for i in range(len(xyzw_list)):
		len1=len(xyzw_list[i])
		temp1=np.zeros((len1))
		xyzw1=xyzw_list[i].ctypes.data_as(POINTER(c_double))
		temp21=temp1.ctypes.data_as(POINTER(c_double))
		cfun.density_pot_list_contra(ba1,p_mat1,s_mat_dia1,xyzw1,temp21,basis_num,c_int(len1))

		density_list.append(temp1)
	return(density_list)



def xc_int_c(ba,ii,jj,s_mat_dia,p_mat,xyzw_list,density_list,p=1):

	sum1=0.0	
	ba=np.array(ba)
	ga=ba[ii].ctypes.data_as(POINTER(c_double))
	gb=ba[jj].ctypes.data_as(POINTER(c_double))

	for k in range(len(density_list)):
		len1=len(xyzw_list[k])
		xyzw1=xyzw_list[k].ctypes.data_as(POINTER(c_double))
		density_list1=density_list[k].ctypes.data_as(POINTER(c_double))
		sum1+=cfun.grid_int(ga,gb,xyzw1,density_list1,c_int(len1))
		
	return(sum1*math.pi*4/(s_mat_dia[ii]*s_mat_dia[jj])**0.5)



