#coding:utf-8
import numpy as np
import math,copy,os
import  lebedev_laikov
from ctypes import cdll,c_double,POINTER,c_int

if os.path.exists("./grid_int.so"):
    cfun=cdll.LoadLibrary("./grid_int.so")
cfun.grid_int.restype=c_double


# becke's weight scheme  (10.1063/1.454033)
def niuw(ra,rb,rab):
    # eq(11)
    temp=(ra-rb)/rab
    # eq (20)
    temp=temp*1.5-temp**3*0.5
    temp=temp*1.5-temp**3*0.5
    temp=temp*1.5-temp**3*0.5
    # eq (21)
    temp=(1-temp)*0.5
    return(temp)

# get weight for atoms
def beck_w(point,ri,atom_i,atom_xyz,atom_distan):
    # eq (22)

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

# Partitioning parameters , table 4 in paper "A standard grid for density functional calculations"(10.1016/0009-2614(93)80125-9)
sg0_a_list=[[0.2500,0.5000,1.0000,4.5000],\
            [0.1667,0.5000,0.9000,3.5000],\
            [0.1000,0.4000,0.8000,2.5000]]

# atom radius and radial extent
radii={"H" :[1.0000, 0.6614 ],"He":[0.5882, 0.5291 ],\
       "Li":[3.0769, 2.4188 ],"Be":[2.0513, 1.8141 ],\
       "B" :[1.5385, 1.5874 ],"C" :[1.2308, 1.3795 ],\
       "N" :[1.0256, 1.3417 ],"O" :[0.8791, 1.2472 ],\
       "F" :[0.7692, 1.0771 ],"Ne":[0.6838, 1.0960 ],\
       "Na":[4.0909, 3.1369 ],"Mg":[3.1579, 2.6645 ],\
       "Al":[2.5714, 2.2866 ],"Si":[2.1687, 2.0976 ],\
       "P" :[1.8750, 2.0220 ],"S" :[1.6514, 1.9842 ],\
       "Cl":[1.4754, 1.9275 ],"Ar":[1.3333, 2.0031 ]}

# generate sphere grid point with lebedev_laikov scheme 
# V. I. Lebedev, and D. N. Laikov, “A quadrature formula for the sphere of the 131st 
# algebraic order of accuracy,” Doklady Mathematics, 59 (3), 477-481 (1999). http://rad.chem.msu.ru/~laikov/ru/DAN_366_741.pdf
# using the code from Rufflewind (https://github.com/Rufflewind/lebedev_laikov)

# the generation of grid points for atoms using the SG0 method 
# reference : "SG-0: A Small Standard Grid for DFT Quadrature on Large Systems“  10.1002/jcc.20383
# and " a standard grid for density functional calculations"  10.1016/0009-2614(93)80125-9
def xyzw2(atom_i,atom_xyz,atom_distan,w=1):
    xyzw_list=[]
    # how many radial points used 
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
        # the weight of every point in integral
        # eq (7) in paper "A standard grid for density functional calculations"(10.1016/0009-2614(93)80125-9)
        ri=p*xi**2/(1-xi)**2
        # eq (6)
        wi=2*p**3/(radial_n+1.0)*xi**5/(1-xi)**7

        # select grid point acoording to the radial Partitioning parameters and atom radii
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
                # the total weight of grid point
                temp[j,3]*=beck_w(temp[j,:],ri,atom_i,atom_xyz,atom_distan)*wi
                temp[j,0]=atom_xyz_i[1]+temp[j,0]*ri
                temp[j,1]=atom_xyz_i[2]+temp[j,1]*ri
                temp[j,2]=atom_xyz_i[3]+temp[j,2]*ri

            xyzw_list+=temp.tolist()


    return(np.array(xyzw_list))

# generate grid point and weight for all atoms
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

# calculate the electron density of grid points for basis function
def density_pot(ga,gb,pot):
    x1,y1,z1,i1,j1,k1,a1,_=ga
    x2,y2,z2,i2,j2,k2,a2,_=gb
    x1,x2=pot[0]-x1,pot[0]-x2
    y1,y2=pot[1]-y1,pot[1]-y2
    z1,z2=pot[2]-z1,pot[2]-z2

    density=x1**i1*y1**j1*z1**k1*math.exp(-(x1**2+y1**2+z1**2)*a1)
    density*=x2**i2*y2**j2*z2**k2*math.exp(-(x2**2+y2**2+z2**2)*a2)
    return(density)

# calculate the electron density at one point from all basis function
def density_pot_t(ba,p_mat,s_mat_dia,pot):
    density=0.0
    for i in range(len(ba)):
        for j in range(len(ba)):
            if i>j:
                density+=2*density_pot(ba[i],ba[j],pot)/(s_mat_dia[i]*s_mat_dia[j])**0.5*p_mat[i,j]
            elif i==j:
                density+=  density_pot(ba[i],ba[i],pot)/s_mat_dia[i]*p_mat[i,i]

    return(density)

# calculate the electron density at one point for contracted basis set
def density_pot_list_t_contra(ba,p_mat,s_mat_dia,pot):
    density=0.0
    for i in range(len(ba)):
        for j in range(len(ba)):
            if i>j:
                density+=2*density_pot(ba[i],ba[j],pot)/(s_mat_dia[i]*s_mat_dia[j])**0.5*p_mat[i,j]*ba[i][-1]*ba[j][-1]
            elif i==j:
                density+=  density_pot(ba[i],ba[i],pot)/s_mat_dia[i]*p_mat[i,i]*ba[i][-1]**2

    return(density)


# calculate the electron density for primitive basis set
def density_pot_list(ba,xyzw_list,s_mat_dia,p_mat):
    density_list=[]
    for i in range(len(xyzw_list)):
        temp1=[]
        xyzw1=xyzw_list[i]
        for j in range(len(xyzw1)):
            density=density_pot_t(ba,p_mat,s_mat_dia,xyzw1[j][:3])
            temp1.append(density)
        density_list.append(temp1)
    return(density_list)

# calculate the electron density for contracted basis set
def density_pot_list_contra(ba,xyzw_list,s_mat_dia,p_mat):
    density_list=[]
    for i in range(len(xyzw_list)):
        temp1=[]
        xyzw1=xyzw_list[i]
        for j in range(len(xyzw1)):
            density=density_pot_list_t_contra(ba,p_mat,s_mat_dia,xyzw1[j][:3])
            temp1.append(density)
        density_list.append(temp1)
    return(density_list)

# density integral of grid point 
def grid_int(ga, gb, xyzw, density_list):
    density=0.0;
    for n,i in enumerate(xyzw):
        density+=density_pot(ga,gb, i[:3])*i[3]*density_list[n]**(1.0/3.0)
    return density;

# integral of grid point 
def xc_int(ba,ii,jj,s_mat_dia,p_mat,xyzw_list,density_list,p=1):

    sum1=0.0    
    ba=np.array(ba)
    for k in range(len(density_list)):
        sum1+=grid_int(ba[ii],ba[jj],xyzw_list[k], density_list[k])
        
    return(sum1*math.pi*4/(s_mat_dia[ii]*s_mat_dia[jj])**0.5)

#
#  c++ lib version functions 
#

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



