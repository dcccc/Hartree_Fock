#coding=utf-8
import math,os
import numpy as np
import scipy.special as spec

# calculate the integral of gauss function g(x1)*g(x2),  overlap integral
# using obara-saikai scheme 
def s_ab(x1,i1,a1,x2,i2,a2):

    # eq (9.2.10)-(9.2.20) in reference book molecular electronic structure theory, and all equations
    # comes from the same book
    p=a1+a2
    niu=a1*a2/p
    px=(a1*x1+a2*x2)/p
    xab=x1-x2
    kab=math.exp(-niu*xab*xab)
    xa=px-x1
    xb=px-x2

    if i1<0 or i2<0:
        s=0.0
    elif i1==0 and i2==0 :
        # eq (9.3.10)  
        s=(math.pi/p)**0.5*kab
    elif i1>0 and i2==0:
        # eq (9.3.8)
        s=xa*s_ab(x1,i1-1,a1,x2,i2,a2)+0.5/p*((i1-1)*s_ab(x1,i1-2,a1,x2,i2,a2))
    else:
        # eq (9.3.9)
        s=xb*s_ab(x1,i1,a1,x2,i2-1,a2)+0.5/p*(i1*s_ab(x1,i1-1,a1,x2,i2-1,a2)+
            (i2-1)*s_ab(x1,i1,a1,x2,i2-2,a2))

    return(s)

# calculate the integration of basis function g(x1,y1,z1)*g(x2,y2,z2)
def sab(ga,gb):
    x1,y1,z1,i1,j1,k1,a1,c1=ga
    x2,y2,z2,i2,j2,k2,a2,c2=gb
    s1=s_ab(x1,i1,a1,x2,i2,a2)
    s2=s_ab(y1,j1,a1,y2,j2,a2)
    s3=s_ab(z1,k1,a1,z2,k2,a2)
    # print(s1,s2,s3)
    return(s1*s2*s3)

# S matrix (contracted basis set)
def get_S_mat_con(ba):
    basis_num = len(ba)
    s_mat=np.zeros((basis_num,basis_num))
    for i in range(basis_num):
        for j in range(basis_num):
            temp=0.0
            for ii in ba[i]:
                for jj in ba[j]:
                    temp+=sab(ii,jj)*ii[-1]*jj[-1] \
                          /(sab(ii,ii)*sab(jj,jj))**0.5
            s_mat[i,j]=temp
            
    return(s_mat)

# S matrix (primitive basis set)
def get_S_mat_pri(ba):
    basis_num = len(ba)
    s_mat_dia=[]
    s_mat=np.zeros((basis_num,basis_num))
    for i in range(basis_num):
        for j in range(basis_num):
            if i>=j:
                s_mat[i,j]=sab(ba[i],ba[j])
                # s_mat[j,i]=s_mat[i,j]
            if i==j:
                s_mat_dia.append(s_mat[i,j])

    for i in range(basis_num):
        for j in range(basis_num):
            if i>=j:
                s_mat[i,j]=s_mat[i,j]/(s_mat_dia[i]*s_mat_dia[j])**0.5
                s_mat[j,i]=s_mat[i,j]
    return(s_mat, s_mat_dia)

# kinetic integral of gauss function
def t_ab(x1,i1,a1,x2,i2,a2):
    p=a1+a2
    niu=a1*a2/p
    px=(a1*x1+a2*x2)/p
    xa=px-x1
    xb=px-x2

    if i1<0 or i2<0:
        t=0.0
    elif i1==0 and i2==0:
        # eq (9.3.42)
        t=(a1-2*a1**2*(xa**2+0.5/p))*s_ab(x1,0,a1,x2,0,a2)
    elif i1>0 and i2==0:
        # eq (9.3.40)
        t=xa*t_ab(x1,i1-1,a1,x2,0,a2)+\
        0.5/p*((i1-1)*t_ab(x1,i1-2,a1,x2,0,a2))+\
        a2/p*(2.0*a1*s_ab(x1,i1,a1,x2,0,a2)-(i1-1)*s_ab(x1,i1-2,a1,x2,0,a2))
    else:
        # eq (9.3.41)
        t=xb*t_ab(x1,i1,a1,x2,i2-1,a2)+\
        0.5/p*(i1*t_ab(x1,i1-1,a1,x2,i2-1,a2)+(i2-1)*t_ab(x1,i1,a1,x2,i2-2,a2))+\
        a1/p*(2.0*a2*s_ab(x1,i1,a1,x2,i2,a2)-(i2-1)*s_ab(x1,i1,a1,x2,i2-2,a2))

    return(t)


def tab(ga,gb):
    x1,y1,z1,i1,j1,k1,a1,c1=ga
    x2,y2,z2,i2,j2,k2,a2,c2=gb
    # eq (9.3.38)
    t1=t_ab(x1,i1,a1,x2,i2,a2)*s_ab(y1,j1,a1,y2,j2,a2)*s_ab(z1,k1,a1,z2,k2,a2)
    t2=s_ab(x1,i1,a1,x2,i2,a2)*t_ab(y1,j1,a1,y2,j2,a2)*s_ab(z1,k1,a1,z2,k2,a2)
    t3=s_ab(x1,i1,a1,x2,i2,a2)*s_ab(y1,j1,a1,y2,j2,a2)*t_ab(z1,k1,a1,z2,k2,a2)
    return(t1+t2+t3)

# T matrix (contracted basis set)
def get_T_mat_con(ba):
    basis_num = len(ba)
    t_mat=np.zeros((basis_num,basis_num))
    for i in range(basis_num):
        for j in range(basis_num):
            temp=0.0
            for ii in ba[i]:
                for jj in ba[j]:
                    temp+=tab(ii,jj)*ii[-1]*jj[-1]\
                    /(sab(ii,ii)*sab(jj,jj))**0.5
            t_mat[i,j]=temp
            
    return(t_mat)

# S matrix (primitive basis set)
def get_T_mat_pri(ba, s_mat_dia):
    basis_num = len(ba)
    t_mat=np.zeros((basis_num,basis_num))
    for i in range(basis_num):
        for j in range(basis_num):
            if i>=j:
                t_mat[i,j]=tab(ba[i],ba[j])/(s_mat_dia[i]*s_mat_dia[j])**0.5
                t_mat[j,i]=t_mat[i,j]
            
    return(t_mat)

# integral of hermite gaussian 
def hermite(x1,i1,a1,x2,i2,a2,t):
    p=a1+a2

    if t<0 or t > i1+i2 or i1<0 or i2<0:
        # eq (9.5.5)
        e=0.0
    elif t==0:
        # eq (9.5.8)
        e=(p/math.pi)**0.5*s_ab(x1,i1,a1,x2,i2,a2)
    else:
        # eq (9.5.17)
        e=(i1*hermite(x1,i1-1,a1,x2,i2,a2,t-1)+i2*hermite(x1,i1,a1,x2,i2-1,a2,t-1))/2.0/p/t
    return(e)

# boys function
def boys(n,x):

    if x==0.0:
        # eq (9.8.1)
        boy=1.0/(2*n+1)
    elif n==0 and x !=0.0:
        # eq (9.8.1)
        boy=(math.pi)**0.5*math.erf(x**0.5)/(x**0.5)/2
    elif x<0.18:
        boy=0.0
        # eq (9.8.11), calculated numeracially directly when x is small 
        for k in range(7):
            boy+=(-x)**k/(math.factorial(k)*(2.0*k+2.0*n+1))

    else:
        # eq (9.8.21), when x is large
        boy=0.5/x**(n+0.5)*spec.gamma(n+0.5)*spec.gammainc(n+0.5,x)
    return(boy)


def ruvt(u,v,t,ga,gb,n,nu):

    # eq (9.2.10)-(9.2.15)
    x1,y1,z1,i1,j1,k1,a1,c1=ga
    x2,y2,z2,i2,j2,k2,a2,c2=gb
    p=a1+a2
    px=(a1*x1+a2*x2)/p
    py=(a1*y1+a2*y2)/p
    pz=(a1*z1+a2*z2)/p
    xa=px-nu[0]
    ya=py-nu[1]
    za=pz-nu[2]

    rpc2=xa**2+ya**2+za**2

    if u<0 or v<0 or t<0 or n<0:
        r=0.0
    elif u==0 and v==0 and t==0 :
        # eq (9.9.14)
        r=(-2*p)**n*boys(n,p*rpc2)
    elif u>0 and v==0 and t==0:
        # eq (9.9.19)
        r=(u-1)*ruvt(u-2,0,0,ga,gb,n+1,nu)+xa*ruvt(u-1,0,0,ga,gb,n+1,nu)
    elif v>0 and t==0:
        # eq (9.9.20)
        r=(v-1)*ruvt(u,v-2,0,ga,gb,n+1,nu)+ya*ruvt(u,v-1,0,ga,gb,n+1,nu)
    else:
        # eq (9.9.18)
        r=(t-1)*ruvt(u,v,t-2,ga,gb,n+1,nu)+za*ruvt(u,v,t-1,ga,gb,n+1,nu)

    return(r)


# Mcmurchie-Davidson scheme for coulomb itegral
def vab(ga,gb,nu):
    x1,y1,z1,i1,j1,k1,a1,c1=ga
    x2,y2,z2,i2,j2,k2,a2,c2=gb
    p=a1+a2
    vv=0.0
    for u in range(i1+i2+1):
        for v in range(j1+j2+1):
            for t in range(k1+k2+1):
                # eq (9.9.40)
                vv+=hermite(x1,i1,a1,x2,i2,a2,u)*hermite(y1,j1,a1,y2,j2,a2,v) \
                *hermite(z1,k1,a1,z2,k2,a2,t)*ruvt(u,v,t,ga,gb,0,nu)

    return(math.pi*vv*2*nu[3]/p)


def get_V_mat_con(ba, nu):
    basis_num = len(ba)
    v_mat=np.zeros((basis_num,basis_num))
    for i in range(basis_num):
        for j in range(basis_num):
            temp=0.0
            for k in nu:
                for ii in ba[i]:
                    for jj in ba[j]:
                        temp+=vab(ii,jj,k)*ii[-1]*jj[-1]\
                        /(sab(ii,ii)*sab(jj,jj))**0.5
        
            v_mat[i,j]=temp
    return(v_mat)

def get_V_mat_pri(ba, nu,s_mat_dia):
    basis_num = len(ba)
    v_mat=np.zeros((basis_num,basis_num))
    for i in range(basis_num):
        for j in range(basis_num):
            if i>=j:
                temp=0.0
                for k in nu:
                    temp+=vab(ba[i],ba[j],k)
                v_mat[i,j]=temp/(s_mat_dia[i]*s_mat_dia[j])**0.5
                v_mat[j,i]=v_mat[i,j]
    return(v_mat)


def ruvt_2(u,v,t,ga,gb,gc,gd,n):

    x1,y1,z1,i1,j1,k1,a1,c1=ga
    x2,y2,z2,i2,j2,k2,a2,c2=gb
    x3,y3,z3,i3,j3,k3,a3,c3=gc
    x4,y4,z4,i4,j4,k4,a4,c4=gd
    p1=a1+a2
    p2=a3+a4
    pp=p1*p2/(p1+p2)
    px=(a1*x1+a2*x2)/p1-(a3*x3+a4*x4)/p2
    py=(a1*y1+a2*y2)/p1-(a3*y3+a4*y4)/p2
    pz=(a1*z1+a2*z2)/p1-(a3*z3+a4*z4)/p2

    rpc2=px**2+py**2+pz**2

    if u<0 or v<0 or t<0 or n<0:
        r=0.0
    elif u==0 and v==0 and t==0 :
        r=(-2*pp)**n*boys(n,pp*rpc2)

    elif u>0 :
        r=(u-1)*ruvt_2(u-2,v,t,ga,gb,gc,gd,n+1)+px*ruvt_2(u-1,v,t,ga,gb,gc,gd,n+1)
    elif v>0:
        r=(v-1)*ruvt_2(u,v-2,t,ga,gb,gc,gd,n+1)+py*ruvt_2(u,v-1,t,ga,gb,gc,gd,n+1)
    elif t>0:
        r=(t-1)*ruvt_2(u,v,t-2,ga,gb,gc,gd,n+1)+pz*ruvt_2(u,v,t-1,ga,gb,gc,gd,n+1)

    return(r)

def gabcd(ga,gb,gc,gd):
    x1,y1,z1,i1,j1,k1,a1,c1=ga
    x2,y2,z2,i2,j2,k2,a2,c2=gb
    x3,y3,z3,i3,j3,k3,a3,c3=gc
    x4,y4,z4,i4,j4,k4,a4,c4=gd
    p1=a1+a2
    p2=a3+a4
    pp=p1*p2/(p1+p2)
    gv=0.0
    # eq (9.9.47)-(9.9.48)
    for u in range(i1+i2+1):
        for v in range(j1+j2+1):
            for t in range(k1+k2+1):
                e1=hermite(x1,i1,a1,x2,i2,a2,u)*hermite(y1,j1,a1,y2,j2,a2,v) \
                *hermite(z1,k1,a1,z2,k2,a2,t)
                for u1 in range(i3+i4+1):
                    for v1 in range(j3+j4+1):
                        for t1 in range(k3+k4+1):
                            e2=hermite(x3,i3,a3,x4,i4,a4,u1)*hermite(y3,j3,a3,y4,j4,a4,v1) \
                            *hermite(z3,k3,a3,z4,k4,a4,t1)
                            gv+=e1*e2*ruvt_2(u+u1,v+v1,t+t1,ga,gb,gc,gd,0)*(-1)**(u1+v1+t1)
    return(2*(math.pi)**2.5*gv/(p1*p2*(p1+p2)**0.5))

def get_Gabcd_mat_con(ba):
    basis_num = len(ba)
    g_mat=np.zeros((basis_num,basis_num,basis_num,basis_num))
    for i in range(basis_num):
        for j in range(basis_num):
            for k in range(basis_num):
                for l in range(basis_num):                
                    if i>=j and k>=l and i>=k  :
                        temp=0.0
                        for ii in ba[i]:
                            for jj in ba[j]:
                                for kk in ba[k]:
                                    for ll in ba[l]:
                            
                                        temp+=gabcd(ii,jj,kk,ll)*ii[-1]*jj[-1]*kk[-1]*ll[-1]/ \
                                        (sab(ii,ii)*sab(jj,jj)*sab(kk,kk)*sab(ll,ll))**0.5
                        
                        g_mat[i,j,k,l]=temp
                        g_mat[i,j,l,k]=temp
                        g_mat[j,i,k,l]=temp
                        g_mat[j,i,l,k]=temp
                        g_mat[k,l,i,j]=temp
                        g_mat[k,l,j,i]=temp
                        g_mat[l,k,i,j]=temp
                        g_mat[l,k,j,i]=temp
    return(g_mat)

def get_Gabcd_mat_pri(ba,s_mat_dia):
    basis_num = len(ba)
    g_mat=np.zeros((basis_num,basis_num,basis_num,basis_num))
    for i in range(basis_num):
        for j in range(basis_num):
            for k in range(basis_num):
                for l in range(basis_num):
                    if i>=j and k>=l and i>=k :
                        
                        temp=gabcd(ba[i],ba[j],ba[k],ba[l])
                        temp=temp/(s_mat_dia[i]*s_mat_dia[j]*s_mat_dia[k]*s_mat_dia[l])**0.5

                        g_mat[i,j,k,l]=temp
                        g_mat[i,j,l,k]=temp
                        g_mat[j,i,k,l]=temp
                        g_mat[j,i,l,k]=temp
                        g_mat[k,l,i,j]=temp
                        g_mat[k,l,j,i]=temp
                        g_mat[l,k,i,j]=temp
                        g_mat[l,k,j,i]=temp
    return(g_mat)


def gabdc_index(i,j,k,l):
    if i<=j:
        g_index1="%d-%d" %(j,i)
    else:
        g_index1="%d-%d" %(i,j)
    if k<=l:
        g_index2="%d-%d" %(l,k)
    else:
        g_index2="%d-%d" %(k,l)
    if (i>=k and i>=l) or (j>i and j>=k and j>=l) :
        g_index=g_index1+"-"+g_index2
    else:
        g_index=g_index2+"-"+g_index1
    return(g_index)

# density matrix
def pab(C_mat,n_electron):
    shape=C_mat.shape
    p_mat=np.zeros(shape)
    for i in range(shape[0]):
        for j in range(shape[0]):
            temp=0.0
            # print(shape)
            for k in range(n_electron//2):
                temp+=C_mat[i,k]*C_mat[j,k]
            p_mat[i,j]=2*temp
    return(p_mat)

# two electron coulumb matrix for coulomb part
def g_ab_j(p_mat,gabcd_mat,i,j,basis_num):
    gv=0.0
    # 
    for k in range(basis_num):
        for l in range(basis_num):
            gv+=p_mat[k,l]*(gabcd_mat[i,j,k,l])  
    return(gv)

# two electron coulumb matrix for coulomb and exchange part
def g_ab_jk(p_mat,gabcd_mat,i,j,basis_num):
    gv=0.0
    # eq (10.3.21)
    for k in range(basis_num):
        for l in range(basis_num):
            gv+=p_mat[k,l]*(gabcd_mat[i,j,k,l]-gabcd_mat[i,l,k,j]/2)
    return(gv)



def c_normarize(c_mat,e_mat,s_mat,basis_num):

    ce_mat=np.vstack((e_mat,c_mat)).T.tolist()
    ce_mat=np.array(sorted(ce_mat,key=lambda x: x[0]))
    e_mat,c_mat=ce_mat.T[0,:],ce_mat.T[1:,:]
    for i in range(basis_num):
        N=0.0
        for j in range(basis_num):
            for k in range(basis_num):
                N+=c_mat[j,i]*c_mat[k,i]*s_mat[j,k]
        c_mat[:,i]=c_mat[:,i]/N**0.5
    return(e_mat,c_mat)

def all_energy(p_mat,F_mat,nu,basis_num):
    v_nn=0.0
    e_e=0.0
    # eq (10.4.15) (10.4.20) 
    for i in range(len(nu)):
        for j in range(len(nu)):
            if i>j:
                v_nn+=nu[i,3]*nu[j,3]/np.sum((nu[i,:3]-nu[j,:3])**2)**0.5
    for k in range(basis_num):
        for l in range(basis_num):
            e_e+=p_mat[k,l]*F_mat[k,l]
    return(v_nn+e_e/2)


def contract_int(ba,contract_cum,int_mat,contract_mat,gabcd=0):
    contract_cum=np.insert(contract_cum,0,0)
    # print(contract_cum)
    if gabcd != 1:
        for i in range(len(contract_cum)-1):
            for j in range(len(contract_cum)-1):
                if i>=j:
                    temp=0.0

                    for ii in range(contract_cum[i],contract_cum[i+1]):
                        for jj in range(contract_cum[j],contract_cum[j+1]):
                            temp+= int_mat[ii,jj]*ba[ii][7]*ba[jj][7]

                    contract_mat[i,j]=temp
                    contract_mat[j,i]=temp

    else:
        for i in range(len(contract_cum)-1):
            for j in range(len(contract_cum)-1):
                for k in range(len(contract_cum)-1):
                    for l in range(len(contract_cum)-1):
                        if i>=j and k>=l and i>=k :
                            temp=0.0
                            for ii in range(contract_cum[i],contract_cum[i+1]):
                                for jj in range(contract_cum[j],contract_cum[j+1]):                        
                                    for kk in range(contract_cum[k],contract_cum[k+1]):
                                        for ll in range(contract_cum[l],contract_cum[l+1]):
                                            temp+=int_mat[ii,jj,kk,ll]*ba[ii][7]*ba[jj][7]*ba[kk][7]*ba[ll][7]
                            contract_mat[i,j,k,l]=temp
                            contract_mat[i,j,l,k]=temp
                            contract_mat[j,i,k,l]=temp
                            contract_mat[j,i,l,k]=temp
                            contract_mat[k,l,i,j]=temp
                            contract_mat[k,l,j,i]=temp
                            contract_mat[l,k,i,j]=temp
                            contract_mat[l,k,j,i]=temp
    

    return(contract_mat)



# 
#   以下是调用编译的库文件中的函数计算需要的各种矩阵元
# 

from ctypes import cdll,c_double,POINTER,c_int
if os.path.exists("./analy_int.so"):
    cfun=cdll.LoadLibrary("./analy_int.so")


def smat_pri(ba,smat,smat_dia):
    smat1=smat.ctypes.data_as(POINTER(c_double))
    ba1=ba.ctypes.data_as(POINTER(c_double))
    smat_dia1=smat_dia.ctypes.data_as(POINTER(c_double))
    len1=c_int(len(ba))

    cfun.smat_pri(ba1,smat1,smat_dia1,len1)

def smat_con(ba,smat,smat_dia,contract_cum):
    smat1=smat.ctypes.data_as(POINTER(c_double))
    ba1=ba.ctypes.data_as(POINTER(c_double))
    smat_dia1=smat_dia.ctypes.data_as(POINTER(c_double))
    contract_cum1=contract_cum.ctypes.data_as(POINTER(c_int))
    len1=c_int(len(contract_cum)-1)
    len2=c_int(len(ba))

    cfun.smat_con(ba1,smat1,smat_dia1,contract_cum1,len1,len2)


def tmat_pri(ba,tmat,smat_dia):
    tmat1=tmat.ctypes.data_as(POINTER(c_double))
    ba1=ba.ctypes.data_as(POINTER(c_double))
    smat_dia1=smat_dia.ctypes.data_as(POINTER(c_double))
    len1=c_int(len(ba))

    cfun.tmat_pri(ba1,tmat1,smat_dia1,len1)

def tmat_con(ba,tmat,smat_dia,contract_cum):
    tmat1=tmat.ctypes.data_as(POINTER(c_double))
    ba1=ba.ctypes.data_as(POINTER(c_double))
    smat_dia1=smat_dia.ctypes.data_as(POINTER(c_double))
    contract_cum1=contract_cum.ctypes.data_as(POINTER(c_int))
    len1=c_int(len(contract_cum)-1)

    cfun.tmat_con(ba1,tmat1,smat_dia1,contract_cum1,len1)


def vmat_pri(ba,vmat,smat_dia,nu):
    vmat1=vmat.ctypes.data_as(POINTER(c_double))
    ba1=ba.ctypes.data_as(POINTER(c_double))
    smat_dia1=smat_dia.ctypes.data_as(POINTER(c_double))
    nu1=nu.ctypes.data_as(POINTER(c_double))
    len1=c_int(len(ba))
    len2=c_int(len(nu))

    cfun.vmat_pri(ba1,vmat1,smat_dia1,nu1,len1,len2)

def vmat_con(ba,vmat,smat_dia,nu,contract_cum):
    vmat1=vmat.ctypes.data_as(POINTER(c_double))
    ba1=ba.ctypes.data_as(POINTER(c_double))
    smat_dia1=smat_dia.ctypes.data_as(POINTER(c_double))
    contract_cum1=contract_cum.ctypes.data_as(POINTER(c_int))
    nu1=nu.ctypes.data_as(POINTER(c_double))
    len1=c_int(len(contract_cum)-1)
    len2=c_int(len(nu))

    cfun.vmat_con(ba1,vmat1,smat_dia1,contract_cum1,nu1,len1,len2)


def gmat_pri(ba,gmat,smat_dia):
    gmat1=gmat.ctypes.data_as(POINTER(c_double))
    ba1=ba.ctypes.data_as(POINTER(c_double))
    smat_dia1=smat_dia.ctypes.data_as(POINTER(c_double))
    len1=c_int(len(ba))

    cfun.gab_pri(ba1,gmat1,smat_dia1,len1)


def gmat_con(ba,gmat,smat_dia,contract_cum):
    gmat1=gmat.ctypes.data_as(POINTER(c_double))
    ba1=ba.ctypes.data_as(POINTER(c_double))
    smat_dia1=smat_dia.ctypes.data_as(POINTER(c_double))
    contract_cum1=contract_cum.ctypes.data_as(POINTER(c_int))
    len1=c_int(len(contract_cum)-1)


    cfun.gab_con(ba1,gmat1,smat_dia1,contract_cum1,len1)


