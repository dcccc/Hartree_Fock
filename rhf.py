#coding=utf-8
import math
import numpy as np
import scipy.special as spec

def s_ab(x1,i1,a1,x2,i2,a2):
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
		s=(math.pi/p)**0.5*kab
	elif i1>0 and i2==0:
		s=xa*s_ab(x1,i1-1,a1,x2,i2,a2)+0.5/p*((i1-1)*s_ab(x1,i1-2,a1,x2,i2,a2))
	else:
		s=xb*s_ab(x1,i1,a1,x2,i2-1,a2)+0.5/p*(i1*s_ab(x1,i1-1,a1,x2,i2-1,a2)+
			(i2-1)*s_ab(x1,i1,a1,x2,i2-2,a2))

	return(s)

def  sab(ga,gb):
	x1,y1,z1,i1,j1,k1,a1=tuple(ga)
	x2,y2,z2,i2,j2,k2,a2=tuple(gb)
	s1=s_ab(x1,i1,a1,x2,i2,a2)
	s2=s_ab(y1,j1,a1,y2,j2,a2)
	s3=s_ab(z1,k1,a1,z2,k2,a2)
	# print(s1,s2,s3)
	return(s1*s2*s3)

def t_ab(x1,i1,a1,x2,i2,a2):
	p=a1+a2
	niu=a1*a2/p
	px=(a1*x1+a2*x2)/p
	xa=px-x1
	xb=px-x2

	if i1<0 or i2<0:
		t=0.0
	elif i1==0 and i2==0:
		t=(a1-2*a1**2*(xa**2+0.5/p))*s_ab(x1,0,a1,x2,0,a2)
	elif i1>0 and i2==0:
		t=xa*t_ab(x1,i1-1,a1,x2,0,a2)+\
		0.5/p*((i1-1)*t_ab(x1,i1-2,a1,x2,0,a2))+\
		a2/p*(2.0*a1*s_ab(x1,i1,a1,x2,0,a2)-(i1-1)*s_ab(x1,i1-2,a1,x2,0,a2))
	else:
		t=xb*t_ab(x1,i1,a1,x2,i2-1,a2)+\
		0.5/p*(i1*t_ab(x1,i1-1,a1,x2,i2-1,a2)+(i2-1)*t_ab(x1,i1,a1,x2,i2-2,a2))+\
		a1/p*(2.0*a2*s_ab(x1,i1,a1,x2,i2,a2)-(i2-1)*s_ab(x1,i1,a1,x2,i2-2,a2))

	return(t)


def tab(ga,gb):
	x1,y1,z1,i1,j1,k1,a1=tuple(ga)
	x2,y2,z2,i2,j2,k2,a2=tuple(gb)
	t1=t_ab(x1,i1,a1,x2,i2,a2)*s_ab(y1,j1,a1,y2,j2,a2)*s_ab(z1,k1,a1,z2,k2,a2)
	t2=s_ab(x1,i1,a1,x2,i2,a2)*t_ab(y1,j1,a1,y2,j2,a2)*s_ab(z1,k1,a1,z2,k2,a2)
	t3=s_ab(x1,i1,a1,x2,i2,a2)*s_ab(y1,j1,a1,y2,j2,a2)*t_ab(z1,k1,a1,z2,k2,a2)
	return(t1+t2+t3)


def hermite(x1,i1,a1,x2,i2,a2,t):
	p=a1+a2

	if t<0 or t > i1+i2 or i1<0 or i2<0:
		e=0.0
	elif t==0:
		e=(p/math.pi)**0.5*s_ab(x1,i1,a1,x2,i2,a2)
	else:
		e=(i1*hermite(x1,i1-1,a1,x2,i2,a2,t-1)+i2*hermite(x1,i1,a1,x2,i2-1,a2,t-1))/2.0/p/t
	return(e)


def boys(n,x):
	# if n<0:
	# 	boy=0.0
	# print(n,x)
	if x==0.0:
		boy=1.0/(2*n+1)
	elif n==0 and x !=0.0:
		boy=(math.pi)**0.5*math.erf(x**0.5)/(x**0.5)/2
	elif x<0.18:
		boy=0.0
		for k in range(7):
			boy+=(-x)**k/(math.factorial(k)*(2.0*k+2.0*n+1))
	# elif x<=19.35 and x > 0.18:
	# 	boy=0.0
	# 	for i in range(12):
	# 		factorial=1
	# 		for j in range(i):
	# 			factorial*=(2*n+2*j+1)
	# 			boy+=math.exp(-x)*(2*x)**i/factorial

	else:
		boy=0.5/x**(n+0.5)*spec.gamma(n+0.5)*spec.gammainc(n+0.5,x)
	# elif x>=35 and n == 0:
	# 	boy=0.5*(math.pi/x)**0.2
	# elif x
	# # else:
	# # 	boy=0.5/x*((2*n-1)*boys(n-1,x)-math.exp(-x))/2.0/x
	# print(boy)
	return(boy)


def ruvt(u,v,t,ga,gb,n,nu):

	x1,y1,z1,i1,j1,k1,a1=tuple(ga)
	x2,y2,z2,i2,j2,k2,a2=tuple(gb)
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
		r=(-2*p)**n*boys(n,p*rpc2)
		# print(r,n,p,rpc2)
	elif u>0 and v==0 and t==0:
		r=(u-1)*ruvt(u-2,0,0,ga,gb,n+1,nu)+xa*ruvt(u-1,0,0,ga,gb,n+1,nu)
	elif v>0 and t==0:
		r=(v-1)*ruvt(u,v-2,0,ga,gb,n+1,nu)+ya*ruvt(u,v-1,0,ga,gb,n+1,nu)
	else:
		r=(t-1)*ruvt(u,v,t-2,ga,gb,n+1,nu)+za*ruvt(u,v,t-1,ga,gb,n+1,nu)

	return(r)



def vab(ga,gb,nu):
	x1,y1,z1,i1,j1,k1,a1=tuple(ga)
	x2,y2,z2,i2,j2,k2,a2=tuple(gb)
	p=a1+a2
	vv=0.0
	for u in range(i1+i2+1):
		for v in range(j1+j2+1):
			for t in range(k1+k2+1):
				vv+=hermite(x1,i1,a1,x2,i2,a2,u)*hermite(y1,j1,a1,y2,j2,a2,v) \
				*hermite(z1,k1,a1,z2,k2,a2,t)*ruvt(u,v,t,ga,gb,0,nu)

	return(math.pi*vv*2*nu[3]/p)


def ruvt_2(u,v,t,ga,gb,gc,gd,n):

	x1,y1,z1,i1,j1,k1,a1=tuple(ga)
	x2,y2,z2,i2,j2,k2,a2=tuple(gb)
	x3,y3,z3,i3,j3,k3,a3=tuple(gc)
	x4,y4,z4,i4,j4,k4,a4=tuple(gd)
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
	x1,y1,z1,i1,j1,k1,a1=tuple(ga)
	x2,y2,z2,i2,j2,k2,a2=tuple(gb)
	x3,y3,z3,i3,j3,k3,a3=tuple(gc)
	x4,y4,z4,i4,j4,k4,a4=tuple(gd)
	p1=a1+a2
	p2=a3+a4
	pp=p1*p2/(p1+p2)
	gv=0.0
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


def g_ab(p_mat,gabcd_mat,i,j,basis_num):
	gv=0.0
	for k in range(basis_num):
		for l in range(basis_num):
			gv+=p_mat[k,l]*(gabcd_mat[gabdc_index(i,j,k,l)]-gabcd_mat[gabdc_index(i,l,k,j)]/2)
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
	for i in range(len(nu)):
		for j in range(len(nu)):
			if i>j:
				v_nn+=nu[i,3]*nu[j,3]/np.sum((nu[i,:3]-nu[j,:3])**2)**0.5
	for k in range(basis_num):
		for l in range(basis_num):
			e_e+=p_mat[k,l]*F_mat[k,l]
	return(v_nn+e_e/2)



#输入参数


delta_e=1.0
delta_p=1.0
n_electron=10
n_c=0

coor=[["C",  0.         ,  0.          ,   0.           ],\
      ["H" , 0.         ,  0.          ,   2.0220696   ],\
      ["H" , 0.         , -1.90636644 ,  -0.67400233  ],\
      ["H" ,-1.65096177,  0.95318323 ,  -0.67400233  ],\
      ["H" , 1.65096177,  0.95318323 ,  -0.67400233  ]]

from pp import basis_data_p

nu=[]
ba=[]
basis_data=basis_data_p("basis_file.txt")

for i in coor :
	if i[0]=="C":
		temp1=[i[1],i[2],i[3],6.0]
		nu.append(temp1)
		basis_data_t=basis_data["C"]["S"]
		for j in range(len(basis_data_t)):
			temp=[i[1],i[2],i[3],0,0,0,float(basis_data_t[j][0])]
			ba.append(temp)
		basis_data_t=basis_data["C"]["SP"]
		for j in range(len(basis_data_t)):
			temp=[i[1],i[2],i[3],0,0,0,float(basis_data_t[j][0])]
			ba.append(temp)
			for k in [[1,0,0],[0,1,0],[0,0,1]]:
				temp=[i[1],i[2],i[3],k[0],k[1],k[2],float(basis_data_t[j][0])]
				ba.append(temp)
	if i[0]=="H":
		temp1=[i[1],i[2],i[3],1.0]
		nu.append(temp1)
		basis_data_t=basis_data["H"]["S"]
		for j in range(len(basis_data_t)):
			temp=[i[1],i[2],i[3],0,0,0,float(basis_data_t[j][0])]
			ba.append(temp)

basis_num=len(ba)

print(basis_data)
# print(ruvt(0,2,0,ba[7],ba[7],0,nu[0]))
# print(vab(ba[7],ba[7],nu[2]))
# print(ba[7],ba[7],nu[0])

# print(ruvt(0,2,0,ba[7],ba[7],0,nu[2]))


# exit()

#重叠矩阵
s_mat=np.zeros((basis_num,basis_num))
for i in range(basis_num):
	for j in range(basis_num):
		s_mat[i,j]=sab(ba[i],ba[j])/(sab(ba[i],ba[i])*sab(ba[j],ba[j]))**0.5

#动能矩阵
t_mat=np.zeros((basis_num,basis_num))
for i in range(basis_num):
	for j in range(basis_num):
		t_mat[i,j]=tab(ba[i],ba[j])/(sab(ba[i],ba[i])*sab(ba[j],ba[j]))**0.5

#势能矩阵
v_mat=np.zeros((basis_num,basis_num))
for i in range(basis_num):
	for j in range(basis_num):
		temp=0.0
		for k in nu:
			temp+=vab(ba[i],ba[j],k)
		v_mat[i,j]=temp/(sab(ba[i],ba[i])*sab(ba[j],ba[j]))**0.5





#双电子积分
gabcd_mat={}

for i in range(basis_num):
	for j in range(basis_num):
		for k in range(basis_num):
			for l in range(basis_num):
				if i>=j and k>=l and i>=k  :
					
					temp=gabcd(ba[i],ba[j],ba[k],ba[l])
					temp=temp/(sab(ba[i],ba[i])*sab(ba[j],ba[j])*sab(ba[k],ba[k])*sab(ba[l],ba[l]))**0.5
					gabcd_index="%d-%d-%d-%d"  %(i,j,k,l)
					gabcd_mat[gabcd_index]=temp



np.savetxt("s_mat.csv",s_mat,delimiter=',')
np.savetxt("t_mat.csv",t_mat,delimiter=',')
np.savetxt("v_mat.csv",v_mat,delimiter=',')


#S^-0.5的获取
d_mat,l_mat=np.linalg.eig(s_mat)
d_mat=np.diag(tuple(d_mat**-0.5))
s_mat_05=np.dot(np.dot(l_mat,d_mat),l_mat.T)



H_mat=t_mat-v_mat
F_mat_1=np.dot(np.dot(s_mat_05,H_mat),s_mat_05.T)


#FOCK的对角化
e_mat,C_mat=np.linalg.eig(F_mat_1)

#参数矩阵的排序归一化
e_mat,C_mat=c_normarize(np.dot(s_mat_05,C_mat),e_mat,s_mat,basis_num)

#密度矩阵
p_mat=pab(C_mat,n_electron)


ene=0.0


#SCF迭代
while(delta_e>1*10**-8 or delta_p>1*10**-6):


	#G矩阵
	G_mat=np.zeros((basis_num,basis_num))
	for i in range(basis_num):
		for j in range(basis_num):
			G_mat[i,j]=g_ab(p_mat,gabcd_mat,i,j,basis_num)

	F_mat=H_mat+G_mat
	F_mat_1=np.dot(np.dot(s_mat_05,F_mat),s_mat_05.T)
	e_mat,C_mat=np.linalg.eig(F_mat_1)
	e_mat,C_mat=c_normarize(np.dot(s_mat_05,C_mat),e_mat,s_mat,basis_num)
	p_mat1=pab(C_mat,n_electron)

	#获得体系能量
	ene1= all_energy(p_mat,F_mat+H_mat,np.array(nu),basis_num)

	#能量差
	delta_e=abs(ene1-ene)

	#密度矩阵差
	delta_p=0.0
	for i in range(basis_num):
		for j in range(basis_num):
			delta_p+=abs(p_mat1[i,j]-p_mat[i,j])/basis_num**2

	ene=ene1
	p_mat=p_mat1[:,:]
	n_c=n_c+1


	#输出能量，能量的变化，以及密度矩阵的变化
	print("%-4d     %10.8f     %10.8f      %10.8f"  %(n_c, ene, delta_e, delta_p))

