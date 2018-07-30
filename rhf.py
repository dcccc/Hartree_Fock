#coding=utf-8
import math
import numpy as np


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
	elif i1>i2:
		s=xa*s_ab(x1,i1-1,a1,x2,i2,a2)+1.0/2.0/p*(i1*s_ab(x1,i1-2,a1,x2,i2,a2)+
			i2*s_ab(x1,i1-1,a1,x2,i2-1,a2))
	else:
		s=xa*s_ab(x1,i1,a1,x2,i2-1,a2)+1.0/2.0/p*(i1*s_ab(x1,i1-1,a1,x2,i2-1,a2)+
			i2*s_ab(x1,i1,a1,x2,i2-2,a2))

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
		t=(a1-2*a1**2*(xa**2+1.0/2.0/p))*s_ab(x1,0,a1,x2,0,a2)
	elif i1>=i2:
		t=xa*t_ab(x1,i1-1,a1,x2,i2,a2)+\
		0.5/p*((i1-1)*t_ab(x1,i1-2,a1,x2,i2,a2)+i2*t_ab(x1,i1-1,a1,x2,i2-1,a2))+\
		a2/p*(2.0*a1*s_ab(x1,i1,a1,x2,i2,a2)-i1*s_ab(x1,i1-2,a1,x2,i2,a2))
	else:
		t=xb*t_ab(x1,i1,a1,x2,i2-1,a2)+\
		0.5/p*(i1*t_ab(x1,i1-1,a1,x2,i2-1,a2)+i2*t_ab(x1,i1,a1,x2,i2-2,a2))+\
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

	if t<0:
		e=0.0
	elif t==0:
		e=(p/math.pi)**0.5*s_ab(x1,i1,a1,x2,i2,a2)
	else:
		e=(i1*hermite(x1-1,i1,a1,x2,i2,a2,t-1)+i2*hermite(x1,i1,a1,x2,i2-1,a2,t-1))/2.0/p/t
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
	elif x<=19.35 and x > 0.18:
		boy=0.0
		for i in range(12):
			factorial=1
			for j in range(i):
				factorial*=(2*n+2*j+1)
			boy+=math.exp(-x)*(2*x)**i/factorial

	elif x>19.35:
		boy=math.factorial(math.factorial(2*n))/2**(n+1)*(math.pi/2**(2*n+1))**0.5
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
	elif u==0 or v==0 or t==0 :
		r=(-2*p)**n*boys(n,p*rpc2)
	elif u>0:
		r=(u-1)*ruvt(u-2,v,t,ga,gb,n+1,nu)+xa*ruvt(u-1,v,t,ga,gb,n+1,nu)
	elif v>0:
		r=(v-1)*ruvt(u,v-2,t,ga,gb,n+1,nu)+ya*ruvt(u,v-1,t,ga,gb,n+1,nu)
	elif t>0:
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

	elif u>0:
		r=(u-1)*ruvt(u-2,v,t,ga,gb,gc,gd,n+1,nu)+px*ruvt(u-1,v,t,ga,gb,gc,gd,n+1,nu)
	elif v>0:
		r=(v-1)*ruvt(u,v-2,t,ga,gb,gc,gd,n+1,nu)+py*ruvt(u,v-1,t,ga,gb,gc,gd,n+1,nu)
	elif t>0:
		r=(t-1)*ruvt(u,v,t-2,ga,gb,gc,gd,n+1,nu)+pz*ruvt(u,v,t-1,ga,gb,gc,gd,n+1,nu)

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


def g_ab(p_mat,gabcd_mat,i,j):
	gv=0.0
	for k in range(6):
		for l in range(6):
			gv+=p_mat[k,l]*(gabcd_mat[gabdc_index(i,j,k,l)]-gabcd_mat[gabdc_index(i,l,k,j)]/2)
	return(gv)



def c_normarize(c_mat,e_mat,s_mat):

	ce_mat=np.vstack((e_mat,c_mat)).T.tolist()
	ce_mat=np.array(sorted(ce_mat,key=lambda x: x[0]))
	e_mat,c_mat=ce_mat.T[0,:],ce_mat.T[1:,:]
	for i in range(6):
		N=0.0
		for j in range(6):
			for k in range(6):
				N+=c_mat[j,i]*c_mat[k,i]*s_mat[j,k]
		c_mat[:,i]=c_mat[:,i]/N**0.5
	return(e_mat,c_mat)

def all_energy(p_mat,F_mat,nu):
	v_nn=0.0
	e_e=0.0
	for i in range(len(nu)):
		for j in range(len(nu)):
			if i>j:
				v_nn+=nu[i,3]*nu[j,3]/np.sum((nu[i,:3]-nu[j,:3])**2)**0.5
	for k in range(6):
		for l in range(6):
			e_e+=p_mat[k,l]*F_mat[k,l]
	return(v_nn+e_e/2)

#########  sto-3g
######### a [ 0.1688554040, 0.6239137298, 3.4252509140 ]

# j=0.370065/0.5291772

# jj=-0.370065/0.5291772


ba=[[0,0,j,0,0,0,i]  for j in [0.370065/0.5291772,-0.370065/0.5291772] for i in [3.4252509140, 0.6239137298,0.1688554040 ]]

# sss=np.zeros((6,6))

# for i in range(6):
# 	for j in range(6):
# 		for k in range(6):
# 			for l in range(6):
# 				if i>=j and k>=l and i>=k  :
# 					ij=i*(i+1)/2+j
# 					sss=g(ba[i],ba[j],ba[k],ba[l])
# 					sss=sss/(sab(ba[i],ba[i])*sab(ba[j],ba[j])*sab(ba[k],ba[k])*sab(ba[l],ba[l]))**0.5
# 					print("%d-%d-%d-%d     %6.5f"  %(i+1,j+1,k+1,l+1,sss))

# i=6-1
# j=6-1
# k=6-1
# l=1-1

# print(ba[i])
# print(ba[l])

# sss=g(ba[i],ba[j],ba[k],ba[l])
# sss=sss/(sab(ba[i],ba[i])*sab(ba[j],ba[j])*sab(ba[k],ba[k])*sab(ba[l],ba[l]))**0.5

# print(sss)


# ijkl=[[6,1,4,2],[6,1,5,2],[6,1,5,4],[6,2,4,1],[6,2,4,3],[6,2,5,1],[6,2,5,4],[6,3,5,1],[6,3,5,4],\
# [6,4,2,1],[6,4,4,3],[6,4,4,4],[6,4,5,1],[6,4,5,2],[6,4,5,3],[6,4,5,4],[6,4,5,5],[6,4,6,1],\
# [6,4,6,2],[6,4,6,3],[6,4,6,4],[6,5,1,1],[6,5,2,1],[6,5,2,2],[6,5,3,1],[6,5,4,1],[6,5,4,2],\
# [6,5,4,3],[6,5,4,4],[6,5,5,1],[6,5,5,2],[6,5,5,3],[6,5,5,4],[6,5,5,5],[6,5,6,1],[6,5,6,2],\
# [6,5,6,3],[6,5,6,4],[6,5,6,5],[6,6,1,1],[6,6,2,1],[6,6,2,2],[6,6,3,1],[6,6,3,2],[6,6,4,1],\
# [6,6,4,2],[6,6,4,3],[6,6,4,4],[6,6,5,1],[6,6,5,2],[6,6,5,3],[6,6,5,4],[6,6,5,5],[6,6,6,1],\
# [6,6,6,2],[6,6,6,3],[6,6,6,4],[6,6,6,5],[6,6,6,6]]

# for ii in ijkl:
# 	i,j,k,l=tuple(ii)
# 	sss=g(ba[i-1],ba[j-1],ba[k-1],ba[l-1])
# 	sss=sss/(sab(ba[i-1],ba[i-1])*sab(ba[j-1],ba[j-1])*sab(ba[k-1],ba[k-1])*sab(ba[l-1],ba[l-1]))**0.5
# 	print("%d-%d-%d-%d     %10.9f"  %(i,j,k,l,sss))



ba=[[0,0,j,0,0,0,i]  for j in [0.370065/0.5291772,-0.370065/0.5291772] for i in [3.4252509140, 0.6239137298,0.1688554040 ]]

nu=[[0,0,j,1]  for j in [0.370065/0.5291772,-0.370065/0.5291772]]
delta_e=1.0
delta_p=1.0
n_electron=2
n_c=0


s_mat=np.zeros((6,6))
for i in range(6):
	for j in range(6):
		s_mat[i,j]=sab(ba[i],ba[j])/(sab(ba[i],ba[i])*sab(ba[j],ba[j]))**0.5


t_mat=np.zeros((6,6))
for i in range(6):
	for j in range(6):
		t_mat[i,j]=tab(ba[i],ba[j])/(sab(ba[i],ba[i])*sab(ba[j],ba[j]))**0.5

v_mat=np.zeros((6,6))
for i in range(6):
	for j in range(6):
		temp=0.0
		temp=vab(ba[i],ba[j],nu[0])+vab(ba[i],ba[j],nu[1])
		v_mat[i,j]=temp/(sab(ba[i],ba[i])*sab(ba[j],ba[j]))**0.5

gabcd_mat={}

for i in range(6):
	for j in range(6):
		for k in range(6):
			for l in range(6):
				if i>=j and k>=l and i>=k  :
					
					temp=gabcd(ba[i],ba[j],ba[k],ba[l])
					temp=temp/(sab(ba[i],ba[i])*sab(ba[j],ba[j])*sab(ba[k],ba[k])*sab(ba[l],ba[l]))**0.5
					gabcd_index="%d-%d-%d-%d"  %(i,j,k,l)
					gabcd_mat[gabcd_index]=temp




d_mat,l_mat=np.linalg.eig(s_mat)
d_mat=np.diag(tuple(d_mat**-0.5))



s_mat_05=np.dot(np.dot(l_mat,d_mat),l_mat.T)

F_mat=t_mat+v_mat

F_mat_1=np.dot(np.dot(s_mat_05,F_mat),s_mat_05.T)
e_mat,C_mat=np.linalg.eig(F_mat_1)
e_mat,C_mat=c_normarize(np.dot(s_mat_05,C_mat),e_mat,s_mat)




p_mat=pab(C_mat,n_electron)


G_mat=np.zeros((6,6))
for i in range(6):
	for j in range(6):
		G_mat[i,j]=g_ab(p_mat,gabcd_mat,i,j)

F_mat=t_mat+v_mat+G_mat
F_mat_1=np.dot(np.dot(s_mat_05,F_mat),s_mat_05.T)
e_mat,C_mat=np.linalg.eig(F_mat_1)
e_mat,C_mat=c_normarize(np.dot(s_mat_05,C_mat),e_mat,s_mat)
p_mat=pab(C_mat,n_electron)
ene= all_energy(p_mat,F_mat+t_mat+v_mat+G_mat,np.array(nu))



# def scf(t_mat,v_mat,s_mat,gabcd_mat,C_mat,s_mat_05,n_electron)

# 	p_mat=pab(C_mat,n_electron)
# 	G_mat=np.zeros((6,6))
# 	for i in range(6):
# 		for j in range(6):
# 			G_mat[i,j]=g_ab(p_mat,gabcd_mat,n_electron,i,j)

# 	F_mat=t_mat+v_mat+G_mat
# 	F_mat_1=np.dot(np.dot(s_mat_05,F_mat),s_mat_05.T)
# 	e_mat,C_mat=np.linalg.eig(F_mat_1)
# 	ene= all_energy(p_mat,F_mat,np.array(nu))


while(delta_e>1*10**-8 or delta_p>1*10**-6):


	
	G_mat=np.zeros((6,6))
	for i in range(6):
		for j in range(6):
			G_mat[i,j]=g_ab(p_mat,gabcd_mat,i,j)

	F_mat=t_mat+v_mat+G_mat
	F_mat_1=np.dot(np.dot(s_mat_05,F_mat),s_mat_05.T)
	e_mat,C_mat=np.linalg.eig(F_mat_1)
	e_mat,C_mat=c_normarize(np.dot(s_mat_05,C_mat),e_mat,s_mat)
	p_mat1=pab(C_mat,n_electron)
	ene1= all_energy(p_mat,F_mat+t_mat+v_mat,np.array(nu))

	delta_e=abs(ene1-ene)
	delta_p=0.0
	for i in range(6):
		for j in range(6):
			delta_p+=abs(p_mat1[i,j]-p_mat[i,j])/36
	ene=ene1
	p_mat=p_mat1[:,:]
	n_c=n_c+1
	# print(n_c)

	print("%-4d     %10.8f     %10.8f     %10.8f      %10.8f"  %(n_c, ene, e_mat[0],delta_e, delta_p))

