#coding=utf-8
import math,time,sys,os
import numpy as np
from grid_int import  xyzw_list_gen,xc_int_c, density_pot_list_c,density_pot_list_c_contra
from pre_pro import basis_list, xyz_file_pp,basis_data_p,input_pp
from analy_int import c_normarize,all_energy,pab,g_ab_j,g_ab_jk,contract_int
from analy_int import smat_pri,tmat_pri,vmat_pri,gmat_pri,smat_con,tmat_con,vmat_con,gmat_con
from ctypes import cdll,c_double,POINTER,c_int
cfun=cdll.LoadLibrary("./grid_int.so")



def scf(ba,nu,contract_list,dft=0,mp2=0,contracted=1):

	basis_num=len(ba)

	atom_num=len(nu)


	delta_e=1.0
	delta_p=1.0
	n_c=0
	
	time0=time.time()

	print("\n单双电子积分计算中....")
	if contracted==1:    #使用收缩基计算时

		contract_cum=np.cumsum(np.array(np.array([y for x in contract_list for y in x])))
		contract_cum1=np.insert(contract_cum,0,0)
		basis_num=len(contract_cum)

		#重叠矩阵
		s_mat_dia=np.zeros((len(ba)))
		s_mat=np.zeros((basis_num,basis_num))
		smat_con(ba,s_mat,s_mat_dia,contract_cum1)

		#动能矩阵
		t_mat=np.zeros((basis_num,basis_num))
		tmat_con(ba,t_mat,s_mat_dia,contract_cum1)
		
		#势能矩阵
		v_mat=np.zeros((basis_num,basis_num))
		vmat_con(ba,v_mat,s_mat_dia,nu,contract_cum1)
		
		#双电子积分
		gabcd_mat=np.zeros((basis_num,basis_num,basis_num,basis_num))
		gmat_con(ba,gabcd_mat,s_mat_dia,contract_cum1)


	else:      #使用原始基计算时
		#重叠矩阵
		s_mat_dia=np.zeros((len(ba)))
		s_mat=np.zeros((basis_num,basis_num))
		smat_pri(ba,s_mat,s_mat_dia)
		# print(s_mat)
		# exit()
		#动能矩阵
		t_mat=np.zeros((basis_num,basis_num))
		tmat_pri(ba,t_mat,s_mat_dia)
		
		#势能矩阵
		v_mat=np.zeros((basis_num,basis_num))
		vmat_pri(ba,v_mat,s_mat_dia,nu)
		
		#双电子积分
		gabcd_mat=np.zeros((basis_num,basis_num,basis_num,basis_num))
		gmat_pri(ba,gabcd_mat,s_mat_dia)
		
		
	#S^-0.5的获取
	d_mat,l_mat=np.linalg.eig(s_mat)
	d_mat=np.diag(tuple(d_mat**-0.5))
	s_mat_05=np.dot(np.dot(l_mat,d_mat),l_mat.T)
	
	
	
	H_mat=t_mat-v_mat
	F_mat_1=np.dot(np.dot(s_mat_05,H_mat),s_mat_05.T)
	
	
	#FOCK的对角化
	e_mat,C_mat=np.linalg.eig(F_mat_1)
	
	#系数矩阵的排序归一化
	e_mat,C_mat=c_normarize(np.dot(s_mat_05,C_mat),e_mat,s_mat,basis_num)
	
	#密度矩阵
	p_mat=pab(C_mat,n_electron)
	

	ene=0.0
	
	print('积分计算消耗时间： %fs'   %(time.time()-time0))
	print("\n进入scf计算")
	time1=time.time()
	print("%-s            %-s            %-s      %-s"  %("cycle", "E(hartree)","delta_e(hartree)","delta_p"))
	#SCF迭代
	while(delta_e>1*10**-8 or delta_p>1*10**-6):
	
	
		#G矩阵
		G_mat=np.zeros((basis_num,basis_num))		
		# print(p_mat)

		#xc矩阵
		if dft==1:

			xyzw_list=[]

			# 生成积分网格点(DFT)
			xyzw_list=np.array(xyzw_list_gen(atom_xyz))
			contract_cum=np.cumsum(np.array(contract_list).flatten())
			basis_num_con=len(contract_list)

			basis_num1=len(ba)
			ks_mat=np.zeros((basis_num1,basis_num1))

			for i in range(basis_num):
				for j in range(basis_num):
					G_mat[i,j]=g_ab_j(p_mat,gabcd_mat,i,j,basis_num)			
			
			 #使用收缩基计算时
			if contracted==1:
				p_mat_uncon=np.zeros((basis_num1,basis_num1))
				contract_cum1=np.insert(contract_cum,0,0)
				for i in range(basis_num):
					for j in range(basis_num):
						if i>=j:
							p_mat_uncon[contract_cum1[i]:contract_cum1[i+1],contract_cum1[j]:contract_cum1[j+1],]=p_mat[i,j]
							p_mat_uncon[contract_cum1[j]:contract_cum1[j+1],contract_cum1[i]:contract_cum1[i+1],]=p_mat[j,i]


				# 收缩基格点电子密度
				density_list=np.array(density_pot_list_c_contra(ba,xyzw_list,s_mat_dia,p_mat_uncon))
			

			#  使用原始基计算时
			else:
				# 原始基格点电子密度
				density_list=np.array(density_pot_list_c(ba,xyzw_list,s_mat_dia,p_mat))
				p_mat_uncon=p_mat

			###  收缩系数对总电子密度的影响

			for i in range(basis_num1):
				for j in range(basis_num1):
					if i>=j :
						
						#  使用原始基计算得到的xc矩阵
						ks_mat[i,j]=xc_int_c(ba,i,j,s_mat_dia,p_mat_uncon,xyzw_list,density_list)*\
						-(9.0/8)*(3.0/math.pi)**(1.0/3)*0.7				
	
						ks_mat[j,i]=ks_mat[i,j]

		else:
			#  计算G矩阵
			for i in range(basis_num):
				for j in range(basis_num):
					G_mat[i,j]=g_ab_jk(p_mat,gabcd_mat,i,j,basis_num)

		if contracted==1 and dft==1:
			#  使用原始基计算得到的xc矩阵 计算收缩的xc矩阵
			ks_mat=contract_int(ba,contract_cum,ks_mat,np.zeros((basis_num,basis_num)))

		elif dft==0 :
			ks_mat=np.zeros((basis_num,basis_num))


		#Fock矩阵 及其对角化等
		F_mat=H_mat+G_mat+ks_mat
		F_mat_1=np.dot(np.dot(s_mat_05,F_mat),s_mat_05.T)
		e_mat,C_mat=np.linalg.eig(F_mat_1)
		e_mat,C_mat=c_normarize(np.dot(s_mat_05,C_mat),e_mat,s_mat,basis_num)
		p_mat1=pab(C_mat,n_electron)
	
		#获得体系能量
		ene1= all_energy(p_mat,F_mat+H_mat+ks_mat,np.array(nu),basis_num)
	
		#能量差
		delta_e=abs(ene1-ene)
	
		#密度矩阵差
		delta_p=0.0
		for i in range(basis_num):
			for j in range(basis_num):
				delta_p+=abs(p_mat1[i,j]-p_mat[i,j])/basis_num**2
	
		ene=ene1

		# 密度矩阵混合
		p_mat=p_mat1
		n_c=n_c+1
	
	
		#输出能量，能量的变化，以及密度矩阵的变化
		print("%-4d            %10.8f           %10.8f            %10.8f"  %(n_c, ene, delta_e, delta_p))


	print("scf计算总时间为%fs"  %(time.time()-time1))


	time2=time.time()
	if mp2==1 :#and dft !=1 :
		print("\nmp2校正能计算中")
		go_mat=np.zeros((basis_num,basis_num,basis_num,basis_num))
		go_mat = np.asarray(go_mat, order='C')
		C_mat = np.asarray(C_mat, order='C')



		go_mat = np.asarray(go_mat, order='C')
		go_mat1=go_mat.ctypes.data_as(POINTER(c_double))
		gabcd_mat = np.asarray(gabcd_mat, order='C')
		gabcd_mat1 =gabcd_mat.ctypes.data_as(POINTER(c_double))
		C_mat = np.asarray(C_mat, order='C')
		C_mat1 =C_mat.ctypes.data_as(POINTER(c_double))

		# 由原子轨道双电子积分计算分子轨道算电子积分   非常耗时
		cfun.aomo1(gabcd_mat1,go_mat1,C_mat1,c_int(basis_num))
		
		
		
		emp2=0.0

		# 计算mp2校正能
		for i in range(n_electron//2):
			for j in range(n_electron//2):
				for k in range(n_electron//2,basis_num):
					for l in range(n_electron//2,basis_num):
						emp2+= go_mat[i,k,j,l]*(2* go_mat[k,i,l,j]- go_mat[l,i,k,j])/\
						(e_mat[i]+e_mat[j]-e_mat[k]-e_mat[l])
		print("mp2校正能：%10.8f               mp2计算时间为%fs"  %(emp2,time.time()-time2))


	print("\n计算总时间为%fs" %(time.time()-time0))
	print("计算完成")




print('''
###############################################
#                                             # 
#     a python code for ristricted hf/ks      #  
#                                             #
###############################################
''')




if len(sys.argv)==1:
	print("没有输入文件！")
	exit()
elif len(sys.argv)==2:
	if os.path.isfile(sys.argv[1]):
		file = sys.argv[1]
		atom_xyz,dft,mp2,basis_type,contracted= input_pp(file)
else:
	print("输入文件不存在！")
	exit()

#输出部分计算设置信息
if len(atom_xyz)>1:
	print("输入结构如下(bohr)")
	for i in atom_xyz:
		print("%s      %-6.4f     %-6.4f    %-6.4f"  %tuple(i))
else:
	print("单个原子或没有原子！")




basis_data=basis_data_p(basis_type)

# 获得基组信息
ba,nu,contract_list=basis_list(atom_xyz, basis_data)
n_electron=int(np.sum(np.array(nu)[:,3]))

if n_electron//2 ==1:
	print("输入构型总电子数为奇数，rhf或rks方法无法计算")
	exit()
else:
	print("输入构型总电子数为%d"  %(n_electron))


print("\n原始基函数数目为 %d " %(len(ba)))
if contracted==1:
	print("使用收缩基进行计算, 收缩后为%d"  %(len(contract_list)))


if dft==0:
	print("\n进入hf计算...")
else:
	print("\n进入DFT计算中(泛函为最简单的XAplha)")


ba=np.array(ba)
nu=np.array(nu)


scf(ba,nu,contract_list,dft=dft,mp2=mp2,contracted=contracted)



