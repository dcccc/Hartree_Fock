#coding=utf-8
import math,time,sys,os
import numpy as np
from grid_int import   *  # xyzw_list_gen,xc_int_c, density_pot_list_c,density_pot_list_c_contra,
from pre_pro import basis_list, xyz_file_pp,basis_data_p,input_pp
from analy_int import *  #c_normarize,all_energy,pab,g_ab_j,g_ab_jk,contract_int,
from analy_int import smat_pri,tmat_pri,vmat_pri,gmat_pri,smat_con,tmat_con,vmat_con,gmat_con
from ctypes import cdll,c_double,POINTER,c_int


def scf(ba,nu,contract_list,dft=0,mp2=0,contracted=1):

    # number of primitive basis function
    basis_num=len(ba)
    atom_num=len(nu)


    delta_e=1.0
    delta_p=1.0
    n_c=0
    
    time0=time.time()
    
    print("\n")
    print("eri calculating....")


    # when contracted basis is used
    if contracted==1:    
        # the index of contracted basis function
        contract_cum=np.cumsum(np.array(np.array([y for x in contract_list for y in x])))
        contract_cum1=np.insert(contract_cum,0,0)
        # number of contracted basis function 
        basis_num=len(contract_cum)
        # contracted basis function list 
        ba_con = [ba[i:contract_cum1[n+1]] for n,i in enumerate(contract_cum1[:-1])]

        # funtions in c++ libs used 
        if func_c :
            # calculaiton of S matrix
            s_mat_dia=np.zeros((len(ba)))
            s_mat=np.zeros((basis_num,basis_num))
            smat_con(ba,s_mat,s_mat_dia,contract_cum1)
    
            # calculaiton of T matrix
            t_mat=np.zeros((basis_num,basis_num))
            tmat_con(ba,t_mat,s_mat_dia,contract_cum1)
            
            # calculaiton of V matrix
            v_mat=np.zeros((basis_num,basis_num))
            vmat_con(ba,v_mat,s_mat_dia,nu,contract_cum1)
        
            # calculaiton of eri
            gabcd_mat=np.zeros((basis_num,basis_num,basis_num,basis_num))
            gmat_con(ba,gabcd_mat,s_mat_dia,contract_cum1)
        
        # pure python funtions used 
        else:
            # calculaiton of S matrix
            s_mat = get_S_mat_con(ba_con)
            s_mat_dia = [sab(i,i) for i in ba]

            # calculaiton of T matrix
            t_mat = get_T_mat_con(ba_con)
            
            # calculaiton of V matrix
            v_mat = get_V_mat_con(ba_con, nu)
        
            # calculaiton of eri
            gabcd_mat = get_Gabcd_mat_con(ba_con)

    # when primitive basis is used
    else: 
        # funtions in c++ libs used 
        if func_c:
            # calculaiton of S matrix
            s_mat_dia=np.zeros((len(ba)))
            s_mat=np.zeros((basis_num,basis_num))
            smat_pri(ba,s_mat,s_mat_dia)

            # calculaiton of T matrix
            t_mat=np.zeros((basis_num,basis_num))
            tmat_pri(ba,t_mat,s_mat_dia)
            
            # calculaiton of V matrix
            v_mat=np.zeros((basis_num,basis_num))
            vmat_pri(ba,v_mat,s_mat_dia,nu)
            
            # calculaiton of eri
            gabcd_mat=np.zeros((basis_num,basis_num,basis_num,basis_num))
            gmat_pri(ba,gabcd_mat,s_mat_dia)

        # pure python funtions used
        else:
            # calculaiton of S matrix
            s_mat, s_mat_dia = get_S_mat_pri(ba)
            
            # calculaiton of T matrix
            t_mat = get_T_mat_pri(ba,s_mat_dia)
            
            # calculaiton of V matrix
            v_mat = get_V_mat_pri(ba, nu, s_mat_dia)
        
            # calculaiton of eri
            gabcd_mat = get_Gabcd_mat_pri(ba, s_mat_dia)
        
        
    #  calculaiton of S^-0.5
    d_mat,l_mat=np.linalg.eig(s_mat)
    d_mat=np.diag(tuple(d_mat**-0.5))
    s_mat_05=np.dot(np.dot(l_mat,d_mat),l_mat.T)
    
    
    # fock matrix 
    H_mat   = t_mat - v_mat
    F_mat_1 = np.dot(np.dot(s_mat_05 ,H_mat), s_mat_05.T)
    
    
    # diagonalization of fock matrix
    e_mat, C_mat=np.linalg.eig(F_mat_1)
    
    # normalization and sortion of coefficient matirx
    e_mat, C_mat = c_normarize(np.dot(s_mat_05, C_mat), e_mat, s_mat, basis_num)
    
    # density matrix
    p_mat = pab(C_mat, n_electron)
    

    ene=0.0
    
    print('eri calcualtion time ： % 10.3f s'   %(time.time()-time0))
    print("\n")
    print("begin scf calcualtion")
    time1=time.time()
    print("%-s            %-s            %-s      %-s"  %("cycle", "E(hartree)","delta_e(hartree)","delta_p"))
   
    # begain scf calculation 
    while(delta_e>1*10**-8 or delta_p>1*10**-6):
    
    
        # initialize a zero G matrix
        G_mat=np.zeros((basis_num, basis_num))        


        # do the DFT calculation 
        if dft==1:
            
            # generate grid point for atoms
            xyzw_list=[]            
            xyzw_list=np.array(xyzw_list_gen(atom_xyz))
            contract_cum=np.cumsum(np.array([y for x in contract_list for y in x]))
            basis_num_con=len(contract_list)

            basis_num1=len(ba)
            ks_mat=np.zeros((basis_num1,basis_num1))

            # calculation of coulomb matrix
            for i in range(basis_num):
                for j in range(basis_num):
                    G_mat[i,j]=g_ab_j(p_mat,gabcd_mat,i,j,basis_num)            
            
            # when contracted basis used 
            if contracted==1:
                # initialize a zero density matrix under primitive basis condition, 
                # and map the density matrix under the contracted basis condition to it
                p_mat_uncon=np.zeros((basis_num1,basis_num1))
                contract_cum1=np.insert(contract_cum,0,0)
                for i in range(basis_num):
                    for j in range(basis_num):
                        if i>=j:
                            p_mat_uncon[contract_cum1[i]:contract_cum1[i+1],contract_cum1[j]:contract_cum1[j+1],]=p_mat[i,j]
                            p_mat_uncon[contract_cum1[j]:contract_cum1[j+1],contract_cum1[i]:contract_cum1[i+1],]=p_mat[j,i]


                # calculation of electron density at grid points
                if func_c:
                    density_list=np.array(density_pot_list_c_contra(ba,xyzw_list,s_mat_dia,p_mat_uncon))
                else:
                    density_list=np.array(density_pot_list_contra(ba,xyzw_list,s_mat_dia,p_mat_uncon))

            # when primitive basis used 
            else:
                # calculation of electron density at grid points
                if func_c:
                    density_list=np.array(density_pot_list_c(ba,xyzw_list,s_mat_dia,p_mat))
                else:
                    density_list=np.array(density_pot_list(ba,xyzw_list,s_mat_dia,p_mat))
                p_mat_uncon=p_mat

            #  calculaiton of XC functional potential(Xalpha functional)
            #  page 562 in book of "quantum chemistry" 7th edition by levine 
            for i in range(basis_num1):
                for j in range(basis_num1):
                    if i>=j :                        
                        if func_c:
                            ks_mat[i,j]=xc_int_c(ba,i,j,s_mat_dia,p_mat_uncon,xyzw_list,density_list)\
                            *-(3.0/2)*(3.0/math.pi)**(1.0/3)*0.7            
                        else:
                            ks_mat[i,j]=xc_int(ba,i,j,s_mat_dia,p_mat_uncon,xyzw_list,density_list)\
                            *-(3.0/2)*(3.0/math.pi)**(1.0/3)*0.7
                        ks_mat[j,i]=ks_mat[i,j]
        
        # do HF calculation
        else:
            #  calculation of G matrix
            for i in range(basis_num):
                for j in range(basis_num):
                    G_mat[i,j]=g_ab_jk(p_mat, gabcd_mat, i, j, basis_num)

        if contracted==1 and dft==1:
            #  calculate XC matrix of contracted basis form XC matrix of primitive basis 
            ks_mat=contract_int(ba,contract_cum,ks_mat,np.zeros((basis_num,basis_num)))

        elif dft==0 :            
            ks_mat=np.zeros((basis_num,basis_num))

        
        # fock matrix and its diagonalization
        F_mat=H_mat+G_mat+ks_mat
        F_mat_1=np.dot(np.dot(s_mat_05,F_mat),s_mat_05.T)
        e_mat,C_mat=np.linalg.eig(F_mat_1)
        e_mat,C_mat=c_normarize(np.dot(s_mat_05,C_mat),e_mat,s_mat,basis_num)
        p_mat1=pab(C_mat,n_electron)
    
        # energy calculaiotn of XC functional (Xalpha functional)     
        # eq (16.61) in book of "quantum chemistry" 7th edition by levine 
        if dft==1:
            E_xc=0.0
            for i in range(len(density_list)):
                for j in range(len(density_list[i])):
                    E_xc+=density_list[i][j]**(4.0/3.0)*xyzw_list[i][j][3]

            E_xc*=4*math.pi*-(9.0/8)*(3.0/math.pi)**(1.0/3)*0.7

            # DFT energy of mol
            ene1= all_energy(p_mat,H_mat*2+G_mat,np.array(nu),basis_num)+E_xc
        else:
            # HF  energy of mol
            ene1= all_energy(p_mat,F_mat+H_mat+ks_mat,np.array(nu),basis_num)
    
    
    
        # energy difference between current scf step and before
        delta_e=abs(ene1-ene)
    
        # difference between old and new density matrixs
        delta_p=0.0
        for i in range(basis_num):
            for j in range(basis_num):
                delta_p+=abs(p_mat1[i,j]-p_mat[i,j])/basis_num**2

        # update system energy
        ene=ene1

        # mix of old and new density 
        p_mat=p_mat1*0.5+p_mat*0.5
        n_c=n_c+1
    
    
        # print the system energy and the changes of system energy, density matrix
        print("% 4d            % 10.8f           % 10.8f            % 10.8f"  %(n_c, ene, delta_e, delta_p))

    print("scf total time is % 10.3f s"  %(time.time()-time1))


    # mp2 calculation
    time2=time.time()
    if mp2==1 :#and dft !=1 :
        print("\n")
        print("mp2 calculating ... ")

        # calculate mol orbit eri from mol orbit overlap, quite time comsuming
        # the calculation may take hours, and calculation is feasible only 
        # using c++ lib function 
        # eq (16.6) in book of "quantum chemistry" 7th edition by levine 
        go_mat=np.zeros((basis_num,basis_num,basis_num,basis_num))
        go_mat = np.asarray(go_mat, order='C')
        C_mat = np.asarray(C_mat, order='C')

        go_mat = np.asarray(go_mat, order='C')
        go_mat1=go_mat.ctypes.data_as(POINTER(c_double))
        gabcd_mat = np.asarray(gabcd_mat, order='C')
        gabcd_mat1 =gabcd_mat.ctypes.data_as(POINTER(c_double))
        C_mat = np.asarray(C_mat, order='C')
        C_mat1 =C_mat.ctypes.data_as(POINTER(c_double))


        cfun.aomo1(gabcd_mat1,go_mat1,C_mat1,c_int(basis_num))
        
        
        
        emp2=0.0

        # calculaiton of mp2 correction energy
        # eq (16.13) in book of "quantum chemistry" 7th edition by levine 
        for i in range(n_electron//2):
            for j in range(n_electron//2):
                for k in range(n_electron//2,basis_num):
                    for l in range(n_electron//2,basis_num):
                        emp2+= go_mat[i,k,j,l]*(2* go_mat[k,i,l,j]- go_mat[l,i,k,j])/\
                        (e_mat[i]+e_mat[j]-e_mat[k]-e_mat[l])
        print("mp2 correction is % 10.8f              mp2 calcualting time is % 10.3f s"  %(emp2,time.time()-time2))

    print("\n")
    print("total calcualting time is % 10.3f s" %(time.time()-time0))
    print("finished")




print('''
###############################################
#                                             # 
#     a python code for ristricted hf/ks      #  
#                 by dcccc                    #
#                                             #
###############################################
''')




if len(sys.argv)==1:
    print("No input file ！")
    exit()
elif len(sys.argv)==2:
    if os.path.isfile(sys.argv[1]):
        file = sys.argv[1]
        # read input content
        #
        # atom_xyz   : the mol atom type and coordinate 
        # dft or mp2 : do dft or mp2 calculation or not
        # basis_type : basis set type used in calculation
        # contracted : whether used the contracted basis set 
        atom_xyz,dft,mp2,basis_type,contracted= input_pp(file)
else:
    print("Input file doesn't exist！")
    exit()

# output the coordinate of mol 
if len(atom_xyz)>1:
    print(" structure coordinate is ( unit is bohr)")
    for i in atom_xyz:
        print("%s      % 6.4f     % 6.4f    % 6.4f"  %tuple(i))
else:
    print("no atom or only one atom in the input file ！")



# read corresponding basis set file and store basis set data in a dict
basis_data=basis_data_p(basis_type)

# generate basis function for the mol 
# ba is the basis function list of mol
# nu is the nuclear charge list of atoms in mol
# contract_list, the number of primitive basis function for contracted basis function
ba,nu,contract_list=basis_list(atom_xyz, basis_data)
nu=np.array(nu)
n_electron=int(np.sum(np.array(nu)[:,3]))

# check the electron number of mol, only closed shell system supported 
if n_electron%2 ==1:
    print("number of electrons is odd ，calcualtion stops")
    exit()
else:
    print("number of electrons is % 5d"  %(n_electron))

print("\n")
# output the information about the basis set 
print("number of original basis function is % 4d " %(len(ba)))
if contracted==1:
    print("use contracted basis set, number of contracted basis function is %d"  \
    %(len([y for x in contract_list for y in x])))

print("\n")
# do a HF or DFT calculaiton
if dft==0:    
    print("begin HF scf calcualtion")
else:
    print("begin DFT scf calcualtion(functional is xaplha)")

# whether c++ lib exist or not 
# functions in C++ libs will be used when C++ libs are compiled, the calculation 
# will much faster then using the pure python fuctions when C++ libs missed
func_c = os.path.exists("./analy_int.so1")
if func_c:
    cfun=cdll.LoadLibrary("./grid_int.so")
    ba=np.array(ba)

# do the scf calculation
scf(ba,nu,contract_list,dft=dft,mp2=mp2,contracted=contracted)


