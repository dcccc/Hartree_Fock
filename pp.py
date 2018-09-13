#coding:utf-8


def basis_data_p(basis_file):

	basis_file=open(basis_file,"r").readlines()
	line_num=len(basis_file)

	basis_data={}
	n=0
	while n<line_num :
		line=basis_file[n].strip().split()
		if len(line) ==2 and len(line[0])<=2:
			temp={}
			element_m=line[0]
			n+=1
			line=basis_file[n].strip().split()
			while len(line)==3 and len(line[0])<=2 :
				obital_m=line[0]
				temp1=[]
				n_tem=int(line[1])
				for i in range(int(line[1])):
					line=basis_file[n+i+1].strip().split()
					temp1.append(line)
				n=n+n_tem+1
				temp[obital_m]=temp1

				if n<line_num:
					line=basis_file[n].strip().split()
			basis_data[element_m]=temp
		else:
			n+=1

	return(basis_data)
