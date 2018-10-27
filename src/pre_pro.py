#coding=utf-8

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
				if obital_m in temp.keys():
					temp[obital_m]=temp[obital_m]+[temp1]
				else:
					temp[obital_m]=[temp1]

				if n<line_num:
					line=basis_file[n].strip().split()
			basis_data[element_m]=temp
		else:
			n+=1

	return(basis_data)












def basis_list(coor,basis_data):
	atom_in={"H":1.0,"He":2.0,"Li":3.0,"Be":4.0,"B":5.0,"C":6.0,"N":7.0,"O":8.0,"F":9.0,"Ne":10.0}
	ba=[]
	nu=[]
	contract_list=[]

	for i in coor :

		if i[0] in ["Li","Be","B","C","N","O","F","Ne"]:
			temp1=[i[1],i[2],i[3],atom_in[i[0]]]
			nu.append(temp1)
			temp2=[]
			basis_data_t=basis_data["C"]["S"]			
			for ii in basis_data_t:
				temp2.append(len(ii))
				for j in ii:
					temp=[i[1],i[2],i[3],0,0,0,float(j[0]),float(j[1])]		
					ba.append(temp)

			basis_data_t=basis_data["C"]["SP"]

			for ii in basis_data_t:
				temp2.append(len(ii))
				for j in ii:
					temp=[i[1],i[2],i[3],0,0,0,float(j[0]),float(j[1])]
					ba.append(temp)
				
			for k in [[1,0,0],[0,1,0],[0,0,1]]:
				for ii in basis_data_t:	
					temp2.append(len(ii))
					for j in ii:
						temp=[i[1],i[2],i[3],k[0],k[1],k[2],float(j[0]),float(j[2])]
						ba.append(temp)							

			contract_list.append(temp2)

		if i[0] in ["H","He"]:
			temp2=[]
			temp1=[i[1],i[2],i[3],atom_in[i[0]]]
			nu.append(temp1)
			basis_data_t=basis_data["H"]["S"]
			for ii in basis_data_t:
				temp2.append(len(ii))
				for j in ii :					
					temp=[i[1],i[2],i[3],0,0,0,float(j[0]),float(j[1])]
					ba.append(temp)

			contract_list.append(temp2)

	return(ba,nu,contract_list)


def xyz_file_pp(xyz_file):

	xyz_file=open(xyz_file).readlines()

	atom_num=int(xyz_file[0].strip())
	
	atom_xyz=[]
	bohr=1.8897261635610906894703702275911
	for line in xyz_file[2:]:
		line=line.strip().split()
		if len(line)==4:
			temp=[line[0],float(line[1])*bohr,float(line[2])*bohr,float(line[3])*bohr]
			atom_xyz.append(temp)
	return(atom_xyz)



def  input_pp(file):
	file=open(file,"r").readlines()
	keywordline = file[0]

	contracted =0
	mp2=0
	dft=0

	if "mp2" in keywordline:
		mp2=1
	elif "dft" in keywordline:
		dft=1
	elif "hf" in keywordline :
		dft=0
	if "contracted" in keywordline :
		contracted=1

	if "dft" in keywordline and "hf" in keywordline:
		print(" 'hf' 和 'dft' 关键词只能同时存在一个!")
		exit()
	elif "dft" not in keywordline and "hf" not in keywordline:
		print(" 'hf' 和 'dft' 关键词必须指定一个!")
		exit()

	if "dft" in keywordline and "mp2" in keywordline:
		print("####")
		print("####在dft计算得到的分子轨道基础三做mp2计算，其校正能意义？")
		print("####")

	if "631g" in keywordline or "6-31g" in keywordline :
		basis_type="6-31g.txt"
	elif "sto3g" in keywordline or "sto-3g" in keywordline :
		basis_type="sto3g.txt"
	elif "3-21g" in keywordline or "321g" in keywordline :
		basis_type="3-21g.txt"

	atom_xyz=[]
	bohr=1.8897261635610906894703702275911
	for line in file[2:]:
		line=line.strip().split()
		if len(line)==4:
			temp=[line[0],float(line[1])*bohr,float(line[2])*bohr,float(line[3])*bohr]
			atom_xyz.append(temp)

	return(atom_xyz,dft,mp2,basis_type,contracted)
