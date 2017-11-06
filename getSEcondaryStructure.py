import os

infoTAG = 'ASG'
ss = 'Coil'


path = "/home/leonardo/Documents/Estagio/PDB/"

files = os.listdir(path)

for file in files:
	content = []
	with open(path+file+"/"+file.lower()+".txt","r") as pdb:
		for line in pdb:
			if (line.split()[0] == infoTAG):
				if line.split()[6] == ss:
					infoAA = []
					infoAA.append(line.split()[1])
					infoAA.append(line.split()[2])
					infoAA.append(line.split()[3])
					infoAA.append(line.split()[6])
					content.append(infoAA)
					print infoAA
	with open(path+file+"/"+file.lower()+ss+".txt","w") as f:
		for aa in content:
			f.write("\t".join(aa)+"\n")