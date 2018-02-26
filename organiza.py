import os



aa_VdW = {'H':1.11,'C':1.68,'N':1.53,'O':1.50,'S':1.82, 'Z': 2.09}
path = '/home/leonardo/Documents/Estagio/PDB/'
one = os.listdir(path)
files= []
critico = {}
dic = {}

with open('crystapack.txt','r') as f:
	for line in f:
		if line != ('\n'):
			files.append(line)
atoms = []

for x in xrange(0,len(files),3):
	atoms.append(files[x:x+3])
teste =[]
for atom in atoms:
	if atom[0].split()[0] in dic:
		dic[(atom[0].split()[0])].append(atom)
	else:
		dic[(atom[0].split()[0])] = []
		dic[(atom[0].split()[0])].append(atom)
	if atom [1][13] in aa_VdW:
		if eval(atom[0].split()[2]) < aa_VdW[atom[1][13]]:	
			if atom[0].split()[0] in critico:
                		critico[(atom[0].split()[0])].append(atom)
        		else:
                		critico[(atom[0].split()[0])] = []
                		critico[(atom[0].split()[0])].append(atom)
	else:
		if eval(atom[0].split()[2]) < aa_VdW[atom[1][12]]:
                        if atom[0].split()[0] in critico:
                                critico[(atom[0].split()[0])].append(atom)
                        else:
                                critico[(atom[0].split()[0])] = []
                                critico[(atom[0].split()[0])].append(atom)
with open('probCryPack.txt','w') as f:
	for item in dic:
		f.write(item+'\n')
		with open(path+item+'/crypack.txt','w') as f1:
			for atom in dic[item]:
				for info in atom:
					f1.write(info)
					f1.write('\n')
with open('critico.txt','w') as f:
        for item in critico:
                f.write(item+'\n')
                with open(path+item+'/critico.txt','w') as f1:
                        for atom in critico[item]:
                                for info in atom:
                                        f1.write(info)
                                        f1.write('\n')

print len(critico),len(dic)
