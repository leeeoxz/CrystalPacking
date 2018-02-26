import os
import json

nonbonddict = {}

with open('nonbondedgromos54a7.txt','r') as f:
	for line in f:
		if (line[0] != '[') & (line[0] != ';'):
			os.system('/usr/local/gromacs/bin/gmx sigeps -c6 '+line.split()[3]+' -cn '+line.split()[4]+" > sigma.txt")
			os.system('rm potje.xvg')
			with open('sigma.txt','r') as f1:
				for linha in f1:
					if len(linha.split()) > 0:
						if linha.split()[0] == 'sigma':
							#print linha.replace(","," ").split()[2],line.split()[0],line.split()[1]
							nonbonddict[line.split()[0]+' '+line.split()[1]]=linha.replace(","," ").split()[2]
							nonbonddict[line.split()[1]+' '+line.split()[0]]=linha.replace(","," ").split()[2]
			os.system('rm sigma.txt')

dic_string = json.dumps(nonbonddict)
print dic_string

with open('dicionario1.json',"w") as f:
	f.write(dic_string)