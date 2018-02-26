import os

path = "/home/leonardo/Documents/Estagio/Teste/"
files = os.listdir(path)

for item in files:
	content = os.listdir(path+item)
	for file in content:
		if file == item.lower()+".pdb":
			pass
		else:
			os.system("rm "+path+item+"/"+file)
