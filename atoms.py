#import json
import os
#from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import math

infoTAGstride = 'ASG'
ss = 'Coil'
cont = 0
count = 0
path = "/home/leonardo/Documents/Estagio/Teste/"
stride = "stride "
dssp =  "dssp "
files = os.listdir(path)
dsspabrv = {'H':'AlphaHelix','B':'Bridge', 'E':'Strand', 'G':'310helix', 'I':'PIHelix','T':'Turn','S':'Bend'}
stridereplace = []
dsspreplace =[]

aa_VdW = {'H':1.11,'C':1.68,'N':1.53,'O':1.50,'S':1.82, 'Z': 2.09} # VdW distance


def euclidian(x1,y1,z1,x2,y2,z2): #calculate the euclidian distance between two points in 3D space
	return math.sqrt(((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))

f = open("1.","w")

for file in files:
	residues = [] #all target residues
	x_max,y_max,z_max = -1000,-1000,-1000 #set the max value to a lower value, so we make sure that none atom is left aside
	x_min,y_min,z_min = 1000,1000,1000 #set the mmin value to a bigger value, so we make sure that none atom is left aside
	aux = []
	pdb_original = []
	pdb_crystal = []

	with open(path+file+"/"+file.lower()+"_residuos.txt") as f:
		for line in f:
			residues.append(line.split())

	with open(path+file+"/"+file.lower()+"hbond.pdb") as f:
		for item in f:
			if 'ATOM' in item[0:6]: #if its atom, then
				if [item[21],item[22:26].split()[0]] in residues:
					if (eval(item[30:38]) > x_max):
						x_max = eval(item[30:38])
					if (eval(item[30:38]) < x_min):
						x_min = eval(item[30:38])
					if (eval(item[38:46]) > y_max):
						y_max = eval(item[38:46])
					if (eval(item[38:46]) < y_min):
						y_min = eval(item[38:46])
					if (eval(item[46:54]) > z_max):
						z_max = eval(item[46:54])
					if (eval(item[46:54]) < z_min):
						z_min = eval(item[46:54])
					pdb_original.append(item)
	
	
	with open(path+file+"/"+file.lower()+"Crystal.pdb") as f:
		for line in f:
			if ('ATOM' in line[0:6]) | ('HETATM' in line[0:6]):
				if (eval(line[30:38]) > x_min-5) & (eval(line[30:38]) < x_max+5):
					if (eval(line[38:46]) > y_min-5) & (eval(line[38:46]) < y_max+5):
						if (eval(line[46:54]) > z_min-5) & (eval(line[46:54]) < z_max+5):
							pdb_crystal.append(line)
	print file
	
	for item in pdb_original:
		for item2 in pdb_crystal:
			try:
				omega = aa_VdW[item[13]]*2
			except:
				omega = 4
				#print item[12:16]
			distance = euclidian(eval(item[30:38]),eval(item[38:46]),eval(item[46:54]),eval(item2[30:38]),eval(item2[38:46]),eval(item2[46:54]))
			if (omega >= distance):
				#print distance,omega
				#print "original",item,"\n\n","crystal",item2,"\n\n\n\n"	
				dsspreplace.append(file)			
	cont += 1
	print file,cont
	
	"""for item1 in pdb_original:
		for item2 in pdb_crystal:
			dist = euclidian(eval(item1[30:38]),eval(item1[38:46]),eval(item1[46:54]),eval(item2[30:38]),eval(item2[38:46]),eval(item2[46:54]))
			if dist < 1:
				print dist"""
for item in set(dsspreplace):
	print item

with open("files_crystal.txt","w") as f:
	for item in set(dsspreplace):
		f.write(item)
		f.write("\n")
