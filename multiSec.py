#LEONARDO ALVES SANTOS
#OCTOBER 2017

from joblib import Parallel,delayed
import multiprocessing
import time
import os
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import math



infoTAGstride = 'ASG'
ss = 'Coil'
num_cores = multiprocessing.cpu_count()
count = 0
path = "/home/leonardo/Documents/Estagio/PDB/"
stride = "stride "
dssp =  "dssp "
files = os.listdir(path)
dsspabrv = {'H':'AlphaHelix','B':'Bridge', 'E':'Strand', 'G':'310helix', 'I':'PIHelix','T':'Turn','S':'Bend'}
stridereplace = []
dsspreplace = []
aa_VdW = {'H':1.11,'C':1.68,'N':1.53,'O':1.50,'S':1.82, 'Z': 2.09} # VdW distance
crystalpacking = []

##########################################
###   CREATES THE ASSYMETRICAL UNITY   ###
##########################################

def assy(file):
	cont = 0 
	os.system("pymol "+path+file+"/"+file.lower()+".pdb -d 'symexp sym, "+file.lower()+", ("+file.lower()+"), 5; delete "+file.lower()+"; save "+path+file+"/"+file.lower()+"Crystal.pdb' -c")
	cont +=1
	print cont

Parallel (n_jobs=num_cores) (delayed(assy) (f) for f in files)
#################################
###   GET THE COIL RESIDUES   ###
#################################

##################
###   STRIDE   ###
##################

def getSSstride(std):
	content = []
	dic = {}
	dictot = {}
	residues = [] #get the residues number under the ss assignment
	with open(path+std+"/"+std.lower()+".stride","r") as pdb:
		for line in pdb:
			if (line.split()[0] == infoTAGstride):
				dictot[line.split()[2],line.split()[3]]=line
				residue = line.split()[3]
				chain = line.split()[2]
				if line.split()[6] == ss: #if the residue is under the ss
					dic[line.split()[2],line.split()[3]]=line
					residues.append([residue,chain]) #save the residue name
					infoAA = []
					infoAA.append(line.split()[1]) # GET THE RESIDUE 
					infoAA.append(line.split()[2]) # GET THE CHAIN
					infoAA.append(line.split()[3]) # GET THE POSITION NUMBER
					infoAA.append(line.split()[6]) # GET THE SECONDARY STRUCTURE 
					content.append(infoAA)
	with open(path+std+"/"+std.lower()+ss+"_STRIDE.txt","w") as f:
		f.write("RESIDUE\tCHAIN\tNUMBER\tSS\n")
		for aa in content:
			f.write("\t".join(aa)+"\n")
	return residues,dic,dictot


################
###   DSSP   ###
################

def getSSdssp(dsp):
	with open(path+dsp+"/"+dsp.lower()+".dssp",'r') as file:
		start = False # filter the file
		content = []
		dic = {} #creates a dictionary with the line corresponded to the aa
		dictot = {}
		residues = [] #get the residues number under the ss assignment
		for line in file:
			if start:
				content.append(line)
			if line.split()[0] == '#':
				start = True #start the ss 

		#CREATES THE OUTPUT #

		with open(path+dsp+"/"+dsp.lower()+ss+"_DSSP.txt","w") as f:
			f.write('POSITION CHAIN RESIDUE SS\n')
			for item in content:
				try:
					dictot[item[11],str(eval(item[5:10]))]=item
				except:
					pass
				residue = item.split()[1]
				chain = item.split()[2]
				if item[16] == " " and item[13] != "!":
					aux = list(item)
					aux[16] = 'C'
					item = "".join(aux)
					dic[item[11],str(eval(item[5:10]))]=item
					residues.append([residue,chain]) #get the residue
					f.write(item[6:14]+" C\n")
		return residues,dic,dictot


#############################################
###   SETS THE DIFFERENCE BETWEEN LISTS   ###
#############################################

def diffbetweenlist(d1,d2,l1,l2):
	notin1 = []
	notin2 = []
	both = []
	for item in d1:
		if item in d2:
			both.append([d1[item],d2[item]])
		else:
			try:
				notin2.append([d1[item],l2[item]])
			except:
				notin2.append([d1[item],'not in STRIDE                                         '])
	for item in d2:
		if item in d1:
			both.append([d1[item],d2[item]])
		else:
			try:	
				notin1.append([d2[item],l1[item]])
			except:
				notin2.append([d2[item],'not in DSSP                                           '])
	return notin1,notin2,both


#############################################
###   ASSIGNMENT OF SECONDARY STRUCTURE   ###
#############################################

def secondary(file):
	
	residues_dssp = None #Zera a lista de residuos do dssp
	residues_stride = None #Zera a lista de residuos do stride
	aa_number = []
	#count += 1 
	#print count,file
	pdb = file.lower()+'.pdb'
	os.system(stride+path+file+"/"+pdb+" > "+path+file+"/"+file.lower()+".stride") #SET SS WITH STRIDE
	os.system(dssp+" -i "+path+file+"/"+pdb+" -o "+path+file+"/"+file.lower()+".dssp") #SET THE SS WITH DSSP
	residues_dssp,dic_dssp,dic_all_dssp = getSSdssp(file) # RETURNS A LIST OF RESIDUES AND A DICTIONARY WITH ALL RESIDUES PER CHAIN INFORMATION
	residues_stride,dic_stride,dic_all_stride = getSSstride(file)# RETURNS A LIST OF RESIDUES AND A DICTIONARY WITH ALL RESIDUES PER CHAIN INFORMATION
	notindssp,notinstride,both = diffbetweenlist(dic_dssp,dic_stride,dic_all_dssp,dic_all_stride)
	with open (path+file+"/"+file.lower()+'_notinSTRIDE.txt',"w") as f:
		f.write("DSSP/STRIDE\n\n")
		for item in notinstride:
			stridereplace.append(item[1][24])
			#print item[1]
			aa_number.append(item[1][9:15])
			f.write(item[0])
			f.write(item[1])
			f.write(item[0][16]+" > "+item[1][24])
			f.write("\n")
			f.write("\n")
	with open (path+file+"/"+file.lower()+'_notinDSSP.txt',"w") as f:
		f.write("STRIDE/DSSP\n\n")
		for item in notindssp:
			dsspreplace.append(item[1][16])
			#print item[0]
			aa_number.append(item[0][9:15])
			f.write(item[0])
			f.write(item[1])
			f.write(item[0][24]+" > "+item[1][16])
			f.write("\n")
			f.write("\n")
	with open(path+file+"/"+file.lower()+'_residuos.txt',"w") as f:
		for item in aa_number:
			f.write(str(item)+"\n")
	return stridereplace,dsspreplace

r = Parallel (n_jobs=num_cores) (delayed(secondary) (f) for f in files)

for item in r:
	for i in item[0]:
		stridereplace.append(i)
	for i in item[1]:
		dsspreplace.append(i)

c = Counter(stridereplace)

s = Counter(dsspreplace)

with open ('/home/leonardo/Documents/Estagio/DSSPTROCOU.txt',"w") as f:
	for item in s:
		f.write(item+"\t"+str(s[item])+"\n")	

with open ('/home/leonardo/Documents/Estagio/STRIDETROCOU.txt',"w") as f:
	for item in c:
		f.write(item+"\t"+str(c[item])+"\n")


##############################
###   EUCLIDIAN DISTANCE   ###
##############################

def euclidian(x1,y1,z1,x2,y2,z2): #calculate the euclidian distance between two points in 3D space
	return math.sqrt(((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))

#################################
###   FIND CRYSTAL PACKING   ####
#################################

def cryst(file):
	cont = 0
	residues = [] #all target residues
	x_max,y_max,z_max = -1000,-1000,-1000 #set the max value to a lower value, so we make sure that none atom is left aside
	x_min,y_min,z_min = 1000,1000,1000 #set the mmin value to a bigger value, so we make sure that none atom is left aside
	aux = []
	pdb_original = []
	pdb_crystal = []
	a = []
	with open(path+file+"/"+file.lower()+"_residuos.txt") as f:
		for line in f:
			residues.append(line.split())

	with open(path+file+"/"+file.lower()+".pdb") as f:
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
	
	for item in pdb_original:
		for item2 in pdb_crystal:
			try:
				omega = aa_VdW[item[13]]*2
			except:
				omega = 2.22
			distance = euclidian(eval(item[30:38]),eval(item[38:46]),eval(item[46:54]),eval(item2[30:38]),eval(item2[38:46]),eval(item2[46:54]))
			if (omega >= distance):
				cry = file+" distance: "+str(distance)+"\n"+item+"\n"+item2+"\n\n\n"
				crystalpacking.append(cry)			
	
	print file,cont
	return crystalpacking

results = Parallel (n_jobs=num_cores) (delayed(cryst) (f) for f in files)

with open("crystapack.txt","w") as f:
	for items in results:
		for item in items:
			f.write(item)

print 'Terminou'
print time.clock()
