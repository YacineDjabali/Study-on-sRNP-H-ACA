#!/usr/bin/env python3

"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

Propos :
Donne la frequence et temps de contact des residus prensent a l'interface du complexes dans un fichier texte.

Pour cela on parse le PDB dans un dictionnaire et on le parcours à l'aide de boucles.
L'occurence des residu appartenant à l'interface est ensuite relevée dans un dictionnaire.
Un histogramme est construit à l'aide d'un dictionnaire des frequences.

"""

################# GESTION DES ARGUMENTS ######################

import argparse # Module de gestion des arguments

arg = argparse.ArgumentParser()
arg.add_argument("i", help="Chemin du fichier d'input") # ajout de l'argument input 
arg.add_argument("outb", help = "Nom du fichier d'output des resultats bruts (de preference .txt)" ) # ajout de l'argument output pour les resultats bruts
arg.add_argument("outh", help = "Nom du fichier d'output de l'histogramme (de preference .png)") # ajout de l'argument output pour les resultats de l'histogrammes 
args = arg.parse_args() # recuperation des arguments
##############################################################


######################## MAIN ################################

import sys, time, func_parse, func_distMD, func_interface

start_time = time.time() # start du timer (pour le temps d'execution)

fichierINTERFACE = open(args.outb, "w") # ouverture du fichier de resultat brut
infile = args.i # recuperation de l'input
dPDB = func_parse.ParsePDBDynamique(infile) # creation du dictionnaire PDB

tmpDist = 0 # variable qui contiendra les distance entre residues

nbConformation = 0 # variable qui determinera le nombre de conformation dans le dictionnaire
interface = {} # dictionnaire qui contiendra les residues faisant partie de l'interface
interface["idProt"] = []

seuil = 9 # Seuil a partir du quel on considere un residu comme faisant partie de l'interface


for i in dPDB["models"]:   #PARCOURS DU DICO PDB
	nbConformation += 1
	for j in dPDB[i]["chains"]:
		if (j == "B" or j == "C"): # si on se trouve dans la chaine C (cation Metallique Zn) ou la chaine B (ARN) on break la boucle
			break
		else:
			for k in dPDB[i][j]["resNumber"]:
				for l in dPDB[i]["B"]["resNumber"]:
					tmpDist =  func_distMD.distRes(dPDB[i][j][k], dPDB[i]["B"][l], "barycentre") # recuperation de la distance entre le residue k(proteine) et le residue l(ARN) par la methode du barycentre
					if(func_interface.isResInterface(tmpDist,seuil)): #si tmpDist < seuil
						idProt = "{}{}_{}".format(dPDB[i][j][k]["resName"], k, j)  # Recuperation de l'ID du residue
						if not idProt in interface["idProt"]:
                             #On initialise le residiue si celui n'etais pas present dans le dico et on lui attribue la valeur de 1
							interface["idProt"].append(idProt)
							interface[idProt] = 1
						else :
							interface[idProt] += 1 # Sinon on incremente de 1 sa presence a l'interface
						break #Interrompt la boucle k car le residue j se trouve à l'interface, rien ne sert de continuer car cela va apporter des doublons



Freq = {} #creation du dictionnaire frequence
Freq["idProt"] = [] 
for i in interface["idProt"]:
    Freq["idProt"].append(i)
    Freq[i] = func_interface.frequency(interface[i],nbConformation) # ajout de la frequence d'appartenance a l'interface du residue

##################################################################


###################### PREPARATION DES OUTPUTS ###################
count = 0 # var qui nous servira pour l'axe des abscisses de l'histogramme
tabLabel = [] # recuperation des noms des residues
frequency = [] # tableau des frequence

#residus cle dont on doit calculer le temps de contact
resCle = ["ARG34_A3", "ASP63_A4", "ALA102_A4", "THR41_A4", "ALA100_A4", "GLU98_A4", "GLU43_A4", "LYS39_A4", "ARG38_A4", "LYS46_A4", "ARG50_A4", "ARG6_23", "LYS26_23", "GLU65_A4", "LYS28_A3", "TYR44_A3", "ALA77_A3", "ARG47_A3", "GLU66_A4", "TYR41_A3"]

for i in Freq["idProt"]:
    count += 1
    tabLabel.append(i)
    frequency.append(Freq[i]/100)
    if i in resCle: # si le res "i" est un residu cle
    	fichierINTERFACE.write("{} : {:.2f}%\ttemps de contact : {:.2f}ns\n".format(i,Freq[i],Freq[i]*10/100))  #
    else: #sinon
    	fichierINTERFACE.write("{} : {:.2f}%\n".format(i,Freq[i]))

fichierINTERFACE.write("\n")

fichierINTERFACE.close()#fermeture de l'output

import numpy as np
import matplotlib.pyplot as plt


x = np.arange(count)
fig = plt.figure(figsize=(6, 3), dpi = 200)
plt.title("Frequence relative des residues appartenant à l'interface",
  fontsize=10)

plt.bar(x, frequency)
																								#Creation de l'histogramme
plt.xticks(x, tabLabel, rotation='vertical', fontsize= 4)
plt.yticks()

plt.savefig(args.outh,figsize=(6, 3), dpi = 200)


print("Les resultats sont enregistrés dans le fichier \"{}\"".format(args.outb))
print("L'histogramme est enregistrés dans le fichier \"{}\"".format(args.outh))
print("\n")

print("Temps d'execution : {:.1f} secondes ---" .format(time.time() - start_time))
###################################################################