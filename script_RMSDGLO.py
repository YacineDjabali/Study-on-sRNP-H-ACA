#!/usr/bin/env python3

"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

Propos :
Donne la RMSD global de la partie "proteique" du complexe pour toutes les conformations presentes
dans le pdb dans un fichier texte. 
Et construit a partir de ces RMSD une carte de covariances 

Pour cela on parse le PDB dans un dictionnaire et on le parcours à l'aide de boucles.
Les RMSD pour le fichier texte seront stockées dans un dictionaire dont les clé seront des tuble (confo1,confo2).
Pour la carte de covariance elles seront stockées dans une matrice

"""

################# GESTION DES ARGUMENTS ######################

import argparse # Module de gestion des arguments

arg = argparse.ArgumentParser()
arg.add_argument("i", help="Chemin du fichier d'input") # ajout de l'argument input 
arg.add_argument("outb", help = "Nom du fichier d'output des resultats bruts (de preference .txt)" ) # ajout de l'argument output pour les resultats bruts
arg.add_argument("outg", help = "Nom du fichier d'output des resultats graphiques (de preference .png)") # ajout de l'argument output pour les resultats graphiques 
args = arg.parse_args()# recuperation des arguments

##############################################################


######################## MAIN ################################

import sys, time, func_parse, func_RMSD, func_distMD

start_time = time.time() # start du timer (pour le temps d'execution)

fichier = open(args.outb, "w") # ouverture du fichier de resultat brut
infile = args.i # recuperation de l'input
dPDB = func_parse.ParsePDBDynamique(infile) # creation du dictionnaire PDB

r = {} # dictionnaire qui contiendra les RMSD pour les resultat brut
r["modName"] = []
heatmapMatrice = [] # matrice qui detiendra les RMSD pour la heatmap

for ConfoI in dPDB["models"]:  #PARCOURS DU DICO PDB POUR LA PREMIERE CONFORMATION
	tmpTab = [] #sous tableau pour la matrice
	for ConfoII in dPDB["models"]: #PARCOURS DU DICO PDB POUR LA DEUXIEME CONFORMATION
		s = 0 #variable qui contiendra la somme des distances au carré 
		N = 0 #variable qui contiendra le nombre de pair d'atome dont la RMSD a ete calculé
		dist = [] #tableau qui contiendra toutes les distance entre les conformations
		for k in dPDB[ConfoII]["chains"]:
			if(k == "B" or k == "C"):
				break # On break la boucle si on est dans la chaine B ou C car ces chaines ne correspondent pas a la proteine mais a l'ARN
			else:
				indice = 0 # variable qui permettra de recuperer l'indice du tableau des distance pour le calcule de s
				for l in dPDB[ConfoII][k]["resNumber"]:
					dist.append(func_distMD.distanceCA(dPDB[ConfoI][k][l],dPDB[ConfoII][k][l])) # recuperation de la distance entre les residue l des deux conformations
					s += dist[indice]**2 # somme de la distance precedement calculé au carré
					indice += 1
			
		N = len(dist)
		modtuple = (ConfoI,ConfoII) # creation d'un tuple des conformations traitée qui sera clé du dictionnaire des RMSD
		r["modName"].append(modtuple) 
		r[modtuple] = func_RMSD.RMSD(s,N) # recuperation de la RMSD dans le dictionnaire
		tmpTab.append(func_RMSD.RMSD(s,N)) # recuperation de la RMSD dans le tableau

	heatmapMatrice.append(tmpTab) # recuperation du tableau

##################################################################


###################### PREPARATION DES OUTPUTS ###################

for i in r["modName"]:
    fichier.write("{} : {:.2f}\n".format(i, r[i]))  # Mise en page du fichier des resultats brut

fichier.close()

import numpy as np
import matplotlib.pyplot as plt

data = np.array(heatmapMatrice)
fig, axis = plt.subplots()										#Creation de la matrice de covariance

plt.title("Global RMSD Heatmap")
heatmap = axis.pcolor(data, cmap=plt.cm.Blues)
plt.savefig(args.outg)

##################################################################

print("Les resultats sont enregistrés dans le fichier \"{}\"".format(args.outb))
print("La matrice de covariance est enregistrés dans le fichier \"{}\"".format(args.outg))
print("\n")
print("Temps d'execution : {:.1f} secondes ---" .format(time.time() - start_time))
