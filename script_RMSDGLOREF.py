#!/usr/bin/env python3

"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

Propos :
Donne la RMSD global de la partie "proteique" du complexe contre une conformations de references (model 10) dans un fichier texte. 
Et construit a partir de ces RMSD un graphique, "conformations en fonctions des RMSD"

Pour cela on parse le PDB dans un dictionnaire et on le parcours à l'aide de boucles.
Les RMSD pour le fichier texte seront stockées dans un dictionaire dont les clé seront des tuples (reference,confo).
La construction du graphique se fais par recuperation des donnée en abscisse (RMSD) et en ordonnée(Conformations)
à l'aide du dictionnaire precedemment obtenu
"""

################# GESTION DES ARGUMENTS ######################

import argparse

arg = argparse.ArgumentParser()
arg.add_argument("i", help="Chemin du fichier d'input") # ajout de l'argument input 
arg.add_argument("outb", help = "Nom du fichier d'output des resultats bruts (de preference .txt)" )  # ajout de l'argument output pour les resultats bruts
arg.add_argument("outg", help = "Nom du fichier d'output des resultats graphiques (de preference .png)")# ajout de l'argument output pour les resultats graphiques
args = arg.parse_args()# recuperation des arguments

##############################################################


######################## MAIN ################################

import sys, time, func_parse, func_RMSD, func_distMD

start_time = time.time() # start du timer (pour le temps d'execution)

fichier = open(args.outb, "w")# ouverture du fichier de resultat brut
infile = args.i# recuperation de l'input
dPDB = func_parse.ParsePDBDynamique(infile)# creation du dictionnaire PDB

r = {}# dictionnaire qui contiendra les RMSD pour les resultat brut
r["modName"] = []


for Confo in dPDB["models"]:#PARCOURS DU DICO PDB POUR LES CONFORMATIONS A COMPARER AVEC LA REFERENCE
	N = 0 #variable qui contiendra le nombre de pair d'atome dont la RMSD a ete calculé
	s = 0 #variable qui contiendra la somme des distances au carré 
	dist = [] #tableau qui contiendra toutes les distance entre les conformations
	for k in dPDB[Confo]["chains"]:
		if(k == "B" or k == "C"): 
			break# On break la boucle si on est dans la chaine B ou C car ces chaines ne correspondent pas a la proteine mais a l'ARN
		else:
			indice = 0 # variable qui permettra de recuperer l'indice du tableau des distance pour le calcule de s
			for l in dPDB[Confo][k]["resNumber"]:
				dist.append(func_distMD.distanceCA(dPDB["10"][k][l],dPDB[Confo][k][l])) # recuperation de la distance entre les residues l de la conformations 10 et de la conformation qui boucle
				s += dist[indice]**2 # somme de la distance precedement calculé au carré
				indice += 1
	
	N = len(dist)
	modtuple = ("10",Confo) # creation d'un tuple des conformations traitée qui sera clé du dictionnaire des RMSD
	r["modName"].append(modtuple)
	r[modtuple] = func_RMSD.RMSD(s,N) # recuperation de la RMSD dans le dictionnaire

##################################################################


###################### PREPARATION DES OUTPUTS ###################

abscisse = [] #tableau d'abscisse pour le graphique
ordonne = [] #tableau d'ordonnée pour le graphique
xlabels = 0 #recuperation du nombre de conformation
for i in r["modName"]:
	fichier.write("{} : {}\n".format(i, r[i])) # Mise en page du fichier des resultats brut
	xlabels += 1
	abscisse.append(xlabels)
	ordonne.append(r[i])

fichier.close()

import numpy as np
import matplotlib.pyplot as plt

x = np.array(abscisse)
y = np.array(ordonne)

plt.plot(x,y)
plt.title("RMSD global contre la conformation de reference")
plt.xlabel("Conformations")                                                #Creation du graphique
plt.ylabel("RMSD")
plt.savefig(args.outg)

##################################################################

print("Les resultats sont enregistrés dans le fichier \"{}\"".format(args.outb))
print("Le plot est enregistrés dans le fichier \"{}\"".format(args.outg))
print("\n")
print("Temps d'execution : {:.1f} secondes ---" .format(time.time() - start_time))
