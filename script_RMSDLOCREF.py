#!/usr/bin/env python3

"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

Propos :
Donne la RMSD local de la partie "proteique" du complexe pour toutes les chaines des conformations presentes
dans le pdb contre la reference (model 10) dans un fichier texte.
Et construit a partir de ces RMSD un graphique pour chaques chaines, "conformations en fonctions des RMSD"


Pour cela on parse le PDB dans un dictionnaire et on le parcours à l'aide de boucles.
Les RMSD pour le fichier texte seront stockées dans un dictionaire dont les clé seront des tuples (confo1,confo2).
La construction du graphique se fais par recuperation des donnée en abscisse (RMSD) et en ordonnée(Conformations)
à l'aide du dictionnaire precedemment obtenu

"""

################# GESTION DES ARGUMENTS ######################

import argparse# Module de gestion des arguments

arg = argparse.ArgumentParser()
arg.add_argument("i", help="Chemin du fichier d'input") # ajout de l'argument input 
arg.add_argument("outb", help = "Nom du fichier d'output des resultats bruts (de preference .txt)" ) # ajout de l'argument output pour les resultats bruts
arg.add_argument("outg", help = "Nom du fichier d'output des resultats graphiques (de preference .png)") # ajout de l'argument output pour les resultats graphiques
args = arg.parse_args() # recuperation des arguments

##############################################################


######################## MAIN ################################

import sys, time, func_parse, func_RMSD, func_distMD

start_time = time.time() # start du timer (pour le temps d'execution)

fichier = open(args.outb, "w") # ouverture du fichier de resultat brut
infile = args.i # recuperation de l'input
dPDB = func_parse.ParsePDBDynamique(infile) # creation du dictionnaire PDB

r = {} # dictionnaire qui contiendra les RMSD pour les resultat brut
r["modName"] = []

for Confo in dPDB["models"]:#PARCOURS DU DICO PDB POUR LES CONFORMATIONS A COMPARER AVEC LA REFERENCE
	modtuple = ("10", Confo) # creation d'un tuple des conformations traitée qui sera clé du dictionnaire des RMSD
	r["modName"].append(modtuple)
	r[modtuple] = {}
	r[modtuple]["chains"] = []
	for k in dPDB[Confo]["chains"]:
		N = 0
		s = 0
		dist = []
		if(k == "B" or k == "C"):
			break  # On break la boucle si on est dans la chaine B ou C car ces chaines ne correspondent pas a la proteine mais a l'ARN
		else:
			r[modtuple]["chains"].append(k) # recuperation de la chaine dans le dictionnaire
			indice = 0 # variable qui permettra de recuperer l'indice du tableau des distance pour le calcule de s
			for l in dPDB[Confo][k]["resNumber"]:
				dist.append(func_distMD.distanceCA(dPDB["10"][k][l],dPDB[Confo][k][l])) # recuperation de la distance entre les residues l de la conformations 10 et de la conformation qui boucle
				s += dist[indice]**2 # somme de la distance precedement calculé au carré
				indice += 1

		N = len(dist)
		r[modtuple][k] = func_RMSD.RMSD(s,N) # recuperation de la RMSD dans le dictionnaire

##################################################################

###################### PREPARATION DES OUTPUTS ###################

fichier.write("RMSD LOCAL : CHAINE A1\n")
for i in r["modName"]:
    fichier.write("{} : {:.2f}\n".format(i, r[i]["A1"]))


fichier.write("\n")


fichier.write("RMSD LOCAL : CHAINE A2\n")
for i in r["modName"]:
    fichier.write("{} : {:.2f}\n".format(i, r[i]["A2"]))

fichier.write("\n")

fichier.write("RMSD LOCAL : CHAINE A3\n")
for i in r["modName"]:
    fichier.write("{} : {:.2f}\n".format(i, r[i]["A3"]))

fichier.write("\n")

abscisse = []
xlabels = 0
ordonnee_A1 = []
ordonnee_A2 = []
ordonnee_A3 = []
ordonnee_A4 = []

fichier.write("RMSD LOCAL : CHAINE A4\n")
for i in r["modName"]:
	fichier.write("{} : {:.2f}\n".format(i, r[i]["A4"]))

	xlabels += 1
	abscisse.append(xlabels)
	ordonnee_A1.append(r[i]["A1"])
	ordonnee_A2.append(r[i]["A2"])
	ordonnee_A3.append(r[i]["A3"])                                # Mise en page du fichier des resultats brut + Creation du graphique
	ordonnee_A4.append(r[i]["A4"])

fichier.close()

import numpy as np
import matplotlib.pyplot as plt

x = np.array(abscisse)
y_A1 = np.array(ordonnee_A1)
y_A2 = np.array(ordonnee_A2)
y_A3 = np.array(ordonnee_A3)
y_A4 = np.array(ordonnee_A4)

plt.plot(x,y_A1, label = "RMSD local chaine A1")
plt.plot(x,y_A2, label = "RMSD local chaine A2")
plt.plot(x,y_A3, label = "RMSD local chaine A3")
plt.plot(x,y_A4, label = "RMSD local chaine A4")
plt.title("RMSD local contre la conformation de reference")
plt.xlabel("Conformations")
plt.ylabel("RMSD")
plt.legend(loc = 4)
plt.savefig(args.outg)

##################################################################

print("Les resultats sont enregistrés dans le fichier \"{}\"".format(args.outb))
print("Le plot est enregistrés dans le fichier \"{}\"".format(args.outg))
print("\n")
print("Temps d'execution : {:.1f} secondes ---" .format(time.time() - start_time))
