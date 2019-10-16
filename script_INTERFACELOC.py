#!/usr/bin/env python3

"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

Propos :
Recupere le nombre de residues à l'interface par chaine,
et sort un graphique Nombre de residues à l'interface
en fonction des conformations.

"""

################# GESTION DES ARGUMENTS ######################

import argparse # Module de gestion des arguments

arg = argparse.ArgumentParser()
arg.add_argument("i", help="Chemin du fichier d'input") # ajout de l'argument input
arg.add_argument("outg", help = "Nom du fichier d'output des resultats graphiques (de preference .png)") # ajout de l'argument output pour les resultats graphiques

args = arg.parse_args() # recuperation des arguments
##############################################################


######################## MAIN ################################

import sys, time, func_parse, func_distMD, func_interface

start_time = time.time() # start du timer (pour le temps d'execution)

infile = args.i # recuperation de l'input
dPDB = func_parse.ParsePDBDynamique(infile) # creation du dictionnaire PDB

tmpDist = 0 # variable qui contiendra les distance entre residues

nbConformation = 0 # variable qui determinera le nombre de conformation dans le dictionnaire
interface = {} # dictionnaire qui contiendra les residues faisant partie de l'interface
interface["idProt"] = []

seuil = 9 # Seuil a partir du quel on considere un residu comme faisant partie de l'interface

nbResIntA1 = []
nbResIntA2 = []
nbResIntA3 = []        #Liste des points pour le plot
nbResIntA4 = []
abscisses = []
for i in dPDB["models"]:   #PARCOURS DU DICO PDB
    abscisses.append(nbConformation) # recuperation de liste de conformations
    nbConformation += 1
    resIntA1 = 0
    resIntA2 = 0  #compteurs de residues à l'interface
    resIntA3 = 0
    resIntA4 = 0
    for j in dPDB[i]["chains"]:
        if (j == "B" or j == "C"): # si on se trouve dans la chaine C (cation Metallique Zn) ou la chaine B (ARN) on break la boucle
            break
        else:
            for k in dPDB[i][j]["resNumber"]:
                for l in dPDB[i]["B"]["resNumber"]:
                    tmpDist =  func_distMD.distRes(dPDB[i][j][k], dPDB[i]["B"][l], "barycentre") # recuperation de la distance entre le residue k(proteine) et le residue l(ARN) par la methode du barycentre
                    if(func_interface.isResInterface(tmpDist,seuil)): #si tmpDist < seuil
                        if(j == "A1"):
                            resIntA1 += 1
                        elif(j == "A2"):
                            resIntA2 += 1                   # recupere le nombre de residues a l'interface suivant la chaine
                        elif(j == "A3"):
                            resIntA3 += 1
                        elif(j == "A4"):
                            resIntA4 += 1


    nbResIntA1.append(resIntA1)
    nbResIntA2.append(resIntA2)
    nbResIntA3.append(resIntA3)  # recuperation de nombre de residues suivant la conformation
    nbResIntA4.append(resIntA4)

##################################################################


###################### PREPARATION DES OUTPUTS ###################
import numpy as np
import matplotlib.pyplot as plt


x = np.array(abscisses)
y_A1 = np.array(nbResIntA1)
y_A2 = np.array(nbResIntA2)
y_A3 = np.array(nbResIntA3)
y_A4 = np.array(nbResIntA4)

plt.plot(x,y_A1, label = "Residue interface chaine A1")
plt.plot(x,y_A2, label = "Residue interface chaine A2")
plt.plot(x,y_A3, label = "Residue interface chaine A3")
plt.plot(x,y_A4, label = "Residue interface chaine A4")

plt.xlabel("Conformations")                                                #Creation du graphique
plt.ylabel("nombre de Residues a l'interface")
plt.legend(loc = 4)
plt.savefig(args.outg)

print("Le plot est enregistrés dans le fichier \"{}\"".format(args.outg))
print("\n")

print("Temps d'execution : {:.1f} secondes ---" .format(time.time() - start_time))
###################################################################
