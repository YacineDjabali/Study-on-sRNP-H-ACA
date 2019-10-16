#!/usr/bin/env python3

"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

"""

import func_distMD

def isResInterface(distance, seuil):
	"""
		Fonction qui donne l'information sur la presence d'un residu a l'interface ou non, a partir d'un seuil donné.
		Prend en parametre :
			- distance : la distance entre residues
			- seuil : seuil limite en dessous du quel on caracterise le residue comme faisant partie de l'interface
		Retourne : 
		 	- 1/0 (booleen) : 1 si fait partie de l'interface, 0 sinon 
	"""
	if(distance<seuil):
		return 1
	else:
		return 0
                
def frequency(occurence, NbConformation):
	"""
		Fonction qui donne la frequence d'appartenance a l'interface a partir du nombre d'occurences d'un residue obtenue sur plusieurs conformations.
		Prend en parametre :
			- occurence : Nombre de fois dont le residu se trouve a l'interface
			- seuil : Nombre de conformation
		Retourne : 
		 	-  (occurence / NbConformation) * 100 : Frequence calculé en %
	"""
	return  (occurence / NbConformation) * 100