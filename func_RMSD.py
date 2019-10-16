#!/usr/bin/env python3

"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

"""

import func_distMD, math

def RMSD(sommeDistances,nombreDePaire): 
	"""
		Fonction qui calcule le RMSD entre deux conformations.
		Prend en parametre :
			- sommeDistance : la somme des distances au carré
			- nombreDePaire : le nombre de pairs d'atome pour lesquels la distance a été calculé
		Retourne : 
		 	- math.sqrt((1/nombreDePaire)*sommeDistances) : La RMSD
	"""

	return math.sqrt((1/nombreDePaire)*sommeDistances)




