#!/usr/bin/env python3

"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

"""

import string, sys


def ParsePDBDynamique(infile): 
	"""
		Fonction dont le role est de parser la contenance d'un fichier pdb de plusieurs conformations.
		La methode consiste a traiter le fichier ligne par ligne et a ranger les informations utiles,
		dans un dicitonnaire.
		Prend en parametre : 
			- infile : le chemin du fichier pdb
		Retourne :
			- dicoPDB : Dictionnaire contenant le pdb parsé
	"""


	file = open(infile, "r")
	lines = file.readlines()
	file.close()


	dicoPDB = {}
	dicoPDB["models"] = []


	for line in lines :

		if line[0:5] == "MODEL":
			mod = line[9:14].strip()
			if not mod in dicoPDB["models"] : 
				dicoPDB["models"].append(mod)
				dicoPDB[mod] = {}
				dicoPDB[mod]["chains"] = []

		if line[0:4] == "ATOM" :
			chain = line[72:76].strip()
			if not chain in dicoPDB[mod]["chains"] :
				dicoPDB[mod]["chains"].append(chain)
				dicoPDB[mod][chain] = {}
				dicoPDB[mod][chain]["resNumber"] = []

			res = "%s"%(line[22:27]).strip()
			if not res in dicoPDB[mod][chain]["resNumber"] :
				dicoPDB[mod][chain]["resNumber"].append(res)
				dicoPDB[mod][chain][res] = {}
				dicoPDB[mod][chain][res]["resName"] = line[16:21].strip()

				dicoPDB[mod][chain][res]["atomName"] = []

			atom = line[12:16].strip()
			dicoPDB[mod][chain][res]["atomName"].append(atom)
			dicoPDB[mod][chain][res][atom] = {}
			dicoPDB[mod][chain][res][atom]["x"] = float(line[30:38])
			dicoPDB[mod][chain][res][atom]["y"] = float(line[38:46])
			dicoPDB[mod][chain][res][atom]["z"] = float(line[46:54])
			dicoPDB[mod][chain][res][atom]["id"] = line[5:11].strip()

	return dicoPDB



