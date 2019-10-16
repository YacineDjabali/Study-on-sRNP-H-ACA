#!/usr/bin/env python3
"""
Auteurs : DJABALI Yacine & PARMENTIER Raphael
Sujet : ETUDE DU COMPLEXE sRNP H /ACA

"""

import math

def distanceCA(dicoResI, dicoResII): 
    """
        Fonction qui calcule la distance entre les Carbones-Alpha de deux residues.
        Pour cela les coordonnées x,y,z de l'atome sont recuperer et la distance calculée.
        Prend en parametre :
            - dicoResI : le dictionnaire du premier residue
            - dicoResII : le dicionnaire du deuxieme residue
        Retourne :
            - math.sqrt((X**2 + Y**2 + Z**2)): la distance entre les deux residues
    """
    X = dicoResI["CA"]["x"] - dicoResII["CA"]["x"]
    Y = dicoResI["CA"]["y"] -  dicoResII["CA"]["y"]
    Z = dicoResI["CA"]["z"] - dicoResII["CA"]["z"]
	
    return math.sqrt((X**2 + Y**2 + Z**2))


def barycentreRes(dicoRes): 
    """
        Fonction qui calcule le barycentre d'un residue
        Pour cela les coordonnées x,y,z de tout les atomes du residue sont sommée un a un et divisé par le nombre de coordonnées.
        Prend en parametre :
            - dicoRes : Le dictionnaire du residue
        Retourne :
            - barycentre : Une liste qui contient les coordonnée x,y,z du barycentre du residue
    """
    nb_point = 0
    x_bar = 0
    y_bar = 0
    z_bar = 0
    for i in dicoRes["atomName"]:
        x_bar += dicoRes[i]["x"]
        y_bar += dicoRes[i]["y"] 
        z_bar += dicoRes[i]["z"]
        nb_point += 1
    baryCentre = [x_bar/nb_point, y_bar/nb_point, z_bar/nb_point]

    return baryCentre


def distRes(dicoResI, dicoResII, mode):


    """
        Fonction qui calcule la distance entre deux residue a partir du barycentre ou de la distance minimale
        Prend en parametre :
            - dicoResI : le dictionnaire du premier residue
            - dicoResII : le dicionnaire du deuxieme residue
            - mode : la methode de calcule de distance
        Retourne :
            - dist : la distance qui a été calculé (suivant le mode choisie)
    """

    import math

    if mode == "min" :
        dist = []
        for i in dicoResI["atomName"]:
            x_resI = dicoResI[i]["x"]
            y_resI = dicoResI[i]["y"]
            z_resI = dicoResI[i]["z"]
            for j in dicoResII["atomName"]:
                x_resII = dicoResII[j]["x"]
                y_resII = dicoResII[j]["y"]
                z_resII = dicoResII[j]["z"]
                dist.append(math.sqrt((x_resII-x_resI)**2 + (y_resII-y_resI)**2 + (z_resII-z_resI)**2))
        return min(dist)

    if mode == "barycentre":
        dist = 0
        barycentreResI = barycentreRes(dicoResI)
        barycentreResII = barycentreRes(dicoResII)
        dist = math.sqrt((barycentreResII[0]-barycentreResI[0])**2+(barycentreResII[1]-barycentreResI[1])**2+(barycentreResII[2]-barycentreResI[2])**2)
        return dist
