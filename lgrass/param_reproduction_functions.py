# coding: utf8
# '''
# Created on 05/04/2020
#
# @author: modelisation - TR
# '''

import random
import numpy as np
import pandas as pd
import math
import os
import subprocess
from lgrass import flowering_functions as flowering_functions


# Créer la matrice de croisement des plantes pour établir une nouvelle génération de nb_plantes
def create_seeds(lstring, nb_plantes, talles, opt_repro, cutting_freq, ParamP):
    matrix = np.zeros((nb_plantes, nb_plantes))
    seeds = []
    elected_seeds = []

    # Méthode de calcul du nombre de graines via les regressions graines/tiges puis tiges/talles du rapport de Pierre Guinard (2012)
    if opt_repro == "SPPR_2012":
        nb_talles = sum(talles) + len(talles)  # talles contient le nombre de talles - 1 par plante
        if cutting_freq < 21:  # Coupes fréquentes
            nb_tiges = 0.32 * nb_talles - 9.70
            nb_seeds = int(9.73 * nb_tiges - 38.19)
        else:   # Coupes peu fréquentes
            nb_tiges = 0.22 * nb_talles - 4.73
            nb_seeds = int(6.35 * nb_tiges - 13.79)
        if nb_seeds >= 1:
            mothers = []
            for ind, i in enumerate(talles):
                for j in range(i+1):
                    mothers.append(ind)
            for i in range(nb_seeds):
                id_mother = random.Random().choice(mothers)
                fathers = [j for j in range(nb_plantes)]
                fathers.remove(id_mother)
                if fathers is not None:
                    id_father = random.Random().choice(fathers)  # séléction aléatoire du père
                    seeds.append((id_mother, id_father))
                else:
                    raise NameError("La génération ne comporte qu'une seule plante, la reproduction est impossible.")
        else:
            raise NameError("Il n'y a pas eu suffisamment de graines produites pour établir une nouvelle génération.")

    # Méthode de calcul via un nombre de graines par épillet
    elif opt_repro == "spikelets":
        # construction des graines et de leurs parents
        mothers = []
        for mod in lstring:
            if mod.name in ('apex',):
                if mod[0].final_spikelet_number is not None:
                    mothers.append((mod[0].id_plante, mod[0].final_spikelet_number))

        for k in range(len(mothers)):
            for i in range(int(mothers[k][1] * int(ParamP[mothers[k][0]]['seeds_by_spikelet']))):   # nb_epillets de la plante * nb_graines/épillet
                id_mother = mothers[k][0]  # séléction de la mère
                fathers = [j for j in range(nb_plantes)]
                fathers.remove(id_mother)
                if fathers is not None:
                    id_father = random.Random().choice(fathers)  # séléction aléatoire du père
                    seeds.append((id_mother, id_father))
                else:
                    raise NameError("La génération ne comporte qu'une seule plante, la reproduction est impossible.")

    # sélection aléatoire des graines pour la génération suivante et création de la matrice de croisement
    if len(seeds) < nb_plantes:
        raise NameError("Il n'y a pas eu suffisamment de graines produites pour établir une nouvelle génération.")

    for rand in range(nb_plantes):
        seed = random.Random().choice(seeds)
        elected_seeds.append(seed)
        seeds.remove(seed)
        # id_mere/ligne et id_pere/colonne
        matrix[elected_seeds[rand][0], elected_seeds[rand][1]] += 1
    return matrix


# Lire le fichier de données génétique (.r) et le formater en dataframe
def get_genet_file(in_genet_file=None):
    if in_genet_file is None:
        in_genet_file = 'modelgenet/ped.r'
    infile = pd.read_csv(in_genet_file, header=None)
    df = pd.DataFrame(columns=['geno', 'B', 'G', 'D', 'E', 'F', 'C', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])
    data = []
    for i in range(len(infile)):
        line = []
        k = 0
        while k < len(infile[0][i]):
            if infile[0][i][k] != ' ':
                chain = ''
                while infile[0][i][k] != ' ':
                    chain += infile[0][i][k]
                    k += 1
                    if k == len(infile[0][i]):
                        break
                line.append(chain)
            else:
                k += 1
        data.append(line)
    for i in range(len(data)):
        df.loc[i] = data[i]
    return df


# La formule pour convertir la donnée génétique en paramètre C
def calculate_C(x):
    # La version de cette fonction est incomplète (Thibault Raquet), voir Leopoldo Sanchez Rodriguez pour compléter la formule
    # se base sur une observation des sorties génétiques dont le paramètre variait entre -4 et 4, à convertir en [0.7, 1.6]
    return 0.1125 * x + 1.15


def calculate_premierCroiss(x):
    # Obtenir une plage de valeurs entre [80;120] à partir de [-4;4]
    pass


# Création du fichier de paramètres d'entrée pour chaque plante et configuration du planteur lgrass
def define_param(in_param_file='inputs/liste_plantes.csv', in_genet_file=None,
                 out_param_file='outputs/Simulation1.csv', id_gener=1, opt_repro=None, number_of_plants=None):
    if opt_repro is not None and opt_repro != "False":
        infile = get_genet_file(in_genet_file=in_genet_file)
        genet_data = infile.loc[infile['D'] == str(id_gener), :]
        if number_of_plants is not None and number_of_plants <= len(genet_data):
            genet_data = genet_data.iloc[:number_of_plants, :]
    data = pd.read_csv(in_param_file)
    # Création du fichier de paramètres d'entrée de chaque plante
    param_init = open(out_param_file, 'w')
    param_name = list(data.columns)
    for par in range(len(param_name)):
        L = [param_name[par]]
        for i in range(len(data)):
            L.append(str(data[param_name[par]].iloc[i]))
        if opt_repro is not None and opt_repro != "False":
            for _ in range(len(data), len(genet_data)):
                L.append(str(data[param_name[par]].iloc[len(data)-1]))
        param_init.write(";".join(L) + "\n")
    param_init.close()
    # lecture du fichier, creation de ParamP et des parametres de floraison
    param_plante = pd.read_csv(out_param_file, sep=";", header=None, index_col=0)
    # Conversion du paramètre génétique en valeur de C
    if opt_repro is not None and opt_repro != 'False':
        for i in range(len(genet_data)):
            param_plante.loc['C'].iloc[i] = calculate_C(float(genet_data['C'].iloc[i]))
            param_plante.loc['geno'].iloc[i] = int(genet_data['geno'].iloc[i])
        print('Infos valeurs de C générées :')
        print(param_plante.loc['C'].describe())
        param_plante.to_csv(out_param_file, sep=';')

    flowering_param = flowering_functions.FloweringFunctions()
    flowering_param.param.__dict__.update(
        (k, list(param_plante.loc[k, :])) for k in param_plante.index.intersection(flowering_param.param.__dict__))

    ParamP = list(dict(zip(param_plante.index, param_plante.iloc[:, col])) for col in range(len(param_plante.columns)))
    # Creation des matrices d'identifiant des plantes et de leur genotype
    nb_plantes = len(ParamP)
    NBlignes = int(math.ceil(np.sqrt(nb_plantes)))
    NBcolonnes = int(math.floor(np.sqrt(nb_plantes)))
    posPlante = [[i, j] for i, j in zip(sorted(list(range(NBlignes)) * NBcolonnes), list(range(NBcolonnes)) * NBlignes)] 
    # Plantes = np.arange(nb_plantes).reshape(NBlignes, NBcolonnes)
    # Genotypes = np.array([i for i in param_plante.loc['geno']]).reshape(NBlignes, NBcolonnes)
    Plantes = np.arange(nb_plantes)
    Genotypes = np.array([i for i in param_plante.loc['geno']])
    return ParamP, nb_plantes, NBlignes, NBcolonnes, posPlante, Plantes, Genotypes, flowering_param


# Remplir le fichier d'entree du modèle génétique avec la matrice de croisement et executer le modèle
def rungenet(src, dst, exe, mat, status):
    nb_plantes = 0
    nb_founder = 0

    # Bidouillage du fichier ped.r pour que le modèle génétique puisse supporter le changement du nombre de fondateurs
    # entre 2 générations
    if status == 1:
        ped = open(os.path.join(dst, 'ped.r'), 'r').readlines()
        new_ped = open(os.path.join(dst, 'ped.r'), 'w')
        for line in ped:
            row = list(line)
            row[47] = '1'
            new_ped.write(''.join(row))
        new_ped.close()

    with open(src, "r") as file:
        source = file.readlines()
        with open(os.path.join(dst, 'insim.txt'), "w") as destination:
            b = False
            for line in range(len(source)):
                if source[line] == "*num_cand \n":
                    nb_plantes = int(source[line+1])
                if b:
                    b = False
                    continue
                else:
                    destination.write(source[line])
                    if source[line] == "*status_gener \n":
                        destination.write(f"{status}\n")
                        b = True
                    elif source[line] == "*mnum&fnum \n":
                        if status == 1:
                            destination.write(f"{mat.shape[0] // 2},{mat.shape[1] // 2}\n")
                            b = True
                        else:
                            nb_founder = source[line+1].split(',')
                            nb_founder = int(nb_founder[0]) + int(nb_founder[1])
                    if source[line] == "*mating_design \n":
                        break
            if mat is None:  # Création aléatoire de la matrice de descendance
                mat = np.zeros((nb_founder, nb_founder))
                for _ in range(nb_plantes):
                    plantes = list(range(nb_founder))
                    mother = random.Random().choice(plantes)
                    plantes.remove(mother)
                    father = random.Random().choice(plantes)
                    # id_mere/ligne et id_pere/colonne
                    mat[mother, father] += 1
            for i in mat:
                for j in i:
                    # écriture de la matrice de croisement ligne par ligne
                    destination.write(str(int(j)) + '\n')
            destination.close()
            file.close()
    working_directory = os.getcwd()
    os.chdir(dst)
    if os.name == "posix":
        subprocess.call([os.path.join(dst, exe)])  
    else:  
        os.startfile(exe)
    os.chdir(working_directory)
