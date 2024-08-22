# coding=utf-8
# '''
# Created on 19/03/2020
#
# @author: modelisation - TR
# '''

# Import the modules necessary to initiate the L-systems
import os
import openalea.lpy as opy
import openalea.plantgl as opal
from lgrass import meteo_ephem as meteo_ephem
from lgrass import param_reproduction_functions as prf
from lgrass import cuts as cuts
from lgrass import run_caribu_lgrass as run_caribu_lgrass
from lgrass import gen_lstring as gen_lstring
import pandas as pd
import numpy as np
import time


# préparation de l'ensemble des conditions/paramètres et exécution de lgrass
def runlsystem(plan_sim=None, id_scenario=0, id_gener=1):
    if plan_sim is None:
        raise NameError('Pas de plan de simulation chargé.')

    # Fichiers d'entrée
    genet_file = 'ped.r'
    # un fichier de remplacement du modèle génétique qui génère une population de C et détermine le nombre de plantes du couvert
    param_plant_file = 'liste_plantes.csv'

    # Répertoires de lecture/écriture
    INPUTS_DIRPATH = 'inputs'
    OUTPUTS_DIRPATH = 'outputs'
    GENET_DIRPATH = 'modelgenet'

    # Charger le plan de simulation et le lsystem
    row = plan_sim.iloc[id_scenario]
    name = str(row["name"])
    lpy_filename = './lgrass.lpy'
    lsystem = opy.Lsystem(lpy_filename)
    lsystem.name_sim = name

    # Choix du fichier de lecture du C en fonction de l'option de reproduction des plantes
    opt_repro = row["option_reproduction"]
    if opt_repro != "False":
        in_genet_file = os.path.join(GENET_DIRPATH, genet_file)
    else:
        in_genet_file = None
    in_param_file = os.path.join(INPUTS_DIRPATH, param_plant_file)

    # Parametres des plantes
    lsystem.ParamP, lsystem.nb_plantes, lsystem.NBlignes, lsystem.NBcolonnes, lsystem.posPlante, lsystem.Plantes, \
    lsystem.Genotypes, lsystem.flowering_model = prf.define_param(
        in_param_file=in_param_file, in_genet_file=in_genet_file, out_param_file=os.path.join(OUTPUTS_DIRPATH, name + '.csv'),
        id_gener=id_gener, opt_repro=opt_repro)


    # Parametres de simulation
    lsystem.option_tallage = row["option_tallage"]
    lsystem.option_senescence = row["option_senescence"]
    lsystem.option_floraison = row["option_floraison"]
    lsystem.option_tiller_regression = row["option_tiller_regression"]
    lsystem.option_morphogenetic_regulation_by_carbone = row["option_morphogenetic_regulation_by_carbone"]
    lsystem.derivationLength = int(row["derivationLength"])
    lsystem.sowing_date = row["sowing_date"]
    lsystem.site = row["site"]
    lsystem.meteo = meteo_ephem.import_meteo_data(row["meteo_path"], row['sowing_date'], row['site'])
    lsystem.output_induction_file_name = name + '_' + 'induction'
    lsystem.output_organ_lengths_file_name = name + '_' + 'organ_lengths'

    # Gestion des tontes
    opt_tontes = row["option_tontes"]
    if opt_tontes:
        lsystem.cutting_dates, lsystem.derivationLength = cuts.define_cutting_dates(lsystem.meteo,
                                                                                    int(row["derivationLength"]),
                                                                                    row["cutting_freq"])
    else:
        lsystem.cutting_dates = []

    # Gestion caribu
    opt_caribu = row["option_caribu"]
    if opt_caribu:
        dico_caribu = run_caribu_lgrass.init(meteo=lsystem.meteo, nb_plantes=lsystem.nb_plantes, scenario=row)
        lsystem.BiomProd = [0.] * lsystem.nb_plantes
        # Rédaction d'un fichier de sortie
        path_out = os.path.join(OUTPUTS_DIRPATH, name + '_caribu.csv')
        output = open(path_out, 'w')
        output.write("GDD;Date;Day;nb_talles;biomasse_aerienne;surface_foliaire;lstring" + "\n")

    # Lancement du lsystem
    lsystem.current_day = 1
    lstring = lsystem.axiom
    for dd in range(0, lsystem.derivationLength):
        day = lsystem.current_day
        lstring = lsystem.derive(lstring, dd, 1)
        lscene = lsystem.sceneInterpretation(lstring)
        if opt_caribu:
            # on exécute caribu une fois par jour
            if lsystem.current_day > day:
                # fonction d'application de caribu
                lsystem.BiomProd, dico_caribu['radiation_interception'], dico_caribu[
                    'Ray'] = run_caribu_lgrass.runcaribu(lstring, lscene, lsystem.current_day,
                                                         lsystem.tiller_appearance,
                                                         lsystem.nb_plantes, dico_caribu)
                # fichier de sortie de caribu
                output.write(";".join(
                    [str(lsystem.TPS), str(lsystem.sowing_date), str(lsystem.current_day), str(lsystem.nb_talle[0]),
                     str(lsystem.BiomProd[0]), str(lsystem.rapportS9_SSol_dict[0])]) + "\n")
        opal.all.Viewer.display(lscene)

    # Matrice de croisement des plantes
    if opt_repro != "False":
        mat = prf.create_seeds(lstring, lsystem.nb_plantes, opt_repro, row["cutting_freq"], lsystem.ParamP)
        np.savetxt(os.path.join(OUTPUTS_DIRPATH, name + "_mat.csv"), mat)
    else:
        mat = 0

    # Sauvegarder la lstring dans un répertoire pour pouvoir la charger dans une prochaine simulation
    if row['option_sauvegarde']:
        gen_lstring.save_lstring(lstring, lsystem)

    # Vider le lsystem
    lsystem.clear()
    print(''.join((name, " - done")))
    return mat


# Algorithme de reproduction des générations via le modèle génétique
def simpraise(plan_sim=None, id_scenario=0):
    if plan_sim is None:
        raise NameError('Pas de plan de simulation chargé.')
    row = plan_sim.iloc[id_scenario]

    # Config des fichiers d'entrée
    INPUTS_DIRPATH = 'inputs'
    src = os.path.join(INPUTS_DIRPATH, 'insim.txt')
    dst = './modelgenet'
    exe = 'simpraise.exe'

    # Génération des fondateurs, première exécution du modèle génétique
    prf.rungenet(src, dst, exe, None)

    # Boucle des générations
    for i in range(1, row['num_gener'] + 1):
        # fichiers de sortie associés à la ième génération
        plan_sim["name"] = row["name"] + "_G" + str(i)
        # modèle morpho et matrice de croisement
        mat = runlsystem(plan_sim=plan_sim, id_scenario=id_scenario, id_gener=i)
        # modèle génétique et paramètre C
        prf.rungenet(src, dst, exe, mat)
    return 0


if __name__ == '__main__':
    timing = time.time()
    plan = pd.read_csv("inputs/plan_simulation.csv", sep=',')

    # runlsystem(plan, 3, 1)
    simpraise(plan_sim=plan, id_scenario=0)
    print('Global execution time : ', time.time() - timing)
