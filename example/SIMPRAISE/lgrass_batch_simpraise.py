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
from lgrass.output_data import CsvGenerator
from lgrass import cuts as cuts
from lgrass import run_caribu_lgrass as run_caribu_lgrass
from lgrass import gen_lstring as gen_lstring

import pandas as pd
import numpy as np
import time


# préparation de l'ensemble des conditions/paramètres et exécution de lgrass
def runlsystem(plan_sim=None, id_scenario=0, genet_path='modelgenet', id_gener=1, display=False):
    if plan_sim is None:
        raise NameError('Pas de plan de simulation chargé.')

    # Fichiers d'entrée
    genet_file = 'ped.r'
    # un fichier de remplacement du modèle génétique qui génère une population de C et détermine le nombre de plantes du couvert
    param_plant_file = 'liste_plantes.csv'

    # Répertoires de lecture/écriture
    current = os.path.dirname(__file__)
    INPUTS_DIRPATH = os.path.join(current, 'inputs')
    OUTPUTS_DIRPATH = os.path.join(current, 'outputs')

    # Charger le plan de simulation et le lsystem
    row = plan_sim.iloc[id_scenario]
    name = str(row["name"])
    lpy_filename = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),"lgrass", "lgrass.lpy")
    lsystem = opy.Lsystem(lpy_filename)
    lsystem.name_sim = name

    # Choix du fichier de lecture du C en fonction de l'option de reproduction des plantes
    opt_repro = row["option_reproduction"]
    in_genet_file = os.path.join(genet_path, genet_file) if opt_repro in ('spikelets', 'SPPR_2012') else None
    in_param_file = os.path.join(INPUTS_DIRPATH, param_plant_file)

    # Parametres des plantes
    lsystem.ParamP, lsystem.nb_plantes, lsystem.NBlignes, lsystem.NBcolonnes, lsystem.posPlante, \
    lsystem.Plantes, lsystem.Genotypes, lsystem.flowering_model = prf.define_param(
        in_param_file=in_param_file, in_genet_file=in_genet_file,
        out_param_file=os.path.join(OUTPUTS_DIRPATH, name + '.csv'), id_gener=id_gener, opt_repro=opt_repro)

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
                                                                                    row["cutting_freq"],
                                                                                    row["cutting_start"])
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
        if display:
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
        mat = prf.create_seeds(lstring, lsystem.nb_plantes, lsystem.nb_talle, opt_repro, row["cutting_freq"], lsystem.ParamP)
        np.savetxt(os.path.join(OUTPUTS_DIRPATH, name + "_mat.csv"), mat)
    else:
        mat = 0

    # Sauvegarder la lstring dans un répertoire pour pouvoir la charger dans une prochaine simulation
    if row['option_sauvegarde']:
        gen_lstring.save_lstring(lstring, lsystem)

    csv_generator = CsvGenerator(lstring, name, OUTPUTS_DIRPATH)
    csv_generator.metadata_to_csv(lsystem)
    csv_generator.leaves_to_csv()
    csv_generator.internodes_to_csv()
    csv_generator.apex_to_csv()

    # Vider le lsystem
    lsystem.clear()
    print(''.join((name, " - done")))
    return mat


# Algorithme de reproduction des générations via le modèle génétique
def simpraise(plan_sim=None, id_scenario=0, display_morpho=False):
    if plan_sim is None:
        raise NameError('Pas de plan de simulation chargé.')
    row = plan_sim.iloc[id_scenario]
    # Config des fichiers d'entrée
    INPUTS_DIRPATH = 'example/SIMPRAISE/inputs'
    src = os.path.join(INPUTS_DIRPATH, 'insim.txt')
    dst = f"modelgenet_{row['name']}"
    if os.name == "posix" :
        exe = 'simpraise'    
    else:
        exe = 'simpraise.exe'

    # TODO : Valable sous windows
    dst = os.path.join(os.path.dirname(__file__), dst)
    genetpath = os.path.join(os.path.dirname(__file__), 'modelgenet')
    os.system(f"mkdir {dst}")
    for f in os.listdir(genetpath):
        if os.name == "posix" :
            os.system(f'cp {os.path.join(genetpath, f)} {os.path.join(dst, f)}')
        else:
            os.system(f'copy {os.path.join(genetpath, f)} {os.path.join(dst, f)}')

    # Génération des fondateurs, première exécution du modèle génétique
    prf.rungenet(src, dst, exe, None, 0)

    # Boucle des générations
    for i in range(1, row['num_gener'] + 1):
        time.sleep(2)  # Laisser le temps au modèle génétique de mettre à jour ses fichiers
        # fichiers de sortie associés à la ième génération
        plan_sim.loc[id_scenario, "name"] = row['name'] + "_G" + str(i)
        # modèle morpho et matrice de croisement
        mat = runlsystem(plan_sim=plan_sim, id_scenario=id_scenario, genet_path=dst, id_gener=i, display=display_morpho)
        # modèle génétique et paramètre C
        prf.rungenet(src, dst, exe, mat, 1)

    return 0


if __name__ == '__main__':
    timing = time.time()
    plan = pd.read_csv("example/SIMPRAISE/inputs/plan_simulation.csv", sep=',')

    # runlsystem(plan_sim=plan, id_scenario=4, id_gener=1, display=False)
    for i in range(2, 4):
        simpraise(plan_sim=plan, id_scenario=i, display_morpho=False)
    print('Global execution time : ', time.time() - timing)
