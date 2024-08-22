# coding: utf8
# '''
# Created on 20/04/2020
#
# @author: modelisation - TR
# '''
import pandas as pd


# Fonction déterminant les jours de coupes programmés en fonction de la fréquence des tontes en jour
def define_cutting_dates(weather, max_degrees, cutting_period, start_date=None):
    cutting_dates = []
    derivation_length = max_degrees
    degrees = 0
    do_cut = False
    if isinstance(start_date, str):
        i = weather.index[weather['date'] == pd.to_datetime(start_date, format='%Y_%m_%d')][0] - cutting_period
    else:
        i = 0
    cut = i
    while degrees <= max_degrees:
        degrees += weather['mean_temperature'].iloc[i]
        i += 1
        if i == cut + cutting_period:
            do_cut = True
        if do_cut & (degrees <= max_degrees):
            cut = i
            cutting_dates.append(cut)
            derivation_length += 3  # une tonte nécessite 3 itérations supplémentaires (x+0.25,x+0.5,x+0.75)
            do_cut = False
    print('Jours de tonte :', cutting_dates)
    print('Mise à jour de la durée de simulation (°j) :', derivation_length)
    return cutting_dates, derivation_length
