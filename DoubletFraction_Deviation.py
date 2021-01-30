# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/29/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction
import os

"""
This code is to calculate the deviation of each ensembled average doublet fraction time series.
"""

start_time = time.time()

path = "/raid6/ctliao/Data/HI_ordering/"

# parser
def parser(string):
    # output [phi, Ca, ensemble id]
    string = string.split("-")
    ensemble_id = int(string[1])
    # format: h24_phi4.4989_Re0.1_Ca0.06_WCA1_zero0.8-8
    if "_" in string[0]:
        string = string[0].split("_")
        phi = float(string[1][3:])
        Ca = float(string[3][2:])
    # format: h24phi4.9488Re0.1Ca0.06WCA1zero0.8-4
    else:
        phi = float(string[0][string[0].find('i')+1:string[0].find('R')])
        Ca = float(string[0][string[0].find('a')+1:string[0].find('W')])
    return [phi, Ca, ensemble_id]


parameter_set = {}
# create parameter_set (two layer dict)
for fn in os.listdir(path):
    if (fn[0] == "h") and (os.path.isdir(path+fn)):
        result = parser(fn)
        [phi, Ca, ensemble_id] = result
        if phi in parameter_set.keys():
            if Ca in parameter_set[phi].keys():
                parameter_set[phi][Ca].append(ensemble_id)
            else:
                parameter_set[phi][Ca] = [ensemble_id]
        else:
            parameter_set[phi] = {}
            parameter_set[phi][Ca] = [ensemble_id]


phis = []

number_of_parameter_set = 0
for phi in parameter_set.keys():
    if (len(parameter_set[phi].keys()) > 5):
        for Ca in parameter_set[phi].keys():
            number_of_parameter_set += len(parameter_set[phi][Ca])
        phis.append(phi)

# calculate deviation
r = 1

deviation = []
ps = []
max_timesteps = 3000

for phi in parameter_set.keys():
    if phi in phis:
        for Ca in parameter_set[phi].keys():
            # initialize
            ensemble_count = 0
            df = np.zeros((len(parameter_set[phi][Ca]), max_timesteps))
            min_timesteps = max_timesteps

            # run over ensemble
            for ensemble_id in parameter_set[phi][Ca]:
                try:
                    result = calcDoubletFraction(phi, Ca, 1, r, 0, ensemble_id, 1)
                    min_timesteps = min(min_timesteps, result[1])
                    df[ensemble_count, :result[1]] = result[0][0][:result[1]]
                    ensemble_count += 1
                except:
                    print("Error: no preprocessed data in (phi = {}, Ca = {}, ensemble id = {})".format(phi, Ca, ensemble_id))

            if (ensemble_count > 0):
                avg_df = np.mean(df[:ensemble_count, :min_timesteps], axis=0)
                first = np.mean(avg_df[500:1200])
                second = np.mean(avg_df[1200:min_timesteps])
                deviation.append(abs((second-first)/first))
                ps.append([phi, Ca])
            else:
                print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))
            

a = sorted(range(len(deviation)), key=lambda k: deviation[k], reverse = True)    
for i in a:
    print("({}, {}) deviation = {}".format(ps[i][0], ps[i][1], deviation[i]))     

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))