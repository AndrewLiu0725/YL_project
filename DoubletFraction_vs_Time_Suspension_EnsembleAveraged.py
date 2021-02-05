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
This code is to plot ensemble averaged doublet fraction time series for suspension system.
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

for phi in parameter_set.keys():
    if (len(parameter_set[phi].keys()) > 5):
        phis.append(phi)

# make plots
r = 1
max_timesteps = 3000

for phi in phis:
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
                pass
                #print("Error: no preprocessed data in (phi = {}, Ca = {}, ensemble id = {})".format(phi, Ca, ensemble_id))

        # make ensemble average time series
        if ensemble_count > 0:
            avg_df = np.mean(df[:ensemble_count, :min_timesteps], axis=0)
            std_df = np.std(df[:ensemble_count, :min_timesteps], axis=0)
            slicing = 10 # python slicing, used to make plot more readable
            plt.figure(figsize = (16,12))
            plt.errorbar(list(range(min_timesteps))[::slicing], avg_df[::slicing], yerr = std_df[::slicing])
            plt.xticks(fontsize = 20)
            plt.xlabel("time", fontsize = 30)
            plt.yticks(fontsize = 20)
            plt.ylabel("doublet fraction", fontsize = 30)
            plt.title("Doublet Fraction vs Time\nphi = {}, Ca = {}, {} ensembles".format(phi, Ca, ensemble_count), fontsize = 30)
            plt.savefig("./Pictures/SuspensionSystem_DoubletFraction_vs_Time_EnsembleAveraged_phi_{}_Ca_{}.png".format(phi, Ca), dpi = 300)
            plt.close()
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))

print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))