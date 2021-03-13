# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/13/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getSuspensionParameterSets
import os
import pickle

"""
This code is to plot ensemble averaged doublet fraction time series for suspension system.
"""

SAVE = 1
PLOT = 0
print("FLAG: SAVE = {}, PLOT = {}".format(SAVE, PLOT))

start_time = time.time()


# make plots
r = 1
max_timesteps = 10000

if SAVE:
    output_dict = {} # dictionary to store data. Format: data[phi][Ca] = avg_df(t)

[phis, parameter_set] = getSuspensionParameterSets()

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
            if SAVE:
                if phi in output_dict.keys():
                    output_dict[phi][Ca] = avg_df
                else:
                    output_dict[phi] = {Ca: avg_df}
            if PLOT:
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

if SAVE:
    with open("AverageDF_Suspension.pickle", 'wb') as handle:
        pickle.dump(output_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))