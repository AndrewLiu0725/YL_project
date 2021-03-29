# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/29/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getSuspensionParameterSets
import os
import pickle

"""
This code is to plot ensemble averaged intrinsic viscosity time series for suspension system.
It also provides flags to save the ensemble averaged and individual time series.
"""

SAVE_ENSEMBLE_AVERAGED = 1
SAVE_INDIVIDUAL = 1
PLOT = 1
print("FLAGS: SAVE_ENSEMBLE_AVERAGED = {}, SAVE_INDIVIDUAL = {}, PLOT = {}".format(SAVE_ENSEMBLE_AVERAGED, SAVE_INDIVIDUAL, PLOT))

start_time = time.time()


# make plots
r = 1
max_timesteps = 10000

if SAVE_ENSEMBLE_AVERAGED:
    output_dict_EA = {} # dictionary to store data. Format: data[phi][Ca] = avg_df(t)

if SAVE_INDIVIDUAL:
    output_dict_I = {} # dictionary to store data. Format: data[phi][Ca] = [df(t)]

[phis, parameter_set] = getSuspensionParameterSets()

for phi in phis:
    for Ca in parameter_set[phi].keys():
        # initialize
        ensemble_count = 0
        iv_data = np.zeros((len(parameter_set[phi][Ca]), max_timesteps))
        min_timesteps = max_timesteps

        # run over ensemble
        for ensemble_id in parameter_set[phi][Ca]:
            try:
                iv_t = getIntrinsicViscosity(phi, Ca, max_timesteps, ensemble_id, 1)
                min_timesteps = min(min_timesteps, len(iv_t))
                iv_data[ensemble_count, :len(iv_t)] = iv_t[:]
                ensemble_count += 1

                if SAVE_INDIVIDUAL:
                    if phi in output_dict_I.keys():
                        if Ca in output_dict_I[phi].keys():
                            output_dict_I[phi][Ca].append(iv_t)
                        else:
                            output_dict_I[phi][Ca] = [iv_t]
                    else:
                        output_dict_I[phi] = {Ca: [iv_t]}

            except:
                pass
                #print("Error: no preprocessed data in (phi = {}, Ca = {}, ensemble id = {})".format(phi, Ca, ensemble_id))

        # make ensemble average time series
        if ensemble_count > 0:
            avg_iv_t = np.mean(df[:ensemble_count, :min_timesteps], axis=0)

            if SAVE_ENSEMBLE_AVERAGED:
                if phi in output_dict_EA.keys():
                    output_dict_EA[phi][Ca] = [avg_iv_t]
                else:
                    output_dict_EA[phi] = {Ca: [avg_iv_t]}

            if PLOT:
                std_iv_t = np.std(df[:ensemble_count, :min_timesteps], axis=0)
                slicing = 10 # python slicing, used to make plot more readable
                plt.figure(figsize = (16,12))
                plt.errorbar(list(range(min_timesteps))[::slicing], avg_iv_t[::slicing], yerr = std_iv_t[::slicing])
                plt.xticks(fontsize = 20)
                plt.xlabel("time", fontsize = 30)
                plt.yticks(fontsize = 20)
                plt.ylabel(r'$\eta _{int}$', fontsize = 30)
                plt.title("Intrinsic viscosity vs Time\nphi = {}, Ca = {}, {} ensembles".format(phi, Ca, ensemble_count), fontsize = 30)
                plt.savefig("./Pictures/SuspensionSystem_IntrinsicViscosity_vs_Time_EnsembleAveraged_phi_{}_Ca_{}.png".format(phi, Ca), dpi = 200)
                plt.close()
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))

if SAVE_ENSEMBLE_AVERAGED:
    with open("IV_Suspension_EnsembleAveraged.pickle", 'wb') as handle:
        pickle.dump(output_dict_EA, handle, protocol=pickle.HIGHEST_PROTOCOL)

if SAVE_INDIVIDUAL:
    with open("IV_Suspension_Indivisual.pickle", 'wb') as handle:
        pickle.dump(output_dict_I, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))