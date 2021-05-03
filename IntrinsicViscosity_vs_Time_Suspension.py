# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 05/03/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import getIntrinsicViscosity, getSuspensionParameterSets
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
                iv_t = getIntrinsicViscosity(phi, Ca, max_timesteps, ensemble_id, 1)[1:]
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
            avg_iv_t = np.mean(iv_data[:ensemble_count, :min_timesteps], axis=0)

            if SAVE_ENSEMBLE_AVERAGED:
                if phi in output_dict_EA.keys():
                    output_dict_EA[phi][Ca] = [avg_iv_t]
                else:
                    output_dict_EA[phi] = {Ca: [avg_iv_t]}

            if PLOT:
                std_iv_t = np.std(iv_data[:ensemble_count, :min_timesteps], axis=0)
                slicing = 20 # python slicing, used to make plot more readable

                fig, ax1 = plt.subplots(figsize = (8, 6))
                ax1.errorbar((np.arange(min_timesteps)[::slicing])*(4000/3669), avg_iv_t[::slicing], yerr = std_iv_t[::slicing])
                ax1.set_xlabel("time ({})".format(r'$\dot \gamma t$'), fontsize = 15)
                ax1.set_ylabel(r'$\eta _{int}$', fontsize = 15)
                ax1.set_title("Intrinsic Viscosity vs Time\nphi = {}, Ca = {}, {} ensembles".format(phi, Ca, ensemble_count), fontsize = 20)

                '''
                # add second x-axis
                ax2 = ax1.twiny() # ax1 and ax2 share y-axis
                ax2.set_xlabel("non-dimensional time scale ({})".format(r'$\dot \gamma t$'), fontsize = 15)
                ax2.set_xticks(ax1.get_xticks())
                ax2.set_xbound(ax1.get_xbound())
                ax2.set_xticklabels([int(x * shear_p) for x in ax1.get_xticks()])
                '''
                
                fig.subplots_adjust(top=0.85)
                fig.tight_layout()
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