# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/10/2022
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getSuspensionParameterSets
import pickle

"""
This code is to plot ensemble averaged doublet fraction time series for suspension system.
It also provides flags to save the ensemble averaged and individual time series.
"""

SAVE_ENSEMBLE_AVERAGED = 1
SAVE_INDIVIDUAL = 0
PLOT = 0 # make ensemble averaged doublet fraction vs time plot for each set of parameters
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
        df = np.zeros((len(parameter_set[phi][Ca]), max_timesteps))
        min_timesteps = max_timesteps

        # run over ensemble
        for ensemble_id in parameter_set[phi][Ca]:
            try:
                result = calcDoubletFraction(phi, Ca, 1, r, 0, ensemble_id, 1)
                min_timesteps = min(min_timesteps, result[1])
                df[ensemble_count, :result[1]] = result[0][0][:result[1]]
                ensemble_count += 1

                df_t = result[0][0][:result[1]]
                if SAVE_INDIVIDUAL:
                    if phi in output_dict_I.keys():
                        if Ca in output_dict_I[phi].keys():
                            output_dict_I[phi][Ca].append(df_t)
                        else:
                            output_dict_I[phi][Ca] = [df_t]
                    else:
                        output_dict_I[phi] = {Ca: [df_t]}

            except:
                pass
                #print("Error: no preprocessed data in (phi = {}, Ca = {}, ensemble id = {})".format(phi, Ca, ensemble_id))

        # make ensemble average time series
        if ensemble_count > 0:
            avg_df = np.mean(df[:ensemble_count, :min_timesteps], axis=0)

            if SAVE_ENSEMBLE_AVERAGED:
                if phi in output_dict_EA.keys():
                    output_dict_EA[phi][Ca] = [avg_df]
                else:
                    output_dict_EA[phi] = {Ca: [avg_df]}

            if PLOT:
                std_df = np.std(df[:ensemble_count, :min_timesteps], axis=0)
                slicing = 20 # python slicing, used to make plot more readable
                fig, ax1 = plt.subplots(figsize = (8, 6))
                ax1.errorbar((np.arange(min_timesteps)[::slicing])*(4000/3669), avg_df[::slicing], yerr = std_df[::slicing])
                ax1.set_xlabel("time ({})".format(r'$\dot \gamma t$'), fontsize = 15)
                ax1.set_ylabel("doublet fraction", fontsize = 15)
                ax1.set_title("Doublet Fraction vs Time\nphi = {}, Ca = {}, {} ensembles".format(phi, Ca, ensemble_count), fontsize = 20)

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
                plt.savefig("./Pictures/SuspensionSystem_DoubletFraction_vs_Time_EnsembleAveraged_phi_{}_Ca_{}.png".format(phi, Ca), dpi = 200)
                plt.close()
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))

if SAVE_ENSEMBLE_AVERAGED:
    with open("AverageDF_Suspension_EnsembleAveraged.pickle", 'wb') as handle:
        pickle.dump(output_dict_EA, handle, protocol=pickle.HIGHEST_PROTOCOL)

if SAVE_INDIVIDUAL:
    with open("AverageDF_Suspension_Indivisual.pickle", 'wb') as handle:
        pickle.dump(output_dict_I, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))