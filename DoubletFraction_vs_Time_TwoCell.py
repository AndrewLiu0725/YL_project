# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/10/2022
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction
import os
import sys
import pickle

"""
This code is to plot angle averaged doublet fraction time series for two-cell system.
It also provides flags to save the ensemble averaged and individual time series.
"""

SAVE_ENSEMBLE_AVERAGED = 1
SAVE_INDIVIDUAL = 0
PLOT = 0 # make ensemble averaged doublet fraction vs time plot for each set of parameters
print("FLAGS: SAVE_ENSEMBLE_AVERAGED = {}, SAVE_INDIVIDUAL = {}, PLOT = {}".format(SAVE_ENSEMBLE_AVERAGED, SAVE_INDIVIDUAL, PLOT))

start_time = time.time()

root_folder = "/userdata4/ajliu/RBC_doublet/"

# setup
ncycle = 4000
vol = 746.3163
phi_range =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phi_range.append(float(str(phi*100)[:3]))
Ca_range = [i*0.01 for i in range(1, 21)]
angle_range = [90 - 10*j for j in range(18)]

if SAVE_ENSEMBLE_AVERAGED:
    output_dict_EA = {} # dictionary to store data. Format: data[phi][Ca] = avg_df(t)

if SAVE_INDIVIDUAL:
    output_dict_I = {} # dictionary to store data. Format: data[phi][Ca] = [df(t)]


# make plot
r = 1

for phi in phi_range:
    for Ca in Ca_range:
        # initialize
        ensemble_count = 0
        data = np.zeros((len(angle_range), ncycle))
        min_timesteps = ncycle

        # run over angle distribution
        for angle in angle_range:
            job = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
            if os.path.isfile(root_folder+"Data/"+job+"/data/bond0_t{}.vtk".format(int(3669*2*(int(ncycle/2)-1)))):
                try:
                    result = calcDoubletFraction(phi, Ca, 1, r, 0, angle, 0)
                    min_timesteps = min(min_timesteps, result[1])
                    data[ensemble_count, :result[1]] = result[0][0][:result[1]]
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

                except KeyboardInterrupt:
                    print("Pressed ctrl+C\n")
                    sys.exit()

                except Exception as e:
                    print("Error: Case (phi = {}, Ca = {}, angle = {}): \n{}".format(phi, Ca, angle, e))

        # make angle averaged time series
        if ensemble_count > 0:
            avg_df = np.mean(data[:ensemble_count, :min_timesteps], axis=0)

            if SAVE_ENSEMBLE_AVERAGED:
                if phi in output_dict_EA.keys():
                    output_dict_EA[phi][Ca] = [avg_df]
                else:
                    output_dict_EA[phi] = {Ca: [avg_df]}

            if PLOT:
                std_df = np.std(data[:ensemble_count, :min_timesteps], axis=0)
                slicing = 20 # python slicing, used to make plot more readable. recommend value: 0.005*ncycle
                
                fig, ax1 = plt.subplots(figsize = (8, 6))
                ax1.errorbar(list(range(min_timesteps))[::slicing], avg_df[::slicing], yerr = std_df[::slicing])
                ax1.set_xlabel("time ({})".format(r'$\dot \gamma t$'), fontsize = 15)
                ax1.set_ylabel("doublet fraction", fontsize = 15)
                ax1.set_title("Doublet Fraction vs Time\nphi = {}, Ca = {}, {} angles".format(phi, Ca, ensemble_count), fontsize = 20)

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
                plt.savefig("./Pictures/TwoCellSystem_DoubletFraction_vs_Time_AngleAveraged_phi_{}_Ca_{}.png".format(phi, Ca), dpi = 200)
                plt.close()
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))

if SAVE_ENSEMBLE_AVERAGED:
    with open("AverageDF_TwoCell_EnsembleAveraged.pickle", 'wb') as handle:
        pickle.dump(output_dict_EA, handle, protocol=pickle.HIGHEST_PROTOCOL)

if SAVE_INDIVIDUAL:
    with open("AverageDF_TwoCell_Indivisual.pickle", 'wb') as handle:
        pickle.dump(output_dict_I, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))