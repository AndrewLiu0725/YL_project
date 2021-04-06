# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 04/06/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import getIntrinsicViscosity
import os
import sys
import pickle

"""
This code is to plot angle averaged intrinsic viscosity time series for two-cell system.
It also provides flags to save the ensemble averaged and individual time series.
"""

SAVE_ENSEMBLE_AVERAGED = 1
SAVE_INDIVIDUAL = 1
PLOT = 1
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
st = 1

for phi in phi_range:
    for Ca in Ca_range:
        # initialize
        ensemble_count = 0
        data = np.zeros((len(angle_range), ncycle-st))

        # run over angle distribution
        for angle in angle_range:
            job = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
            if os.path.isfile(root_folder+"Data/"+job+"/data/bond0_t{}.vtk".format(int(3669*2*(int(ncycle/2)-1)))):
                try:
                    iv_t = getIntrinsicViscosity(phi, Ca, ncycle, angle, 0)[st:]
                    data[ensemble_count, :] = iv_t[:]
                    ensemble_count += 1

                    if SAVE_INDIVIDUAL:
                        if phi in output_dict_I.keys():
                            if Ca in output_dict_I[phi].keys():
                                output_dict_I[phi][Ca].append(iv_t)
                            else:
                                output_dict_I[phi][Ca] = [iv_t]
                        else:
                            output_dict_I[phi] = {Ca: [iv_t]}

                except KeyboardInterrupt:
                    print("Pressed ctrl+C\n")
                    sys.exit()

                except Exception as e:
                    print("Error: Case (phi = {}, Ca = {}, angle = {}): \n{}".format(phi, Ca, angle, e))

        # make angle averaged time series
        if ensemble_count > 0:
            avg_iv_t = np.mean(data[:ensemble_count, :], axis=0)

            if SAVE_ENSEMBLE_AVERAGED:
                if phi in output_dict_EA.keys():
                    output_dict_EA[phi][Ca] = [avg_iv_t]
                else:
                    output_dict_EA[phi] = {Ca: [avg_iv_t]}

            if PLOT:
                shear_p = 1026.4*Ca
                std_iv_t = np.std(data[:ensemble_count, :], axis=0)
                slicing = 40 # python slicing, used to make plot more readable. recommend value: 0.005*ncycle
                
                fig, ax1 = plt.subplots(figsize = (8, 6))
                ax1.errorbar(list(range(ncycle))[::slicing], avg_iv_t[::slicing], yerr = std_iv_t[::slicing])
                ax1.set_xlabel("time", fontsize = 15)
                ax1.set_ylabel(r'$\eta _{int}$', fontsize = 15)
                ax1.set_title("Intrinsic Viscosity vs Time\nphi = {}, Ca = {}, {} angles".format(phi, Ca, ensemble_count), fontsize = 20)

                # add second x-axis
                ax2 = ax1.twiny() # ax1 and ax2 share y-axis
                ax2.set_xlabel("non-dimensional time scale ({})".format(r'$\dot \gamma t$'), fontsize = 15)
                ax2.set_xticks(ax1.get_xticks())
                ax2.set_xbound(ax1.get_xbound())
                ax2.set_xticklabels([int(x * shear_p) for x in ax1.get_xticks()])
                
                fig.subplots_adjust(top=0.85)
                fig.tight_layout()
                plt.savefig("./Pictures/TwoCellSystem_IntrinsicViscosity_vs_Time_AngleAveraged_phi_{}_Ca_{}.png".format(phi, Ca), dpi = 200)
                plt.close()
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))

if SAVE_ENSEMBLE_AVERAGED:
    with open("IV_TwoCell_EnsembleAveraged.pickle", 'wb') as handle:
        pickle.dump(output_dict_EA, handle, protocol=pickle.HIGHEST_PROTOCOL)

if SAVE_INDIVIDUAL:
    with open("IV_TwoCell_Indivisual.pickle", 'wb') as handle:
        pickle.dump(output_dict_I, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))