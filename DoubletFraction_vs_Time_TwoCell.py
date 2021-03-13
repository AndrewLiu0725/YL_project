# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/13/2021
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
"""

SAVE = 1
PLOT = 0
print("FLAG: SAVE = {}, PLOT = {}".format(SAVE, PLOT))

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

if SAVE:
    output_dict = {} # dictionary to store data. Format: data[phi][Ca] = avg_df(t)


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

                except KeyboardInterrupt:
                    print("Pressed ctrl+C\n")
                    sys.exit()

                except Exception as e:
                    print("Error: Case (phi = {}, Ca = {}, angle = {}): \n{}".format(phi, Ca, angle, e))

        # make angle averaged time series
        if ensemble_count > 0:
            avg_df = np.mean(data[:ensemble_count, :min_timesteps], axis=0)
            if SAVE:
                if phi in output_dict.keys():
                    output_dict[phi][Ca] = avg_df
                else:
                    output_dict[phi] = {Ca: avg_df}
            if PLOT:
                std_df = np.std(data[:ensemble_count, :min_timesteps], axis=0)
                slicing = 20 # python slicing, used to make plot more readable. recommend value: 0.005*ncycle
                plt.figure(figsize = (16,12))
                plt.errorbar(list(range(min_timesteps))[::slicing], avg_df[::slicing], yerr = std_df[::slicing])
                plt.xticks(fontsize = 20)
                plt.xlabel("time", fontsize = 30)
                plt.yticks(fontsize = 20)
                plt.ylabel("doublet fraction", fontsize = 30)
                plt.title("Doublet Fraction vs Time\nphi = {}, Ca = {}, {} angles".format(phi, Ca, ensemble_count), fontsize = 30)
                plt.savefig("./Pictures/TwoCellSystem_DoubletFraction_vs_Time_AngleAveraged_phi_{}_Ca_{}.png".format(phi, Ca), dpi = 300)
                plt.close()
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))

if SAVE:
    with open("AverageDF_TwoCell.pickle", 'wb') as handle:
        pickle.dump(output_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))