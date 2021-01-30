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
This code is to make two plots: doublet fraction vs Ca plot and its phasediagram (in colormap).
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
phis.sort()

# make plots
r = 1
plt.figure(figsize = (16,12))
for phi in phis:
    Cas = list(parameter_set[phi].keys())
    Cas.sort()
    #if phi == 5.9986: tmp.remove(0.04)
    Cas = np.array(Cas)
    avg_df = np.zeros(len(Cas))

    for Ca_index, Ca in enumerate(Cas):
        ensemble_count = 0
        sum_df = 0

        # run over ensemble
        for ensemble_id in parameter_set[phi][Ca]:
            try:
                result = calcDoubletFraction(phi, Ca, 1, r, 0, ensemble_id, 1)
                sum_df += np.mean(result[0][0][-int(result[1]/2):]) # last half
                ensemble_count += 1
            except:
                print("Error: no preprocessed data in (phi = {}, Ca = {}, ensemble_id = {})".format(phi, Ca, ensemble_id))

        if (ensemble_count > 0):
            avg_df[Ca_index] = sum_df/ensemble_count
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))
    plt.plot(Cas[np.nonzero(avg_df)], avg_df[np.nonzero(avg_df)], label = phi)

plt.xticks(fontsize = 20)
plt.xlabel("Ca", fontsize = 30)
plt.yticks(fontsize = 20)
plt.ylabel("doublet fraction", fontsize = 30)
plt.title("Doublet Fraction vs Ca (Suspension system)\ncriteria_r = {}Dm, criteria_T = 1t_rot, second half time series".format(r), fontsize = 30)
plt.legend(fontsize = 20)
plt.savefig("./Pictures/SuspensionSystem_DoubletFraction_vs_Ca_EnsembleAveraged_SecondHalf_r_{}.png".format(r), dpi = 300)
plt.close()
print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))