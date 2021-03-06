# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/13/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity, getSuspensionParameterSets
import os
from scipy.stats import linregress
import sys

"""
This code is to plot slope of doublet fraction vs intrinsic viscosity vs Ca for suspension system
"""

start_time = time.time()

[phis, parameter_set] = getSuspensionParameterSets()

phis.sort()

Cas = set(list(parameter_set[phis[0]].keys()))
for phi in phis:
    Cas = Cas.intersection(list(parameter_set[phi].keys()))
Cas = list(Cas)
Cas.sort()


# get the data matrix
df = np.ones((len(phis), len(Cas)))*(-1)
viscosity = np.ones((len(phis), len(Cas)))*(-1)
r = 1

for phi_index, phi in enumerate(phis):
    for Ca_index, Ca in enumerate(Cas):
        ensemble_count = 0
        sum_df, sum_iv = 0, 0
        for ensemble_id in parameter_set[phi][Ca]:
            try:
                result = calcDoubletFraction(phi, Ca, 1, r, 0, ensemble_id, 1)
                sum_df += np.mean(result[0][0][-int(result[1]/2):]) # second half
                iv = getIntrinsicViscosity(phi, Ca, None, ensemble_id, 1)[result[1]/2 : result[1]]
                sum_iv += np.mean(iv)
                ensemble_count += 1

            except KeyboardInterrupt:
                print("Pressed ctrl+C\n")
                sys.exit()

            except:
                pass
                #print("Error occured: (phi = {}, Ca = {}, ensemble_id = {})\n".format(phi, Ca, ensemble_id))
        
        if ensemble_count > 0:
            df[phi_index, Ca_index] = sum_df/ensemble_count
            viscosity[phi_index, Ca_index] = sum_iv/ensemble_count
        else: # the value of this cell remains -1
            print("Error: no data in (phi = {}, Ca = {})\n".format(phi, Ca))
print("Calculated the data matrices!\n")

# make plot
slope_list = []
for i, Ca in enumerate(Cas):
    x = viscosity[np.where(df[:, i] >= 0)[0], i]
    y = df[np.where(df[:, i] >= 0)[0], i]
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    slope_list.append(slope)
    fig = plt.figure(figsize = (8,6))
    plt.plot(x, y, 'o')
    plt.xlabel("Intrinsic viscosity")
    plt.ylabel("Doublet fraction")
    plt.title("Suspension System\nCa = {}".format(Ca))
    plt.savefig("./Pictures/SuspensionSystem_DoubletFraction_vs_IntrinsicViscosity_r_{}_Ca_{}.png".format(r, Ca), dpi = 300)
    plt.close()


fig = plt.figure(figsize = (8,6))
plt.plot(Cas, slope_list)
plt.xlabel("Ca", fontsize = 20)
plt.ylabel("Slope", fontsize = 20)
plt.title(r"$Slope\left( \frac{doublet\;fraction}{intrinsic\;viscosity}\right)$"+"vs Ca", fontsize = 25)
fig.tight_layout()
plt.savefig("./Pictures/SuspensionSystem_Slope_vs_Ca.png", dpi=300)
plt.close()

print('Total time elapsed = {}\n'.format(str(datetime.timedelta(seconds=time.time()-start_time))))