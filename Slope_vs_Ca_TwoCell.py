# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/30/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity, getRelativeViscosity
from scipy.stats import linregress

"""
This code is to plot slope of doublet fraction vs intrinsic viscosity vs Ca for two-cell system
"""

start_time = time.time()

# set up
ncycle = 2000
vol = 746.3163
phis  =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phis.append(float(str(phi*100)[:3]))
Ca_list = [(i+1)*0.01 for i in range(20)]
angles = [90-10*i for i in range(18)]
r = 1

# get doublet fraction
df = np.zeros((len(phis), len(Ca_list)))
viscosity = np.zeros((len(phis), len(Ca_list)))

for phi_index, phi in enumerate(phis):
    print('phi = ', phi)
    for Ca_index, Ca in enumerate(Ca_list):
        for angle_index, angle in enumerate(angles):
            result = calcDoubletFraction(phi, Ca, 1, r, 0, angle, 0)
            iv = getIntrinsicViscosity(phi, Ca, int(ncycle/2), ncycle, angle, 0)
            df[phi_index, Ca_index] += np.mean(result[0][0][-int(result[1]/2):])
            viscosity[phi_index, Ca_index] += np.mean(iv)
df = df/len(angles) # average over the ensemble
viscosity = viscosity/len(angles)
print("Calculated the data matrices!\n")


# plot slope here
slope_list = []
for i, Ca in enumerate(Ca_list):
    x = viscosity[:, i]
    y = df[:, i]
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    slope_list.append(slope)
    fig = plt.figure(figsize = (8,6))
    plt.plot(x, y, 'o')
    plt.xlabel("Intrinsic viscosity")
    plt.ylabel("Doublet fraction")
    plt.title("Two Cell System\nCa = {}".format(Ca))
    plt.savefig("./Pictures/TwoCellSystem_DoubletFraction_vs_IntrinsicViscosity_r_{}_Ca_{}.png".format(r, Ca), dpi = 300)


fig = plt.figure(figsize = (8,6))
plt.plot(Ca_list, slope_list)
plt.xlabel("Ca", fontsize = 20)
plt.ylabel("Slope", fontsize = 20)
plt.title(r"$Slope\left( \frac{doublet\;fraction}{intrinsic\;viscosity}\right)$"+"vs Ca", fontsize = 25)
fig.tight_layout()
plt.savefig("./Pictures/TwoCellSystem_Slope_vs_Ca.png", dpi=300)

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))