# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/30/2021
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt
from RBC_Utilities import calcDoubletFraction, getInstrinsicViscosity, getRelativeViscosity
import time as time
import datetime as datetime

"""
This code is to plot intrinsic and relative viscosity vs Ca 
for both doublet and two-singlet state with all volume fraction.
Note that states are not seperated!
"""
start_time = time.time()

# Parameters
r = 1.0
ncycle = 2000

vol = 746.3163
phis  =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phis.append(float(str(phi*100)[:3]))
phi_range = np.array(phis)
Ca_range = np.array([i*0.01 for i in range(1, 21)])


# Store the calculated data
# (iv/rv, mean/std, phi, Ca)
data = np.zeros((2, 2, len(phi_range), len(Ca_range)))

for phi_index, phi in enumerate(phi_range):
    print(phi)
    for Ca_index, Ca in enumerate(Ca_range):
        ## Store stat for each set of phi and Ca: x/x^2/count, d/k, iv/rv
        for j in range(18):
            angle = 90-10*j
            iv = getInstrinsicViscosity(phi, Ca, 0, ncycle, angle, 0)
            rv = getRelativeViscosity(phi, Ca, 0, ncycle, angle, 0)

            ## store the value       
            data[0, 0, phi_index, Ca_index] = np.mean(iv)
            data[0, 1, phi_index, Ca_index] = np.std(iv)
            data[1, 0, phi_index, Ca_index] = np.mean(rv)
            data[1, 1, phi_index, Ca_index] = np.std(rv)

### Plot here

### vs Ca
import matplotlib.colors as mcolors
import random
cmap = random.sample(list(mcolors.CSS4_COLORS.keys()), len(phi_range))

for i in range(2): # run over iv and rv
    plt.figure(figsize=(12,9))
    for phi_index, phi in enumerate(phi_range):
        plt.plot(Ca_range, data[i, 0, phi_index, :], label = 'phi = {}'.format(phi), color = cmap[phi_index])   
            
    plt.xlabel('Ca', fontsize = 20)
    if i == 0:
        text = "Intrinsic Viscosity"
    elif i == 1:
        text = "Relative viscosity"
    plt.ylabel(text, fontsize = 20)
    plt.title(text+' vs Ca', fontsize=30)
    plt.legend()
    plt.savefig('./Pictures/TwoCellSystem_{}_Ca_summary_r_{}_unseperated.png'.format("".join(text.split()), r), dpi = 300)


### vs phi
for Ca_index, Ca in enumerate(Ca_range):
    fig, ax1 = plt.subplots(figsize = (9, 9))
    color = 'tab:red'
    ax1.set_xlabel(r'$\phi$', fontsize = 20)
    ax1.set_ylabel("Intrinsic Viscosity", color = color, fontsize = 20)
    ax1.set_ylim([np.min(data[0, 0, :, :])*0.9, np.max(data[0, 0, :, :])*1.1])
    ax1.set_title("Ca = {}".format(Ca), fontsize = 30)
    ax1.tick_params(axis='y', labelcolor = color, labelsize = 20)
    ax1.tick_params(axis='x', labelsize = 20)
    p1, = ax1.plot(0.01*phi_range, data[0, 0, :, Ca_index], color = color, label = "Intrinsic Viscosity")

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel("Relative Viscosity", color = color, fontsize = 20)
    ax2.set_ylim([np.min(data[1, 0, :, :])*0.9, np.max(data[1, 0, :, :])*1.1])
    p2, = ax2.plot(0.01*phi_range, data[1, 0, :, Ca_index], color = color, label = "Relative Viscosity")
    ax2.tick_params(axis = 'y', labelcolor = color, labelsize = 20)

    lines = [p1, p2]
    ax1.legend(lines, [l.get_label() for l in lines], bbox_to_anchor=(1.0,1.2))
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.savefig('./Pictures/TwoCellSystem_Viscosity_vs_phi_Ca_{}_fixed_ylim_r_{}_unseperated.png'.format(Ca, r), dpi = 100)
    ax2.set_ylim([np.min(data[1, 0, :, Ca_index]), np.max(data[1, 0, :, Ca_index])])
    ax1.set_ylim([np.min(data[0, 0, :, Ca_index]), np.max(data[0, 0, :, Ca_index])])
    fig.savefig('./Pictures/TwoCellSystem_Viscosity_vs_phi_Ca_{}_r_{}__unseperated.png'.format(Ca, r), dpi = 300)

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
